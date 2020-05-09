/*
 Calculate histogram of insert sizes
 c shofmeyr@lbl.gov
 March 2020
*/

#include <iostream>
#include <math.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"

#include "utils.hpp"
#include "alignments.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;


bool bad_alignment(Aln *aln) {
  if (aln->rstart > KLIGN_UNALIGNED_THRES) return true;
  if (aln->rlen - aln->rstop > KLIGN_UNALIGNED_THRES) return true;
  if (aln->score1 < aln->rlen - HISTOGRAMMER_BAD_ALN) return true;
  return false;
}

int calculate_unmerged_rlen(Alns &alns) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  // get the unmerged read length - most common read length
  HASH_TABLE<int, int64_t> rlens;
  int64_t sum_rlens = 0;
  for (auto &aln : alns) {
    rlens[aln.rlen]++;
    sum_rlens += aln.rlen;
  }
  auto all_sum_rlens = upcxx::reduce_all(sum_rlens, op_fast_add).wait();
  auto all_nalns = upcxx::reduce_all(alns.size(), op_fast_add).wait();
  auto avg_rlen = all_sum_rlens / all_nalns;
  int most_common_rlen = avg_rlen;
  int64_t max_count = 0;
  for (auto &rlen : rlens) {
    if (rlen.second > max_count) {
      max_count = rlen.second;
      most_common_rlen = rlen.first;
    }
  }
  SLOG_VERBOSE("Computed unmerged read length as ", most_common_rlen, " with a count of ", max_count, " and average of ",
               avg_rlen, "\n");
  return most_common_rlen;
}

pair<int, int> calculate_insert_size(Alns &alns, int expected_ins_avg, int expected_ins_stddev, int max_expected_ins_size,
                                     const string &dump_msa_fname="") {
  BarrierTimer timer(__FILEFUNC__, false, true);
  auto unmerged_rlen = calculate_unmerged_rlen(alns);
  ProgressBar progbar(alns.size(), "Processing alignments to compute insert size");
  Aln *prev_aln = nullptr;
  int64_t prev_cid = -1;
  int64_t num_valid_pairs = 0;
  int64_t num_overlap_rejected = 0;
  int64_t sum_insert_size = 0;
  int64_t num_large = 0;
  int64_t num_small = 0;
  int min_insert_size = 1000000;
  int max_insert_size = 0;
  vector<string> misassembled_ctgs;
  vector<int> insert_sizes;
  insert_sizes.reserve(alns.size());
  for (auto &aln : alns) {
    progbar.update();
    assert(aln.cstop <= aln.clen);
    // if the read length is greater than the median read length, this must have been merged
    if (aln.rlen > unmerged_rlen) {
      num_valid_pairs++;
      // for merged reads, the merged length is the insert size
      auto insert_size = aln.rlen;
      min_insert_size = min(min_insert_size, insert_size);
      max_insert_size = max(max_insert_size, insert_size);
      sum_insert_size += insert_size;
      insert_sizes.push_back(insert_size);
    } else if (prev_aln) {
      auto read_id = substr_view(aln.read_id, 0, aln.read_id.length() - 2);
      char pair_num = aln.read_id[aln.read_id.length() - 1];
      auto prev_read_id = substr_view(prev_aln->read_id, 0, prev_aln->read_id.length() - 2);
      char prev_pair_num = prev_aln->read_id[prev_aln->read_id.length() - 1];
      if (read_id == prev_read_id && prev_aln->cid == aln.cid) {
        assert(prev_pair_num == '1' && pair_num == '2');
        if (bad_alignment(prev_aln) || bad_alignment(&aln)) {
          num_overlap_rejected++;
        } else {
          int prev_cstart = prev_aln->cstart;
          int prev_cstop = prev_aln->cstop;
          // if both alns have the same orientation, flip one of them
          if (prev_aln->orient == aln.orient) {
            prev_cstart = prev_aln->clen - prev_aln->cstop;
            prev_cstop = prev_aln->clen - prev_aln->cstart;
          }
          auto insert_size = max(prev_cstop, aln.cstop) - min(prev_cstart, aln.cstart);
          if (max_expected_ins_size && insert_size >= max_expected_ins_size) {
            //WARN("large insert size: ", insert_size, " prev: ", prev_aln->clen, " ",
            //     prev_aln->orient, " ", prev_aln->cstart, " ", prev_aln->cstop, " ", prev_cstart, " ", prev_cstop,
            //     " current ", aln.clen, " ", aln.orient, " ", aln.cstart, " ", aln.cstop);
            if (!dump_msa_fname.empty())
              misassembled_ctgs.push_back(to_string(prev_aln->cid) + " " + to_string(min(prev_cstart, aln.cstart)) + " " +
                                          to_string(max(prev_cstop, aln.cstop)));
            num_large++;
          } else if (insert_size < 10) {
            if (!dump_msa_fname.empty())
              misassembled_ctgs.push_back(to_string(prev_aln->cid) + " " + to_string(min(prev_cstart, aln.cstart)) + " " +
                                          to_string(max(prev_cstop, aln.cstop)));
            num_small++;
          } else {
            num_valid_pairs++;
            min_insert_size = min(min_insert_size, insert_size);
            max_insert_size = max(max_insert_size, insert_size);
            sum_insert_size += insert_size;
            insert_sizes.push_back(insert_size);
          }
        }
      }
    }
    prev_aln = &aln;
  }
  progbar.done();

  if (!dump_msa_fname.empty()) {
    string out_str = "";
    for (auto msa_ctg : misassembled_ctgs) {
      out_str += msa_ctg + "\n";
    }
    upcxx::atomic_domain<size_t> ad({upcxx::atomic_op::fetch_add, upcxx::atomic_op::load});
    upcxx::global_ptr<size_t> fpos = nullptr;
    if (!upcxx::rank_me()) fpos = upcxx::new_<size_t>(0);
    fpos = upcxx::broadcast(fpos, 0).wait();
    auto sz = out_str.length();
    size_t my_fpos = ad.fetch_add(fpos, sz, std::memory_order_relaxed).wait();
    // wait until all ranks have updated the global counter
    upcxx::barrier();
    int fileno = -1;
    size_t fsize = 0;
    if (!upcxx::rank_me()) {
      fsize = ad.load(fpos, std::memory_order_relaxed).wait();
      // rank 0 creates the file and truncates it to the correct length
      fileno = open(dump_msa_fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
      if (fileno == -1) WARN("Error trying to create file ", dump_msa_fname, ": ", strerror(errno), "\n");
      if (ftruncate(fileno, fsize) == -1) WARN("Could not truncate ", dump_msa_fname, " to ", fsize, " bytes\n");
    }
    upcxx::barrier();
    ad.destroy();
    // wait until rank 0 has finished setting up the file
    if (rank_me()) fileno = open(dump_msa_fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) WARN("Error trying to open file ", dump_msa_fname, ": ", strerror(errno), "\n");
    auto bytes_written = pwrite(fileno, out_str.c_str(), sz, my_fpos);
    close(fileno);
    if (bytes_written != sz) DIE("Could not write all ", sz, " bytes; only wrote ", bytes_written, "\n");
    upcxx::barrier();
    auto tot_bytes_written = upcxx::reduce_one(bytes_written, upcxx::op_fast_add, 0).wait();
    upcxx::barrier();
    SLOG_VERBOSE("Successfully wrote ", get_size_str(tot_bytes_written), " bytes to ", dump_msa_fname, "\n");
  }

  auto all_num_overlap_rejected = reduce_one(num_overlap_rejected, op_fast_add, 0).wait();
  auto all_num_valid_pairs = reduce_all(num_valid_pairs, op_fast_add).wait();
  auto all_num_alns = reduce_one(alns.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(all_num_valid_pairs, all_num_alns), " pairs with valid alignments to the same contig\n");
  SLOG_VERBOSE("Rejected ", perc_str(all_num_overlap_rejected, all_num_alns), " possible candidates with bad overlaps\n");
  auto all_num_large = reduce_one(num_large, op_fast_add, 0).wait();
  auto all_num_small = reduce_one(num_small, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(all_num_large, all_num_alns), " large (> ", max_expected_ins_size, ") and ",
               perc_str(all_num_small, all_num_alns), " small outliers\n");

  auto all_sum_insert_size = reduce_all(sum_insert_size, op_fast_add).wait();
  if (!all_num_valid_pairs) {
    if (expected_ins_avg) {
      WARN("Could not find any suitable alignments for calculating the insert size. Using the parameters ",
           expected_ins_avg, " ", expected_ins_stddev);
      return {expected_ins_avg, expected_ins_stddev};
    } else {
      DIE("Could not find any suitable alignments for calculating the insert size and no parameters are set.");
    }
  }
  auto insert_avg = all_sum_insert_size / all_num_valid_pairs;
  double sum_sqs = 0;
  for (auto insert_size : insert_sizes) {
    sum_sqs += pow((double)insert_size - insert_avg, 2.0);
  }
  auto all_min_insert_size = reduce_one(min_insert_size, op_fast_min, 0).wait();
  auto all_max_insert_size = reduce_one(max_insert_size, op_fast_max, 0).wait();
  auto all_sum_sqs = reduce_all(sum_sqs, op_fast_add).wait();
  auto insert_stddev = sqrt(all_sum_sqs / all_num_valid_pairs);
  SLOG_VERBOSE("Calculated insert size: average ", insert_avg, " stddev ", insert_stddev,
               " min ", all_min_insert_size, " max ", all_max_insert_size, "\n");
  if (expected_ins_avg) {
    if (abs(insert_avg - expected_ins_avg) > 100)
      SWARN("Large difference in calculated vs expected insert average sizes (expected ", expected_ins_avg, ")");
    if (abs(insert_stddev - expected_ins_stddev) > 100)
      SWARN("Large difference in calculated vs expected insert std dev (expected ", expected_ins_stddev, ")");
    SLOG_VERBOSE("Using specified ", expected_ins_avg, " avg insert size and ", expected_ins_stddev, " stddev\n");
    return {expected_ins_avg, expected_ins_stddev};
  }
  return {insert_avg, insert_stddev};
}

