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
  return false;
}

pair<int, int> calculate_insert_size(Alns &alns, int expected_ins_avg, int expected_ins_stddev, int max_expected_ins_size) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  ProgressBar progbar(alns.size(), "Processing alignments to compute insert size");
  Aln *prev_aln = nullptr;
  int64_t prev_cid = -1;
  int64_t num_valid_pairs = 0;
  int64_t num_overlap_rejected = 0;
  int64_t sum_insert_size = 0;
  int64_t num_large = 0;
  int min_insert_size = 1000000;
  int max_insert_size = 0;
  vector<int> insert_sizes;
  insert_sizes.reserve(alns.size());
  for (auto &aln : alns) {
    progbar.update();
    assert(aln.cstop <= aln.clen);
    if (prev_aln) {
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
            num_large++;
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
  // FIXME: now process the reads and find all merged reads and use their lengths as insert sizes. Actually, this could
  // also be done in merge reads, I guess... but simpler and cleaner to do it here, rather than passing that vector of
  // insert sizes through main()

  auto all_num_overlap_rejected = reduce_one(num_overlap_rejected, op_fast_add, 0).wait();
  auto all_num_valid_pairs = reduce_all(num_valid_pairs, op_fast_add).wait();
  auto all_num_alns = reduce_one(alns.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(all_num_valid_pairs, all_num_alns), " pairs with valid alignments to the same contig\n");
  SLOG_VERBOSE("Rejected ", perc_str(all_num_overlap_rejected, all_num_alns), " possible candidates with bad overlaps\n");
  auto all_num_large = reduce_one(num_large, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(all_num_large, all_num_alns), " outliers\n");

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

