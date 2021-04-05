/*
 HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from the U.S. Dept. of Energy).  All rights reserved."

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 (1) Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 (2) Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 (3) Neither the name of the University of California, Lawrence Berkeley National
 Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 endorse or promote products derived from this software without specific prior
 written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 to the features, functionality or performance of the source code ("Enhancements") to
 anyone; however, if you choose to make your Enhancements available either publicly,
 or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 written license agreement for such Enhancements, then you hereby grant the following
 license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 prepare derivative works, incorporate into other computer software, distribute, and
 sublicense such enhancements or derivative works thereof, in binary and source code
 form.
*/

#include <math.h>

#include <iostream>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

bool bad_alignment(Aln *aln) {
  if (aln->rstart > KLIGN_UNALIGNED_THRES) return true;
  if (aln->rlen - aln->rstop > KLIGN_UNALIGNED_THRES) return true;
  // with match and mismatch scores of 1 this means at most two errors
  if (aln->score1 < aln->rlen - 2) return true;
  return false;
}

pair<int, int> calculate_insert_size(Alns &alns, int expected_ins_avg, int expected_ins_stddev, int max_expected_ins_size,
                                     const string &dump_large_alns_fname = "") {
  BarrierTimer timer(__FILEFUNC__);
  auto unmerged_rlen = alns.calculate_unmerged_rlen();
  ProgressBar progbar(alns.size(), "Processing alignments to compute insert size");
  Aln *prev_aln = nullptr;
  int64_t prev_cid = -1;
  int64_t num_valid_pairs = 0;
  int64_t num_overlap_rejected = 0;
  int64_t sum_insert_size = 0;
  int64_t num_large = 0;
  int64_t num_small = 0;
  int64_t num_repetitive_conflicts = 0;
  int min_insert_size = 1000000;
  int max_insert_size = 0;
  vector<string> large_alns;
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
        if (prev_pair_num != '1' || pair_num != '2') {
#ifdef DEBUG
          // i.e. when one read mate has two mappings to the same contig (i.e. repetitive region)
          // OR this can happen when the mapping is wildly wrong.  Keep a count an report it.
          LOG("pair nums wrong: ", prev_pair_num, " ", pair_num, ", aln:", aln.to_string(), "\n");
#endif
          prev_aln = nullptr;
          num_repetitive_conflicts++;
          continue;
        }
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
            // WARN("large insert size: ", insert_size, " prev: ", prev_aln->clen, " ",
            //     prev_aln->orient, " ", prev_aln->cstart, " ", prev_aln->cstop, " ", prev_cstart, " ", prev_cstop,
            //     " current ", aln.clen, " ", aln.orient, " ", aln.cstart, " ", aln.cstop);
            // the metaquast extensive misassembly minimum size is 1000
            if (!dump_large_alns_fname.empty())
              large_alns.push_back(to_string(prev_aln->cid) + " " + to_string(prev_aln->clen) + " " +
                                   to_string(min(prev_cstart, aln.cstart)) + " " + to_string(max(prev_cstop, aln.cstop)));
            num_large++;
          } else if (insert_size < 10) {
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

  if (!dump_large_alns_fname.empty()) {
    upcxx_utils::dist_ofstream outf(dump_large_alns_fname);
    if (!upcxx::rank_me()) {
      outf << "#cid clen aln_start aln_stop\n";
    }
    for (auto large_alns_ctg : large_alns) {
      outf << large_alns_ctg + "\n";
    }
    outf.close();  // syncs and write stats
  }

  auto all_num_repetitive_conflicts = reduce_one(num_repetitive_conflicts, op_fast_add, 0).wait();
  auto all_num_overlap_rejected = reduce_one(num_overlap_rejected, op_fast_add, 0).wait();
  auto all_num_valid_pairs = reduce_all(num_valid_pairs, op_fast_add).wait();
  auto all_num_alns = reduce_one(alns.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(all_num_valid_pairs, all_num_alns), " pairs with valid alignments to the same contig\n");
  SLOG_VERBOSE("Rejected ", perc_str(all_num_overlap_rejected, all_num_alns), " possible candidates with bad overlaps\n");
  SLOG_VERBOSE("Repetitive Conflicts: ", perc_str(all_num_repetitive_conflicts, all_num_alns),
               " pairs with multiple mappings of the same mate to the same contig\n");
  auto all_num_large = reduce_one(num_large, op_fast_add, 0).wait();
  auto all_num_small = reduce_one(num_small, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(all_num_large, all_num_alns), " large (> ", max_expected_ins_size, ") and ",
               perc_str(all_num_small, all_num_alns), " small outliers\n");

  auto all_sum_insert_size = reduce_all(sum_insert_size, op_fast_add).wait();
  if (!all_num_valid_pairs) {
    if (expected_ins_avg) {
      WARN("Could not find any suitable alignments for calculating the insert size. Using the parameters ", expected_ins_avg, " ",
           expected_ins_stddev);
      return {expected_ins_avg, expected_ins_stddev};
    } else {
      SWARN("Could not find any suitable alignments for calculating the insert size and no parameters are set.");
      return {0, 0};
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
  SLOG_VERBOSE("Calculated insert size: average ", insert_avg, " stddev ", insert_stddev, " min ", all_min_insert_size, " max ",
               all_max_insert_size, "\n");
  if (expected_ins_avg) {
    if (abs(insert_avg - expected_ins_avg) > 100)
      SWARN("Large difference in calculated (", insert_avg, ") vs expected (", expected_ins_avg, ") insert average sizes");
    if (abs(insert_stddev - expected_ins_stddev) > 100)
      SWARN("Large difference in calculated (", insert_stddev, ") vs expected (", expected_ins_stddev,
            ") insert standard deviation");
    SLOG_VERBOSE("Using specified ", expected_ins_avg, " avg insert size and ", expected_ins_stddev, " stddev\n");
    return {expected_ins_avg, expected_ins_stddev};
  }
  return {insert_avg, insert_stddev};
}
