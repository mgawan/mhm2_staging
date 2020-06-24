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

#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <string>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/flat_aggr_store.hpp"
#include "upcxx_utils/timers.hpp"

#include "utils.hpp"
#include "contigs.hpp"
#include "alignments.hpp"

using namespace std;

struct CtgBaseDepths {
  int64_t cid;
  vector<int> base_counts;

  UPCXX_SERIALIZED_FIELDS(cid, base_counts);
};

class CtgsDepths {
 private:
  using ctgs_depths_map_t = upcxx::dist_object<HASH_TABLE<int64_t, CtgBaseDepths>>;
  ctgs_depths_map_t ctgs_depths;
  HASH_TABLE<int64_t, CtgBaseDepths>::iterator ctgs_depths_iter;

  size_t get_target_rank(int64_t cid) {
    return std::hash<int64_t>{}(cid) % upcxx::rank_n();
  }

 public:
  CtgsDepths() : ctgs_depths({}) {}

  int64_t get_num_ctgs() {
    return upcxx::reduce_one(ctgs_depths->size(), upcxx::op_fast_add, 0).wait();
  }

  void add_ctg(CtgBaseDepths &ctg) {
    upcxx::rpc(
        get_target_rank(ctg.cid),
        [](ctgs_depths_map_t &ctgs_depths, CtgBaseDepths ctg) {
          ctgs_depths->insert({ctg.cid, ctg});
        },
        ctgs_depths, ctg)
        .wait();
  }

  void update_ctg_aln_depth(int64_t cid, int aln_start, int aln_stop) {
    upcxx::rpc(
        get_target_rank(cid),
        [](ctgs_depths_map_t &ctgs_depths, int64_t cid, int aln_start, int aln_stop) {
          assert(aln_start < aln_stop);
          const auto it = ctgs_depths->find(cid);
          if (it == ctgs_depths->end()) DIE("could not fetch vertex ", cid, "\n");
          auto ctg = &it->second;
          for (int i = aln_start; i < aln_stop; i++) {
            ctg->base_counts[i]++;
          }
        },
        ctgs_depths, cid, aln_start, aln_stop)
        .wait();
  }

  CtgBaseDepths *get_first_local_ctg() {
    ctgs_depths_iter = ctgs_depths->begin();
    if (ctgs_depths_iter == ctgs_depths->end()) return nullptr;
    auto ctg = &ctgs_depths_iter->second;
    ctgs_depths_iter++;
    return ctg;
  }

  CtgBaseDepths *get_next_local_ctg() {
    if (ctgs_depths_iter == ctgs_depths->end()) return nullptr;
    auto ctg = &ctgs_depths_iter->second;
    ctgs_depths_iter++;
    return ctg;
  }

};

void compute_aln_depths(const string &fname, Contigs &ctgs, Alns &alns, int kmer_len, int min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  CtgsDepths ctgs_depths;
  SLOG_VERBOSE("Loading contigs\n");
  for (auto &ctg : ctgs) {
    int clen = ctg.seq.length();
    if (clen < min_ctg_len) continue;
    CtgBaseDepths ctg_base_depths = {.cid = ctg.id, .base_counts = vector<int>(clen, 0)};
    ctgs_depths.add_ctg(ctg_base_depths);
  }
  barrier();
  auto num_ctgs = ctgs_depths.get_num_ctgs();
  SLOG_VERBOSE("Computing aln depths for ", num_ctgs, " ctgs\n");
  ProgressBar progbar(alns.size(), "Processing alignments");
  for (auto &aln : alns) {
    progbar.update();
    // convert to coords for use here
    auto cstart = aln.cstart;
    auto cstop = aln.cstop;
    if (aln.orient == '-') {
      int tmp = cstart;
      cstart = aln.clen - cstop;
      cstop = aln.clen - tmp;
    }
    int unaligned_left = min(aln.rstart, cstart);
    int unaligned_right = min(aln.rlen - aln.rstop, aln.clen - cstop);
    if (unaligned_left <= KLIGN_UNALIGNED_THRES && unaligned_right <= KLIGN_UNALIGNED_THRES)
      ctgs_depths.update_ctg_aln_depth(aln.cid, aln.cstart, aln.cstop);
  }
  progbar.done();
  barrier();
  // get string to dump
  string out_str = "";
  for (auto ctg = ctgs_depths.get_first_local_ctg(); ctg != nullptr; ctg = ctgs_depths.get_next_local_ctg()) {
    double avg_depth = 0;
    for (auto base_count : ctg->base_counts) {
      avg_depth += base_count;
    }
    avg_depth /= ctg->base_counts.size();
    double sum_sqs = 0;
    for (auto base_count : ctg->base_counts) {
      sum_sqs += pow((double)base_count - avg_depth, 2.0);
    }
    double stddev_depth = sqrt(sum_sqs / ctg->base_counts.size());
    if (avg_depth < 2) avg_depth = 2;
    out_str += "Contig" + to_string(ctg->cid) + "\t" + to_string(avg_depth) + "\t" + to_string(stddev_depth) + "\n";
  }
  dump_single_file(fname, out_str);
  barrier();
}
