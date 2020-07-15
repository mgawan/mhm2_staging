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
  int edge_base_len;
  HASH_TABLE<int64_t, CtgBaseDepths>::iterator ctgs_depths_iter;

  size_t get_target_rank(int64_t cid) {
    return std::hash<int64_t>{}(cid) % upcxx::rank_n();
  }

 public:
  CtgsDepths(int edge_base_len) : ctgs_depths({}), edge_base_len(edge_base_len) {}

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

  void update_ctg_aln_depth(int64_t cid, int aln_start, int aln_stop, int aln_merge_start, int aln_merge_stop) {
    if (aln_start >= aln_stop) return;
    upcxx::rpc(
        get_target_rank(cid),
        [](ctgs_depths_map_t &ctgs_depths, int64_t cid, int aln_start, int aln_stop, int aln_merge_start, int aln_merge_stop) {
          const auto it = ctgs_depths->find(cid);
          if (it == ctgs_depths->end()) DIE("could not fetch vertex ", cid, "\n");
          auto ctg = &it->second;
          for (int i = aln_start; i < aln_stop; i++) {
            ctg->base_counts[i]++;
          }
          // this is the merged region in a merged read - should be counted double
          if (aln_merge_start != -1 && aln_merge_stop != -1) {
            for (int i = aln_merge_start; i < aln_merge_stop; i++) {
              ctg->base_counts[i]++;
            }
          }
        },
        ctgs_depths, cid, aln_start, aln_stop, aln_merge_start, aln_merge_stop)
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

  std::pair<double, double> get_depth(int64_t cid) {
    return upcxx::rpc(
               get_target_rank(cid),
               [](ctgs_depths_map_t &ctgs_depths, int64_t cid, int edge_base_len) -> pair<double, double> {
                 const auto it = ctgs_depths->find(cid);
                 if (it == ctgs_depths->end()) DIE("could not fetch vertex ", cid, "\n");
                 auto ctg_base_depths = &it->second;
                 double avg_depth = 0;
                 for (int i = edge_base_len; i < (int) ctg_base_depths->base_counts.size() - edge_base_len; i++) {
                   avg_depth += ctg_base_depths->base_counts[i];
                 }
                 size_t clen = ctg_base_depths->base_counts.size() - 2 * edge_base_len;
                 avg_depth /= clen;
                 double sum_sqs = 0;
                 for (int i = edge_base_len; i < (int) ctg_base_depths->base_counts.size() - edge_base_len; i++) {
                   sum_sqs += pow((double)ctg_base_depths->base_counts[i] - avg_depth, 2.0);
                 }
                 double var_depth = sum_sqs / clen;
                 if (avg_depth < 2) avg_depth = 2;
                 return {avg_depth, var_depth};
               },
               ctgs_depths, cid, edge_base_len)
        .wait();
  }

};

void compute_aln_depths(const string &fname, Contigs &ctgs, Alns &alns, int kmer_len, int min_ctg_len, bool use_kmer_depths) {
  BarrierTimer timer(__FILEFUNC__);
  int edge_base_len = (min_ctg_len >= 75 ? 75 : 0);
  CtgsDepths ctgs_depths(edge_base_len);
  SLOG_VERBOSE("Processing contigs, using an edge base length of ", edge_base_len, " and a min ctg len of ", min_ctg_len, "\n");
  for (auto &ctg : ctgs) {
    int clen = ctg.seq.length();
    if (clen < min_ctg_len) continue;
    CtgBaseDepths ctg_base_depths = {.cid = ctg.id, .base_counts = vector<int>(clen, 0)};
    ctgs_depths.add_ctg(ctg_base_depths);
    upcxx::progress();
  }
  barrier();
  auto unmerged_rlen = alns.calculate_unmerged_rlen();
  int64_t num_bad_overlaps = 0;
  int64_t num_bad_alns = 0;
  auto num_ctgs = ctgs_depths.get_num_ctgs();
  SLOG_VERBOSE("Computing aln depths for ", num_ctgs, " ctgs\n");
  ProgressBar progbar(alns.size(), "Processing alignments");
  for (auto &aln : alns) {
    progbar.update();
    // require at least this much overlap with the read
    // what this does is drop alns that hang too much over the ends of contigs
    // this gives abundances more in line with what we see in MetaBAT, although that uses 97% identity as the cutoff and we're using
    // 85% here (our alns differ somewhat because of different seed lengths, etc) .
    // In practice, when using aln depths for scaffolding, this tends to reduce msa without any benefits so we only use it in the
    // final round, i.e. if min_ctg_len > 0
    /*
    if (min_ctg_len && aln.identity < 85) {
      num_bad_alns++;
      continue;
    }
    */
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
    if (unaligned_left <= KLIGN_UNALIGNED_THRES && unaligned_right <= KLIGN_UNALIGNED_THRES) {
      // set to -1 if this read is not merged
      int aln_cstart_merge = -1, aln_cstop_merge = -1;
      // FIXME: need to somehow communicate to the update func the range of double counting for a merged read.
      // This is the area > length of read pair that is in the middle of the read
      if (aln.rlen > unmerged_rlen) {
        // merged read
        int merge_offset = (aln.rlen - unmerged_rlen) / 2;
        aln_cstart_merge = (merge_offset > aln.rstart ? merge_offset - aln.rstart : 0) + aln.cstart;
        int stop_merge = aln.rlen - merge_offset;
        aln_cstop_merge = aln.cstop - (stop_merge < aln.rstop ? aln.rstop - stop_merge : 0);
        // the aln may not include the merged region
        if (aln_cstart_merge >= aln_cstop_merge) aln_cstart_merge = -1;
      }
      // as per MetaBAT analysis, ignore the 75 bases at either end because they are likely to be in error
      ctgs_depths.update_ctg_aln_depth(aln.cid, std::max(aln.cstart, edge_base_len), std::min(aln.cstop, aln.clen - edge_base_len),
                                       aln_cstart_merge, aln_cstop_merge);
    } else {
      num_bad_overlaps++;
    }
    upcxx::progress();
  }
  progbar.done();
  barrier();
  auto all_num_alns = reduce_one(alns.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Dropped ", perc_str(reduce_one(num_bad_overlaps, op_fast_add, 0).wait(), all_num_alns), " bad overlaps and ",
               perc_str(reduce_one(num_bad_alns, op_fast_add, 0).wait(), all_num_alns), " low quality alns\n");
  // get string to dump
  string out_str = "";
  if (!upcxx::rank_me()) out_str = "contigName\tcontigLen\ttotalAvgDepth\tavg_depth\tvar_depth\n";
  // FIXME: the depths need to be in the same order as the contigs in the final_assembly.fasta file. This is an inefficient
  // way of ensuring that
  for (auto &ctg : ctgs) {
    if ((int) ctg.seq.length() < min_ctg_len) continue;
    auto [avg_depth, var_depth] = ctgs_depths.get_depth(ctg.id);
    ostringstream oss;
    oss << "Contig" << ctg.id << "\t" << ctg.seq.length() << "\t" << avg_depth << "\t" << avg_depth << "\t" << var_depth << "\n";
    out_str += oss.str();
    // it seems that using aln depths improves the ctgy at the cost of an increase in msa
    if (!use_kmer_depths) ctg.depth = avg_depth;
    upcxx::progress();
  }
  dump_single_file(fname, out_str);
}
