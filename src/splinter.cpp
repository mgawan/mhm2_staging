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

#include <fstream>
#include <iostream>
#include <regex>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "contigs.hpp"
#include "ctg_graph.hpp"
#include "fastq.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "utils.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

static CtgGraph *_graph = nullptr;

struct AlnStats {
  int64_t nalns, unaligned, short_alns, containments, circular, bad_overlaps;

  void print() {
    int64_t tot_nalns = reduce_one(nalns, op_fast_add, 0).wait();
    SLOG_VERBOSE(setprecision(2), fixed, "Processed ", tot_nalns, " alignments:\n");
    SLOG_VERBOSE(setprecision(2), fixed,
                 "    unaligned:          ", perc_str(reduce_one(unaligned, op_fast_add, 0).wait(), tot_nalns), "\n");
    SLOG_VERBOSE(setprecision(2), fixed,
                 "    containments:       ", perc_str(reduce_one(containments, op_fast_add, 0).wait(), tot_nalns), "\n");
    SLOG_VERBOSE(setprecision(2), fixed,
                 "    short alns:         ", perc_str(reduce_one(short_alns, op_fast_add, 0).wait(), tot_nalns), "\n");
    SLOG_VERBOSE(setprecision(2), fixed,
                 "    circular:           ", perc_str(reduce_one(circular, op_fast_add, 0).wait(), tot_nalns), "\n");
    SLOG_VERBOSE(setprecision(2), fixed,
                 "    bad overlaps:       ", perc_str(reduce_one(bad_overlaps, op_fast_add, 0).wait(), tot_nalns), "\n");
  }
};

// gets all the alns for a single read
static void get_alns_for_read(Alns &alns, int64_t &i, vector<Aln> &alns_for_read, int64_t *nalns, int64_t *unaligned) {
  string start_read_id = "";
  for (; i < (int64_t)alns.size(); i++) {
    Aln aln = alns.get_aln(i);
    // alns for a new read
    if (start_read_id != "" && aln.read_id != start_read_id) return;
    (*nalns)++;
    // convert to coords for use here
    if (aln.orient == '-') {
      int tmp = aln.cstart;
      aln.cstart = aln.clen - aln.cstop;
      aln.cstop = aln.clen - tmp;
    }
    int unaligned_left = min(aln.rstart, aln.cstart);
    int unaligned_right = min(aln.rlen - aln.rstop, aln.clen - aln.cstop);
    if (unaligned_left > KLIGN_UNALIGNED_THRES || unaligned_right > KLIGN_UNALIGNED_THRES) {
      (*unaligned)++;
      //      DBG("unaligned ", aln.read_id, " ", aln.rstart, " ", aln.rstop, " ", aln.rlen, " ",
      //          aln.cid, " ", aln.cstart, " ", aln.cstop, " ", aln.clen, " ", aln.orient, " ", aln.score1, " ", aln.score2, "\n");
    } else {
      alns_for_read.push_back(aln);
    }
    start_read_id = aln.read_id;
    progress();
  }
}

static bool add_splint(const Aln *aln1, const Aln *aln2, AlnStats &stats) {
  struct AlnCoords {
    int start, stop;
  };

  auto get_aln_coords = [](const Aln *aln) -> AlnCoords {
    // get contig start and end positions in read coordinate set
    return {.start = aln->rstart - aln->cstart, .stop = aln->rstop + (aln->clen - aln->cstop)};
  };

  auto is_contained = [](const AlnCoords &ctg1, const AlnCoords &ctg2) -> bool {
    if (ctg1.start >= ctg2.start && ctg1.stop <= ctg2.stop)
      return true;
    else
      return false;
  };

  AlnCoords ctg1 = get_aln_coords(aln1);
  AlnCoords ctg2 = get_aln_coords(aln2);

  if (is_contained(ctg1, ctg2) || is_contained(ctg2, ctg1)) {
    stats.containments++;
    return false;
  }
  if (aln1->cid == aln2->cid) {
    stats.circular++;
    return false;
  }
  int end1, end2;
  int gap;
  int gap_start;
  if (ctg1.start <= ctg2.start) {
    end1 = aln1->orient == '+' ? 3 : 5;
    end2 = aln2->orient == '+' ? 5 : 3;
    gap = ctg2.start - ctg1.stop;
    gap_start = ctg1.stop;
  } else {
    end1 = aln1->orient == '+' ? 5 : 3;
    end2 = aln2->orient == '+' ? 3 : 5;
    gap = ctg1.start - ctg2.stop;
    gap_start = ctg2.stop;
  }
  int min_aln_len = min(aln1->rstop - aln1->rstart, aln2->rstop - aln2->rstart);
  int min_aln_score = min(aln1->score1, aln2->score1);
  if (gap < -min_aln_len) {
    stats.short_alns++;
    return false;
  }
  if (gap < -aln1->rlen) DIE("Gap is too small: ", gap, ", read length ", aln1->rlen, "\n");
  // check for bad overlaps
  if (gap < 0 && (aln1->clen < -gap || aln2->clen < -gap)) {
    stats.bad_overlaps++;
    return false;
  }
  char orient1 = aln1->orient;
  CidPair cids = {.cid1 = aln1->cid, .cid2 = aln2->cid};
  if (aln1->cid < aln2->cid) {
    swap(end1, end2);
    swap(cids.cid1, cids.cid2);
    orient1 = aln2->orient;
  }
  Edge edge = {.cids = cids,
               .end1 = end1,
               .end2 = end2,
               .gap = gap,
               .support = 1,
               .aln_len = min_aln_len,
               .aln_score = min_aln_score,
               .edge_type = EdgeType::SPLINT,
               .seq = "",
               .mismatch_error = false,
               .conflict_error = false,
               .excess_error = false,
               .short_aln = false,
               .gap_reads = {}};
  if (edge.gap > 0) {
    edge.gap_reads = vector<GapRead>{GapRead(aln1->read_id, gap_start, orient1, cids.cid1)};
    _graph->add_pos_gap_read(aln1->read_id);
  }
  _graph->add_or_update_edge(edge);
  return true;
}

void get_splints_from_alns(Alns &alns, CtgGraph *graph) {
  BarrierTimer timer(__FILEFUNC__);
  _graph = graph;
  AlnStats stats = {0};
  IntermittentTimer t_get_alns(__FILENAME__ + string(":") + "get alns splints");
  ProgressBar progbar(alns.size(), "Adding edges to graph from splints");
  int64_t aln_i = 0;
  int64_t num_splints = 0;
  while (aln_i < (int64_t)alns.size()) {
    vector<Aln> alns_for_read;
    t_get_alns.start();
    get_alns_for_read(alns, aln_i, alns_for_read, &stats.nalns, &stats.unaligned);
    t_get_alns.stop();
    progbar.update(aln_i);
    for (int i = 0; i < (int)alns_for_read.size(); i++) {
      auto aln = &alns_for_read[i];
      for (int j = i + 1; j < (int)alns_for_read.size(); j++) {
        progress();
        auto other_aln = &alns_for_read[j];
        if (other_aln->read_id != aln->read_id) DIE("Mismatched read ids: ", other_aln->read_id, " != ", aln->read_id, "\n");
        if (add_splint(other_aln, aln, stats)) num_splints++;
      }
    }
  }
  progbar.done();
  barrier();
  t_get_alns.done_all();
  stats.print();
  SLOG_VERBOSE("Found ", reduce_one(num_splints, op_fast_add, 0).wait(), " splints\n");
}
