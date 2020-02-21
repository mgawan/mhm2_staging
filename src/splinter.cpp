#include <iostream>
#include <fstream>
#include <regex>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "progressbar.hpp"
#include "zstr.hpp"
#include "ctg_graph.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "fastq.hpp"

using namespace std;
using namespace upcxx;


static CtgGraph *_graph = nullptr;



struct AlnStats {
  int64_t nalns, unaligned, short_alns, containments, circular, bad_overlaps;

  void print() {
    int64_t tot_nalns = reduce_one(nalns, op_fast_add, 0).wait();
    SLOG_VERBOSE(setprecision(2), fixed,
                 "Processed ", tot_nalns, " alignments:\n",
                 "    unaligned:          ", perc_str(reduce_one(unaligned, op_fast_add, 0).wait(), tot_nalns), "\n",
                 "    containments:       ", perc_str(reduce_one(containments, op_fast_add, 0).wait(), tot_nalns), "\n",
                 "    short alns:         ", perc_str(reduce_one(short_alns, op_fast_add, 0).wait(), tot_nalns), "\n",
                 "    circular:           ", perc_str(reduce_one(circular, op_fast_add, 0).wait(), tot_nalns), "\n", 
                 "    bad overlaps:       ", perc_str(reduce_one(bad_overlaps, op_fast_add, 0).wait(), tot_nalns), "\n");
  }
};


// gets all the alns for a single read
static void get_alns_for_read(Alns &alns, int64_t &i, vector<Aln> &alns_for_read, int64_t *nalns, int64_t *unaligned) {
  string start_read_id = "";
  for (; i < alns.size(); i++) {
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
    if (unaligned_left > ALN_WIGGLE || unaligned_right > ALN_WIGGLE) {
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
    return { .start = aln->rstart - aln->cstart, .stop = aln->rstop + (aln->clen - aln->cstop) };
  };

  auto is_contained = [](const AlnCoords &ctg1, const AlnCoords &ctg2) -> bool {
    if (ctg1.start >= ctg2.start && ctg1.stop <= ctg2.stop) return true;
    else return false;
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
  CidPair cids = { .cid1 = aln1->cid, .cid2 = aln2->cid };
  if (aln1->cid < aln2->cid) {
    swap(end1, end2);
    swap(cids.cid1, cids.cid2);
    orient1 = aln2->orient;
  }
  Edge edge = { .cids = cids, .end1 = end1, .end2 = end2, .gap = gap, .support = 1, .aln_len = min_aln_len,
                .aln_score = min_aln_score, .edge_type = EdgeType::SPLINT, .seq = "",
                .mismatch_error = false, .conflict_error = false, .excess_error = false, .short_aln = false, .gap_reads = {}};
  if (edge.gap > 0) {
    edge.gap_reads = vector<GapRead>{GapRead(aln1->read_id, gap_start, orient1, cids.cid1)};
    _graph->add_pos_gap_read(aln1->read_id);
  }
  _graph->add_or_update_edge(edge);
  return true;
}


void get_splints_from_alns(Alns &alns, CtgGraph *graph) {
  Timer timer(__FILEFUNC__);
  _graph = graph;
  AlnStats stats = {0};
  IntermittentTimer t_get_alns(__FILENAME__ + string(":") + "get alns splints");
  ProgressBar progbar(alns.size(), "Adding edges to graph from splints");
  int64_t aln_i = 0;
  int64_t num_splints = 0;
  while (aln_i < alns.size()) {
    vector<Aln> alns_for_read;
    t_get_alns.start();
    get_alns_for_read(alns, aln_i, alns_for_read, &stats.nalns, &stats.unaligned);
    t_get_alns.stop();
    progbar.update(aln_i);
    for (int i = 0; i < alns_for_read.size(); i++) {
      auto aln = &alns_for_read[i];
      for (int j = i + 1; j < alns_for_read.size(); j++) {
        progress();
        auto other_aln = &alns_for_read[j];
        if (other_aln->read_id != aln->read_id) DIE("Mismatched read ids: ", other_aln->read_id, " != ", aln->read_id, "\n");
        if (add_splint(other_aln, aln, stats)) num_splints++;
      }
    }
  }
  progbar.done();
  barrier();
  t_get_alns.done_barrier();
  stats.print();
  SLOG_VERBOSE("Found ", reduce_one(num_splints, op_fast_add, 0).wait(), " splints\n");
}


