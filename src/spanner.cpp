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
#include "zstr.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

//#define DUMP_LINKS

static CtgGraph *_graph = nullptr;

// functions for evaluating the gap correction analytically assuming Gaussian insert size distribution
// this is taken from meraculous bmaToLinks.pl

double erf(int x) {
  double absX = (x < 0) ? -x : x;
  double t = 1 / (1 + 0.5 * absX);
  double t2 = t * t;
  double t3 = t * t2;
  double t4 = t * t3;
  double t5 = t * t4;
  double t6 = t * t5;
  double t7 = t * t6;
  double t8 = t * t7;
  double t9 = t * t8;
  double a0 = -1.26551223;
  double a1 = 1.00002368;
  double a2 = 0.37409196;
  double a3 = 0.09678418;
  double a4 = -0.18628806;
  double a5 = 0.27886807;
  double a6 = -1.13520398;
  double a7 = 1.48851587;
  double a8 = -0.82215223;
  double a9 = 0.17087277;
  double tau = t * exp(-x * x + a0 + a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5 + a6 * t6 + a7 * t7 + a8 * t8 + a9 * t9);
  if (x < 0)
    return (tau - 1);
  else
    return (1 - tau);
}

double G0(double a, double b, double mu, double sigma) {
  const double sqrtTwo = sqrt(2);
  const double pi = 3.14159265359;
  const double sqrtPi = sqrt(pi);

  double rt2sig = sqrtTwo * sigma;
  double erfa = erf((a - mu) / rt2sig);
  double erfb = erf((b - mu) / rt2sig);
  return (sqrtPi / sqrtTwo) * sigma * (erfb - erfa);
}

double G1(double a, double b, double mu, double sigma) {
  double za = (a - mu) / sigma;
  double zb = (b - mu) / sigma;
  double expa = exp(-0.5 * za * za);
  double expb = exp(-0.5 * zb * zb);
  double g0 = G0(a, b, mu, sigma);
  return (sigma * sigma * (expa - expb) + mu * g0);
}

double G2(double a, double b, double mu, double sigma) {
  double za = (a - mu) / sigma;
  double zb = (b - mu) / sigma;
  double expa = exp(-0.5 * za * za);
  double expb = exp(-0.5 * zb * zb);
  double g0 = G0(a, b, mu, sigma);
  double g1 = G1(a, b, mu, sigma);
  double sigma2 = sigma * sigma;
  return (sigma2 * g0 + mu * g1 + sigma2 * (a * expa - b * expb));
}

double mean_spanning_clone(double g, double k, double l, double c1, double c2, double mu, double sigma) {
  double x1 = g + 2 * k - 1;
  double x2 = g + c1 + l;
  double alpha = x2 - x1;
  double x3 = g + c2 + l;
  double x4 = x3 + alpha;
  double num = 0;
  double den = 0;
  double N1 = G2(x1, x2, mu, sigma) - x1 * G1(x1, x2, mu, sigma);
  double N2 = (x2 - x1) * G1(x2, x3, mu, sigma);
  double N3 = x4 * G1(x3, x4, mu, sigma) - G2(x3, x4, mu, sigma);
  double D1 = G1(x1, x2, mu, sigma) - x1 * G0(x1, x2, mu, sigma);
  double D2 = (x2 - x1) * G0(x2, x3, mu, sigma);
  double D3 = x4 * G0(x3, x4, mu, sigma) - G1(x3, x4, mu, sigma);
  num = N1 + N2 + N3;
  den = D1 + D2 + D3;
  if (den) {
    return num / den;
  } else {
    // WARN("mean_spanning_clone failed for (g,k,l,c1,c2,mu,sigma)");
    return 0;
  }
}

double estimate_gap_size(double meanAnchor, double k, double l, double c1, double c2, double mu, double sigma) {
  double gMax = mu + 3 * sigma - 2 * k;
  double gMin = -(k - 2);
  double gMid = mu - meanAnchor;
  // Negative gap size padding disabled for metagenomes
  // if (gMid < gMin) gMid = gMin + 1;
  if (gMid > gMax) gMid = gMax - 1;
  double aMax = mean_spanning_clone(gMax, k, l, c1, c2, mu, sigma) - gMax;
  double aMin = mean_spanning_clone(gMin, k, l, c1, c2, mu, sigma) - gMin;
  double aMid = mean_spanning_clone(gMid, k, l, c1, c2, mu, sigma) - gMid;
  double deltaG = gMax - gMin;
  double iterations = 0;
  while (deltaG > 10) {
    iterations++;
    if (meanAnchor > aMid) {
      gMax = gMid;
      aMax = aMid;
      gMid = (gMid + gMin) / 2;
      aMid = mean_spanning_clone(gMid, k, l, c1, c2, mu, sigma) - gMid;
    } else if (meanAnchor < aMid) {
      gMin = gMid;
      aMin = aMid;
      gMid = (gMid + gMax) / 2;
      aMid = mean_spanning_clone(gMid, k, l, c1, c2, mu, sigma) - gMid;
    } else {
      break;
    }
    deltaG = gMax - gMin;
  }
  return gMid;
}

static bool get_best_span_aln(int insert_avg, int insert_stddev, vector<Aln> &alns, Aln &best_aln, string &read_status,
                              string &type_status, int64_t *reject_5_trunc, int64_t *reject_3_trunc, int64_t *reject_uninf) {
  if (alns.size() == 0) return false;
  int min_rstart = -1;
  read_status = "";
  // sort in order of highest to lowest aln scores
  sort(alns.begin(), alns.end(), [](auto &aln1, auto &aln2) { return aln1.score1 > aln2.score1; });
  for (const auto &aln : alns) {
    // Assess alignment for completeness (do this before scaffold coordinate conversion!)
    // Important: Truncations are performed before reverse complementation
    // and apply to the end of the actual read
    int rstart = aln.rstart + 1;
    int cstart = aln.cstart + 1;
    int unaligned_start = rstart - 1;
    int projected_start = (aln.orient == '+' ? cstart - unaligned_start : aln.cstop + unaligned_start);
    string start_status = "";
    int projected_off = 0;
    if (projected_start < 1)
      projected_off = 1 - projected_start;
    else if (projected_start > aln.clen)
      projected_off = projected_start - aln.clen;
    int missing_start_bases = unaligned_start - projected_off;
    if (unaligned_start == 0) {
      start_status = "FULL";
    } else if (projected_off > 0 && missing_start_bases < KLIGN_UNALIGNED_THRES) {
      start_status = "GAP";
    } else if (unaligned_start < KLIGN_UNALIGNED_THRES) {
      start_status = "INC";
    } else {
      (*reject_5_trunc)++;
      continue;
    }

    int unaligned_end = aln.rlen - aln.rstop;
    int projected_end = (aln.orient == '+' ? aln.cstop + unaligned_end : cstart - unaligned_end);
    string end_status = "";
    projected_off = 0;
    if (projected_end < 1)
      projected_off = 1 - projected_end;
    else if (projected_end > aln.clen)
      projected_off = projected_end - aln.clen;
    int missing_end_bases = unaligned_end - projected_off;
    if (unaligned_end == 0) {
      end_status = "FULL";
    } else if (projected_off > 0 && missing_end_bases < KLIGN_UNALIGNED_THRES) {
      end_status = "GAP";
    } else if (unaligned_end < KLIGN_UNALIGNED_THRES) {
      end_status = "INC";
    } else {
      (*reject_3_trunc)++;
      continue;
    }

    int endDistance = insert_avg + 3 * insert_stddev;
    // Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle)
    if ((aln.orient == '+' && projected_start > aln.clen - endDistance) || (aln.orient == '-' && projected_start < endDistance)) {
      if (min_rstart == -1 || rstart < min_rstart) {
        read_status = start_status + "." + end_status + ".OUT";
        if (start_status == "GAP")
          type_status = "OUTGAP";
        else if (end_status == "GAP")
          type_status = "INTGAP";
        else
          type_status = "ANCHOR";
        best_aln = aln;
        min_rstart = rstart;
      } else {
        (*reject_uninf)++;
      }
    } else {
      (*reject_uninf)++;
    }
  }
  if (min_rstart == -1) return false;
  return true;
}

// gets all the alns for a single read
static void get_all_alns_for_read(Alns &alns, int64_t &i, vector<Aln> &alns_for_read) {
  string start_read_id = "";
  for (; i < (int64_t)alns.size(); i++) {
    Aln aln = alns.get_aln(i);
    // alns for a new read
    if (start_read_id != "" && aln.read_id != start_read_id) return;
    alns_for_read.push_back(aln);
    start_read_id = aln.read_id;
  }
}

static void add_span_pos_gap_read(Edge *edge, Aln &aln) {
  int gap_start = 0;
  if (aln.cid == edge->cids.cid1) {
    if ((edge->end1 == 3 && aln.orient == '+') || (edge->end1 == 5 && aln.orient == '-'))
      gap_start = aln.rstop;
    else
      gap_start = aln.rlen - aln.rstart;
    if (gap_start == aln.rlen) return;
  } else if (aln.cid == edge->cids.cid2) {
    // head
    if ((edge->end2 == 5 && aln.orient == '+') || (edge->end2 == 3 && aln.orient == '-'))
      gap_start = aln.rstart;
    else
      gap_start = aln.rlen - aln.rstop;
    if (gap_start == 0) return;
  } else {
    DIE("cid doesn't match in pos edge\n");
  }
  edge->gap_reads.emplace_back(aln.read_id, gap_start, aln.orient, aln.cid);
  _graph->add_pos_gap_read(aln.read_id);
}

enum class ProcessPairResult { FAIL_SMALL, FAIL_SELF_LINK, FAIL_EDIST, FAIL_MIN_GAP, SUCCESS };

ProcessPairResult process_pair(int insert_avg, int insert_stddev, Aln &aln1, Aln &aln2, const string &type_status1,
                               const string &type_status2, const string &read_status1, const string &read_status2) {
  auto get_dist = [=](int &d, Aln &aln) -> bool {
    aln.rstart++;
    aln.cstart++;
    // Omit alignments to very short (potentially "popped-out") scaffolds
    if (aln.clen < insert_avg / 2) return false;
    d = (aln.orient == '+' ? (aln.clen - aln.cstart + 1) + (aln.rstart - 1) : aln.cstop + (aln.rstart - 1));
    return true;
  };

  DBG_VERBOSE("process_pair: ", aln1.to_string(), " vs ", aln2.to_string(), "\n");
  assert(aln1.read_id.size() == aln2.read_id.size());
  assert(aln1.read_id.substr(0, aln1.read_id.size() - 1).compare(aln2.read_id.substr(0, aln2.read_id.size() - 1)) == 0);
  assert(aln1.read_id[aln1.read_id.size() - 1] != aln2.read_id[aln2.read_id.size() - 1]);
  assert(aln1.read_id != aln2.read_id);

  // check for inappropriate self-linkage
  if (aln1.cid == aln2.cid) return ProcessPairResult::FAIL_SELF_LINK;
  int end_distance = insert_avg + 3 * insert_stddev;
  int d1, d2;
  if (!get_dist(d1, aln1)) return ProcessPairResult::FAIL_SMALL;
  if (!get_dist(d2, aln2)) return ProcessPairResult::FAIL_SMALL;
  if (d1 >= end_distance || d2 >= end_distance) return ProcessPairResult::FAIL_EDIST;

  int end_separation = insert_avg - (d1 + d2);
  CidPair cids = {.cid1 = aln1.cid, .cid2 = aln2.cid};
  int end1 = aln1.orient == '+' ? 3 : 5;
  int end2 = aln2.orient == '+' ? 3 : 5;
  if (cids.cid1 < cids.cid2) {
    swap(end1, end2);
    swap(cids.cid1, cids.cid2);
  }
  Edge edge = {.cids = cids,
               .end1 = end1,
               .end2 = end2,
               .gap = end_separation,
               .support = 1,
               .aln_len = min(aln1.rstop - aln1.rstart, aln2.rstop - aln2.rstart),
               .aln_score = min(aln1.score1, aln2.score1),
               .edge_type = EdgeType::SPAN,
               .seq = "",
               .mismatch_error = false,
               .conflict_error = false,
               .excess_error = false,
               .short_aln = false,
               .gap_reads = {}};
  if (edge.gap > 0) {
    add_span_pos_gap_read(&edge, aln1);
    add_span_pos_gap_read(&edge, aln2);
  }
  _graph->add_or_update_edge(edge);
  return ProcessPairResult::SUCCESS;
}

// so the way meraculous spanner seems to work is that it only makes a pair out of the alignments with the shortest rstart.
void get_spans_from_alns(int insert_avg, int insert_stddev, int kmer_len, Alns &alns, CtgGraph *graph) {
  _graph = graph;
  BarrierTimer timer(__FILEFUNC__);
  IntermittentTimer t_get_alns(__FILENAME__ + string(":") + "get alns spans");
  ProgressBar progbar(alns.size(), "Adding edges to graph from spans");
  int64_t reject_5_trunc = 0;
  int64_t reject_3_trunc = 0;
  int64_t reject_uninf = 0;
  int64_t result_counts[(int)ProcessPairResult::SUCCESS + 1] = {0};
  int64_t num_pairs = 0;
  int64_t aln_i = 0;
  int read_len = 0;
  Aln prev_best_aln;
  string read_status = "", type_status = "", prev_read_status = "", prev_type_status = "";
  while (aln_i < (int64_t)alns.size()) {
    vector<Aln> alns_for_read;
    t_get_alns.start();
    get_all_alns_for_read(alns, aln_i, alns_for_read);
    t_get_alns.stop();
    progbar.update(aln_i);
    Aln best_aln;
    if (get_best_span_aln(insert_avg, insert_stddev, alns_for_read, best_aln, read_status, type_status, &reject_5_trunc,
                          &reject_3_trunc, &reject_uninf)) {
      if (!prev_best_aln.read_id.empty()) {
        read_len = best_aln.rlen;
        auto best_read_id_len = best_aln.read_id.length();
        auto prev_read_id_len = prev_best_aln.read_id.length();
        if (best_read_id_len == prev_read_id_len &&
            best_aln.read_id.compare(0, best_read_id_len - 1, prev_best_aln.read_id, 0, prev_read_id_len - 1) == 0) {
          if (best_aln.cid == prev_best_aln.cid) {
            result_counts[(int)ProcessPairResult::FAIL_SELF_LINK]++;
          } else {
            num_pairs++;
            auto res = process_pair(insert_avg, insert_stddev, prev_best_aln, best_aln, prev_type_status, type_status,
                                    prev_read_status, read_status);
            result_counts[(int)res]++;
          }
          // there will be no previous one next time
          prev_best_aln.read_id = "";
          prev_read_status = "";
          prev_type_status = "";
          continue;
        }
      }
    }
    prev_best_aln = best_aln;
    prev_read_status = read_status;
    prev_type_status = type_status;
  }
  progbar.done();
  barrier();
  t_get_alns.done_all();

  auto tot_alignments = reduce_one(alns.size(), op_fast_add, 0).wait();
  auto tot_num_pairs = reduce_one(num_pairs, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed ", tot_num_pairs, " pairs (out of ", tot_alignments, " alignments)\n");
  SLOG_VERBOSE("Rejected pairs:\n");
  SLOG_VERBOSE("  small:     ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_SMALL], op_fast_add, 0).wait(), "\n");
  SLOG_VERBOSE("  self link: ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_SELF_LINK], op_fast_add, 0).wait(), "\n");
  SLOG_VERBOSE("  edist:     ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_EDIST], op_fast_add, 0).wait(), "\n");
  SLOG_VERBOSE("  min gap:   ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_MIN_GAP], op_fast_add, 0).wait(), "\n");
  SLOG_VERBOSE("  uninf.:    ", reduce_one(reject_uninf, op_fast_add, 0).wait(), "\n");
#ifdef DUMP_LINKS
  string links_fname = "links-" + to_string(kmer_len) + ".spans.gz";
  get_rank_path(links_fname, rank_me());
  zstr::ofstream links_file(links_fname);
#endif
  int64_t num_spans_only = 0;
  int64_t num_pos_spans = 0;
  for (auto edge = _graph->get_first_local_edge(); edge != nullptr; edge = _graph->get_next_local_edge()) {
    if (edge->edge_type == EdgeType::SPAN) {
      num_spans_only++;
      int mean_gap_estimate = edge->gap / edge->support;
      int mean_offset = insert_avg - mean_gap_estimate;
      int clen1 = _graph->get_vertex(edge->cids.cid1)->clen;
      int clen2 = _graph->get_vertex(edge->cids.cid2)->clen;
      if (clen1 >= clen2) swap(clen1, clen2);
      // it appears that the full complex calculation (taken from meraculous) doesn't actually improve anything
      // compared to a simple setting based on the insert average
      edge->gap = estimate_gap_size(mean_offset, kmer_len, read_len, clen1, clen2, insert_avg, insert_stddev);
      //      edge->gap = mean_gap_estimate;
      if (edge->gap > 0) num_pos_spans++;
      // debug print in form comparable to mhm
      string ctg1 = "Contig" + to_string(edge->cids.cid1) + "." + to_string(edge->end1);
      string ctg2 = "Contig" + to_string(edge->cids.cid2) + "." + to_string(edge->end2);
      string link = (ctg1 < ctg2 ? ctg1 + "<=>" + ctg2 : ctg2 + "<=>" + ctg1);
#ifdef DUMP_LINKS
      links_file << "SPAN\t" << link << "\t0|" << edge->support << "\t" << edge->gap << "\t" << mean_gap_estimate << "\n";
#endif
    }
  }
#ifdef DUMP_LINKS
  links_file.close();
#endif
  barrier();
  auto tot_num_spans_only = reduce_one(num_spans_only, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(tot_num_spans_only, _graph->get_num_edges()), " spans\n");
  SLOG_VERBOSE("Found ", perc_str(reduce_one(num_pos_spans, op_fast_add, 0).wait(), tot_num_spans_only), " pos gap spans\n");
}
