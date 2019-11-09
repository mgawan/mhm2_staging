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

#define DUMP_LINKS

void add_pos_gap_read(Edge* edge, Aln &aln);

static CtgGraph *_graph = nullptr;

// functions for evaluating the gap correction analytically assuming Gaussian insert size distribution
// this is taken from meraculaus bmaToLinks.pl 

int erf(int x) {
  int absX = (x < 0) ? -x : x;
  int t = 1 / (1 + 0.5 * absX);
  int t2 = t * t;
  int t3 = t * t2;
  int t4 = t * t3;
  int t5 = t * t4;
  int t6 = t * t5;
  int t7 = t * t6;
  int t8 = t * t7;
  int t9 = t * t8;
  int a0 = -1.26551223;
  int a1 = 1.00002368;
  int a2 = 0.37409196;
  int a3 = 0.09678418;
  int a4 = -0.18628806;
  int a5 = 0.27886807;
  int a6 = -1.13520398;
  int a7 = 1.48851587;
  int a8 = -0.82215223;
  int a9 = 0.17087277;
  int tau = t * exp(-x * x + a0 + a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5 + a6 * t6 + a7 * t7 + a8 * t8 + a9 * t9);
  if (x < 0) return (tau - 1);
  else return (1 - tau);
}


int G0(int a, int b, int mu, int sigma) {
  const int sqrtTwo = sqrt(2);
  const int pi = 3.14159265359;
  const int sqrtPi = sqrt(pi);

  int rt2sig = sqrtTwo * sigma;
  int erfa = erf((a - mu)/rt2sig);
  int erfb = erf((b - mu)/rt2sig);
  return (sqrtPi/sqrtTwo) * sigma * (erfb - erfa);
}


int G1(int a, int b, int mu, int sigma) {
  int za = (a - mu) / sigma;
  int zb = (b - mu) / sigma;
  int expa = exp(-0.5 * za * za);
  int expb = exp(-0.5 * zb * zb);
  int g0 = G0(a, b, mu, sigma);
  return (sigma * sigma * (expa-expb) + mu * g0);
}


int G2(int a, int b, int mu, int sigma) {
  int za = (a - mu) / sigma;
  int zb = (b - mu) / sigma;
  int expa = exp(-0.5 * za * za);
  int expb = exp(-0.5 * zb * zb);
  int g0 = G0(a, b, mu, sigma);
  int g1 = G1(a, b, mu, sigma);
  int sigma2 = sigma * sigma;
  return (sigma2 * g0 + mu * g1 + sigma2 * (a * expa - b * expb));
}


int mean_spanning_clone(int g, int k, int l, int c1, int c2, int mu, int sigma) {
  int x1 = g + 2 * k - 1;
  int x2 = g + c1 + l;
  int alpha = x2 - x1;
  int x3 = g + c2 + l;
  int x4 = x3 + alpha;
  int num = 0;
  int den = 0;
  int N1 = G2(x1, x2, mu, sigma) - x1 * G1(x1, x2, mu, sigma);
  int N2 = (x2 - x1) * G1(x2, x3, mu, sigma);
  int N3 = x4 * G1(x3, x4, mu, sigma) - G2(x3, x4, mu, sigma);
  int D1 = G1(x1, x2, mu, sigma) - x1 * G0(x1, x2, mu, sigma);
  int D2 = (x2 - x1) * G0(x2, x3, mu, sigma);
  int D3 = x4 * G0(x3, x4, mu, sigma) - G1(x3, x4, mu, sigma);
  num = N1 + N2 + N3;
  den = D1 + D2 + D3;
  if (den) {
    return num/den;
  } else {
    //WARN("mean_spanning_clone failed for (g,k,l,c1,c2,mu,sigma)");
    return 0;
  }
}


int estimate_gap_size(int meanAnchor, int k, int l, int c1, int c2, int mu, int sigma) {
  int gMax = mu+3 * sigma - 2 * k;
  int gMin = -(k - 2);
  int gMid = mu-meanAnchor;
  // Negative gap size padding disabled for metagenomes
  //if (gMid < gMin) gMid = gMin + 1;
  if (gMid > gMax) gMid = gMax - 1;
  int aMax = mean_spanning_clone(gMax, k, l, c1, c2, mu, sigma) - gMax;
  int aMin = mean_spanning_clone(gMin, k, l, c1, c2, mu, sigma) - gMin;
  int aMid = mean_spanning_clone(gMid, k, l, c1, c2, mu, sigma) - gMid;
  int deltaG = gMax-gMin;
  int iterations = 0;
  while (deltaG > 10) {
    iterations++;
    if (meanAnchor > aMid) {
      gMax = gMid;
      aMax = aMid;
      gMid = (gMid+gMin) / 2;
      aMid = mean_spanning_clone(gMid, k, l, c1, c2, mu, sigma) - gMid;
    } else if (meanAnchor < aMid) {
      gMin = gMid;
      aMin = aMid;
      gMid = (gMid+gMax) / 2;
      aMid = mean_spanning_clone(gMid, k, l, c1, c2, mu, sigma) - gMid;
    } else {
      break;
    }
    deltaG = gMax-gMin;
  }
  return gMid;
}


static bool get_best_span_aln(int insert_avg, int insert_stddev, vector<Aln> &alns, Aln &best_aln, string &read_status,
                              string &type_status, int64_t *reject_5_trunc, int64_t *reject_3_trunc, int64_t *reject_uninf) {
  if (alns.size() == 0) return false;
  int min_rstart = -1;
  read_status = "";
  // sort in order of highest to lowest aln scores
  sort(alns.begin(), alns.end(),
       [](auto &aln1, auto &aln2) {
         return aln1.score1 > aln2.score1;
       });
  for (auto aln : alns) {
    // Assess alignment for completeness (do this before scaffold coordinate conversion!)
    // Important: Truncations are performed before reverse complementation
    // and apply to the end of the actual read
    int rstart = aln.rstart + 1;
    int cstart = aln.cstart + 1;
    int unaligned_start = rstart - 1;
    int projected_start = (aln.orient == '+' ? cstart - unaligned_start : aln.cstop + unaligned_start);
    string start_status = "";
    int projected_off = 0;
    if (projected_start < 1) projected_off = 1 - projected_start;
    else if (projected_start > aln.clen) projected_off = projected_start - aln.clen;
    int missing_start_bases = unaligned_start - projected_off;
    if (unaligned_start == 0) {
      start_status = "FULL";
    } else if (projected_off > 0 && missing_start_bases < FIVE_PRIME_WIGGLE_ROOM) {
      start_status = "GAP";
    } else if (unaligned_start < FIVE_PRIME_WIGGLE_ROOM) {
      start_status = "INC";
    } else {
      (*reject_5_trunc)++;
      continue;
    }

    int unaligned_end = aln.rlen - aln.rstop;
    int projected_end = (aln.orient == '+' ? aln.cstop + unaligned_end : cstart - unaligned_end);
    string end_status = "";
    projected_off = 0;
    if (projected_end < 1) projected_off = 1 - projected_end;
    else if (projected_end > aln.clen) projected_off = projected_end - aln.clen;
    int missing_end_bases = unaligned_end - projected_off;
    if (unaligned_end == 0) {
      end_status = "FULL";
    } else if (projected_off > 0 && missing_end_bases < THREE_PRIME_WIGGLE_ROOM) {
      end_status = "GAP";
    } else if (unaligned_end < THREE_PRIME_WIGGLE_ROOM) {
      end_status = "INC";
    } else {
      (*reject_3_trunc)++;
      continue;
    }

    int endDistance = insert_avg + 3 * insert_stddev;
    // Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle)
    if ((aln.orient == '+' && projected_start > aln.clen - endDistance) ||
        (aln.orient == '-' && projected_start < endDistance)) {
      if (min_rstart == -1 || rstart < min_rstart) {
        read_status = start_status + "." + end_status + ".OUT";
        if (start_status == "GAP") type_status = "OUTGAP";
        else if (end_status == "GAP") type_status = "INTGAP";
        else type_status = "ANCHOR";
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


static string get_ctg_aln_str(Aln &aln, const string &read_status) {
  ostringstream os;
  os << "Contig" << aln.cid << "." << (aln.orient == '+' ? 3 : 5) << "\t["
     << read_status << " " << aln.read_id << " " << (aln.rstart + 1) << " " << aln.rstop << " " << aln.rlen
     << " Contig" << aln.cid << " " << (aln.cstart + 1) << " " << aln.cstop << " " << aln.clen << " "
     << (aln.orient == '+' ? "Plus" : "Minus") << "]";
  return os.str();
}


// gets all the alns for a single read, and returns true if there are more alns
static bool get_all_alns_for_read(Alns &alns, int64_t &i, vector<Aln> &alns_for_read) {
  string start_read_id = "";
  for (; i < alns.size(); i++) {
    Aln aln = alns.get_aln(i);
    // alns for a new read
    if (start_read_id != "" && aln.read_id != start_read_id) return true;
    alns_for_read.push_back(aln);
    start_read_id = aln.read_id;
  }
  return false;
}

enum class ProcessPairResult { FAIL_SMALL, FAIL_SELF_LINK, FAIL_EDIST, FAIL_MIN_GAP, SUCCESS };
 
ProcessPairResult process_pair(int insert_avg, int insert_stddev, Aln &aln1, Aln &aln2, const string &type_status1,
                               const string &type_status2, const string &read_status1, const string &read_status2,
                               zstr::ofstream &spans_file, int max_kmer_len) {
#ifdef DUMP_SPANS
  spans_file << "PAIR\t" << type_status1 << "." << type_status2 << "\t"
             << get_ctg_aln_str(aln1, read_status1) << "\t" << get_ctg_aln_str(aln2, read_status2) << endl;
#endif
  auto get_dist = [=](int &d, Aln &aln) -> bool {
    aln.rstart++;
    aln.cstart++;
    // Omit alignments to very short (potentially "popped-out") scaffolds
    if (aln.clen < insert_avg / 2) return false;
    d = (aln.orient == '+' ? (aln.clen - aln.cstart + 1) + (aln.rstart - 1) : aln.cstop + (aln.rstart - 1));
    return true;
  };

  // check for inappropriate self-linkage
  if (aln1.cid == aln2.cid) return ProcessPairResult::FAIL_SELF_LINK;
  int end_distance = insert_avg + 3 * insert_stddev;
  int d1, d2;
  if (!get_dist(d1, aln1)) return ProcessPairResult::FAIL_SMALL;
  if (!get_dist(d2, aln2)) return ProcessPairResult::FAIL_SMALL;
  if (d1 >= end_distance || d2 >= end_distance) return ProcessPairResult::FAIL_EDIST;
  
  int end_separation = insert_avg - (d1 + d2);
  //if (end_separation < -max_kmer_len * 2) return ProcessPairResult::FAIL_MIN_GAP;
  CidPair cids = { .cid1 = aln1.cid, .cid2 = aln2.cid };
  int end1 = aln1.orient == '+' ? 3 : 5;
  int end2 = aln2.orient == '+' ? 3 : 5;
  if (cids.cid1 < cids.cid2) {
    swap(end1, end2);
    swap(cids.cid1, cids.cid2);
  }
  Edge edge = { .cids = cids, .end1 = end1, .end2 = end2, .gap = end_separation, .support = 1,
                .aln_len = min(aln1.rstop - aln1.rstart, aln2.rstop - aln2.rstart), .aln_score = min(aln1.score1, aln2.score1),
                .edge_type = EdgeType::SPAN, .seq = "",
                .mismatch_error = false, .conflict_error = false, .excess_error = false, .short_aln = false, .gap_reads = {}};
  if (edge.gap > 0) {
    add_pos_gap_read(&edge, aln1);
    add_pos_gap_read(&edge, aln2);
  }
  _graph->add_or_update_edge(edge);
  return ProcessPairResult::SUCCESS;
}


// so the way meraculous spanner seems to work is that it only makes a pair out of the alignments with the shortest rstart. 
void run_spanner(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, Alns &alns, CtgGraph *graph) {
  _graph = graph;
  Timer timer(__func__);
  IntermittentTimer t_get_alns("get alns spanner");
  ProgressBar progbar(alns.size(), "Running spanner");
  int64_t reject_5_trunc = 0;
  int64_t reject_3_trunc = 0;
  int64_t reject_uninf = 0;
  int64_t result_counts[(int)ProcessPairResult::SUCCESS + 1] = {0};
  int64_t num_pairs = 0;
  int64_t aln_i = 0;
  Aln prev_best_aln = { .read_id = "" };
  string read_status = "", type_status = "", prev_read_status = "", prev_type_status = "";
  string spans_fname = "spanner-" + to_string(kmer_len) + ".spans.gz";
  get_rank_path(spans_fname, rank_me());
  zstr::ofstream spans_file(spans_fname);
  while (true) {
    vector<Aln> alns_for_read;
    t_get_alns.start();
    if (!get_all_alns_for_read(alns, aln_i, alns_for_read)) break;
    t_get_alns.stop();
    progbar.update(aln_i);
    Aln best_aln = { .read_id = "" };
    if (get_best_span_aln(insert_avg, insert_stddev, alns_for_read, best_aln, read_status, type_status, &reject_5_trunc,
                          &reject_3_trunc, &reject_uninf)) {
      if (!prev_best_aln.read_id.empty()) {
        auto read_id_len = best_aln.read_id.length();
        if (best_aln.read_id.compare(0, read_id_len - 2, prev_best_aln.read_id, 0, read_id_len - 2) == 0) {
          if (best_aln.cid == prev_best_aln.cid) {
            reject_uninf++;
          } else {
            num_pairs++;
            auto res = process_pair(insert_avg, insert_stddev, prev_best_aln, best_aln, prev_type_status, type_status,
                                    prev_read_status, read_status, spans_file, max_kmer_len);
            result_counts[(int)res]++;
            // there will be no previous one next time 
            prev_best_aln.read_id = "";
            prev_read_status = "";
            prev_type_status = "";
            continue;
          }
        }
      }
    }
    prev_best_aln = best_aln;
    prev_read_status = read_status;
    prev_type_status = type_status;
  }
  progbar.done();
  spans_file.close();
  barrier();

  auto tot_num_pairs = reduce_one(num_pairs, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed ", tot_num_pairs, " pairs and rejected:\n");
  SLOG_VERBOSE("  small:     ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_SMALL], op_fast_add, 0).wait(), "\n");
  SLOG_VERBOSE("  self link: ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_SELF_LINK], op_fast_add, 0).wait(), "\n");
  SLOG_VERBOSE("  edist:     ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_EDIST], op_fast_add, 0).wait(), "\n");
  SLOG_VERBOSE("  min gap:   ", reduce_one(result_counts[(int)ProcessPairResult::FAIL_MIN_GAP], op_fast_add, 0).wait(), "\n");
  string links_fname = "links-" + to_string(kmer_len) + ".spans.gz";
  get_rank_path(links_fname, rank_me());
  zstr::ofstream links_file(links_fname);
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
      edge->gap = estimate_gap_size(mean_offset, kmer_len, TMP_READ_LEN, clen1, clen2, insert_avg, insert_stddev);
      
      if (edge->gap > 0) num_pos_spans++;
      // debug print in form comparable to mhm
      string ctg1 = "Contig" + to_string(edge->cids.cid1) + "." + to_string(edge->end1);
      string ctg2 = "Contig" + to_string(edge->cids.cid2) + "." + to_string(edge->end2);
      string link = (ctg1 < ctg2 ? ctg1 + "<=>" + ctg2 : ctg2 + "<=>" + ctg1);
#ifdef DUMP_LINKS
      links_file << "SPAN\t" << link << "\t0|" << edge->support << "\t" << edge->gap << "\t" << mean_gap_estimate << endl;
#endif
    }
  }
  links_file.close();
  barrier();
  auto tot_num_spans_only = reduce_one(num_spans_only, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(tot_num_spans_only, _graph->get_num_edges()), " spans\n");
  SLOG_VERBOSE("Found ", perc_str(reduce_one(num_pos_spans, op_fast_add, 0).wait(), tot_num_spans_only), " pos gap spans\n");
}
