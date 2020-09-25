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
#include "packed_reads.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

#define SIZE_TN_LOOKUP 257

void get_spans_from_alns(int insert_avg, int insert_stddev, int kmer_len, Alns &alns, CtgGraph *graph);
void get_splints_from_alns(Alns &alns, CtgGraph *graph);

static CtgGraph *_graph = nullptr;

#ifdef CALC_TNF
// return A C G T (0-3) or 4 if not a base
static int base_to_number(const char base) {
  int num = 0;
  switch (base) {
    case ('A'):
    case ('a'): num = 0; break;
    case ('C'):
    case ('c'): num = 1; break;
    case ('G'):
    case ('g'): num = 2; break;
    case ('T'):
    case ('t'): num = 3; break;
    default: num = 4;
  }
  return num;
}

// returns 0-255 for any valid tetamer with A C G or T and 256 for an invalid one
static int tn_to_number(const string &seq, int pos) {
  int num0 = base_to_number(seq[pos + 0]);
  int num1 = base_to_number(seq[pos + 1]);
  int num2 = base_to_number(seq[pos + 2]);
  int num3 = base_to_number(seq[pos + 3]);
  if (num0 == 4 || num1 == 4 || num2 == 4 || num3 == 4) return 256;
  return num0 + 4 * num1 + 16 * num2 + 64 * num3;
}

static void compute_tnfs(Contigs &ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  const string TN[] = {"GGTA", "AGCC", "AAAA", "ACAT", "AGTC", "ACGA", "CATA", "CGAA", "AAGT", "CAAA", "CCAG", "GGAC", "ATTA",
                       "GATC", "CCTC", "CTAA", "ACTA", "AGGC", "GCAA", "CCGC", "CGCC", "AAAC", "ACTC", "ATCC", "GACC", "GAGA",
                       "ATAG", "ATCA", "CAGA", "AGTA", "ATGA", "AAAT", "TTAA", "TATA", "AGTG", "AGCT", "CCAC", "GGCC", "ACCC",
                       "GGGA", "GCGC", "ATAC", "CTGA", "TAGA", "ATAT", "GTCA", "CTCC", "ACAA", "ACCT", "TAAA", "AACG", "CGAG",
                       "AGGG", "ATCG", "ACGC", "TCAA", "CTAC", "CTCA", "GACA", "GGAA", "CTTC", "GCCC", "CTGC", "TGCA", "GGCA",
                       "CACG", "GAGC", "AACT", "CATG", "AATT", "ACAG", "AGAT", "ATAA", "CATC", "GCCA", "TCGA", "CACA", "CAAC",
                       "AAGG", "AGCA", "ATGG", "ATTC", "GTGA", "ACCG", "GATA", "GCTA", "CGTC", "CCCG", "AAGC", "CGTA", "GTAC",
                       "AGGA", "AATG", "CACC", "CAGC", "CGGC", "ACAC", "CCGG", "CCGA", "CCCC", "TGAA", "AACA", "AGAG", "CCCA",
                       "CGGA", "TACA", "ACCA", "ACGT", "GAAC", "GTAA", "ATGC", "GTTA", "TCCA", "CAGG", "ACTG", "AAAG", "AAGA",
                       "CAAG", "GCGA", "AACC", "ACGG", "CCAA", "CTTA", "AGAC", "AGCG", "GAAA", "AATC", "ATTG", "GCAC", "CCTA",
                       "CGAC", "CTAG", "AGAA", "CGCA", "CGCG", "AATA"};
  // Palindromic sequences
  const string TNP[] = {"ACGT", "AGCT", "TCGA", "TGCA", "CATG", "CTAG", "GATC", "GTAC",
                        "ATAT", "TATA", "CGCG", "GCGC", "AATT", "TTAA", "CCGG", "GGCC"};
  const char bases[]{'A', 'C', 'G', 'T'};
  unordered_map<std::string, int> TN_map;
  unordered_set<std::string> TNP_map;
  // lookup table 0 - 255 of raw 4-mer to tetramer index in TNF
  array<int, SIZE_TN_LOOKUP> TN_lookup;
  // initialize the TN data structures
  for (size_t i = 0; i < nTNF; ++i) {
    TN_map[TN[i]] = i;
  }
  for (size_t i = 0; i < 16; ++i) {
    TNP_map.insert(TNP[i]);
  }
  TN_lookup[256] = nTNF;  // any non-base in the kmer
  string tnf_seq(4, 0);
  for (int i0 = 0; i0 < 4; i0++) {
    tnf_seq[0] = bases[i0];
    for (int i1 = 0; i1 < 4; i1++) {
      tnf_seq[1] = bases[i1];
      for (int i2 = 0; i2 < 4; i2++) {
        tnf_seq[2] = bases[i2];
        for (int i3 = 0; i3 < 4; i3++) {
          tnf_seq[3] = bases[i3];
          string tn = tnf_seq;
          int tn_number = tn_to_number(tnf_seq, 0);
          assert(tn_number <= 255);
          auto it = TN_map.find(tn);
          if (it != TN_map.end()) {
            TN_lookup[tn_number] = it->second;
            continue;
          }
          tn = revcomp(tn);
          if (TNP_map.find(tn) == TNP_map.end()) {  // if it is palindromic, then skip
            it = TN_map.find(tn);
            if (it != TN_map.end()) {
              TN_lookup[tn_number] = it->second;
            } else {
              WARN("Unknown TNF ", tn);
              continue;
            }
          } else {
            TN_lookup[tn_number] = nTNF;  // skip
          }
        }
      }
    }
  }

  ProgressBar progbar(ctgs.size(), "Computing TNFs for contigs");
  for (auto &ctg : ctgs) {
    for (int i = 0; i < nTNF; ++i) {
      ctg.tnf[i] = 0;
    }
    for (size_t i = 0; i < ctg.seq.length() - 3; ++i) {
      int tn_num = tn_to_number(ctg.seq, i);
      if (tn_num < 0 || tn_num >= SIZE_TN_LOOKUP) DIE("out of range ", tn_num);
      int tn_idx = TN_lookup[tn_num];
      if (tn_idx < nTNF) ++ctg.tnf[tn_idx];
    }
    // normalize to unit size (L2 norm)
    double rsum = 0;
    for (size_t c = 0; c < ctg.tnf.size(); ++c) {
      rsum += ctg.tnf[c] * ctg.tnf[c];
    }
    rsum = sqrt(rsum);
    for (size_t c = 0; c < ctg.tnf.size(); ++c) {
      ctg.tnf[c] /= rsum;
    }
    progbar.update();
  }
  progbar.done();
  barrier();
}
#endif

static void set_nbs() {
  BarrierTimer timer(__FILEFUNC__);
  {
    int64_t clen_excess = 0;
    int64_t num_excess_ctgs = 0;
    int max_excess_degree = 0;
    ProgressBar progbar(_graph->get_local_num_edges(), "Updating edges");
    // add edges to vertices
    for (auto edge = _graph->get_first_local_edge(); edge != nullptr; edge = _graph->get_next_local_edge()) {
      progbar.update();
      for (cid_t cid : {edge->cids.cid1, edge->cids.cid2}) {
        auto v = _graph->get_vertex(cid);
        assert(!edge->excess_error);
        if (v->end5.size() + v->end3.size() > CGRAPH_MAX_DEGREE) {
          edge->excess_error = true;
          clen_excess += v->clen;
          num_excess_ctgs++;
          max_excess_degree = std::max(max_excess_degree, v->clen);
          break;
        }
      }
      if (edge->excess_error) continue;
      // minor race condition here but it shouldn't lead to very high degrees
      _graph->add_vertex_nb(edge->cids.cid1, edge->cids.cid2, edge->end1);
      _graph->add_vertex_nb(edge->cids.cid2, edge->cids.cid1, edge->end2);
    }
    progbar.done();
    barrier();
    auto tot_clen_excess = reduce_one(clen_excess, op_fast_add, 0).wait();
    auto tot_num_excess_ctgs = reduce_one(num_excess_ctgs, op_fast_add, 0).wait();
    int all_max_excess_clen = reduce_one(max_excess_degree, op_fast_max, 0).wait();
    auto num_edges = _graph->get_num_edges();
    if (tot_num_excess_ctgs)
      SLOG_VERBOSE("Average excess clen ", ((double)tot_clen_excess / tot_num_excess_ctgs), " max ", all_max_excess_clen, "\n");
    int64_t num_excess_edges = reduce_one(_graph->purge_excess_edges(), op_fast_add, 0).wait();
    if (num_excess_edges) {
      if (rank_me() == 0 && tot_num_excess_ctgs == 0) SDIE("We have no excess ctgs, but ", num_excess_edges, " excess edges");
      SLOG_VERBOSE("Purged ", perc_str(num_excess_edges, num_edges), " excess degree edges\n");
    }
  }
  barrier();
}

static void add_vertices_from_ctgs(Contigs &ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  ProgressBar progbar(ctgs.size(), "Adding contig vertices to graph");
  for (auto &ctg : ctgs) {
    Vertex v = {.cid = ctg.id, .clen = (int)ctg.seq.length(), .depth = ctg.depth};
#ifdef TNF_PATH_RESOLUTION
    v.tnf = ctg.tnf;
#endif
    _graph->add_vertex(v, ctg.seq);
    progbar.update();
  }
  progbar.done();
  barrier();
  auto num_vertices = _graph->get_num_vertices();
  SLOG_VERBOSE("Added ", num_vertices, " vertices\n");
}

string get_consensus_seq(const vector<string> &seqs, int max_len) {
  static char bases[5] = {'A', 'C', 'G', 'T', 'N'};
  auto base_freqs = new int[max_len][4]();
  for (auto seq : seqs) {
    for (size_t i = 0; i < seq.size(); i++) {
      switch (seq[i]) {
        case 'A': base_freqs[i][0]++; break;
        case 'C': base_freqs[i][1]++; break;
        case 'G': base_freqs[i][2]++; break;
        case 'T': base_freqs[i][3]++; break;
        case 'N': break;
        default: WARN("Unknown base at pos ", i, " (", seq.size(), "): ", seq[i], "\n", seq);
      }
    }
  }
  string consensus_seq = "";
  for (int i = 0; i < max_len; i++) {
    int max_freq = 0;
    // start off with N
    int highest_idx = 4;
    for (int j = 0; j < 4; j++) {
      if (base_freqs[i][j] > max_freq) {
        max_freq = base_freqs[i][j];
        highest_idx = j;
      }
    }
    consensus_seq += bases[highest_idx];
  }
  delete[] base_freqs;
  // trim any Ns off the front
  consensus_seq.erase(0, consensus_seq.find_first_not_of('N'));
  return consensus_seq;
}

static string get_splint_edge_seq(int kmer_len, Edge *edge) {
  vector<string> seqs;
  // tail and end for checking primer matches
  int gap_size = edge->gap + 2 * (kmer_len - 1);
  for (auto gap_read : edge->gap_reads) {
    auto seq = _graph->get_read_seq(gap_read.read_name);
    if (seq == "") DIE("Could not find read seq for read ", gap_read.read_name, "\n");
    if (gap_read.gap_start < kmer_len) {
      // WARN("Positive gap overlap is less than kmer length, ", gap_read.gap_start, " < ", kmer_len, "\n");
      continue;
    }
    int rstart = gap_read.gap_start - (kmer_len - 1);
    if ((int)seq.length() < gap_size + rstart) {
      // WARN("seq length is less than sub string access with rstart ", rstart, ", gap ", gap_size, "\n",
      //     gap_read.read_name, " ", seq);
      continue;
    }
    string gap_seq = seq.substr(rstart, gap_size);
    if (edge->gap_reads.size() > 1) {
      string gap_seq_rc = revcomp(gap_seq);
      // these sequences should all be similar, so choosing the lexicographically least should ensure they have the same orientation
      if (gap_seq > gap_seq_rc) gap_seq = gap_seq_rc;
    }
    seqs.push_back(gap_seq);
  }
  return get_consensus_seq(seqs, gap_size);
}

static string get_span_edge_seq(int kmer_len, Edge *edge, bool tail) {
  string ctg_seq = "";
  cid_t cid = (tail ? edge->cids.cid1 : edge->cids.cid2);
  vector<string> seqs;
  int max_len = 0;
  for (auto gap_read : edge->gap_reads) {
    if (gap_read.cid != cid) continue;
    // sprintf(buf, "gs %3d rs %3d rp %3d %c ", gap_read.gap_start, gap_read.rstart, gap_read.rstop, gap_read.orient);
    auto gap_seq = _graph->get_read_seq(gap_read.read_name);
    if (gap_seq == "") DIE("Could not find read seq for read ", gap_read.read_name, "\n");
    if (tail) {
      if ((edge->end1 == 5 && gap_read.orient == '+') || (edge->end1 == 3 && gap_read.orient == '-')) gap_seq = revcomp(gap_seq);
      if (gap_read.gap_start > kmer_len) {
        gap_seq.erase(0, gap_read.gap_start - kmer_len);
        gap_read.gap_start = kmer_len;
      }
      // DBG_SPANS(buf, gap_seq, "\n");
    } else {
      if ((edge->end2 == 3 && gap_read.orient == '+') || (edge->end2 == 5 && gap_read.orient == '-')) gap_seq = revcomp(gap_seq);
      if (gap_read.gap_start + kmer_len < (int)gap_seq.size()) gap_seq.erase(gap_read.gap_start + kmer_len);
      // pad the front of the gap sequence with Ns to make them all the same length
      string offset_padding(1 + _graph->max_read_len - kmer_len - gap_read.gap_start, 'N');
      gap_seq = offset_padding + gap_seq;
      // DBG_SPANS(buf, gap_seq, "\n");
    }
    seqs.push_back(gap_seq);
    if ((int)gap_seq.size() > max_len) max_len = gap_seq.size();
  }
  if (seqs.empty()) {
    if (tail) {
      auto vertex = _graph->get_vertex(edge->cids.cid1);
      ctg_seq = _graph->get_vertex_seq(vertex->seq_gptr, vertex->clen);
      if (edge->end1 == 5) ctg_seq = revcomp(ctg_seq);
      int tail_len = vertex->clen - kmer_len;
      if (tail_len < 0) tail_len = 0;
      // ctg_seq = ctg_seq.substr(tail_len);
      ctg_seq.erase(0, tail_len);
      DBG_SPANS("TAIL contig", vertex->cid, ".", edge->end1, "\t", ctg_seq, "\n");
    } else {
      auto vertex = _graph->get_vertex(edge->cids.cid2);
      ctg_seq = _graph->get_vertex_seq(vertex->seq_gptr, vertex->clen);
      if (edge->end2 == 3) ctg_seq = revcomp(ctg_seq);
      // ctg_seq = ctg_seq.substr(0, kmer_len);
      ctg_seq.erase(kmer_len);
      DBG_SPANS("HEAD contig", vertex->cid, ".", edge->end2, "\t", ctg_seq, "\n");
    }
    return ctg_seq;
  }
  return get_consensus_seq(seqs, max_len);
}

static bool is_overlap_mismatch(int dist, int overlap) {
  return (dist > CGRAPH_SPAN_OVERLAP_MISMATCH_THRES || dist > overlap / 10);
}

static std::pair<int, int> min_hamming_dist(const string &s1, const string &s2, int max_overlap, int expected_overlap = -1) {
  int min_dist = max_overlap;
  if (expected_overlap != -1) {
    int min_dist = hamming_dist(tail(s1, expected_overlap), head(s2, expected_overlap));
    if (!is_overlap_mismatch(min_dist, expected_overlap)) return {min_dist, expected_overlap};
  }
  for (int d = std::min(max_overlap, (int)std::min(s1.size(), s2.size())); d >= 10; d--) {
    int dist = hamming_dist(tail(s1, d), head(s2, d));
    if (dist < min_dist) {
      min_dist = dist;
      expected_overlap = d;
      if (dist == 0) break;
    }
  }
  return {min_dist, expected_overlap};
}

static void parse_reads(unsigned kmer_len, const vector<PackedReads *> &packed_reads_list) {
  BarrierTimer timer(__FILEFUNC__);

  int64_t num_seqs_added = 0;
  int64_t num_reads = 0;
  unsigned max_read_len = 0;
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    string id, seq, quals;
    ProgressBar progbar(packed_reads->get_local_num_reads(), "Processing reads for gap sequences");
    while (true) {
      progress();
      if (!packed_reads->get_next_read(id, seq, quals)) break;
      progbar.update();
      // this happens when we have a placeholder entry because reads merged
      if (kmer_len > seq.length()) continue;
      if (_graph->update_read_seq(id, seq)) num_seqs_added++;
      if (seq.length() > max_read_len) max_read_len = seq.length();
      num_reads++;
    }
    progbar.done();
    barrier();
  }
  _graph->max_read_len = reduce_all(max_read_len, op_fast_max).wait();
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed a total of ", tot_num_reads, " reads, found max read length ", _graph->max_read_len, "\n");
  SLOG_VERBOSE("Extracted ", perc_str(reduce_one(num_seqs_added, op_fast_add, 0).wait(), tot_num_reads),
               " read sequences for pos gaps\n");

  int64_t num_pos_gaps = 0;
  int64_t num_pos_spans = 0;
  int64_t num_pos_spans_closed = 0;
  int64_t num_pos_spans_w_ns = 0;
  {
    ProgressBar progbar(_graph->get_local_num_edges(), "Fill pos gaps");
    for (auto edge = _graph->get_first_local_edge(); edge != nullptr; edge = _graph->get_next_local_edge()) {
      progbar.update();
      if (edge->gap <= 0) continue;
      num_pos_gaps++;
      if (edge->edge_type == EdgeType::SPAN) {
        num_pos_spans++;
        DBG_SPANS("SPAN pos gap ", edge->gap, " with ", edge->gap_reads.size(), " reads\n");
        string tail_seq = get_span_edge_seq(kmer_len, edge, true);
        if (tail_seq == "") continue;
        string head_seq = get_span_edge_seq(kmer_len, edge, false);
        if (head_seq == "") continue;
        DBG_SPANS("tail consensus         ", tail_seq, "\n");
        int offset = kmer_len + edge->gap - head_seq.size() + kmer_len;
        if (offset < 1) offset = 1;
        DBG_SPANS("head consensus         ", string(offset, ' '), head_seq, "\n");
        // now try to merge tail_seq and head_seq using the best (lowest hamming dist) overlap
        int min_len = min(tail_seq.size(), head_seq.size());
        int max_len = max(tail_seq.size(), head_seq.size());
        auto [min_dist, expected_overlap] = min_hamming_dist(tail_seq, head_seq, min_len);
        if (is_overlap_mismatch(min_dist, expected_overlap)) {
          min_dist = min_len;
          expected_overlap = -1;
          for (int i = 0; i < max_len - min_len; i++) {
            int dist = hamming_dist(substr_view(tail_seq, 0, min_len), substr_view(head_seq, 0, min_len));
            if (dist < min_dist) {
              min_dist = dist;
              expected_overlap = i;
            }
          }
          if (is_overlap_mismatch(min_dist, expected_overlap)) {
            DBG_SPANS("overlap mismatch: hdist ", min_dist, " best overlap ", expected_overlap, " original gap ", edge->gap, "\n");
            if (tail_seq.size() + head_seq.size() < 2 * kmer_len + edge->gap) {
              // the gap is not closed - fill with Ns
              int num_ns = 2 * kmer_len + edge->gap - tail_seq.size() - head_seq.size();
              edge->seq = tail_seq + string(num_ns, 'N') + head_seq;
              // remove one of either end because the later checking will use kmer_len - 1
              edge->seq = edge->seq.substr(1, edge->seq.size() - 2);
              DBG_SPANS("using orig gap ", edge->gap, " with ", num_ns, " Ns: ", edge->seq, "\n");
              num_pos_spans_w_ns++;
              num_pos_spans_closed++;
            }
            continue;
          }
        }
        num_pos_spans_closed++;
        int gap_size = tail_seq.size() + head_seq.size() - expected_overlap - 2 * kmer_len;
        DBG_SPANS("overlap is ", expected_overlap, " original gap is ", edge->gap, " corrected gap is ", gap_size, "\n");
        edge->gap = gap_size;
        DBG_SPANS(tail_seq, "\n");
        DBG_SPANS(string(tail_seq.size() - expected_overlap, ' '), head_seq, "\n");
        tail_seq.erase(tail_seq.size() - expected_overlap);
        edge->seq = tail_seq + head_seq;
        edge->seq = edge->seq.substr(1, edge->seq.size() - 2);
        DBG_SPANS(edge->seq, "\n");
        if (edge->seq.size() != 2 * (kmer_len - 1) + gap_size)
          WARN("fill mismatch ", edge->seq.size(), " != ", 2 * kmer_len + gap_size);
      } else {
        edge->seq = get_splint_edge_seq(kmer_len, edge);
      }
    }
    progbar.done();
    barrier();
  }
  auto tot_pos_gaps = reduce_one(num_pos_gaps, op_fast_add, 0).wait();
  SLOG_VERBOSE("Filled ", tot_pos_gaps, " positive gaps with ", _graph->get_num_read_seqs(), " reads\n");
  auto tot_pos_spans = reduce_one(num_pos_spans, op_fast_add, 0).wait();
  auto tot_pos_spans_closed = reduce_one(num_pos_spans_closed, op_fast_add, 0).wait();
  auto tot_pos_spans_w_ns = reduce_one(num_pos_spans_w_ns, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(tot_pos_spans, tot_pos_gaps), " positive spans, of which ",
               perc_str(tot_pos_spans_closed, tot_pos_spans), " were closed - ", perc_str(tot_pos_spans_w_ns, tot_pos_spans_closed),
               " with Ns\n");
}

static bool merge_end(Vertex *curr_v, const vector<cid_t> &nb_cids, vector<vector<cid_t> > &nb_cids_merged,
                      IntermittentTimer &t_merge_get_nbs, IntermittentTimer &t_merge_sort_nbs,
                      IntermittentTimer &t_merge_output_nbs) {
  nb_cids_merged.clear();

  // dead-end
  if (!nb_cids.size()) {
    DBG_BUILD("\tdead end\n");
    progress();
    return false;
  }

  struct NbPair {
    shared_ptr<Vertex> vertex;
    shared_ptr<Edge> edge;
  };

  t_merge_get_nbs.start();
  vector<NbPair> nbs;
  // select the neighbors that are valid for connections
  for (auto nb_cid : nb_cids) {
    auto nb_edge = _graph->get_edge_cached(curr_v->cid, nb_cid);
    // drop low support for high degree nodes - this reduces computational overhead and load imbalance, which can get severe
    if (nb_cids.size() > 50 && nb_edge->support < 2) continue;
    if (nb_cids.size() > 100 && nb_edge->support < 3) continue;
    if (nb_cids.size() > 150 && nb_edge->support < 4) continue;
    nbs.push_back({_graph->get_vertex_cached(nb_cid), nb_edge});
  }
  t_merge_get_nbs.stop();
  if (nbs.empty()) return false;
  // just one nb found
  if (nbs.size() == 1) {
    nb_cids_merged.push_back({nbs[0].vertex->cid});
    DBG_BUILD("\t", nbs[0].vertex->cid, " gap ", nbs[0].edge->gap, " len ", nbs[0].vertex->clen, " depth ", nbs[0].vertex->depth,
              "\n");
    return false;
  }
  t_merge_sort_nbs.start();
  // found multiple nbs, check for overlaps that can be merged
  // first, sort nbs by gap size
  sort(nbs.begin(), nbs.end(), [](const auto &elem1, const auto &elem2) { return elem1.edge->gap < elem2.edge->gap; });
  t_merge_sort_nbs.stop();

  // gather a vector of merged paths (there can be more than one because of forks)
  vector<vector<NbPair *> > all_next_nbs = {};
  // attempt to merge all neighbors as overlaps
  for (size_t i = 0; i < nbs.size(); i++) {
    NbPair *nb = &nbs[i];
    DBG_BUILD("\t", nb->vertex->cid, " gap ", nb->edge->gap, " len ", nb->vertex->clen, " depth ", nb->vertex->depth, "\n");
    if (i == 0) {
      all_next_nbs.push_back({nb});
      continue;
    }
    bool found_merge = false;
    for (auto &next_nbs : all_next_nbs) {
      progress();
      auto prev_nb = next_nbs.back();
      int g1 = -prev_nb->edge->gap;
      int g2 = -nb->edge->gap;
      int gdiff = g1 - g2;
      // check gap spacing to see if current nb overlaps previous nb
      if ((prev_nb->edge->gap == nb->edge->gap) || (gdiff >= prev_nb->vertex->clen - 10) ||
          (gdiff + nb->vertex->clen <= prev_nb->vertex->clen)) {
        DBG_BUILD("\tMerge conflict ", prev_nb->vertex->cid, " ", nb->vertex->cid, "\n");
        continue;
      }
      // now check that there is an edge from this vertex to the previous one
      auto intermediate_edge = _graph->get_edge_cached(nb->vertex->cid, prev_nb->vertex->cid);
      if (!intermediate_edge) {
        DBG_BUILD("\tNo edge found between ", prev_nb->vertex->cid, " and ", nb->vertex->cid, "\n");
        continue;
      }
      // now check the overlaps are correct
      DBG_BUILD("\tMERGE ", prev_nb->vertex->cid, " ", nb->vertex->cid, "\n");
      next_nbs.push_back(nb);
      found_merge = true;
      break;
    }
    if (!found_merge) all_next_nbs.push_back({nb});
  }
  if (all_next_nbs.empty()) DIE("empty all_next_nbs. How?\n");

  t_merge_output_nbs.start();
  // update the merged nbs and write out the merged paths found
  for (auto &next_nbs : all_next_nbs) {
    progress();
    nb_cids_merged.push_back({});
    DBG_BUILD("\tmerged path (len ", next_nbs.back()->vertex->clen + next_nbs.back()->edge->gap, ", depth ",
              next_nbs.back()->vertex->depth, "): ");
    for (auto &next_nb : next_nbs) {
      DBG_BUILD(next_nb->vertex->cid, " ");
      nb_cids_merged.back().push_back(next_nb->vertex->cid);
    }
    DBG_BUILD("\n");
  }
  t_merge_output_nbs.stop();

  return true;
}

static void merge_nbs() {
  barrier();
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_merges = 0;
  int64_t num_orphans = 0;
  int max_orphan_len = 0;
  int64_t num_nbs = 0, num_merged_nbs = 0, max_nbs = 0;
  double max_orphan_depth = 0;
  {
    IntermittentTimer t_merge_ends(__FILENAME__ + string(":") + "merge ends"),
        t_merge_get_nbs(__FILENAME__ + string(":") + "merge get nbs"),
        t_merge_sort_nbs(__FILENAME__ + string(":") + "merge sort nbs"),
        t_merge_output_nbs(__FILENAME__ + string(":") + "merge output nbs");
    ProgressBar progbar(_graph->get_local_num_vertices(), "Merge nbs");
    // mark all the vertices that have forks and the side of the forks. Note that in many cases what look like forks are
    // actually vertices that should be merged into a single neighbor
    for (auto v = _graph->get_first_local_vertex(); v != nullptr; v = _graph->get_next_local_vertex()) {
      DBG_BUILD(v->cid, " len ", v->clen, " depth ", v->depth, ":\n");
      if (v->end5.empty() && v->end3.empty()) {
        num_orphans++;
        if (v->clen > max_orphan_len) max_orphan_len = v->clen;
        if (v->depth > max_orphan_depth) max_orphan_depth = v->depth;
        DBG_BUILD("\torphan: ", v->cid, " len ", v->clen, " depth ", v->depth, "\n");
      }
      t_merge_ends.start();
      DBG_BUILD("  5-end:\n");
      if (merge_end(v, v->end5, v->end5_merged, t_merge_get_nbs, t_merge_sort_nbs, t_merge_output_nbs)) num_merges++;
      DBG_BUILD("  3-end:\n");
      if (merge_end(v, v->end3, v->end3_merged, t_merge_get_nbs, t_merge_sort_nbs, t_merge_output_nbs)) num_merges++;
      t_merge_ends.stop();
      int v_num_nbs = v->end5.size() + v->end3.size();
      num_nbs += v_num_nbs;
      num_merged_nbs += v->end5_merged.size() + v->end3_merged.size();
      if (v_num_nbs > max_nbs) max_nbs = v_num_nbs;
      progbar.update();
    }
    DBG("Number of nbs ", num_nbs, " avg degree ", ((double)num_nbs / _graph->get_local_num_vertices()), " merged ", num_merged_nbs,
        " avg degree ", ((double)num_merged_nbs / _graph->get_local_num_vertices()), " max degree ", max_nbs, "\n");
    progbar.done();
    t_merge_ends.done_all();
    t_merge_get_nbs.done_all();
    t_merge_sort_nbs.done_all();
    t_merge_output_nbs.done_all();
  }
  barrier();
  auto tot_merges = reduce_one(num_merges, op_fast_add, 0).wait();
  SLOG_VERBOSE("Merged ", perc_str(tot_merges, 2 * _graph->get_num_vertices()), " vertices\n");
  auto tot_orphans = reduce_one(num_orphans, op_fast_add, 0).wait();
  auto all_max_orphan_len = reduce_one(max_orphan_len, op_fast_max, 0).wait();
  auto all_max_orphan_depth = reduce_one(max_orphan_depth, op_fast_max, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(tot_orphans, _graph->get_num_vertices()), " orphaned vertices (no edges), max length ",
               all_max_orphan_len, ", max depth ", all_max_orphan_depth, "\n");
}

void build_ctg_graph(CtgGraph *graph, int insert_avg, int insert_stddev, int kmer_len, vector<PackedReads *> &packed_reads_list,
                     Contigs &ctgs, Alns &alns) {
  BarrierTimer timer(__FILEFUNC__);
  _graph = graph;
#ifdef TNF_PATH_RESOLUTION
  compute_tnfs(ctgs);
#endif
  add_vertices_from_ctgs(ctgs);
  get_splints_from_alns(alns, graph);
  get_spans_from_alns(insert_avg, insert_stddev, kmer_len, alns, graph);
#ifdef TNF_PATH_RESOLUTION
  _graph->compute_edge_tnfs();
#endif
  int64_t mismatched = 0, conflicts = 0, empty_spans = 0;
  _graph->purge_error_edges(&mismatched, &conflicts, &empty_spans);
  auto num_edges = _graph->get_num_edges();
  SLOG_VERBOSE("Purged edges:\n");
  SLOG_VERBOSE("  mismatched:  ", perc_str(reduce_one(mismatched, op_fast_add, 0).wait(), num_edges), "\n");
  SLOG_VERBOSE("  conflicts:   ", perc_str(reduce_one(conflicts, op_fast_add, 0).wait(), num_edges), "\n");
  SLOG_VERBOSE("  empty spans: ", perc_str(reduce_one(empty_spans, op_fast_add, 0).wait(), num_edges), "\n");
  barrier();
  set_nbs();
  // mark_short_aln_edges(max_kmer_len);
  parse_reads(kmer_len, packed_reads_list);
  merge_nbs();
}
