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

#include <fcntl.h>
#include <unistd.h>

#include <iostream>
#include <queue>
#include <upcxx/upcxx.hpp>

#include "contigs.hpp"
#include "ctg_graph.hpp"
#include "ssw.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/reduce_prefix.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

static CtgGraph *_graph = nullptr;

// used for building up scaffolds
struct Walk {
  int64_t len;
  int64_t start_clen;
  double depth;
  vector<pair<cid_t, Orient>> vertices;
};

struct ScaffVertex {
  cid_t cid;
  Orient orient;
  int depth;
  int len;
};

struct Scaffold {
  int64_t id;
  string seq;
  vector<ScaffVertex> vertices;
  vector<int> gaps;
  double depth;
};

struct WalkStats {
  int64_t num_steps, dead_ends, term_visited, term_no_candidate, term_multi_candidates;

  void print() {
    int64_t tot_dead_ends = reduce_one(dead_ends, op_fast_add, 0).wait();
    int64_t tot_steps = reduce_one(num_steps, op_fast_add, 0).wait();

    int64_t tot_visited = reduce_one(term_visited, op_fast_add, 0).wait();
    int64_t tot_no_candidate = reduce_one(term_no_candidate, op_fast_add, 0).wait();
    int64_t tot_multi_candidates = reduce_one(term_multi_candidates, op_fast_add, 0).wait();
    int64_t tot_terms = tot_dead_ends + tot_visited + tot_no_candidate + tot_multi_candidates;

    cout << setprecision(2) << fixed;
    SLOG_VERBOSE("Walks statistics:\n");
    SLOG_VERBOSE("  total walk steps:           ", tot_steps, "\n");
    SLOG_VERBOSE("  walk terminations:          ", perc_str(tot_terms, tot_steps), "\n");
    SLOG_VERBOSE("    dead ends:                ", perc_str(tot_dead_ends, tot_terms), "\n");
    SLOG_VERBOSE("    no viable candidates:     ", perc_str(tot_no_candidate, tot_terms), "\n");
    SLOG_VERBOSE("    multiple candidates:      ", perc_str(tot_multi_candidates, tot_terms), "\n");
    SLOG_VERBOSE("    already visited:          ", perc_str(tot_visited, tot_terms), "\n");
  }
};

struct GapStats {
  int64_t mismatched_splints, mismatched_spans, gaps, positive, unclosed, corrected_splints, corrected_spans, num_break_scaffs,
      num_excess_breaks, num_tolerance_breaks, num_Ns_breaks;

  void print() {
    int64_t tot_mismatched_splints = reduce_one(mismatched_splints, op_fast_add, 0).wait();
    int64_t tot_mismatched_spans = reduce_one(mismatched_spans, op_fast_add, 0).wait();
    int64_t tot_gaps = reduce_one(gaps, op_fast_add, 0).wait();
    int64_t tot_positive = reduce_one(positive, op_fast_add, 0).wait();
    int64_t tot_unclosed = reduce_one(unclosed, op_fast_add, 0).wait();
    int64_t tot_break_scaffs = reduce_one(num_break_scaffs, op_fast_add, 0).wait();
    int64_t tot_excess_breaks = reduce_one(num_excess_breaks, op_fast_add, 0).wait();
    int64_t tot_tolerance_breaks = reduce_one(num_tolerance_breaks, op_fast_add, 0).wait();
    int64_t tot_Ns_breaks = reduce_one(num_Ns_breaks, op_fast_add, 0).wait();

    cout << setprecision(2) << fixed;
    SLOG_VERBOSE("Gaps statistics:\n");
    SLOG_VERBOSE("  total:                    ", tot_gaps, "\n");
    SLOG_VERBOSE("  positive:                 ", perc_str(tot_positive, tot_gaps), "\n");
    SLOG_VERBOSE("    unclosed:               ", perc_str(tot_unclosed, tot_positive), "\n");
    SLOG_VERBOSE("  mismatched splints:       ", perc_str(tot_mismatched_splints, tot_gaps), "\n");
    SLOG_VERBOSE("  mismatched spans:         ", perc_str(tot_mismatched_spans, tot_gaps), "\n");
    //    SLOG_VERBOSE("  corrected splints:        ", perc_str(tot_corrected_splints, tot_gaps), "\n");
    //    SLOG_VERBOSE("  corrected spans:          ", perc_str(tot_corrected_spans, tot_gaps), "\n");
    SLOG_VERBOSE("  scaffold breaks:          ", perc_str(tot_break_scaffs, tot_gaps), "\n");
    SLOG_VERBOSE("    excess:                 ", perc_str(tot_excess_breaks, tot_break_scaffs), "\n");
    SLOG_VERBOSE("    tolerance:              ", perc_str(tot_tolerance_breaks, tot_break_scaffs), "\n");
    SLOG_VERBOSE("    too many Ns:            ", perc_str(tot_Ns_breaks, tot_break_scaffs), "\n");
  }
};

static bool is_overlap_mismatch(int dist, int overlap) {
  return (dist > CGRAPH_GAP_CLOSING_OVERLAP_MISMATCH_THRES || dist > overlap / 10);
}

static void get_ctgs_from_walks(int max_kmer_len, int kmer_len, int break_scaff_Ns, vector<Walk> &walks, Contigs &ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  // match, mismatch, gap opening, gap extending, ambiguious
  StripedSmithWaterman::Aligner ssw_aligner(ALN_MATCH_SCORE, ALN_MISMATCH_COST, ALN_GAP_OPENING_COST, ALN_GAP_EXTENDING_COST,
                                            ALN_AMBIGUITY_COST);
  StripedSmithWaterman::Filter ssw_filter;
  ssw_filter.report_cigar = false;
  GapStats gap_stats = {0};
  for (auto &walk : walks) {
    shared_ptr<Vertex> prev_v(nullptr);
    Contig ctg = {0};
    ctg.depth = walk.depth;
    for (size_t i = 0; i < walk.vertices.size(); i++) {
      bool break_scaffold = false;
      bool break_for_Ns = false;
      auto v = _graph->get_vertex(walk.vertices[i].first);
      auto orient = walk.vertices[i].second;
      auto seq = _graph->get_vertex_seq(v->seq_gptr, v->clen);
      if (orient == Orient::REVCOMP) seq = revcomp(seq);
      if (!prev_v) {
        // no previous vertex - the start of the scaffold
        ctg.seq = seq;
      } else {
        gap_stats.gaps++;
        auto edge = _graph->get_edge(v->cid, prev_v->cid);
        if (edge->gap > 0) {
          gap_stats.positive++;
          string gap_seq;
          auto gap_edge_seq = edge->seq;
          if (gap_edge_seq.size()) {
            // can only close when the gap fill sequence exists
            // check that gap filling seq matches properly, kmer_len - 1 on each side
            int len = kmer_len - 1;
            // first check the left
            int dist = hamming_dist(tail(ctg.seq, len), head(gap_edge_seq, len));
            if (is_overlap_mismatch(dist, len)) {
              auto orig_gap_seq = gap_edge_seq;
              // norm is bad, the revcomp may be good
              gap_edge_seq = revcomp(gap_edge_seq);
              dist = hamming_dist(tail(ctg.seq, len), head(gap_edge_seq, len));
              if (is_overlap_mismatch(dist, len)) {
                DBG_WALK("breaking pos gap after revcomp, orig:\n", orig_gap_seq, "\n");
                // the revcomp is also bad, break the scaffold
                break_scaffold = true;
              }
            }
            if (!break_scaffold) {
              // now check the right
              dist = hamming_dist(tail(gap_edge_seq, len), head(seq, len));
              if (is_overlap_mismatch(dist, len)) {
                // for debugging, check non-revcomp
                auto rc_gap_edge_seq = revcomp(gap_edge_seq);
                DBG_WALK("dist to revcomp ", dist, " and to non-revcomp ", hamming_dist(tail(rc_gap_edge_seq, len), head(seq, len)),
                         ":\n  ", tail(rc_gap_edge_seq, len), "\n  ", tail(gap_edge_seq, len), "\n  ", head(seq, len), "\n");
                // the left was fine, but the right is not, so break
                break_scaffold = true;
              } else {
                gap_seq = gap_edge_seq.substr(len, edge->gap);
              }
              // now break if too many Ns
              if (break_scaff_Ns > 0 && !break_scaffold && gap_seq.find(string(break_scaff_Ns, 'N')) != string::npos) {
                break_scaffold = true;
                break_for_Ns = true;
                gap_stats.num_Ns_breaks++;
              }
            }
          } else {  // gap sequence does not exist
            gap_stats.unclosed++;
            // fill with Ns - we could also break the scaffold here
            DBG_WALK("SPAN (", edge->cids.cid1, ", ", edge->cids.cid2, ") gap size ", edge->gap, "\n");
            // gap_seq = string(edge->gap, 'N');
            break_scaffold = true;
          }
          if (!break_scaffold) ctg.seq += gap_seq + seq;
        } else if (edge->gap < 0) {
          int gap_excess = -edge->gap;
          if (gap_excess > (int)ctg.seq.size() || gap_excess > (int)seq.size()) gap_excess = min(ctg.seq.size(), seq.size()) - 5;
          // FIXME: this should be a full SSW alignment check to deal with indels
          if (ctg.seq.compare(ctg.seq.length() - gap_excess, gap_excess, seq, 0, gap_excess) == 0) {
            DBG_WALK("perfect neg overlap ", gap_excess, "\n");
            ctg.seq += tail(seq, seq.size() - gap_excess);
          } else {
            int max_overlap = max(gap_excess + 20, max_kmer_len + 2);
            if (max_overlap > (int)ctg.seq.size()) max_overlap = ctg.seq.size();
            if (max_overlap > (int)seq.size()) max_overlap = seq.size();
            StripedSmithWaterman::Alignment ssw_aln;
            // make sure upcxx progress is done before starting alignment
            progress();
            discharge();
            ssw_aligner.Align(tail(ctg.seq, max_overlap).c_str(), head(seq, max_overlap).c_str(), max_overlap, ssw_filter, &ssw_aln,
                              max((int)(max_overlap / 2), 15));
            int left_excess = max_overlap - (ssw_aln.query_end + 1);
            int right_excess = ssw_aln.ref_begin;
            int aln_len = ssw_aln.query_end - ssw_aln.query_begin;
            DBG_WALK("SSW aln for gap ", gap_excess, " left clen ", ctg.seq.length(), " right clen ", seq.length(), " left begin ",
                     ssw_aln.query_begin, " left end ", ssw_aln.query_end + 1, " right begin ", ssw_aln.ref_begin, " right end ",
                     ssw_aln.ref_end + 1, " left excess ", left_excess, " right excess ", right_excess, " score ", ssw_aln.sw_score,
                     " score next best ", ssw_aln.sw_score_next_best, " aln len ", aln_len, "\n", tail(ctg.seq, max_overlap), "\n",
                     head(seq, max_overlap), "\n");
            if (left_excess > KLIGN_UNALIGNED_THRES || right_excess > KLIGN_UNALIGNED_THRES) {
              break_scaffold = true;
              gap_stats.num_excess_breaks++;
              DBG_WALK("break neg gap\n");
            } else if (ssw_aln.sw_score < aln_len - ALN_MISMATCH_COST * CGRAPH_MAX_MISMATCHES_THRES) {
              break_scaffold = true;
              gap_stats.num_tolerance_breaks++;
              DBG_WALK("break poor aln: score ", ssw_aln.sw_score, " < ", aln_len - ALN_MISMATCH_COST * CGRAPH_MAX_MISMATCHES_THRES,
                       "\n");
            } else {
              DBG_WALK("close neg gap trunc left at ", ctg.seq.length() - max_overlap + ssw_aln.query_end + 1,
                       " and from right at ", ssw_aln.ref_end + 1, "\n");
              ctg.seq.erase(ctg.seq.length() - max_overlap + ssw_aln.query_end + 1);
              ctg.seq += tail(seq, seq.size() - (ssw_aln.ref_end + 1));
            }
          }
        } else {
          // gap is exactly 0
          ctg.seq += seq;
        }
        if (break_scaffold) {
          gap_stats.num_break_scaffs++;
          DBG_WALK("break scaffold from ", prev_v->cid, " to ", v->cid, " gap ", edge->gap, " type ",
                   edge_type_str(edge->edge_type), " prev_v clen ", prev_v->clen, " curr_v clen ", v->clen, "\n");
          if (!break_for_Ns) {
            if (edge->edge_type == EdgeType::SPLINT)
              gap_stats.mismatched_splints++;
            else
              gap_stats.mismatched_spans++;
          }
          // save current scaffold
          ctgs.add_contig(ctg);
          // start new scaffold
          ctg.depth = v->depth;
          ctg.seq = seq;
        }
      }
      prev_v = v;
    }
    // done with all walk vertices
    if (ctg.seq != "") ctgs.add_contig(ctg);
  }
  barrier();
  gap_stats.print();
  // now get unique ids for all the contigs
  size_t num_ctgs = ctgs.size();
  auto fut = upcxx_utils::reduce_prefix(num_ctgs, upcxx::op_fast_add).then([num_ctgs, &ctgs](size_t my_prefix) {
    auto my_counter = my_prefix - num_ctgs;  // get my start
    for (auto it = ctgs.begin(); it != ctgs.end(); it++) it->id = my_counter++;
  });
  fut.wait();
  barrier();
}

static bool depth_match(double depth, double walk_depth) {
  double depth_diff = fabs(depth - walk_depth);
  double allowable_diff = CGRAPH_DEPTH_DIFF_THRES * walk_depth;
  if (allowable_diff > CGRAPH_MAX_DEPTH_DIFF) allowable_diff = CGRAPH_MAX_DEPTH_DIFF;
  if (allowable_diff < CGRAPH_MIN_DEPTH_DIFF) allowable_diff = CGRAPH_MIN_DEPTH_DIFF;
  return (depth_diff <= allowable_diff);
}

static vector<shared_ptr<Vertex>> get_vertex_list(vector<cid_t> &cids) {
  vector<shared_ptr<Vertex>> vertices;
  for (auto &cid : cids) vertices.push_back(_graph->get_vertex_cached(cid));
  return vertices;
}

static string vertex_list_to_cid_string(vector<shared_ptr<Vertex>> &vertices) {
  string s;
  for (auto &v : vertices) s += to_string(v->cid) + " ";
  return s;
}

static cid_t bfs_branch(shared_ptr<Vertex> curr_v, int end, double walk_depth) {
  queue<pair<shared_ptr<Vertex>, int>> q;
  HASH_TABLE<cid_t, bool> visited;

  vector<shared_ptr<Vertex>> frontier = {};

  q.push({curr_v, end});
  // nullptr is a level marker
  q.push({nullptr, 0});
  int search_level = 0;
  cid_t candidate = -1;
  while (!q.empty()) {
    progress();
    tie(curr_v, end) = q.front();
    q.pop();
    if (!curr_v) {
      if (q.empty()) return -1;
      search_level++;
      string offset(6 + search_level * 2, ' ');
      // break if the search level is too high, or if the queue size is too big
      if (search_level >= CGRAPH_MAX_SEARCH_LEVEL) {
        DBG_WALK(offset, "Reached max search level ", CGRAPH_MAX_SEARCH_LEVEL, ", stopping...\n");
        break;
      }
      if (q.size() >= CGRAPH_MAX_QUEUE_SIZE) {
        DBG_WALK(offset, "Reached max queue size ", q.size(), " > ", CGRAPH_MAX_QUEUE_SIZE, " stopping...\n");
        break;
      }
      q.push({nullptr, 0});
      DBG_WALK(offset, "* level ", search_level, "\n");
      continue;
    }
    string offset(6 + search_level * 2, ' ');
    auto nb_cids = (end == 5 ? curr_v->end5 : curr_v->end3);
    DBG_WALK(offset, curr_v->cid, " depth ", curr_v->depth, " length ", curr_v->clen, " num nbs ", nb_cids.size(), "\n");
    if (!depth_match(curr_v->depth, walk_depth) && curr_v->depth < walk_depth / 2) {
      DBG_WALK(offset, "-> vertex ", curr_v->cid, " depth is too low ", curr_v->depth, "\n");
      continue;
    }
    // terminate if visited before
    if (visited.find(curr_v->cid) != visited.end()) {
      DBG_WALK(offset, "-> vertex ", curr_v->cid, " is already visited\n");
      continue;
    }
    visited[curr_v->cid] = true;
    // terminate if a suitable vertex is found
    // if (depth_match(curr_v->depth, walk_depth) && curr_v->clen > 200) {
    if (depth_match(curr_v->depth, walk_depth)) {
      DBG_WALK(offset, "-> found candidate vertex ", curr_v->cid, "\n");
      frontier.push_back(curr_v);
      candidate = curr_v->cid;
      // short circuit
      break;
    }
    DBG_WALK(offset, "adding ", nb_cids.size(), " vertices to the queue\n");
    for (auto &nb_cid : nb_cids) {
      auto nb = _graph->get_vertex_cached(nb_cid);
      auto edge = _graph->get_edge_cached(curr_v->cid, nb->cid);
      DBG_WALK(offset, "added ", nb->cid, " to the queue\n");
      auto nb_end = (_graph->get_other_end(curr_v, nb, edge) == 3 ? 5 : 3);
      q.push({nb, nb_end});
    }
  }
  if (!frontier.empty() || !q.empty()) {
    while (!q.empty()) {
      auto elem = q.front();
      q.pop();
      if (elem.first) frontier.push_back(elem.first);
    }
#ifdef DEBUG
    string cids_str = "      frontier consists of: ";
    for (auto v : frontier) cids_str += to_string(v->cid) + " ";
    DBG_WALK(cids_str, "\n");
#endif
  }
  return candidate;
}

static vector<shared_ptr<Vertex>> search_for_next_nbs(int max_kmer_len, int kmer_len, shared_ptr<Vertex> curr_v, int end,
                                                      double walk_depth, WalkStats &stats, cid_t fwd_cid = -1) {
  stats.num_steps++;
  // get the nbs from the correct end
  auto nbs_cids = (end == 5 ? curr_v->end5_merged : curr_v->end3_merged);
  DBG_WALK_CONT("curr_v ", curr_v->cid, " depth ", curr_v->depth, " length ", curr_v->clen, " nbs ", nbs_cids.size(),
                " walk_depth ", walk_depth, "\n");
  if (nbs_cids.empty()) {
    stats.dead_ends++;
    DBG_WALK("    -> terminate: dead end\n");
    return {};
  }

  vector<shared_ptr<Vertex>> nb_vertices;
  for (auto nb_cids : nbs_cids) nb_vertices.push_back(_graph->get_vertex_cached(nb_cids.back()));
  vector<shared_ptr<Edge>> nb_edges;
  for (auto nb_cids : nbs_cids) nb_edges.push_back(_graph->get_edge_cached(curr_v->cid, nb_cids.back()));

  vector<pair<int, int>> candidate_branches;
  HASH_TABLE<cid_t, int> candidates;
  bool bulge = false;
  // candidate first search from each of the neighbors (branches)
  for (size_t i = 0; i < nb_vertices.size(); i++) {
    auto nb = nb_vertices[i];
    auto edge = nb_edges[i];
    cid_t candidate = -1;
    if (fwd_cid == nb->cid) {
      candidate = fwd_cid;
      DBG_WALK("      ", i, ". branch ", nb->cid, " edge support ", edge->support, ", fwd cid, accept candidate\n");
    } else {
      auto nb_end = _graph->get_other_end(curr_v, nb, edge);
      DBG_WALK("      ", i, ". branch ", nb->cid, " edge support ", edge->support, ", searching...\n");
      candidate = bfs_branch(nb, nb_end == 5 ? 3 : 5, walk_depth);
    }
    progress();
    if (candidate != -1) {
      auto it = candidates.find(candidate);
      if (it != candidates.end()) {
        bulge = true;
        DBG_WALK("      -> ", nb->cid, " (", i, ") is a bulge with branch ", it->second, "\n");
      }
      candidates[candidate] = i;
      candidate_branches.push_back({i, nb->clen});
      DBG_WALK("      -> ", nb->cid, " (", i, ") found candidate ", candidate, "\n");
    }
  }

  int branch_chosen = -1;

  if (candidate_branches.size() == 1) {
    branch_chosen = candidate_branches[0].first;
    DBG_WALK("      -> viable candidate found on branch ", candidate_branches[0].first, "\n");
  } else if (candidate_branches.size() > 1) {
    stats.term_multi_candidates++;
    DBG_WALK("      -> found ", candidate_branches.size(), " viable candidates\n");
    if (max_kmer_len > kmer_len) {
      // if one branch has much better aln len than the others, choose it
      int num_max_kmer_alns = 0;
      vector<pair<int, int>> nbs_aln_scores;
      for (auto candidate : candidate_branches) {
        auto edge = nb_edges[candidate.first];
        if (edge->aln_score >= max_kmer_len) {
          num_max_kmer_alns++;
          if (num_max_kmer_alns > 1) break;
        }
        nbs_aln_scores.push_back({edge->aln_score, candidate.first});
      }
      if (num_max_kmer_alns < 2) {
        sort(nbs_aln_scores.begin(), nbs_aln_scores.end(), [](auto &a, auto &b) { return a.first > b.first; });
        DBG_WALK("    -> best aln score branch is ", nbs_aln_scores[0].second, " (", nbs_aln_scores[0].first, "), next best is ",
                 nbs_aln_scores[1].second, " (", nbs_aln_scores[1].first, ")\n");
        if (num_max_kmer_alns == 1) {
          branch_chosen = nbs_aln_scores[0].second;
          DBG_WALK("    -> resolve only max aln len ", branch_chosen, "\n");
        } else {
          if (nbs_aln_scores[0].first >= 2 * nbs_aln_scores[1].first) {
            branch_chosen = nbs_aln_scores[0].second;
            DBG_WALK("    -> resolve best aln len ", branch_chosen, "\n");
          }
        }
      }
    }
#ifdef TNF_PATH_RESOLUTION
    if (branch_chosen == -1) {
      // if one branch has much higher TNF than another, choose it
      vector<pair<int, int>> nbs_tnf;
      for (auto candidate : candidate_branches) {
        auto edge = nb_edges[candidate.first];
        nbs_tnf.push_back({edge->tnf_prob, candidate.first});
      }
      sort(nbs_tnf.begin(), nbs_tnf.end(), [](auto &a, auto &b) { return a.first > b.first; });
      DBG_WALK("    -> highest TNF branch is ", nbs_tnf[0].second, " (", nbs_tnf[0].first, "), next best is ", nbs_tnf[1].second,
               " (", nbs_tnf[1].first, ")\n");
      if (nbs_tnf[0].first >= CGRAPH_MIN_TNF_CLEN && nbs_tnf[1].first < CGRAPH_MIN_TNF_CLEN) {
        branch_chosen = nbs_tnf[0].second;
        DBG_WALK("    -> resolve highest TNF ", branch_chosen, "\n");
      }
    }
#endif
    if (branch_chosen == -1) {
      // if one branch has much higher edge support than the others, choose it
      vector<pair<int, int>> nbs_support;
      for (auto candidate : candidate_branches) {
        auto edge = nb_edges[candidate.first];
        nbs_support.push_back({edge->support, candidate.first});
      }
      sort(nbs_support.begin(), nbs_support.end(), [](auto &a, auto &b) { return a.first > b.first; });
      DBG_WALK("    -> most supported branch is ", nbs_support[0].second, " (", nbs_support[0].first, "), next best is ",
               nbs_support[1].second, " (", nbs_support[1].first, ")\n");
      if (nbs_support[0].first >= CGRAPH_WALK_SUPPORT_THRES * nbs_support[1].first && nbs_support[0].first > 2) {
        branch_chosen = nbs_support[0].second;
        DBG_WALK("    -> resolve most supported ", branch_chosen, "\n");
      }
    }
  }

  vector<shared_ptr<Vertex>> next_nbs = {};
  if (branch_chosen != -1) {
    //    DBG_WALK("Branch chosen has TNF of ", nb_edges[branch_chosen]->tnf_prob, "\n");
    next_nbs = get_vertex_list(nbs_cids[branch_chosen]);
    // make sure that this list contains the fwd cid, if it is specified
    if (fwd_cid != -1) {
      bool found = false;
      for (auto &next_nb : next_nbs) {
        // check to see if the merged list for the back candidate contains the current vertex
        if (next_nb->cid == fwd_cid) {
          found = true;
          break;
        }
      }
      if (!found) {
        DBG_WALK("    -> chosen bwd path does not include fwd vertex ", fwd_cid, "\n");
        next_nbs = {};
      }
    }
  }
  if (next_nbs.empty())
    stats.term_no_candidate++;
  else
    DBG_WALK("    -> resolved: ", vertex_list_to_cid_string(next_nbs), "\n");
  return next_nbs;
}

static vector<Walk> do_walks(int max_kmer_len, int kmer_len, vector<pair<cid_t, int32_t>> &sorted_ctgs, WalkStats &walk_stats,
                             IntermittentTimer &next_nbs_timer) {
  auto is_visited = [](HASH_TABLE<cid_t, bool> &visited, shared_ptr<Vertex> v) {
    if (v->visited) return true;
    if (visited.find(v->cid) == visited.end()) return false;
    return true;
  };

  auto get_start_vertex = [&](vector<pair<cid_t, int32_t>> &sorted_ctgs, int64_t *ctg_pos) -> shared_ptr<Vertex> {
    while (*ctg_pos < (int64_t)sorted_ctgs.size()) {
      auto v = _graph->get_local_vertex(sorted_ctgs[*ctg_pos].first);
      (*ctg_pos)++;
      if (!v->visited) return v;
    }
    return nullptr;
  };

  int64_t sum_scaff_lens = 0;
  int64_t ctg_pos = 0;
  shared_ptr<Vertex> start_v;
  // temporarily store the scaffolds from this rank - some may get discarded due to conflicts
  vector<Walk> tmp_walks;
  _graph->clear_caches();
  // each rank does an independent set of walks over the graph until it has run out of start vertices
  while ((start_v = get_start_vertex(sorted_ctgs, &ctg_pos)) != nullptr) {
    DBG_WALK("start ", start_v->cid, " len ", start_v->clen, " depth ", start_v->depth, "\n");
    // store the walk in a double ended queue because we could start walking in the middle and so need to add to
    // either the front or back
    deque<pair<shared_ptr<Vertex>, Orient>> walk_vertices;
    // start walk backwards, going from 5 to 3 ends
    Dirn dirn = Dirn::BACKWARD;
    int end = 5;
    Orient orient = Orient::NORMAL;

    double walk_depth = start_v->depth;
    HASH_TABLE<cid_t, bool> visited;
    visited[start_v->cid] = true;

    walk_vertices.push_front({start_v, orient});
    auto curr_v = start_v;
    int64_t scaff_len = curr_v->clen;
    // do walk
    while (curr_v) {
      DBG_WALK("    search fwd: ");
      next_nbs_timer.start();
      auto next_nbs = search_for_next_nbs(max_kmer_len, kmer_len, curr_v, end, walk_depth, walk_stats);
      next_nbs_timer.stop();
      if (!next_nbs.empty()) {
        // we have possibly multiple next nbs in sequence
        // update the last one and reject if it has insufficient depth remaining
        if (is_visited(visited, next_nbs.back())) {
          walk_stats.term_visited++;
          DBG_WALK_CONT("    -> terminate: ", next_nbs.back()->cid, " is already visited\n");
          curr_v = nullptr;
        }
      } else {
        curr_v = nullptr;
      }
      if (curr_v) {
        bool join_resolved = false;
        auto next_nb = next_nbs.back();
        // resolve joins (backward forks) from next_nb
        int next_nb_end = _graph->get_other_end(curr_v, next_nb);
        if ((next_nb_end == 5 && next_nb->end5_merged.size() > 1) || (next_nb_end == 3 && next_nb->end3_merged.size() > 1)) {
          DBG_WALK("    search bwd: ");
          next_nbs_timer.start();
          auto back_nbs = search_for_next_nbs(max_kmer_len, kmer_len, next_nb, next_nb_end, walk_depth, walk_stats, curr_v->cid);
          next_nbs_timer.stop();
          if (!back_nbs.empty()) join_resolved = true;
        } else {
          DBG_WALK("    accept single bwd path\n");
          join_resolved = true;
        }
        if (join_resolved) {
          if (depth_match(next_nbs.back()->depth, walk_depth)) walk_depth = next_nbs.back()->depth;
          // add all vertices in merged path to the walk
          for (auto next_nb : next_nbs) {
            visited[next_nb->cid] = true;
            auto edge = _graph->get_edge_cached(curr_v->cid, next_nb->cid);
            if (!edge) DIE("edge not found\n");
            scaff_len += next_nb->clen + edge->gap;
            auto next_nb_end = _graph->get_other_end(curr_v, next_nb, edge);
            curr_v = next_nb;
            // if the ends are the same, change the orientation of the next vertex
            if (end == next_nb_end) orient = flip_orient(orient);
            end = (next_nb_end == 3 ? 5 : 3);
            if (dirn == Dirn::BACKWARD)
              walk_vertices.push_front({curr_v, orient});
            else
              walk_vertices.push_back({curr_v, orient});
          }
        } else {
          DBG_WALK("    -> terminate: join not resolved\n");
          curr_v = nullptr;
        }
      }
      if (!curr_v && dirn == Dirn::BACKWARD) {
        // backward walk terminated, change walk direction
        dirn = Dirn::FORWARD;
        end = 3;
        curr_v = start_v;
        walk_depth = curr_v->depth;
        orient = Orient::NORMAL;
        DBG_WALK("  switch to dirn FORWARD\n");
      }
    }  // walk loop
    sum_scaff_lens += scaff_len;
    // unique ids are generated later
    Walk walk = {.len = scaff_len, .start_clen = start_v->clen, .depth = walk_depth, .vertices = {}};
    for (auto &w : walk_vertices) walk.vertices.push_back({w.first->cid, w.second});
    tmp_walks.push_back(walk);
  }
  return tmp_walks;
}

static vector<pair<cid_t, int32_t>> sort_ctgs(int min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__);
  vector<pair<cid_t, int32_t>> sorted_ctgs;
  for (auto v = _graph->get_first_local_vertex(); v != nullptr; v = _graph->get_next_local_vertex()) {
    // don't start on a contig that has already been used
    if (v->visited) continue;
    // don't start walks on short contigs
    if (v->clen < min_ctg_len) continue;
    // don't start walks on a high depth contig, which is a contig that has at least 2x depth higher than its nb average
    // and doesn't have any nbs of higher depth

    if (v->depth > 10) {
      auto all_nb_cids = v->end5;
      all_nb_cids.insert(all_nb_cids.end(), v->end3.begin(), v->end3.end());
      double sum_nb_depths = 0;
      bool found_higher_depth = false;
      for (auto nb_cid : all_nb_cids) {
        auto nb = _graph->get_vertex_cached(nb_cid);
        if (nb->depth > v->depth) {
          found_higher_depth = true;
          break;
        }
        sum_nb_depths += nb->depth;
      }
      if (!found_higher_depth) {
        if (v->depth > 2.0 * sum_nb_depths / all_nb_cids.size()) continue;
      }
    }

    sorted_ctgs.push_back({v->cid, v->clen});
  }
  sort(sorted_ctgs.begin(), sorted_ctgs.end(), [](auto &a, auto &b) { return a.second > b.second; });
  return sorted_ctgs;
}

void walk_graph(CtgGraph *graph, int max_kmer_len, int kmer_len, int break_scaff_Ns, Contigs &ctgs) {
  // The general approach is to have each rank do walks starting from its local vertices only.
  // First, to prevent loops within a walk, the vertices visited locally are kept track of using a visited hash table.
  // Once walks starting from all local vertices have been completed, any conflicts (overlaps) between walks are resolved.
  // These are resolved in favor of the longest walks (and for ties, the highest numbered rank). The walks
  // that lose are discarded, and the whole process is repeated, since there are potentially left-over vertices freed
  // up when walks are dropped.
  // The vertices in winning walks are marked as visited in the vertex structure.
  // This is repeated until there are no more walks found.
  _graph = graph;

  BarrierTimer timer(__FILEFUNC__);
  vector<Walk> walks;
  WalkStats walk_stats = {0};
  int64_t num_rounds = 0;
  auto sorted_ctgs = sort_ctgs(2 * max_kmer_len);
  barrier();
  int64_t num_start_ctgs = reduce_one(sorted_ctgs.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Number of starting contigs: ", perc_str(num_start_ctgs, _graph->get_num_vertices()), "\n");
  // need to repeat the sets of walks because there may be a conflicts between walks across ranks, which results
  // in one of the walks being dropped. So this loop will repeat until no more scaffolds can be built.
  {
    IntermittentTimer next_nbs_timer(__FILENAME__ + string(":") + "next_nbs");
    IntermittentTimer walks_timer(__FILENAME__ + string(":") + "walks");
    while (true) {
      walks_timer.start();
      ProgressBar progbar(sorted_ctgs.size(), "Walking graph round " + to_string(num_rounds));
      auto tmp_walks = do_walks(max_kmer_len, kmer_len, sorted_ctgs, walk_stats, next_nbs_timer);
      progbar.done();
      walks_timer.stop();
      barrier();
      // now eliminate duplicate walks. Each vertex will get labeled with the rank that has the longest walk,
      // first set all the vertex fields to empty
      for (auto v = _graph->get_first_local_vertex(); v != nullptr; v = _graph->get_next_local_vertex()) {
        v->walk_rank = -1;
        v->walk_score = 0;
      }
      barrier();
      for (size_t walk_i = 0; walk_i < tmp_walks.size(); walk_i++) {
        auto walk = &tmp_walks[walk_i];
        // resolve conflict in favor of longest walk - this marks the walk the vertex belongs to
        for (auto &w : walk->vertices) _graph->update_vertex_walk(w.first, walk->len, walk_i);
        // resolve in favor of longest starting ctg - reduces the msa on synth64d around 19%, with a reduction in ctgy too)
        // for (auto &w : walk->vertices) _graph->update_vertex_walk(w.first, walk->start_clen, walk_i);
      }
      barrier();
      int num_walks_added = 0;
      // now drop all walks where any vertex's rank does not match this rank
      for (int walk_i = 0; walk_i < (int)tmp_walks.size(); walk_i++) {
        auto walk = &tmp_walks[walk_i];
        bool add_walk = true;
        for (auto &w : walk->vertices) {
          auto v = _graph->get_vertex(w.first);
          if (v->walk_rank != rank_me() || v->walk_i != walk_i) {
            add_walk = false;
            break;
          }
        }
        if (add_walk) {
          num_walks_added++;
          // update depth remaining
          for (auto &w : walk->vertices) _graph->set_vertex_visited(w.first);
          walks.push_back(*walk);
        }
      }
      barrier();
      auto tot_walks_added = reduce_all(num_walks_added, op_fast_add).wait();
      if (tot_walks_added == 0) break;
      // SLOG_VERBOSE("Walk round ", num_rounds, " found ", tot_walks_added, " new walks\n");
      num_rounds++;
      if (num_rounds > rank_n() * 5) {
        SWARN("breaking contig graph walk on high count\n");
        break;
      }
    }  // loop until no more walks are found
    next_nbs_timer.done_all();
    walks_timer.done_all();
  }
  barrier();
  // now add any unvisited to the walks
  int64_t num_unvisited = 0;
  int64_t max_unvisited_len = 0;
  int64_t unvisited_len = 0;
  for (auto v = _graph->get_first_local_vertex(); v != nullptr; v = _graph->get_next_local_vertex()) {
    if (!v->visited) {
      num_unvisited++;
      if (v->clen > max_unvisited_len) max_unvisited_len = v->clen;
      unvisited_len += v->clen;
      Walk walk = {.len = v->clen, .start_clen = v->clen, .depth = v->depth, .vertices = {{v->cid, Orient::NORMAL}}};
      walks.push_back(walk);
      // DBG_WALK("unvisited ", v->cid, "\n");
    }
  }
  barrier();
  auto tot_unvisited = reduce_all(num_unvisited, op_fast_add).wait();
  auto tot_max_unvisited_len = reduce_all(max_unvisited_len, op_fast_max).wait();
  auto tot_unvisited_len = reduce_all(unvisited_len, op_fast_add).wait();
  if (tot_unvisited)
    SLOG_VERBOSE("Didn't visit ", tot_unvisited, " vertices, max len ", tot_max_unvisited_len, " total length ", tot_unvisited_len,
                 "\n");
  walk_stats.print();
  get_ctgs_from_walks(max_kmer_len, kmer_len, break_scaff_Ns, walks, ctgs);
  barrier();
}
