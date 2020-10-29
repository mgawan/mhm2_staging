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

#include "ctg_graph.hpp"

#include <stdarg.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <upcxx/upcxx.hpp>

#include "contigs.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"

using std::endl;
using std::get;
using std::istream;
using std::istringstream;
using std::make_shared;
using std::make_tuple;
using std::max;
using std::min;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::tuple;
using std::vector;

using namespace upcxx_utils;

bool GapRead::operator==(const GapRead &other) const { return (read_name == other.read_name); }

bool CidPair::operator==(const CidPair &other) const { return (cid1 == other.cid1 && cid2 == other.cid2); }

bool CidPair::operator!=(const CidPair &other) const { return (cid1 != other.cid1 || cid2 != other.cid2); }

ostream &operator<<(ostream &os, const CidPair &cids) {
  os << "(" << cids.cid1 << ", " << cids.cid2 << ")";
  return os;
}

string edge_type_str(EdgeType edge_type) {
  switch (edge_type) {
    case EdgeType::SPLINT: return "splint";
    case EdgeType::SPAN: return "span";
    default: return "unknown";
  }
}

size_t CtgGraph::get_vertex_target_rank(cid_t cid) { return std::hash<cid_t>{}(cid) % upcxx::rank_n(); }

size_t CtgGraph::get_edge_target_rank(CidPair &cids) { return std::hash<CidPair>{}(cids) % upcxx::rank_n(); }

size_t CtgGraph::get_read_target_rank(const string &r) { return std::hash<string>{}(r) % upcxx::rank_n(); }

#ifdef TNF_PATH_RESOLUTION
double CtgGraph::calc_tnf_dist(shared_ptr<Vertex> v1, shared_ptr<Vertex> v2) {
  auto tnf1 = &v1->tnf;
  auto tnf2 = &v2->tnf;
  // EXP(pre_prob) <= 9 yields prob >= 0.1, so pre_prob <= LOG(9.0);
  const double floor_prob = 0.1;
  const double floor_pre_prob = log((1.0 / floor_prob) - 1.0);
  double d = 0;
  for (size_t i = 0; i < nTNF; ++i) {
    d += pow((*tnf1)[i] - (*tnf2)[i], 2);  // euclidean double
  }
  d = sqrt(d);
  double log_c1 = log10(min(v1->clen, 500000));
  double log_c2 = log10(min(v2->clen, 500000));
  double lw11 = std::min(log_c1, log_c2);
  double lw21 = std::max(log_c1, log_c2);
  double lw12 = lw11 * lw11;
  double lw13 = lw12 * lw11;
  double lw14 = lw13 * lw11;
  double lw15 = lw14 * lw11;
  double lw16 = lw15 * lw11;
  double lw17 = lw16 * lw11;
  double lw22 = lw21 * lw21;
  double lw23 = lw22 * lw21;
  double lw24 = lw23 * lw21;
  double lw25 = lw24 * lw21;
  double lw26 = lw25 * lw21;
  double b = 46349.1624324381 + -76092.3748553155 * lw11 + -639.918334183 * lw21 + 53873.3933743949 * lw12 +
             -156.6547554844 * lw22 + -21263.6010657275 * lw13 + 64.7719132839 * lw23 + 5003.2646455284 * lw14 +
             -8.5014386744 * lw24 + -700.5825500292 * lw15 + 0.3968284526 * lw25 + 54.037542743 * lw16 + -1.7713972342 * lw17 +
             474.0850141891 * lw11 * lw21 + -23.966597785 * lw12 * lw22 + 0.7800219061 * lw13 * lw23 + -0.0138723693 * lw14 * lw24 +
             0.0001027543 * lw15 * lw25;
  double c = -443565.465710869 + 718862.10804858 * lw11 + 5114.1630934534 * lw21 + -501588.206183097 * lw12 +
             784.4442123743 * lw22 + 194712.394138513 * lw13 + -377.9645994741 * lw23 + -45088.7863182741 * lw14 +
             50.5960513287 * lw24 + 6220.3310639927 * lw15 + -2.3670776453 * lw25 + -473.269785487 * lw16 + 15.3213264134 * lw17 +
             -3282.8510348085 * lw11 * lw21 + 164.0438603974 * lw12 * lw22 + -5.2778800755 * lw13 * lw23 +
             0.0929379305 * lw14 * lw24 + -0.0006826817 * lw15 * lw25;
  // logistic model
  // prob = 1.0 / (1 + EXP(-(b + c * d)));
  // if (prob >= .1)  //second logistic model
  double pre_prob = -(b + c * d);
  // pre_prob <= LOG(9.0) yields prob > 0.1, so use second logistic model
  double prob = (pre_prob <= floor_pre_prob ? floor_prob : 1.0 / (1 + exp(pre_prob)));
  if (prob >= floor_prob) {  // second logistic model
    b = 6770.9351457442 + -5933.7589419767 * lw11 + -2976.2879986855 * lw21 + 3279.7524685865 * lw12 + 1602.7544794819 * lw22 +
        -967.2906583423 * lw13 + -462.0149190219 * lw23 + 159.8317289682 * lw14 + 74.4884405822 * lw24 + -14.0267151808 * lw15 +
        -6.3644917671 * lw25 + 0.5108811613 * lw16 + 0.2252455343 * lw26 + 0.965040193 * lw12 * lw22 + -0.0546309127 * lw13 * lw23 +
        0.0012917084 * lw14 * lw24 + -1.14383e-05 * lw15 * lw25;
    c = 39406.5712626297 + -77863.1741143294 * lw11 + 9586.8761567725 * lw21 + 55360.1701572325 * lw12 + -5825.2491611377 * lw22 +
        -21887.8400068324 * lw13 + 1751.6803621934 * lw23 + 5158.3764225203 * lw14 + -290.1765894829 * lw24 +
        -724.0348081819 * lw15 + 25.364646181 * lw25 + 56.0522105105 * lw16 + -0.9172073892 * lw26 + -1.8470088417 * lw17 +
        449.4660736502 * lw11 * lw21 + -24.4141920625 * lw12 * lw22 + 0.8465834103 * lw13 * lw23 + -0.0158943762 * lw14 * lw24 +
        0.0001235384 * lw15 * lw25;
    // prob = 1.0 / (1 + EXP(-(b + c * d)));
    // prob = prob < .1 ? .1 : prob;
    // EXP(pre_prob) <= 9 yields prob >= 0.1, so pre_prob <= LOG(9.0) to calculate, otherwise use the floor
    pre_prob = -(b + c * d);
    prob = pre_prob <= floor_pre_prob ? 1.0 / (1 + exp(pre_prob)) : floor_prob;
  }
  return prob;
}
#endif

CtgGraph::CtgGraph()
    : vertices({})
    , edges({})
    , read_seqs({})
    , vertex_cache({})
    , edge_cache({}) {
  vertex_cache.reserve(CGRAPH_MAX_CACHE_SIZE);
  edge_cache.reserve(CGRAPH_MAX_CACHE_SIZE);
}

void CtgGraph::clear() {
  for (auto it = vertices->begin(); it != vertices->end();) {
    it = vertices->erase(it);
  }
  for (auto it = edges->begin(); it != edges->end();) {
    it = edges->erase(it);
  }
}

CtgGraph::~CtgGraph() { clear(); }

int64_t CtgGraph::get_num_vertices(bool all) {
  if (!all)
    return upcxx::reduce_one(vertices->size(), upcxx::op_fast_add, 0).wait();
  else
    return upcxx::reduce_all(vertices->size(), upcxx::op_fast_add).wait();
}

int64_t CtgGraph::get_local_num_vertices(void) { return vertices->size(); }

int64_t CtgGraph::get_num_edges(bool all) {
  if (!all)
    return upcxx::reduce_one(edges->size(), upcxx::op_fast_add, 0).wait();
  else
    return upcxx::reduce_all(edges->size(), upcxx::op_fast_add).wait();
}

int64_t CtgGraph::get_local_num_edges(void) { return edges->size(); }

shared_ptr<Vertex> CtgGraph::get_vertex(cid_t cid) {
  size_t target_rank = get_vertex_target_rank(cid);
  if (target_rank == (size_t)upcxx::rank_me()) {
    upcxx::progress();
    const auto it = vertices->find(cid);
    if (it == vertices->end()) return nullptr;
    return make_shared<Vertex>(it->second);
  }
  return upcxx::rpc(
             target_rank,
             [](vertex_map_t &vertices, cid_t cid) {
               const auto it = vertices->find(cid);
               if (it == vertices->end()) return Vertex({.cid = -1});
               return it->second;
             },
             vertices, cid)
      .then([](Vertex v) -> shared_ptr<Vertex> {
        if (v.cid == -1) return nullptr;
        return make_shared<Vertex>(v);
      })
      .wait();
}

int CtgGraph::get_vertex_clen(cid_t cid) {
  return upcxx::rpc(
             get_vertex_target_rank(cid),
             [](vertex_map_t &vertices, cid_t cid) {
               const auto it = vertices->find(cid);
               if (it == vertices->end()) return -1;
               return it->second.clen;
             },
             vertices, cid)
      .wait();
}

shared_ptr<Vertex> CtgGraph::get_local_vertex(cid_t cid) {
  const auto it = vertices->find(cid);
  if (it == vertices->end()) return nullptr;
  return make_shared<Vertex>(it->second);
}

void CtgGraph::set_vertex_visited(cid_t cid) {
  upcxx::rpc(
      get_vertex_target_rank(cid),
      [](vertex_map_t &vertices, cid_t cid) {
        const auto it = vertices->find(cid);
        if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
        auto v = &it->second;
        v->visited = true;
      },
      vertices, cid)
      .wait();
}

void CtgGraph::update_vertex_walk(cid_t cid, int walk_score, int walk_i) {
  upcxx::rpc(
      get_vertex_target_rank(cid),
      [](vertex_map_t &vertices, cid_t cid, int walk_score, int walk_i, int myrank) {
        const auto it = vertices->find(cid);
        if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
        auto v = &it->second;
        if (myrank == v->walk_rank) {
          // same rank, select in favor of highest score
          if (walk_score > v->walk_score) {
            v->walk_score = walk_score;
            v->walk_i = walk_i;
          }
        } else {
          // different rank, select highest score, break ties with highest rank
          if ((walk_score == v->walk_score && myrank > v->walk_rank) || walk_score > v->walk_score) {
            v->walk_score = walk_score;
            v->walk_rank = myrank;
            v->walk_i = walk_i;
          }
        }
      },
      vertices, cid, walk_score, walk_i, upcxx::rank_me())
      .wait();
}

Vertex *CtgGraph::get_first_local_vertex() {
  vertex_iter = vertices->begin();
  if (vertex_iter == vertices->end()) return nullptr;
  auto v = &vertex_iter->second;
  vertex_iter++;
  return v;
}

Vertex *CtgGraph::get_next_local_vertex() {
  if (vertex_iter == vertices->end()) return nullptr;
  auto v = &vertex_iter->second;
  vertex_iter++;
  return v;
}

void CtgGraph::add_vertex(Vertex &v, const string &seq) {
  v.clen = seq.length();
  v.seq_gptr = upcxx::allocate<char>(v.clen + 1);
  strcpy(v.seq_gptr.local(), seq.c_str());
  upcxx::rpc(
      get_vertex_target_rank(v.cid),
      [](vertex_map_t &vertices, Vertex v) {
        v.visited = false;
        vertices->insert({v.cid, v});
      },
      vertices, v)
      .wait();
}

void CtgGraph::add_vertex_nb(cid_t cid, cid_t nb, char end) {
  upcxx::rpc(
      get_vertex_target_rank(cid),
      [](vertex_map_t &vertices, cid_t cid, cid_t nb, int end) {
        const auto it = vertices->find(cid);
        if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
        auto v = &it->second;
#ifdef DEBUG
        // sanity checks
        for (auto &prev_v : v->end5) {
          if (prev_v == nb) {
            WARN("end5 already includes nb ", nb);
            return;
          }
        }
        for (auto &prev_v : v->end3) {
          if (prev_v == nb) {
            WARN("end3 already includes nb ", nb);
            return;
          }
        }
#endif
        if (end == 5)
          v->end5.push_back(nb);
        else
          v->end3.push_back(nb);
      },
      vertices, cid, nb, end)
      .wait();
}

string CtgGraph::get_vertex_seq(upcxx::global_ptr<char> seq_gptr, int64_t seq_len) {
  char buf[seq_len + 1];
  upcxx::rget(seq_gptr, buf, seq_len + 1).wait();
  string s(buf);
  return s;
}

shared_ptr<Edge> CtgGraph::get_edge(cid_t cid1, cid_t cid2) {
  CidPair cids = {.cid1 = cid1, .cid2 = cid2};
  if (cid1 < cid2) std::swap(cids.cid1, cids.cid2);
  size_t target_rank = get_edge_target_rank(cids);
  if (target_rank == (size_t)upcxx::rank_me()) {
    upcxx::progress();
    const auto it = edges->find(cids);
    if (it == edges->end()) return nullptr;
    return make_shared<Edge>(it->second);
  }
  return upcxx::rpc(
             target_rank,
             [](edge_map_t &edges, CidPair cids) -> Edge {
               const auto it = edges->find(cids);
               if (it == edges->end()) return Edge({.cids = {-1, -1}});
               return it->second;
             },
             edges, cids)
      .then([](Edge edge) -> shared_ptr<Edge> {
        if (edge.cids.cid1 == -1 && edge.cids.cid2 == -1) return nullptr;
        return make_shared<Edge>(edge);
      })
      .wait();
}

Edge *CtgGraph::get_first_local_edge() {
  edge_iter = edges->begin();
  if (edge_iter == edges->end()) return nullptr;
  auto edge = &edge_iter->second;
  edge_iter++;
  return edge;
}

Edge *CtgGraph::get_next_local_edge() {
  if (edge_iter == edges->end()) return nullptr;
  auto edge = &edge_iter->second;
  edge_iter++;
  return edge;
}

void CtgGraph::add_or_update_edge(Edge &edge) {
  upcxx::rpc(
      get_edge_target_rank(edge.cids),
      [](edge_map_t &edges, Edge new_edge) {
        const auto it = edges->find(new_edge.cids);
        if (it == edges->end()) {
          // not found, always insert
          edges->insert({new_edge.cids, new_edge});
        } else {
          auto edge = &it->second;
          // always a failure
          if (edge->mismatch_error) return;
          if (edge->edge_type == EdgeType::SPLINT && new_edge.edge_type == EdgeType::SPAN) {
            DBG_BUILD("span confirms splint: ", edge->cids, "\n");
            edge->support++;
          } else {
            if (edge->edge_type == EdgeType::SPLINT && new_edge.edge_type == EdgeType::SPLINT) {
              // check for mismatches in gap size if they're both splints
              if (abs(new_edge.gap - edge->gap) > 2) {
                DBG_BUILD("gap mismatch for ", new_edge.cids, " ", new_edge.gap, " != ", edge->gap, "\n");
                edge->mismatch_error = true;
                // this edge will be dropped
                return;
              }
              edge->gap = min(new_edge.gap, edge->gap);
            }
            if (edge->edge_type == EdgeType::SPAN && new_edge.edge_type == EdgeType::SPAN) edge->gap += new_edge.gap;
            edge->support++;
            edge->aln_len = max(edge->aln_len, new_edge.aln_len);
            edge->aln_score = max(edge->aln_score, new_edge.aln_score);
            // count conflicts
            if (edge->end1 != new_edge.end1 || edge->end2 != new_edge.end2) edge->conflict_error = true;
            if (new_edge.gap > 0) {
              // add reads to positive gap for splints
              edge->gap_reads.insert(edge->gap_reads.end(), new_edge.gap_reads.begin(), new_edge.gap_reads.end());
            }
          }
        }
      },
      edges, edge)
      .wait();
}

void CtgGraph::purge_error_edges(int64_t *mismatched, int64_t *conflicts, int64_t *empty_spans) {
  for (auto it = edges->begin(); it != edges->end();) {
    auto edge = make_shared<Edge>(it->second);
    if (edge->mismatch_error) {
      (*mismatched)++;
      it = edges->erase(it);
    } else if (edge->edge_type == EdgeType::SPAN && edge->gap > 0 && !edge->gap_reads.size()) {
      // don't use positive span gaps without filler
      (*empty_spans)++;
      it = edges->erase(it);
    } else {
      // retain the conflicts to prevent falsely choosing a path when it should be a fork
      if (edge->conflict_error) (*conflicts)++;
      it++;
    }
  }
}

int64_t CtgGraph::purge_excess_edges() {
  int64_t excess = 0;
  for (auto it = edges->begin(); it != edges->end();) {
    auto edge = make_shared<Edge>(it->second);
    if (edge->excess_error) {
      excess++;
      it = edges->erase(it);
    } else {
      it++;
    }
  }
  return excess;
}

void CtgGraph::remove_nb(cid_t cid, int end, cid_t nb) {
  upcxx::rpc(
      get_vertex_target_rank(cid),
      [](vertex_map_t &vertices, cid_t cid, int end, cid_t nb) {
        const auto it = vertices->find(cid);
        if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
        auto v = &it->second;
        if (end == 3) {
          for (auto it = v->end3.begin(); it != v->end3.end(); it++) {
            if (*it == nb) {
              v->end3.erase(it);
              return;
            }
          }
        } else {
          for (auto it = v->end5.begin(); it != v->end5.end(); it++) {
            if (*it == nb) {
              v->end5.erase(it);
              return;
            }
          }
        }
        DIE("Could not find the nb to remove");
      },
      vertices, cid, end, nb)
      .wait();
}

int64_t CtgGraph::purge_short_aln_edges() {
  int64_t num_short = 0;
  for (auto it = edges->begin(); it != edges->end();) {
    auto edge = make_shared<Edge>(it->second);
    if (edge->short_aln) {
      num_short++;
      remove_nb(edge->cids.cid1, edge->end1, edge->cids.cid2);
      remove_nb(edge->cids.cid2, edge->end2, edge->cids.cid1);
      it = edges->erase(it);
    } else {
      it++;
    }
  }
  return num_short;
}

void CtgGraph::add_pos_gap_read(const string &read_name) {
  upcxx::rpc(
      get_read_target_rank(read_name),
      [](reads_map_t &read_seqs, string read_name) {
        read_seqs->insert({read_name, ""});
      },
      read_seqs, read_name)
      .wait();
}

bool CtgGraph::update_read_seq(const string &read_name, const string &seq) {
  return upcxx::rpc(
             get_read_target_rank(read_name),
             [](reads_map_t &read_seqs, string read_name, string seq) {
               auto it = read_seqs->find(read_name);
               if (it == read_seqs->end()) return false;
               (*read_seqs)[read_name] = seq;
               return true;
             },
             read_seqs, read_name, seq)
      .wait();
}

string CtgGraph::get_read_seq(const string &read_name) {
  return upcxx::rpc(
             get_read_target_rank(read_name),
             [](reads_map_t &read_seqs, string read_name) -> string {
               const auto it = read_seqs->find(read_name);
               if (it == read_seqs->end()) return string("");
               return it->second;
             },
             read_seqs, read_name)
      .wait();
}

size_t CtgGraph::get_num_read_seqs(bool all) {
  if (!all)
    return upcxx::reduce_one(read_seqs->size(), upcxx::op_fast_add, 0).wait();
  else
    return upcxx::reduce_all(read_seqs->size(), upcxx::op_fast_add).wait();
}

int CtgGraph::get_other_end(shared_ptr<Vertex> v1, shared_ptr<Vertex> v2, shared_ptr<Edge> edge) {
  if (!edge) edge = get_edge_cached(v1->cid, v2->cid);
  // the cids in the edge are organized as (largest, smallest), so the edges are determined that way too
  return v1->cid >= v2->cid ? edge->end2 : edge->end1;
}

int CtgGraph::get_other_end_local(Vertex *v1, shared_ptr<Vertex> v2, shared_ptr<Edge> edge) {
  if (!edge) edge = get_edge(v1->cid, v2->cid);
  // the cids in the edge are organized as (largest, smallest), so the edges are determined that way too
  return v1->cid >= v2->cid ? edge->end2 : edge->end1;
}

void CtgGraph::clear_caches() {
  vertex_cache.clear();
  edge_cache.clear();
}

shared_ptr<Vertex> CtgGraph::get_vertex_cached(cid_t cid) {
  auto it = vertex_cache.find(cid);
  if (it != vertex_cache.end()) {
    upcxx::progress();
    return it->second;
  }
  auto v = get_vertex(cid);
  // load factor around .5
  if (vertex_cache.size() < CGRAPH_MAX_CACHE_SIZE / 2) vertex_cache[cid] = v;
  return v;
}

shared_ptr<Edge> CtgGraph::get_edge_cached(cid_t cid1, cid_t cid2) {
  CidPair cids = {.cid1 = cid1, .cid2 = cid2};
  if (cid1 < cid2) std::swap(cids.cid1, cids.cid2);
  auto it = edge_cache.find(cids);
  if (it != edge_cache.end()) {
    upcxx::progress();
    return it->second;
  }
  auto edge = get_edge(cids.cid1, cids.cid2);
  if (edge_cache.size() < CGRAPH_MAX_CACHE_SIZE / 2) edge_cache[cids] = edge;
  return edge;
}

void CtgGraph::print_stats(int min_ctg_print_len, string graph_fname) {
  BarrierTimer timer(__FILEFUNC__);
  auto get_avg_min_max = [](vector<int64_t> &vals) -> string {
    int64_t total = (!vals.size() ? 0 : std::accumulate(vals.begin(), vals.end(), 0));
    int64_t max_val = (!vals.size() ? 0 : *std::max_element(vals.begin(), vals.end()));
    int64_t min_val = (!vals.size() ? 0 : *std::min_element(vals.begin(), vals.end()));
    int64_t all_min_val = upcxx::reduce_one(min_val, upcxx::op_fast_min, 0).wait();
    int64_t all_max_val = upcxx::reduce_one(max_val, upcxx::op_fast_max, 0).wait();
    double all_total = upcxx::reduce_one(total, upcxx::op_fast_add, 0).wait();
    size_t all_nvals = upcxx::reduce_one(vals.size(), upcxx::op_fast_add, 0).wait();
    ostringstream os;
    os.precision(2);
    os << std::fixed;
    os << (all_total / all_nvals) << " [" << all_min_val << ", " << all_max_val << "]";
    return os.str();
  };

  vector<int64_t> depths;
  depths.reserve(get_local_num_vertices());
  vector<int64_t> clens;
  clens.reserve(get_local_num_vertices());
  for (auto v = get_first_local_vertex(); v != nullptr; v = get_next_local_vertex()) {
    if (v->clen >= min_ctg_print_len) {
      depths.push_back(round(v->depth));
      clens.push_back(v->clen);
    }
  }
  vector<int64_t> supports;
  supports.reserve(get_local_num_edges());
  vector<int64_t> aln_lens;
  aln_lens.reserve(get_local_num_edges());
  vector<int64_t> aln_scores;
  aln_scores.reserve(get_local_num_edges());
  vector<int64_t> gaps;
  gaps.reserve(get_local_num_edges());
  {
    ProgressBar progbar(get_local_num_edges(), "Compute graph stats");
    for (auto edge = get_first_local_edge(); edge != nullptr; edge = get_next_local_edge()) {
      aln_lens.push_back(edge->aln_len);
      aln_scores.push_back(edge->aln_score);
      auto clen1 = get_vertex_clen(edge->cids.cid1);
      auto clen2 = get_vertex_clen(edge->cids.cid2);
      if (clen1 >= min_ctg_print_len || clen2 >= min_ctg_print_len) {
        supports.push_back(edge->support);
        gaps.push_back(edge->gap);
      }
      progbar.update();
    }
    progbar.done();
  }

  auto num_vertices = get_num_vertices();
  auto num_edges = get_num_edges();
  SLOG_VERBOSE("Graph statistics:\n");
  SLOG_VERBOSE("    vertices:  ", num_vertices, "\n");
  SLOG_VERBOSE("    edges:     ", num_edges, "\n");
  SLOG_VERBOSE("    degree:    ", (double)num_edges / num_vertices, "\n");
  SLOG_VERBOSE("    aln_len:   ", get_avg_min_max(aln_lens), "\n");
  SLOG_VERBOSE("    aln_score: ", get_avg_min_max(aln_scores), "\n");
  SLOG_VERBOSE("  for contigs >= ", min_ctg_print_len, " length:\n");
  SLOG_VERBOSE("    depth:     ", get_avg_min_max(depths), "\n");
  SLOG_VERBOSE("    clen:      ", get_avg_min_max(clens), "\n");
  SLOG_VERBOSE("    support:   ", get_avg_min_max(supports), "\n");
  SLOG_VERBOSE("    gap:       ", get_avg_min_max(gaps), "\n");
}

#ifdef TNF_PATH_RESOLUTION
void CtgGraph::compute_edge_tnfs() {
  BarrierTimer timer(__FILEFUNC__);
  int num_short_ctgs = 0;
  double sum_tnf_prob = 0;
  double min_tnf_prob = 1;
  double max_tnf_prob = 0;
  for (auto edge = get_first_local_edge(); edge != nullptr; edge = get_next_local_edge()) {
    auto v1 = get_vertex(edge->cids.cid1);
    auto v2 = get_vertex(edge->cids.cid2);
    if (v1->clen < CGRAPH_MIN_TNF_CLEN || v2->clen < CGRAPH_MIN_TNF_CLEN) {
      edge->tnf_prob = 0;
      num_short_ctgs++;
      // DBG_BUILD("TNF: short ctgs in edge ", v1->cid, "->", v2->cid, " lengths ", v1->clen, " ", v2->clen, "\n");
      continue;
    }
    edge->tnf_prob = calc_tnf_dist(v1, v2);
    DBG_BUILD("TNF: edge ", v1->cid, "->", v2->cid, " ", edge->tnf_prob, "\n");
    max_tnf_prob = std::max(edge->tnf_prob, max_tnf_prob);
    min_tnf_prob = std::min(edge->tnf_prob, min_tnf_prob);
    sum_tnf_prob += edge->tnf_prob;
  }
  upcxx::barrier();
  int all_num_short_ctgs = upcxx::reduce_one(num_short_ctgs, upcxx::op_fast_add, 0).wait();
  double all_sum_tnf_prob = upcxx::reduce_one(sum_tnf_prob, upcxx::op_fast_add, 0).wait();
  double all_max_tnf_prob = upcxx::reduce_one(max_tnf_prob, upcxx::op_fast_max, 0).wait();
  double all_min_tnf_prob = upcxx::reduce_one(min_tnf_prob, upcxx::op_fast_min, 0).wait();
  auto num_edges = get_num_edges();
  SLOG_VERBOSE("Computed TNF probabilities for edges:\n");
  SLOG_VERBOSE("  average: ", all_sum_tnf_prob / num_edges, "\n");
  SLOG_VERBOSE("  min: ", all_min_tnf_prob, "\n");
  SLOG_VERBOSE("  max: ", all_max_tnf_prob, "\n");
  SLOG_VERBOSE("  short ctgs: ", perc_str(all_num_short_ctgs, num_edges), "\n");
}
#endif

void CtgGraph::print_gfa2(const string &gfa_fname, int min_ctg_print_len) {
  BarrierTimer timer(__FILEFUNC__);
  dist_ofstream of(gfa_fname + ".gfa");
  for (auto v = get_first_local_vertex(); v != nullptr; v = get_next_local_vertex()) {
    if (v->clen < min_ctg_print_len) continue;
    // don't include the sequence, and have a user tag 'kd' for kmer depth, and 'ad' for alignment depth
    of << "S\t" << to_string(v->cid) << "\t" << to_string(v->clen) << "\t*\tkd: " << to_string(v->depth)
       << " ad: " << to_string(v->aln_depth) << "\n";
  }
  for (auto edge = get_first_local_edge(); edge != nullptr; edge = get_next_local_edge()) {
    int clen1 = get_vertex_clen(edge->cids.cid1);
    if (clen1 < min_ctg_print_len) continue;
    int clen2 = get_vertex_clen(edge->cids.cid2);
    if (clen2 < min_ctg_print_len) continue;
    auto cid_str = to_string(edge->cids.cid1) + (edge->end1 == 5 ? "+" : "-") + "\t" + to_string(edge->cids.cid2) +
                   (edge->end2 == 3 ? "+" : "-");
    if (edge->gap >= 0) {
      // this is a gap, not an edge, according to the GFA terminology
      of << "G\t*\t" << cid_str << "\t" << to_string(edge->gap) << "\t*";
    } else {
      // this is a negative gap - set as an alignment overlap in the GFA output
      // we don't have the actual alignments, but we know that for this to be valid, the alignment must be almost
      // perfect over the tail and front of the two contigs, so we can give the positions based on the gap size
      of << "E\t*\t" << cid_str << "\t";
      int overlap = -edge->gap;
      // positions are specified *before* revcomp
      int begin_pos1 = (edge->end1 == 5 ? clen1 - overlap : 0);
      int end_pos1 = (edge->end1 == 5 ? clen1 : overlap);
      int begin_pos2 = (edge->end1 == 3 ? clen2 - overlap : 0);
      int end_pos2 = (edge->end1 == 3 ? clen2 : overlap);
      of << to_string(begin_pos1) << "\t" << to_string(end_pos1);
      if (end_pos1 == clen1) of << "$";
      if (end_pos2 == clen2) of << "$";
    }
    // add MHM specific tags
    of << "\tsp: " << to_string(edge->support) << "\ttp: " << edge_type_str(edge->edge_type) << "\n";
  }
  of.close();  // sync and prints stats
}
