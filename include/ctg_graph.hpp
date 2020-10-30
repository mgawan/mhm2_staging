#pragma once

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

#include <memory>
#include <ostream>
#include <string>
#include <vector>

using std::make_shared;
using std::ostream;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::vector;

#include "utils.hpp"

//#define DBG_BUILD DBG
#define DBG_BUILD(...)

//#define DBG_WALK DBG
#define DBG_WALK(...)

//#define DBG_WALK_CONT DBG_CONT
#define DBG_WALK_CONT(...)

//#define DBG_SPANS DBG
#define DBG_SPANS(...)

using cid_t = int64_t;

enum class Dirn { FORWARD, BACKWARD };
enum class Orient { NORMAL, REVCOMP };

inline string dirn_str(Dirn dirn) { return (dirn == Dirn::FORWARD ? "forward" : "backward"); }

inline string orient_str(Orient orient) { return (orient == Orient::NORMAL ? "+" : "-"); }

inline Orient flip_orient(Orient orient) { return (orient == Orient::NORMAL ? Orient::REVCOMP : Orient::NORMAL); }

struct GapRead {
  string read_name;
  // used for resolving positive splints
  int gap_start;
  // used for resolving positive spans
  // int rstart, rstop;
  char orient;
  cid_t cid;

  UPCXX_SERIALIZED_FIELDS(read_name, gap_start, orient, cid);

  GapRead() {}

  // GapRead(const string &read_name, int gap_start, int rstart, int rstop, int orient, cid_t cid) {
  GapRead(const string &read_name, int gap_start, int orient, cid_t cid)
      : read_name(read_name)
      , gap_start(gap_start)
      , orient(orient)
      , cid(cid) {}

  bool operator==(const GapRead &other) const;
};

struct CidPair {
  cid_t cid1, cid2;

  bool operator==(const CidPair &other) const;

  bool operator!=(const CidPair &other) const;

  friend ostream &operator<<(ostream &os, const CidPair &cids); /*{
    os << "(" << cids.cid1 << ", " << cids.cid2 << ")";
    return os;
  }*/
};

inline const CidPair NULL_CIDS = CidPair{-1, -1};

namespace std {
template <>
struct hash<CidPair> {
  std::size_t operator()(const CidPair &cids) const {
    return (std::hash<cid_t>()(cids.cid1) ^ (std::hash<cid_t>()(cids.cid2) << 1));
  }
};
}  // namespace std

enum class EdgeType { SPLINT, SPAN };

string edge_type_str(EdgeType edge_type);

struct Edge {
  // the cids of the vertices connected by this edge. Note that the largest number cid is always first
  CidPair cids;
  // the ends correspond to the cids above.
  int end1, end2;
  int gap;
  int support;
  int aln_len;
  // the best aln score that generated this edge
  int aln_score;
  EdgeType edge_type;
  // the sequence filling a positive gap - if the gap is non-positive, this is empty
  string seq;
  // these flags are set during graph construction to keep track of errors in edge construction
  bool mismatch_error, conflict_error, excess_error, short_aln;
  // contains information of reads that map to a positive gap - used for filling the gap
  vector<GapRead> gap_reads;
#ifdef TNF_PATH_RESOLUTION
  // TNF probability for this link
  double tnf_prob;
#endif

  UPCXX_SERIALIZED_FIELDS(cids, end1, end2, gap, support, aln_len, aln_score, edge_type, seq, mismatch_error, conflict_error,
                          excess_error, short_aln, gap_reads
#ifdef TNF_PATH_RESOLUTION
                          ,
                          tnf_prob
#endif
  );
};

struct Vertex {
  cid_t cid;
  int clen;
  double depth;
  // track depth from alignments
  double aln_depth;
  // set to true if visited in previous round
  bool visited;
  upcxx::global_ptr<char> seq_gptr;
  // the neighbors at the different ends
  vector<cid_t> end5;
  vector<cid_t> end3;
  // the merged series of nbs
  // FIXME: when using spans make sure these are valid
  vector<vector<cid_t>> end5_merged;
  vector<vector<cid_t>> end3_merged;
#ifdef TNF_PATH_RESOLUTION
  // the tnf distribution for this vertex
  tnf_t tnf;
#endif
  // book-keeping fields for resolving walk conflicts between ranks -
  // choose the walk with the longest scaffold, and if there is a tie, choose the highest rank
  int walk_score;
  int walk_rank;
  int walk_i;
  UPCXX_SERIALIZED_FIELDS(cid, clen, depth, aln_depth, visited, seq_gptr, end5, end3, end5_merged, end3_merged,
#ifdef TNF_PATH_RESOLUTION
                          tnf,
#endif
                          walk_score, walk_rank, walk_i);
};

class CtgGraph {
 private:
  using vertex_map_t = upcxx::dist_object<HASH_TABLE<cid_t, Vertex>>;
  using edge_map_t = upcxx::dist_object<HASH_TABLE<CidPair, Edge>>;
  using reads_map_t = upcxx::dist_object<HASH_TABLE<string, string>>;
  vertex_map_t vertices;
  edge_map_t edges;
  reads_map_t read_seqs;
  HASH_TABLE<cid_t, shared_ptr<Vertex>> vertex_cache;
  HASH_TABLE<CidPair, shared_ptr<Edge>> edge_cache;

  struct VertexDepthInfo {
    cid_t cid;
    double depth;
    int clen;
  };

  HASH_TABLE<cid_t, Vertex>::iterator vertex_iter;
  HASH_TABLE<CidPair, Edge>::iterator edge_iter;

  size_t get_vertex_target_rank(cid_t cid);

  size_t get_edge_target_rank(CidPair &cids);

  size_t get_read_target_rank(const string &r);

#ifdef TNF_PATH_RESOLUTION
  double calc_tnf_dist(shared_ptr<Vertex> v1, shared_ptr<Vertex> v2);
#endif

 public:
  int max_read_len;

  CtgGraph();

  void clear();

  ~CtgGraph();

  int64_t get_num_vertices(bool all = false);

  int64_t get_local_num_vertices(void);

  int64_t get_num_edges(bool all = false);

  int64_t get_local_num_edges(void);

  shared_ptr<Vertex> get_vertex(cid_t cid);

  int get_vertex_clen(cid_t cid);

  shared_ptr<Vertex> get_local_vertex(cid_t cid);

  void set_vertex_visited(cid_t cid);

  void update_vertex_walk(cid_t cid, int walk_score, int walk_i);

  Vertex *get_first_local_vertex();

  Vertex *get_next_local_vertex();

  void add_vertex(Vertex &v, const string &seq);

  void add_vertex_nb(cid_t cid, cid_t nb, char end);

  string get_vertex_seq(upcxx::global_ptr<char> seq_gptr, int64_t seq_len);

  shared_ptr<Edge> get_edge(cid_t cid1, cid_t cid2);

  Edge *get_first_local_edge();

  Edge *get_next_local_edge();

  void add_or_update_edge(Edge &edge);

  void purge_error_edges(int64_t *mismatched, int64_t *conflicts, int64_t *empty_spans);

  int64_t purge_excess_edges();

  void remove_nb(cid_t cid, int end, cid_t nb);

  int64_t purge_short_aln_edges();

  void add_pos_gap_read(const string &read_name);

  bool update_read_seq(const string &read_name, const string &seq);

  string get_read_seq(const string &read_name);

  size_t get_num_read_seqs(bool all = false);

  int get_other_end(shared_ptr<Vertex> v1, shared_ptr<Vertex> v2, shared_ptr<Edge> edge = nullptr);

  int get_other_end_local(Vertex *v1, shared_ptr<Vertex> v2, shared_ptr<Edge> edge = nullptr);

  void clear_caches();

  shared_ptr<Vertex> get_vertex_cached(cid_t cid);

  shared_ptr<Edge> get_edge_cached(cid_t cid1, cid_t cid2);

  void print_stats(int min_ctg_print_len, string graph_fname = "");
#ifdef TNF_PATH_RESOLUTION
  void compute_edge_tnfs();
#endif

  void print_gfa2(const string &gfa_fname, int min_ctg_print_len);
};
