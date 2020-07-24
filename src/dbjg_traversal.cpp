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
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"

#include "utils.hpp"
#include "kmer_dht.hpp"
#include "contigs.hpp"

#define DBG_TRAVERSE DBG
//#define DBG_TRAVERSE(...)

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

extern ofstream _dbgstream;

enum class Dirn { LEFT, RIGHT, NONE };
#define DIRN_STR(d) ((d) == Dirn::LEFT ? "left" : (d) == Dirn::RIGHT ? "right" : "none")

enum class WalkStatus { RUNNING = '-', DEADEND = 'X', FORK = 'F', CONFLICT = 'O', REPEAT = 'R', VISITED = 'V'};

struct FragElem {
  global_ptr<FragElem> left_gptr, right_gptr;
  bool left_is_rc, right_is_rc;
  global_ptr<char> frag_seq;
  int frag_len;
  int64_t sum_depths;
  bool visited;

  FragElem() : left_gptr(nullptr), right_gptr(nullptr), left_is_rc(false), right_is_rc(false), frag_seq(nullptr),
               frag_len(0), sum_depths(0), visited(false) {}
};

struct StepInfo {
  WalkStatus walk_status;
  ext_count_t count;
  char left, right;
  global_ptr<FragElem> visited_frag_elem_gptr;
};

struct WalkTermStats {
  int64_t num_deadends, num_forks, num_conflicts, num_repeats, num_visited;

  void update(WalkStatus walk_status) {
    switch (walk_status) {
      case WalkStatus::DEADEND: num_deadends++; break;
      case WalkStatus::FORK: num_forks++; break;
      case WalkStatus::CONFLICT: num_conflicts++; break;
      case WalkStatus::REPEAT: num_repeats++; break;
      case WalkStatus::VISITED: num_visited++; break;
      default: DIE("Should never get here\n"); break;
    }
  }

  void print() {
    auto all_num_deadends = reduce_one(num_deadends, op_fast_add, 0).wait();
    auto all_num_forks = reduce_one(num_forks, op_fast_add, 0).wait();
    auto all_num_conflicts = reduce_one(num_conflicts, op_fast_add, 0).wait();
    auto all_num_repeats = reduce_one(num_repeats, op_fast_add, 0).wait();
    auto all_num_visited = reduce_one(num_visited, op_fast_add, 0).wait();
    auto tot_ends = all_num_forks + all_num_deadends + all_num_conflicts + all_num_repeats + all_num_visited;
    SLOG_VERBOSE("Walk statistics:\n");
    SLOG_VERBOSE("  deadends:  ", perc_str(all_num_deadends, tot_ends), "\n");
    SLOG_VERBOSE("  forks:     ", perc_str(all_num_forks, tot_ends), "\n");
    SLOG_VERBOSE("  conflicts: ", perc_str(all_num_conflicts, tot_ends), "\n");
    SLOG_VERBOSE("  repeats:   ", perc_str(all_num_repeats, tot_ends), "\n");
    SLOG_VERBOSE("  visited:   ", perc_str(all_num_visited, tot_ends), "\n");
  }
};

static string gptr_str(global_ptr<FragElem> gptr) {
  if (!gptr) return string(10, '0');
  ostringstream oss;
  oss << setw(11);
  oss << gptr.raw_ptr_;
  string s = oss.str();
  s.erase(0, s.length() - 6);
  return to_string(gptr.where()) + ":" + s;
}

#ifdef DEBUG
template<int MAX_K>
static bool check_kmers(const string &seq, dist_object<KmerDHT<MAX_K>> &kmer_dht, int kmer_len) {
  vector<Kmer<MAX_K>> kmers;
  Kmer<MAX_K>::get_kmers(kmer_len, seq, kmers);
  for (auto kmer : kmers) {
    if (!kmer_dht->kmer_exists(kmer)) return false;
  }
  return true;
}
#endif

template<int MAX_K>
static future<StepInfo> get_next_step(dist_object<KmerDHT<MAX_K>> &kmer_dht, Kmer<MAX_K> kmer, Dirn dirn, char prev_ext,
                                      global_ptr<FragElem> frag_elem_gptr, bool revisit_allowed) {
  auto kmer_rc = kmer.revcomp();
  bool is_rc = false;
  if (kmer_rc < kmer) {
    kmer = kmer_rc;
    is_rc = true;
  }
  return rpc(kmer_dht->get_kmer_target_rank(kmer),
             [](dist_object<KmerDHT<MAX_K>> &kmer_dht, Kmer<MAX_K> kmer, Dirn dirn, char prev_ext, bool revisit_allowed,
                bool is_rc, global_ptr<FragElem> frag_elem_gptr) -> StepInfo {
               KmerCounts *kmer_counts = kmer_dht->get_local_kmer_counts(kmer);
               // this kmer doesn't exist, abort
               if (!kmer_counts) return {.walk_status = WalkStatus::DEADEND,};
               char left = kmer_counts->left;
               char right = kmer_counts->right;
               if (left == 'X' || right == 'X') return {.walk_status = WalkStatus::DEADEND};
               if (left == 'F' || right == 'F') return {.walk_status = WalkStatus::FORK};
               if (is_rc) {
                 left = comp_nucleotide(left);
                 right = comp_nucleotide(right);
                 swap(left, right);
               }
               // check for conflicts
               if (prev_ext && ((dirn == Dirn::LEFT && prev_ext != right) || (dirn == Dirn::RIGHT && prev_ext != left)))
                 return {.walk_status = WalkStatus::CONFLICT};
               // if visited by another rank first
               if (kmer_counts->uutig_frag && kmer_counts->uutig_frag != frag_elem_gptr)
                 return {.walk_status = WalkStatus::VISITED, .count = 0, .left = 0, .right = 0,
                         .visited_frag_elem_gptr = kmer_counts->uutig_frag};
               // a repeat, abort (but allowed if traversing right after adding start kmer previously)
               if (kmer_counts->uutig_frag == frag_elem_gptr && !revisit_allowed)
                 return {.walk_status = WalkStatus::REPEAT};
               // mark as visited
               kmer_counts->uutig_frag = frag_elem_gptr;
               return {.walk_status = WalkStatus::RUNNING, .count = kmer_counts->count, .left = left, .right = right};
             }, kmer_dht, kmer, dirn, prev_ext, revisit_allowed, is_rc, frag_elem_gptr);
}

template<int MAX_K>
static global_ptr<FragElem> traverse_dirn(dist_object<KmerDHT<MAX_K>> &kmer_dht, Kmer<MAX_K> kmer,
                                          global_ptr<FragElem> frag_elem_gptr, Dirn dirn, string &uutig, int64_t &sum_depths,
                                          WalkTermStats &walk_term_stats) {
  auto kmer_str = kmer.to_string();
  char prev_ext = 0;
  char next_ext = (dirn == Dirn::LEFT ? kmer_str.front() : kmer_str.back());
  bool revisit_allowed = (dirn == Dirn::LEFT ? false : true);
  if (dirn == Dirn::RIGHT) uutig += substr_view(kmer_str, 1, kmer_str.length() - 2);
  while (true) {
    //progress();
    auto step_info = get_next_step(kmer_dht, kmer, dirn, prev_ext, frag_elem_gptr, revisit_allowed).wait();
    revisit_allowed = false;
    if (step_info.walk_status != WalkStatus::RUNNING) {
      walk_term_stats.update(step_info.walk_status);
      // reverse it because we were walking backwards
      if (dirn == Dirn::LEFT) reverse(uutig.begin(), uutig.end());
      return step_info.visited_frag_elem_gptr;
    }
    sum_depths += step_info.count;
    // add the last nucleotide to the uutig
    uutig += next_ext;
    // now attempt to walk to next kmer
    kmer_str = kmer.to_string();
    // get next extension
    if (dirn == Dirn::LEFT) {
      next_ext = step_info.left;
      prev_ext = kmer_str.back();
      kmer = kmer.backward_base(next_ext);
    } else {
      next_ext = step_info.right;
      prev_ext = kmer_str.front();
      kmer = kmer.forward_base(next_ext);
    }
    if (next_ext == 'X' || next_ext == 'F') DIE("Found X or F");
  }
}

template<int MAX_K>
static void construct_frags(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht,
                            vector<global_ptr<FragElem>> &frag_elems) {
  BarrierTimer timer(__FILEFUNC__);
  // allocate space for biggest possible uutig in global storage
  WalkTermStats walk_term_stats = {0};
  int64_t num_walks = 0;
  barrier();
  ProgressBar progbar(kmer_dht->get_local_num_kmers(), "DeBruijn graph traversal to construct uutig fragments");
  for (auto it = kmer_dht->local_kmers_begin(); it != kmer_dht->local_kmers_end(); it++) {
    progress();
    progbar.update();
    auto kmer = it->first;
    auto kmer_counts = &it->second;
    // don't start any new walk if this kmer has already been visited
    if (kmer_counts->uutig_frag) continue;
    // don't start walks on kmers without extensions on both sides
    if (kmer_counts->left == 'X' || kmer_counts->left == 'F' || kmer_counts->right == 'X' || kmer_counts->right == 'F')
      continue;
    string uutig;
    int64_t sum_depths = 0;
    global_ptr<FragElem> frag_elem_gptr = new_<FragElem>();
    auto left_gptr = traverse_dirn(kmer_dht, kmer, frag_elem_gptr, Dirn::LEFT, uutig, sum_depths, walk_term_stats);
    auto right_gptr = traverse_dirn(kmer_dht, kmer, frag_elem_gptr, Dirn::RIGHT, uutig, sum_depths, walk_term_stats);
    FragElem *frag_elem = frag_elem_gptr.local();
    frag_elem->frag_seq = new_array<char>(uutig.length() + 1);
    strcpy(frag_elem->frag_seq.local(), uutig.c_str());
    frag_elem->frag_seq.local()[uutig.length()] = 0;
    frag_elem->frag_len = uutig.length();
    frag_elem->sum_depths = sum_depths;
    frag_elem->left_gptr = left_gptr;
    frag_elem->right_gptr = right_gptr;
    frag_elems.push_back(frag_elem_gptr);
    num_walks++;
  }
  progbar.done();
  barrier();
  walk_term_stats.print();
}

static int64_t print_link_stats(int64_t num_links, int64_t num_overlaps, int64_t num_overlaps_rc, const string &dirn_str) {
  auto all_num_links = reduce_one(num_links, op_fast_add, 0).wait();
  auto all_num_overlaps = reduce_one(num_overlaps, op_fast_add, 0).wait();
  auto all_num_overlaps_rc = reduce_one(num_overlaps_rc, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", all_num_links, " ", dirn_str, " links with ",
               perc_str(all_num_overlaps, all_num_links), " overlaps and ",
               perc_str(all_num_overlaps_rc, all_num_links), " revcomped overlaps\n");
  return all_num_links;
}

static bool is_overlap(const string &left_seq, const string &right_seq, int overlap_len) {
  return (left_seq.compare(left_seq.length() - overlap_len, overlap_len, right_seq, 0, overlap_len) == 0);
}

static string get_frag_seq(FragElem &frag_elem) {
#ifdef ENABLE_GPUS
  // FIXME: this hack is to avoid the crazy slowdown we see with IBV initialization
  global_ptr<char> gbuf = new_array<char>(frag_elem.frag_len + 1);
  char *buf = gbuf.local();
#else
  char *buf = new char[frag_elem.frag_len + 1];
#endif
  rget(frag_elem.frag_seq, buf, frag_elem.frag_len + 1).wait();
  string frag_seq(buf);
  assert(frag_seq.length() == frag_elem.frag_len);
#ifdef ENABLE_GPUS
  delete_array(gbuf);
#else
  delete[] buf;
#endif
  return frag_seq;
}

static void set_link_status(Dirn dirn, global_ptr<FragElem> &nb_gptr, bool &is_rc, string &uutig, int kmer_len,
                            int64_t &num_overlaps, int64_t &num_overlaps_rc, int64_t &num_non_recip) {
  if (nb_gptr) {
    FragElem nb_frag_elem = rget(nb_gptr).wait();
    string nb_frag_seq = get_frag_seq(nb_frag_elem);
    string *s1 = (dirn == Dirn::LEFT ? &nb_frag_seq : &uutig);
    string *s2 = (dirn == Dirn::LEFT ? &uutig : &nb_frag_seq);
    if (is_overlap(*s1, *s2, kmer_len - 1)) {
      if ((dirn == Dirn::LEFT ? nb_frag_elem.right_gptr : nb_frag_elem.left_gptr) == nb_gptr) {
        num_non_recip++;
        nb_gptr = nullptr;
        return;
      }
      num_overlaps++;
      return;
    }
    auto nb_frag_seq_rc = revcomp(nb_frag_seq);
    s1 = (dirn == Dirn::LEFT ? &nb_frag_seq_rc : &uutig);
    s2 = (dirn == Dirn::LEFT ? &uutig : &nb_frag_seq_rc);
    if (is_overlap(*s1, *s2, kmer_len - 1)) {
      if ((dirn == Dirn::LEFT ? nb_frag_elem.left_gptr : nb_frag_elem.right_gptr) == nb_gptr) {
        num_non_recip++;
        nb_gptr = nullptr;
        return;
      }
      num_overlaps_rc++;
      is_rc = true;
      return;
    }
    DBG_TRAVERSE("No ", DIRN_STR(dirn), " overlap:\n", uutig, "\n", nb_frag_seq, "\n", nb_frag_seq_rc, "\n");
  }
}

template<int MAX_K>
static void clean_frag_links(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht,
                             vector<global_ptr<FragElem>> &frag_elems) {
  BarrierTimer timer(__FILEFUNC__);
  // put all the uutigs found by this rank into my_uutigs
  int64_t num_equal_links = 0, num_non_recip = 0, num_short = 0,
          num_left_links = 0, num_left_overlaps = 0, num_left_overlaps_rc = 0,
          num_right_links = 0, num_right_overlaps = 0, num_right_overlaps_rc = 0;
  ProgressBar progbar(frag_elems.size(), "Cleaning fragment links");
  for (auto frag_elem_gptr : frag_elems) {
    progbar.update();
    FragElem *frag_elem = frag_elem_gptr.local();
    if (frag_elem->frag_len < kmer_len) {
      num_short++;
      continue;
    }
    if (frag_elem->left_gptr) num_left_links++;
    if (frag_elem->right_gptr) num_right_links++;
    string uutig(frag_elem->frag_seq.local());
    if (frag_elem->left_gptr && frag_elem->left_gptr == frag_elem->right_gptr) {
      num_equal_links++;
      frag_elem->left_gptr = nullptr;
      frag_elem->right_gptr = nullptr;
      continue;
    }
    set_link_status(Dirn::LEFT, frag_elem->left_gptr, frag_elem->left_is_rc, uutig, kmer_len,
                    num_left_overlaps, num_left_overlaps_rc, num_non_recip);
    set_link_status(Dirn::RIGHT, frag_elem->right_gptr, frag_elem->right_is_rc, uutig, kmer_len,
                    num_right_overlaps, num_right_overlaps_rc, num_non_recip);
  }
  progbar.done();
  barrier();
  auto all_num_frags = reduce_one(frag_elems.size(), op_fast_add, 0).wait();
  auto all_num_short = reduce_one(num_short, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", all_num_frags, " uutig fragments of which ", perc_str(all_num_short, all_num_frags), " are short\n");
  auto all_num_left = print_link_stats(num_left_links, num_left_overlaps, num_left_overlaps_rc, "left");
  auto all_num_right = print_link_stats(num_right_links, num_right_overlaps, num_right_overlaps_rc, "right");
  auto all_num_equal_links = reduce_one(num_equal_links, op_fast_add, 0).wait();
  auto all_num_non_recip = reduce_one(num_non_recip, op_fast_add, 0).wait();
  SLOG_VERBOSE("There were ", perc_str(all_num_equal_links, all_num_left + all_num_right), " equal left and right links\n");
  SLOG_VERBOSE("There were ", perc_str(all_num_non_recip, all_num_left + all_num_right), " non-reciprocating links\n");
}

static global_ptr<FragElem> get_other_side_gptr(const FragElem &frag_elem, global_ptr<FragElem> frag_elem_gptr) {
  if (frag_elem.left_gptr == frag_elem_gptr) return frag_elem.right_gptr;
  return frag_elem.left_gptr;
}

static bool walk_frags_dirn(unsigned kmer_len, global_ptr<FragElem> frag_elem_gptr, global_ptr<FragElem> next_gptr,
                            string &uutig, int64_t &depths, int64_t &walk_steps, int64_t &num_repeats,
                            vector<FragElem*> &my_frag_elems_visited) {
  if (!next_gptr) return true;
  global_ptr<FragElem> prev_gptr = frag_elem_gptr;
  FragElem prev_frag_elem = *frag_elem_gptr.local();
  // for checking that we haven't got a bug - frags should never be revisited in a walk
  HASH_TABLE<global_ptr<FragElem>, bool> visited;
  visited[frag_elem_gptr] = true;
#ifdef DEBUG
  string padding;
  DBG_TRAVERSE(uutig, "\n");
#endif
  Dirn dirn = Dirn::NONE;
  while (next_gptr) {
    DBG_TRAVERSE(padding, gptr_str(get_other_side_gptr(prev_frag_elem, next_gptr)), " <-- ", gptr_str(prev_gptr), " ==> ",
                 gptr_str(next_gptr), "\n");
    if (next_gptr.where() > rank_me()) {
      DBG_TRAVERSE(padding, "DROP: owner ", next_gptr.where(), " > ", rank_me(), "\n");
      return false;
    }
    if (visited.find(next_gptr) != visited.end()) {
      DBG_TRAVERSE(padding, "REPEAT: ", gptr_str(next_gptr), "\n");
      num_repeats++;
      return true;
    }
    visited[next_gptr] = true;
    FragElem next_frag_elem = rget(next_gptr).wait();
    if (next_gptr.where() == rank_me()) {
      if (next_frag_elem.visited) DIE("gptr ", next_gptr, " should not be already visited");
      my_frag_elems_visited.push_back(next_gptr.local());
    }
    string next_frag_seq = get_frag_seq(next_frag_elem);
    string next_frag_seq_rc = revcomp(next_frag_seq);
    if (dirn == Dirn::NONE) {
      if (is_overlap(uutig, next_frag_seq, kmer_len - 1)) dirn = Dirn::RIGHT;
      else if (is_overlap(next_frag_seq, uutig, kmer_len - 1)) dirn = Dirn::LEFT;
      if (dirn == Dirn::NONE) {
        if (is_overlap(uutig, next_frag_seq_rc, kmer_len - 1)) dirn = Dirn::RIGHT;
        else if (is_overlap(next_frag_seq_rc, uutig, kmer_len - 1)) dirn = Dirn::LEFT;
        else DIE("No overlap");
      }
      DBG_TRAVERSE(padding, "Direction is set to ", DIRN_STR(dirn), "\n");
    }
    if (dirn == Dirn::LEFT) {
      int slen = next_frag_seq.length() - kmer_len + 1;
      DBG_TRAVERSE(string(slen, ' '), uutig, "\n");
      if (is_overlap(next_frag_seq, uutig, kmer_len - 1)) uutig.insert(0, substr_view(next_frag_seq, 0, slen));
      else if (is_overlap(next_frag_seq_rc, uutig, kmer_len - 1)) uutig.insert(0, substr_view(next_frag_seq_rc, 0, slen));
      else DIE("No valid overlap in dirn ", DIRN_STR(dirn));
    } else {
      if (is_overlap(uutig, next_frag_seq, kmer_len - 1)) uutig += substr_view(next_frag_seq, kmer_len - 1);
      else if (is_overlap(uutig, next_frag_seq_rc, kmer_len - 1)) uutig += substr_view(next_frag_seq_rc, kmer_len - 1);
      else DIE("No valid overlap in dirn ", DIRN_STR(dirn));
    }
    DBG_TRAVERSE(uutig, "\n");
    depths += (next_frag_elem.sum_depths * (1.0 - (kmer_len - 1) / next_frag_elem.frag_len));
    auto other_side_gptr = get_other_side_gptr(next_frag_elem, prev_gptr);
    prev_frag_elem = next_frag_elem;
    prev_gptr = next_gptr;
    next_gptr = other_side_gptr;
#ifdef DEBUG
    padding += string(4, ' ');
#endif
    walk_steps++;
  }
  DBG_TRAVERSE(padding, "DEADEND\n");
  return true;
}

template<int MAX_K>
static void connect_frags(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, vector<global_ptr<FragElem>> &frag_elems,
                          Contigs &my_uutigs) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_steps = 0, max_steps = 0, num_drops = 0, num_prev_visited = 0, num_repeats = 0;
  ProgressBar progbar(frag_elems.size(), "Connecting fragments");
  for (auto frag_elem_gptr : frag_elems) {
    progbar.update();
    FragElem *frag_elem = frag_elem_gptr.local();
    if (frag_elem->frag_len < kmer_len) continue;
    if (frag_elem->visited) {
      num_prev_visited++;
      continue;
    }
    vector<FragElem*> my_frag_elems_visited;
    string uutig(frag_elem->frag_seq.local());
    int64_t depths = frag_elem->sum_depths;
    int64_t walk_steps = 1;
    bool walk_ok = walk_frags_dirn(kmer_len, frag_elem_gptr, frag_elem->left_gptr, uutig, depths, walk_steps, num_repeats,
                                   my_frag_elems_visited);
    if (walk_ok)
      walk_ok = walk_frags_dirn(kmer_len, frag_elem_gptr, frag_elem->right_gptr, uutig, depths, walk_steps, num_repeats,
                        my_frag_elems_visited);
    if (walk_ok) {
      num_steps += walk_steps;
      max_steps = max(walk_steps, max_steps);
      Contig contig = {0, uutig, (double)depths / (uutig.length() - kmer_len + 2)};
      my_uutigs.add_contig(contig);
      // the walk is successful, so set the visited for all the local elems
      for (auto &elem : my_frag_elems_visited) elem->visited = true;
    } else {
      num_drops++;
    }
  }
  progbar.done();
  auto all_num_steps = reduce_one(num_steps, op_fast_add, 0).wait();
  auto all_max_steps = reduce_one(max_steps, op_fast_max, 0).wait();
  auto all_num_drops = reduce_one(num_drops, op_fast_add, 0).wait();
  auto all_num_repeats = reduce_one(num_repeats, op_fast_add, 0).wait();
  auto all_num_uutigs = reduce_one(my_uutigs.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Constructed ", all_num_uutigs, " uutigs with ", (double)all_num_steps / all_num_uutigs,
               " avg path length (max ", all_max_steps, "), dropped ", perc_str(all_num_drops, all_num_uutigs), " paths\n");
  auto all_num_prev_visited = reduce_one(num_prev_visited, op_fast_add, 0).wait();
  auto all_num_frags = reduce_one(frag_elems.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Skipped ", perc_str(all_num_prev_visited, all_num_frags), " already visited fragments, and found ",
               perc_str(all_num_repeats, all_num_frags), " repeats\n");
  barrier();
  for (auto frag_elem_gptr : frag_elems) {
    FragElem *frag_elem = frag_elem_gptr.local();
    delete_array(frag_elem->frag_seq);
    delete_(frag_elem_gptr);
  }
}

template<int MAX_K>
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, Contigs &my_uutigs) {
  BarrierTimer timer(__FILEFUNC__);
  {
    // scope for frag_elems
    vector<global_ptr<FragElem>> frag_elems;
    construct_frags(kmer_len, kmer_dht, frag_elems);
    clean_frag_links(kmer_len, kmer_dht, frag_elems);
    // put all the uutigs found by this rank into my_uutigs
    my_uutigs.clear();
    connect_frags(kmer_len, kmer_dht, frag_elems, my_uutigs);
  }
  // now get unique ids for the uutigs
  atomic_domain<size_t> ad({atomic_op::fetch_add, atomic_op::load});
  global_ptr<size_t> counter = nullptr;
  if (!rank_me()) counter = new_<size_t>(0);
  counter = broadcast(counter, 0).wait();
  size_t my_counter = ad.fetch_add(counter, my_uutigs.size(), memory_order_relaxed).wait();
  // wait until all ranks have updated the global counter
  barrier();
  if (!rank_me()) upcxx::delete_(counter);
  ad.destroy();
  // set the unique ids
  int64_t cid = my_counter;
  for (auto it = my_uutigs.begin(); it != my_uutigs.end(); ++it) {
    it->id = cid;
    cid++;
  }
  barrier();
#ifdef DEBUG
  ProgressBar progbar(my_uutigs.size(), "Checking kmers in uutigs");
  for (auto uutig : my_uutigs) {
    progbar.update();
    if (!check_kmers(uutig.seq, kmer_dht, kmer_len)) DIE("kmer not found in uutig");
  }
  progbar.done();
#endif
}

#define TDG_K(KMER_LEN) \
    template \
    void traverse_debruijn_graph<KMER_LEN>(unsigned kmer_len, dist_object<KmerDHT<KMER_LEN>> &kmer_dht, Contigs &my_uutigs)

TDG_K(32);
#if MAX_BUILD_KMER >= 64
TDG_K(64);
#endif
#if MAX_BUILD_KMER >= 96
TDG_K(96);
#endif
#if MAX_BUILD_KMER >= 128
TDG_K(128);
#endif
#if MAX_BUILD_KMER >= 160
TDG_K(160);
#endif

#undef TDG_K
