// kcount - kmer counting
// Steven Hofmeyr, LBNL, June 2019

#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "progressbar.hpp"
#include "kmer_dht.hpp"
#include "contigs.hpp"

#define DBG_TRAVERSE DBG
//#define DBG_TRAVERSE(...)

using namespace std;
using namespace upcxx;

extern ofstream _dbgstream;

// A problem with the DBJG traversal is that we can extreme load imbalance at scale in the later rounds of contigging, for 
// example, an average that is 0.01 compared to the maximum. This is a consequence of some paths being very long, and ideally
// the work building those paths should be distributed across multiple ranks.
//
// The way this is done is by having each rank create a doubly-linked list of uutig fragments for each dbjg traversal. Each 
// element is contains a global ptr to the left and right elements, and a global ptr to fragment sequence buffer (which is a 
// char array of fixed maximum size). 
//
// As the traversal proceeds, a rank builds a uutig fragment until it reaches a termination in the graph, or encounters a kmer
// already visited, or fills up the fragment sequence buffer. Note that the traversal only needs to proceed in a single
// direction (left to right), because all kmers to the left of the starting kmer will be visited at some point by some rank.
// If the fragment buffer is full, it creates a new element and links that to the current one (setting the current element's
// right gptr to the new element, and the new element's left gptr to the current element) and continues the traverse. If the 
// kmer is already visited, it links the current element to that fragment (which is recorded with a global ptr in the kmer 
// hash table), setting the right gptr of the current element, and the left gptr of the already visited element. If the 
// traversal terminates on an 'F' or 'X', then there are no more added fragments. At the end of the dbjg traversal, there 
// will be a set of doubly-linked lists of fragments representing walks through the graph.
// 
// In the next stage, each rank traverses each list with a local element, and rgets each next element in turn to construct the
// full uutig. If the gptr for the next element is owned by a higher rank, then the link up is abandoned and the rank moves on
// to the next local element. Otherwise the traverse continues. In this way, no duplicate uutigs are built, and there is no 
// need to store additional information with the elements about visited vertices. Repeats are taken care of by keeping a local
// hash table of all fragments currently in the path. In this stage, the traverse will have to follow each path in both left 
// and right directions.


struct FragElem {
  global_ptr<FragElem> left_gptr, right_gptr;
  char frag_seq[MAX_FRAG_LEN + 1]; // extra for null char
  int64_t sum_depths;
  bool visited;
  
  FragElem() : left_gptr(nullptr), right_gptr(nullptr), sum_depths(0), visited(false) {
    frag_seq[0] = 0;
  }
};

enum class WalkStatus { RUNNING = '-', DEADEND = 'X', FORK = 'F', CONFLICT = 'O', VISITED = 'V', REPEAT = 'R'};

enum class WalkDirn {LEFT = 0, RIGHT = 1};
#define DIRN_STR(d) ((d) == WalkDirn::LEFT ? "left" : "right")

struct StepInfo {
  WalkStatus walk_status;
  ext_count_t count;
  char left, right;
  global_ptr<FragElem> frag_elem;
};

struct WalkTermStats {
  int64_t num_deadends, num_forks, num_conflicts, num_visited, num_repeats;
  
  void update(WalkStatus walk_status) {
    switch (walk_status) {
      case WalkStatus::DEADEND: num_deadends++; break;
      case WalkStatus::FORK: num_forks++; break;
      case WalkStatus::CONFLICT: num_conflicts++; break;
      case WalkStatus::VISITED: num_visited++; break;
      case WalkStatus::REPEAT: num_repeats++; break;
    }
  }  
  
  void print() {
    auto all_num_deadends = reduce_one(num_deadends, op_fast_add, 0).wait();
    auto all_num_forks = reduce_one(num_forks, op_fast_add, 0).wait();
    auto all_num_conflicts = reduce_one(num_conflicts, op_fast_add, 0).wait();
    auto all_num_visited = reduce_one(num_visited, op_fast_add, 0).wait();
    auto all_num_repeats = reduce_one(num_repeats, op_fast_add, 0).wait();
    auto tot_ends = all_num_forks + all_num_deadends + all_num_conflicts + all_num_visited + all_num_repeats;
    SLOG_VERBOSE("Walk statistics:\n");
    SLOG_VERBOSE("  deadends:  ", perc_str(all_num_deadends, tot_ends), "\n");
    SLOG_VERBOSE("  forks:     ", perc_str(all_num_forks, tot_ends), "\n");
    SLOG_VERBOSE("  conflicts: ", perc_str(all_num_conflicts, tot_ends), "\n");
    SLOG_VERBOSE("  visited:   ", perc_str(all_num_visited, tot_ends), "\n");
    SLOG_VERBOSE("  repeats:   ", perc_str(all_num_repeats, tot_ends), "\n");
  }
};

static StepInfo get_next_step(dist_object<KmerDHT> &kmer_dht, Kmer kmer, WalkDirn dirn, char prev_ext, bool revisit_allowed,
                              global_ptr<FragElem> frag_elem) {
  auto kmer_rc = kmer.revcomp();
  bool is_rc = false;
  if (kmer_rc < kmer) {
    kmer = kmer_rc;
    is_rc = true;
  }
  return rpc(kmer_dht->get_kmer_target_rank(kmer), 
             [](dist_object<KmerDHT> &kmer_dht, MerArray merarr, WalkDirn dirn, char prev_ext, bool revisit_allowed, 
                bool is_rc, global_ptr<FragElem> frag_elem) -> StepInfo {
               Kmer kmer(merarr);
               KmerCounts *kmer_counts = kmer_dht->get_local_kmer_counts(kmer);
               // this kmer doesn't exist, abort
               if (!kmer_counts) return {.walk_status = WalkStatus::DEADEND};
               if (kmer_counts->uutig_frag == frag_elem) {
                 if (!revisit_allowed) return {.walk_status = WalkStatus::REPEAT};
                 kmer_counts->uutig_frag = nullptr;
               }
               char left = kmer_counts->left;
               char right = kmer_counts->right;
               if (left == 'X' || right == 'X') return {.walk_status = WalkStatus::DEADEND};
               if (left == 'F' || right == 'F') return {.walk_status = WalkStatus::FORK};
               if (is_rc) {
                 if (left != 'X' && left != 'F') left = comp_nucleotide(left);
                 if (right != 'X' && right != 'F') right = comp_nucleotide(right);
                 swap(left, right);
               }
               if (prev_ext && ((dirn == WalkDirn::LEFT && prev_ext != right) || 
                                (dirn == WalkDirn::RIGHT && prev_ext != left))) {
                 return {.walk_status = WalkStatus::CONFLICT};
               }
               // already visited
               if (kmer_counts->uutig_frag) {
                 DBG_TRAVERSE("  walk in dirn ", DIRN_STR(dirn), " term visited frag ", kmer_counts->uutig_frag, 
                              " rc ", is_rc, "\n");
                 return {.walk_status = WalkStatus::VISITED, .count = 0, .left = 0, .right = 0, 
                         .frag_elem = kmer_counts->uutig_frag};
               }
               // mark as visited
               kmer_counts->uutig_frag = frag_elem;
               return {.walk_status = WalkStatus::RUNNING, .count = kmer_counts->count, .left = left, .right = right};
             }, kmer_dht, kmer.get_array(), dirn, prev_ext, revisit_allowed, is_rc, frag_elem).wait();
}

WalkStatus walk_dbjg_dirn(dist_object<KmerDHT> &kmer_dht, Kmer kmer, WalkDirn dirn, string &uutig, int64_t &sum_depths,
                          global_ptr<FragElem> frag_elem) {
  auto kmer_str = kmer.to_string();
  char prev_ext = 0;
  char next_ext = (dirn == WalkDirn::LEFT ? kmer_str.front() : kmer_str.back());
  bool revisit_allowed = (dirn == WalkDirn::LEFT ? false : true);
  if (dirn == WalkDirn::RIGHT) uutig += kmer_str.substr(1, kmer_str.length() - 2);
  while (true) {
    auto step_info = get_next_step(kmer_dht, kmer, dirn, prev_ext, revisit_allowed, frag_elem);
    revisit_allowed = false;
    if (step_info.walk_status != WalkStatus::RUNNING) {
      // reverse it because we were walking backwards
      if (dirn == WalkDirn::LEFT) reverse(uutig.begin(), uutig.end());
      return step_info.walk_status;
    }
    sum_depths += step_info.count;
    // add the last nucleotide to the uutig
    uutig += next_ext;
    // now attempt to walk to next kmer
    kmer_str = kmer.to_string();
    // get next extension
    if (dirn == WalkDirn::LEFT) {
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

/*
  while (true) {
    char ext = (dirn == WalkDirn::RIGHT ? right_ext : left_ext);
    char other_ext = (dirn == WalkDirn::RIGHT ? left_ext : right_ext);
    if (ext == 'X') return WalkStatus::DEADEND;
    if (ext == 'F') return WalkStatus::FORK;
    // construct the next kmer
    if (dirn == WalkDirn::RIGHT) kmer = kmer.forward_base(ext);
    else kmer = kmer.backward_base(ext);
    // try to add that kmer and get the next step
    auto step_info = get_next_step(kmer_dht, kmer, other_ext, frag_elem, dirn);
    if (step_info.walk_status != WalkStatus::RUNNING) {
      if (step_info.walk_status == WalkStatus::VISITED) {
        if (dirn == WalkDirn::RIGHT) frag_elem.local()->right_gptr = step_info.frag_elem;
        else frag_elem.local()->left_gptr = step_info.frag_elem;
      }
      return step_info.walk_status;
    }
    // still running
    sum_depths += step_info.count;
    if (uutig.length() >= MAX_FRAG_LEN) {
      // FIXME: allocate new fragment element here
      DIE("Fragment sequence too long");
    }
    uutig += ext;
    // add the next nucleotide to the fragment
    if (dirn == WalkDirn::RIGHT) {
      right_ext = step_info.right;
      left_ext = kmer.to_string().front();
    } else {
      left_ext = step_info.left;
      right_ext = kmer.to_string().back();
    }
  }
}
*/

vector<global_ptr<FragElem>> construct_fragments(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht) {
  Timer timer(__FILEFUNC__);
  WalkTermStats walk_term_stats = {0};
  vector<global_ptr<FragElem>> frag_elems;
  ProgressBar progbar(kmer_dht->get_local_num_kmers(), "Traversing deBruijn graph");
  for (auto it = kmer_dht->local_kmers_begin(); it != kmer_dht->local_kmers_end(); it++) {
    progress();
    progbar.update();
    auto kmer = it->first;
    auto kmer_counts = &it->second;
    // don't start any new walk if this kmer has already been visited
    if (kmer_counts->uutig_frag) continue;
    // don't start walks on kmers without extensions in both directions
    if (kmer_counts->left == 'X' || kmer_counts->right == 'X' || kmer_counts->left == 'F' || kmer_counts->right == 'F')
      continue;
    // add this kmer to the frag elem
    global_ptr<FragElem> frag_elem = new_<FragElem>();
    DBG_TRAVERSE("Frag ", frag_elem, " start dbjg traverse at kmer ", kmer, "\n");
    frag_elems.push_back(frag_elem);
    FragElem *lfrag_elem = frag_elem.local();
    // claim this kmer for this frag elem
    kmer_counts->uutig_frag = frag_elem;
    string uutig;
    int64_t sum_depths = 0;
    WalkStatus walk_status = walk_dbjg_dirn(kmer_dht, kmer, WalkDirn::LEFT, uutig, sum_depths, frag_elem);
    DBG_TRAVERSE("  left walk terminates with '", (char)walk_status, "'\n");
    walk_term_stats.update(walk_status);
    walk_status = walk_dbjg_dirn(kmer_dht, kmer, WalkDirn::RIGHT, uutig, sum_depths, frag_elem);
    DBG_TRAVERSE("  right walk terminates with '", (char)walk_status, "'\n");
    walk_term_stats.update(walk_status);
    strcpy(lfrag_elem->frag_seq, uutig.c_str());
    lfrag_elem->sum_depths = sum_depths;
  }
  progbar.done();
  barrier();
  walk_term_stats.print();
  return frag_elems;
}

bool is_overlap(const string &left_seq, const string &right_seq, int overlap_len) {
  return (left_seq.compare(left_seq.length() - overlap_len, overlap_len, right_seq, 0, overlap_len) == 0);
}

bool walk_frags(WalkDirn dirn, int kmer_len, global_ptr<FragElem> start_frag_elem, string &uutig, int64_t &walk_steps,
                vector<FragElem*> &lfrag_elems_visited, int64_t &num_overlaps_rc) {
  auto lfrag_elem = start_frag_elem.local();
  lfrag_elems_visited.push_back(lfrag_elem);
  walk_steps = 1;
  global_ptr<FragElem> next_gptr = (dirn == WalkDirn::RIGHT ? lfrag_elem->right_gptr : lfrag_elem->left_gptr);
#ifdef DEBUG
  // for checking that we haven't got a bug - frags should never be revisited in a walk
  HASH_TABLE<global_ptr<FragElem>, bool> visited;
#endif
  DBG_TRAVERSE("Walking ", DIRN_STR(dirn), " from ", start_frag_elem, "\n");
  while (next_gptr) {
    DBG_TRAVERSE("  next_gptr ", next_gptr, "\n");
    if (next_gptr.where() > rank_me()) {
      DBG_TRAVERSE("    -> DROP: owner ", next_gptr.where(), " > ", rank_me(), "\n");
      return false;
    }
#ifdef DEBUG
    if (visited.find(next_gptr) != visited.end()) {
      DBG_TRAVERSE("    DIE: repeat - ", next_gptr, "\n");
      DIE("repeat: ", next_gptr);
    }
    visited[next_gptr] = true;
#endif 
    if (next_gptr.where() == rank_me()) {
      if (next_gptr.local()->visited) DIE("Should not be already visited");
      lfrag_elems_visited.push_back(next_gptr.local());
    } 
    FragElem next_frag_elem = rget(next_gptr).wait();
    DBG_TRAVERSE("    left gptr ", next_frag_elem.left_gptr, " right gptr ", next_frag_elem.right_gptr, "\n");
    string next_frag_seq(next_frag_elem.frag_seq);
    if (dirn == WalkDirn::RIGHT) {
      if (is_overlap(uutig, next_frag_seq, kmer_len - 1)) {
        uutig += next_frag_seq.substr(kmer_len - 1);
        DBG_TRAVERSE("    GOOD right overlap\n");
      } else {
        num_overlaps_rc++;      
        string next_frag_seq_rc = revcomp(next_frag_seq);
        if (!is_overlap(uutig, next_frag_seq_rc, kmer_len - 1)) {
          DBG_TRAVERSE("    DIE: right fragments don't overlap\n");
          DIE("right fragments RC don't overlap:",
              "\nnext_gptr  ", next_gptr, 
              "\nleft_gptr  ", next_frag_elem.left_gptr,
              "\nright_gptr ", next_frag_elem.right_gptr,
              "\nuutig:  ", uutig, 
              "\nnextrc: ", next_frag_seq_rc, 
              "\nnext:   ", next_frag_seq);
        }
        DBG_TRAVERSE("    REVCOMP right overlap\n");
        uutig += next_frag_seq_rc.substr(kmer_len - 1);
        // now flip dirn since it was revcomped
        dirn = WalkDirn::LEFT;
        DBG_TRAVERSE("    flipping from RIGHT to LEFT\n");
      }
    } else {
      if (is_overlap(next_frag_seq, uutig, kmer_len - 1)) {
        DBG_TRAVERSE("    GOOD left overlap\n");
        uutig.insert(0, next_frag_seq.substr(0, next_frag_seq.length() - (kmer_len - 1)));
      } else {
        num_overlaps_rc++;      
        string next_frag_seq_rc = revcomp(next_frag_seq);
        if (!is_overlap(uutig, next_frag_seq_rc, kmer_len - 1)) {
          DBG_TRAVERSE("    DIE: left fragments don't overlap\n");
          DIE("left fragments RC don't overlap ", next_frag_seq_rc.length(), 
              "\nuutig:  ", uutig, 
              "\nnextrc: ", next_frag_seq_rc, 
              "\nnext:   ", next_frag_seq);
        }
        DBG_TRAVERSE("    REVCOMP left overlap\n");
        uutig += next_frag_seq_rc.substr(kmer_len - 1);
        dirn = WalkDirn::RIGHT;
        DBG_TRAVERSE("    flipping from LEFT to RIGHT\n");
      }
    }
    walk_steps++;
    next_gptr = (dirn == WalkDirn::RIGHT ? next_frag_elem.right_gptr : next_frag_elem.left_gptr);
  }
  return true;
}
        
void connect_fragments(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, vector<global_ptr<FragElem>> &frag_elems, 
                       Contigs &my_uutigs) {
  Timer timer(__FILEFUNC__);
  // connect the fragments and put into my_uutigs
  int64_t num_drops = 0, num_steps = 0, max_steps = 0, num_overlaps_rc = 0;
  my_uutigs.clear();
  ProgressBar progbar(frag_elems.size(), "Connecting fragments");
  int64_t num_frags = 0, num_left_links = 0, num_right_links = 0;
  for (auto it = frag_elems.begin(); it != frag_elems.end(); ++it) {
    progbar.update();
    auto frag_elem = *it;
    auto lfrag_elem = frag_elem.local();
    DBG_TRAVERSE("Fragment: ", frag_elem, " left ", lfrag_elem->left_gptr, " right ", lfrag_elem->right_gptr, "\n");
    if (lfrag_elem->left_gptr) {
      num_left_links++;
      FragElem other_frag_elem = rget(lfrag_elem->left_gptr).wait();
      DBG_TRAVERSE("  on left ", lfrag_elem->left_gptr, " left ", other_frag_elem.left_gptr, 
                   " right ", other_frag_elem.right_gptr, "\n");
    }
    if (lfrag_elem->right_gptr) {
      num_right_links++;
      FragElem other_frag_elem = rget(lfrag_elem->right_gptr).wait();
      DBG_TRAVERSE("  on right ", lfrag_elem->right_gptr, " left ", other_frag_elem.left_gptr, 
                   " right ", other_frag_elem.right_gptr, "\n");
    }
    if (lfrag_elem->visited) continue;
    string uutig(lfrag_elem->frag_seq);
    int64_t sum_depths = lfrag_elem->sum_depths;
#ifdef DEBUG      
    if (lfrag_elem->left_gptr && lfrag_elem->right_gptr && lfrag_elem->left_gptr == lfrag_elem->right_gptr) {
      DIE("left gptr is equal to right gptr ", lfrag_elem->left_gptr, " == ", lfrag_elem->right_gptr);
    }
    // check for errors in kmers
    auto kmers = Kmer::get_kmers(kmer_len, uutig);
    for (auto kmer : kmers) {
      assert(kmer_dht->get_kmer_uutig_frag(kmer) == frag_elem);
    }
#endif
    Contig contig = {0, uutig, (double)sum_depths / (uutig.length() - kmer_len + 2)};
    my_uutigs.add_contig(contig);
    
    /*
    int64_t walk_steps;
    vector<FragElem*> lfrag_elems_visited;
    if (walk_frags(WalkDirn::RIGHT, kmer_len, frag_elem, uutig, walk_steps, lfrag_elems_visited, num_overlaps_rc) &&
            walk_frags(WalkDirn::LEFT, kmer_len, frag_elem, uutig, walk_steps, lfrag_elems_visited, num_overlaps_rc)) {
      num_steps += walk_steps;
      max_steps = max(walk_steps, max_steps);
      Contig contig = {0, uutig, (double)sum_depths / (uutig.length() - kmer_len + 2)};
      my_uutigs.add_contig(contig);
      // the walk is successful, so set the visited for all the local elems
      for (auto &elem : lfrag_elems_visited) {
        elem->visited = true;
      }
    } else {
      num_drops++;
    }
    */
    num_frags++;
  }
  progbar.done();
  barrier();
  auto all_num_frags = reduce_one(frag_elems.size(), op_fast_add, 0).wait();
  auto all_num_left_links = reduce_one(num_left_links, op_fast_add, 0).wait();
  auto all_num_right_links = reduce_one(num_right_links, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", all_num_frags, " uutig fragments with ", all_num_left_links, " left links and ",
               all_num_right_links, " right links\n");
  auto all_num_drops = reduce_one(num_drops, op_fast_add, 0).wait();
  auto all_num_steps = reduce_one(num_steps, op_fast_add, 0).wait();
  auto all_max_steps = reduce_one(max_steps, op_fast_max, 0).wait();
  auto all_num_overlaps_rc = reduce_one(num_overlaps_rc, op_fast_add, 0).wait();
  auto all_num_uutigs = reduce_one(my_uutigs.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Completed ", all_num_uutigs, " fragment walks:\n");
  SLOG_VERBOSE("  drops:     ", perc_str(all_num_drops, all_num_frags), "\n");
  SLOG_VERBOSE("  steps:     ", all_num_steps, "\n");
  SLOG_VERBOSE("  average walk length: ", (double)all_num_steps / all_num_uutigs, "\n");
  SLOG_VERBOSE("  max walk length:     ", all_max_steps, "\n");
  SLOG_VERBOSE("  RC overlaps:         ", all_num_overlaps_rc, "\n");
}

void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs) {
  Timer timer(__FILEFUNC__);
  barrier();
  vector<global_ptr<FragElem>> frag_elems = construct_fragments(kmer_len, kmer_dht);
  connect_fragments(kmer_len, kmer_dht, frag_elems, my_uutigs);
  // now get unique ids for the uutigs
  atomic_domain<size_t> ad({atomic_op::fetch_add, atomic_op::load});
  global_ptr<size_t> counter = nullptr;
  if (!rank_me()) counter = new_<size_t>(0);
  counter = broadcast(counter, 0).wait();
  size_t my_counter = ad.fetch_add(counter, my_uutigs.size(), memory_order_relaxed).wait();
  // wait until all ranks have updated the global counter
  barrier();
  ad.destroy();
  // set the unique ids
  int64_t cid = my_counter;
  for (auto it = my_uutigs.begin(); it != my_uutigs.end(); ++it) {
    it->id = cid;
    cid++;
#ifdef DEBUG
    // check that all kmers in sequence actually exist
    auto kmers = Kmer::get_kmers(kmer_len, it->seq);
    for (auto kmer : kmers) {
      if (!kmer_dht->kmer_exists(kmer)) DIE("kmer not found in dht");
    }
#endif
  }
  SLOG_VERBOSE("Constructed ", reduce_one(my_uutigs.size(), op_fast_add, 0).wait(), " UU-tigs\n");
}

