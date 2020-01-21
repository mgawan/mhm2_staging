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

enum class Dirn { LEFT, RIGHT };
#define DIRN_STR(d) ((d) == Dirn::LEFT ? "left" : "right")

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

static bool check_kmers(const string &seq, dist_object<KmerDHT> &kmer_dht, int kmer_len) {
  auto kmers = Kmer::get_kmers(kmer_len, seq);
  for (auto kmer : kmers) {
    if (!kmer_dht->kmer_exists(kmer)) return false;
  }
  return true;
}

static future<StepInfo> get_next_step(dist_object<KmerDHT> &kmer_dht, Kmer kmer, Dirn dirn, char prev_ext, 
                                      global_ptr<FragElem> frag_elem_gptr, bool revisit_allowed) {
  auto kmer_rc = kmer.revcomp();
  bool is_rc = false;
  if (kmer_rc < kmer) {
    kmer = kmer_rc;
    is_rc = true;
  }
  return rpc(kmer_dht->get_kmer_target_rank(kmer), 
             [](dist_object<KmerDHT> &kmer_dht, MerArray merarr, Dirn dirn, char prev_ext, bool revisit_allowed, 
                bool is_rc, global_ptr<FragElem> frag_elem_gptr) -> StepInfo {
               Kmer kmer(merarr);
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
             }, kmer_dht, kmer.get_array(), dirn, prev_ext, revisit_allowed, is_rc, frag_elem_gptr);
}

static global_ptr<FragElem> traverse_dirn(dist_object<KmerDHT> &kmer_dht, Kmer kmer, global_ptr<FragElem> frag_elem_gptr, 
                                          Dirn dirn, string &uutig, int64_t &sum_depths, WalkTermStats &walk_term_stats) {
  auto kmer_str = kmer.to_string();
  char prev_ext = 0;
  char next_ext = (dirn == Dirn::LEFT ? kmer_str.front() : kmer_str.back());
  bool revisit_allowed = (dirn == Dirn::LEFT ? false : true);
  if (dirn == Dirn::RIGHT) uutig += kmer_str.substr(1, kmer_str.length() - 2);
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

static void construct_frags(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, vector<global_ptr<FragElem>> &frag_elems) {
  Timer timer(__FILEFUNC__);
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
    frag_elem_gptr.local()->frag_seq = new_array<char>(uutig.length() + 1);
    strcpy(frag_elem_gptr.local()->frag_seq.local(), uutig.c_str());
    frag_elem_gptr.local()->frag_seq.local()[uutig.length()] = 0;
    frag_elem_gptr.local()->frag_len = uutig.length();
    frag_elem_gptr.local()->sum_depths = sum_depths;      
    frag_elem_gptr.local()->left_gptr = left_gptr;
    frag_elem_gptr.local()->right_gptr = right_gptr;
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
  char *buf = new char[frag_elem.frag_len + 1];
  rget(frag_elem.frag_seq, buf, frag_elem.frag_len + 1);
  string frag_seq(buf);
  delete[] buf;
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

static void clean_frag_links(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, vector<global_ptr<FragElem>> &frag_elems) {
  Timer timer(__FILEFUNC__);
  // put all the uutigs found by this rank into my_uutigs
  int64_t num_uutigs = 0, num_equal_links = 0, num_non_recip = 0, num_short = 0,
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

static bool walk_frags_dirn(Dirn dirn, unsigned kmer_len, global_ptr<FragElem> frag_elem_gptr, string &uutig,
                            int64_t &walk_steps, vector<FragElem*> my_frag_elems_visited) {
//  auto frag_elem = frag_elem_gptr.local();
//  uutig = string(frag_elem->frag_seq.local());

  auto start_frag_elem = frag_elem_gptr.local();
  my_frag_elems_visited.push_back(start_frag_elem);
  if (dirn == Dirn::LEFT) uutig = string(frag_elem_gptr.local()->frag_seq.local());
  walk_steps = 1;
  global_ptr<FragElem> prev_gptr = frag_elem_gptr;
  FragElem prev_frag_elem = *start_frag_elem;
  global_ptr<FragElem> next_gptr = (dirn == Dirn::RIGHT ? start_frag_elem->right_gptr : start_frag_elem->left_gptr);
#ifdef DEBUG
  // for checking that we haven't got a bug - frags should never be revisited in a walk
  HASH_TABLE<global_ptr<FragElem>, bool> visited;
#endif
  if (next_gptr) 
    DBG_TRAVERSE("Walking ", DIRN_STR(dirn), " from ", frag_elem_gptr, " (length ", uutig.length(), ") to ", next_gptr, "\n");
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
    string uutig_prev = uutig;
#endif 
    if (next_gptr.where() == rank_me()) {
      if (next_gptr.local()->visited) {
        DBG_TRAVERSE("    DIE: gptr", next_gptr, " is already visited\n");
        DIE("gptr ", next_gptr, " should not be already visited");
      }
      my_frag_elems_visited.push_back(next_gptr.local());
    } 
    FragElem next_frag_elem = rget(next_gptr).wait();
    DBG_TRAVERSE("    left gptr ", next_frag_elem.left_gptr, " right gptr ", next_frag_elem.right_gptr, "\n");
    string next_frag_seq = get_frag_seq(next_frag_elem);
    if (dirn == Dirn::LEFT) {
      if (next_frag_elem.left_gptr == prev_gptr) DBG_TRAVERSE("  Need to flip direction\n");
      int ext = next_frag_seq.length() - kmer_len + 1;
      if (!prev_frag_elem.left_is_rc) {
        uutig.insert(0, next_frag_seq.substr(0, next_frag_seq.length() - kmer_len + 1));
        DBG_TRAVERSE("    extending to the left by ", ext,
                     "\n  ", string(ext, ' '), uutig_prev,
                     "\n  ", next_frag_seq, 
                     "\n  ", uutig, "\n");
      } else {
        string next_frag_seq_rc = revcomp(next_frag_seq);
        uutig.insert(0, next_frag_seq_rc.substr(0, next_frag_seq.length() - (kmer_len - 1)));
        DBG_TRAVERSE("    extending to the left RC by ", (next_frag_seq_rc.length() - kmer_len + 1), 
                     "\n  ", string(ext, ' '), uutig_prev,
                     "\n  ", next_frag_seq_rc, 
                     "\n  ", uutig, "\n");
#ifdef DEBUG        
        if (!is_overlap(next_frag_seq_rc, uutig_prev, kmer_len - 1)) 
          DBG_TRAVERSE("  INCORRECT left overlap RC\n", next_frag_seq, "\n");
#endif
        return true;
        
        // now flip dirn since it was revcomped
        dirn = Dirn::RIGHT;
        DBG_TRAVERSE("    flipping from LEFT to RIGHT\n");
      }
    } else {
      if (next_frag_elem.right_gptr == prev_gptr) DBG_TRAVERSE("  Need to flip direction\n");
      int ext = next_frag_seq.length() - kmer_len + 1;
      if (!prev_frag_elem.right_is_rc) {
        uutig += next_frag_seq.substr(kmer_len - 1);
        DBG_TRAVERSE("    extending to the right by ", ext,
                     "\n  ", uutig_prev, 
                     "\n  ", string(uutig_prev.length() - kmer_len + 1, ' '), next_frag_seq, 
                     "\n  ", uutig, "\n");
      } else {
        string next_frag_seq_rc = revcomp(next_frag_seq);
        uutig += next_frag_seq_rc.substr(kmer_len - 1);
        DBG_TRAVERSE("    extending to the right RC by ", (next_frag_seq.length() - kmer_len + 1), 
                     "\n  ", uutig_prev, 
                     "\n  ", string(uutig_prev.length() - kmer_len + 1, ' '), next_frag_seq_rc, 
                     "\n  ", uutig, "\n");
#ifdef DEBUG        
        if (!is_overlap(uutig_prev, next_frag_seq_rc, kmer_len - 1)) 
          DBG_TRAVERSE("  INCORRECT right overlap RC\n", next_frag_seq, "\n");
#endif
        
        return true;
        
        dirn = Dirn::LEFT;
        DBG_TRAVERSE("    flipping from RIGHT to LEFT\n");
      }
    }
    walk_steps++;
    prev_gptr = next_gptr;
    prev_frag_elem = next_frag_elem;
    next_gptr = (dirn == Dirn::RIGHT ? next_frag_elem.right_gptr : next_frag_elem.left_gptr);
  }

  return true;
}


static void connect_frags(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, vector<global_ptr<FragElem>> &frag_elems, 
                          Contigs &my_uutigs) {
  Timer timer(__FILEFUNC__);
  int64_t num_steps = 0, max_steps = 0, num_drops = 0;
  ProgressBar progbar(frag_elems.size(), "Connecting fragments");
  for (auto frag_elem_gptr : frag_elems) {
    progbar.update();
    FragElem *frag_elem = frag_elem_gptr.local();
    if (frag_elem->frag_len < kmer_len) continue;
    int64_t walk_steps;
    vector<FragElem*> my_frag_elems_visited;
    string uutig;
    int64_t sum_depths = 0;
    if (walk_frags_dirn(Dirn::LEFT, kmer_len, frag_elem_gptr, uutig, walk_steps, my_frag_elems_visited) &&
        walk_frags_dirn(Dirn::RIGHT, kmer_len, frag_elem_gptr, uutig, walk_steps, my_frag_elems_visited)) {
      num_steps += walk_steps;
      if (walk_steps > max_steps) DBG_TRAVERSE("MAX path length ", walk_steps, "\n");
      max_steps = max(walk_steps, max_steps);
      Contig contig = {0, uutig, (double)sum_depths / (uutig.length() - kmer_len + 2)};
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
  auto all_num_uutigs = reduce_one(my_uutigs.size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Constructed ", all_num_uutigs, " uutigs with ", (double)all_num_steps / all_num_uutigs, 
               " avg path length (max ", all_max_steps, "), dropped ", perc_str(all_num_drops, all_num_uutigs), " paths\n");
}

void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs) {
  Timer timer(__FILEFUNC__);
  vector<global_ptr<FragElem>> frag_elems;
  construct_frags(kmer_len, kmer_dht, frag_elems);
  clean_frag_links(kmer_len, kmer_dht, frag_elems);
  // put all the uutigs found by this rank into my_uutigs
  my_uutigs.clear();
  connect_frags(kmer_len, kmer_dht, frag_elems, my_uutigs);
  barrier();
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
  }
#ifdef DEBUG
  ProgressBar progbar(frag_elems.size(), "Checking kmers in uutigs");
  for (auto uutig : my_uutigs) {
    progbar.update();
    if (!check_kmers(uutig.seq, kmer_dht, kmer_len)) DIE("kmer not found in uutig");
  }
  progbar.done();
#endif
}
