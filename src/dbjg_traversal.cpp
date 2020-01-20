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

enum class TraverseDirn { LEFT, RIGHT };

enum class WalkStatus { RUNNING = '-', DEADEND = 'X', FORK = 'F', CONFLICT = 'O', REPEAT = 'R', VISITED = 'V'};

struct FragElem {
  global_ptr<FragElem> left_gptr, right_gptr;
  global_ptr<char> frag_seq;
  int frag_len;
  int64_t sum_depths;
  bool visited;
  
  FragElem() : left_gptr(nullptr), right_gptr(nullptr), frag_seq(nullptr), frag_len(0), sum_depths(0), visited(false) {}
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

static future<StepInfo> get_next_step(dist_object<KmerDHT> &kmer_dht, Kmer kmer, TraverseDirn dirn, char prev_ext, 
                                      global_ptr<FragElem> frag_elem_gptr, bool revisit_allowed) {
  auto kmer_rc = kmer.revcomp();
  bool is_rc = false;
  if (kmer_rc < kmer) {
    kmer = kmer_rc;
    is_rc = true;
  }
  return rpc(kmer_dht->get_kmer_target_rank(kmer), 
             [](dist_object<KmerDHT> &kmer_dht, MerArray merarr, TraverseDirn dirn, char prev_ext, bool revisit_allowed, 
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
               if (prev_ext && ((dirn == TraverseDirn::LEFT && prev_ext != right) || 
                                (dirn == TraverseDirn::RIGHT && prev_ext != left))) 
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
                                          TraverseDirn dirn, string &uutig, int64_t &sum_depths, WalkTermStats &walk_term_stats) {
  auto kmer_str = kmer.to_string();
  char prev_ext = 0;
  char next_ext = (dirn == TraverseDirn::LEFT ? kmer_str.front() : kmer_str.back());
  bool revisit_allowed = (dirn == TraverseDirn::LEFT ? false : true);
  if (dirn == TraverseDirn::RIGHT) uutig += kmer_str.substr(1, kmer_str.length() - 2);
  while (true) {
    //progress();
    auto step_info = get_next_step(kmer_dht, kmer, dirn, prev_ext, frag_elem_gptr, revisit_allowed).wait();
    revisit_allowed = false;
    if (step_info.walk_status != WalkStatus::RUNNING) {
      walk_term_stats.update(step_info.walk_status);
      // reverse it because we were walking backwards
      if (dirn == TraverseDirn::LEFT) reverse(uutig.begin(), uutig.end());
      return step_info.visited_frag_elem_gptr;  
    }
    sum_depths += step_info.count;
    // add the last nucleotide to the uutig
    uutig += next_ext;
    // now attempt to walk to next kmer
    kmer_str = kmer.to_string();
    // get next extension
    if (dirn == TraverseDirn::LEFT) {
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

void construct_frags(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, vector<global_ptr<FragElem>> &frag_elems) {
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
    auto left_gptr = traverse_dirn(kmer_dht, kmer, frag_elem_gptr, TraverseDirn::LEFT, uutig, sum_depths, walk_term_stats);
    auto right_gptr = traverse_dirn(kmer_dht, kmer, frag_elem_gptr, TraverseDirn::RIGHT, uutig, sum_depths, walk_term_stats);
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

bool is_overlap(const string &left_seq, const string &right_seq, int overlap_len) {
  return (left_seq.compare(left_seq.length() - overlap_len, overlap_len, right_seq, 0, overlap_len) == 0);
}

bool check_nb_gptr(TraverseDirn dirn, global_ptr<FragElem> nb_gptr, string &uutig, int kmer_len, 
                   int64_t &num_overlaps, int64_t &num_overlaps_rc, int64_t num_non_recip) {
  if (nb_gptr) {
    FragElem nb_frag_elem = rget(nb_gptr).wait();
    char *buf = new char[nb_frag_elem.frag_len + 1];
    rget(nb_frag_elem.frag_seq, buf, nb_frag_elem.frag_len + 1);
    string nb_frag_seq(buf);
    delete[] buf;
    string *s1 = (dirn == TraverseDirn::LEFT ? &nb_frag_seq : &uutig);
    string *s2 = (dirn == TraverseDirn::LEFT ? &uutig : &nb_frag_seq);
    if (is_overlap(*s1, *s2, kmer_len - 1)) {
      num_overlaps++;
      if ((dirn == TraverseDirn::LEFT ? nb_frag_elem.right_gptr : nb_frag_elem.left_gptr) == nb_gptr) return true;
      num_non_recip++;
      return false;
    }
    auto nb_frag_seq_rc = revcomp(nb_frag_seq);
    s1 = (dirn == TraverseDirn::LEFT ? &nb_frag_seq_rc : &uutig);
    s2 = (dirn == TraverseDirn::LEFT ? &uutig : &nb_frag_seq_rc);
    if (is_overlap(*s1, *s2, kmer_len - 1)) {
      num_overlaps_rc++;
      if ((dirn == TraverseDirn::LEFT ? nb_frag_elem.left_gptr : nb_frag_elem.right_gptr) == nb_gptr) return true;
      num_non_recip++;
      return true;
    }
    DBG_TRAVERSE("No ", (dirn == TraverseDirn::LEFT ? "left" : "right"), " overlap:\n", uutig, 
                 "\n", nb_frag_seq, "\n", nb_frag_seq_rc, "\n");
    return false;
  }
  return false;
}

int64_t print_link_stats(int64_t num_links, int64_t num_overlaps, int64_t num_overlaps_rc, const string &dirn_str) {
  auto all_num_links = reduce_one(num_links, op_fast_add, 0).wait();
  auto all_num_overlaps = reduce_one(num_overlaps, op_fast_add, 0).wait();
  auto all_num_overlaps_rc = reduce_one(num_overlaps_rc, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", all_num_links, " ", dirn_str, " links with ",
               perc_str(all_num_overlaps, all_num_links), " overlaps and ", 
               perc_str(all_num_overlaps_rc, all_num_links), " revcomped overlaps\n");
  return all_num_links;
}
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs) {
  Timer timer(__FILEFUNC__);
  vector<global_ptr<FragElem>> frag_elems;
  construct_frags(kmer_len, kmer_dht, frag_elems);
  // put all the uutigs found by this rank into my_uutigs
  int64_t num_uutigs = 0, num_equal_links = 0, num_non_recip = 0, num_short = 0,
          num_left_links = 0, num_left_overlaps = 0, num_left_overlaps_rc = 0, 
          num_right_links = 0, num_right_overlaps = 0, num_right_overlaps_rc = 0;
  my_uutigs.clear();
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
      continue;      
    }
    if (!check_nb_gptr(TraverseDirn::LEFT, frag_elem->left_gptr, uutig, kmer_len, num_left_overlaps, 
                       num_left_overlaps_rc, num_non_recip))
      frag_elem->left_gptr = nullptr;
    if (!check_nb_gptr(TraverseDirn::RIGHT, frag_elem->right_gptr, uutig, kmer_len, num_right_overlaps, 
                       num_right_overlaps_rc, num_non_recip))
      frag_elem->right_gptr = nullptr;
#ifdef DEBUG
    // check that all kmers in sequence actually exist
    auto kmers = Kmer::get_kmers(kmer_len, uutig);
    for (auto kmer : kmers) {
      if (!kmer_dht->kmer_exists(kmer)) DIE("kmer not found in dht");
    }
#endif
    num_uutigs++;
    Contig contig = {0, uutig, (double)frag_elem->sum_depths / (uutig.length() - kmer_len + 2)};
    my_uutigs.add_contig(contig);
  }
  progbar.done();
  barrier();
  auto all_num_uutigs = reduce_one(num_uutigs, op_fast_add, 0).wait();
  auto all_num_short = reduce_one(num_short, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", all_num_uutigs, " uutigs of which ", perc_str(all_num_short, all_num_uutigs), " are short\n");
  auto all_num_left = print_link_stats(num_left_links, num_left_overlaps, num_left_overlaps_rc, "left");
  auto all_num_right = print_link_stats(num_right_links, num_right_overlaps, num_right_overlaps_rc, "right");
  auto all_num_equal_links = reduce_one(num_equal_links, op_fast_add, 0).wait();
  auto all_num_non_recip = reduce_one(num_non_recip, op_fast_add, 0).wait();
  SLOG_VERBOSE("There are ", perc_str(all_num_equal_links, all_num_left + all_num_right), " equal left and right links\n");
  SLOG_VERBOSE("There are ", perc_str(all_num_non_recip, all_num_left + all_num_right), " non-reciprocating links\n");
  // now get unique ids for the uutigs
  atomic_domain<size_t> ad({atomic_op::fetch_add, atomic_op::load});
  global_ptr<size_t> counter = nullptr;
  if (!rank_me()) counter = new_<size_t>(0);
  counter = broadcast(counter, 0).wait();
  size_t my_counter = ad.fetch_add(counter, num_uutigs, memory_order_relaxed).wait();
  // wait until all ranks have updated the global counter
  barrier();
  ad.destroy();
  // set the unique ids
  int64_t cid = my_counter;
  for (auto it = my_uutigs.begin(); it != my_uutigs.end(); ++it) {
    it->id = cid;
    cid++;
  }
  SLOG_VERBOSE("Constructed ", reduce_one(my_uutigs.size(), op_fast_add, 0).wait(), " UU-tigs\n");
}
