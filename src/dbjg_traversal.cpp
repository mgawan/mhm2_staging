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

//#define DBG_TRAVERSE DBG
#define DBG_TRAVERSE(...)

using namespace std;
using namespace upcxx;

extern ofstream _dbgstream;

enum class TraverseDirn { LEFT, RIGHT };

enum class WalkStatus { RUNNING = '-', DEADEND = 'X', FORK = 'F', CONFLICT = 'O', REPEAT = 'R'};

struct StepInfo {
  WalkStatus walk_status;
  bool drop_walk;
  ext_count_t count;
  char left, right;
};

struct WalkTermStats {
  int64_t num_deadends, num_forks, num_conflicts, num_repeats;
  
  void update(WalkStatus walk_status) {
    switch (walk_status) {
      case WalkStatus::DEADEND: num_deadends++; break;
      case WalkStatus::FORK: num_forks++; break;
      case WalkStatus::CONFLICT: num_conflicts++; break;
      case WalkStatus::REPEAT: num_repeats++; break;
    }
  }  
  
  void print() {
    auto all_num_deadends = reduce_one(num_deadends, op_fast_add, 0).wait();
    auto all_num_forks = reduce_one(num_forks, op_fast_add, 0).wait();
    auto all_num_conflicts = reduce_one(num_conflicts, op_fast_add, 0).wait();
    auto all_num_repeats = reduce_one(num_repeats, op_fast_add, 0).wait();
    auto tot_ends = all_num_forks + all_num_deadends + all_num_conflicts + all_num_repeats;
    SLOG_VERBOSE("Walk statistics:\n");
    SLOG_VERBOSE("  deadends:  ", perc_str(all_num_deadends, tot_ends), "\n");
    SLOG_VERBOSE("  forks:     ", perc_str(all_num_forks, tot_ends), "\n");
    SLOG_VERBOSE("  conflicts: ", perc_str(all_num_conflicts, tot_ends), "\n");
    SLOG_VERBOSE("  repeats:   ", perc_str(all_num_repeats, tot_ends), "\n");
  }
};

static future<StepInfo> get_next_step(dist_object<KmerDHT> &kmer_dht, Kmer kmer, TraverseDirn dirn, char prev_ext, 
                                        int32_t start_walk_us, bool revisit_allowed) {
  auto kmer_rc = kmer.revcomp();
  bool is_rc = false;
  if (kmer_rc < kmer) {
    kmer = kmer_rc;
    is_rc = true;
  }
  return rpc(kmer_dht->get_kmer_target_rank(kmer), 
             [](dist_object<KmerDHT> &kmer_dht, MerArray merarr, TraverseDirn dirn, intrank_t start_rank, char prev_ext, 
                int32_t start_walk_us, bool revisit_allowed, bool is_rc) -> StepInfo {
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
                 if (kmer_counts->visited >= 0 && 
                         (kmer_counts->visited != start_rank || kmer_counts->start_walk_us != start_walk_us)) {
                   // resolve in favor of earliest starting time
                   if (kmer_counts->start_walk_us < start_walk_us || 
                           (kmer_counts->start_walk_us == start_walk_us && kmer_counts->visited > start_rank)) 
                     return {.walk_status = WalkStatus::CONFLICT, .drop_walk = true};
                 }
                 // a repeat, abort (but allowed if traversing right after adding start kmer previously)
                 if (kmer_counts->visited == start_rank && kmer_counts->start_walk_us == start_walk_us && !revisit_allowed) 
                   return {.walk_status = WalkStatus::REPEAT};
                   // mark as visited
                   kmer_counts->visited = start_rank;
                   // set the walk time
                   kmer_counts->start_walk_us = start_walk_us;
                   return {.walk_status = WalkStatus::RUNNING, .drop_walk = false, .count = kmer_counts->count, 
                           .left = left, .right = right};
             }, kmer_dht, kmer.get_array(), dirn, rank_me(), prev_ext, start_walk_us, revisit_allowed, is_rc);
}

static bool traverse_dirn(dist_object<KmerDHT> &kmer_dht, Kmer kmer, int32_t start_walk_us, TraverseDirn dirn, string &uutig, 
                            int64_t &sum_depths, WalkTermStats &walk_term_stats) {
  auto kmer_str = kmer.to_string();
  char prev_ext = 0;
  char next_ext = (dirn == TraverseDirn::LEFT ? kmer_str.front() : kmer_str.back());
  bool revisit_allowed = (dirn == TraverseDirn::LEFT ? false : true);
  if (dirn == TraverseDirn::RIGHT) uutig += kmer_str.substr(1, kmer_str.length() - 2);
  while (true) {
    //progress();
    auto step_info = get_next_step(kmer_dht, kmer, dirn, prev_ext, start_walk_us, revisit_allowed).wait();
    revisit_allowed = false;
    if (step_info.drop_walk) return false;
    if (step_info.walk_status != WalkStatus::RUNNING) {
      walk_term_stats.update(step_info.walk_status);
      // reverse it because we were walking backwards
      if (dirn == TraverseDirn::LEFT) reverse(uutig.begin(), uutig.end());
      return true;
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
  return true;
}

void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs) {
  Timer timer(__FILEFUNC__);
  // allocate space for biggest possible uutig in global storage
  WalkTermStats walk_term_stats = {0};
  int64_t num_drops = 0, num_walks = 0;
  Contigs uutigs;
  barrier();
  {
    ProgressBar progbar(kmer_dht->get_local_num_kmers(), "Traversing deBruijn graph");
    for (auto it = kmer_dht->local_kmers_begin(); it != kmer_dht->local_kmers_end(); it++) {
      progress();
      progbar.update();
      auto kmer = it->first;
      auto kmer_counts = &it->second;
      // don't start any new walk if this kmer has already been visited
      if (kmer_counts->visited > -1) continue;
      // don't start walks on kmers without extensions on both sides
      if (kmer_counts->left == 'X' || kmer_counts->left == 'F' || kmer_counts->right == 'X' || kmer_counts->right == 'F') 
        continue;
      int32_t start_walk_us = kmer_dht->get_time_offset_us();
      string uutig;
      int64_t sum_depths = 0;
      if (!traverse_dirn(kmer_dht, kmer, start_walk_us, TraverseDirn::LEFT, uutig, sum_depths, walk_term_stats)) {
        num_drops++;
        continue;
      }
      if (!traverse_dirn(kmer_dht, kmer, start_walk_us, TraverseDirn::RIGHT, uutig, sum_depths, walk_term_stats)) {
        num_drops++;
        continue;
      }
      Contig contig = {0, uutig, (double)sum_depths / (uutig.length() - kmer_len + 2)};
      uutigs.add_contig(contig);
      num_walks++;
    }
    progbar.done();
    // FIXME: now steal kmers from others???
  }
  barrier();
  walk_term_stats.print();

  // put all the uutigs found by this rank into my_uutigs
  int64_t num_uutigs = 0;
  my_uutigs.clear();
  {
    ProgressBar progbar(uutigs.size(), "Dropping duplicate walks");
    for (auto it = uutigs.begin(); it != uutigs.end(); ++it) {
      auto uutig = it;
      progbar.update();
      if (uutig->seq.length() < kmer_len)
        DIE("uutig length ", uutig->seq.length(), " less than kmer len, seq:\n", uutig->seq);
      // now check to make sure we're the owner of this one - this is after all ranks have finished traversals
      Kmer start_kmer(uutig->seq.substr(0, kmer_len).c_str());
      auto start_kmer_rc = start_kmer.revcomp();
      if (start_kmer_rc < start_kmer) start_kmer = start_kmer_rc;
      auto visited = kmer_dht->get_visited(start_kmer);
      if (visited == -2) DIE("Could not find kmer ", start_kmer, "\n");
      if (visited == -1) DIE("unitialized visited for kmer ", start_kmer, "\n");
      if (visited != rank_me()) {
        num_drops++;
        DBG_TRAVERSE("Drop kmer ", start_kmer, " visited is ", visited, " not my rank ", rank_me(), "\n");
      } else {
        num_uutigs++;
        my_uutigs.add_contig(*uutig);
      }
    }
    progbar.done();
  }
  barrier();
  DBG_TRAVERSE("dropped ", num_drops, " out of ", num_walks, " walks\n");
  auto all_num_drops = reduce_one(num_drops, op_fast_add, 0).wait();
  auto all_num_walks = reduce_one(num_walks, op_fast_add, 0).wait();
  SLOG_VERBOSE("Dropped ", perc_str(all_num_drops, all_num_walks), " duplicate traversals\n");
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


/*
void compute_kmer_ctg_depths(int kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &ctgs) {
  Timer timer(__FILEFUNC__);
  ProgressBar progbar(ctgs.size(), "Computing contig kmer depths");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() >= kmer_len + 2) {
      int num_kmers = ctg->seq.length() - kmer_len;
      ctg->kmer_depths.reserve(num_kmers);
      auto kmers = Kmer::get_kmers(kmer_len, ctg->seq);
      for (auto kmer : kmers) {
        auto kmer_rc = kmer.revcomp();
        if (kmer_rc < kmer) kmer = kmer_rc;
        uint16_t count = kmer_dht->get_kmer_count(kmer);
        assert(count != 0);
        ctg->kmer_depths.push_back(count);
        progress();
      }
    }
  }
  progbar.done();
  barrier();
}
*/
