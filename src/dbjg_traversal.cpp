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

enum class WalkStatus { RUNNING, DEADEND, FORK, CONFLICT, REPEAT };

struct WalkInfo {
  WalkStatus status;
  bool drop;
  int32_t len;
  int64_t sum_depths;
};

static void abort_walk(dist_object<WalkInfo> &walk_info, bool walk_drop, int32_t walk_len, int64_t sum_depths, intrank_t start_rank,
                       WalkStatus walk_status) {
  rpc_ff(start_rank,
         [](dist_object<WalkInfo> &walk_info, bool walk_drop, int32_t walk_len, int64_t sum_depths, WalkStatus walk_status) {
           if (walk_info->status != WalkStatus::RUNNING) DIE("walk status is already done");
           walk_info->status = walk_status;
           walk_info->drop = walk_drop;
           walk_info->len = walk_len;
           walk_info->sum_depths = sum_depths;
         }, walk_info, walk_drop, walk_len, sum_depths, walk_status);
}

static void traverse_step(dist_object<KmerDHT> &kmer_dht, const MerArray &merarr, char ext,
                          TraverseDirn dirn, dist_object<WalkInfo> &walk_info, int32_t walk_len, int64_t sum_depths,
                          intrank_t start_rank, global_ptr<char> uutig, bool revisit_allowed, char prev_ext, int32_t start_walk_us) {
  Kmer kmer(merarr);
  auto kmer_rc = kmer.revcomp();
  auto rc = false;
  KmerCounts *kmer_counts = nullptr;
  if (kmer_rc < kmer) {
    rc = true;
    kmer_counts = kmer_dht->get_local_kmer_counts(kmer_rc);
  } else {
    kmer_counts = kmer_dht->get_local_kmer_counts(kmer);
  }
  // this kmer doesn't exist, abort
  if (!kmer_counts) {
    DBG_TRAVERSE("TERM: could not find kmer ", kmer, " with extension ", ext, " or revcomp kmer ", kmer_rc, "\n");
    abort_walk(walk_info, false, walk_len, sum_depths, start_rank, WalkStatus::DEADEND);
    return;
  }
  char left = kmer_counts->left;
  char right = kmer_counts->right;
  if (left == 'X' || right == 'X') {
    abort_walk(walk_info, false, walk_len, sum_depths, start_rank, WalkStatus::DEADEND);
    return;
  }    
  if (left == 'F' || right == 'F') {
    abort_walk(walk_info, false, walk_len, sum_depths, start_rank, WalkStatus::FORK);
    return;
  }    
  if (rc) {
    DBG_TRAVERSE("using rc (", kmer_rc, ")\n");
    left = comp_nucleotide(left);
    right = comp_nucleotide(right);
    swap(left, right);
  }
  // check for conflicts
  if (prev_ext && ((dirn == TraverseDirn::LEFT && prev_ext != right) || (dirn == TraverseDirn::RIGHT && prev_ext != left))) {
    DBG_TRAVERSE("TERM: conflicting kmers ", kmer, ", ext ", left, " ", right, " prev ", prev_ext, "\n");
    //WARN("TERM: conflicting kmers ", kmer, ", ext ", left, " ", right, " prev ", prev_ext);
    abort_walk(walk_info, false, walk_len, sum_depths, start_rank, WalkStatus::CONFLICT);
    return;
  }
  
  DBG_TRAVERSE("kmer ", kmer, " with ext ", ext, " left ", left, " right ", right, "\n");
  // if visited by another rank first
  if (kmer_counts->visited >= 0 && kmer_counts->visited != start_rank) {
    // resolve in favor of earliest starting time
    if (kmer_counts->start_walk_us < start_walk_us || (kmer_counts->start_walk_us == start_walk_us && kmer_counts->visited > start_rank)) {
      DBG_TRAVERSE("TERM: kmer ", kmer, " already visited ", kmer_counts->visited, " start walk ", kmer_counts->start_walk_us,
                   " new start walk ", start_walk_us, " start rank ", start_rank, "\n");
      abort_walk(walk_info, true, walk_len, sum_depths, start_rank, WalkStatus::CONFLICT);
      return;
    }
  }
  // already visited by start rank, repeat, abort (but allowed if traversing right after adding start kmer previously)
  if (kmer_counts->visited == start_rank && !revisit_allowed) {
    DBG_TRAVERSE("TERM: kmer ", kmer, " already visited ", kmer_counts->visited, " start rank ", start_rank, "\n");
    abort_walk(walk_info, false, walk_len, sum_depths, start_rank, WalkStatus::REPEAT);
    return;
  }
  // mark as visited
  kmer_counts->visited = start_rank;
  // set the walk time
  kmer_counts->start_walk_us = start_walk_us;
  // add the last nucleotide to the uutig
  DBG_TRAVERSE("wrote ext ", ext, " at position ", walk_len, "\n");
  walk_len++;
  sum_depths += kmer_counts->count;
  // now attempt to walk to next kmer
  //Kmer next_kmer();
  Kmer next_kmer;
  string kmer_str = kmer.to_string();
  // get next extension
  char next_ext;
  if (dirn == TraverseDirn::LEFT) {
    next_ext = left;
    prev_ext = kmer_str.back();
    next_kmer = kmer.backward_base(next_ext);
  } else {
    next_ext = right;
    prev_ext = kmer_str.front();
    next_kmer = kmer.forward_base(next_ext);
  }
  if (next_ext == 'X' || next_ext == 'F') DIE("Found X or F");
  auto next_kmer_rc = next_kmer.revcomp();
  auto target_rank = (next_kmer_rc < next_kmer ? kmer_dht->get_kmer_target_rank(next_kmer_rc) :
                      kmer_dht->get_kmer_target_rank(next_kmer));
  // valid new kmer can be constructed
  DBG_TRAVERSE("next kmer ", next_kmer, "\n");
  rput(ext, uutig + walk_len - 1).then(
    [=, &kmer_dht, &walk_info]() {
      // only do the rpc to the next kmer vertex once the rput has completed
      rpc_ff(target_rank, traverse_step, kmer_dht, next_kmer.get_array(), next_ext, dirn, walk_info, walk_len,
             sum_depths, start_rank, uutig, false, prev_ext, start_walk_us);
    });
}

static bool traverse_left(dist_object<KmerDHT> &kmer_dht, Kmer &kmer, global_ptr<char> &uutig_gptr, 
                          string &uutig_str, dist_object<WalkInfo> &walk_info, int32_t start_walk_us) {
  walk_info->status = WalkStatus::RUNNING;
  walk_info->drop = false;
  walk_info->len = 0;
  walk_info->sum_depths = 0;
  int walk_len = 0;
  int64_t sum_depths = 0;
  char prev_ext = 0;
  string kmer_str = kmer.to_string();
  traverse_step(kmer_dht, kmer.get_array(), kmer_str.front(), TraverseDirn::LEFT, walk_info, walk_len, sum_depths,
                rank_me(), uutig_gptr, false, prev_ext, start_walk_us);
  while (walk_info->status == WalkStatus::RUNNING) progress();
  if (walk_info->drop) return false;
  char *local_uutig = uutig_gptr.local();
  // reverse it because we were walking backwards
  for (int i = walk_info->len - 1; i >= 0; i--) uutig_str += local_uutig[i];
  return true;
}

static bool traverse_right(dist_object<KmerDHT> &kmer_dht, Kmer &kmer, global_ptr<char> &uutig_gptr, 
                           string &uutig_str, dist_object<WalkInfo> &walk_info, int32_t start_walk_us) {
  walk_info->status = WalkStatus::RUNNING;
  walk_info->drop = false;
  walk_info->len = 0;
  int walk_len = 0;
  int64_t sum_depths = 0;
  string kmer_str = kmer.to_string();
  char prev_ext = 0;
  traverse_step(kmer_dht, kmer.get_array(), kmer_str.back(), TraverseDirn::RIGHT, walk_info, walk_len, sum_depths,
                rank_me(), uutig_gptr, true, prev_ext, start_walk_us);
  while (walk_info->status == WalkStatus::RUNNING) progress(); 
  if (walk_info->drop) return false;
  char *local_uutig = uutig_gptr.local();
  local_uutig[walk_info->len] = 0;
  uutig_str += local_uutig;
  return true;
}

void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs) {
  Timer timer(__func__, true);
  // allocate space for biggest possible uutig in global storage
  const int MAX_UUTIG_LEN = 10000000;
  global_ptr<char> uutig_gptr = new_array<char>(MAX_UUTIG_LEN);
  dist_object<WalkInfo> walk_info({WalkStatus::RUNNING, false, 0});
  int64_t num_drops = 0, num_walks = 0, num_deadends = 0, num_forks = 0, num_conflicts = 0, num_repeats = 0;
  Contigs uutigs;
  barrier();
  {
    ProgressBar progbar(kmer_dht->get_local_num_kmers(), "Traversing deBruijn graph");
    for (auto it = kmer_dht->local_kmers_begin(); it != kmer_dht->local_kmers_end(); it++) {
      progbar.update();
      auto kmer = it->first;
      auto kmer_counts = &it->second;
      // don't start any new walk if this kmer has already been visited
      if (kmer_counts->visited > -1) continue;
      // don't start walks on kmers without extensions on both sides
      if (kmer_counts->left == 'X' || kmer_counts->left == 'F' || kmer_counts->right == 'X' || kmer_counts->right == 'F') continue;
      int32_t start_walk_us = kmer_dht->get_time_offset_us();
      // walk left first
      string uutig_str = "";
      if (!traverse_left(kmer_dht, kmer, uutig_gptr, uutig_str, walk_info, start_walk_us)) {
        num_drops++;
        continue;
      }
      auto sum_depths = walk_info->sum_depths;
      auto kmer_str = kmer.to_string();
      uutig_str += kmer_str.substr(1, kmer_len - 2);
      if (!traverse_right(kmer_dht, kmer, uutig_gptr, uutig_str, walk_info, start_walk_us)) {
        num_drops++;
        continue;
      }
      if (walk_info->status == WalkStatus::DEADEND) num_deadends++;
      if (walk_info->status == WalkStatus::FORK) num_forks++;
      if (walk_info->status == WalkStatus::CONFLICT) num_conflicts++;
      if (walk_info->status == WalkStatus::REPEAT) num_repeats++;
      sum_depths += walk_info->sum_depths;
      Contig contig = {0, uutig_str, (double)sum_depths / (uutig_str.length() - kmer_len + 2)};
      uutigs.add_contig(contig);
      num_walks++;
    }
    progbar.done();
    // FIXME: now steal kmers from others???
  }
  delete_array(uutig_gptr);
  barrier();
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
  Timer timer(__func__);
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
