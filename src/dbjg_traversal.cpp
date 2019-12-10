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

enum class WalkStatus { DEADEND = 'X', FORK = 'F', CONFLICT = 'O', REPEAT = 'R'};

struct WalkInfo {
  bool drop;
  int32_t len;
  int64_t sum_depths;
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

static void abort_walk(intrank_t target_rank, global_ptr<char> uutig, WalkStatus walk_status, WalkInfo walk_info) {
  // we use a rpc here to ensure synchronization
  rpc_ff(target_rank,
         [](global_ptr<char> uutig, WalkStatus walk_status, WalkInfo walk_info){
           char *local_uutig = uutig.local();
           local_uutig[0] = (char)walk_status;
           memcpy(local_uutig + 1, (char*)&walk_info, sizeof(walk_info));
         }, uutig, walk_status, walk_info);
}

static void traverse_step(dist_object<KmerDHT> &kmer_dht, const MerArray &merarr, char ext, TraverseDirn dirn, 
                          int32_t walk_len, int64_t sum_depths, intrank_t start_rank, global_ptr<char> uutig, 
                          bool revisit_allowed, char prev_ext, int32_t start_walk_us) {
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
    abort_walk(start_rank, uutig, WalkStatus::DEADEND, {false, walk_len, sum_depths});
    return;
  }
  char left = kmer_counts->left;
  char right = kmer_counts->right;
  if (left == 'X' || right == 'X') {
    DBG_TRAVERSE("TERM: deadend for ", kmer, " left is ", left, " and right is ", right, "\n");
    abort_walk(start_rank, uutig, WalkStatus::DEADEND, {false, walk_len, sum_depths}); 
    return;
  }    
  if (left == 'F' || right == 'F') {
    DBG_TRAVERSE("TERM: fork for ", kmer, " left is ", left, " and right is ", right, "\n");
    abort_walk(start_rank, uutig, WalkStatus::FORK, {false, walk_len, sum_depths});
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
    abort_walk(start_rank, uutig, WalkStatus::CONFLICT, {false, walk_len, sum_depths});
    return;
  }
  
  DBG_TRAVERSE("kmer ", kmer, " with ext ", ext, " left ", left, " right ", right, "\n");
  // if visited by another rank first
  if (kmer_counts->visited >= 0 && (kmer_counts->visited != start_rank || kmer_counts->start_walk_us != start_walk_us)) {
    // resolve in favor of earliest starting time
    if (kmer_counts->start_walk_us < start_walk_us || 
        (kmer_counts->start_walk_us == start_walk_us && kmer_counts->visited > start_rank)) {
      DBG_TRAVERSE("TERM: kmer ", kmer, " already visited ", kmer_counts->visited, " start walk ", kmer_counts->start_walk_us,
                   " new start walk ", start_walk_us, " start rank ", start_rank, "\n");
      abort_walk(start_rank, uutig, WalkStatus::CONFLICT, {true, walk_len, sum_depths});
      return;
    }
  }
  // already visited by start rank, repeat, abort (but allowed if traversing right after adding start kmer previously)
  if (kmer_counts->visited == start_rank && kmer_counts->start_walk_us == start_walk_us && !revisit_allowed) {
    DBG_TRAVERSE("TERM: kmer ", kmer, " already visited ", kmer_counts->visited, " start rank ", start_rank, "\n");
    abort_walk(start_rank, uutig, WalkStatus::REPEAT, {false, walk_len, sum_depths});
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
  if (walk_len >= MAX_UUTIG_BUF_LEN - 2 - sizeof(WalkInfo)) DIE("walk is greater than buf size ", MAX_UUTIG_BUF_LEN);
  rput(ext, uutig + 1 + sizeof(WalkInfo) + walk_len - 1).then(
    [=, &kmer_dht]() {
      // only do the rpc to the next kmer vertex once the rput has completed
      rpc_ff(target_rank, traverse_step, kmer_dht, next_kmer.get_array(), next_ext, dirn, walk_len,
             sum_depths, start_rank, uutig, false, prev_ext, start_walk_us);
    });
}

static bool traverse_start(dist_object<KmerDHT> &kmer_dht, Kmer &kmer, int32_t start_walk_us, 
                           bool revisit_allowed, TraverseDirn dirn, global_ptr<char> uutig_gptr, char next_ext) {
  char *local_uutig = uutig_gptr.local();
  local_uutig[0] = 0;
  int walk_len = 0;
  int64_t sum_depths = 0;
  char prev_ext = 0;
  traverse_step(kmer_dht, kmer.get_array(), next_ext, dirn, walk_len, sum_depths, rank_me(), uutig_gptr,
                revisit_allowed, prev_ext, start_walk_us);
}

static bool traverse_complete(string &uutig_str, WalkTermStats &walk_term_stats, global_ptr<char> uutig_gptr, 
                              TraverseDirn dirn, int64_t &sum_depths) {
  char *local_uutig = uutig_gptr.local();
  while (local_uutig[0] == 0) progress();
  assert(local_uutig[0] == 'X' || local_uutig[0] == 'F' || local_uutig[0] == 'O' || local_uutig[0] == 'R');
  WalkStatus walk_status = (WalkStatus)local_uutig[0];
  WalkInfo walk_info({false, 0, 0});
  memcpy((char*)&walk_info, local_uutig + 1, sizeof(WalkInfo));
  if (!walk_info.drop) {
    walk_term_stats.update(walk_status);
    char *uutig = local_uutig + 1 + sizeof(WalkInfo);
    uutig[walk_info.len] = 0;
    sum_depths += walk_info.sum_depths;
#ifdef DEBUG
    for (int i = 0; i < walk_info.len; i++) {
      if (uutig[i] != 'A' && uutig[i] != 'C' && uutig[i] != 'G' && uutig[i] != 'T' && uutig[i] != 'N')
        DIE("left: bad uutig character '", uutig[i], "' (", (int)uutig[i], ") ", uutig, "\n");
    }
#endif
    if (dirn == TraverseDirn::LEFT) {
      // reverse it because we were walking backwards
      for (int i = walk_info.len - 1; i >= 0; i--) uutig_str += uutig[i];
    } else {
      uutig_str += uutig;
    }
  }
  delete_array(uutig_gptr);
  return (!walk_info.drop);
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
    deque<global_ptr<char>> uutig_gptrs;
    deque<string> kmer_strs;
    int64_t num_kmers = 0;
    for (auto it = kmer_dht->local_kmers_begin(); it != kmer_dht->local_kmers_end(); it++) {
      progbar.update();
      progress();
      auto kmer = it->first;
      auto kmer_counts = &it->second;
      // don't start any new walk if this kmer has already been visited
      if (kmer_counts->visited > -1) continue;
      // don't start walks on kmers without extensions on both sides
      if (kmer_counts->left == 'X' || kmer_counts->left == 'F' || kmer_counts->right == 'X' || kmer_counts->right == 'F') 
        continue;
      int32_t start_walk_us = kmer_dht->get_time_offset_us();
      auto kmer_str = kmer.to_string();
      // walk left first
      kmer_strs.push_back(kmer_str);
      global_ptr<char> left_uutig_gptr = new_array<char>(MAX_UUTIG_BUF_LEN);
      uutig_gptrs.push_back(left_uutig_gptr);
      traverse_start(kmer_dht, kmer, start_walk_us, false, TraverseDirn::LEFT, left_uutig_gptr, kmer_str.front());
      global_ptr<char> right_uutig_gptr = new_array<char>(MAX_UUTIG_BUF_LEN);
      uutig_gptrs.push_back(right_uutig_gptr);
      traverse_start(kmer_dht, kmer, start_walk_us, true, TraverseDirn::RIGHT, right_uutig_gptr, kmer_str.back());
      num_kmers++;
      // limit the number of walks in flight
      while (uutig_gptrs.size() > MAX_DBJG_WALKS_IN_FLIGHT) {
        auto kmer_str = kmer_strs.front();
        kmer_strs.pop_front();
        auto uutig_gptr_left = uutig_gptrs.front();
        uutig_gptrs.pop_front();
        auto uutig_gptr_right = uutig_gptrs.front();
        uutig_gptrs.pop_front();
        string uutig_str = "";
        int64_t sum_depths = 0;
        if (!traverse_complete(uutig_str, walk_term_stats, uutig_gptr_left, TraverseDirn::LEFT, sum_depths)) {
          num_drops++;
          continue;
        }
        uutig_str += kmer_str.substr(1, kmer_len - 2);
        if (!traverse_complete(uutig_str, walk_term_stats, uutig_gptr_right, TraverseDirn::RIGHT, sum_depths)) {
          num_drops++;
          continue;
        }
        Contig contig = {0, uutig_str, (double)sum_depths / (uutig_str.length() - kmer_len + 2)};
        uutigs.add_contig(contig);
        num_walks++;
      }
    }
    for (int i = 0; i < uutig_gptrs.size(); i += 2) {
      auto &kmer_str = kmer_strs[i / 2];
      string uutig_str = "";
      int64_t sum_depths = 0;
      if (!traverse_complete(uutig_str, walk_term_stats, uutig_gptrs[i], TraverseDirn::LEFT, sum_depths)) {
        num_drops++;
        continue;
      }
      uutig_str += kmer_str.substr(1, kmer_len - 2);
      if (!traverse_complete(uutig_str, walk_term_stats, uutig_gptrs[i + 1], TraverseDirn::RIGHT, sum_depths)) {
        num_drops++;
        continue;
      }
      Contig contig = {0, uutig_str, (double)sum_depths / (uutig_str.length() - kmer_len + 2)};
      uutigs.add_contig(contig);
      num_walks++;
    }
    uutig_gptrs.clear();
    kmer_strs.clear();
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
