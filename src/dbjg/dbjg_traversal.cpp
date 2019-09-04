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
#include "options_kcount.hpp"
#include "progressbar.hpp"
#include "kmer_dht.hpp"
#include "dbjg_traversal.hpp"

//#define DBG_TRAVERSE DBG
#define DBG_TRAVERSE(...)

using namespace std;
using namespace upcxx;

extern ofstream _dbgstream;

enum class TraverseDirn { LEFT, RIGHT };

struct WalkStatus {
  bool done;
  bool drop;
  int32_t len;
  int64_t sum_depths;
};

static void abort_walk(dist_object<WalkStatus> &walk_status, bool walk_drop, int32_t walk_len, int64_t sum_depths, intrank_t start_rank)
{
  rpc_ff(start_rank,
         [](dist_object<WalkStatus> &walk_status, bool walk_drop, int32_t walk_len, int64_t sum_depths) {
           if (walk_status->done) DIE("walk status is already done");
           walk_status->done = true;
           walk_status->drop = walk_drop;
           walk_status->len = walk_len;
           walk_status->sum_depths = sum_depths;
         }, walk_status, walk_drop, walk_len, sum_depths);
}

static int num_conflicts = 0;

static void traverse_step(dist_object<KmerDHT> &kmer_dht, const Kmer::MERARR &merarr, char ext, TraverseDirn dirn,
                          dist_object<WalkStatus> &walk_status, int32_t walk_len, int64_t sum_depths, intrank_t start_rank,
                          global_ptr<char> uutig, bool revisit_allowed, char prev_ext, int32_t start_walk_us)
{
  Kmer kmer(merarr);
  auto kmer_rc = kmer.twin();
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
    abort_walk(walk_status, false, walk_len, sum_depths, start_rank);
    return;
  }
  char left = kmer_counts->left;
  char right = kmer_counts->right;
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
    abort_walk(walk_status, false, walk_len, sum_depths, start_rank);
    num_conflicts++;
    return;
  }
  
  DBG_TRAVERSE("kmer ", kmer, " with ext ", ext, " left ", left, " right ", right, "\n");
  // if visited by another rank first
  if (kmer_counts->visited >= 0 && kmer_counts->visited != start_rank) {
    // resolve in favor of earliest starting time
    if (kmer_counts->start_walk_us < start_walk_us || (kmer_counts->start_walk_us == start_walk_us && kmer_counts->visited > start_rank)) {
      DBG_TRAVERSE("TERM: kmer ", kmer, " already visited ", kmer_counts->visited, " start walk ", kmer_counts->start_walk_us,
                   " new start walk ", start_walk_us, " start rank ", start_rank, "\n");
      abort_walk(walk_status, true, walk_len, sum_depths, start_rank);
      return;
    }
  }
  // already visited by start rank, circular, abort (but allowed if traversing right after adding start kmer previously)
  if (kmer_counts->visited == start_rank && !revisit_allowed) {
    DBG_TRAVERSE("TERM: kmer ", kmer, " already visited ", kmer_counts->visited, " start rank ", start_rank, "\n");
    abort_walk(walk_status, false, walk_len, sum_depths, start_rank);
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
  Kmer next_kmer;
  string kmer_str = kmer.toString();
  // get next extension
  char next_ext;
  if (dirn == TraverseDirn::LEFT) {
    next_ext = left;
    prev_ext = kmer_str.back();
    next_kmer = kmer.backwardBase(next_ext);
  } else {
    next_ext = right;
    prev_ext = kmer_str.front();
    next_kmer = kmer.forwardBase(next_ext);
  }
  if (next_ext == 'X' || next_ext == 'F') DIE("Found X or F");
  auto next_kmer_rc = next_kmer.twin();
  auto target_rank = (next_kmer_rc < next_kmer ? kmer_dht->get_kmer_target_rank(next_kmer_rc) :
                      kmer_dht->get_kmer_target_rank(next_kmer));
  // valid new kmer can be constructed
  DBG_TRAVERSE("next kmer ", next_kmer, "\n");
  rput(ext, uutig + walk_len - 1).then(
    [=, &kmer_dht, &walk_status]() {
      // only do the rpc to the next kmer vertex once the rput has completed
      rpc_ff(target_rank, traverse_step, kmer_dht, next_kmer.getArray(), next_ext, dirn, walk_status, walk_len, sum_depths, start_rank,
             uutig, false, prev_ext, start_walk_us);
    });
}

bool traverse_left(int kmer_len, dist_object<KmerDHT> &kmer_dht, Kmer &kmer, global_ptr<char> &uutig_gptr, 
                   string &uutig_str, dist_object<WalkStatus> &walk_status, int32_t start_walk_us)
{
  walk_status->done = false;
  walk_status->drop = false;
  walk_status->len = 0;
  walk_status->sum_depths = 0;
  int walk_len = 0;
  int64_t sum_depths = 0;
  char prev_ext = 0;
  string kmer_str = kmer.toString();
  traverse_step(kmer_dht, kmer.getArray(), kmer_str.front(), TraverseDirn::LEFT, walk_status, walk_len, sum_depths, rank_me(),
                uutig_gptr, false, prev_ext, start_walk_us);
  while (!walk_status->done) progress();
  if (walk_status->drop) return false;
  char *local_uutig = uutig_gptr.local();
  // reverse it because we were walking backwards
  for (int i = walk_status->len - 1; i >= 0; i--) uutig_str += local_uutig[i];
  return true;
}

bool traverse_right(int kmer_len, dist_object<KmerDHT> &kmer_dht, Kmer &kmer, global_ptr<char> &uutig_gptr, 
                    string &uutig_str, dist_object<WalkStatus> &walk_status, int32_t start_walk_us)
{
  walk_status->done = false;
  walk_status->drop = false;
  walk_status->len = 0;
  int walk_len = 0;
  int64_t sum_depths = 0;
  string kmer_str = kmer.toString();
  char prev_ext = 0;
  traverse_step(kmer_dht, kmer.getArray(), kmer_str.back(), TraverseDirn::RIGHT, walk_status, walk_len, sum_depths, rank_me(),
                uutig_gptr, true, prev_ext, start_walk_us);
  while (!walk_status->done) progress(); 
  if (walk_status->drop) return false;
  char *local_uutig = uutig_gptr.local();
  local_uutig[walk_status->len] = 0;
  uutig_str += local_uutig;
  return true;
}

void traverse_debruijn_graph(shared_ptr<OptionsKcount> options, dist_object<KmerDHT> &kmer_dht)
{
  Timer timer(__func__);
  // allocate space for biggest possible uutig in global storage
  const int MAX_UUTIG_LEN = 10000000;
  global_ptr<char> uutig_gptr = new_array<char>(MAX_UUTIG_LEN);
  dist_object<WalkStatus> walk_status({false, false, 0});
  int64_t num_drops = 0;
  int64_t num_walks = 0;
  vector<string> uutigs;
  vector<double> depths;
  num_conflicts = 0;
  barrier();
  {
    ProgressBar progbar(kmer_dht->get_local_num_kmers(), "Traversing deBruijn graph");
    for (auto it = kmer_dht->local_kmers_begin(); it != kmer_dht->local_kmers_end(); it++) {
      progbar.update();
      auto kmer = it->first;
      auto kmer_counts = &it->second;
      // don't start any new walk if this kmer has already been visited
      if (kmer_counts->visited > -1) continue;
      int32_t start_walk_us = kmer_dht->get_time_offset_us();
      // walk left first
      string uutig_str = "";
      if (!traverse_left(options->kmer_len, kmer_dht, kmer, uutig_gptr, uutig_str, walk_status, start_walk_us)) {
        num_drops++;
        continue;
      }
      auto sum_depths = walk_status->sum_depths;
      auto kmer_str = kmer.toString();
      uutig_str += kmer_str.substr(1, options->kmer_len - 2);
      if (!traverse_right(options->kmer_len, kmer_dht, kmer, uutig_gptr, uutig_str, walk_status, start_walk_us)) {
        num_drops++;
        continue;
      }
      sum_depths += walk_status->sum_depths;
      depths.push_back((double)sum_depths / (uutig_str.length() - options->kmer_len + 2));
      uutigs.push_back(uutig_str);
      num_walks++;
    }
    progbar.done();
    // now steal kmers from others
    
  }
  delete_array(uutig_gptr);
  barrier();
  auto all_num_conflicts = reduce_one(num_conflicts, op_fast_add, 0).wait();
  SOUT("Found ", perc_str(all_num_conflicts, kmer_dht->get_num_kmers()), " conflicting kmers\n");
  // count number of uutigs found by this rank, and save a vector of pointers to them
  int64_t num_uutigs = 0;
  vector<string*> my_uutigs;
  vector<double> my_depths;
  {
    ProgressBar progbar(uutigs.size(), "Dropping duplicate walks");
    int i = 0;
    for (auto &uutig : uutigs) {
      progbar.update();
      if (uutig.length() < options->kmer_len) DIE("uutig length ", uutig.length(), " less than kmer len");
      // now check to make sure we're the owner of this one - this is after all ranks have finished traversals
      Kmer start_kmer(uutig.substr(0, options->kmer_len).c_str());
      auto start_kmer_rc = start_kmer.twin();
      if (start_kmer_rc < start_kmer) start_kmer = start_kmer_rc;
      auto visited = kmer_dht->get_visited(start_kmer);
      if (visited == -2) DIE("Could not find kmer ", start_kmer, "\n");
      if (visited == -1) DIE("unitialized visited for kmer ", start_kmer, "\n");
      if (visited != rank_me()) {
        num_drops++;
        DBG_TRAVERSE("Drop kmer ", start_kmer, " visited is ", visited, " not my rank ", rank_me(), "\n");
      } else {
        num_uutigs++;
        my_uutigs.push_back(&uutig);
        my_depths.push_back(depths[i]);
      }
      i++;
    }
    progbar.done();
  }
  barrier();
  DBG_TRAVERSE("dropped ", num_drops, " out of ", num_walks, " walks\n");
  auto all_num_drops = reduce_one(num_drops, op_fast_add, 0).wait();
  auto all_num_walks = reduce_one(num_walks, op_fast_add, 0).wait();
  SOUT("Dropped ", perc_str(all_num_drops, all_num_walks), " duplicate traversals\n");
  // now get unique ids for the uutigs
  atomic_domain<size_t> ad({atomic_op::fetch_add, atomic_op::load});
  global_ptr<size_t> counter = nullptr;
  if (!rank_me()) counter = new_<size_t>(0);
  counter = broadcast(counter, 0).wait();
  size_t my_counter = ad.fetch_add(counter, num_uutigs, memory_order_relaxed).wait();
  // wait until all ranks have updated the global counter
  barrier();
  ad.destroy();
  {
    string uutigs_fname = string(options->cached_io ? "/dev/shm/" : "./");
    //uutigs_fname += "UUtigs-" + to_string(options->kmer_len) + ".fasta.gz";
    uutigs_fname += "UUtigs_contigs-" + to_string(options->kmer_len) + ".fasta.gz";
    get_rank_path(uutigs_fname, rank_me());
    zstr::ofstream uutigs_file(uutigs_fname);
    ostringstream uutigs_out_buf;
    
    string depths_fname = string(options->cached_io ? "/dev/shm/" : "./");
    //depths_fname += "UUtigs-" + to_string(options->kmer_len) + ".depths.gz";
    depths_fname += "merDepth_UUtigs_contigs-" + to_string(options->kmer_len) + ".txt.gz";
    get_rank_path(depths_fname, rank_me());
    zstr::ofstream depths_file(depths_fname);
    ostringstream depths_out_buf;
    
    ProgressBar progbar(my_uutigs.size(), "Writing uutigs");
    int64_t cid = my_counter;
    for (auto uutig : my_uutigs) {
      uutigs_out_buf << ">Contig_" << cid << " " << (uutig->length() - options->kmer_len + 1) << " " << my_depths[cid - my_counter] << endl;
      string rc_uutig = revcomp(*uutig);
      if (rc_uutig < *uutig) *uutig = rc_uutig;
      // fold output
      for (int p = 0; p < uutig->length(); p += 50)
        uutigs_out_buf << uutig->substr(p, 50) << endl;
      depths_out_buf << "Contig" << cid << "\t" << (uutig->length() - options->kmer_len + 1) << "\t" << my_depths[cid - my_counter] << endl;
      progbar.update();
      cid++;
      if (!(cid % 1000)) {
        uutigs_file << uutigs_out_buf.str();
        uutigs_out_buf = ostringstream();
        depths_file << depths_out_buf.str();
        depths_out_buf = ostringstream();
      }
    }
    if (!uutigs_out_buf.str().empty()) uutigs_file << uutigs_out_buf.str();
    uutigs_file.close();
    if (!depths_out_buf.str().empty()) depths_file << depths_out_buf.str();
    depths_file.close();
    progbar.done();
    SOUT("Wrote ", reduce_one(my_uutigs.size(), op_fast_add, 0).wait(), " uutigs to ", uutigs_fname, "\n");

    // print some metafiles used in downstream steps
    auto tot_num_kmers = kmer_dht->get_num_kmers();
    if (!rank_me()) {
      string count_fname = (options->cached_io ? "/dev/shm/" : "./");
      count_fname += "nUUtigs_contigs-" + to_string(options->kmer_len) + ".txt";
      get_rank_path(count_fname, -1);
      ofstream count_file(count_fname);
      count_file << tot_num_kmers << endl;
    }
    {
      string my_count_fname = (options->cached_io ? "/dev/shm/" : "./");
      my_count_fname += "myUUtigs_contigs-" + to_string(options->kmer_len) + ".txt";
      get_rank_path(my_count_fname, rank_me());
      ofstream my_count_file(my_count_fname);
      my_count_file << kmer_dht->get_local_num_kmers() << endl;
    }
  }
  barrier();
}

