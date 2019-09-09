#ifndef __KMER_DHT
#define __KMER_DHT

#include <limits>
#include <iostream>
#include <map>
#include <fstream>
#include <chrono>
#include <stdarg.h>
#include <upcxx/upcxx.hpp>

#include "bytell_hash_map.hpp"
#include "progressbar.hpp"
#include "utils.hpp"
#include "kmer.hpp"
#include "bloom.hpp"
#include "zstr.hpp"
#include "histogram_dht.hpp"
#include "aggr_store.hpp"

using std::vector;
using std::pair;
using std::ostream;
using std::ostringstream;
using std::sort;
using std::numeric_limits;
using std::make_shared;
using std::make_pair;
using std::shared_ptr;
using std::swap;
using std::array;
using std::endl;

using upcxx::intrank_t;
using upcxx::rank_me;
using upcxx::rank_n;
using upcxx::barrier;
using upcxx::dist_object;
using upcxx::reduce_one;
using upcxx::reduce_all;
using upcxx::op_fast_add;
using upcxx::op_fast_max;
using upcxx::progress;
using upcxx::global_ptr;
using upcxx::new_array;
using upcxx::delete_array;

#define USE_BYTELL

#define DBG_ON

//#define DBG_INS_CTG_KMER DBG
#define DBG_INS_CTG_KMER(...)


#ifndef BLOOM_FP
#define BLOOM_FP 0.05
#endif

enum PASS_TYPE { BLOOM_SET_PASS, BLOOM_COUNT_PASS, NO_BLOOM_PASS, CTG_BLOOM_SET_PASS, CTG_KMERS_PASS };

using ext_count_t = uint16_t;
    
struct ExtCounts {
  ext_count_t count_A;
  ext_count_t count_C;
  ext_count_t count_G;
  ext_count_t count_T;

  array<pair<char, int>, 4> get_sorted() {
    array<pair<char, int>, 4> counts = {make_pair('A', (int)count_A), make_pair('C', (int)count_C),
                                        make_pair('G', (int)count_G), make_pair('T', (int)count_T)};
    sort(std::begin(counts), std::end(counts),
         [](const auto &elem1, const auto &elem2) {
           if (elem1.second == elem2.second) return elem1.first > elem2.first;
           else return elem1.second > elem2.second;
         });
    return counts;
  }

  bool is_zero() {
    if (count_A + count_C + count_G + count_T == 0) return true;
    return false;
  }
};

// total bytes: 2+8+8=18
struct KmerCounts {
  // how many times this kmer has occurred: don't need to count beyond 65536
  // count of high quality forward and backward exts
  ExtCounts left_exts;
  ExtCounts right_exts;
  // the final extensions chosen - A,C,G,T, or F,X
  char left;
  char right;
  uint16_t count;
  bool from_ctg;
  int32_t visited;
  int32_t start_walk_us;

  pair<char, char> get_exts(int min_depth_cutoff, double dynamic_min_depth) {
    auto sorted_lefts = left_exts.get_sorted();
    auto sorted_rights = right_exts.get_sorted();
    char left = sorted_lefts[0].first;
    int leftmax = sorted_lefts[0].second;
    int leftmin = sorted_lefts[1].second;
    char right = sorted_rights[0].first;
    int rightmax = sorted_rights[0].second;
    int rightmin = sorted_rights[1].second;
    int dmin_dyn = (1.0 - dynamic_min_depth) * count;      // integer floor
    if (dmin_dyn < min_depth_cutoff) dmin_dyn = min_depth_cutoff;
    if (leftmax < dmin_dyn) left = 'X';
    else if (leftmin >= dmin_dyn) left = 'F';
    if (rightmax < dmin_dyn) right = 'X';
    else if (rightmin >= dmin_dyn) right = 'F';
    return {left, right};
  }
};
  
#ifdef USE_BYTELL
using kmer_map_t = ska::bytell_hash_map<Kmer, KmerCounts, KmerHash, KmerEqual>;
#else
using kmer_map_t = std::unordered_map<Kmer, KmerCounts, KmerHash, KmerEqual>;
#endif

class KmerDHT {
private:

  // total bytes for k = 51: 16+18+18=52
  struct MerarrAndExt {
    Kmer::MERARR merarr;
    char left_ext, right_ext;
    uint16_t count;
  };

  dist_object<kmer_map_t> kmers;
  AggrStore<Kmer::MERARR> kmer_store_bloom;
  AggrStore<MerarrAndExt> kmer_store;
  // The first bloom filter stores all kmers and is used to check for single occurrences to filter out
  dist_object<BloomFilter> bloom_filter1;
  // the second bloom filer stores only kmers that are above the repeat depth, and is used for correctly sizing the kmer hash table
  dist_object<BloomFilter> bloom_filter2;
  dist_object<int> min_depth_cutoff;
  dist_object<double> dynamic_min_depth;
  int64_t max_kmer_store_bytes;
  int64_t initial_kmer_dht_reservation;
  int64_t bloom1_cardinality;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_t;
  
public:
  int64_t num_prev_mers_from_ctgs;
private:

  struct InsertKmer {
    void operator()(MerarrAndExt &merarr_and_ext, dist_object<kmer_map_t> &kmers) {
      Kmer new_kmer(merarr_and_ext.merarr);
      // find it - if it isn't found then insert it, otherwise increment the counts
      const auto it = kmers->find(new_kmer);
      if (it == kmers->end()) {
        KmerCounts kmer_counts = { .left_exts = {0}, .right_exts = {0}, .left = 'X', .right = 'X', .count = merarr_and_ext.count,
                                   .from_ctg = false, .visited = -1, .start_walk_us = 0 };
        inc_ext(kmer_counts.left_exts, merarr_and_ext.left_ext);
        inc_ext(kmer_counts.right_exts, merarr_and_ext.right_ext);
        auto prev_bucket_count = kmers->bucket_count();
        kmers->insert({new_kmer, kmer_counts});
        if (prev_bucket_count < kmers->bucket_count())
          SOUT("*** Hash table on rank 0 was resized from ", prev_bucket_count, " to ", kmers->bucket_count(), "***\n");
      } else {
        auto kmer = &it->second;
        int newCount = kmer->count + merarr_and_ext.count;
        kmer->count = (newCount < numeric_limits<uint16_t>::max()) ? newCount : numeric_limits<uint16_t>::max();
        inc_ext(kmer->left_exts, merarr_and_ext.left_ext);
        inc_ext(kmer->right_exts, merarr_and_ext.right_ext);
      }
    }
  };
  dist_object<InsertKmer> insert_kmer;

  struct BloomSet {
    void operator()(Kmer::MERARR &merarr, dist_object<BloomFilter> &bloom_filter1, dist_object<BloomFilter> &bloom_filter2) {
      // look for it in the first bloom filter - if not found, add it just to the first bloom filter
      // if found, add it to the second bloom filter
      Kmer new_kmer(merarr);
      if (!bloom_filter1->possibly_contains(new_kmer.getBytes(), new_kmer.getNumBytes())) 
        bloom_filter1->add(new_kmer.getBytes(), new_kmer.getNumBytes());
      else 
        bloom_filter2->add(new_kmer.getBytes(), new_kmer.getNumBytes());
    }
  };
  dist_object<BloomSet> bloom_set;

  struct CtgBloomSet {
    void operator()(Kmer::MERARR &merarr, dist_object<BloomFilter> &bloom_filter2) {
      // only add to bloom_filter2
      Kmer new_kmer(merarr);
      bloom_filter2->add(new_kmer.getBytes(), new_kmer.getNumBytes());
    }
  };
  dist_object<CtgBloomSet> ctg_bloom_set;

  struct BloomCount {
    void operator()(MerarrAndExt &merarr_and_ext, dist_object<kmer_map_t> &kmers, dist_object<BloomFilter> &bloom_filter) {
      Kmer new_kmer(merarr_and_ext.merarr);
      // if the kmer is not found in the bloom filter, skip it
      if (!bloom_filter->possibly_contains(new_kmer.getBytes(), new_kmer.getNumBytes())) return;
      // add or update the kmer count
      const auto it = kmers->find(new_kmer);
      if (it == kmers->end()) {
        KmerCounts kmer_counts = { .left_exts = {0}, .right_exts = {0}, .left = 'X', .right = 'X', .count = merarr_and_ext.count,
                                   .from_ctg = false, .visited = -1, .start_walk_us = 0 };
        inc_ext(kmer_counts.left_exts, merarr_and_ext.left_ext);
        inc_ext(kmer_counts.right_exts, merarr_and_ext.right_ext);
        auto prev_bucket_count = kmers->bucket_count();
        kmers->insert({new_kmer, kmer_counts});
        // this shouldn't happen 
        if (prev_bucket_count < kmers->bucket_count())
          WARN("Hash table on rank 0 was resized from ", prev_bucket_count, " to ", kmers->bucket_count(), "\n");
      } else {
        auto kmer = &it->second;
        int newCount = kmer->count + merarr_and_ext.count;
        kmer->count = (newCount < numeric_limits<uint16_t>::max()) ? newCount : numeric_limits<uint16_t>::max();
        inc_ext(kmer->left_exts, merarr_and_ext.left_ext);
        inc_ext(kmer->right_exts, merarr_and_ext.right_ext);
      }
    }
  };
  dist_object<BloomCount> bloom_count;

  struct InsertCtgKmer {
    void operator()(MerarrAndExt &merarr_and_ext, dist_object<kmer_map_t> &kmers, dist_object<int> &min_depth_cutoff,
                    dist_object<double> &dynamic_min_depth) {
      // insert a new kmer derived from the previous round's contigs
      Kmer new_kmer(merarr_and_ext.merarr);
      const auto it = kmers->find(new_kmer);
      bool insert = false;
      if (it == kmers->end()) {
        // if it isn't found then insert it
        insert = true;
        DBG_INS_CTG_KMER("new ", merarr_and_ext.count, " ", merarr_and_ext.left_ext, " ", merarr_and_ext.right_ext, "\n");
      } else {
        auto kmer_counts = &it->second;
        DBG_INS_CTG_KMER(new_kmer, " old/new ", kmer_counts->count, " ", merarr_and_ext.count, " ",
                         merarr_and_ext.left_ext, " ", merarr_and_ext.right_ext, " ",
                         "A", kmer_counts->left_exts.count_A, " C", kmer_counts->left_exts.count_C, " ", 
                         "G", kmer_counts->left_exts.count_G, " T", kmer_counts->left_exts.count_T, " ", 
                         "A", kmer_counts->right_exts.count_A, " C", kmer_counts->right_exts.count_C, " ", 
                         "G", kmer_counts->right_exts.count_G, " T", kmer_counts->right_exts.count_T, "\n");
        if (!kmer_counts->from_ctg) {
          // existing kmer is from a read, only replace with new contig kmer if the existing kmer is not UU
          auto exts = kmer_counts->get_exts(*min_depth_cutoff, *dynamic_min_depth);
          if (exts.first == 'X' || exts.first == 'F' || exts.second == 'X' || exts.second == 'F') {
            // non-UU, replace
            insert = true;
            // but keep the count from the read kmer
            //if (kmer_counts->count > *min_depth_cutoff) merarr_and_ext.count = kmer_counts->count;
            DBG_INS_CTG_KMER("replace non-UU read kmer\n");
          }
        } else {
          // existing kmer is from the previous round's contigs
          // if this was previously a conflict and the depth was set to zero, do nothing
          if (!kmer_counts->count) {
            DBG_INS_CTG_KMER("skip conflicted kmer, depth 0\n");
          } else {
            // if the two contig kmers disagree on exts, set to purge this one by setting the count to 0
            auto exts = kmer_counts->get_exts(*min_depth_cutoff, *dynamic_min_depth);
            if (exts.first != merarr_and_ext.left_ext || exts.second != merarr_and_ext.right_ext) {
              merarr_and_ext.count = 0;
              insert = true;
              DBG_INS_CTG_KMER("purge mismatch\n");
            } else {
              // we have multiple occurrences of the same kmer derived from different contigs or
              // parts of contigs - sum the depths
              int sum_counts = (int)merarr_and_ext.count + kmer_counts->count;
              if (sum_counts > numeric_limits<uint16_t>::max()) sum_counts = numeric_limits<uint16_t>::max();
              merarr_and_ext.count = sum_counts;
              insert = true;
              DBG_INS_CTG_KMER("increase count\n");
            }
          }
        }
      }
      if (insert) {
        uint16_t count = merarr_and_ext.count;
        KmerCounts kmer_counts = { .left_exts = {0}, .right_exts = {0}, .left = 'X', .right = 'X', .count = count,
                                   .from_ctg = true, .visited = -1, .start_walk_us = 0 };
        inc_ext(kmer_counts.left_exts, merarr_and_ext.left_ext, count);
        inc_ext(kmer_counts.right_exts, merarr_and_ext.right_ext, count);
        (*kmers)[new_kmer] = kmer_counts;
      }
    }
  };
  dist_object<InsertCtgKmer> insert_ctg_kmer;
  
  static void inc_ext(ExtCounts &ext_counts, char ext, int count=1) {
    switch (ext) {
      case 'A': 
        count += ext_counts.count_A;
        ext_counts.count_A = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
      case 'C': 
        count += ext_counts.count_C;
        ext_counts.count_C = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
      case 'G':
        count += ext_counts.count_G;
        ext_counts.count_G = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
      case 'T': 
        count += ext_counts.count_T;
        ext_counts.count_T = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
    }
  }

public:

  KmerDHT(uint64_t cardinality, int max_kmer_store_bytes, int min_depth_cutoff, double dynamic_min_depth,
          bool use_bloom) : kmers({}), min_depth_cutoff(min_depth_cutoff), dynamic_min_depth(dynamic_min_depth), bloom_filter1({}),
                            bloom_filter2({}), kmer_store({}), kmer_store_bloom({}),
                            insert_kmer({}), bloom_set({}), ctg_bloom_set({}), bloom_count({}), insert_ctg_kmer({}),
                            max_kmer_store_bytes(max_kmer_store_bytes), initial_kmer_dht_reservation(0), bloom1_cardinality(0),
                            num_prev_mers_from_ctgs(0) {
    if (use_bloom) {
        kmer_store_bloom.set_size(max_kmer_store_bytes);
    } else {
        kmer_store.set_size(max_kmer_store_bytes);
    }
#ifdef USE_BYTELL
    SOUT("Using bytell hash map\n");
#else
    SOUT("Using std::unordered_map\n");
#endif
    if (use_bloom) {
      // in this case we get an accurate estimate of the hash table size after the first bloom round, so the hash table space is
      // reserved then
      double init_mem_free = get_free_mem_gb();
      bloom_filter1->init(cardinality, BLOOM_FP);
      bloom_filter2->init(cardinality/4, BLOOM_FP); // second bloom will have far fewer entries - assume 75% are filtered out
      SOUT("Bloom filters used ", (init_mem_free - get_free_mem_gb()), "GB memory on node 0\n");
    } else {
      double init_mem_free = get_free_mem_gb();
      barrier();
      // in this case we have to roughly estimate the hash table size - we do the allocation here
      // rough estimate at 5x depth
      cardinality /= 5;
      initial_kmer_dht_reservation = cardinality;    
      kmers->reserve(cardinality);
      double kmers_space_reserved = cardinality * (sizeof(Kmer) + sizeof(KmerCounts));
      SOUT("Rank 0 is reserving ", get_size_str(kmers_space_reserved), " for kmer hash table with ", cardinality, " entries (",
           kmers->bucket_count(), " buckets)\n");
      barrier();
      SOUT("After reserving space for local hash tables, ", (init_mem_free - get_free_mem_gb()), "GB memory was used on node 0\n");
    }
    start_t = std::chrono::high_resolution_clock::now();
  }
  
  void clear() {
    for (auto it = kmers->begin(); it != kmers->end(); ) {
      it = kmers->erase(it);
    }
  }
  
  ~KmerDHT() {
    kmer_store_bloom.clear();
    kmer_store.clear();
    clear();
  }

  int64_t get_num_kmers(bool all = false) {
    if (!all) return reduce_one(kmers->size(), op_fast_add, 0).wait();
    else return reduce_all(kmers->size(), op_fast_add).wait();
  }

  float max_load_factor() {
    return reduce_one(kmers->max_load_factor(), op_fast_max, 0).wait();
  }
  
  float load_factor() {
    int64_t cardinality = initial_kmer_dht_reservation * rank_n();
    int64_t num_kmers = get_num_kmers();
    SOUT("Originally reserved ", cardinality, " and now have ", num_kmers, " elements\n");
    return reduce_one(kmers->load_factor(), op_fast_add, 0).wait() / upcxx::rank_n();
  }
  
  int64_t get_local_num_kmers(void) {
    return kmers->size();
  }
  
  size_t get_kmer_target_rank(Kmer &kmer) {
    return std::hash<Kmer>{}(kmer) % rank_n();
  }

  KmerCounts *get_local_kmer_counts(Kmer &kmer) {
    const auto it = kmers->find(kmer);
    if (it == kmers->end()) return nullptr;
    return &it->second;
  }

  int32_t get_visited(Kmer &kmer) {
    return rpc(get_kmer_target_rank(kmer),
               [](Kmer::MERARR merarr, dist_object<kmer_map_t> &kmers) -> int32_t {
                 Kmer kmer(merarr);
                 const auto it = kmers->find(kmer);
                 if (it == kmers->end()) return -2;
                 else return it->second.visited;
               }, kmer.getArray(), kmers).wait();
  }
  
  void add_kmer(Kmer kmer, char left_ext, char right_ext, uint16_t count, PASS_TYPE pass_type) {
    // get the lexicographically smallest
    Kmer kmer_twin = kmer.twin();
    if (kmer_twin < kmer) {
      kmer = kmer_twin;
      swap(left_ext, right_ext);
      left_ext = comp_nucleotide(left_ext);
      right_ext = comp_nucleotide(right_ext);
    }
    auto target_rank = get_kmer_target_rank(kmer);
    MerarrAndExt merarr_and_ext = { kmer.getArray(), left_ext, right_ext, count };
    switch (pass_type) {
      case BLOOM_SET_PASS:
        if (count != 0) kmer_store_bloom.update(target_rank, merarr_and_ext.merarr, bloom_set, bloom_filter1, bloom_filter2);
        break;
      case CTG_BLOOM_SET_PASS:
        kmer_store_bloom.update(target_rank, merarr_and_ext.merarr, ctg_bloom_set, bloom_filter2);
        break;
      case BLOOM_COUNT_PASS:
        kmer_store.update(target_rank, merarr_and_ext, bloom_count, kmers, bloom_filter2);
        break;
      case NO_BLOOM_PASS:
        kmer_store.update(target_rank, merarr_and_ext, insert_kmer, kmers);
        break;
      case CTG_KMERS_PASS:
        kmer_store.update(target_rank, merarr_and_ext, insert_ctg_kmer, kmers, min_depth_cutoff, dynamic_min_depth);
        break;
    };
  }

  void reserve_space_and_clear_bloom1() {
    Timer timer(__func__);
    // at this point we're done with generating the bloom filters, so we can drop the first bloom filter and
    // allocate the kmer hash table

    // purge the kmer store and prep the kmer + count
    kmer_store_bloom.clear();
    kmer_store.set_size(max_kmer_store_bytes);

    int64_t cardinality1 = bloom_filter1->estimate_num_items();
    int64_t cardinality2 = bloom_filter2->estimate_num_items();
    bloom1_cardinality = cardinality1;
    SOUT("Rank 0: first bloom filter size estimate is ", cardinality1, " and second size estimate is ", cardinality2,
         " ratio is ", (double)cardinality2 / cardinality1, "\n");
    bloom_filter1->clear(); // no longer need it

    double init_mem_free = get_free_mem_gb();
    barrier();
    initial_kmer_dht_reservation = (int64_t) (cardinality2 * (1+BLOOM_FP) * (1+BLOOM_FP) + 1000); // two bloom false positive rates applied
    kmers->reserve( initial_kmer_dht_reservation );
    double kmers_space_reserved = initial_kmer_dht_reservation * (sizeof(Kmer) + sizeof(KmerCounts));
    SOUT("Rank 0 is reserving ", get_size_str(kmers_space_reserved), " for kmer hash table with ", initial_kmer_dht_reservation, " entries (",
         kmers->bucket_count(), " buckets)\n");
    barrier();
  }

  void flush_updates(PASS_TYPE pass_type) {
    Timer timer(__func__);
    barrier();
    if (pass_type == BLOOM_SET_PASS) {
        kmer_store_bloom.flush_updates(bloom_set, bloom_filter1, bloom_filter2);
    } else if (pass_type == CTG_BLOOM_SET_PASS) {
        kmer_store_bloom.flush_updates(ctg_bloom_set, bloom_filter2);
    } else if (pass_type == BLOOM_COUNT_PASS) {
        kmer_store.flush_updates(bloom_count, kmers, bloom_filter2);
    } else if (pass_type == NO_BLOOM_PASS) {
        kmer_store.flush_updates(insert_kmer, kmers);
    } else if (pass_type == CTG_KMERS_PASS) {
        kmer_store.flush_updates(insert_ctg_kmer, kmers, min_depth_cutoff, dynamic_min_depth);
    } else {
      DIE("bad pass type\n");
    }
    barrier();
  }

  void purge_kmers(int threshold) {
    Timer timer(__func__);
    auto num_prior_kmers = get_num_kmers();
    int64_t num_purged = 0;
    for (auto it = kmers->begin(); it != kmers->end(); ) {
      auto kmer_counts = make_shared<KmerCounts>(it->second);
      if ((kmer_counts->count < threshold) || (kmer_counts->left_exts.is_zero() && kmer_counts->right_exts.is_zero())) {
        num_purged++;
        it = kmers->erase(it);
      } else {
        ++it;
      }
    }
    auto all_num_purged = reduce_one(num_purged, op_fast_add, 0).wait();
    SOUT("Purged ", perc_str(all_num_purged, num_prior_kmers), " kmers below frequency threshold of ", threshold, "\n");
  }

  void purge_fx_kmers() {
    Timer timer(__func__);
    auto num_prior_kmers = get_num_kmers();
    int64_t num_purged = 0;
    for (auto it = kmers->begin(); it != kmers->end(); ) {
      auto kmer_counts = make_shared<KmerCounts>(it->second);
      if (kmer_counts->left == 'X' || kmer_counts->left == 'F' || kmer_counts->right == 'X' || kmer_counts->right == 'F') {
        num_purged++;
        it = kmers->erase(it);
      } else {
        ++it;
      }
    }
    auto all_num_purged = reduce_one(num_purged, op_fast_add, 0).wait();
    SOUT("Purged ", perc_str(all_num_purged, num_prior_kmers), " kmers with F or X extensions\n");
  }

  void compute_kmer_exts() {
    Timer timer(__func__);
    for (auto &elem : *kmers) {
      auto kmer_counts = &elem.second;
      auto exts = kmer_counts->get_exts(*min_depth_cutoff, *dynamic_min_depth);
      kmer_counts->left = exts.first;
      kmer_counts->right = exts.second;
    }
 }
  
  // one line per kmer, format:
  // KMERCHARS LR N
  // where L is left extension and R is right extension, one char, either X, F or A, C, G, T
  // where N is the count of the kmer frequency
  void dump_kmers(int k) {
    Timer timer(__func__);
    string dump_fname = string("./");
    dump_fname += "ALL_INPUTS.fofn-" + to_string(k) + ".ufx.bin";
    dump_fname += ".gz";
    get_rank_path(dump_fname, rank_me());
    zstr::ofstream dump_file(dump_fname);
    ostringstream out_buf;
    ProgressBar progbar(kmers->size(), "Dumping kmers to " + dump_fname);
    int64_t i = 0;
    for (auto &elem : *kmers) {
      out_buf << elem.first << " " << elem.second.count << " " << elem.second.left << " " << elem.second.right << endl;
      i++;
      if (!(i % 1000)) {
        dump_file << out_buf.str();
        out_buf = ostringstream();
      }
      progbar.update();
    }
    if (!out_buf.str().empty()) dump_file << out_buf.str();
    dump_file.close();
    progbar.done();
    SOUT("Dumped ", this->get_num_kmers(), " kmers\n");
    string entries_fname = dump_fname + ".entries";
    std::ofstream entries_file(entries_fname);
    entries_file << kmers->size() << endl;
    entries_file.close();
  }

  // build and output histogram
  void write_histogram() {
    Timer timer(__func__);
    SOUT("Generating histogram\n");
    shared_ptr<HistogramDHT> histogram = make_shared<HistogramDHT>();
    size_t numKmers = 0;
    for (auto &elem : *kmers) {
      size_t count = elem.second.count;
      histogram->add_count(count,1);
      numKmers++;
    }
    if (bloom1_cardinality > 0) {
        int64_t num_one_counts = bloom1_cardinality - numKmers;
        if (num_one_counts > 0) {
            histogram->add_count(1, num_one_counts);
        }
    }
    histogram->flush_updates();
    string hist_fname = "per_thread/histogram_k" + std::to_string(Kmer::k) + ".txt";
    histogram->write_histogram(hist_fname);
  }

  kmer_map_t::const_iterator local_kmers_begin() {
    return kmers->begin();
  }

  kmer_map_t::const_iterator local_kmers_end() {
    return kmers->end();
  }

  int32_t get_time_offset_us() {
    std::chrono::duration<double> t_elapsed = std::chrono::high_resolution_clock::now() - start_t;
    return std::chrono::duration_cast<std::chrono::microseconds>(t_elapsed).count();
  }      
    
};

#endif