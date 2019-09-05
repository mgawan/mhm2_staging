#ifndef __HISTOGRAM_DHT
#define __HISTOGRAM_DHT

#include <limits>
#include <iostream>
#include <map>
#include <deque>
#include <fstream>
#include <stdarg.h>
#include <upcxx/upcxx.hpp>

#include "bytell_hash_map.hpp"
#include "progressbar.hpp"
#include "utils.hpp"
#include "zstr.hpp"
#include "aggr_store.hpp"

using std::vector;
using std::deque;
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
using upcxx::make_view;
using upcxx::view;

class HistogramDHT {
private:
    
  //using hist_map_t = ska::bytell_hash_map<size_t, size_t>; // distributed hash to store all counts
  using hist_map_t = std::unordered_map<size_t, size_t>; // distributed hash to store all counts
  using bucket_entry_t = pair<size_t, size_t>;
  using hist_store_t = AggrStore< bucket_entry_t >;

  static const size_t HIGH_COUNT = 512;
  using cache_vector_t = vector<size_t>; // local cache to store the low order counts < HIGH_COUT before a global reduction

  dist_object<hist_map_t> hist_map;
  cache_vector_t cache_vector;
  hist_store_t hist_store;

  intrank_t get_hist_target_rank(size_t bucket) {
    return std::hash<size_t>{}(bucket) % rank_n();
  }

  struct InsertOrAddHistEntry {
    void operator()(bucket_entry_t &elem, dist_object<hist_map_t> &hist_map) {
        if (elem.second == 0) return; // do nothing
        // find it - if it isn't found then insert it, otherwise increment the counts
        const auto it = hist_map->find(elem.first);
        if (it == hist_map->end()) {
            hist_map->insert(elem);
        } else {
            assert(it->first == elem.first);
            it->second += elem.second;
        }
    }
  };
  dist_object<InsertOrAddHistEntry> insert_or_add_hist_entry;
  
  void add_low_counts() {
    Timer timer(__func__);
    for(int i = 0 ; i < HIGH_COUNT; i++) {
      if (cache_vector[i] > 0) {
        add_count(i, cache_vector[i], true);
        cache_vector[i] = 0;
      }
    }
  }

public:

  HistogramDHT() : hist_map({}), insert_or_add_hist_entry({}), hist_store({}) {
#ifdef DBG_ON
    hist_store.set_size(0); // force no store for testing
#else
    hist_store.set_size(1024*1024);
#endif
    cache_vector.resize(HIGH_COUNT, 0);
  }

  void clear() {
    for (auto it = hist_map->begin(); it != hist_map->end(); ) {
      it = hist_map->erase(it);
    }
    hist_map->clear();
    cache_vector.clear();
    cache_vector.resize(HIGH_COUNT, 0);
  }
  
  ~HistogramDHT() {
    clear();
  }

  void add_count(size_t bucket, size_t count = 1, bool forcePush = false) {
    // get the lexicographically smallest
    if (bucket < HIGH_COUNT && !forcePush) {
        cache_vector[bucket] += count;
        return;
    }
    auto target_rank = get_hist_target_rank(bucket);
    bucket_entry_t entry = {bucket, count};
    hist_store.update(target_rank, entry, insert_or_add_hist_entry, hist_map);
  }

  void flush_updates() {
    Timer timer(__func__);
    add_low_counts();
    hist_store.flush_updates(insert_or_add_hist_entry, hist_map);
  }

  // one file output by process 0
  // one line per bucket\tcount
  void write_histogram(string hist_fname) {
    Timer timer(__func__);
    size_t myMaxBucket = 0;
    for (auto elem : *hist_map) {
        if (elem.first > myMaxBucket) myMaxBucket = elem.first;
    }
    DBG("myMaxBucket=", myMaxBucket, " ", (*hist_map)[myMaxBucket], "\n");
    size_t maxBucket = reduce_one(myMaxBucket, op_fast_max, 0).wait();
    upcxx::global_ptr<size_t> sortedHistogram;
    if (rank_me() == 0) {
        sortedHistogram = upcxx::new_array<size_t>(maxBucket+1);
        size_t *h = sortedHistogram.local();
        memset(h, 0, sizeof(size_t) * (maxBucket+1));
    }
    sortedHistogram = upcxx::broadcast(sortedHistogram, 0).wait();

    // serialize local maps into a vector of non-zero map entries
    vector< bucket_entry_t > entries;
    entries.reserve( hist_map->size() );
    for (auto elem : *hist_map) {
        if (elem.second > 0) {
            entries.push_back( elem );
        }
    }

    // set the values it the array on rank 0
    intrank_t root(0);
    auto fut = rpc(
        root,
        [](upcxx::global_ptr<size_t> sh, view< bucket_entry_t > entries) {
           if (!sh.is_local()) DIE("!\n");
           size_t *hist = sh.local();
           for( auto entry : entries ) {
              hist[ entry.first ] = entry.second;
           }
        }, sortedHistogram, make_view(entries.begin(), entries.end()));
    fut.wait();
    upcxx::barrier();

    if (rank_me() == 0) {
        ofstream hist_file(hist_fname);
        size_t *sh = sortedHistogram.local();
        for( size_t i = 1; i < maxBucket+1; i++) {
            hist_file << i << "\t" << sh[i] << "\n";
        }
        hist_file.close();
        upcxx::delete_array(sortedHistogram);
    }
    sortedHistogram = NULL;
  }
};

#endif
