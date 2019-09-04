// mhm - UPC++ version
// Steven Hofmeyr, LBNL, June 2019

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <chrono>
#include <upcxx/upcxx.hpp>

using namespace std;

#ifndef NDEBUG
#define DEBUG
ofstream _dbgstream;
#endif

#include "common/utils.hpp"
#include "options.hpp"

using namespace std;
using namespace upcxx;


int main(int argc, char **argv)
{
  upcxx::init();
  auto start_t = chrono::high_resolution_clock::now();

#ifdef DEBUG
  time_t curr_t = std::time(nullptr);
  string dbg_fname = "dbg" + to_string(curr_t) + ".log"; 
  get_rank_path(dbg_fname, rank_me());
  _dbgstream.open(dbg_fname);
  SOUT(KRED, "DEBUG mode - expect low performance\n", KNORM);
#endif

  SOUT("MHM version XXX, date\n");
  double start_mem_free = get_free_mem_gb();
  SOUT("Initial free memory on node 0: ", start_mem_free, "GB\n");

  auto options = make_shared<Options>();
  options->load(argc, argv);

  // first merge reads - the results will go in the per_rank directory
  
  /*
  auto my_cardinality = estimate_cardinality(options);
  Kmer::set_k(options->kmer_len);
  dist_object<KmerDHT> kmer_dht(world(), my_cardinality, options->max_kmer_store, options->min_depth_cutoff,
                                options->dynamic_min_depth, options->use_bloom);
  barrier();
  if (options->use_bloom) {
    count_kmers(options, kmer_dht, BLOOM_SET_PASS);
    if (options->ctgs_fname != "") {
      SOUT("Scanning contigs file to populate bloom2\n");
      add_ctg_kmers(options, kmer_dht, true, 1);
    }
    kmer_dht->reserve_space_and_clear_bloom1();
    count_kmers(options, kmer_dht, BLOOM_COUNT_PASS);
  } else {
    count_kmers(options, kmer_dht, NO_BLOOM_PASS);
  }
  barrier();
  SOUT("kmer DHT load factor: ", kmer_dht->load_factor(), "\n");
  barrier();
  kmer_dht->write_histogram();
  barrier();
  kmer_dht->purge_kmers(options->min_depth_cutoff);
  int64_t newCount = kmer_dht->get_num_kmers();
  SOUT("After purge of kmers <", options->min_depth_cutoff, " there are ", newCount, " unique kmers\n");
  barrier();
  if (options->ctgs_fname != "") {
    add_ctg_kmers(options, kmer_dht, options->use_bloom, options->use_bloom ? 2 : 3);
    kmer_dht->purge_kmers(1);
  }
  barrier();
  kmer_dht->compute_kmer_exts();
  kmer_dht->dump_kmers(options->kmer_len, options->cached_io);
  barrier();
  kmer_dht->purge_fx_kmers();
  traverse_debruijn_graph(options, kmer_dht);
  double end_mem_free = get_free_mem_gb();
  SOUT("Final free memory on node 0: ", end_mem_free, "GB, used ", (start_mem_free - end_mem_free), "GB\n");
  barrier();

  Timer lastly("Reductions");
  auto tot_upc_mem_leak = reduce_one(upc_mem_alloced - upc_mem_freed, op_fast_add, 0).wait();
  auto tot_upc_mem_peak = reduce_one(upc_mem_peak, op_fast_add, 0).wait();
  SOUT("Peak UPC memory ", (double)tot_upc_mem_peak / ONE_GB / rank_n(), "GB per rank\n");
  if (tot_upc_mem_leak) SOUT("Apparent memory leak of ", tot_upc_mem_leak, " across all ranks\n");

  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SOUT("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(), "\n"); 
  barrier();
*/
  DBG("test debugging\n");
#ifdef DEBUG
  _dbgstream.flush();
  _dbgstream.close();
#endif
  barrier();
  upcxx::finalize();
  return 0;
}

