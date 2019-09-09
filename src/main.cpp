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

#include "utils.hpp"
#include "options.hpp"
#include "merge_reads.hpp"
#include "kcount.hpp"
#include "dbjg_traversal.hpp"

using namespace std;
using namespace upcxx;

unsigned int Kmer::k = 0;
unsigned int Kmer::max_k = MAX_KMER_SIZE;


int main(int argc, char **argv)
{
  upcxx::init();
  auto start_t = chrono::high_resolution_clock::now();

#ifdef DEBUG
  time_t curr_t = std::time(nullptr);
  string dbg_fname = "dbg" + to_string(curr_t) + ".log"; 
  get_rank_path(dbg_fname, rank_me());
  _dbgstream.open(dbg_fname);
#endif

  SOUT("MHM version ", MHM_VERSION, "\n");
  double start_mem_free = get_free_mem_gb();
  SOUT("Initial free memory on node 0: ", start_mem_free, "GB\n");
  SOUT("Running on ", rank_n(), " ranks\n");

  auto options = make_shared<Options>();
  options->load(argc, argv);

  // first merge reads - the results will go in the per_rank directory
  merge_reads(options->reads_fname_list, options->qual_offset);

  Kmer::k = options->kmer_len;
  
  auto my_cardinality = estimate_cardinality(options);
  dist_object<KmerDHT> kmer_dht(world(), my_cardinality, options->max_kmer_store, options->min_depth_cutoff,
                                options->dynamic_min_depth, options->use_bloom);
  barrier();

  if (options->use_bloom) {
    count_kmers(options, kmer_dht, BLOOM_SET_PASS);
    /*
    if (options->ctgs_fname != "") {
      SOUT("Scanning contigs file to populate bloom2\n");
      add_ctg_kmers(options, kmer_dht, true, 1);
    }
    */
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
  /*
  if (options->ctgs_fname != "") {
    add_ctg_kmers(options, kmer_dht, options->use_bloom, options->use_bloom ? 2 : 3);
    kmer_dht->purge_kmers(1);
  }
  */
  barrier();
  kmer_dht->compute_kmer_exts();
  kmer_dht->dump_kmers(options->kmer_len);
  barrier();
  kmer_dht->purge_fx_kmers();
  traverse_debruijn_graph(options, kmer_dht);
  double end_mem_free = get_free_mem_gb();
  SOUT("Final free memory on node 0: ", end_mem_free, "GB, used ", (start_mem_free - end_mem_free), "GB\n");
  barrier();

  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SOUT("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(), "\n"); 
  barrier();

#ifdef DEBUG
  _dbgstream.flush();
  _dbgstream.close();
#endif
  barrier();
  upcxx::finalize();
  return 0;
}

