//mhm - UPC++ version
//Steven Hofmeyr, LBNL, June 2019

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
#include "kmer.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "kmer_dht.hpp"

using namespace std;
using namespace upcxx;

unsigned int Kmer::k = 0;

// Implementations in various .cpp files. Declarations here to prevent explosion of header files with one function in each one
uint64_t estimate_num_kmers(unsigned kmer_len, vector<string> &reads_fname_list);
void analyze_kmers(unsigned int kmer_len, int qual_offset, vector<string> &reads_fname_list, bool use_bloom, int min_depth_cutoff,
                   double dynamic_min_depth, Contigs &ctgs, dist_object<KmerDHT> &kmer_dht);
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs);
void find_alignments(unsigned kmer_len, unsigned seed_space, vector<string> &reads_fname_list, int max_store_size,
                     int max_ctg_cache, Contigs &ctgs, Alns *alns);
void run_scaffolding(int max_kmer_len, int kmer_len, vector<string> &reads_fname_list, Contigs &ctgs, Alns &alns);


int main(int argc, char **argv) {
  upcxx::init();
  auto start_t = chrono::high_resolution_clock::now();
  double start_mem_free = get_free_mem_gb();

#ifdef DEBUG
  //time_t curr_t = std::time(nullptr);
  //string dbg_fname = "debug" + to_string(curr_t) + ".log";
  string dbg_fname = "debug.log";
  get_rank_path(dbg_fname, rank_me());
  _dbgstream.open(dbg_fname);
#endif

  auto options = make_shared<Options>();
  options->load(argc, argv);

  // get total file size across all libraries
  double tot_file_size = 0;
  if (!rank_me()) {
    for (auto const &reads_fname : options->reads_fname_list) tot_file_size += get_file_size(reads_fname);
    SOUT("Total size of ", options->reads_fname_list.size(), " input file", (options->reads_fname_list.size() > 1 ? "s" : ""),
         " is ", (tot_file_size / ONE_GB), " GB\n");
  }

  // first merge reads - the results will go in the per_rank directory
  merge_reads(options->reads_fname_list, options->qual_offset);

  Contigs ctgs;
  int max_kmer_len = options->kmer_lens.back();
  for (auto kmer_len : options->kmer_lens) {
    auto loop_start_t = chrono::high_resolution_clock::now();
    SOUT(KBLUE "_________________________\nContig generation k = ", kmer_len, "\n\n", KNORM);
    Kmer::k = kmer_len;
    {
      // scope is to ensure that kmer_dht is freed by destructor
      auto my_num_kmers = estimate_num_kmers(kmer_len, options->reads_fname_list);
      dist_object<KmerDHT> kmer_dht(world(), my_num_kmers, options->max_kmer_store, options->use_bloom);
      barrier();
      analyze_kmers(kmer_len, options->qual_offset, options->reads_fname_list, options->use_bloom, options->min_depth_cutoff,
                    options->dynamic_min_depth, ctgs, kmer_dht);
      barrier();
      traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
      // FIXME: dump as single file if checkpoint option is specified
      ctgs.dump_contigs("uutigs", kmer_len);
    }
    {
      Alns alns;
      find_alignments(kmer_len, options->seed_space, options->reads_fname_list, options->max_kmer_store, options->max_ctg_cache,
                      ctgs, &alns);
      run_scaffolding(max_kmer_len, kmer_len, options->reads_fname_list, ctgs, alns);
    }
    
    chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
    SOUT("Completed contig round k = ", kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(), " s at ",
         get_current_time(), ", free memory ", get_free_mem_gb(), " GB\n");
    barrier();
  }
  SOUT(KBLUE "_________________________\nScaffolding\n\n", KNORM);
  SOUT(KBLUE "_________________________\n\n", KNORM);
  double end_mem_free = get_free_mem_gb();
  SOUT("Final free memory on node 0: ", setprecision(3), fixed, end_mem_free,
       " GB (unreclaimed ", (start_mem_free - end_mem_free), " GB)\n");
  
  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SOUT("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
       " for MHM version ", MHM_VERSION, "\n"); 
  barrier();

#ifdef DEBUG
  _dbgstream.flush();
  _dbgstream.close();
#endif
  barrier();
  upcxx::finalize();
  return 0;
}

