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
#include "kmer.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "kmer_dht.hpp"

using namespace std;
using namespace upcxx;

ofstream _logstream;
bool _verbose = false;
bool _show_progress = false;


unsigned int Kmer::k = 0;

// Implementations in various .cpp files. Declarations here to prevent explosion of header files with one function in each one
void merge_reads(vector<string> reads_fname_list, int qual_offset);
uint64_t estimate_num_kmers(unsigned kmer_len, vector<string> &reads_fname_list);
void analyze_kmers(unsigned kmer_len, int qual_offset, vector<string> &reads_fname_list, bool use_bloom,
                   double dynamic_min_depth, int dmin_thres, Contigs &ctgs, dist_object<KmerDHT> &kmer_dht);
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs);
void find_alignments(unsigned kmer_len, vector<string> &reads_fname_list, 
                     int max_store_size, int max_ctg_cache, Contigs &ctgs, Alns &alns);
void localassm(int max_kmer_len, int kmer_len, vector<string> &reads_fname_list, int insert_avg, int insert_stddev,
               int qual_offset, double dynamic_min_depth, Contigs &ctgs, Alns &alns);
void traverse_ctg_graph(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, vector<string> &reads_fname_list,
                        int break_scaffolds, QualityLevel quality_level, Contigs &ctgs, Alns &alns);


int main(int argc, char **argv) {
  upcxx::init();
  init_logger();
  IntermittentTimer merge_reads_dt(__FILENAME__ + string(":") + "Merge reads", "Merging reads"),
    analyze_kmers_dt(__FILENAME__ + string(":") + "Analyze kmers", "Analyzing kmers"),
    dbjg_traversal_dt(__FILENAME__ + string(":") + "Traverse deBruijn graph", "Traversing deBruijn graph"),
    alignments_dt(__FILENAME__ + string(":") + "Alignments", "Aligning reads to contigs"),
    localassm_dt(__FILENAME__ + string(":") + "Local assembly", "Locally extending ends of contigs"),
    cgraph_dt(__FILENAME__ + string(":") + "Traverse contig graph", "Traversing contig graph");
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
  if (!options->load(argc, argv)) return 0;
  set_logger_verbose(options->verbose);
  _show_progress = options->show_progress;
  // get total file size across all libraries
  double tot_file_size = 0;
  if (!rank_me()) {
    for (auto const &reads_fname : options->reads_fname_list) tot_file_size += get_file_size(reads_fname);
    SLOG("Total size of ", options->reads_fname_list.size(), " input file", (options->reads_fname_list.size() > 1 ? "s" : ""),
         " is ", (tot_file_size / ONE_GB), " GB\n");
  }
  // first merge reads - the results will go in the per_rank directory
  {
    merge_reads_dt.start();
    merge_reads(options->reads_fname_list, options->qual_offset);
    merge_reads_dt.stop();
  }
  Contigs ctgs;
  if (!options->ctgs_fname.empty()) ctgs.load_contigs(options->ctgs_fname);
  int max_kmer_len = 0;
  if (options->kmer_lens.size()) {
    max_kmer_len = options->kmer_lens.back();
    for (auto kmer_len : options->kmer_lens) {
      auto loop_start_t = chrono::high_resolution_clock::now();
      auto free_mem = get_free_mem_gb();
      SLOG(KBLUE "_________________________\nContig generation k = ", kmer_len, "\n\n", KNORM);
      Kmer::k = kmer_len;
      // duration of kmer_dht
      analyze_kmers_dt.start();
      auto my_num_kmers = estimate_num_kmers(kmer_len, options->reads_fname_list);
      dist_object<KmerDHT> kmer_dht(world(), my_num_kmers, options->max_kmer_store, options->use_bloom);
      barrier();
      analyze_kmers(kmer_len, options->qual_offset, options->reads_fname_list, options->use_bloom,
                    options->dynamic_min_depth, options->dmin_thres, ctgs, kmer_dht);
      analyze_kmers_dt.stop();
      barrier();
      dbjg_traversal_dt.start();
      traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
      dbjg_traversal_dt.stop();
      barrier();
      //if (kmer_len < max_kmer_len) compute_kmer_ctg_depths(kmer_len, kmer_dht, ctgs);
#ifdef DEBUG
      ctgs.dump_contigs("uutigs-" + to_string(kmer_len), 0);
#endif
      if (kmer_len < options->kmer_lens.back()) {
        Alns alns;
        alignments_dt.start();
        find_alignments(kmer_len, options->reads_fname_list, options->max_kmer_store, options->max_ctg_cache, ctgs, alns);
        alignments_dt.stop();
        barrier();
        localassm_dt.start();
        localassm(LASSM_MAX_KMER_LEN, kmer_len, options->reads_fname_list, options->insert_avg, options->insert_stddev,
                  options->qual_offset, options->dynamic_min_depth, ctgs, alns);
        localassm_dt.stop();
      }
      if (options->checkpoint) ctgs.dump_contigs("contigs-" + to_string(kmer_len), 0);
      SLOG(KBLUE "_________________________\n", KNORM);
      ctgs.print_stats(500);
      chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
      SLOG(KBLUE, "\nCompleted contig round k = ", kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(), " s at ",
           get_current_time(), ", used ", (free_mem - get_free_mem_gb()), " GB memory\n", KNORM);
      barrier();
    }
  }
  if (options->scaff_kmer_lens.size()) {
    if (!max_kmer_len) {
      if (options->max_kmer_len) max_kmer_len = options->max_kmer_len;
      else max_kmer_len = options->scaff_kmer_lens.front();
    }
    for (auto scaff_kmer_len : options->scaff_kmer_lens) {
      auto loop_start_t = chrono::high_resolution_clock::now();
      auto free_mem = get_free_mem_gb();
      Kmer::k = scaff_kmer_len;
      SLOG(KBLUE "_________________________\nScaffolding k = ", scaff_kmer_len, "\n\n", KNORM);
      Alns alns;
      alignments_dt.start();
      find_alignments(scaff_kmer_len, options->reads_fname_list, options->max_kmer_store, options->max_ctg_cache, ctgs, alns);
#ifdef DEBUG      
      alns.dump_alns("scaff-" + to_string(scaff_kmer_len) + ".alns.gz");
#endif
      alignments_dt.stop();
      int break_scaff_Ns = (scaff_kmer_len == options->scaff_kmer_lens.back() ? BREAK_SCAFF_NS : 1);
      cgraph_dt.start();
      traverse_ctg_graph(options->insert_avg, options->insert_stddev, max_kmer_len, scaff_kmer_len, options->reads_fname_list,
                         break_scaff_Ns, QualityLevel::ALL, ctgs, alns);
      cgraph_dt.stop();
      if (scaff_kmer_len != options->scaff_kmer_lens.back()) {
        if (options->checkpoint) ctgs.dump_contigs("scaff-contigs-" + to_string(scaff_kmer_len), 0);
        SLOG(KBLUE "_________________________\n", KNORM);
        ctgs.print_stats(ASSM_CLEN_THRES);
      }
      chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
      SLOG(KBLUE, "\nCompleted scaffolding round k = ", scaff_kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(), " s at ",
           get_current_time(), ", used ", (free_mem - get_free_mem_gb()), " GB memory\n", KNORM);
      barrier();
    }
  }
  SLOG(KBLUE "_________________________\n", KNORM);
  ctgs.dump_contigs("final_assembly", MIN_CTG_PRINT_LEN);
  SLOG(KBLUE "_________________________\n", KNORM);
  ctgs.print_stats(ASSM_CLEN_THRES);
  SLOG(KBLUE "_________________________\n", KNORM);
  SLOG("Stage timing:\n");
  SLOG("    " + merge_reads_dt.get_final() + "\n");
  SLOG("    " + analyze_kmers_dt.get_final() + "\n");
  SLOG("    " + dbjg_traversal_dt.get_final() + "\n");
  SLOG("    " + alignments_dt.get_final() + "\n");
  SLOG("    " + localassm_dt.get_final() + "\n");
  SLOG("    " + cgraph_dt.get_final() + "\n");
  SLOG(KBLUE "_________________________\n", KNORM);
  double end_mem_free = get_free_mem_gb();
  SLOG("Final free memory on node 0: ", setprecision(3), fixed, end_mem_free,
       " GB (unreclaimed ", (start_mem_free - end_mem_free), " GB)\n");
  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
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


