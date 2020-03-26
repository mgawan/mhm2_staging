//mhm - UPC++ version
//Steven Hofmeyr, LBNL, June 2019

#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <string>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

using namespace std;

#ifndef NDEBUG
#define DEBUG
ofstream _dbgstream;
#endif

#include "utils.hpp"
#include "options.hpp"
#include "fastq.hpp"
#include "kmer.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "kmer_dht.hpp"

using namespace std;
using namespace upcxx;

ofstream _logstream;
bool _verbose = false;
bool _show_progress = false;

// Implementations in various .cpp files. Declarations here to prevent explosion of header files with one function in each one
int merge_reads(vector<string> reads_fname_list, int qual_offset, double &elapsed_write_io_t);
uint64_t estimate_num_kmers(unsigned kmer_len, vector<FastqReader*> &fqr_list);
template<int MAX_K>
void analyze_kmers(unsigned kmer_len, unsigned prev_kmer_len, int qual_offset, vector<FastqReader*> &fqr_list, 
                          bool use_bloom, double dynamic_min_depth, int dmin_thres, Contigs &ctgs, 
                          dist_object<KmerDHT<MAX_K>> &kmer_dht);
template<int MAX_K>
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, Contigs &my_uutigs);
template<int MAX_K> 
void find_alignments(unsigned kmer_len, vector<FastqReader*> &fqr_list, int max_store_size, int max_ctg_cache, Contigs &ctgs, 
                     Alns &alns);
void localassm(int max_kmer_len, int kmer_len, vector<FastqReader*> &fqr_list, int insert_avg, int insert_stddev,
               int qual_offset, Contigs &ctgs, Alns &alns);
void traverse_ctg_graph(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, int read_len, int min_ctg_print_len,
                        vector<FastqReader*> &fqr_list, int break_scaffolds, QualityLevel quality_level, 
                        Contigs &ctgs, Alns &alns);


#ifdef USE_KMER_DEPTH
static void compute_depths_of_kmers_in_ctgs(int kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, Contigs &ctgs) {
  Timer timer(__func__);
  int64_t num_zero_count = 0;
  int64_t tot_num_kmers = 0;
  vector<Kmer<MAX_K>> kmers;
  ProgressBar progbar(ctgs.size(), "Computing depths of kmers in contigs");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    //if (ctg->seq.length() >= kmer_len + 2) {
      int num_kmers = ctg->seq.length() - kmer_len;
      ctg->kmer_depths.reserve(num_kmers);
      Kmer<MAX_K>::get_kmers(kmer_len, ctg->seq, kmers);
      tot_num_kmers += kmers.size();
      for (auto kmer : kmers) {
        auto kmer_rc = kmer.revcomp();
        if (kmer_rc < kmer) kmer = kmer_rc;
        uint16_t count = kmer_dht->get_kmer_count(kmer);
        if (count == 0) {
          num_zero_count++;
          count = ctg->depth;
        }
        ctg->kmer_depths.push_back(count);
        progress();
      }
    //}
  }
  progbar.done();
  barrier();
  auto all_num_zero_count = reduce_one(num_zero_count, op_fast_add, 0).wait();
  auto all_tot_num_kmers = reduce_one(tot_num_kmers, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", all_tot_num_kmers, " kmers with ", perc_str(all_num_zero_count, all_tot_num_kmers), 
               " of zero count\n");
}
#endif

struct StageTimers {
  IntermittentTimer merge_reads, analyze_kmers, dbjg_traversal, alignments, localassm;
};

template<int MAX_K> 
void contigging(int kmer_len, int prev_kmer_len, vector<FastqReader*> fqr_list, Contigs &ctgs, 
                IntermittentTimer &analyze_kmers_dt, IntermittentTimer &dbjg_traversal_dt, IntermittentTimer &alignments_dt,
                IntermittentTimer &localassm_dt, shared_ptr<Options> options) {
  SLOG(KBLUE "_________________________\nContig generation k = ", kmer_len, KNORM, "\n\n");
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
  Kmer<MAX_K>::set_k(kmer_len);
  // duration of kmer_dht
  analyze_kmers_dt.start();
  auto my_num_kmers = estimate_num_kmers(kmer_len, fqr_list);
  dist_object<KmerDHT<MAX_K>> kmer_dht(world(), my_num_kmers, max_kmer_store, options->use_bloom);
  barrier();
  analyze_kmers(kmer_len, prev_kmer_len, options->qual_offset, fqr_list, options->use_bloom,
                options->dynamic_min_depth, options->dmin_thres, ctgs, kmer_dht);
  analyze_kmers_dt.stop();
  barrier();
  dbjg_traversal_dt.start();
  traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
  dbjg_traversal_dt.stop();
#ifdef DEBUG
  ctgs.dump_contigs("uutigs-" + to_string(kmer_len), 0);
#endif
  if (kmer_len < options->kmer_lens.back()) {
    Alns alns;
    alignments_dt.start();
    find_alignments<MAX_K>(kmer_len, fqr_list, max_kmer_store, options->max_ctg_cache, ctgs, alns);
    alignments_dt.stop();
    barrier();
    localassm_dt.start();
    localassm(LASSM_MAX_KMER_LEN, kmer_len, fqr_list, options->insert_size[0], options->insert_size[1],
              options->qual_offset, ctgs, alns);
    localassm_dt.stop();
  }
#ifdef USE_KMER_DEPTHS
  barrier();
  compute_kmer_depths_dt.start();
  compute_depths_of_kmers_in_ctgs(kmer_len, kmer_dht, ctgs);
  compute_kmer_depths_dt.stop();
#endif
  barrier();
}


int main(int argc, char **argv) {
  upcxx::init();
  init_logger();
  IntermittentTimer merge_reads_dt(__FILENAME__ + string(":") + "Merge reads", "Merging reads"),
          load_cache_dt(__FILENAME__ + string(":") + "Load reads into cache", "Loading reads into cache"),
          analyze_kmers_dt(__FILENAME__ + string(":") + "Analyze kmers", "Analyzing kmers"),
          dbjg_traversal_dt(__FILENAME__ + string(":") + "Traverse deBruijn graph", "Traversing deBruijn graph"),
          alignments_dt(__FILENAME__ + string(":") + "Alignments", "Aligning reads to contigs"),
          localassm_dt(__FILENAME__ + string(":") + "Local assembly", "Locally extending ends of contigs"),
          cgraph_dt(__FILENAME__ + string(":") + "Traverse contig graph", "Traversing contig graph"),
          dump_ctgs_dt(__FILENAME__ + string(":") + "Dump contigs", "Dumping contigs"),
          compute_kmer_depths_dt(__FILENAME__ + string(":") + "Compute kmer depths", "Computing kmer depths");
  auto start_t = chrono::high_resolution_clock::now();

#ifdef DEBUG
  //time_t curr_t = std::time(nullptr);
  //string dbg_fname = "debug" + to_string(curr_t) + ".log";
  string dbg_fname = "debug.log";
  get_rank_path(dbg_fname, rank_me());
  _dbgstream.open(dbg_fname);
#endif

  auto options = make_shared<Options>();
  if (!options->load(argc, argv)) return 0;
  _show_progress = options->show_progress;
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;

  MemoryTrackerThread memory_tracker;
  memory_tracker.start();
  // first merge reads - the results will go in the per_rank directory
  double elapsed_write_io_t = 0;
  merge_reads_dt.start();
  int read_len = merge_reads(options->reads_fnames, options->qual_offset, elapsed_write_io_t);
  merge_reads_dt.stop();
  vector<FastqReader*> fqr_list;
  for (auto const &reads_fname : options->reads_fnames) {
    fqr_list.push_back(new FastqReader(get_merged_reads_fname(reads_fname)));
  }
  if (options->cache_reads) {
    load_cache_dt.start();
    double free_mem = (!rank_me() ? get_free_mem() : 0);
    upcxx::barrier();
    for (auto fqr : fqr_list) {
      fqr->load_cache(options->qual_offset);
    }
    load_cache_dt.stop();
    SLOG(KBLUE, "Cache used ", setprecision(2), fixed, get_size_str(free_mem - get_free_mem()), " memory on node 0", 
         KNORM, "\n");
  }
  Contigs ctgs;
  if (!options->ctgs_fname.empty()) ctgs.load_contigs(options->ctgs_fname);
#ifdef USE_KMER_DEPTHS
  if (!options->kmer_depths_fname.empty()) ctgs.load_kmer_depths(options->kmer_depths_fname);
#endif
  int max_kmer_len = 0;
  int prev_kmer_len = options->prev_kmer_len;
  if (options->kmer_lens.size()) {
    max_kmer_len = options->kmer_lens.back();
    for (auto kmer_len : options->kmer_lens) {
      auto loop_start_t = chrono::high_resolution_clock::now();
      auto max_k = (kmer_len / 32 + 1) * 32;
      switch (max_k) {
        case 32:
          contigging<32>(kmer_len, prev_kmer_len, fqr_list, ctgs, analyze_kmers_dt, dbjg_traversal_dt, alignments_dt, 
                  localassm_dt, options);
          break;
        case 64:
          contigging<64>(kmer_len, prev_kmer_len, fqr_list, ctgs, analyze_kmers_dt, dbjg_traversal_dt, alignments_dt, 
                  localassm_dt, options);
          break;
        case 96:
          contigging<96>(kmer_len, prev_kmer_len, fqr_list, ctgs, analyze_kmers_dt, dbjg_traversal_dt, alignments_dt, 
                  localassm_dt, options);
          break;
        case 128:
          contigging<128>(kmer_len, prev_kmer_len, fqr_list, ctgs, analyze_kmers_dt, dbjg_traversal_dt, alignments_dt, 
                  localassm_dt, options);
          break;
        case 160:
          contigging<160>(kmer_len, prev_kmer_len, fqr_list, ctgs, analyze_kmers_dt, dbjg_traversal_dt, alignments_dt, 
                  localassm_dt, options);
          break;
      }                  
      if (options->checkpoint) {
        dump_ctgs_dt.start();
        ctgs.dump_contigs("contigs-" + to_string(kmer_len), 0);
#ifdef USE_KMER_DEPTHS
        ctgs.dump_kmer_depths("kmer_depths-" + to_string(kmer_len));
#endif
        dump_ctgs_dt.stop();
      }        
      SLOG(KBLUE "_________________________", KNORM, "\n");
      ctgs.print_stats(500);
      chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
      SLOG(KBLUE, "\nCompleted contig round k = ", kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(), 
           " s at ", get_current_time(), KNORM, "\n");
      barrier();
      prev_kmer_len = kmer_len;
    }
  }
  if (options->scaff_kmer_lens.size()) {
    if (!max_kmer_len) {
      if (options->max_kmer_len) max_kmer_len = options->max_kmer_len;
      else max_kmer_len = options->scaff_kmer_lens.front();
    }
    for (auto scaff_kmer_len : options->scaff_kmer_lens) {
      auto loop_start_t = chrono::high_resolution_clock::now();
      SLOG(KBLUE "_________________________\nScaffolding k = ", scaff_kmer_len, KNORM, "\n\n");
      Alns alns;
      alignments_dt.start();
#ifdef DEBUG      
      alns.dump_alns("scaff-" + to_string(scaff_kmer_len) + ".alns.gz");
#endif      
      auto max_k = (scaff_kmer_len / 32 + 1) * 32;
      switch (max_k) {
        case 32:
          find_alignments<32>(scaff_kmer_len, fqr_list, max_kmer_store, options->max_ctg_cache, ctgs, alns);
          break;
        case 64:
          find_alignments<64>(scaff_kmer_len, fqr_list, max_kmer_store, options->max_ctg_cache, ctgs, alns);
          break;
        case 96:
          find_alignments<96>(scaff_kmer_len, fqr_list, max_kmer_store, options->max_ctg_cache, ctgs, alns);
          break;
        case 128:
          find_alignments<128>(scaff_kmer_len, fqr_list, max_kmer_store, options->max_ctg_cache, ctgs, alns);
          break;
        case 160:
          find_alignments<160>(scaff_kmer_len, fqr_list, max_kmer_store, options->max_ctg_cache, ctgs, alns);
          break;
      }                  
      alignments_dt.stop();
      int break_scaff_Ns = (scaff_kmer_len == options->scaff_kmer_lens.back() ? options->break_scaff_Ns : 1);
      cgraph_dt.start();
      traverse_ctg_graph(options->insert_size[0], options->insert_size[1], max_kmer_len, scaff_kmer_len, read_len,
                         options->min_ctg_print_len, fqr_list, break_scaff_Ns, QualityLevel::ALL, ctgs, alns);
      cgraph_dt.stop();
      if (scaff_kmer_len != options->scaff_kmer_lens.back()) {
        if (options->checkpoint) {
          dump_ctgs_dt.start();
          ctgs.dump_contigs("scaff-contigs-" + to_string(scaff_kmer_len), 0);
          dump_ctgs_dt.stop();
        }
        SLOG(KBLUE "_________________________", KNORM, "\n");
        ctgs.print_stats(options->min_ctg_print_len);
      }
      chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
      SLOG(KBLUE, "\nCompleted scaffolding round k = ", scaff_kmer_len, " in ", setprecision(2), fixed, 
           loop_t_elapsed.count(), " s at ", get_current_time(), KNORM, "\n");
      barrier();
    }
  }
  for (auto fqr : fqr_list) {
    delete fqr;
  }  
  
  SLOG(KBLUE "_________________________", KNORM, "\n");
  ctgs.dump_contigs("final_assembly", options->min_ctg_print_len);
  SLOG(KBLUE "_________________________", KNORM, "\n");
  ctgs.print_stats(options->min_ctg_print_len);
  SLOG(KBLUE "_________________________", KNORM, "\n");
  SLOG("Stage timing:\n");
  SLOG("    ", merge_reads_dt.get_final(), "\n");
  if (options->cache_reads) SLOG("    ", load_cache_dt.get_final(), "\n");
  SLOG("    ", analyze_kmers_dt.get_final(), "\n");
  SLOG("    ", dbjg_traversal_dt.get_final(), "\n");
  SLOG("    ", alignments_dt.get_final(), "\n");
  SLOG("    ", localassm_dt.get_final(), "\n");
  SLOG("    ", cgraph_dt.get_final(), "\n");
  SLOG("    FASTQ total read time: ", FastqReader::get_io_time(), "\n");
  SLOG("    merged FASTQ write time: ", elapsed_write_io_t, "\n");
  SLOG("    Contigs write time: ", dump_ctgs_dt.get_elapsed(), "\n");
  SLOG(KBLUE "_________________________", KNORM, "\n");
  memory_tracker.stop();
  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
       " for MHMXX version ", MHMXX_VERSION, "\n");
#ifdef DEBUG
  _dbgstream.flush();
  _dbgstream.close();
#endif
  barrier();
  upcxx::finalize();
  return 0;
}


