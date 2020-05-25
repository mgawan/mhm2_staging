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

#include "upcxx_utils.hpp"
#include "options.hpp"
#include "fastq.hpp"
#include "packed_reads.hpp"
#include "kmer.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "kmer_dht.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

bool _verbose = false;

// Implementations in various .cpp files. Declarations here to prevent explosion of header files with one function in each one
void merge_reads(vector<string> reads_fname_list, int qual_offset, double &elapsed_write_io_t,
                 vector<PackedReads*> &packed_reads_list, bool checkpoint);
uint64_t estimate_num_kmers(unsigned kmer_len, vector<PackedReads*> &packed_reads_list);
template<int MAX_K>
void analyze_kmers(unsigned kmer_len, unsigned prev_kmer_len, int qual_offset, vector<PackedReads*> &packed_reads_list,
                   double dynamic_min_depth, int dmin_thres, Contigs &ctgs, dist_object<KmerDHT<MAX_K>> &kmer_dht,
                   double &num_kmers_factor);
template<int MAX_K>
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, Contigs &my_uutigs);
template<int MAX_K>
void find_alignments(unsigned kmer_len, vector<PackedReads*> &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                     Contigs &ctgs, Alns &alns, bool compute_cigar=false, int min_ctg_len=0);
void localassm(int max_kmer_len, int kmer_len, vector<PackedReads*> &packed_reads_list, int insert_avg, int insert_stddev,
               int qual_offset, Contigs &ctgs, Alns &alns);
void traverse_ctg_graph(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, int min_ctg_print_len,
                        vector<PackedReads*> &packed_reads_list, int break_scaffolds, Contigs &ctgs, Alns &alns);
pair<int, int> calculate_insert_size(Alns &alns, int ins_avg, int ins_stddev, int max_expected_ins_size,
                                     const string &dump_large_alns_fname="");

struct StageTimers {
  IntermittentTimer *merge_reads, *cache_reads, *analyze_kmers, *dbjg_traversal, *alignments, *localassm, *cgraph, *dump_ctgs,
    *compute_kmer_depths;
};

static StageTimers stage_timers = {
  .merge_reads = new IntermittentTimer(__FILENAME__ + string(":") + "Merge reads", "Merging reads"),
  .cache_reads = new IntermittentTimer(__FILENAME__ + string(":") + "Load reads into cache", "Loading reads into cache"),
  .analyze_kmers = new IntermittentTimer(__FILENAME__ + string(":") + "Analyze kmers", "Analyzing kmers"),
  .dbjg_traversal = new IntermittentTimer(__FILENAME__ + string(":") + "Traverse deBruijn graph", "Traversing deBruijn graph"),
  .alignments = new IntermittentTimer(__FILENAME__ + string(":") + "Alignments", "Aligning reads to contigs"),
  .localassm = new IntermittentTimer(__FILENAME__ + string(":") + "Local assembly", "Locally extending ends of contigs"),
  .cgraph = new IntermittentTimer(__FILENAME__ + string(":") + "Traverse contig graph", "Traversing contig graph"),
  .dump_ctgs = new IntermittentTimer(__FILENAME__ + string(":") + "Dump contigs", "Dumping contigs"),
  .compute_kmer_depths = new IntermittentTimer(__FILENAME__ + string(":") + "Compute kmer depths", "Computing kmer depths")
};

template<int MAX_K>
void contigging(int kmer_len, int prev_kmer_len, vector<PackedReads*> packed_reads_list, Contigs &ctgs, double &num_kmers_factor,
                int &max_expected_ins_size, int &ins_avg, int &ins_stddev, shared_ptr<Options> options) {
  auto loop_start_t = chrono::high_resolution_clock::now();
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Contig generation k = ", kmer_len, KNORM, "\n");
  SLOG("\n");
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
  {
    Kmer<MAX_K>::set_k(kmer_len);
    // duration of kmer_dht
    stage_timers.analyze_kmers->start();
    int64_t my_num_kmers = estimate_num_kmers(kmer_len, packed_reads_list);
    dist_object<KmerDHT<MAX_K>> kmer_dht(world(), my_num_kmers, num_kmers_factor, max_kmer_store, options->max_rpcs_in_flight, 
                                         options->force_bloom, options->use_heavy_hitters);
    barrier();
    analyze_kmers(kmer_len, prev_kmer_len, options->qual_offset, packed_reads_list, options->dynamic_min_depth, options->dmin_thres,
                  ctgs, kmer_dht, num_kmers_factor);
    stage_timers.analyze_kmers->stop();
    barrier();
    stage_timers.dbjg_traversal->start();
    traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
    stage_timers.dbjg_traversal->stop();
  }
#ifdef DEBUG
  ctgs.dump_contigs("uutigs-" + to_string(kmer_len), 0);
#endif
  if (kmer_len < options->kmer_lens.back()) {
    Alns alns;
    stage_timers.alignments->start();
    find_alignments<MAX_K>(kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight, ctgs, alns);
    stage_timers.alignments->stop();
    barrier();
    tie(ins_avg, ins_stddev) = calculate_insert_size(alns, options->insert_size[0], options->insert_size[1],
                                                     max_expected_ins_size);
    // insert size should never be larger than this; if it is that signals some error in the assembly
    max_expected_ins_size = ins_avg + 8 * ins_stddev;
    barrier();
    stage_timers.localassm->start();
    localassm(LASSM_MAX_KMER_LEN, kmer_len, packed_reads_list, ins_avg, ins_stddev, options->qual_offset, ctgs, alns);
    stage_timers.localassm->stop();
  }
  barrier();
  if (options->checkpoint) {
    stage_timers.dump_ctgs->start();
    ctgs.dump_contigs("contigs-" + to_string(kmer_len), 0);
    stage_timers.dump_ctgs->stop();
  }
  SLOG(KBLUE "_________________________", KNORM, "\n");
  ctgs.print_stats(500);
  chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
  SLOG("\n");
  SLOG(KBLUE, "Completed contig round k = ", kmer_len, " in ", setprecision(2),
       fixed, loop_t_elapsed.count(), " s at ", get_current_time(), " (",
       get_size_str(get_free_mem()), " free memory on node 0)", KNORM, "\n");
  barrier();
}

template <int MAX_K>
void scaffolding(int scaff_kmer_len, int max_kmer_len, vector<PackedReads *> packed_reads_list, Contigs &ctgs,
                 int &max_expected_ins_size, int &ins_avg, int &ins_stddev, shared_ptr<Options> options) {
  auto loop_start_t = chrono::high_resolution_clock::now();
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Scaffolding k = ", scaff_kmer_len, KNORM, "\n");
  SLOG("\n");
  Alns alns;
  stage_timers.alignments->start();
#ifdef DEBUG
  alns.dump_alns("scaff-" + to_string(scaff_kmer_len) + ".alns.gz");
#endif
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
  find_alignments<MAX_K>(scaff_kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight, ctgs, alns);
  stage_timers.alignments->stop();
  // always recalculate the insert size because we may need it for resumes of
  // Failed runs
  tie(ins_avg, ins_stddev) = calculate_insert_size(alns, options->insert_size[0], options->insert_size[1], max_expected_ins_size);
  // insert size should never be larger than this; if it is that signals some
  // error in the assembly
  max_expected_ins_size = ins_avg + 8 * ins_stddev;
  int break_scaff_Ns = (scaff_kmer_len == options->scaff_kmer_lens.back() ? options->break_scaff_Ns : 1);
  stage_timers.cgraph->start();
  traverse_ctg_graph(ins_avg, ins_stddev, max_kmer_len, scaff_kmer_len, options->min_ctg_print_len, packed_reads_list,
                     break_scaff_Ns, ctgs, alns);
  stage_timers.cgraph->stop();
  if (scaff_kmer_len != options->scaff_kmer_lens.back()) {
    if (options->checkpoint) {
      stage_timers.dump_ctgs->start();
      ctgs.dump_contigs("scaff-contigs-" + to_string(scaff_kmer_len), 0);
      stage_timers.dump_ctgs->stop();
    }
    SLOG(KBLUE "_________________________", KNORM, "\n");
    ctgs.print_stats(options->min_ctg_print_len);
  }
  chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
  SLOG("\n");
  SLOG(KBLUE, "Completed scaffolding round k = ", scaff_kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(),
       " s at ", get_current_time(), " (", get_size_str(get_free_mem()), " free memory on node 0)", KNORM, "\n");
  barrier();
}

template <int MAX_K>
void post_assembly(int max_kmer_len, Contigs &ctgs, shared_ptr<Options> options, int max_expected_ins_size) {
  auto loop_start_t = chrono::high_resolution_clock::now();
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Post processing", KNORM, "\n\n");
  vector<PackedReads*> packed_reads_list;
  for (auto const &reads_fname : options->reads_fnames) {
    packed_reads_list.push_back(new PackedReads(options->qual_offset, reads_fname, true));
  }
  stage_timers.cache_reads->start();
  double free_mem = (!rank_me() ? get_free_mem() : 0);
  upcxx::barrier();
  for (auto packed_reads : packed_reads_list) {
    packed_reads->load_reads();
  }
  stage_timers.cache_reads->stop();
  Alns alns;
  stage_timers.alignments->start();
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
  find_alignments<MAX_K>(max_kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight, ctgs, alns, true,
                         options->min_ctg_print_len);
  stage_timers.alignments->stop();
  for (auto packed_reads : packed_reads_list) {
    delete packed_reads;
  }
  packed_reads_list.clear();
  alns.dump_single_file_alns("final_assembly.sam", true);
  calculate_insert_size(alns, options->insert_size[0], options->insert_size[1], max_expected_ins_size, "large_alns_ctgs.txt");
  SLOG("\n", KBLUE, "Aligned unmerged reads to final assembly: SAM file can be found at ", options->output_dir,
       "/final_assembly.sam", KNORM, "\n");
  SLOG(KBLUE, "_________________________", KNORM, "\n");
}

int main(int argc, char **argv) {
  upcxx::init();
  auto start_t = chrono::high_resolution_clock::now();

  auto init_start_t = chrono::high_resolution_clock::now();
  auto options = make_shared<Options>();
  if (!options->load(argc, argv)) return 0;
  ProgressBar::SHOW_PROGRESS = options->show_progress;
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;

//#ifndef DEBUG
  // pin ranks to a single core in production

  if (options->pin_by == "thread") {
    if (pin_thread(getpid(), local_team().rank_me()) == -1) SWARN("Could not pin process ", getpid(), " to logical core ", rank_me());
    else SLOG_VERBOSE("Pinned processes, with process 0 (pid ", getpid(), ") pinned to a logical core ", local_team().rank_me(), "\n");
  } else if (options->pin_by == "socket") {
    if (pin_socket() < 0) SWARN("Could not pin processes by socket\n");
    else SLOG_VERBOSE("Pinned processes by socket\n");
  } else if (options->pin_by == "core") {
    if (pin_core() < 0) SWARN("Could not pin processes by physical core\n");
    else SLOG_VERBOSE("Pinned processes by physical core\n");
  } else if (options->pin_by == "clear") {
    if (pin_clear() < 0) SWARN("Could not clear pinning of proccesses\n");
    else SLOG_VERBOSE("Cleared any pinning of processes\n");
  } else {
    assert(options->pin_by == "none");
    SLOG_VERBOSE("No process pinning enabled\n");
  }
  
//#endif
  
  if (!upcxx::rank_me()) {
    // get total file size across all libraries
    double tot_file_size = 0;
    for (auto const &reads_fname : options->reads_fnames) {
      tot_file_size += get_file_size(reads_fname);
    }
    SLOG("Total size of ", options->reads_fnames.size(), " input file", (options->reads_fnames.size() > 1 ? "s" : ""),
         " is ", get_size_str(tot_file_size), "\n");
  }

  Contigs ctgs;
  int max_kmer_len = 0;
  int max_expected_ins_size = 0;
  if (!options->post_assm_only) {
    MemoryTrackerThread memory_tracker;
    memory_tracker.start();
    SLOG(KBLUE, "Starting with ", get_size_str(get_free_mem()), " free on node 0", KNORM, "\n");
    vector<PackedReads*> packed_reads_list;
    for (auto const &reads_fname : options->reads_fnames) {
      packed_reads_list.push_back(new PackedReads(options->qual_offset, get_merged_reads_fname(reads_fname)));
    }
    double elapsed_write_io_t = 0;
    if (!options->restart) {
      // merge the reads and insert into the packed reads memory cache
      stage_timers.merge_reads->start();
      merge_reads(options->reads_fnames, options->qual_offset, elapsed_write_io_t, packed_reads_list, options->checkpoint);
      stage_timers.merge_reads->stop();
    } else {
      // since this is a restart, the merged reads should be on disk already
      stage_timers.cache_reads->start();
      double free_mem = (!rank_me() ? get_free_mem() : 0);
      upcxx::barrier();
      for (auto packed_reads : packed_reads_list) {
        packed_reads->load_reads();
      }
      stage_timers.cache_reads->stop();
      SLOG_VERBOSE(KBLUE, "Cache used ", setprecision(2), fixed, get_size_str(free_mem - get_free_mem()), " memory on node 0",
                  KNORM, "\n");
    }
    if (!options->ctgs_fname.empty()) ctgs.load_contigs(options->ctgs_fname);
    chrono::duration<double> init_t_elapsed = chrono::high_resolution_clock::now() - init_start_t;
    SLOG("\n");
    SLOG(KBLUE, "Completed initialization in ", setprecision(2), fixed, init_t_elapsed.count(), " s at ",
        get_current_time(), " (", get_size_str(get_free_mem()), " free memory on node 0)", KNORM, "\n");
    int prev_kmer_len = options->prev_kmer_len;
    double num_kmers_factor = 1.0 / 3;
    int ins_avg = 0;
    int ins_stddev = 0;

    // contigging loops
    if (options->kmer_lens.size()) {
      max_kmer_len = options->kmer_lens.back();
      for (auto kmer_len : options->kmer_lens) {
        auto max_k = (kmer_len / 32 + 1) * 32;

  #define CONTIG_K(KMER_LEN) \
        case KMER_LEN: \
          contigging<KMER_LEN>(kmer_len, prev_kmer_len, packed_reads_list, ctgs, num_kmers_factor, max_expected_ins_size, \
                               ins_avg, ins_stddev, options); \
          break

        switch (max_k) {

          CONTIG_K(32);
  #if MAX_BUILD_KMER >= 64
          CONTIG_K(64);
  #endif
  #if MAX_BUILD_KMER >= 96
          CONTIG_K(96);
  #endif
  #if MAX_BUILD_KMER >= 128
          CONTIG_K(128);
  #endif
  #if MAX_BUILD_KMER >= 160
          CONTIG_K(160);
  #endif
          default: DIE("Built for max k = ", MAX_BUILD_KMER, " not k = ", max_k);
        }
  #undef CONTIG_K

        prev_kmer_len = kmer_len;
      }
    }

    // scaffolding loops
    if (options->scaff_kmer_lens.size()) {
      if (!max_kmer_len) {
        if (options->max_kmer_len) max_kmer_len = options->max_kmer_len;
        else max_kmer_len = options->scaff_kmer_lens.front();
      }
      for (auto scaff_kmer_len : options->scaff_kmer_lens) {
        auto max_k = (scaff_kmer_len / 32 + 1) * 32;

  #define SCAFFOLD_K(KMER_LEN) \
        case KMER_LEN: \
          scaffolding<KMER_LEN>(scaff_kmer_len, max_kmer_len, packed_reads_list, ctgs, max_expected_ins_size, ins_avg,\
                                ins_stddev, options); \
          break

        switch (max_k) {

          SCAFFOLD_K(32);
  #if MAX_BUILD_KMER >= 64
          SCAFFOLD_K(64);
  #endif
  #if MAX_BUILD_KMER >= 96
          SCAFFOLD_K(96);
  #endif
  #if MAX_BUILD_KMER >= 128
          SCAFFOLD_K(128);
  #endif
  #if MAX_BUILD_KMER >= 160
          SCAFFOLD_K(160);
  #endif
          default: DIE("Built for max k = ", MAX_BUILD_KMER, " not k = ", max_k);
        }
  #undef SCAFFOLD_K

      }
    }

    // cleanup
    auto fin_start_t = chrono::high_resolution_clock::now();
    for (auto packed_reads : packed_reads_list) {
      delete packed_reads;
    }
    packed_reads_list.clear();

    // output final assembly
    SLOG(KBLUE "_________________________", KNORM, "\n");
    stage_timers.dump_ctgs->start();
    ctgs.dump_contigs("final_assembly", options->min_ctg_print_len);
    stage_timers.dump_ctgs->stop();

    SLOG(KBLUE "_________________________", KNORM, "\n");
    ctgs.print_stats(options->min_ctg_print_len);
    chrono::duration<double> fin_t_elapsed = chrono::high_resolution_clock::now() - fin_start_t;
    SLOG("\n");
    SLOG(KBLUE, "Completed finalization in ", setprecision(2), fixed, fin_t_elapsed.count(), " s at ", get_current_time(),
        " (", get_size_str(get_free_mem()), " free memory on node 0)", KNORM, "\n");

    SLOG(KBLUE "_________________________", KNORM, "\n");
    SLOG("Stage timing:\n");
    if (!options->restart) SLOG("    ", stage_timers.merge_reads->get_final(), "\n");
    else SLOG("    ", stage_timers.cache_reads->get_final(), "\n");
    SLOG("    ", stage_timers.analyze_kmers->get_final(), "\n");
    SLOG("    ", stage_timers.dbjg_traversal->get_final(), "\n");
    SLOG("    ", stage_timers.alignments->get_final(), "\n");
    SLOG("    ", stage_timers.localassm->get_final(), "\n");
    SLOG("    ", stage_timers.cgraph->get_final(), "\n");
    SLOG("    FASTQ total read time: ", FastqReader::get_io_time(), "\n");
    SLOG("    merged FASTQ write time: ", elapsed_write_io_t, "\n");
    SLOG("    Contigs write time: ", stage_timers.dump_ctgs->get_elapsed(), "\n");
    SLOG(KBLUE "_________________________", KNORM, "\n");
    memory_tracker.stop();
    chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
    SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
        " for MHMXX version ", MHMXX_VERSION, "\n");
  }
  // post processing
  if (options->post_assm_aln || options->post_assm_only) {
    if (options->post_assm_only) {
      if (!options->ctgs_fname.empty()) ctgs.load_contigs(options->ctgs_fname);
      if (!max_kmer_len) max_kmer_len = options->max_kmer_len;
    }
    auto max_k = (max_kmer_len / 32 + 1) * 32;

#define POST_ASSEMBLY(KMER_LEN)                           \
    case KMER_LEN:                                            \
      post_assembly<KMER_LEN>(max_kmer_len, ctgs, options, max_expected_ins_size); \
      break

    switch (max_k) {
      POST_ASSEMBLY(32);
  #if MAX_BUILD_KMER >= 64
      POST_ASSEMBLY(64);
  #endif
  #if MAX_BUILD_KMER >= 96
      POST_ASSEMBLY(96);
  #endif
  #if MAX_BUILD_KMER >= 128
      POST_ASSEMBLY(128);
  #endif
  #if MAX_BUILD_KMER >= 160
      POST_ASSEMBLY(160);
  #endif
      default:
        DIE("Built for maximum kmer of ", MAX_BUILD_KMER, " not ", max_k);
        break;
    }
  #undef POST_ASSEMBLY
  }

  barrier();
  upcxx::finalize();
  return 0;
}


