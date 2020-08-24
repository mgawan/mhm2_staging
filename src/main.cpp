/*
 HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from the U.S. Dept. of Energy).  All rights reserved."

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 (1) Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 (2) Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 (3) Neither the name of the University of California, Lawrence Berkeley National
 Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 endorse or promote products derived from this software without specific prior
 written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 to the features, functionality or performance of the source code ("Enhancements") to
 anyone; however, if you choose to make your Enhancements available either publicly,
 or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 written license agreement for such Enhancements, then you hereby grant the following
 license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 prepare derivative works, incorporate into other computer software, distribute, and
 sublicense such enhancements or derivative works thereof, in binary and source code
 form.
*/


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

#ifdef ENABLE_GPUS
#include <thread>
#include "adept-sw/driver.hpp"
#endif

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

ofstream _logstream;
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
double find_alignments(unsigned kmer_len, vector<PackedReads*> &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                       Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit, bool compute_cigar=false, int min_ctg_len=0, int ranks_per_gpu=0);
void localassm(int max_kmer_len, int kmer_len, vector<PackedReads*> &packed_reads_list, int insert_avg, int insert_stddev,
               int qual_offset, Contigs &ctgs, Alns &alns);
void traverse_ctg_graph(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, int min_ctg_print_len,
                        vector<PackedReads *> &packed_reads_list, int break_scaffolds, Contigs &ctgs, Alns &alns,
                        const string &graph_fname);
pair<int, int> calculate_insert_size(Alns &alns, int ins_avg, int ins_stddev, int max_expected_ins_size,
                                     const string &dump_large_alns_fname="");
void compute_aln_depths(const string &fname, Contigs &ctgs, Alns &alns, int kmer_len, int min_ctg_len, bool use_kmer_depths);


struct StageTimers {
  IntermittentTimer *merge_reads, *cache_reads, *analyze_kmers, *dbjg_traversal, *alignments, *kernel_alns, *localassm, *cgraph,
    *dump_ctgs, *compute_kmer_depths;
};

static StageTimers stage_timers = {
  .merge_reads = new IntermittentTimer(__FILENAME__ + string(":") + "Merge reads", "Merging reads"),
  .cache_reads = new IntermittentTimer(__FILENAME__ + string(":") + "Load reads into cache", "Loading reads into cache"),
  .analyze_kmers = new IntermittentTimer(__FILENAME__ + string(":") + "Analyze kmers", "Analyzing kmers"),
  .dbjg_traversal = new IntermittentTimer(__FILENAME__ + string(":") + "Traverse deBruijn graph", "Traversing deBruijn graph"),
  .alignments = new IntermittentTimer(__FILENAME__ + string(":") + "Alignments", "Aligning reads to contigs"),
  .kernel_alns = new IntermittentTimer(__FILENAME__ + string(":") + "Kernel alignments", ""),
  .localassm = new IntermittentTimer(__FILENAME__ + string(":") + "Local assembly", "Locally extending ends of contigs"),
  .cgraph = new IntermittentTimer(__FILENAME__ + string(":") + "Traverse contig graph", "Traversing contig graph"),
  .dump_ctgs = new IntermittentTimer(__FILENAME__ + string(":") + "Dump contigs", "Dumping contigs"),
  .compute_kmer_depths = new IntermittentTimer(__FILENAME__ + string(":") + "Compute kmer depths", "Computing kmer depths")
};

template<int MAX_K>
void contigging(int kmer_len, int prev_kmer_len, int rlen_limit, vector<PackedReads*> packed_reads_list, Contigs &ctgs,
                double &num_kmers_factor, int &max_expected_ins_size, int &ins_avg, int &ins_stddev, shared_ptr<Options> options) {
  auto loop_start_t = chrono::high_resolution_clock::now();
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Contig generation k = ", kmer_len, KNORM, "\n");
  SLOG("\n");
  bool is_debug = false;
#ifdef DEBUG
  is_debug = true;
#endif

  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
  string contigs_fname("uutigs-" + to_string(kmer_len) + ".fasta");
  if (options->restart && file_exists(contigs_fname)) {
    ctgs.load_contigs(contigs_fname);
  } else {
    Kmer<MAX_K>::set_k(kmer_len);
    // duration of kmer_dht
    stage_timers.analyze_kmers->start();
    int64_t my_num_kmers = estimate_num_kmers(kmer_len, packed_reads_list);
    my_num_kmers = upcxx::reduce_all(my_num_kmers, upcxx::op_max).wait();
    dist_object<KmerDHT<MAX_K>> kmer_dht(world(), my_num_kmers, num_kmers_factor, max_kmer_store, options->max_rpcs_in_flight,
                                         options->force_bloom);
    barrier();
    analyze_kmers(kmer_len, prev_kmer_len, options->qual_offset, packed_reads_list, options->dynamic_min_depth, options->dmin_thres,
                  ctgs, kmer_dht, num_kmers_factor);
    stage_timers.analyze_kmers->stop();
    barrier();
    stage_timers.dbjg_traversal->start();
    traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
    stage_timers.dbjg_traversal->stop();
    if (is_debug || options->checkpoint) {
      stage_timers.dump_ctgs->start();
      ctgs.dump_contigs(contigs_fname, 0);
      stage_timers.dump_ctgs->stop();
    }
  }
  
  if (kmer_len < options->kmer_lens.back()) {
    Alns alns;
    stage_timers.alignments->start();
    double kernel_elapsed = find_alignments<MAX_K>(kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight, ctgs,
                                                   alns, KLIGN_SEED_SPACE, rlen_limit, false, 0, options->ranks_per_gpu);
    stage_timers.kernel_alns->inc_elapsed(kernel_elapsed);
    stage_timers.alignments->stop();
    barrier();
#ifdef DEBUG
    alns.dump_alns("ctg-" + to_string(kmer_len) + ".alns.gz");
#endif
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
  if (is_debug || options->checkpoint) {
    stage_timers.dump_ctgs->start();
    string contigs_fname("contigs-" + to_string(kmer_len) + ".fasta");
    ctgs.dump_contigs(contigs_fname, 0);
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
void scaffolding(int scaff_i, int max_kmer_len, int rlen_limit, vector<PackedReads *> packed_reads_list, Contigs &ctgs,
                 int &max_expected_ins_size, int &ins_avg, int &ins_stddev, shared_ptr<Options> options) {
  auto loop_start_t = chrono::high_resolution_clock::now();
  int scaff_kmer_len = options->scaff_kmer_lens[scaff_i];
  bool gfa_iter = (options->dump_gfa && scaff_i == options->scaff_kmer_lens.size() - 1) ? true : false;
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  if (gfa_iter) SLOG(KBLUE, "Computing contig graph for GFA output, k = ", scaff_kmer_len, KNORM, "\n\n");
  else SLOG(KBLUE, "Scaffolding k = ", scaff_kmer_len, KNORM, "\n\n");
  bool is_debug = false;
#ifdef DEBUG
  is_debug = true;
#endif
  string scaff_contigs_fname("scaff-contigs-" + to_string(scaff_kmer_len) + ".fasta");
  if ((options->restart || is_debug) && file_exists(scaff_contigs_fname)) {
    SLOG_VERBOSE("(Re)loading scaffold contigs ", scaff_contigs_fname, "\n");
    ctgs.load_contigs(scaff_contigs_fname);
  } else {
    Alns alns;
    stage_timers.alignments->start();
    auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
    int seed_space = KLIGN_SEED_SPACE;
    if (options->dump_gfa && scaff_i == options->scaff_kmer_lens.size() - 1) seed_space = 4;
    double kernel_elapsed = find_alignments<MAX_K>(scaff_kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight,
                                                 ctgs, alns, seed_space, rlen_limit, false, 0, options->ranks_per_gpu);
    stage_timers.kernel_alns->inc_elapsed(kernel_elapsed);
    stage_timers.alignments->stop();
#ifdef DEBUG
    alns.dump_alns("scaff-" + to_string(scaff_kmer_len) + ".alns.gz");
#endif
    compute_aln_depths("", ctgs, alns, max_kmer_len, 0, options->use_kmer_depths);
    // always recalculate the insert size because we may need it for resumes of failed runs
    tie(ins_avg, ins_stddev) = calculate_insert_size(alns, options->insert_size[0], options->insert_size[1], max_expected_ins_size);
    // insert size should never be larger than this; if it is that signals some
    // error in the assembly
    max_expected_ins_size = ins_avg + 8 * ins_stddev;
    int break_scaff_Ns = (scaff_kmer_len == options->scaff_kmer_lens.back() ? options->break_scaff_Ns : 1);
    stage_timers.cgraph->start();
    traverse_ctg_graph(ins_avg, ins_stddev, max_kmer_len, scaff_kmer_len, options->min_ctg_print_len, packed_reads_list,
                       break_scaff_Ns, ctgs, alns, (gfa_iter ? "final_assembly" : ""));
    stage_timers.cgraph->stop();
    ctgs.print_stats(options->min_ctg_print_len);
    int max_scaff_i = (options->dump_gfa ? options->scaff_kmer_lens.size() - 2 : options->scaff_kmer_lens.size() - 1);
    if ((is_debug || options->checkpoint) && scaff_i < max_scaff_i) {
      SLOG_VERBOSE("Saving scaffold contigs ", scaff_contigs_fname, "\n");
      stage_timers.dump_ctgs->start();
      ctgs.dump_contigs(scaff_contigs_fname, 0);
      stage_timers.dump_ctgs->stop();
    }
  }

  chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
  SLOG("\n");
  SLOG(KBLUE, "Completed ", (gfa_iter ? "GFA output" : "scaffolding"), " round k = ", scaff_kmer_len, " in ", setprecision(2),
       fixed, loop_t_elapsed.count(), " s at ", get_current_time(), " (", get_size_str(get_free_mem()), " free memory on node 0)",
       KNORM, "\n");
  barrier();
}

template <int MAX_K>
void post_assembly(int kmer_len, Contigs &ctgs, shared_ptr<Options> options, int max_expected_ins_size) {
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
  int rlen_limit = 0;
  for (auto packed_reads : packed_reads_list) {
    rlen_limit = max(rlen_limit, packed_reads->get_max_read_len());
  }
  Alns alns;
  stage_timers.alignments->start();
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
  double kernel_elapsed = find_alignments<MAX_K>(kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight, ctgs,
                                                 alns, 4, rlen_limit, true, options->min_ctg_print_len, options->ranks_per_gpu);
  stage_timers.kernel_alns->inc_elapsed(kernel_elapsed);
  stage_timers.alignments->stop();
  for (auto packed_reads : packed_reads_list) {
    delete packed_reads;
  }
  packed_reads_list.clear();
  calculate_insert_size(alns, options->insert_size[0], options->insert_size[1], max_expected_ins_size);
  if (options->post_assm_aln) {
    alns.dump_single_file_alns("final_assembly.sam", true, &ctgs);
    SLOG("\n", KBLUE, "Aligned unmerged reads to final assembly: SAM file can be found at ", options->output_dir,
         "/final_assembly.sam", KNORM, "\n");
  }
  if (options->post_assm_abundances) {
    compute_aln_depths("final_assembly_depths.txt", ctgs, alns, kmer_len, options->min_ctg_print_len, options->use_kmer_depths);
    SLOG(KBLUE, "Contig depths (abundances) can be found at ", options->output_dir, "/final_assembly_depths.txt", KNORM, "\n");
  }
  SLOG(KBLUE, "_________________________", KNORM, "\n");
}

int main(int argc, char **argv) {
  upcxx::init();
  // we wish to have all ranks start at the same time to determine actual timing
  barrier();
  auto start_t = chrono::high_resolution_clock::now();
  auto init_start_t = start_t;
  auto options = make_shared<Options>();
  if (!options->load(argc, argv)) return 0;
  ProgressBar::SHOW_PROGRESS = options->show_progress;
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;

  SLOG_VERBOSE("Process 0 on node 0 is initially pinned to ", get_proc_pin(), "\n");
  // pin ranks only in production
  if (options->pin_by == "cpu") pin_cpu();
  else if (options->pin_by == "core") pin_core();
  else if (options->pin_by == "numa") pin_numa();
  
  if (!upcxx::rank_me()) {
    // get total file size across all libraries
    double tot_file_size = 0;
    for (auto const &reads_fname : options->reads_fnames) {
      tot_file_size += get_file_size(reads_fname);
    }
    SLOG("Total size of ", options->reads_fnames.size(), " input file", (options->reads_fnames.size() > 1 ? "s" : ""),
         " is ", get_size_str(tot_file_size), "\n");
  }
#ifdef ENABLE_GPUS
  // initialize the GPU and first-touch memory and functions in a new thread as this can take many seconds to complete
  double gpu_startup_duration = 0;
  int num_gpus = -1;
  size_t gpu_mem = 0;
  bool init_gpu_thread = true;
  SLOG_VERBOSE("Detecting GPUs\n");
  auto detect_gpu_fut = execute_in_new_thread(
    [&gpu_startup_duration, &num_gpus, &gpu_mem]() {
      adept_sw::initialize_gpu(gpu_startup_duration, num_gpus, gpu_mem);
  }).then(
    [&gpu_startup_duration, &num_gpus, &gpu_mem]() {
        if (num_gpus>0) {
            SLOG_VERBOSE("Using ", num_gpus, " GPUs on node 0, with ", get_size_str(gpu_mem), " available memory. Detected in ", gpu_startup_duration, " s.\n");
        } else {
            SWARN("Compiled for GPUs but no GPUs available...");
        }
    });
#endif

  Contigs ctgs;
  int max_kmer_len = 0;
  int max_expected_ins_size = 0;
  if (!options->post_assm_only) {
    MemoryTrackerThread memory_tracker("memory_tracker.log");
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
    int rlen_limit = 0;
    for (auto packed_reads : packed_reads_list) {
      rlen_limit = max(rlen_limit, packed_reads->get_max_read_len());
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
          contigging<KMER_LEN>(kmer_len, prev_kmer_len, rlen_limit, packed_reads_list, ctgs, num_kmers_factor, \
                               max_expected_ins_size, ins_avg, ins_stddev, options); \
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
    
#ifdef ENABLE_GPUS
    if (init_gpu_thread) {
      Timer t("Waiting for GPU to be initialized (should be noop)");
      init_gpu_thread = false;
      detect_gpu_fut.wait();
    }
#endif
    
    // scaffolding loops
    if (options->scaff_kmer_lens.size()) {
      if (!max_kmer_len) {
        if (options->max_kmer_len) max_kmer_len = options->max_kmer_len;
        else max_kmer_len = options->scaff_kmer_lens.front();
      }
      if (options->dump_gfa) options->scaff_kmer_lens.push_back(options->scaff_kmer_lens.back());
      for (int i = 0; i < options->scaff_kmer_lens.size(); ++i) {
        auto scaff_kmer_len = options->scaff_kmer_lens[i];
        auto max_k = (scaff_kmer_len / 32 + 1) * 32;

  #define SCAFFOLD_K(KMER_LEN) \
        case KMER_LEN: \
          scaffolding<KMER_LEN>(i, max_kmer_len, rlen_limit, packed_reads_list, ctgs, max_expected_ins_size, ins_avg, \
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
    } else {
        SLOG_VERBOSE("Skipping scaffolding stage - no scaff_kmer_lens specified\n");
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
    ctgs.dump_contigs("final_assembly.fasta", options->min_ctg_print_len);
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
    SLOG("      -> ", stage_timers.kernel_alns->get_final(), "\n");
    SLOG("    ", stage_timers.localassm->get_final(), "\n");
    SLOG("    ", stage_timers.cgraph->get_final(), "\n");
    SLOG("    FASTQ total read time: ", FastqReader::get_io_time(), "\n");
    SLOG("    merged FASTQ write time: ", elapsed_write_io_t, "\n");
    SLOG("    Contigs write time: ", stage_timers.dump_ctgs->get_elapsed(), "\n");
    SLOG(KBLUE "_________________________", KNORM, "\n");
    memory_tracker.stop();
    chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
    SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
        " for ", MHMXX_VERSION, "\n");
  }
  // post processing
  if (options->post_assm_aln || options->post_assm_only || options->post_assm_abundances) {
    int kmer_len = 33;
    if (options->post_assm_only && !options->ctgs_fname.empty()) ctgs.load_contigs(options->ctgs_fname);
    auto max_k = (kmer_len / 32 + 1) * 32;

#define POST_ASSEMBLY(KMER_LEN)                                         \
    case KMER_LEN:                                                      \
      post_assembly<KMER_LEN>(kmer_len, ctgs, options, max_expected_ins_size); \
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

#ifdef DEBUG
  _dbgstream.flush();
  while(close_dbg());
#endif
  barrier();
  upcxx::finalize();
  return 0;
}


