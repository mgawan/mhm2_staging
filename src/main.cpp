#include "main.hpp"

#include "contigging.hpp"
#include "scaffolding.hpp"
#include "post_assembly.hpp"

bool _verbose = false;

StageTimers stage_timers = {
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


int main(int argc, char **argv) {
  upcxx::init();
  // we wish to have all ranks start at the same time to determine actual timing
  barrier();
  auto start_t = std::chrono::high_resolution_clock::now();
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
  std::thread *init_gpu_thread = nullptr;
  double gpu_startup_duration = 0;
  auto num_gpus_per_node = (rank_me() == 0 ? adept_sw::get_num_node_gpus() : 0);
  if (num_gpus_per_node) {
    SLOG("Using ", num_gpus_per_node, " GPUs on node 0, with ", get_size_str(adept_sw::get_tot_gpu_mem()), " available memory\n");
    init_gpu_thread = adept_sw::initialize_gpu(gpu_startup_duration);
  } else {
    SWARN("Compiled for GPUs but no GPUs available...");
  }
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
    unsigned rlen_limit = 0;
    for (auto packed_reads : packed_reads_list) {
      rlen_limit = max(rlen_limit, packed_reads->get_max_read_len());
    }
    
    if (!options->ctgs_fname.empty()) ctgs.load_contigs(options->ctgs_fname);
    std::chrono::duration<double> init_t_elapsed = std::chrono::high_resolution_clock::now() - init_start_t;
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
      SLOG_VERBOSE("Waiting for GPU to be initialized (should be noop)\n");
      init_gpu_thread->join();
      delete init_gpu_thread;
      LOG("GPU took ", gpu_startup_duration, " seconds to initialize\n");
    }
#endif
    
    // scaffolding loops
    if (options->dump_gfa) {
      if (options->scaff_kmer_lens.size()) options->scaff_kmer_lens.push_back(options->scaff_kmer_lens.back());
      else options->scaff_kmer_lens.push_back(options->kmer_lens[0]);
    }
    if (options->scaff_kmer_lens.size()) {
      if (!max_kmer_len) {
        if (options->max_kmer_len) max_kmer_len = options->max_kmer_len;
        else max_kmer_len = options->scaff_kmer_lens.front();
      }
      for (unsigned i = 0; i < options->scaff_kmer_lens.size(); ++i) {
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
    FastqReaders::close_all(); // needed to cleanup any open files in this singleton
    auto fin_start_t = std::chrono::high_resolution_clock::now();
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
    std::chrono::duration<double> fin_t_elapsed = std::chrono::high_resolution_clock::now() - fin_start_t;
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
    std::chrono::duration<double> t_elapsed = std::chrono::high_resolution_clock::now() - start_t;
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


