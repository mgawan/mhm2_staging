#ifndef __OPTIONS_H
#define __OPTIONS_H

#include <iostream>
#include <regex>
#include <sys/stat.h>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "CLI11.hpp"

using std::cout;
using std::endl;
using std::vector;

#define YES_NO(X) ((X) ? "YES" : "NO")

class Options {
  vector<string> splitter(string in_pattern, string& content) {
    vector<string> split_content;
    std::regex pattern(in_pattern);
    copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1),
         std::sregex_token_iterator(), back_inserter(split_content));
    return split_content;
  }

  template <typename T>
  string vec_to_str(vector<T> &vec, const string &delimiter) {
    std::ostringstream oss;
    for (auto elem : vec) {
      oss << elem;
      if (elem != vec.back()) oss << delimiter;
    }
    return oss.str();
  }

public:

  vector<string> reads_fnames;
  vector<unsigned> kmer_lens = {21, 33, 55, 77, 99};
  int max_kmer_len = 0;
  int prev_kmer_len = 0;
  vector<unsigned> scaff_kmer_lens = {99, 33};
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store_mb = 200;
  int max_ctg_cache = 0;
  bool use_bloom = true;
  bool cache_reads = true;
  double dynamic_min_depth = 0.9;
  int dmin_thres = 2.0;
  bool checkpoint = true;
  bool show_progress = false;
  string ctgs_fname;
#ifdef USE_KMER_DEPTHS
  string kmer_depths_fname;
#endif
  vector<int> insert_size = {0, 0};
  int min_ctg_print_len = 500;
  int break_scaff_Ns = 10;

  bool load(int argc, char **argv) {
    CLI::App app("MHMXX (" + string(MHMXX_VERSION) + ")");

    app.add_option("-r, --reads", reads_fnames,
                   "Files containing merged and unmerged reads in FASTQ format (comma separated)")
                   ->required() ->delimiter(',');
    app.add_option("-i, --insert", insert_size,
                   "Insert size (average:stddev)")
                   ->required() ->delimiter(':') ->expected(2);
    auto *kmer_lens_opt = app.add_option("-k, --kmer-lens", kmer_lens,
                   "kmer lengths (comma separated)\n")
                   ->delimiter(',') ->capture_default_str();
    app.add_option("--max-kmer-len", max_kmer_len,
                   "Maximum kmer length (need to specify if only scaffolding)");
    app.add_option("--prev-kmer-len", prev_kmer_len,
                   "Previous kmer length (need to specify if contigging and contig file is specified)");
    auto *scaff_kmer_lens_opt = app.add_option("-s, --scaff-kmer-lens", scaff_kmer_lens,
                   "kmer lengths for scaffolding (comma separated)")
                   ->delimiter(',') ->capture_default_str();
    app.add_option("-Q, --quality-offset", qual_offset,
                   "Phred encoding offset")
                   ->capture_default_str();
    app.add_option("-c, --contigs", ctgs_fname,
                   "File with contigs used for restart");
#ifdef USE_KMER_DEPTHS
    app.add_option("-d, --kmer-depths", kmer_depths_fname, "File with kmer depths for restart");
#endif
    app.add_option("--dynamic-min-depth", dynamic_min_depth,
                   "Dynamic min. depth for DeBruijn graph traversal - set to 1.0 for a single genome")
                   ->capture_default_str();
    app.add_option("--min-depth-thres", dmin_thres,
                   "Absolute mininimum depth threshold for DeBruijn graph traversal")
                   ->capture_default_str();
    app.add_option("--max-kmer-store", max_kmer_store_mb,
                   "Maximum size for kmer store in MB")
                   ->capture_default_str();
    app.add_option("--max-ctg-cache", max_ctg_cache,
                   "Maximum entries for alignment contig cache")
                   ->capture_default_str();
    app.add_option("--min-ctg-print-len", min_ctg_print_len,
                   "Minimum length required for printing a contig in the final assembly")
                   ->capture_default_str();
    app.add_option("--break-scaff-Ns", break_scaff_Ns,
                   "Number of Ns allowed before a scaffold is broken")
                   ->capture_default_str();
    app.add_flag("--use-bloom", use_bloom,
                 "Use bloom filter to reduce memory at the increase of runtime")
                 ->capture_default_str();
    app.add_flag("--cache-reads", cache_reads,
                 "Cache reads in memory")
                 ->capture_default_str();
    app.add_flag("--checkpoint", checkpoint,
                 "Checkpoint after each contig round")
                 ->capture_default_str();
    app.add_flag("--progress", show_progress,
                 "Show progress")
                 ->capture_default_str();
    app.add_flag("-v, --verbose", verbose,
                 "Verbose output")
                 ->capture_default_str();

    app.set_config("--config");

    try {
      app.parse(argc, argv);
    } catch(const CLI::ParseError &e) {
      if (upcxx::rank_me() == 0) app.exit(e);
      return false;
    }

    // make sure we only use defaults for kmer lens if none of them were set by the user
    if (*kmer_lens_opt && !*scaff_kmer_lens_opt) scaff_kmer_lens = {};
    if (*scaff_kmer_lens_opt && !*kmer_lens_opt) kmer_lens_opt = {};

    if (show_progress) verbose = true;
    set_logger_verbose(verbose);

    if (upcxx::rank_me() == 0) {
      SLOG(KLBLUE, "MHMXX version ", MHMXX_VERSION, KNORM, "\n");
      // print out all compiler definitions
      SLOG_VERBOSE(KLBLUE, "_________________________", KNORM, "\n");
      SLOG_VERBOSE(KLBLUE, "Compiler definitions:", KNORM, "\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
      for (auto &def : all_defs) SLOG_VERBOSE("  ", def, "\n");
      SLOG_VERBOSE("_________________________\n");
      SLOG("Options:\n");
      SLOG("  reads files:           ", vec_to_str(reads_fnames, ","), "\n");
      SLOG("  kmer lengths:          ", vec_to_str(kmer_lens, ","), "\n");
      if (max_kmer_len) SLOG("  max kmer length:       ", max_kmer_len, "\n");
      if (prev_kmer_len) SLOG("  prev kmer length:      ", prev_kmer_len, "\n");
      SLOG("  scaffold kmer lengths: ", vec_to_str(scaff_kmer_lens, ","), "\n");
      SLOG("  quality offset:        ", qual_offset, "\n");
      SLOG("  max kmer store:        ", max_kmer_store_mb, "MB\n");
      SLOG("  max ctg cache:         ", max_ctg_cache, "\n");
      SLOG("  dynamic min depth:     ", dynamic_min_depth, "\n");
      SLOG("  min depth threshold:   ", dmin_thres, "\n");
      if (!ctgs_fname.empty()) SLOG("  contig file name:      ", ctgs_fname, "\n");
#ifdef USE_KMER_DEPTHS
      if (!kmer_depths_fname.empty()) SLOG("  kmer depths file name: ", kmer_depths_fname, "\n");
#endif
      SLOG("  insert sizes:          ", vec_to_str(insert_size, ":"), "\n");
      SLOG("  min ctg print length:  ", min_ctg_print_len, "\n");
      SLOG("  break scaff Ns:        ", break_scaff_Ns, "\n");
      SLOG("  use bloom:             ", YES_NO(use_bloom), "\n");
      SLOG("  cache reads:           ", YES_NO(cache_reads), "\n");
      SLOG("  show progress:         ", YES_NO(show_progress), "\n");
      SLOG("  verbose:               ", YES_NO(verbose), "\n");
      SLOG("_________________________", KNORM, "\n");

    }
    auto num_nodes = upcxx::rank_n() / upcxx::local_team().rank_n();
    SLOG("Starting run with ", upcxx::rank_n(), " processes on ", num_nodes, " node", (num_nodes > 1 ? "s" : ""), " at ",
         get_current_time(), "\n");
#ifdef DEBUG
    SWARN("Running low-performance debug mode");
#endif
    if (!upcxx::rank_me()) {
      // get total file size across all libraries
      double tot_file_size = 0;
      for (auto const &reads_fname : reads_fnames) tot_file_size += get_file_size(reads_fname);
      SLOG("Total size of ", reads_fnames.size(), " input file", (reads_fnames.size() > 1 ? "s" : ""),
           " is ", get_size_str(tot_file_size), "\n");
      SLOG("Starting run at ", get_current_time(), "\n");
      // write out configuration file for restarts
      ofstream ofs("mhmxx.config");
      ofs << app.config_to_str(false, true);
      /*
      // check to see if mhmxx.log exists. If so, and not restarting, rename it
      if (file_exists("mhmxx.log") && !restart) {
        string new_log_fname = "mhmxx-" + get_current_time(true) + ".log";
        WARN("mhmxx.log exists: renaming to ", new_log_fname, "\n");
        if (rename("mhmxx.log", new_log_fname.c_str()) == -1) DIE("Could not rename mhmxx.log: ", strerror(errno));
      } else if (!file_exists("mhmxx.log") && restart) {
        DIE("Could not restart - missing mhmxx.log in this directory");
      }
      */
    }
    upcxx::barrier();
    return true;
  }
};


#endif
