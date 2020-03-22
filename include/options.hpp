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

public:

  vector<string> reads_fname_list;
  vector<unsigned> kmer_lens = {};
  int max_kmer_len = 0;
  int prev_kmer_len = 0;
  vector<unsigned> scaff_kmer_lens = {};
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store_mb = 200;
  int max_ctg_cache = 1000000;
  bool use_bloom = true;
  bool cache_reads = false;
  double dynamic_min_depth = 0.9;
  int dmin_thres = 2.0;
  bool checkpoint = true;
  bool show_progress = false;
  string ctgs_fname;
  int cores_per_node = upcxx::rank_n();
#ifdef USE_KMER_DEPTHS
  string kmer_depths_fname;
#endif
  int insert_avg = 0;
  int insert_stddev = 0;
  int min_ctg_print_len = 500;
  int break_scaff_Ns = 10;

  bool load(int argc, char **argv) {
    CLI::App app("MHM (" + string(MHM_VERSION) + ")");

    string reads_fnames;
    string kmer_lens_str;
    string scaff_kmer_lens_str;
    string ins_size_params;

    app.add_option("-r, --reads", reads_fnames, "Files containing merged and unmerged reads in FASTQ format (comma separated)")
      -> required();
    app.add_option("-i, --insert", ins_size_params, "Average insert length:standard deviation in insert size")
      ->required();
    app.add_option("-k, --kmer-lens", kmer_lens_str, "kmer lengths (comma separated)");
    app.add_option("--max-kmer-len", max_kmer_len, "Maximum kmer length (only need to specify if -k is not specified)");
    app.add_option("--prev-kmer-len", prev_kmer_len, "Previous kmer length (only need to specify if restarting contigging)");
    app.add_option("-s, --scaff_kmer_lens", scaff_kmer_lens_str, "kmer lengths for scaffolding (comma separated)");
    app.add_option("-Q, --quality-offset", qual_offset, "Phred encoding offset (default " + to_string(qual_offset) + ")");
    app.add_option("-c, --contigs", ctgs_fname, "File with contigs used for restart");
#ifdef USE_KMER_DEPTHS
    app.add_option("-d, --kmer-depths", kmer_depths_fname, "File with kmer depths for restart");
#endif
    app.add_option("--dynamic-min-depth", dynamic_min_depth,
                   "Dynamic min. depth for DeBruijn graph traversal - set to 1.0 for a single genome (default " +
                   to_string(dynamic_min_depth) + ")");
    app.add_option("--min-depth-thres", dmin_thres,
                   "Absolute mininimum depth threshold for DeBruijn graph traversal (default " + to_string(dmin_thres) + ")");
    app.add_option("--max-kmer-store", max_kmer_store_mb,
                   "Maximum size for kmer store in MB (default " + to_string(max_kmer_store_mb) + "MB)");
    app.add_option("--max-ctg-cache", max_ctg_cache,
                   "Maximum entries for alignment contig cache (default " + to_string(max_ctg_cache) + ")");
    app.add_option("--cores-per-node", cores_per_node,
                   "Number of cores per node - needed for accurate memory estimate (default " +
                   to_string(cores_per_node) + ")");
    app.add_option("--min-ctg-print-len", min_ctg_print_len,
                   "Minimum length required for printing a contig in the final assembly (default " +
                   to_string(min_ctg_print_len) + ")");
    app.add_option("--break-scaff-Ns", break_scaff_Ns,
                   "Number of Ns allowed before a scaffold is broken (default " + to_string(break_scaff_Ns) + ")");
    app.add_flag("--use-bloom", use_bloom, "Use bloom filter to reduce memory at the increase of runtime");
    app.add_flag("--cache-reads", cache_reads, "Cache reads in memory");
    app.add_flag("--checkpoint", checkpoint, "Checkpoint after each contig round");
    app.add_flag("--progress", show_progress, "Show progress");
    app.add_flag("-v, --verbose", verbose, "Verbose output");

    try {
      app.parse(argc, argv);
    } catch(const CLI::ParseError &e) {
      if (upcxx::rank_me() == 0) app.exit(e);
      return false;
    }

    reads_fname_list = split(reads_fnames, ',');
    for (auto kmer_len : split(kmer_lens_str, ',')) kmer_lens.push_back(std::stoi(kmer_len.c_str()));
    for (auto scaff_kmer_len : split(scaff_kmer_lens_str, ',')) scaff_kmer_lens.push_back(std::stoi(scaff_kmer_len.c_str()));
    auto ins_size_params_list = split(ins_size_params, ':');
    insert_avg = std::stoi(ins_size_params_list[0]);
    insert_stddev = std::stoi(ins_size_params_list[1]);
    if (show_progress) verbose = true;

    set_logger_verbose(verbose);

    if (upcxx::rank_me() == 0) {
      SLOG(KLBLUE, "MHM version ", MHM_VERSION, KNORM, "\n");
      // print out all compiler definitions
      SLOG_VERBOSE(KLBLUE, "_________________________", KNORM, "\n");
      SLOG_VERBOSE(KLBLUE, "Compiler definitions:", KNORM, "\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
      for (auto &def : all_defs) SLOG_VERBOSE("  ", def, "\n");
      SLOG_VERBOSE("_________________________\n");
      SLOG("Options:\n");
      SLOG("  reads files:           ");
      for (auto read_fname : reads_fname_list) SLOG(read_fname, ",");
      SLOG("\n");
      SLOG("  kmer lengths:          ");
      for (auto kmer_len : kmer_lens) SLOG(kmer_len, ",");
      SLOG("\n");
      if (max_kmer_len) SLOG("  max kmer length:       ", max_kmer_len, "\n");
      if (prev_kmer_len) SLOG("  prev kmer length:      ", prev_kmer_len, "\n");
      SLOG("  scaffold kmer lengths: ");
      for (auto scaff_kmer_len : scaff_kmer_lens) SLOG(scaff_kmer_len, ",");
      SLOG("\n");
      SLOG("  quality offset:        ", qual_offset, "\n");
      SLOG("  max kmer store:        ", max_kmer_store_mb, "MB\n");
      SLOG("  max ctg cache:         ", max_ctg_cache, "\n");
      SLOG("  dynamic min depth:     ", dynamic_min_depth, "\n");
      SLOG("  min depth threshold:   ", dmin_thres, "\n");
      if (!ctgs_fname.empty()) SLOG("  contig file name:      ", ctgs_fname, "\n");
#ifdef USE_KMER_DEPTHS
      if (!kmer_depths_fname.empty()) SLOG("  kmer depths file name: ", kmer_depths_fname, "\n");
#endif
      SLOG("  insert sizes:          ", insert_avg, ":", insert_stddev, "\n");
      SLOG("  cores per node:        ", cores_per_node, "\n");
      SLOG("  min ctg print length:  ", min_ctg_print_len, "\n");
      SLOG("  break scaff Ns:        ", break_scaff_Ns, "\n");
      SLOG("  use bloom:             ", YES_NO(use_bloom), "\n");
      SLOG("  cache reads:           ", YES_NO(cache_reads), "\n");
      SLOG("  show progress:         ", YES_NO(show_progress), "\n");
      SLOG("  verbose:               ", YES_NO(verbose), "\n");
      SLOG("_________________________", KNORM, "\n");

    }
    auto num_nodes = upcxx::rank_n() / cores_per_node;
    SLOG("Starting run with ", upcxx::rank_n(), " processes on ", num_nodes, " node", (num_nodes > 1 ? "s" : ""), " at ",
         get_current_time(), "\n");
#ifdef DEBUG
    SWARN("Running low-performance debug mode");
#endif
    if (!upcxx::rank_me()) {
      // get total file size across all libraries
      double tot_file_size = 0;
      for (auto const &reads_fname : reads_fname_list) tot_file_size += get_file_size(reads_fname);
      SLOG("Total size of ", reads_fname_list.size(), " input file", (reads_fname_list.size() > 1 ? "s" : ""),
           " is ", get_size_str(tot_file_size), "\n");
      SLOG("Starting run at ", get_current_time(), "\n");
    }
    upcxx::barrier();
    return true;
  }
};


#endif
