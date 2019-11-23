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
  vector<unsigned> scaff_kmer_lens = {};
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store = ONE_MB;
  int max_ctg_cache = 1000000;
  bool use_bloom = true;
  double dynamic_min_depth = 0.9;
  bool checkpoint = true;
  string ctgs_fname;
  int insert_avg = 0;
  int insert_stddev = 0;
  bool compress_reads = true;

  
  bool load(int argc, char **argv) {
    //SLOG(KBLUE, "MHM version ", MHM_VERSION, KNORM, "\n");
    CLI::App app("MHM (version) " + string(MHM_VERSION));

    string reads_fnames;
    string kmer_lens_str;
    string scaff_kmer_lens_str;
    string ins_size_params;

    app.add_option("-r, --reads", reads_fnames, "Files containing merged and unmerged reads in FASTQ format (comma separated)")
      -> required();
    app.add_option("-i, --insert", ins_size_params, "Average insert length:standard deviation in insert size")
      ->required();
    app.add_option("-k, --kmer-lens", kmer_lens_str, "kmer lengths (comma separated)");
    app.add_option("--max-kmer-len", max_kmer_len, "Maximum kmer length (only need specify if -k is not specified)");
    app.add_option("-s, --scaff_kmer_lens", scaff_kmer_lens_str, "kmer lengths for scaffolding (comma separated)");
    app.add_option("-Q, --quality-offset", qual_offset, "Phred encoding offset (default " + to_string(qual_offset) + ")");
    app.add_option("-c, --contigs", ctgs_fname, "File with contigs used for restart");
    app.add_option("--dynamic-min-depth", dynamic_min_depth,
                   "Dynamic min. depth for DeBruijn graph traversal (default " + to_string(dynamic_min_depth) + ")");
    app.add_option("--max-kmer-store", max_kmer_store,
                   "Maximum size for kmer store (default " + to_string(max_kmer_store) + ")");
    app.add_option("--max-ctg-cache", max_ctg_cache,
                   "Maximum entries for alignment contig cache (default " + to_string(max_ctg_cache) + ")");
    app.add_flag("--use-bloom", use_bloom, "Use bloom filter to reduce memory at the increase of runtime");
    app.add_flag("--checkpoint", checkpoint, "Checkpoint after each contig round");
    app.add_flag("--compress-reads", compress_reads, "Compress read files after merge");
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

    if (upcxx::rank_me() == 0) {
      // print out all compiler definitions
      SLOG(KBLUE "_________________________\nCompiler definitions:\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
      for (auto &def : all_defs) SLOG("  ", def, "\n");
      SLOG(KLBLUE, "_________________________\n");
      SLOG("MHM options:\n");
      SLOG("  reads files:           ");
      for (auto read_fname : reads_fname_list) SLOG(read_fname, ",");
      SLOG("\n");
      SLOG("  kmer lengths:          ");
      for (auto kmer_len : kmer_lens) SLOG(kmer_len, ",");
      SLOG("\n");
      if (max_kmer_len) SLOG("  max kmer length:       ", max_kmer_len, "\n");
      SLOG("  scaffold kmer lengths: ");
      for (auto scaff_kmer_len : scaff_kmer_lens) SLOG(scaff_kmer_len, ",");
      SLOG("\n");
      SLOG("  quality offset:        ", qual_offset, "\n");
      SLOG("  max kmer store:        ", max_kmer_store, "\n");
      SLOG("  max ctg cache:         ", max_ctg_cache, "\n");
      SLOG("  dynamic min depth:     ", dynamic_min_depth, "\n");
      if (!ctgs_fname.empty()) SLOG("  contig file name:      ", ctgs_fname, "\n");
      SLOG("  insert sizes:          ", insert_avg, ":", insert_stddev, "\n");
      SLOG("  use bloom:             ", YES_NO(use_bloom), "\n");
      SLOG("  compress reads:        ", YES_NO(compress_reads), "\n");
      SLOG("  verbose:               ", YES_NO(verbose), "\n");
      SLOG("_________________________", KNORM, "\n");
      
      double start_mem_free = get_free_mem_gb();
      SLOG("Initial free memory on node 0: ", std::setprecision(3), std::fixed, start_mem_free, " GB\n");
      SLOG("Running on ", upcxx::rank_n(), " ranks\n");
#ifdef DEBUG
      SLOG(KLRED "WARNING: Running low-performance debug mode\n", KNORM);
#endif
    }
    upcxx::barrier();
    return true;
  }
};


#endif
