#ifndef __OPTIONS_H
#define __OPTIONS_H

#include <iostream>
#include <regex>
#include <sys/stat.h>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "argh.hpp"

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

  void get_options(argh::parser &parser, string &usage) {
    vector<string> lines = splitter(R"(\n)", usage);
    for (string line: lines) {
      parser.add_param(line.substr(0, 2));
    }
  }
  
public:
  
  vector<string> reads_fname_list;
  vector<unsigned> kmer_lens = {};
  vector<unsigned> scaff_kmer_lens = {};
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store = ONE_MB;
  int max_ctg_cache = 1000000;
  bool use_bloom = false;
  double dynamic_min_depth = 0.9;
  bool checkpoint = false;
  string ctgs_fname;
  int insert_avg = 320;
  int insert_stddev = 30;

  
  void load(int argc, char **argv) {
    string usage = string(argv[0]) + "\n" +
      "-r    readsfile       Files containing merged and unmerged reads in FASTQ format (comma separated)\n" + 
      "-k    kmerlens        kmer lengths\n" +
      "-k    prevkmerlen     prev kmer length used for generating contigs (optional)\n" +
      "-Q    qualoffset      Phred encoding offset\n" +
      "-m    maxkmerstore    Maximum size for kmer store\n" +
      "-C    maxctgcache     Maximum number of entries for contig cache in aligner\n" + 
      "-d    mindepthcutoff  Min. allowable depth\n" +
      "-D    dynamicmindepth Dynamic min depth setting\n" +
      "-c    string          File with contigs for restart\n" +
      "-b    bool            Use bloom filter to reduce memory at the increase of runtime\n" +
      "-x    bool            Checkpoint after each contig round\n" +
      "-s    scaffkmerlens   Kmer lengths for scaffolding rounds\n" +
      "-i    insavg:stddev   Average insert length:standard deviation in insert size\n" + 
      "-v                    Verbose mode\n" + 
      "-h                    Display help message\n";

    SLOG(KBLUE, "MHM version ", MHM_VERSION, KNORM, "\n");

    argh::parser args;
    get_options(args, usage);
    args.parse(argc, argv);

    string reads_fnames;
    if (!(args("-r") >> reads_fnames) || args["-h"]) {
      SOUT(usage);
      upcxx::barrier();
      upcxx::finalize();
      exit(0);
    }
    reads_fname_list = split(reads_fnames, ',');
    
    string kmer_lens_str = "";
    if (args("-k") >> kmer_lens_str) {
      auto kmer_lens_split = split(kmer_lens_str, ',');
      for (auto kmer_len : kmer_lens_split) kmer_lens.push_back(std::stoi(kmer_len.c_str()));
    }

    string scaff_kmer_lens_str = "";
    if (args("-s") >> scaff_kmer_lens_str) {
      auto scaff_kmer_lens_split = split(scaff_kmer_lens_str, ',');
      for (auto scaff_kmer_len : scaff_kmer_lens_split) scaff_kmer_lens.push_back(std::stoi(scaff_kmer_len.c_str()));
    }

    string ins_size_params = "";
    if (args("-i") >> ins_size_params) {
      auto ins_size_params_split = split(ins_size_params, ':');
      insert_avg = std::stoi(ins_size_params_split[0]);
      insert_stddev = std::stoi(ins_size_params_split[1]);
    }
    args("-Q") >> qual_offset;
    args("-m") >> max_kmer_store;
    args("-C") >> max_ctg_cache;
    args("-D") >> dynamic_min_depth;
    args("-c") >> ctgs_fname;
    if (args["-b"]) use_bloom = true;
    if (args["-x"]) checkpoint = true;
    if (args["-v"]) verbose = true;

    if (upcxx::rank_me() == 0) {
      // print out all compiler definitions
      SLOG(KBLUE "_________________________\nCompiler definitions:\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
      for (auto &def : all_defs) SLOG("  ", def, "\n");
      SLOG(KLBLUE, "_________________________\n");
      SLOG("MHM options:\n");
      SLOG("  (-r) reads files:           ", reads_fnames, "\n");
      SLOG("  (-k) kmer lengths:          ", kmer_lens_str, "\n");
      SLOG("  (-s) scaffold kmer lengths: ", scaff_kmer_lens_str, "\n");
      SLOG("  (-Q) quality offset:        ", qual_offset, "\n");
      SLOG("  (-m) max kmer store:        ", max_kmer_store, "\n");
      SLOG("  (-C) max ctg cache:         ", max_ctg_cache, "\n");
      SLOG("  (-D) dynamic min depth:     ", dynamic_min_depth, "\n");
      if (!ctgs_fname.empty()) SLOG("  (-c) contig file name:      ", ctgs_fname, "\n");
      SLOG("  (-i) insert sizes:          ", insert_avg, ":", insert_stddev, "\n");
      SLOG("  (-b) use bloom:             ", YES_NO(use_bloom), "\n");
      SLOG("  (-x) checkpoint:            ", YES_NO(checkpoint), "\n");
      SLOG("  (-v) verbose:               ", YES_NO(verbose), "\n");
      SLOG("_________________________", KNORM, "\n");
      
      double start_mem_free = get_free_mem_gb();
      SLOG("Initial free memory on node 0: ", std::setprecision(3), std::fixed, start_mem_free, " GB\n");
      SLOG("Running on ", upcxx::rank_n(), " ranks\n");
#ifdef DEBUG
      SLOG(KLRED "WARNING: Running low-performance debug mode\n", KNORM);
#endif
    }
    upcxx::barrier();
  }
};


#endif
