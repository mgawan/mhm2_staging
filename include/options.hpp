#ifndef __OPTIONS_H
#define __OPTIONS_H

#include <iostream>
#include <regex>
#include <sys/stat.h>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "argh.h"

using std::cout;
using std::endl;
using std::vector;

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
  vector<unsigned> kmer_lens = {21, 33, 55, 77, 99};
  unsigned int prev_kmer_len = 0;
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store = ONE_MB;
  int max_ctg_cache = 0;
  bool use_bloom = false;
  double dynamic_min_depth = 0.9;
  int min_depth_cutoff = 2;
  
  void load(int argc, char **argv) {
    string usage = string(argv[0]) + "\n" +
      "-r    readsfile       Files containing merged and unmerged reads in FASTQ format (comma separated)\n" + 
      "-k    kmerlens        kmer lengths\n" +
      "-p    prevkmerlen     prev kmer length used for generating contigs (optional)\n" +
      "-Q    qualoffset      Phred encoding offset\n" +
      "-m    maxkmerstore    Maximum size for kmer store\n" +
      "-C    maxctgcache     Maximum number of entries for contig cache in aligner\n" + 
      "-b    usebloom        Use bloom filter to reduce memory at the increase of runtime\n" +
      "-d    mindepthcutoff  Min. allowable depth\n" +
      "-D    dynamicmindepth Dynamic min depth setting\n" +
      "-v                    Verbose mode\n" + 
      "-h                    Display help message\n";

    SOUT(KBLUE, "MHM version ", MHM_VERSION, KNORM, "\n");

    argh::parser args;
    get_options(args, usage);
    args.parse(argc, argv);

    string reads_fnames;
    if (!(args("-r") >> reads_fnames) || args["-h"]) {
      SOUT(usage);
      exit(0);
    }
    reads_fname_list = split(reads_fnames, ',');
    
    string kmer_lens_str = "";
    if (args("-k") >> kmer_lens_str) {
      auto kmer_lens_split = split(kmer_lens_str, ',');
      kmer_lens.clear();
      for (auto kmer_len : kmer_lens_split) kmer_lens.push_back(std::stoi(kmer_len.c_str()));
    } else {
      kmer_lens_str = "";
      for (auto kmer_len : kmer_lens) kmer_lens_str += to_string(kmer_len) + ",";
    }

    args("-p") >> prev_kmer_len;
    args("-Q") >> qual_offset;
    args("-m") >> max_kmer_store;
    args("-C") >> max_ctg_cache;
    args("-d") >> min_depth_cutoff;
    args("-D") >> dynamic_min_depth;
    if (args["-b"]) use_bloom = true;
    if (args["-v"]) verbose = true;
    
    if (upcxx::rank_me() == 0) {
      // print out all compiler definitions
      SOUT(KBLUE "_________________________\nCompiler definitions:\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
      for (auto &def : all_defs) SOUT("  ", def, "\n");
      cout << KLBLUE << "_________________________\n";
      cout << "MHM options:\n";
      cout << "  (-r) reads files:           " << reads_fnames << endl;
      cout << "  (-k) kmer lengths:          " << kmer_lens_str << endl;
      if (prev_kmer_len) cout << "  (-p) prev kmer length:      " << prev_kmer_len << endl;
      cout << "  (-Q) quality offset:        " << qual_offset << endl;
      cout << "  (-m) max kmer store:        " << max_kmer_store << endl;
      cout << "  (-C) max ctg cache:         " << max_ctg_cache << endl;
      cout << "  (-d) min depth cutoff:      " << min_depth_cutoff << endl;
      cout << "  (-D) dynamic min depth:     " << dynamic_min_depth << endl;
      cout << "  (-b) use bloom:             " << use_bloom << endl;
      cout << "  (-v) verbose:               " << (verbose ? "YES" : "NO") << endl;
      cout << "_________________________\n" << KNORM;
      cout << std::flush;
      
      double start_mem_free = get_free_mem_gb();
      SOUT("Initial free memory on node 0: ", std::setprecision(3), std::fixed, start_mem_free, " GB\n");
      SOUT("Running on ", upcxx::rank_n(), " ranks\n");
#ifdef DEBUG
      SOUT(KLRED "WARNING: Running low-performance debug mode\n", KNORM);
#endif
    }
    upcxx::barrier();
  }
};


#endif
