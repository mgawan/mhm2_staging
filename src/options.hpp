#ifndef __OPTIONS_H
#define __OPTIONS_H

#include <iostream>
#include <regex>
#include <sys/stat.h>
#include <upcxx/upcxx.hpp>

#include "common/utils.hpp"
#include "argh.h"

using std::cout;
using std::endl;
using std::vector;

class Options {
private:
  vector<string> splitter(string in_pattern, string& content)
  {
    vector<string> split_content;
    std::regex pattern(in_pattern);
    copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1),
         std::sregex_token_iterator(),back_inserter(split_content));
    return split_content;
  }

  void get_options(argh::parser &parser, string &usage)
  {
    vector<string> lines = splitter(R"(\n)", usage);
    for (string line: lines) {
      parser.add_param(line.substr(0, 2));
    }
  }
  
public:
  vector<string> reads_fname_list;
  string ctgs_fname;
  string ctg_depths_fname;
  int kmer_len = 99;
  int prev_kmer_len = 0;
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store = ONE_MB;
  bool use_bloom = false;
  double dynamic_min_depth = 0.9;
  int min_depth_cutoff = 2;
  
  void load(int argc, char **argv) {
    string usage = string(argv[0]) + "\n" +
      "-r    readsfile       Files containing merged and unmerged reads in FASTQ format (comma separated)\n" + 
      "-c    ctgsfile        Optional file containing contigs in FASTA format\n" +
      "-f    depthsfile      Optional file containing contig depths\n" +
      "-r    readsfile       Files containing merged and unmerged reads in FASTQ format (comma separated)\n" + 
      "-k    kmerlen         kmer length\n" +
      "-p    prevkmerlen     prev kmer length used for generating contigs (optional)\n" +
      "-Q    qualoffset      Phred encoding offset\n" +
      "-m    kmerstore       Maximum size for kmer store\n" +
      "-b    usebloom        Use bloom filter to reduce memory at the increase of runtime\n" +
      "-d    mindepthcutoff  Min. allowable depth\n" +
      "-D    dynamicmindepth Dynamic min depth setting\n" +
      "-v                    Verbose mode\n" + 
      "-h                    Display help message\n";
    argh::parser args;
    get_options(args, usage);
    args.parse(argc, argv);

    string reads_fnames;
    
    if (!(args("-r") >> reads_fnames) || args["-h"]) {
      SOUT(usage);
      exit(0);
    }
    reads_fname_list = find_per_thread_files(reads_fnames, "", false);
    args("-c") >> ctgs_fname;
    args("-f") >> ctg_depths_fname;
    args("-k") >> kmer_len;
    args("-p") >> prev_kmer_len;
    args("-Q") >> qual_offset;
    args("-m") >> max_kmer_store;
    args("-d") >> min_depth_cutoff;
    args("-D") >> dynamic_min_depth;
    if (args["-b"]) use_bloom = true;
    if (args["-v"]) verbose = true;
    if (upcxx::rank_me() == 0) {
      cout << KLBLUE << "----\n";
      cout << "MHM options:\n";
      cout << "  (-r) reads files:           " << reads_fnames << endl;
      if (ctgs_fname != "") cout << "  (-c) contigs file:          " << ctgs_fname << endl;
      if (ctg_depths_fname != "") cout << "  (-f) contig depths file:    " << ctg_depths_fname << endl;
      cout << "  (-k) kmer length:           " << kmer_len << endl;
      if (prev_kmer_len) cout << "  (-p) prev kmer length:      " << prev_kmer_len << endl;
      cout << "  (-Q) quality offset:        " << qual_offset << endl;
      cout << "  (-m) kmer store:            " << max_kmer_store << endl;
      cout << "  (-d) min depth cutoff:      " << min_depth_cutoff << endl;
      cout << "  (-D) dynamic min depth:     " << dynamic_min_depth << endl;
      cout << "  (-b) use bloom:             " << use_bloom << endl;
      cout << "  (-v) verbose:               " << (verbose ? "YES" : "NO") << endl;
      cout << "----\n" << KNORM;
      cout << std::flush;
    }
    upcxx::barrier();
  }
};


#endif
