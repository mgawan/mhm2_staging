#pragma once

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
#include <sstream>
#include <string>
#include <vector>

#include "version.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

class Options {
  vector<string> splitter(string in_pattern, string &content);

  template <typename T>
  string vec_to_str(const vector<T> &vec, const string &delimiter = ",") {
    std::ostringstream oss;
    for (auto elem : vec) {
      oss << elem;
      if (elem != vec.back()) oss << delimiter;
    }
    return oss.str();
  }

  bool extract_previous_lens(vector<unsigned> &lens, unsigned k);

  void get_restart_options();

  void setup_output_dir();

  void setup_log_file();

 public:
  vector<string> reads_fnames;
  vector<string> paired_fnames;
  vector<unsigned> kmer_lens = {};
  int max_kmer_len = 0;
  int prev_kmer_len = 0;
  vector<unsigned> scaff_kmer_lens = {};
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store_mb = 0;  // per rank - default to use 1% of node memory
  int max_rpcs_in_flight = 100;
  bool use_heavy_hitters = false;  // only enable when files are localized
  bool force_bloom = false;
  double dynamic_min_depth = 0.9;
  int dmin_thres = 2.0;
  bool checkpoint = true;
  bool use_kmer_depths = false;
  bool post_assm_aln = false;
  bool post_assm_abundances = false;
  bool post_assm_only = false;
  bool dump_gfa = false;
  bool show_progress = false;
  string pin_by = "core";
  int ranks_per_gpu = 0;  // autodetect
  string ctgs_fname;
#ifdef USE_KMER_DEPTHS
  string kmer_depths_fname;
#endif
  vector<int> insert_size = {0, 0};
  int min_ctg_print_len = 500;
  int break_scaff_Ns = 10;
  string output_dir;
  bool restart = false;

  Options();
  ~Options();

  void cleanup();

  bool load(int argc, char **argv);
};
