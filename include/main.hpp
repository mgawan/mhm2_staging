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

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "contigs.hpp"
#include "fastq.hpp"
#include "kcount.hpp"
#include "klign.hpp"
#include "kmer.hpp"
#include "kmer_dht.hpp"
#include "options.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils.hpp"

#ifdef ENABLE_GPUS
#include <thread>

#include "adept-sw/driver.hpp"
#endif

#if defined(ENABLE_GASNET_STATS)

extern string _gasnet_stats_stage;
extern void mhm2_trace_set_mask(const char *newmask);

#define BEGIN_GASNET_STATS(stage)                        \
  if (_gasnet_stats_stage == stage) {                    \
    mhm2_trace_set_mask("PGA");                        \
    SWARN("Collecting communication stats for ", stage); \
  }
#define END_GASNET_STATS() \
  if (_gasnet_stats_stage != "") mhm2_trace_set_mask("")

#else
#define BEGIN_GASNET_STATS(stage)
#define END_GASNET_STATS()
#endif


using namespace upcxx;
using namespace upcxx_utils;

extern bool _verbose;

// Implementations in various .cpp files. Declarations here to prevent explosion of header files with one function in each one
void merge_reads(vector<string> reads_fname_list, int qual_offset, double &elapsed_write_io_t,
                 vector<PackedReads *> &packed_reads_list, bool checkpoint);
uint64_t estimate_num_kmers(unsigned kmer_len, vector<PackedReads *> &packed_reads_list);
template <int MAX_K>
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, Contigs &my_uutigs);
void localassm(int max_kmer_len, int kmer_len, vector<PackedReads *> &packed_reads_list, int insert_avg, int insert_stddev,
               int qual_offset, Contigs &ctgs, Alns &alns);
void traverse_ctg_graph(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, int min_ctg_print_len,
                        vector<PackedReads *> &packed_reads_list, int break_scaffolds, Contigs &ctgs, Alns &alns,
                        const string &graph_fname);
pair<int, int> calculate_insert_size(Alns &alns, int ins_avg, int ins_stddev, int max_expected_ins_size,
                                     const string &dump_large_alns_fname = "");
void compute_aln_depths(const string &fname, Contigs &ctgs, Alns &alns, int kmer_len, int min_ctg_len, vector<string> &filenames, bool double_count_merged_region);

struct StageTimers {
  IntermittentTimer *merge_reads, *cache_reads, *load_ctgs, *analyze_kmers, *dbjg_traversal, *alignments, *kernel_alns, *localassm,
      *cgraph, *dump_ctgs, *compute_kmer_depths;
};

extern StageTimers stage_timers;

void mhm2_trace_setmask(const char *newmask);
