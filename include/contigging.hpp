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

#include "main.hpp"

using namespace upcxx;
using namespace upcxx_utils;

using std::tie;

template <int MAX_K>

void contigging(int kmer_len, int prev_kmer_len, int rlen_limit, vector<PackedReads *> packed_reads_list, Contigs &ctgs,
                double &num_kmers_factor, int &max_expected_ins_size, int &ins_avg, int &ins_stddev, shared_ptr<Options> options) {
  auto loop_start_t = std::chrono::high_resolution_clock::now();
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Contig generation k = ", kmer_len, KNORM, "\n");
  SLOG("\n");
  bool is_debug = false;
#ifdef DEBUG
  is_debug = true;
#endif

  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;

  string uutigs_fname("uutigs-" + to_string(kmer_len) + ".fasta");
  if (options->ctgs_fname != uutigs_fname) {
    Kmer<MAX_K>::set_k(kmer_len);
    // duration of kmer_dht
    stage_timers.analyze_kmers->start();
    int64_t my_num_kmers = estimate_num_kmers(kmer_len, packed_reads_list);
    // use the max among all ranks
    my_num_kmers = upcxx::reduce_all(my_num_kmers, upcxx::op_max).wait();
    dist_object<KmerDHT<MAX_K>> kmer_dht(world(), my_num_kmers, num_kmers_factor, max_kmer_store, options->max_rpcs_in_flight,
                                         options->force_bloom, options->use_heavy_hitters);
    barrier();
    BEGIN_GASNET_STATS("kmer_analysis");
    analyze_kmers(kmer_len, prev_kmer_len, options->qual_offset, packed_reads_list, options->dmin_thres, ctgs, kmer_dht,
                  num_kmers_factor);
    END_GASNET_STATS();
    stage_timers.analyze_kmers->stop();
    barrier();
    stage_timers.dbjg_traversal->start();
    BEGIN_GASNET_STATS("dbjg_traversal");
    traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
    END_GASNET_STATS();
    stage_timers.dbjg_traversal->stop();
    if (is_debug || options->checkpoint) {
      stage_timers.dump_ctgs->start();
      ctgs.dump_contigs(uutigs_fname, 0);
      stage_timers.dump_ctgs->stop();
    }
  }

  if (kmer_len < options->kmer_lens.back()) {
    Alns alns;
    stage_timers.alignments->start();
    BEGIN_GASNET_STATS("alignment");
    double kernel_elapsed = find_alignments<MAX_K>(kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight, ctgs,
                                                   alns, KLIGN_SEED_SPACE, rlen_limit, false, 0, options->ranks_per_gpu);
    END_GASNET_STATS();
    stage_timers.kernel_alns->inc_elapsed(kernel_elapsed);
    stage_timers.alignments->stop();
    barrier();
#ifdef DEBUG
    alns.dump_rank_file("ctg-" + to_string(kmer_len) + ".alns.gz");
#endif
    tie(ins_avg, ins_stddev) = calculate_insert_size(alns, options->insert_size[0], options->insert_size[1], max_expected_ins_size);
    // insert size should never be larger than this; if it is that signals some error in the assembly
    max_expected_ins_size = ins_avg + 8 * ins_stddev;
    barrier();
    stage_timers.localassm->start();
    BEGIN_GASNET_STATS("local_assembly");
    localassm(LASSM_MAX_KMER_LEN, kmer_len, packed_reads_list, ins_avg, ins_stddev, options->qual_offset, ctgs, alns);
    END_GASNET_STATS();
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
  std::chrono::duration<double> loop_t_elapsed = std::chrono::high_resolution_clock::now() - loop_start_t;
  SLOG("\n");
  SLOG(KBLUE, "Completed contig round k = ", kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(), " s at ",
       get_current_time(), " (", get_size_str(get_free_mem()), " free memory on node 0)", KNORM, "\n");
  barrier();
}

#define __MACRO_CONTIGGING__(KMER_LEN, MODIFIER)                                                                     \
  MODIFIER void contigging<KMER_LEN>(int, int, int, vector<PackedReads *>, Contigs &, double &, int &, int &, int &, \
                                     shared_ptr<Options>);

// Reduce compile time by instantiating templates of common types
// extern template declarations are in contigging.hpp
// template instantiations each happen in src/CMakeLists via contigging-extern-template.in.cpp

__MACRO_CONTIGGING__(32, extern template);

#if MAX_BUILD_KMER >= 64

__MACRO_CONTIGGING__(64, extern template);

#endif
#if MAX_BUILD_KMER >= 96

__MACRO_CONTIGGING__(96, extern template);

#endif
#if MAX_BUILD_KMER >= 128

__MACRO_CONTIGGING__(128, extern template);

#endif
#if MAX_BUILD_KMER >= 160

__MACRO_CONTIGGING__(160, extern template);

#endif
