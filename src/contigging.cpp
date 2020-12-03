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

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORvoid shuffle_reads(int kmer_len, vector<PackedReads *>
&packed_reads_list, Alns &alns) { rant the following license: a  non-exclusive, royalty-free perpetual license to install, use,
modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.
*/

#include "contigging.hpp"

#include "gasnet_stats.hpp"
#include "histogrammer.hpp"
#include "kcount.hpp"
#include "klign.hpp"
#include "kmer_dht.hpp"
#include "stage_timers.hpp"

using namespace upcxx;
using namespace upcxx_utils;

using std::fixed;
using std::setprecision;
using std::shared_ptr;
using std::string;
using std::tie;
using std::vector;

template <int MAX_K>
void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, Contigs &my_uutigs);
void localassm(int max_kmer_len, int kmer_len, vector<PackedReads *> &packed_reads_list, int insert_avg, int insert_stddev,
               int qual_offset, Contigs &ctgs, Alns &alns);
void shuffle_reads(int qual_offset, vector<PackedReads *> &packed_reads_list, Alns &alns);

static uint64_t estimate_num_kmers(unsigned kmer_len, vector<PackedReads *> &packed_reads_list) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  int64_t num_reads = 0;
  int64_t tot_num_reads = 0;
  for (auto packed_reads : packed_reads_list) {
    tot_num_reads += packed_reads->get_local_num_reads();
    packed_reads->reset();
    string id, seq, quals;
    ProgressBar progbar(packed_reads->get_local_num_reads(), "Scanning reads to estimate number of kmers");

    for (int i = 0; i < 100000; i++) {
      if (!packed_reads->get_next_read(id, seq, quals)) break;
      progbar.update();
      // do not read the entire data set for just an estimate
      if (seq.length() < kmer_len) continue;
      num_kmers += seq.length() - kmer_len + 1;
      num_reads++;
    }
    progbar.done();
    barrier();
  }
  DBG("This rank processed ", num_reads, " reads, and found ", num_kmers, " kmers\n");
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_tot_num_reads = reduce_one(tot_num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_all(num_kmers, op_fast_add).wait();

  SLOG_VERBOSE("Processed ", perc_str(all_num_reads, all_tot_num_reads), " reads, and estimated a maximum of ",
               all_num_kmers * (all_tot_num_reads / all_num_reads), " kmers\n");
  return num_reads > 0 ? num_kmers * tot_num_reads / num_reads : 0;
}

template <int MAX_K>
void contigging(int kmer_len, int prev_kmer_len, int rlen_limit, vector<PackedReads *> &packed_reads_list, Contigs &ctgs,
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
    if (kmer_len == options->kmer_lens.front()) {
      if (options->shuffle_reads) {
        stage_timers.shuffle_reads->start();
        shuffle_reads(options->qual_offset, packed_reads_list, alns);
        stage_timers.shuffle_reads->stop();
      }
    }
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
