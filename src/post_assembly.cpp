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

#include "post_assembly.hpp"

#include "aln_depths.hpp"
#include "gasnet_stats.hpp"
#include "histogrammer.hpp"
#include "klign.hpp"
#include "packed_reads.hpp"
#include "stage_timers.hpp"
#include "upcxx_utils.hpp"

using namespace upcxx_utils;

using std::shared_ptr;
using std::string;
using std::vector;

template <int MAX_K>
void post_assembly(int kmer_len, Contigs &ctgs, shared_ptr<Options> options, int max_expected_ins_size) {
  auto loop_start_t = std::chrono::high_resolution_clock::now();
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Post processing", KNORM, "\n\n");
  vector<PackedReads *> packed_reads_list;
  for (auto const &reads_fname : options->reads_fnames) {
    packed_reads_list.push_back(new PackedReads(options->qual_offset, reads_fname, true));
  }
  stage_timers.cache_reads->start();
  double free_mem = (!rank_me() ? get_free_mem() : 0);
  upcxx::barrier();
  for (auto packed_reads : packed_reads_list) {
    packed_reads->load_reads();
  }
  stage_timers.cache_reads->stop();
  unsigned rlen_limit = 0;
  for (auto packed_reads : packed_reads_list) {
    rlen_limit = max(rlen_limit, packed_reads->get_max_read_len());
  }
  Alns alns;
  stage_timers.alignments->start();
  auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
  double kernel_elapsed = find_alignments<MAX_K>(kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight, ctgs,
                                                 alns, 4, rlen_limit, true, options->min_ctg_print_len, options->ranks_per_gpu);
  stage_timers.kernel_alns->inc_elapsed(kernel_elapsed);
  stage_timers.alignments->stop();
  for (auto packed_reads : packed_reads_list) {
    delete packed_reads;
  }
  packed_reads_list.clear();
  calculate_insert_size(alns, options->insert_size[0], options->insert_size[1], max_expected_ins_size);
  if (options->post_assm_aln) {
    alns.dump_sam_file("final_assembly.sam", options->reads_fnames, ctgs, options->min_ctg_print_len);
    SLOG("\n", KBLUE, "Aligned unmerged reads to final assembly: SAM file can be found at ", options->output_dir,
         "/final_assembly.sam", KNORM, "\n");
  }
  if (options->post_assm_abundances) {
    compute_aln_depths("final_assembly_depths.txt", ctgs, alns, kmer_len, options->min_ctg_print_len, options->reads_fnames, false);
    SLOG(KBLUE, "Contig depths (abundances) can be found at ", options->output_dir, "/final_assembly_depths.txt", KNORM, "\n");
  }
  SLOG(KBLUE, "_________________________", KNORM, "\n");
}
