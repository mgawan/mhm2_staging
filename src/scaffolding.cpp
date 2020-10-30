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

#include "scaffolding.hpp"

#include "aln_depths.hpp"
#include "gasnet_stats.hpp"
#include "histogrammer.hpp"
#include "klign.hpp"
#include "stage_timers.hpp"

using namespace upcxx;
using namespace upcxx_utils;

using std::fixed;
using std::setprecision;
using std::shared_ptr;
using std::string;
using std::tie;
using std::vector;

void traverse_ctg_graph(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, int min_ctg_print_len,
                        vector<PackedReads *> &packed_reads_list, int break_scaffolds, Contigs &ctgs, Alns &alns,
                        const string &graph_fname);

template <int MAX_K>
void scaffolding(int scaff_i, int max_kmer_len, int rlen_limit, vector<PackedReads *> packed_reads_list, Contigs &ctgs,
                 int &max_expected_ins_size, int &ins_avg, int &ins_stddev, shared_ptr<Options> options) {
  auto loop_start_t = std::chrono::high_resolution_clock::now();
  unsigned scaff_kmer_len = options->scaff_kmer_lens[scaff_i];
  bool gfa_iter = (options->dump_gfa && scaff_i == options->scaff_kmer_lens.size() - 1) ? true : false;
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  if (gfa_iter)
    SLOG(KBLUE, "Computing contig graph for GFA output, k = ", scaff_kmer_len, KNORM, "\n\n");
  else
    SLOG(KBLUE, "Scaffolding k = ", scaff_kmer_len, KNORM, "\n\n");
  bool is_debug = false;
#ifdef DEBUG
  is_debug = true;
#endif
  string scaff_contigs_fname("scaff-contigs-" + to_string(scaff_kmer_len) + ".fasta");
  if ((options->restart || is_debug) && file_exists(scaff_contigs_fname)) {
    SLOG_VERBOSE("(Re)loading scaffold contigs ", scaff_contigs_fname, "\n");
    ctgs.load_contigs(scaff_contigs_fname);
  } else {
    Alns alns;
    stage_timers.alignments->start();
    auto max_kmer_store = options->max_kmer_store_mb * ONE_MB;
    int seed_space = KLIGN_SEED_SPACE;
    if (options->dump_gfa && scaff_i == options->scaff_kmer_lens.size() - 1) seed_space = 1;
    BEGIN_GASNET_STATS("alignment");
    double kernel_elapsed = find_alignments<MAX_K>(scaff_kmer_len, packed_reads_list, max_kmer_store, options->max_rpcs_in_flight,
                                                   ctgs, alns, seed_space, rlen_limit, false, 0, options->ranks_per_gpu);
    END_GASNET_STATS();
    stage_timers.kernel_alns->inc_elapsed(kernel_elapsed);
    stage_timers.alignments->stop();
#ifdef DEBUG
    alns.dump_rank_file("scaff-" + to_string(scaff_kmer_len) + ".alns.gz");
#endif
    BEGIN_GASNET_STATS("alignment_depths");
    compute_aln_depths("", ctgs, alns, max_kmer_len, 0, options->reads_fnames, true);
    END_GASNET_STATS();
    // always recalculate the insert size because we may need it for resumes of failed runs
    BEGIN_GASNET_STATS("insert_size");
    tie(ins_avg, ins_stddev) = calculate_insert_size(alns, options->insert_size[0], options->insert_size[1], max_expected_ins_size);
    END_GASNET_STATS();
    // insert size should never be larger than this; if it is that signals some
    // error in the assembly
    max_expected_ins_size = ins_avg + 8 * ins_stddev;
    int break_scaff_Ns = (scaff_kmer_len == options->scaff_kmer_lens.back() ? options->break_scaff_Ns : 1);
    stage_timers.cgraph->start();
    BEGIN_GASNET_STATS("traverse_ctg_graph");
    traverse_ctg_graph(ins_avg, ins_stddev, max_kmer_len, scaff_kmer_len, options->min_ctg_print_len, packed_reads_list,
                       break_scaff_Ns, ctgs, alns, (gfa_iter ? "final_assembly" : ""));
    END_GASNET_STATS();
    stage_timers.cgraph->stop();
    ctgs.print_stats(options->min_ctg_print_len);
    int max_scaff_i = (options->dump_gfa ? options->scaff_kmer_lens.size() - 2 : options->scaff_kmer_lens.size() - 1);
    if ((is_debug || options->checkpoint) && scaff_i < max_scaff_i) {
      SLOG_VERBOSE("Saving scaffold contigs ", scaff_contigs_fname, "\n");
      stage_timers.dump_ctgs->start();
      ctgs.dump_contigs(scaff_contigs_fname, 0);
      stage_timers.dump_ctgs->stop();
    }
  }

  std::chrono::duration<double> loop_t_elapsed = std::chrono::high_resolution_clock::now() - loop_start_t;
  SLOG("\n");
  SLOG(KBLUE, "Completed ", (gfa_iter ? "GFA output" : "scaffolding"), " round k = ", scaff_kmer_len, " in ", setprecision(2),
       fixed, loop_t_elapsed.count(), " s at ", get_current_time(), " (", get_size_str(get_free_mem()), " free memory on node 0)",
       KNORM, "\n");
  barrier();
}
