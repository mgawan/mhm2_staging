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
#include <iostream>
#include <upcxx/upcxx.hpp>

#include "klign.hpp"
#include "kmer.hpp"
#include "ssw.hpp"
#include "upcxx_utils/limit_outstanding.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/three_tier_aggr_store.hpp"
#include "utils.hpp"
#include "zstr.hpp"
#ifdef ENABLE_GPUS
#include "adept-sw/driver.hpp"
#endif

#ifdef __PPC64__  // FIXME remove after solving Issues #60 #35 #49
#define NO_KLIGN_CPU_WORK_STEAL
#endif

using namespace upcxx;
using namespace upcxx_utils;

using std::fixed;
using std::max;
using std::min;
using std::pair;
using std::setprecision;

using cid_t = int64_t;

struct CtgLoc {
  cid_t cid;
  global_ptr<char> seq_gptr;
  int clen;
  int pos_in_ctg;
  bool is_rc;  // FIXME pack with bitfields this bool likely adds 8 bytes!
};

struct ReadAndCtgLoc {
  int pos_in_read;
  bool read_is_rc;  // FIXME  pack this bool better, if possible
  CtgLoc ctg_loc;
};

template <int MAX_K>
struct KmerAndCtgLoc {
  Kmer<MAX_K> kmer;
  CtgLoc ctg_loc;
  UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
};

// global variables to avoid passing dist objs to rpcs
static int64_t _num_dropped_seed_to_ctgs = 0;

template <int MAX_K>
class KmerCtgDHT {
  // maps a kmer to a contig - the first part of the pair is set to true if this is a conflict,
  // with a kmer mapping to multiple contigs
  using local_kmer_map_t = HASH_TABLE<Kmer<MAX_K>, pair<bool, CtgLoc>>;
  using kmer_map_t = dist_object<local_kmer_map_t>;
  kmer_map_t kmer_map;

#ifndef FLAT_AGGR_STORE
  ThreeTierAggrStore<KmerAndCtgLoc<MAX_K>, kmer_map_t &> kmer_store;
#else
  FlatAggrStore<KmerAndCtgLoc<MAX_K>, kmer_map_t &> kmer_store;
#endif

  int64_t num_alns;
  int64_t num_perfect_alns;
  int64_t num_overlaps;

  // default aligner and filter
  StripedSmithWaterman::Aligner ssw_aligner;
  StripedSmithWaterman::Filter ssw_filter;
  AlnScoring aln_scoring;

#ifdef ENABLE_GPUS
  adept_sw::GPUDriver gpu_driver;
#endif

  vector<Aln> kernel_alns;
  vector<string> ctg_seqs;
  vector<string> read_seqs;

  future<> active_kernel_fut;
  IntermittentTimer aln_cpu_bypass_timer;

  int64_t max_clen = 0;
  int64_t max_rlen = 0;
  size_t gpu_mem_avail = 0;
  int gpu_devices = 0;

  Alns *alns;

  int64_t ctg_bytes_fetched = 0;
  HASH_TABLE<cid_t, string> ctg_cache;

  bool use_minimizers;

  int get_cigar_length(const string &cigar) {
    // check that cigar string length is the same as the sequence, but only if the sequence is included
    int base_count = 0;
    string num = "";
    for (char c : cigar) {
      switch (c) {
        case 'M':
        case 'S':
        case '=':
        case 'X':
        case 'I':
          base_count += stoi(num);
          num = "";
          break;
        case 'D':
          // base_count -= stoi(num);
          num = "";
          break;
        default:
          if (!isdigit(c)) DIE("Invalid char detected in cigar: '", c, "'");
          num += c;
          break;
      }
    }
    return base_count;
  }

  static void set_sam_string(Aln &aln, string read_seq, string cigar) {
    assert(aln.is_valid());
    aln.sam_string = aln.read_id + "\t";
    if (aln.orient == '-') {
      aln.sam_string += "16\t";
      if (read_seq != "*") read_seq = revcomp(read_seq);
      // reverse(read_quals.begin(), read_quals.end());
    } else {
      aln.sam_string += "0\t";
    }
    aln.sam_string += "Contig" + to_string(aln.cid) + "\t" + to_string(aln.cstart + 1) + "\t";
    uint32_t mapq;
    // for perfect match, set to same maximum as used by minimap or bwa
    if (aln.score2 == 0) {
      mapq = 60;
    } else {
      mapq = -4.343 * log(1 - (double)abs(aln.score1 - aln.score2) / (double)aln.score1);
      mapq = (uint32_t)(mapq + 4.99);
      mapq = mapq < 254 ? mapq : 254;
    }
    aln.sam_string += to_string(mapq) + "\t";
    // aln.sam_string += cigar + "\t*\t0\t0\t" + read_subseq + "\t*\t";
    // Don't output either the read sequence or quals - that causes the SAM file to bloat up hugely, and that info is already
    // available in the read files
    aln.sam_string += cigar + "\t*\t0\t0\t*\t*\t";
    aln.sam_string +=
        "AS:i:" + to_string(aln.score1) + "\tNM:i:" + to_string(aln.mismatches) + "\tRG:Z:" + to_string(aln.read_group_id);
    // for debugging
    // aln.sam_string += " rstart " + to_string(aln.rstart) + " rstop " + to_string(aln.rstop) + " cstop " + to_string(aln.cstop) +
    //                  " clen " + to_string(aln.clen) + " alnlen " + to_string(aln.rstop - aln.rstart);
    /*
#ifdef DEBUG
    // only used if we actually include the read seq and quals in the SAM, which we don't
    int base_count = get_cigar_length(cigar);
    if (base_count != read_seq.length())
      DIE("number of bases in cigar != aln rlen, ", base_count, " != ", read_subseq.length(), "\nsam string ", aln.sam_string);
#endif
    */
  }

  void align_read(const string &rname, int64_t cid, const string &rseq, const string &cseq, int rstart, int rlen, int cstart,
                  int clen, char orient, int overlap_len, int read_group_id, IntermittentTimer &aln_kernel_timer) {
    if (cseq.compare(0, overlap_len, rseq, rstart, overlap_len) == 0) {
      num_perfect_alns++;
      int rstop = rstart + overlap_len;
      int cstop = cstart + overlap_len;
      if (orient == '-') switch_orient(rstart, rstop, rlen);
      int score1 = overlap_len * aln_scoring.match;
      int identity = 100 * score1 / aln_scoring.match / rlen;
      Aln aln(rname, cid, rstart, rstop, rlen, cstart, cstop, clen, orient, score1, 0, identity, 0, read_group_id);
      assert(aln.is_valid());
      if (ssw_filter.report_cigar) set_sam_string(aln, rseq, to_string(overlap_len) + "=");  // exact match '=' not 'M'
      alns->add_aln(aln);
    } else {
      max_clen = max((int64_t)cseq.size(), max_clen);
      max_rlen = max((int64_t)rseq.size(), max_rlen);
      int64_t num_alns = kernel_alns.size() + 1;
      unsigned max_matrix_size = (max_clen + 1) * (max_rlen + 1);
      int64_t tot_mem_est = num_alns * (max_clen + max_rlen + 2 * sizeof(int) + 5 * sizeof(short));
      // contig is the ref, read is the query - done this way so that we can potentially do multiple alns to each read
      // this is also the way it's done in meraligner
      kernel_alns.emplace_back(rname, cid, 0, 0, rlen, cstart, 0, clen, orient);
      ctg_seqs.emplace_back(cseq);
      read_seqs.emplace_back(rseq);
      bool will_run_kernel = (tot_mem_est >= gpu_mem_avail) | (num_alns >= KLIGN_GPU_BLOCK_SIZE);
      if (will_run_kernel) {
        DBG("tot_mem_est (", tot_mem_est, ") >= gpu_mem_avail (", gpu_mem_avail, " - dispatching ", kernel_alns.size(),
            " alignments\n");
        kernel_align_block(aln_kernel_timer, read_group_id);
      }
    }
  }

  static void ssw_align_read(StripedSmithWaterman::Aligner &ssw_aligner, StripedSmithWaterman::Filter &ssw_filter, Alns *alns,
                             AlnScoring &aln_scoring, IntermittentTimer &aln_kernel_timer, Aln &aln, const string &cseq,
                             const string &rseq, int read_group_id) {
    StripedSmithWaterman::Alignment ssw_aln;
    aln_kernel_timer.start();
    // align query, ref, reflen
    ssw_aligner.Align(cseq.c_str(), rseq.c_str(), rseq.length(), ssw_filter, &ssw_aln, max((int)(rseq.length() / 2), 15));
    aln_kernel_timer.stop();

    aln.rstop = aln.rstart + ssw_aln.ref_end + 1;
    aln.rstart += ssw_aln.ref_begin;
    aln.cstop = aln.cstart + ssw_aln.query_end + 1;
    aln.cstart += ssw_aln.query_begin;
    if (aln.orient == '-') switch_orient(aln.rstart, aln.rstop, aln.rlen);

    aln.score1 = ssw_aln.sw_score;
    aln.score2 = ssw_aln.sw_score_next_best;
    aln.mismatches = ssw_aln.mismatches;
    aln.identity = 100 * ssw_aln.sw_score / aln_scoring.match / aln.rlen;
    aln.read_group_id = read_group_id;
    if (ssw_filter.report_cigar) set_sam_string(aln, rseq, ssw_aln.cigar_string);
    alns->add_aln(aln);
  }

  void ssw_align_read(IntermittentTimer &aln_kernel_timer, Aln &aln, const string &cseq, const string &rseq, int read_group_id) {
    ssw_align_read(ssw_aligner, ssw_filter, alns, aln_scoring, aln_kernel_timer, aln, cseq, rseq, read_group_id);
  }

  // encapsulate the data for a kernel to run a block independently
  struct AlignBlockData {
#ifdef ENABLE_GPUS
    adept_sw::GPUDriver &gpu_driver;
#endif
    vector<Aln> kernel_alns;
    vector<string> ctg_seqs;
    vector<string> read_seqs;
    shared_ptr<Alns> alns;
    AlnScoring aln_scoring;
    StripedSmithWaterman::Aligner ssw_aligner;
    StripedSmithWaterman::Filter ssw_filter;
    int64_t max_clen;
    int64_t max_rlen;
    int read_group_id;

    AlignBlockData(KmerCtgDHT &kmer_ctg_dht, int read_group_id)
        : aln_scoring(kmer_ctg_dht.aln_scoring)
        , ssw_aligner(kmer_ctg_dht.ssw_aligner)
        , ssw_filter(kmer_ctg_dht.ssw_filter)
        , max_clen(kmer_ctg_dht.max_clen)
        , max_rlen(kmer_ctg_dht.max_rlen)
        , read_group_id(read_group_id)
#ifdef ENABLE_GPUS
        , gpu_driver(kmer_ctg_dht.gpu_driver)
#endif

    {
      DBG_VERBOSE("Created AlignBlockData for kernel ", kmer_ctg_dht.kernel_alns.size(), "\n");
      // copy/swap/reserve necessary data and configs
      kernel_alns.swap(kmer_ctg_dht.kernel_alns);
      kmer_ctg_dht.kernel_alns.reserve(kernel_alns.size());
      ctg_seqs.swap(kmer_ctg_dht.ctg_seqs);
      kmer_ctg_dht.ctg_seqs.reserve(ctg_seqs.size());
      read_seqs.swap(kmer_ctg_dht.read_seqs);
      kmer_ctg_dht.read_seqs.reserve(read_seqs.size());
      alns = make_shared<Alns>();
      alns->reserve(kernel_alns.size());
    }
  };

  static void _ssw_align_block(AlignBlockData &abd, IntermittentTimer &aln_kernel_timer) {
    DBG_VERBOSE("Starting _ssw_align_block of ", abd.kernel_alns.size(), "\n");
    auto alns_ptr = abd.alns.get();
    for (int i = 0; i < abd.kernel_alns.size(); i++) {
      Aln &aln = abd.kernel_alns[i];
      string &cseq = abd.ctg_seqs[i];
      string &rseq = abd.read_seqs[i];
      ssw_align_read(abd.ssw_aligner, abd.ssw_filter, alns_ptr, abd.aln_scoring, aln_kernel_timer, aln, cseq, rseq,
                     abd.read_group_id);
    }
  }

  upcxx::future<> ssw_align_block(IntermittentTimer &aln_kernel_timer, int read_group_id) {
    AsyncTimer t("ssw_align_block (thread)");
    auto &myself = *this;
    shared_ptr<AlignBlockData> sh_abd = make_shared<AlignBlockData>(myself, read_group_id);
    assert(kernel_alns.empty());

    future<> fut = upcxx_utils::execute_in_thread_pool([sh_abd, t, &aln_kernel_timer]() {
      t.start();
      assert(!sh_abd->kernel_alns.empty());
      _ssw_align_block(*sh_abd, aln_kernel_timer);
      t.stop();
    });
    fut = fut.then([&myself, sh_abd, t]() {
      SLOG_VERBOSE("Finished CPU SSW aligning block of ", sh_abd->kernel_alns.size(), " in ", t.get_elapsed(), " s (",
                   (t.get_elapsed() > 0 ? sh_abd->kernel_alns.size() / t.get_elapsed() : 0.0), " aln/s)\n");
      DBG_VERBOSE("appending and returning ", sh_abd->alns->size(), "\n");
      myself.alns->append(*(sh_abd->alns));
    });

    return fut;
  }

#ifdef ENABLE_GPUS

  static void _gpu_align_block_kernel(AlignBlockData &abd, IntermittentTimer &aln_kernel_timer) {
    DBG_VERBOSE("Starting _gpu_align_block_kernel of ", abd.kernel_alns.size(), "\n");
    aln_kernel_timer.start();

    // align query_seqs, ref_seqs, max_query_size, max_ref_size
    abd.gpu_driver.run_kernel_forwards(abd.read_seqs, abd.ctg_seqs, abd.max_rlen, abd.max_clen);
    abd.gpu_driver.kernel_block();
    abd.gpu_driver.run_kernel_backwards(abd.read_seqs, abd.ctg_seqs, abd.max_rlen, abd.max_clen);
    abd.gpu_driver.kernel_block();
    aln_kernel_timer.stop();

    auto aln_results = abd.gpu_driver.get_aln_results();

    for (int i = 0; i < abd.kernel_alns.size(); i++) {
      // progress();
      Aln &aln = abd.kernel_alns[i];
      aln.rstop = aln.rstart + aln_results.query_end[i];
      aln.rstart += aln_results.query_begin[i];
      aln.cstop = aln.cstart + aln_results.ref_end[i];
      aln.cstart += aln_results.ref_begin[i];
      if (aln.orient == '-') switch_orient(aln.rstart, aln.rstop, aln.rlen);
      aln.score1 = aln_results.top_scores[i];
      // FIXME: needs to be set to the second best
      aln.score2 = 0;
      // FIXME: need to get the mismatches
      aln.mismatches = 0;  // ssw_aln.mismatches;
      aln.identity = 100 * aln.score1 / abd.aln_scoring.match / aln.rlen;
      aln.read_group_id = abd.read_group_id;
      // FIXME: need to get cigar
      if (abd.ssw_filter.report_cigar) set_sam_string(aln, "*", "*");  // FIXME until there is a valid:ssw_aln.cigar_string);
      abd.alns->add_aln(aln);
    }
  }

  upcxx::future<> gpu_align_block(IntermittentTimer &aln_kernel_timer, int read_group_id) {
    AsyncTimer t("gpu_align_block (thread)");
    auto &myself = *this;
    shared_ptr<AlignBlockData> sh_abd = make_shared<AlignBlockData>(myself, read_group_id);
    assert(kernel_alns.empty());

    future<> fut = upcxx_utils::execute_in_thread_pool([&myself, t, sh_abd, &aln_kernel_timer] {
      t.start();
      _gpu_align_block_kernel(*sh_abd, aln_kernel_timer);
      t.stop();
    });
    fut = fut.then([&myself, t, sh_abd]() {
      SLOG_VERBOSE("Finished GPU SSW aligning block of ", sh_abd->kernel_alns.size(), " in ", t.get_elapsed(), "s (",
                   (t.get_elapsed() > 0 ? sh_abd->kernel_alns.size() / t.get_elapsed() : 0.0), " aln/s)\n");
      DBG_VERBOSE("appending and returning ", sh_abd->alns->size(), "\n");
      myself.alns->append(*(sh_abd->alns));
    });

    return fut;
  }
#endif

 public:
  unsigned kmer_len;

  // aligner construction: SSW internal defaults are 2 2 3 1

  KmerCtgDHT(int kmer_len, int max_store_size, int max_rpcs_in_flight, Alns &alns, AlnScoring &aln_scoring, int rlen_limit,
             bool compute_cigar, bool use_minimizers, int all_num_ctgs, int ranks_per_gpu = 0)
      : kmer_map({})
      , kmer_store(kmer_map)
      , num_alns(0)
      , num_perfect_alns(0)
      , num_overlaps(0)
      , ssw_aligner(aln_scoring.match, aln_scoring.mismatch, aln_scoring.gap_opening, aln_scoring.gap_extending,
                    aln_scoring.ambiguity)
      , kernel_alns({})
      , ctg_seqs({})
      , read_seqs({})
      , active_kernel_fut(make_future())
      , aln_cpu_bypass_timer("klign.cpp:CPU_BSW-bypass")
      , alns(&alns)
      , kmer_len(kmer_len)
      , use_minimizers(use_minimizers) {
    this->aln_scoring = aln_scoring;
    ssw_filter.report_cigar = compute_cigar;
    kmer_store.set_size("insert ctg seeds", max_store_size, max_rpcs_in_flight);
    kmer_store.set_update_func([](KmerAndCtgLoc<MAX_K> kmer_and_ctg_loc, kmer_map_t &kmer_map) {
      CtgLoc ctg_loc = kmer_and_ctg_loc.ctg_loc;
      const auto it = kmer_map->find(kmer_and_ctg_loc.kmer);
      if (it == kmer_map->end()) {
        kmer_map->insert({kmer_and_ctg_loc.kmer, {false, ctg_loc}});
      } else {
        // in this case, we have a conflict, i.e. the kmer maps to multiple contigs
        it->second.first = true;
        _num_dropped_seed_to_ctgs++;
      }
    });
    gpu_devices = 0;
#ifdef ENABLE_GPUS
    gpu_devices = adept_sw::get_num_node_gpus();
    if (gpu_devices <= 0) {
      // CPU only
      gpu_devices = 0;
    } else {
      if (ranks_per_gpu == 0) {
        // auto detect
        gpu_mem_avail = adept_sw::get_avail_gpu_mem_per_rank(local_team().rank_n(), gpu_devices);
      } else {
        gpu_mem_avail = adept_sw::get_avail_gpu_mem_per_rank(ranks_per_gpu, 1);
      }
      if (gpu_mem_avail) {
        SLOG_VERBOSE("GPU memory available: ", get_size_str(gpu_mem_avail), "\n");
        auto init_time =
            gpu_driver.init(local_team().rank_me(), local_team().rank_n(), (short)aln_scoring.match, (short)-aln_scoring.mismatch,
                            (short)-aln_scoring.gap_opening, (short)-aln_scoring.gap_extending, rlen_limit);
        SLOG_VERBOSE("Initialized adept_sw driver in ", init_time, " s\n");
      } else {
        gpu_devices = 0;
      }
    }
    if (gpu_devices == 0) {
      SWARN("No GPU will be used for alignments");
      gpu_mem_avail = 32 * 1024 * 1024;  // cpu needs a block of memory
    }
#else
    gpu_mem_avail = 32 * 1024 * 1024;  // cpu needs a block of memory
#endif
    ctg_cache.reserve(all_num_ctgs / rank_n());
  }

  void clear() {
    if (kernel_alns.size() || !active_kernel_fut.ready())
      DIE("clear called with alignments in the buffer or active kernel - was flush_remaining called before destrutor?\n");
    clear_aln_bufs();
    aln_cpu_bypass_timer.print_out();
    local_kmer_map_t().swap(*kmer_map);  // release all memory
    kmer_store.clear();
  }

  ~KmerCtgDHT() { clear(); }

  void reserve(int64_t mysize) { kmer_map->reserve(mysize); }
  int64_t size() const { return kmer_map->size(); }

  intrank_t get_target_rank(const Kmer<MAX_K> &kmer, const Kmer<MAX_K> *kmer_rc = nullptr) const {
    if (use_minimizers)
      return kmer.minimizer_hash_fast(MINIMIZER_LEN, kmer_rc) % rank_n();
    else
      return std::hash<Kmer<MAX_K>>{}(kmer) % rank_n();
  }

  int64_t get_num_kmers(bool all = false) {
    if (!all) return reduce_one(kmer_map->size(), op_fast_add, 0).wait();
    return reduce_all(kmer_map->size(), op_fast_add).wait();
  }

  int64_t get_num_perfect_alns(bool all = false) {
    if (!all) return reduce_one(num_perfect_alns, op_fast_add, 0).wait();
    return reduce_all(num_perfect_alns, op_fast_add).wait();
  }

  int64_t get_num_alns(bool all = false) {
    if (!all) return reduce_one(num_alns, op_fast_add, 0).wait();
    return reduce_all(num_alns, op_fast_add).wait();
  }

  int64_t get_num_overlaps(bool all = false) {
    if (!all) return reduce_one(num_overlaps, op_fast_add, 0).wait();
    return reduce_all(num_overlaps, op_fast_add).wait();
  }

  int64_t get_num_dropped_seed_to_ctgs(bool all = false) {
    if (!all) return reduce_one(_num_dropped_seed_to_ctgs, op_fast_add, 0).wait();
    return reduce_all(_num_dropped_seed_to_ctgs, op_fast_add).wait();
  }

  void add_kmer(Kmer<MAX_K> kmer, CtgLoc &ctg_loc) {
    Kmer<MAX_K> kmer_rc = kmer.revcomp();
    ctg_loc.is_rc = false;
    if (kmer_rc < kmer) {
      kmer.swap(kmer_rc);
      ctg_loc.is_rc = true;
    }
    KmerAndCtgLoc<MAX_K> kmer_and_ctg_loc = {kmer, ctg_loc};
    kmer_store.update(get_target_rank(kmer, &kmer_rc), kmer_and_ctg_loc);
  }

  void flush_add_kmers() {
    BarrierTimer timer(__FILEFUNC__, false);  // barrier on exit, not entrance
    kmer_store.flush_updates();
    kmer_store.clear();
  }

  void clear_aln_bufs() {
    kernel_alns.clear();
    ctg_seqs.clear();
    read_seqs.clear();
    max_clen = 0;
    max_rlen = 0;
  }

  void kernel_align_block(IntermittentTimer &aln_kernel_timer, int read_group_id) {
    BaseTimer steal_t("CPU work steal");
    steal_t.start();
    auto num = kernel_alns.size();

    // steal work from this kernel block if the previous kernel is still active
    // if true, this balances the block size that will be sent to the kernel
    while ((gpu_devices == 0 || !active_kernel_fut.ready()) && !kernel_alns.empty()) {
      assert(!ctg_seqs.empty());
      assert(!read_seqs.empty());
#ifndef NO_KLIGN_CPU_WORK_STEAL
      // steal one from the block
      ssw_align_read(aln_cpu_bypass_timer, kernel_alns.back(), ctg_seqs.back(), read_seqs.back(), read_group_id);
      kernel_alns.pop_back();
      ctg_seqs.pop_back();
      read_seqs.pop_back();
#endif
#ifndef ENABLE_GPUS
      std::this_thread::yield();  // yield if the kernel is CPU based
#endif
      progress();
    }

    steal_t.stop();
    auto steal_secs = steal_t.get_elapsed();
    if (num != kernel_alns.size()) {
      auto num_stole = num - kernel_alns.size();
      LOG("Stole from kernel block ", num_stole, " alignments in ", steal_secs, "s (",
          (steal_secs > 0 ? num_stole / steal_secs : 0.0), " aln/s), while waiting for previous block to complete",
          (kernel_alns.empty() ? " - THE ENTIRE BLOCK" : ""), "\n");
    } else if (steal_secs > 0.01) {
      LOG("Waited ", steal_secs, "s for previous block to complete\n");
    }
    if (!kernel_alns.empty()) {
      assert(active_kernel_fut.ready() && "active_kernel_fut should already be ready");
      active_kernel_fut.wait();  // should be ready already
#ifdef ENABLE_GPUS
      // for now, the GPU alignment doesn't support cigars
      if (!ssw_filter.report_cigar && gpu_devices > 0) {
        active_kernel_fut = gpu_align_block(aln_kernel_timer, read_group_id);
      } else {
#ifdef __PPC64__
        SWARN("FIXME Issue #49,#60 no cigars for gpu alignments\n");
        active_kernel_fut = gpu_align_block(aln_kernel_timer, read_group_id);
#else
        active_kernel_fut = ssw_align_block(aln_kernel_timer, read_group_id);
#endif
      }
#else
      active_kernel_fut = ssw_align_block(aln_kernel_timer, read_group_id);
#endif
    }

    clear_aln_bufs();
  }

  void flush_remaining(IntermittentTimer &aln_kernel_timer, int read_group_id) {
    BaseTimer t(__FILEFUNC__);
    t.start();
    auto num = kernel_alns.size();
    if (num) {
      kernel_align_block(aln_kernel_timer, read_group_id);
    }
    bool is_ready = active_kernel_fut.ready();
    active_kernel_fut.wait();
    t.stop();
    if (num || !is_ready) {
      SLOG_VERBOSE("Aligned and waited for final block with ", num, " alignments in ", t.get_elapsed(), "\n");
    }
  }

  future<vector<KmerAndCtgLoc<MAX_K>>> get_ctgs_with_kmers(int target_rank, vector<Kmer<MAX_K>> &kmers) {
    return rpc(target_rank,
               [](vector<Kmer<MAX_K>> kmers, kmer_map_t &kmer_map) {
                 vector<KmerAndCtgLoc<MAX_K>> kmer_ctg_locs;
                 kmer_ctg_locs.reserve(kmers.size());
                 for (auto &kmer : kmers) {
                   const auto it = kmer_map->find(kmer);
                   if (it == kmer_map->end()) continue;
                   // skip conflicts
                   if (it->second.first) continue;
                   // now add it
                   kmer_ctg_locs.push_back({kmer, it->second.second});
                 }
                 return kmer_ctg_locs;
               },
               kmers, kmer_map);
  }

#ifdef DEBUG
  void dump_ctg_kmers() {
    BarrierTimer timer(__FILEFUNC__, false);  // barrier on exit not entrance
    string dump_fname = "ctg_kmers-" + to_string(kmer_len) + ".txt.gz";
    get_rank_path(dump_fname, rank_me());
    zstr::ofstream dump_file(dump_fname);
    ostringstream out_buf;
    ProgressBar progbar(kmer_map->size(), "Dumping kmers to " + dump_fname);
    int64_t i = 0;
    for (auto &elem : *kmer_map) {
      // FIXME this was broken when I got here.... made my best guess as to what the fields are supposed to be. -Rob
      auto &pair_ctg_loc = elem.second;
      auto &ctg_loc = pair_ctg_loc.second;
      out_buf << elem.first << " " << ctg_loc.cid << " " << ctg_loc.clen << " " << ctg_loc.pos_in_ctg << " " << ctg_loc.is_rc
              << "\n";
      i++;
      if (!(i % 1000)) {
        dump_file << out_buf.str();
        out_buf = ostringstream();
      }
      progbar.update();
    }
    if (!out_buf.str().empty()) dump_file << out_buf.str();
    dump_file.close();
    progbar.done();
    SLOG_VERBOSE("Dumped ", this->get_num_kmers(), " kmers\n");
  }
#endif

  void compute_alns_for_read(HASH_TABLE<cid_t, ReadAndCtgLoc> *aligned_ctgs_map, const string &rname, string rseq,
                             int read_group_id, IntermittentTimer &fetch_ctg_seqs_timer, IntermittentTimer &aln_kernel_timer) {
    int rlen = rseq.length();
    string rseq_rc = revcomp(rseq);
    // make the buffer pretty big, but expand in the loop if it's too small for any one contig
    int buf_len = 1000000;
    char *seq_buf = new char[2 * rlen + 1];
    for (auto &elem : *aligned_ctgs_map) {
      progress();
      int pos_in_read = elem.second.pos_in_read;
      bool read_kmer_is_rc = elem.second.read_is_rc;
      CtgLoc ctg_loc = elem.second.ctg_loc;
      char orient = '+';
      string *rseq_ptr = &rseq;
      if (ctg_loc.is_rc != read_kmer_is_rc) {
        // it's revcomp in either contig or read, but not in both or neither
        orient = '-';
        pos_in_read = rlen - (kmer_len + pos_in_read);
        rseq_ptr = &rseq_rc;
      }
      // calculate available bases before and after the seeded kmer
      int ctg_bases_left_of_kmer = ctg_loc.pos_in_ctg;
      int ctg_bases_right_of_kmer = ctg_loc.clen - ctg_bases_left_of_kmer - kmer_len;
      int read_bases_left_of_kmer = pos_in_read;
      int read_bases_right_of_kmer = rlen - kmer_len - pos_in_read;
      int left_of_kmer = min(read_bases_left_of_kmer, ctg_bases_left_of_kmer);
      int right_of_kmer = min(read_bases_right_of_kmer, ctg_bases_right_of_kmer);

      int cstart = ctg_loc.pos_in_ctg - left_of_kmer;
      int rstart = pos_in_read - left_of_kmer;
      int overlap_len = left_of_kmer + kmer_len + right_of_kmer;

      // use the whole read, to account for possible indels
      string read_subseq = rseq_ptr->substr(0, rlen);

      assert(cstart >= 0 && cstart + overlap_len <= ctg_loc.clen);
      assert(overlap_len <= 2 * rlen);

      string ctg_subseq;
      bool found = false;
      auto it = ctg_cache.find(ctg_loc.cid);
      if (it != ctg_cache.end()) {
        found = true;
        ctg_subseq.resize(overlap_len);
        for (int i = 0; i < overlap_len; i++) {
          if (it->second[i + cstart] == ' ') {
            found = false;
            break;
          }
          ctg_subseq[i] = it->second[i + cstart];
        }
      } else {
        it = ctg_cache.insert({ctg_loc.cid, string(ctg_loc.clen, ' ')}).first;
      }
      if (!found) {
        fetch_ctg_seqs_timer.start();
        rget(ctg_loc.seq_gptr + cstart, seq_buf, overlap_len).wait();
        fetch_ctg_seqs_timer.stop();
        ctg_bytes_fetched += overlap_len;
        ctg_subseq = string(seq_buf, overlap_len);
        for (int i = 0; i < overlap_len; i++) {
          it->second[i + cstart] = ctg_subseq[i];
        }
      }
      align_read(rname, ctg_loc.cid, read_subseq, ctg_subseq, rstart, rlen, cstart, ctg_loc.clen, orient, overlap_len,
                 read_group_id, aln_kernel_timer);
      num_alns++;
    }
    delete[] seq_buf;
  }

  void sort_alns() {
    if (!kernel_alns.empty()) {
      DIE("sort_alns called while alignments are still pending to be processed - ", kernel_alns.size());
    }
    if (!active_kernel_fut.ready()) {
      SWARN("Waiting for active_kernel - has flush_remaining() been called?\n");
    }
    active_kernel_fut.wait();
    alns->sort_alns().wait();
  }

  int get_gpu_mem_avail() { return gpu_mem_avail; }

  void log_ctg_bytes_fetched() {
    auto all_ctg_bytes_fetched = reduce_one(ctg_bytes_fetched, op_fast_add, 0).wait();
    auto max_ctg_bytes_fetched = reduce_one(ctg_bytes_fetched, op_fast_max, 0).wait();
    SLOG_VERBOSE("Contig bytes fetched ", get_size_str(all_ctg_bytes_fetched), " balance ",
         (double)all_ctg_bytes_fetched / (rank_n() * max_ctg_bytes_fetched), "\n");
  }
};

template <int MAX_K>
static void build_alignment_index(KmerCtgDHT<MAX_K> &kmer_ctg_dht, Contigs &ctgs, unsigned min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
  // estimate and reserve room in the local map to avoid excessive reallocations
  int64_t est_num_kmers = 0;
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    auto len = ctg->seq.length();
    if (len < min_ctg_len) continue;
    est_num_kmers += len - kmer_ctg_dht.kmer_len + 1;
  }
  est_num_kmers = upcxx::reduce_all(est_num_kmers, upcxx::op_fast_add).wait();
  auto my_reserve = 1.2 * est_num_kmers / rank_n() + 2000;  // 120% to keep the map fast
  kmer_ctg_dht.reserve(my_reserve);
  vector<Kmer<MAX_K>> kmers;
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() < min_ctg_len) continue;
    global_ptr<char> seq_gptr = allocate<char>(ctg->seq.length() + 1);
    strcpy(seq_gptr.local(), ctg->seq.c_str());
    CtgLoc ctg_loc = {.cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length()};
    Kmer<MAX_K>::get_kmers(kmer_ctg_dht.kmer_len, ctg->seq, kmers);
    num_kmers += kmers.size();
    for (unsigned i = 0; i < kmers.size(); i++) {
      ctg_loc.pos_in_ctg = i;
      kmer_ctg_dht.add_kmer(kmers[i], ctg_loc);
    }
    progress();
  }
  auto fut = progbar.set_done();
  kmer_ctg_dht.flush_add_kmers();
  auto tot_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  fut.wait();
  auto num_kmers_in_ht = kmer_ctg_dht.get_num_kmers();
  LOG("Estimated room for ", my_reserve, " my final count ", kmer_ctg_dht.size(), "\n");
  SLOG_VERBOSE("Processed ", tot_num_kmers, " seeds from contigs, added ", num_kmers_in_ht, "\n");
  auto num_dropped_seed_to_ctgs = kmer_ctg_dht.get_num_dropped_seed_to_ctgs();
  if (num_dropped_seed_to_ctgs)
    SLOG_VERBOSE("Dropped ", num_dropped_seed_to_ctgs, " non-unique seed-to-contig mappings (", setprecision(2), fixed,
                 (100.0 * num_dropped_seed_to_ctgs / tot_num_kmers), "%)\n");
}

struct ReadRecord {
  string id;
  string seq;
  string quals;

  HASH_TABLE<cid_t, ReadAndCtgLoc> aligned_ctgs_map;

  ReadRecord(const string &id, const string &seq, const string &quals)
      : id(id)
      , seq(seq)
      , quals(quals) {}
};

struct KmerToRead {
  ReadRecord *read_record;
  int pos_in_read;
  bool is_rc;
};

template <int MAX_K>
static int align_kmers(KmerCtgDHT<MAX_K> &kmer_ctg_dht, HASH_TABLE<Kmer<MAX_K>, vector<KmerToRead>> &kmer_read_map,
                       vector<ReadRecord *> &read_records, HASH_TABLE<Kmer<MAX_K>, KmerAndCtgLoc<MAX_K>> &kmer_cache,
                       int64_t &kmer_cache_hits, IntermittentTimer &compute_alns_timer, IntermittentTimer &get_ctgs_timer,
                       IntermittentTimer &fetch_ctg_seqs_timer, IntermittentTimer &aln_kernel_timer, int64_t &num_excess_alns_reads,
                       int &read_group_id, int64_t &kmer_bytes_sent, int64_t &kmer_bytes_received) {
  auto process_kmer_ctg_loc = [](HASH_TABLE<Kmer<MAX_K>, vector<KmerToRead>> &kmer_read_map, int64_t &num_excess_alns_reads,
                                 int64_t &kmer_bytes_received, const KmerAndCtgLoc<MAX_K> &kmer_ctg_loc) {
    kmer_bytes_received += sizeof(kmer_ctg_loc);
    // get the reads that this kmer mapped to
    auto kmer_read_map_it = kmer_read_map.find(kmer_ctg_loc.kmer);
    if (kmer_read_map_it == kmer_read_map.end()) DIE("Could not find kmer ", kmer_ctg_loc.kmer);
    // this is a list of the reads
    auto &kmer_to_reads = kmer_read_map_it->second;
    // now add the ctg loc to all the reads
    for (auto &kmer_to_read : kmer_to_reads) {
      auto read_record = kmer_to_read.read_record;
      int pos_in_read = kmer_to_read.pos_in_read;
      bool read_is_rc = kmer_to_read.is_rc;
      if (KLIGN_MAX_ALNS_PER_READ && read_record->aligned_ctgs_map.size() >= KLIGN_MAX_ALNS_PER_READ) {
        // too many mappings for this read, stop adding to it
        num_excess_alns_reads++;
        continue;
      }
      // this here ensures that we don't insert duplicate mappings
      read_record->aligned_ctgs_map.insert({kmer_ctg_loc.ctg_loc.cid, {pos_in_read, read_is_rc, kmer_ctg_loc.ctg_loc}});
    }
  };

  // get the contigs that match one read
  // extract a list of kmers for each target rank
  auto kmer_lists = new vector<Kmer<MAX_K>>[rank_n()];
  for (auto &elem : kmer_read_map) {
    auto kmer = elem.first;
    // FIXME: when creating this list, look for kmers in cache.
    // If found, don't add to the list instead add to read_record
    auto it = kmer_cache.find(kmer);
    if (it == kmer_cache.end()) {
      kmer_lists[kmer_ctg_dht.get_target_rank(kmer)].push_back(kmer);
    } else {
      progress();
      process_kmer_ctg_loc(kmer_read_map, num_excess_alns_reads, kmer_bytes_received, it->second);
      kmer_cache_hits++;
    }
  }
  get_ctgs_timer.start();
  future<> fut_serial_results = make_future();
  // fetch ctgs for each set of kmers from target ranks
  auto lranks = local_team().rank_n();
  auto nnodes = rank_n() / lranks;
  for (auto target_rank : upcxx_utils::foreach_rank_by_node()) {  // stagger by rank_me, round robin by node
    progress();
    // skip targets that have no ctgs - this should reduce communication at scale
    if (kmer_lists[target_rank].empty()) continue;
    kmer_bytes_sent += kmer_lists[target_rank].size() * sizeof(Kmer<MAX_K>);
    auto fut_get_ctgs = kmer_ctg_dht.get_ctgs_with_kmers(target_rank, kmer_lists[target_rank]);
    auto fut_rpc_returned = fut_get_ctgs.then([&](const vector<KmerAndCtgLoc<MAX_K>> kmer_ctg_locs) {
      // iterate through the kmers, each one has an associated ctg location
      for (auto &kmer_ctg_loc : kmer_ctg_locs) {
        process_kmer_ctg_loc(kmer_read_map, num_excess_alns_reads, kmer_bytes_received, kmer_ctg_loc);
        // now cache it
        kmer_cache.insert({kmer_ctg_loc.kmer, kmer_ctg_loc});
      }
    });

    upcxx_utils::limit_outstanding_futures(fut_rpc_returned, std::max(nnodes * 2, lranks * 4)).wait();
  }

  upcxx_utils::flush_outstanding_futures();

  get_ctgs_timer.stop();
  delete[] kmer_lists;
  kmer_read_map.clear();

  compute_alns_timer.start();
  int num_reads_aligned = 0;
  // create a new list of records with all reads having < KLIGN_MAX_ALNS_PER_READ, i.e. those without excessive mappings
  // setting KLIGN_MAX_ALNS_PER_READ to zero means don't drop any
  vector<ReadRecord *> good_read_records;
  good_read_records.reserve(read_records.size());
  for (auto read_record : read_records) {
    if (!KLIGN_MAX_ALNS_PER_READ || read_record->aligned_ctgs_map.size() < KLIGN_MAX_ALNS_PER_READ)
      good_read_records.push_back(read_record);
  }
  // compute alignments for each read
  for (auto read_record : good_read_records) {
    progress();
    // compute alignments
    if (read_record->aligned_ctgs_map.size()) {
      num_reads_aligned++;
      kmer_ctg_dht.compute_alns_for_read(&read_record->aligned_ctgs_map, read_record->id, read_record->seq, read_group_id,
                                         fetch_ctg_seqs_timer, aln_kernel_timer);
    }
    delete read_record;
  }
  read_records.clear();
  compute_alns_timer.stop();
  return num_reads_aligned;
}

template <int MAX_K>
static double do_alignments(KmerCtgDHT<MAX_K> &kmer_ctg_dht, vector<PackedReads *> &packed_reads_list, int seed_space,
                            bool compute_cigar) {
  BarrierTimer timer(__FILEFUNC__);
  SLOG_VERBOSE("Using a seed space of ", seed_space, "\n");
  int64_t tot_num_kmers = 0;
  int64_t num_reads = 0;
  int64_t num_reads_aligned = 0, num_excess_alns_reads = 0;
  IntermittentTimer compute_alns_timer(__FILENAME__ + string(":") + "Compute alns");
  IntermittentTimer get_ctgs_timer(__FILENAME__ + string(":") + "Get ctgs with kmer");
  IntermittentTimer fetch_ctg_seqs_timer(__FILENAME__ + string(":") + "Fetch ctg seqs");
#ifdef ENABLE_GPUS
  IntermittentTimer aln_kernel_timer(__FILENAME__ + string(":") +
                                     ((compute_cigar && kmer_ctg_dht.get_gpu_mem_avail()) ? "SSW" : "GPU_BSW"));
#else
  IntermittentTimer aln_kernel_timer(__FILENAME__ + string(":") + "SSW");
#endif
  kmer_ctg_dht.clear_aln_bufs();
  barrier();
  int64_t kmer_bytes_received = 0;
  int64_t kmer_bytes_sent = 0;
  upcxx::future<> all_done = make_future();
  int read_group_id = 0;
  HASH_TABLE<Kmer<MAX_K>, KmerAndCtgLoc<MAX_K>> kmer_cache;
  int64_t kmer_cache_hits = 0;
  int64_t num_lookups = 0;
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    string read_id, read_seq, quals;
    ProgressBar progbar(packed_reads->get_local_num_reads(), "Aligning reads to contigs");
    vector<ReadRecord *> read_records;
    HASH_TABLE<Kmer<MAX_K>, vector<KmerToRead>> kmer_read_map;
    vector<Kmer<MAX_K>> kmers;
    while (true) {
      progress();
      if (!packed_reads->get_next_read(read_id, read_seq, quals)) break;
      progbar.update();
      // this happens when a placeholder read with just a single N character is added after merging reads
      if (kmer_ctg_dht.kmer_len > read_seq.length()) continue;
      Kmer<MAX_K>::get_kmers(kmer_ctg_dht.kmer_len, read_seq, kmers);
      tot_num_kmers += kmers.size();
      ReadRecord *read_record = new ReadRecord(read_id, read_seq, quals);
      read_records.push_back(read_record);
      bool filled = false;
      for (int i = 0; i < (int)kmers.size(); i += seed_space) {
        Kmer<MAX_K> kmer = kmers[i];
        Kmer<MAX_K> kmer_rc = kmer.revcomp();
        bool is_rc = false;
        if (kmer_rc < kmer) {
          kmer = kmer_rc;
          is_rc = true;
        }
        auto it = kmer_read_map.find(kmer);
        if (it == kmer_read_map.end()) it = kmer_read_map.insert({kmer, {}}).first;
        it->second.push_back({read_record, i, is_rc});
        if (kmer_read_map.size() >= KLIGN_CTG_FETCH_BUF_SIZE) filled = true;
      }
      if (filled) {
        num_lookups += kmer_read_map.size();
        num_reads_aligned += align_kmers(kmer_ctg_dht, kmer_read_map, read_records, kmer_cache, kmer_cache_hits, compute_alns_timer,
                                         get_ctgs_timer, fetch_ctg_seqs_timer, aln_kernel_timer, num_excess_alns_reads,
                                         read_group_id, kmer_bytes_sent, kmer_bytes_received);
      }
      num_reads++;
    }
    if (read_records.size()) {
      num_lookups += kmer_read_map.size();
      num_reads_aligned += align_kmers(kmer_ctg_dht, kmer_read_map, read_records, kmer_cache, kmer_cache_hits, compute_alns_timer,
                                       get_ctgs_timer, fetch_ctg_seqs_timer, aln_kernel_timer, num_excess_alns_reads, read_group_id,
                                       kmer_bytes_sent, kmer_bytes_received);
    }
    kmer_ctg_dht.flush_remaining(aln_kernel_timer, read_group_id);
    read_group_id++;
    all_done = when_all(all_done, progbar.set_done());
  }
  kmer_cache.clear();
  HASH_TABLE<Kmer<MAX_K>, KmerAndCtgLoc<MAX_K>>().swap(kmer_cache);
  auto all_kmer_cache_hits = reduce_one(kmer_cache_hits, op_fast_add, 0).wait();
  auto all_lookups = reduce_one(num_lookups, op_fast_add, 0).wait();

  SLOG_VERBOSE("Got ", perc_str(all_kmer_cache_hits, all_lookups), " cache hits for ", all_lookups, " kmer lookups\n");

  // make sure to do any outstanding kernel block alignments

  auto tot_num_reads_fut = reduce_one(num_reads, op_fast_add, 0);
  auto num_excess_alns_reads_fut = reduce_one(num_excess_alns_reads, op_fast_add, 0);
  auto num_seeds_fut = reduce_one(tot_num_kmers, op_fast_add, 0);
  auto tot_num_reads_aligned_fut = reduce_one(num_reads_aligned, op_fast_add, 0);

  kmer_ctg_dht.sort_alns();
  auto num_overlaps = kmer_ctg_dht.get_num_overlaps();
  all_done.wait();
  barrier();

  auto tot_num_reads = tot_num_reads_fut.wait();
  SLOG_VERBOSE("Parsed ", tot_num_reads, " reads, with ", num_seeds_fut.wait(), " seeds\n");
  auto tot_num_alns = kmer_ctg_dht.get_num_alns();
  SLOG_VERBOSE("Found ", tot_num_alns, " alignments of which ", perc_str(kmer_ctg_dht.get_num_perfect_alns(), tot_num_alns),
               " are perfect\n");
  auto tot_excess_alns_reads = num_excess_alns_reads_fut.wait();
  if (num_excess_alns_reads)
    SLOG_VERBOSE("Dropped ", tot_excess_alns_reads, " reads because of alignments in excess of ", KLIGN_MAX_ALNS_PER_READ, "\n");
  if (num_overlaps) SLOG_VERBOSE("Dropped ", perc_str(num_overlaps, tot_num_alns), " alignments becasue of overlaps\n");
  auto tot_num_reads_aligned = tot_num_reads_aligned_fut.wait();
  SLOG_VERBOSE("Mapped ", perc_str(tot_num_reads_aligned, tot_num_reads), " reads to contigs\n");
  SLOG_VERBOSE("Average mappings per read ", (double)tot_num_alns / tot_num_reads_aligned, "\n");
  auto all_kmer_bytes_sent = reduce_one(kmer_bytes_sent, op_fast_add, 0).wait();
  auto all_kmer_bytes_received = reduce_one(kmer_bytes_received, op_fast_add, 0).wait();

  SLOG_VERBOSE("Sent ", get_size_str(all_kmer_bytes_sent), " and received ", get_size_str(all_kmer_bytes_received), " of kmers\n");

  kmer_ctg_dht.log_ctg_bytes_fetched();

  get_ctgs_timer.done_all();
  fetch_ctg_seqs_timer.done_all();
  compute_alns_timer.done_all();
  double aln_kernel_secs = aln_kernel_timer.get_elapsed();
  aln_kernel_timer.done_all();
  return aln_kernel_secs;
}

template <int MAX_K>
double find_alignments(unsigned kmer_len, vector<PackedReads *> &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                       Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit, bool use_minimizers, bool compute_cigar,
                       int min_ctg_len, int ranks_per_gpu) {
  BarrierTimer timer(__FILEFUNC__);
  _num_dropped_seed_to_ctgs = 0;
  Kmer<MAX_K>::set_k(kmer_len);
  SLOG_VERBOSE("Aligning with seed size of ", kmer_len, "\n");
  // default for normal alignments in the pipeline, but for final alignments, uses minimap2 defaults
  AlnScoring aln_scoring = {.match = ALN_MATCH_SCORE,
                            .mismatch = ALN_MISMATCH_COST,
                            .gap_opening = ALN_GAP_OPENING_COST,
                            .gap_extending = ALN_GAP_EXTENDING_COST,
                            .ambiguity = ALN_AMBIGUITY_COST};
  if (compute_cigar) {
    AlnScoring alt_aln_scoring = {.match = 2, .mismatch = 4, .gap_opening = 4, .gap_extending = 2, .ambiguity = 1};
    aln_scoring = alt_aln_scoring;
  }
  auto all_num_ctgs = reduce_all(ctgs.size(), op_fast_add).wait();
  SLOG_VERBOSE("Alignment scoring parameters: ", aln_scoring.to_string(), "\n");
  KmerCtgDHT<MAX_K> kmer_ctg_dht(kmer_len, max_store_size, max_rpcs_in_flight, alns, aln_scoring, rlen_limit, compute_cigar,
                                 use_minimizers, all_num_ctgs, ranks_per_gpu);
  barrier();
  build_alignment_index(kmer_ctg_dht, ctgs, min_ctg_len);
#ifdef DEBUG
// kmer_ctg_dht.dump_ctg_kmers();
#endif
  double kernel_elapsed = do_alignments(kmer_ctg_dht, packed_reads_list, seed_space, compute_cigar);
  barrier();
  auto num_alns = kmer_ctg_dht.get_num_alns();
  auto num_dups = alns.get_num_dups();
  if (num_dups) SLOG_VERBOSE("Number of duplicate alignments ", perc_str(num_dups, num_alns), "\n");
  barrier();
  return kernel_elapsed;
}
