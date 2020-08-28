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
#include <algorithm>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/flat_aggr_store.hpp"

#include "utils.hpp"
#include "ssw.hpp"
#include "contigs.hpp"
#include "kmer.hpp"
#include "alignments.hpp"
#include "packed_reads.hpp"
#ifdef ENABLE_GPUS
#include "adept-sw/driver.hpp"
#endif

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

struct AlnScoring {
  int match, mismatch, gap_opening, gap_extending, ambiguity;

  string to_string() {
    ostringstream oss;
    oss << "match " << match << " mismatch " << mismatch << " gap open " << gap_opening << " gap extend " << gap_extending
        << " ambiguity " << ambiguity;
    return oss.str();
  }
};

using cid_t = int64_t;

struct CtgLoc {
  cid_t cid;
  global_ptr<char> seq_gptr;
  int clen;
  int pos_in_ctg;
  bool is_rc;
};


struct ReadAndCtgLoc {
  int pos_in_read;
  bool read_is_rc;
  CtgLoc ctg_loc;
};

template<int MAX_K>
struct KmerCtgLoc {
  Kmer<MAX_K> kmer;
  CtgLoc ctg_loc;
  UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
};


// global variables to avoid passing dist objs to rpcs
static int64_t _num_dropped_seed_to_ctgs = 0;

template <int MAX_K>
class KmerCtgDHT {
  struct KmerAndCtgLoc {
    Kmer<MAX_K> kmer;
    CtgLoc ctg_loc;
    UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
  };

  // maps a kmer to a contig - the first part of the pair is set to true if this is a conflict,
  // with a kmer mapping to multiple contigs
  using kmer_map_t = dist_object<HASH_TABLE<Kmer<MAX_K>, pair<bool, CtgLoc>>>;
  kmer_map_t kmer_map;

  FlatAggrStore<KmerAndCtgLoc, kmer_map_t &> kmer_store;

  int64_t num_alns;
  int64_t num_perfect_alns;
  int64_t num_overlaps;

  int64_t ctg_seq_bytes_fetched;

  Alns *alns;

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

  int64_t max_clen = 0;
  int64_t max_rlen = 0;
  size_t gpu_mem_avail = 0;

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

  void set_sam_string(Aln &aln, string read_seq, string cigar) {
    aln.sam_string = aln.read_id + "\t";
    if (aln.orient == '-') {
      aln.sam_string += "16\t";
      read_seq = revcomp(read_seq);
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
    aln.sam_string += "AS:i:" + to_string(aln.score1) + "\tNM:i:" + to_string(aln.mismatches);
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
                  int clen, char orient, int overlap_len, IntermittentTimer &aln_kernel_timer) {
    if (cseq.compare(0, overlap_len, rseq, rstart, overlap_len) == 0) {
      num_perfect_alns++;
      int rstop = rstart + overlap_len;
      int cstop = cstart + overlap_len;
      if (orient == '-') switch_orient(rstart, rstop, rlen);
      int score1 = overlap_len * aln_scoring.match;
      int identity = 100 * score1 / aln_scoring.match / rlen;
      Aln aln = {rname, cid, rstart, rstop, rlen, cstart, cstop, clen, orient, score1, 0, identity, 0};
      if (ssw_filter.report_cigar) set_sam_string(aln, rseq, to_string(overlap_len) + "M");
      alns->add_aln(aln);
    } else {
      max_clen = max((int64_t)cseq.size(), max_clen);
      max_rlen = max((int64_t)rseq.size(), max_rlen);
      int64_t num_alns = kernel_alns.size() + 1;
      unsigned max_matrix_size = (max_clen + 1) * (max_rlen + 1);
      int64_t tot_mem_est = num_alns * (max_clen + max_rlen + 2 * sizeof(int) + 5 * sizeof(short));
      // contig is the ref, read is the query - done this way so that we can potentially do multiple alns to each read
      // this is also the way it's done in meraligner
      kernel_alns.push_back({rname, cid, 0, 0, rlen, cstart, 0, clen, orient});
      ctg_seqs.push_back(cseq);
      read_seqs.push_back(rseq);
      if ((tot_mem_est >= gpu_mem_avail && kernel_alns.size()) || kernel_alns.size() >= KLIGN_GPU_BLOCK_SIZE) {
        DBG("tot_mem_est (", tot_mem_est, ") >= gpu_mem_avail (", gpu_mem_avail, " - dispatching ", kernel_alns.size(),
            " alignments\n");
        kernel_align_block(aln_kernel_timer);
      }
    }
  }

  void ssw_align_block(IntermittentTimer &aln_kernel_timer) {
    StripedSmithWaterman::Alignment ssw_aln;
    for (int i = 0; i < kernel_alns.size(); i++) {
      Aln &aln = kernel_alns[i];
      string &cseq = ctg_seqs[i];
      string &rseq = read_seqs[i];
      // make sure upcxx progress is done before starting alignment
      // discharge is for internal progress;  progress defaults to user level
      progress();
      discharge();
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
      if (ssw_filter.report_cigar) set_sam_string(aln, rseq, ssw_aln.cigar_string);
      alns->add_aln(aln);
    }
  }

#ifdef ENABLE_GPUS
  void gpu_align_block(IntermittentTimer &aln_kernel_timer) {
    progress();
    discharge();

    aln_kernel_timer.start();
    // align query_seqs, ref_seqs, max_query_size, max_ref_size
    gpu_driver.run_kernel_forwards(read_seqs, ctg_seqs, max_rlen, max_clen);
    while (!gpu_driver.kernel_is_done()) progress();
    gpu_driver.run_kernel_backwards(read_seqs, ctg_seqs, max_rlen, max_clen);
    while (!gpu_driver.kernel_is_done()) progress();
    aln_kernel_timer.stop();

    auto aln_results = gpu_driver.get_aln_results();

    for (int i = 0; i < kernel_alns.size(); i++) {
      progress();
      Aln &aln = kernel_alns[i];
      aln.rstop = aln.rstart + aln_results.query_end[i];
      aln.rstart += aln_results.query_begin[i];
      aln.cstop = aln.cstart + aln_results.ref_end[i];
      aln.cstart += aln_results.ref_begin[i];
      if (aln.orient == '-') switch_orient(aln.rstart, aln.rstop, aln.rlen);
      aln.score1 = aln_results.top_scores[i];
      // FIXME: needs to be set to the second best
      aln.score2 = 0;
      // FIXME: need to get the mismatches
      aln.mismatches = 0;//ssw_aln.mismatches;
      aln.identity = 100 * aln.score1 / aln_scoring.match / aln.rlen;
      // FIXME: need to get cigar
      //if (ssw_filter.report_cigar) set_sam_string(aln, rseq, ssw_aln.cigar_string);
      alns->add_aln(aln);
    }
  }
#endif

 public:
  int kmer_len;

  // aligner construction: SSW internal defaults are 2 2 3 1

  KmerCtgDHT(int kmer_len, int max_store_size, int max_rpcs_in_flight, Alns &alns, AlnScoring &aln_scoring, int rlen_limit,
             bool compute_cigar, int ranks_per_gpu = 0)
      : kmer_map({})
      , kmer_store(kmer_map)
      , ssw_aligner(aln_scoring.match, aln_scoring.mismatch, aln_scoring.gap_opening, aln_scoring.gap_extending,
                    aln_scoring.ambiguity)
      , num_alns(0)
      , num_perfect_alns(0)
      , num_overlaps(0)
      , kmer_len(kmer_len)
      , ctg_seq_bytes_fetched(0)
      , alns(&alns)
      , kernel_alns({})
      , ctg_seqs({})
      , read_seqs({}) {
    this->aln_scoring = aln_scoring;
    ssw_filter.report_cigar = compute_cigar;
    kmer_store.set_size("insert ctg seeds", max_store_size, max_rpcs_in_flight);
    kmer_store.set_update_func([](KmerAndCtgLoc kmer_and_ctg_loc, kmer_map_t &kmer_map) {
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
#ifdef ENABLE_GPUS
    auto device_count = adept_sw::get_num_node_gpus();
    if (ranks_per_gpu == 0) {
        // auto detect
        gpu_mem_avail = adept_sw::get_avail_gpu_mem_per_rank(local_team().rank_n(), device_count);
    } else {
        gpu_mem_avail = adept_sw::get_avail_gpu_mem_per_rank(ranks_per_gpu, 1);
    }
    if (gpu_mem_avail) {
      SLOG_VERBOSE("GPU memory available: ", get_size_str(gpu_mem_avail), "\n");
      auto init_time = gpu_driver.init(local_team().rank_me(), local_team().rank_n(), (short)aln_scoring.match, (short)-aln_scoring.mismatch,
                      (short)-aln_scoring.gap_opening, (short)-aln_scoring.gap_extending, rlen_limit);
      SLOG_VERBOSE("Initialized adept_sw driver in ", init_time, " s\n");
      //} else {
      //SWARN("No GPU memory detected!  (", device_count, " devices detected)\n");
    }
#else
    gpu_mem_avail = 32 * 1024 * 1024;
#endif
  }

  void clear() {
    for (auto it = kmer_map->begin(); it != kmer_map->end();) {
      it = kmer_map->erase(it);
    }
    kmer_store.clear();
  }

  ~KmerCtgDHT() { clear(); }

  intrank_t get_target_rank(Kmer<MAX_K> &kmer) { return std::hash<Kmer<MAX_K>>{}(kmer) % rank_n(); }

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

  int64_t get_ctg_seq_bytes_fetched() { return reduce_one(ctg_seq_bytes_fetched, op_fast_add, 0).wait(); }

  void add_kmer(Kmer<MAX_K> kmer, CtgLoc &ctg_loc) {
    Kmer<MAX_K> kmer_rc = kmer.revcomp();
    ctg_loc.is_rc = false;
    if (kmer_rc < kmer) {
      kmer = kmer_rc;
      ctg_loc.is_rc = true;
    }
    KmerAndCtgLoc kmer_and_ctg_loc = {kmer, ctg_loc};
    kmer_store.update(get_target_rank(kmer), kmer_and_ctg_loc);
  }

  void flush_add_kmers() {
    BarrierTimer timer(__FILEFUNC__);
    kmer_store.flush_updates();
  }

  void clear_aln_bufs() {
    kernel_alns.clear();
    ctg_seqs.clear();
    read_seqs.clear();
    max_clen = 0;
    max_rlen = 0;
  }

  void kernel_align_block(IntermittentTimer &aln_kernel_timer) {
    BaseTimer t(__FILEFUNC__);
    t.start();
#ifdef ENABLE_GPUS
    // for now, the GPU alignment doesn't support cigars
    if (!ssw_filter.report_cigar && gpu_mem_avail) gpu_align_block(aln_kernel_timer);
    else ssw_align_block(aln_kernel_timer);
#else
    ssw_align_block(aln_kernel_timer);
#endif
    t.stop();
    // interferes with progress ticker
    DBG("Aligned block with ", kernel_alns.size(), " alignments in ", t.get_elapsed(), "\n");
    clear_aln_bufs();
  }

  future<vector<KmerCtgLoc<MAX_K>>> get_ctgs_with_kmers(int target_rank, vector<Kmer<MAX_K>> &kmers) {
    return rpc(
        target_rank,
        [](vector<Kmer<MAX_K>> kmers, kmer_map_t &kmer_map) {
          vector<KmerCtgLoc<MAX_K>> kmer_ctg_locs;
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
    BarrierTimer timer(__FILEFUNC__);
    string dump_fname = "ctg_kmers-" + to_string(kmer_len) + ".txt.gz";
    get_rank_path(dump_fname, rank_me());
    zstr::ofstream dump_file(dump_fname);
    ostringstream out_buf;
    ProgressBar progbar(kmer_map->size(), "Dumping kmers to " + dump_fname);
    int64_t i = 0;
    for (auto &elem : *kmer_map) {
      auto ctg_loc = &elem.second;
      out_buf << elem.first << " " << elem.second.size() << " " << ctg_loc->clen << " " << ctg_loc->pos_in_ctg << " "
              << ctg_loc->is_rc << endl;
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
                             IntermittentTimer &fetch_ctg_seqs_timer, IntermittentTimer &aln_kernel_timer) {
    int rlen = rseq.length();
    string rseq_rc = revcomp(rseq);
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
      // fetch only the substring
      fetch_ctg_seqs_timer.start();
      rget(ctg_loc.seq_gptr + cstart, seq_buf, overlap_len).wait();
      fetch_ctg_seqs_timer.stop();
      string ctg_subseq(seq_buf, overlap_len);
      ctg_seq_bytes_fetched += overlap_len;
      align_read(rname, ctg_loc.cid, read_subseq, ctg_subseq, rstart, rlen, cstart, ctg_loc.clen, orient, overlap_len,
                 aln_kernel_timer);
      num_alns++;
    }
    delete[] seq_buf;
  }

  void sort_alns() { alns->sort_alns(); }

  int get_gpu_mem_avail() { return gpu_mem_avail; }
};

template<int MAX_K>
static void build_alignment_index(KmerCtgDHT<MAX_K> &kmer_ctg_dht, Contigs &ctgs, int min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
  vector<Kmer<MAX_K>> kmers;
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() < min_ctg_len) continue;
    global_ptr<char> seq_gptr = allocate<char>(ctg->seq.length() + 1);
    strcpy(seq_gptr.local(), ctg->seq.c_str());
    CtgLoc ctg_loc = { .cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length() };
    Kmer<MAX_K>::get_kmers(kmer_ctg_dht.kmer_len, ctg->seq, kmers);
    num_kmers += kmers.size();
    for (int i = 0; i < kmers.size(); i++) {
      ctg_loc.pos_in_ctg = i;
      kmer_ctg_dht.add_kmer(kmers[i], ctg_loc);
      progress();
    }
  }
  progbar.done();
  kmer_ctg_dht.flush_add_kmers();
  auto tot_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto num_kmers_in_ht = kmer_ctg_dht.get_num_kmers();
  SLOG_VERBOSE("Processed ", tot_num_kmers, " seeds from contigs, added ", num_kmers_in_ht, "\n");
  auto num_dropped_seed_to_ctgs = kmer_ctg_dht.get_num_dropped_seed_to_ctgs();
  if (num_dropped_seed_to_ctgs)
    SLOG_VERBOSE("Dropped ", num_dropped_seed_to_ctgs, " non-unique seed-to-contig mappings (",
                 setprecision(2), fixed, (100.0 * num_dropped_seed_to_ctgs / tot_num_kmers), "%)\n");
}


struct ReadRecord {
  string id;
  string seq;
  string quals;

  HASH_TABLE<cid_t, ReadAndCtgLoc> aligned_ctgs_map;

  ReadRecord(const string &id, const string &seq, const string &quals) : id(id), seq(seq), quals(quals) {}
};


struct KmerToRead {
  ReadRecord *read_record;
  int pos_in_read;
  bool is_rc;
};


template<int MAX_K>
static int align_kmers(KmerCtgDHT<MAX_K> &kmer_ctg_dht, HASH_TABLE<Kmer<MAX_K>, vector<KmerToRead>> &kmer_read_map,
                       vector<ReadRecord*> &read_records, IntermittentTimer &compute_alns_timer,
                       IntermittentTimer &get_ctgs_timer, IntermittentTimer &fetch_ctg_seqs_timer,
                       IntermittentTimer &aln_kernel_timer, int64_t &num_excess_alns_reads) {
  // extract a list of kmers for each target rank
  auto kmer_lists = new vector<Kmer<MAX_K>>[rank_n()];
  for (auto &elem : kmer_read_map) {
    auto kmer = elem.first;
    kmer_lists[kmer_ctg_dht.get_target_rank(kmer)].push_back(kmer);
  }
  get_ctgs_timer.start();
  future<> fut_chain = make_future();
  // fetch ctgs for each set of kmers from target ranks
  for (int i = 0; i < rank_n(); i++) {
    progress();
    // skip targets that have no ctgs - this should reduce communication at scale
    if (kmer_lists[i].empty()) continue;
    auto fut = kmer_ctg_dht.get_ctgs_with_kmers(i, kmer_lists[i]).then(
      [&](vector<KmerCtgLoc<MAX_K>> kmer_ctg_locs) {
        // iterate through the kmers, each one has an associated ctg location
        for (auto &kmer_ctg_loc : kmer_ctg_locs) {
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
        }
      });
    fut_chain = when_all(fut_chain, fut);
  }
  fut_chain.wait();
  get_ctgs_timer.stop();
  delete[] kmer_lists;
  kmer_read_map.clear();

  compute_alns_timer.start();
  int num_reads_aligned = 0;
  // create a new list of records with all reads having < KLIGN_MAX_ALNS_PER_READ, i.e. those without excessive mappings
  // setting KLIGN_MAX_ALNS_PER_READ to zero means don't drop any
  vector<ReadRecord*> good_read_records;
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
      kmer_ctg_dht.compute_alns_for_read(&read_record->aligned_ctgs_map, read_record->id, read_record->seq,
                                         fetch_ctg_seqs_timer, aln_kernel_timer);
    }
    delete read_record;
  }
  read_records.clear();
  compute_alns_timer.stop();
  return num_reads_aligned;
}

template<int MAX_K>
static double do_alignments(KmerCtgDHT<MAX_K> &kmer_ctg_dht, vector<PackedReads*> &packed_reads_list, int seed_space,
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
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    string read_id, read_seq, quals;
    ProgressBar progbar(packed_reads->get_local_num_reads(), "Aligning reads to contigs");
    vector<ReadRecord*> read_records;
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
      for (int i = 0; i < kmers.size(); i += seed_space) {
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
      if (filled)
        num_reads_aligned += align_kmers(kmer_ctg_dht, kmer_read_map, read_records, compute_alns_timer, get_ctgs_timer,
                                         fetch_ctg_seqs_timer, aln_kernel_timer, num_excess_alns_reads);
      num_reads++;
    }
    if (read_records.size())
      num_reads_aligned += align_kmers(kmer_ctg_dht, kmer_read_map, read_records, compute_alns_timer, get_ctgs_timer,
                                       fetch_ctg_seqs_timer, aln_kernel_timer, num_excess_alns_reads);
    // make sure to do any outstanding kernel block alignments
    kmer_ctg_dht.kernel_align_block(aln_kernel_timer);
    kmer_ctg_dht.sort_alns();
    progbar.done();
    barrier();
  }
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Parsed ", tot_num_reads, " reads, with ", reduce_one(tot_num_kmers, op_fast_add, 0).wait(), " seeds\n");
  auto tot_num_alns = kmer_ctg_dht.get_num_alns();
  SLOG_VERBOSE("Found ", tot_num_alns, " alignments of which ", perc_str(kmer_ctg_dht.get_num_perfect_alns(), tot_num_alns),
               " are perfect\n");
  auto tot_excess_alns_reads = reduce_one(num_excess_alns_reads, op_fast_add, 0).wait();
  if (num_excess_alns_reads)
    SLOG_VERBOSE("Dropped ", tot_excess_alns_reads, " reads because of alignments in excess of ",
                 KLIGN_MAX_ALNS_PER_READ, "\n");
  auto num_overlaps = kmer_ctg_dht.get_num_overlaps();
  if (num_overlaps) SLOG_VERBOSE("Dropped ", perc_str(num_overlaps, tot_num_alns), " alignments becasue of overlaps\n");
  auto tot_num_reads_aligned = reduce_one(num_reads_aligned, op_fast_add, 0).wait();
  SLOG_VERBOSE("Mapped ", perc_str(tot_num_reads_aligned, tot_num_reads), " reads to contigs\n");
  SLOG_VERBOSE("Average mappings per read ", (double)tot_num_alns / tot_num_reads_aligned, "\n");
  SLOG_VERBOSE("Fetched ", get_size_str(kmer_ctg_dht.get_ctg_seq_bytes_fetched()), " of contig sequences\n");

  get_ctgs_timer.done_all();
  fetch_ctg_seqs_timer.done_all();
  compute_alns_timer.done_all();
  double aln_kernel_secs = aln_kernel_timer.get_elapsed();
  aln_kernel_timer.done_all();
  return aln_kernel_secs;
}

template<int MAX_K>
double find_alignments(unsigned kmer_len, vector<PackedReads*> &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                       Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit, bool compute_cigar=false, int min_ctg_len=0, int ranks_per_gpu = 0) {

  BarrierTimer timer(__FILEFUNC__);
  _num_dropped_seed_to_ctgs = 0;
  Kmer<MAX_K>::set_k(kmer_len);
  SLOG_VERBOSE("Aligning with seed size of ", kmer_len, "\n");
  // default for normal alignments in the pipeline, but for final alignments, uses minimap2 defaults
  AlnScoring aln_scoring = { .match = ALN_MATCH_SCORE, .mismatch = ALN_MISMATCH_COST, .gap_opening = ALN_GAP_OPENING_COST,
                             .gap_extending = ALN_GAP_EXTENDING_COST, .ambiguity = ALN_AMBIGUITY_COST };
  if (compute_cigar) {
    AlnScoring alt_aln_scoring = { .match = 2, .mismatch = 4, .gap_opening = 4, .gap_extending = 2, .ambiguity = 1};
    aln_scoring = alt_aln_scoring;
  }
  SLOG_VERBOSE("Alignment scoring parameters: ", aln_scoring.to_string(), "\n");
  KmerCtgDHT<MAX_K> kmer_ctg_dht(kmer_len, max_store_size, max_rpcs_in_flight, alns, aln_scoring, rlen_limit, compute_cigar, ranks_per_gpu);
  barrier();
  build_alignment_index(kmer_ctg_dht, ctgs, min_ctg_len);
#ifdef DEBUG
  //kmer_ctg_dht.dump_ctg_kmers();
#endif
  double kernel_elapsed = do_alignments(kmer_ctg_dht, packed_reads_list, seed_space, compute_cigar);
  barrier();
  auto num_alns = kmer_ctg_dht.get_num_alns();
  auto num_dups = alns.get_num_dups();
  if (num_dups) SLOG_VERBOSE("Number of duplicate alignments ", perc_str(num_dups, num_alns), "\n");
  barrier();
  return kernel_elapsed;
}

#define FA(KMER_LEN) \
    template \
    double find_alignments<KMER_LEN>(unsigned, vector<PackedReads*>&, int, int, Contigs&, Alns&, int, int, bool, int, int)

FA(32);
#if MAX_BUILD_KMER >= 64
FA(64);
#endif
#if MAX_BUILD_KMER >= 96
FA(96);
#endif
#if MAX_BUILD_KMER >= 128
FA(128);
#endif
#if MAX_BUILD_KMER >= 160
FA(160);
#endif

#undef FA

