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

#include <fstream>
#include <iostream>
#include <regex>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "contigs.hpp"
#include "kmer_dht.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/three_tier_aggr_store.hpp"
#include "utils.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

enum class AlnStatus { NO_ALN, OVERLAPS_CONTIG, EXTENDS_CONTIG };

struct CtgInfo {
  int64_t cid;
  char orient;
  char side;
  UPCXX_SERIALIZED_FIELDS(cid, orient, side);
};
struct ReadCtgInfo {
  string read_id;
  CtgInfo ctg_info;
  UPCXX_SERIALIZED_FIELDS(read_id, ctg_info);
};

class ReadsToCtgsDHT {
 public:
  using reads_to_ctgs_map_t = HASH_TABLE<string, vector<CtgInfo>>;
  using dist_reads_to_ctgs_map_t = dist_object<reads_to_ctgs_map_t>;
  static size_t get_target_rank(const string &read_id) { return std::hash<string>{}(read_id) % rank_n(); }

 private:
  dist_reads_to_ctgs_map_t reads_to_ctgs_map;
  ThreeTierAggrStore<ReadCtgInfo> rtc_store;

 public:
  ReadsToCtgsDHT(int64_t initial_size)
      : reads_to_ctgs_map({})
      , rtc_store() {
    reads_to_ctgs_map->reserve(initial_size);
    rtc_store.set_update_func([&reads_to_ctgs_map = this->reads_to_ctgs_map](const ReadCtgInfo &read_ctg_info) {
      const auto it = reads_to_ctgs_map->find(read_ctg_info.read_id);
      if (it == reads_to_ctgs_map->end())
        reads_to_ctgs_map->insert({read_ctg_info.read_id, {read_ctg_info.ctg_info}});
      else
        it->second.push_back(read_ctg_info.ctg_info);
    });
    // TODO rtc_store.set_size();
  }

  void clear() {
    reads_to_ctgs_map_t().swap(*reads_to_ctgs_map);
    rtc_store.clear();
  }

  void add(const string &read_id, int64_t cid, char orient, char side) {
    ReadCtgInfo read_ctg_info = {.read_id = read_id, .ctg_info = {.cid = cid, .orient = orient, .side = side}};
    rtc_store.update(get_target_rank(read_id), read_ctg_info);
  }

  void flush_updates() {
    rtc_store.flush_updates();
    rtc_store.clear();
  }

  int64_t get_num_mappings() { return reduce_one(reads_to_ctgs_map->size(), op_fast_add, 0).wait(); }

  vector<CtgInfo> get_ctgs(string &read_id) {
    return upcxx::rpc(get_target_rank(read_id),
                      [](upcxx::dist_object<reads_to_ctgs_map_t> &reads_to_ctgs_map, string read_id) -> vector<CtgInfo> {
                        const auto it = reads_to_ctgs_map->find(read_id);
                        if (it == reads_to_ctgs_map->end()) return {};
                        return it->second;
                      },
                      reads_to_ctgs_map, read_id)
        .wait();
  }

  future<vector<vector<CtgInfo>>> get_ctgs(intrank_t target_rank, vector<string> &read_ids) {
    return upcxx::rpc(
        target_rank,
        [](upcxx::dist_object<reads_to_ctgs_map_t> &reads_to_ctgs_map, view<string> read_ids) -> vector<vector<CtgInfo>> {
          vector<vector<CtgInfo>> results(read_ids.size());
          size_t i = 0;
          for (const auto &read_id : read_ids) {
            assert(get_target_rank(read_id) == rank_me());
            const auto it = reads_to_ctgs_map->find(read_id);
            if (it != reads_to_ctgs_map->end()) results[i] = it->second;
            i++;
          }
          return results;
        },
        reads_to_ctgs_map, make_view(read_ids.begin(), read_ids.end()));
  }
};

struct ReadSeq {
  string read_id;
  string seq;
  string quals;
  UPCXX_SERIALIZED_FIELDS(read_id, seq, quals);
};

struct CtgWithReads {
  int64_t cid;
  string seq;
  double depth;
  vector<ReadSeq> reads_left;
  vector<ReadSeq> reads_right;
};

struct CtgData {
  int64_t cid;
  string seq;
  double depth;
  UPCXX_SERIALIZED_FIELDS(cid, seq, depth);
};

struct CtgReadData {
  int64_t cid;
  char side;
  ReadSeq read_seq;
  UPCXX_SERIALIZED_FIELDS(cid, side, read_seq);
};

class CtgsWithReadsDHT {
 public:
  using ctgs_map_t = HASH_TABLE<int64_t, CtgWithReads>;
  static size_t get_target_rank(int64_t cid) { return std::hash<int64_t>{}(cid) % rank_n(); }

 private:
  dist_object<ctgs_map_t> ctgs_map;
  ctgs_map_t::iterator ctgs_map_iter;
  ThreeTierAggrStore<CtgData> ctg_store;
  ThreeTierAggrStore<CtgReadData> ctg_read_store;

 public:
  CtgsWithReadsDHT(int64_t num_ctgs)
      : ctgs_map({})
      , ctg_store()
      , ctg_read_store() {
    // pad the local ctg count a bit for this estimate
    ctgs_map->reserve(num_ctgs * 1.2);

    ctg_store.set_update_func([&ctgs_map = this->ctgs_map](const CtgData &ctg_data) {
      const auto it = ctgs_map->find(ctg_data.cid);
      if (it != ctgs_map->end()) DIE("Found duplicate ctg ", ctg_data.cid);
      CtgWithReads ctg_with_reads = {
          .cid = ctg_data.cid, .seq = ctg_data.seq, .depth = ctg_data.depth, .reads_left = {}, .reads_right = {}};
      ctgs_map->insert({ctg_data.cid, ctg_with_reads});
    });
    // TODO ctg_store.set_size();

    ctg_read_store.set_update_func([&ctgs_map = this->ctgs_map](const CtgReadData &ctg_read_data) {
      const auto it = ctgs_map->find(ctg_read_data.cid);
      if (it == ctgs_map->end()) DIE("Could not find ctg ", ctg_read_data.cid);
      const ReadSeq &read_seq = ctg_read_data.read_seq;
      if (ctg_read_data.side == 'L')
        it->second.reads_left.push_back(read_seq);
      else
        it->second.reads_right.push_back(read_seq);
    });
  }

  void add_ctg(Contig &ctg) {
    CtgData ctg_data = {ctg.id, ctg.seq, ctg.depth};
    ctg_store.update(get_target_rank(ctg.id), ctg_data);
    /*
    rpc(get_target_rank(ctg.id),
        [](dist_object<ctgs_map_t> &ctgs_map, int64_t cid, string seq, double depth) {
          const auto it = ctgs_map->find(cid);
          if (it != ctgs_map->end()) DIE("Found duplicate ctg ", cid);
          CtgWithReads ctg_with_reads = {.cid = cid, .seq = seq, .depth = depth, .reads_left = {}, .reads_right = {}};
          ctgs_map->insert({cid, ctg_with_reads});
        },
        ctgs_map, ctg.id, ctg.seq, ctg.depth)
        .wait();
        */
  }

  void add_read(int64_t cid, char side, ReadSeq read_seq) {
    CtgReadData ctg_read_data = {.cid = cid, .side = side, .read_seq = read_seq};
    add_read(ctg_read_data);
  }
  void add_read(CtgReadData &ctg_read_data) {
    ctg_read_store.update(get_target_rank(ctg_read_data.cid), ctg_read_data);
    /*
    rpc(get_target_rank(cid),
        [](dist_object<ctgs_map_t> &ctgs_map, int64_t cid, char side, string read_id, string seq, string quals) {
          const auto it = ctgs_map->find(cid);
          if (it == ctgs_map->end()) DIE("Could not find ctg ", cid);
          if (side == 'L')
            it->second.reads_left.push_back({read_id, seq, quals});
          else
            it->second.reads_right.push_back({read_id, seq, quals});
        },
        ctgs_map, cid, side, read_seq.read_id, read_seq.seq, read_seq.quals)
        .wait();
        */
  }

  void add_reads(vector<CtgReadData> &ctg_read_datas) {
    for (auto &ctg_read_data : ctg_read_datas) {
      add_read(ctg_read_data);
    }
    ctg_read_datas.clear();
  }

  void flush_ctg_updates() {
    ctg_store.flush_updates();
    ctg_store.clear();
    // TODO ctg_read_store.set_size();
  }

  void flush_read_updates() {
    ctg_read_store.flush_updates();
    ctg_read_store.clear();
  }

  int64_t get_num_ctgs() { return reduce_one(ctgs_map->size(), op_fast_add, 0).wait(); }

  int64_t get_local_num_ctgs() { return ctgs_map->size(); }

  CtgWithReads *get_first_local_ctg() {
    ctgs_map_iter = ctgs_map->begin();
    if (ctgs_map_iter == ctgs_map->end()) return nullptr;
    auto ctg = &ctgs_map_iter->second;
    ctgs_map_iter++;
    return ctg;
  }

  CtgWithReads *get_next_local_ctg() {
    if (ctgs_map_iter == ctgs_map->end()) return nullptr;
    auto ctg = &ctgs_map_iter->second;
    ctgs_map_iter++;
    return ctg;
  }
};

struct MerFreqs {
  // how many times this kmer has occurred: don't need to count beyond 65536
  // count of high quality extensions and low quality extensions - structure comes from kmer_dht.hpp
  ExtCounts hi_q_exts, low_q_exts;
  // the final extensions chosen - A,C,G,T, or F,X
  char ext;
  // the count of the final extension
  int count;

  struct MerBase {
    char base;
    uint16_t nvotes_hi_q, nvotes, rating;

    uint16_t get_base_rating(int depth) {
      double min_viable = max(LASSM_MIN_VIABLE_DEPTH * depth, 2.0);
      double min_expected_depth = max(LASSM_MIN_EXPECTED_DEPTH * depth, 2.0);
      if (nvotes == 0) return 0;
      if (nvotes == 1) return 1;
      if (nvotes < min_viable) return 2;
      if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q < min_viable) return 3;
      if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q >= min_viable) return 4;
      if (nvotes >= min_expected_depth && nvotes_hi_q < min_viable) return 5;
      if (nvotes >= min_expected_depth && min_viable < nvotes_hi_q && nvotes_hi_q < min_expected_depth) return 6;
      return 7;
    }
  };

  void set_ext(int seq_depth) {
    // set extension similarly to how it is done with localassm in mhm
    MerBase mer_bases[4] = {{.base = 'A', .nvotes_hi_q = hi_q_exts.count_A, .nvotes = low_q_exts.count_A},
                            {.base = 'C', .nvotes_hi_q = hi_q_exts.count_C, .nvotes = low_q_exts.count_C},
                            {.base = 'G', .nvotes_hi_q = hi_q_exts.count_G, .nvotes = low_q_exts.count_G},
                            {.base = 'T', .nvotes_hi_q = hi_q_exts.count_T, .nvotes = low_q_exts.count_T}};
    for (int i = 0; i < 4; i++) {
      mer_bases[i].rating = mer_bases[i].get_base_rating(seq_depth);
    }
    // sort bases in descending order of quality
    sort(mer_bases, mer_bases + sizeof(mer_bases) / sizeof(mer_bases[0]), [](const auto &elem1, const auto &elem2) -> bool {
      if (elem1.rating != elem2.rating) return elem1.rating > elem2.rating;
      if (elem1.nvotes_hi_q != elem2.nvotes_hi_q) return elem1.nvotes_hi_q > elem2.nvotes_hi_q;
      if (elem1.nvotes != elem2.nvotes) return elem1.nvotes > elem2.nvotes;
      return true;
    });
    int top_rating = mer_bases[0].rating;
    int runner_up_rating = mer_bases[1].rating;
    if (top_rating < runner_up_rating) DIE("top_rating ", top_rating, " < ", runner_up_rating, "\n");
    assert(top_rating >= runner_up_rating);
    int top_rated_base = mer_bases[0].base;
    ext = 'X';
    count = 0;
    // no extension (base = 0) if the runner up is close to the top rating
    // except, if rating is 7 (best quality), then all bases of rating 7 are forks
    if (top_rating > LASSM_RATING_THRES) {  // must have at least minViable bases
      if (top_rating <= 3) {                // must be uncontested
        if (runner_up_rating == 0) ext = top_rated_base;
      } else if (top_rating < 6) {
        if (runner_up_rating < 3) ext = top_rated_base;
      } else if (top_rating == 6) {  // viable and fair hiQ support
        if (runner_up_rating < 4) ext = top_rated_base;
      } else {  // strongest rating trumps
        if (runner_up_rating < 7) {
          ext = top_rated_base;
        } else {
          if (mer_bases[2].rating == 7 || mer_bases[0].nvotes == mer_bases[1].nvotes)
            ext = 'F';
          else if (mer_bases[0].nvotes > mer_bases[1].nvotes)
            ext = mer_bases[0].base;
          else if (mer_bases[1].nvotes > mer_bases[0].nvotes)
            ext = mer_bases[1].base;
        }
      }
    }
    for (int i = 0; i < 4; i++) {
      if (mer_bases[i].base == ext) {
        count = mer_bases[i].nvotes;
        break;
      }
    }
  }
};

using MerMap = HASH_TABLE<string, MerFreqs>;

static void process_reads(unsigned kmer_len, vector<PackedReads *> &packed_reads_list, ReadsToCtgsDHT &reads_to_ctgs,
                          CtgsWithReadsDHT &ctgs_dht) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_reads = 0;
  int64_t num_read_maps_found = 0;
  future<> all_done = make_future();
  vector<CtgReadData> ctgs_to_add;

  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    string id, seq, quals;
    ProgressBar progbar(packed_reads->get_local_num_reads() * 2, "Processing reads - two stage");
    vector<vector<string>> rank_read_ids;
    vector<vector<uint64_t>> rank_read_idx;
    rank_read_ids.resize(rank_n());
    rank_read_idx.resize(rank_n());
    while (true) {
      progress();
      auto read_idx = packed_reads->get_read_index();
      if (!packed_reads->get_next_read(id, seq, quals)) break;
      progbar.update();
      // this happens when we have a placeholder entry because reads merged
      if (kmer_len > seq.length()) continue;
      num_reads++;
      auto target_rank = ReadsToCtgsDHT::get_target_rank(id);
      rank_read_ids[target_rank].push_back(id);
      rank_read_idx[target_rank].push_back(read_idx);
    }

    for (auto target_rank : foreach_rank_by_node()) {
      progress();
      auto &read_ids = rank_read_ids[target_rank];
      auto &read_idxs = rank_read_idx[target_rank];
      assert(read_ids.size() == read_idxs.size());
      if (read_ids.empty()) continue;
      auto read_ctgs_fut = reads_to_ctgs.get_ctgs(target_rank, read_ids);
      auto fut = read_ctgs_fut.then(
          [&read_ids, &read_idxs, &progbar, &packed_reads, &num_read_maps_found, &ctgs_to_add](vector<vector<CtgInfo>> read_ctgs) {
            assert(read_ctgs.size() == read_ids.size());
            string id, seq, quals, seq_rc, quals_rc;
            for (size_t i = 0; i < read_ctgs.size(); i++) {
              progbar.update();
              auto &read_id = read_ids[i];
              auto &read_idx = read_idxs[i];
              auto &ctgs = read_ctgs[i];
              if (ctgs.size()) {
                num_read_maps_found++;
                packed_reads->get_read(read_idx, id, seq, quals);
                assert(id.compare(read_id) == 0);
                bool was_revcomp = false;
                for (auto &ctg : ctgs) {
                  if ((ctg.orient == '-' && ctg.side == 'R') || (ctg.orient == '+' && ctg.side == 'L')) {
                    if (!was_revcomp) {
                      seq_rc = revcomp(seq);
                      quals_rc = quals;
                      reverse(quals_rc.begin(), quals_rc.end());
                      was_revcomp = true;
                    }
                    ctgs_to_add.push_back({ctg.cid, ctg.side, {id, seq_rc, quals_rc}});
                    // ctgs_dht.add_read(ctg.cid, ctg.side, {id, seq_rc, quals_rc});
                  } else {
                    ctgs_to_add.push_back({ctg.cid, ctg.side, {id, seq, quals}});
                    // ctgs_dht.add_read(ctg.cid, ctg.side, {id, seq, quals});
                  }
                }
              }
              // clear out some memory
              vector<string>().swap(read_ids);
              vector<uint64_t>().swap(read_idxs);
            }
          });
      limit_outstanding_futures(fut).wait();
      ctgs_dht.add_reads(ctgs_to_add);
    }

    auto tot_num_reads_fut = reduce_one(num_reads, op_fast_add, 0);
    auto tot_num_read_maps_found_fut = reduce_one(num_read_maps_found, op_fast_add, 0);
    all_done = when_all(all_done, tot_num_reads_fut, tot_num_read_maps_found_fut)
                   .then([](int64_t tot_num_reads, int64_t tot_num_read_maps_found) {
                     SLOG_VERBOSE("Found ", perc_str(tot_num_read_maps_found, tot_num_reads), " reads that map to contigs\n");
                   });
    all_done = when_all(all_done, progbar.set_done());
  }
  auto all_outstanding = flush_outstanding_futures_async();
  while (!all_outstanding.ready()) {
    ctgs_dht.add_reads(ctgs_to_add);
    progress();
  }
  ctgs_dht.add_reads(ctgs_to_add);
  assert(flush_outstanding_futures_async().ready());
  assert(ctgs_to_add.empty());
  ctgs_dht.flush_read_updates();
  all_done.wait();
  // implicit barrier on exit
}

static void get_best_aln_for_read(const Alns &alns, int64_t &i, Aln &best_aln, AlnStatus &best_start_status,
                                  AlnStatus &best_end_status, int64_t &num_alns_found, int64_t &num_alns_invalid) {
  auto classify_aln = [](int runaligned, int cunaligned) -> AlnStatus {
    if (runaligned > cunaligned && cunaligned < KLIGN_UNALIGNED_THRES) return AlnStatus::EXTENDS_CONTIG;
    if (runaligned <= cunaligned && runaligned < KLIGN_UNALIGNED_THRES) return AlnStatus::OVERLAPS_CONTIG;
    return AlnStatus::NO_ALN;
  };

  // choose the highest scoring aln for this read that is useful
  best_start_status = AlnStatus::NO_ALN;
  best_end_status = AlnStatus::NO_ALN;
  string start_read_id = "";
  int best_aln_score = 0;
  best_aln.read_id = "";
  for (; i < (int64_t)alns.size(); i++) {
    const Aln aln = alns.get_aln(i);
    // alns for a new read
    if (start_read_id != "" && aln.read_id != start_read_id) return;
    num_alns_found++;
    if (aln.score1 < best_aln_score) continue;
    AlnStatus start_status, end_status;
    if (aln.orient == '+') {
      start_status = classify_aln(aln.rstart - 1, aln.cstart - 1);
      end_status = classify_aln(aln.rlen - aln.rstop, aln.clen - aln.cstop);
    } else {
      // for '-' strand, aln is between read and revcomp of contig
      start_status = classify_aln(aln.rstart - 1, aln.clen - aln.cstop);
      end_status = classify_aln(aln.rlen - aln.rstop, aln.cstart - 1);
    }
    if (start_status == AlnStatus::NO_ALN || end_status == AlnStatus::NO_ALN) {
      num_alns_invalid++;
      continue;
    }
    best_aln = aln;
    best_aln_score = aln.score1;
    best_start_status = start_status;
    best_end_status = end_status;
    start_read_id = aln.read_id;
  }
}

void process_alns(const Alns &alns, ReadsToCtgsDHT &reads_to_ctgs, int insert_avg, int insert_stddev) {
  auto pair_overlap = [](Aln &aln, int min_pair_len) -> bool {
    // make sure that the mate won't overlap the same contig
    if (aln.orient == '+') {
      if (min_pair_len - aln.rlen - aln.rstart + 1 <= aln.clen - aln.cstart) return true;
    } else {
      if (min_pair_len - 2 * aln.rlen + aln.rstart - 1 <= aln.cstart) return true;
    }
    return false;
  };

  BarrierTimer timer(__FILEFUNC__);
  int64_t num_alns_found = 0, num_alns_invalid = 0, num_direct = 0, num_proj = 0;
  int min_pair_len = insert_avg + 3 * insert_stddev;
  IntermittentTimer t_get_alns(__FILENAME__ + string(":") + "get alns reads to contigs");
  int64_t aln_i = 0;
  AlnStatus start_status, end_status;
  ProgressBar progbar(alns.size(), "Getting read-to-contig mappings from alignments");
  while (aln_i < (int64_t)alns.size()) {
    progress();
    Aln aln;
    t_get_alns.start();
    get_best_aln_for_read(alns, aln_i, aln, start_status, end_status, num_alns_found, num_alns_invalid);
    t_get_alns.stop();
    progbar.update(aln_i);
    if (aln.read_id.empty()) continue;
    // add a direct extension to the contig, start or end
    if (start_status == AlnStatus::EXTENDS_CONTIG) {
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient, aln.orient == '+' ? 'L' : 'R');
      num_direct++;
    } else if (end_status == AlnStatus::EXTENDS_CONTIG) {
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient, aln.orient == '+' ? 'R' : 'L');
      num_direct++;
    }
    // FIXME: if read is longer than one pair read length, don't look for mate since it is a merged read
    // add mate pair if feasible
    if (!pair_overlap(aln, min_pair_len)) {
      // indicate the other pair number
      int len = aln.read_id.length();
      assert(len > 1);
      if (aln.read_id[len - 1] == '1')
        aln.read_id[len - 1] = '2';
      else if (aln.read_id[len - 1] == '2')
        aln.read_id[len - 1] = '1';
      else
        DIE("Bad pair number ", (int)aln.read_id[len - 1], " in read: ", aln.read_id);
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient == '+' ? '-' : '+', aln.orient == '+' ? 'R' : 'L');
      num_proj++;
    }
  }
  reads_to_ctgs.flush_updates();
  progbar.done();
  barrier();
  t_get_alns.done_all();
  auto tot_alns_found = reduce_one(num_alns_found, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed ", tot_alns_found, " alignments:\n");
  SLOG_VERBOSE("  invalid:   ", perc_str(reduce_one(num_alns_invalid, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("  direct:    ", perc_str(reduce_one(num_direct, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("  projected: ", perc_str(reduce_one(num_proj, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("Added ", reads_to_ctgs.get_num_mappings(), " mappings\n");
}

static void add_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  // process the local ctgs and insert into the distributed hash table
  ProgressBar progbar(ctgs.size(), "Adding contigs to distributed hash table");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    progbar.update();
    ctgs_dht.add_ctg(*it);
    progress();
  }
  ctgs_dht.flush_ctg_updates();
  progbar.done();
  barrier();
  SLOG_VERBOSE("Added ", ctgs_dht.get_num_ctgs(), " contigs\n");
}

static void count_mers(vector<ReadSeq> &reads, MerMap &mers_ht, int seq_depth, int mer_len, int qual_offset,
                       int64_t &excess_reads) {
  int num_reads = 0;
  // split reads into kmers and count frequency of high quality extensions
  for (auto &read_seq : reads) {
    num_reads++;
    if (num_reads > LASSM_MAX_COUNT_MERS_READS) {
      excess_reads += reads.size() - LASSM_MAX_COUNT_MERS_READS;
      break;
    }
    progress();
    if (mer_len >= (int)read_seq.seq.length()) continue;
    int num_mers = read_seq.seq.length() - mer_len;
    for (int start = 0; start < num_mers; start++) {
      // skip mers that contain Ns
      if (read_seq.seq.find("N", start) != string::npos) continue;
      string mer = read_seq.seq.substr(start, mer_len);
      auto it = mers_ht.find(mer);
      if (it == mers_ht.end()) {
        mers_ht.insert({mer, {.hi_q_exts = {0}, .low_q_exts = {0}, .ext = 0, .count = 0}});
        it = mers_ht.find(mer);
      }
      int ext_pos = start + mer_len;
      assert(ext_pos < (int)read_seq.seq.length());
      char ext = read_seq.seq[ext_pos];
      if (ext == 'N') continue;
      int qual = read_seq.quals[ext_pos] - qual_offset;
      if (qual >= LASSM_MIN_QUAL) it->second.low_q_exts.inc(ext, 1);
      if (qual >= LASSM_MIN_HI_QUAL) it->second.hi_q_exts.inc(ext, 1);
    }
  }
  // now set extension choices
  for (auto &elem : mers_ht) {
    elem.second.set_ext(seq_depth);
  }
}

// return the result of the walk (f, r or x)
static char walk_mers(MerMap &mers_ht, string &mer, string &walk, int mer_len, int walk_len_limit) {
  HASH_TABLE<string, bool> loop_check_ht;
  char walk_result = 'X';
  for (int nsteps = 0; nsteps < walk_len_limit; nsteps++) {
    if (!(nsteps % 10)) progress();
    // check for a cycle in the graph
    if (loop_check_ht.find(mer) != loop_check_ht.end()) {
      walk_result = 'R';
      break;
    } else {
      loop_check_ht.insert({mer, true});
    }
    auto it = mers_ht.find(mer);
    if (it == mers_ht.end()) {
      walk_result = 'X';
      break;
    }
    char ext = it->second.ext;
    if (ext == 'F' || ext == 'X') {
      walk_result = ext;
      break;
    }
    mer.erase(0, 1);
    mer += ext;
    walk += ext;
  }
  return walk_result;
}

static string iterative_walks(string &seq, int seq_depth, vector<ReadSeq> &reads, int max_mer_len, int kmer_len, int qual_offset,
                              int walk_len_limit, array<int64_t, 3> &term_counts, int64_t &num_walks, int64_t &max_walk_len,
                              int64_t &sum_ext, IntermittentTimer &count_mers_timer, IntermittentTimer &walk_mers_timer,
                              int64_t &excess_reads) {
  int min_mer_len = LASSM_MIN_KMER_LEN;
  max_mer_len = min(max_mer_len, (int)seq.length());
  // iteratively walk starting from kmer_size, increasing mer size on a fork (F) or repeat (R),
  // and decreasing on an end of path (X)
  // look for the longest walk. Note we have to restart from the beginning for each walk to ensure
  // that all loops will be detected
  string longest_walk = "";
  int shift = 0;
  DBG_VERBOSE("  reads:\n");
#ifdef DEBUG
  for (auto &read_seq : reads) {
    DBG_VERBOSE("    ", read_seq.read_id, "\n", read_seq.seq, "\n");
  }
#endif
  for (int mer_len = kmer_len; mer_len >= min_mer_len && mer_len <= max_mer_len; mer_len += shift) {
    count_mers_timer.start();
    MerMap mers_ht;
    count_mers(reads, mers_ht, seq_depth, mer_len, qual_offset, excess_reads);
    count_mers_timer.stop();
    string mer = seq.substr(seq.length() - mer_len);
    string walk = "";
    walk_mers_timer.start();
    char walk_result = walk_mers(mers_ht, mer, walk, mer_len, walk_len_limit);
    walk_mers_timer.stop();
    int walk_len = walk.length();
    if (walk_len > (int)longest_walk.length()) longest_walk = walk;
    if (walk_result == 'X') {
      term_counts[0]++;
      // walk reaches a dead-end, downshift, unless we were upshifting
      if (shift == LASSM_SHIFT_SIZE) break;
      shift = -LASSM_SHIFT_SIZE;
    } else {
      if (walk_result == 'F')
        term_counts[1]++;
      else
        term_counts[2]++;
      // otherwise walk must end with a fork or repeat, so upshift
      if (shift == -LASSM_SHIFT_SIZE) break;
      if (mer_len > (int)seq.length()) break;
      shift = LASSM_SHIFT_SIZE;
    }
  }
  if (!longest_walk.empty()) {
    num_walks++;
    max_walk_len = max(max_walk_len, (int64_t)longest_walk.length());
    sum_ext += longest_walk.length();
  }
  return longest_walk;
}

struct WalkMetrics {
  int64_t num_walks = 0, sum_clen = 0, sum_ext = 0, max_walk_len = 0, num_reads = 0, num_sides = 0, max_num_reads = 0,
          excess_reads = 0;
  array<int64_t, 3> term_counts = {0};
  WalkMetrics &operator+(const WalkMetrics &wm) {
    num_walks += wm.num_walks;
    sum_clen += wm.sum_clen;
    max_walk_len += wm.max_walk_len;
    num_reads += wm.num_reads;
    num_sides += wm.num_sides;
    max_num_reads = std::max(max_num_reads, wm.max_num_reads);
    excess_reads += wm.excess_reads;
    term_counts[0] += wm.term_counts[0];
    term_counts[1] += wm.term_counts[1];
    term_counts[2] += wm.term_counts[2];
    return *this;
  }
};

static void extend_ctg(CtgWithReads *ctg, WalkMetrics &wm, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len,
                       int qual_offset, int walk_len_limit, IntermittentTimer &count_mers_timer,
                       IntermittentTimer &walk_mers_timer) {
  wm.sum_clen += ctg->seq.length();
  if (ctg->reads_right.size()) {
    wm.num_sides++;
    wm.num_reads += ctg->reads_right.size();
    wm.max_num_reads = max(wm.max_num_reads, (int64_t)ctg->reads_right.size());
    DBG_VERBOSE("walk right ctg ", ctg->cid, " ", ctg->depth, "\n", ctg->seq, "\n");
    // have to do right first because the contig needs to be revcomped for the left
    string right_walk =
        iterative_walks(ctg->seq, ctg->depth, ctg->reads_right, max_kmer_len, kmer_len, qual_offset, walk_len_limit, wm.term_counts,
                        wm.num_walks, wm.max_walk_len, wm.sum_ext, count_mers_timer, walk_mers_timer, wm.excess_reads);
    if (!right_walk.empty()) ctg->seq += right_walk;
  }
  if (ctg->reads_left.size()) {
    wm.num_sides++;
    wm.num_reads += ctg->reads_left.size();
    wm.max_num_reads = max(wm.max_num_reads, (int64_t)ctg->reads_left.size());
    string seq_rc = revcomp(ctg->seq);
    DBG_VERBOSE("walk left ctg ", ctg->cid, " ", ctg->depth, "\n", seq_rc, "\n");
    string left_walk =
        iterative_walks(seq_rc, ctg->depth, ctg->reads_left, max_kmer_len, kmer_len, qual_offset, walk_len_limit, wm.term_counts,
                        wm.num_walks, wm.max_walk_len, wm.sum_ext, count_mers_timer, walk_mers_timer, wm.excess_reads);
    if (!left_walk.empty()) {
      left_walk = revcomp(left_walk);
      ctg->seq.insert(0, left_walk);
    }
  }
}

static void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs, int insert_avg, int insert_stddev, int max_kmer_len,
                        int kmer_len, int qual_offset) {
  BarrierTimer timer(__FILEFUNC__);
  // walk should never be more than this. Note we use the maximum insert size from all libraries
  int walk_len_limit = insert_avg + 2 * insert_stddev;
  WalkMetrics wm;

  IntermittentTimer count_mers_timer(__FILENAME__ + string(":") + "count_mers"),
      walk_mers_timer(__FILENAME__ + string(":") + "walk_mers");
  ProgressBar progbar(ctgs_dht.get_local_num_ctgs(), "Extending contigs");
  for (auto ctg = ctgs_dht.get_first_local_ctg(); ctg != nullptr; ctg = ctgs_dht.get_next_local_ctg()) {
    progress();
    progbar.update();
    Contig ext_contig;
    extend_ctg(ctg, wm, insert_avg, insert_stddev, max_kmer_len, kmer_len, qual_offset, walk_len_limit, count_mers_timer,
               walk_mers_timer);
    ctgs.add_contig({.id = ctg->cid, .seq = ctg->seq, .depth = ctg->depth});
  }
  progbar.done();
  count_mers_timer.done_all();
  walk_mers_timer.done_all();
  barrier();
  SLOG_VERBOSE("Walk terminations: ", reduce_one(wm.term_counts[0], op_fast_add, 0).wait(), " X, ",
               reduce_one(wm.term_counts[1], op_fast_add, 0).wait(), " F, ", reduce_one(wm.term_counts[2], op_fast_add, 0).wait(),
               " R\n");
  auto tot_num_reads = reduce_one(wm.num_reads, op_fast_add, 0).wait();
  auto tot_num_walks = reduce_one(wm.num_walks, op_fast_add, 0).wait();
  auto tot_sum_ext = reduce_one(wm.sum_ext, op_fast_add, 0).wait();
  auto tot_sum_clen = reduce_one(wm.sum_clen, op_fast_add, 0).wait();
  auto tot_max_walk_len = reduce_one(wm.max_walk_len, op_fast_max, 0).wait();
  auto tot_excess_reads = reduce_one(wm.excess_reads, op_fast_add, 0).wait();
  auto num_ctgs = ctgs_dht.get_num_ctgs();
  auto tot_num_sides = reduce_one(wm.num_sides, op_fast_add, 0).wait();
  SLOG_VERBOSE("Used a total of ", tot_num_reads, " reads, max per ctg ", wm.max_num_reads, " avg per ctg ",
               (num_ctgs > 0 ? (tot_num_reads / num_ctgs) : 0), ", dropped ", perc_str(tot_excess_reads, tot_num_reads),
               " excess reads\n");
  SLOG_VERBOSE("Could walk ", perc_str(tot_num_sides, num_ctgs * 2), " contig sides\n");
  if (tot_sum_clen)
    SLOG_VERBOSE("Found ", tot_num_walks, " walks, total extension length ", tot_sum_ext, " extended ",
                 (double)(tot_sum_ext + tot_sum_clen) / tot_sum_clen, "\n");
  if (tot_num_walks)
    SLOG_VERBOSE("Average walk length ", tot_sum_ext / tot_num_walks, ", max walk length ", tot_max_walk_len, "\n");
}

void localassm(int max_kmer_len, int kmer_len, vector<PackedReads *> &packed_reads_list, int insert_avg, int insert_stddev,
               int qual_offset, Contigs &ctgs, const Alns &alns) {
  BarrierTimer timer(__FILEFUNC__);
  CtgsWithReadsDHT ctgs_dht(ctgs.size());
  add_ctgs(ctgs_dht, ctgs);
  ReadsToCtgsDHT reads_to_ctgs(100);
  // extract read id to ctg id mappings from alignments
  process_alns(alns, reads_to_ctgs, insert_avg, insert_stddev);
  // extract read seqs and add to ctgs
  process_reads(max_kmer_len, packed_reads_list, reads_to_ctgs, ctgs_dht);
  // free the reads_to_contigs map
  reads_to_ctgs.clear();
  // clear out the local contigs
  ctgs.clear();
  ctgs.set_capacity(ctgs_dht.get_local_num_ctgs());
  // extend contigs using locally mapped reads
  extend_ctgs(ctgs_dht, ctgs, insert_avg, insert_stddev, max_kmer_len, kmer_len, qual_offset);
}
