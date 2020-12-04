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
#include <experimental/random>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "utils.hpp"
#include "hash_funcs.h"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

intrank_t get_target_rank(int64_t val) {
  return MurmurHash3_x64_64(reinterpret_cast<const void *>(&val), sizeof(int64_t)) % rank_n();
}

using read_to_cid_map_t = HASH_TABLE<int64_t, pair<int, int64_t>>;
using ctg_to_targets_map_t = HASH_TABLE<int64_t, pair<int, int>>;

dist_object<read_to_cid_map_t> generate_read_to_cid_map(vector<PackedReads *> &packed_reads_list, Alns &alns, int64_t tot_reads) {
  dist_object<read_to_cid_map_t> read_to_cid_map({});
  read_to_cid_map->reserve(tot_reads / rank_n() + 1);
  // FIXME: this should use an TTAS
  for (auto &aln : alns) {
    progress();
    // using abs ensures that both reads in a pair are mapped to the same location
    int64_t packed_read_id = abs(PackedRead::to_packed_id(aln.read_id));
    rpc(
        get_target_rank(packed_read_id),
        [](dist_object<read_to_cid_map_t> &read_to_cid_map, int64_t read_id, int score, int64_t cid) {
          auto it = read_to_cid_map->find(read_id);
          if (it == read_to_cid_map->end()) {
            read_to_cid_map->insert({read_id, {score, cid}});
          } else {
            if (it->second.first < score) {
              it->second.first = score;
              it->second.second = cid;
            }
          }
        },
        read_to_cid_map, packed_read_id, aln.score1, aln.cid)
        .wait();
  }
  barrier();
  auto tot_aln_reads = reduce_one(read_to_cid_map->size(), op_fast_add, 0).wait();
  SLOG("Number of read pairs with alignments ", perc_str(tot_aln_reads, tot_reads / 2), "\n");
  return read_to_cid_map;
}

dist_object<ctg_to_targets_map_t> generate_ctg_to_targets_map(dist_object<read_to_cid_map_t> &read_to_cid_map, int64_t num_ctgs,
                                                              int64_t tot_reads) {
  barrier();
  auto tot_num_ctgs = reduce_all(num_ctgs, op_fast_add).wait();
  dist_object<ctg_to_targets_map_t> ctg_to_targets_map({});
  ctg_to_targets_map->reserve((int64_t)ceil((double)tot_num_ctgs / rank_n()));
  for (auto &[read_id, score_cid_pair] : *read_to_cid_map) {
    rpc(
        get_target_rank(score_cid_pair.second),
        [](dist_object<ctg_to_targets_map_t> &ctg_to_targets_map, int64_t cid) {
          auto it = ctg_to_targets_map->find(cid);
          if (it == ctg_to_targets_map->end())
            ctg_to_targets_map->insert({cid, {1, rank_me()}});
          else
            it->second.first++;
        },
        ctg_to_targets_map, score_cid_pair.second)
        .wait();
  }
  barrier();
  atomic_domain<int64_t> fetch_add_domain({atomic_op::fetch_add});
  dist_object<global_ptr<int64_t>> ctg_counter_dobj = (!rank_me() ? new_<int64_t>(0) : nullptr);
  global_ptr<int64_t> ctg_counter = ctg_counter_dobj.fetch(0).wait();
  barrier();
  int64_t num_rank_reads = 0;
  for (auto &[cid, target_pair] : *ctg_to_targets_map) {
    num_rank_reads += target_pair.first * 2;
  }
  auto offset = fetch_add_domain.fetch_add(ctg_counter, num_rank_reads, memory_order_relaxed).wait();
  barrier();
  auto tot_num_rank_reads = reduce_one(num_rank_reads, op_fast_add, 0).wait();
  auto avg_num_rank_reads = tot_num_rank_reads / rank_n();
  auto max_num_rank_reads = reduce_one(num_rank_reads, op_fast_max, 0).wait();
  SLOG("Avg number of reads per rank ", avg_num_rank_reads, " max ", max_num_rank_reads, " balance ",
       (double)avg_num_rank_reads / max_num_rank_reads, "\n");
  int64_t block_size = tot_reads / rank_n();
  for (auto &[cid, target_pair] : *ctg_to_targets_map) {
    auto new_target_rank = offset / block_size;
    DBG("got block ", offset, " size ", target_pair.first * 2, " target is ", new_target_rank, "\n");
    target_pair.second = new_target_rank;
    offset += target_pair.first * 2;
  }
  barrier();
  fetch_add_domain.destroy();
  return ctg_to_targets_map;
}

dist_object<vector<PackedRead>> shuffle_reads_to_targets(vector<PackedReads *> &packed_reads_list,
                                                         dist_object<read_to_cid_map_t> &read_to_cid_map,
                                                         dist_object<ctg_to_targets_map_t> &ctg_to_targets_map, int64_t tot_reads) {
  barrier();
  int64_t num_not_found = 0;
  dist_object<vector<PackedRead>> new_packed_reads({});
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    for (int i = 0; i < packed_reads->get_local_num_reads(); i += 2) {
      auto &packed_read1 = (*packed_reads)[i];
      auto &packed_read2 = (*packed_reads)[i + 1];
      auto read_id = abs(packed_read1.get_id());
      int64_t cid = rpc(
                        get_target_rank(read_id),
                        [](dist_object<read_to_cid_map_t> &read_to_cid_map, int64_t read_id) -> int64_t {
                          const auto it = read_to_cid_map->find(read_id);
                          if (it == read_to_cid_map->end()) return -1;
                          return it->second.second;
                        },
                        read_to_cid_map, read_id)
                        .wait();
      int target_ctg_rank;
      if (cid != -1) {
        target_ctg_rank = rpc(
                              get_target_rank(cid),
                              [](dist_object<ctg_to_targets_map_t> &ctg_to_targets_map, int64_t cid) -> int {
                                const auto it = ctg_to_targets_map->find(cid);
                                if (it == ctg_to_targets_map->end()) return -1;
                                return it->second.second;
                              },
                              ctg_to_targets_map, cid)
                              .wait();
      }
      if (cid == -1) {
        num_not_found++;
        target_ctg_rank = std::experimental::randint(0, rank_n() - 1);
      }
      assert(target_ctg_rank >= 0 && target_ctg_rank < rank_n());
      rpc(
          target_ctg_rank,
          [](dist_object<vector<PackedRead>> &new_packed_reads, PackedRead packed_read1, PackedRead packed_read2) {
            new_packed_reads->push_back(packed_read1);
            new_packed_reads->push_back(packed_read2);
          },
          new_packed_reads, packed_read1, packed_read2)
          .wait();
    }
  }
  barrier();
  auto tot_num_not_found = reduce_one(num_not_found, op_fast_add, 0).wait();
  SLOG("Didn't find contig targets for ", perc_str(tot_num_not_found, tot_reads / 2), " pairs\n");
  return new_packed_reads;
}

void shuffle_reads(int qual_offset, vector<PackedReads *> &packed_reads_list, Alns &alns, size_t num_ctgs) {
  BarrierTimer timer(__FILEFUNC__);

  auto num_reads = 0;
  for (auto packed_reads : packed_reads_list) num_reads += packed_reads->get_local_num_reads();
  auto tot_reads = reduce_all(num_reads, op_fast_add).wait();

  auto read_to_cid_map = generate_read_to_cid_map(packed_reads_list, alns, tot_reads);
  auto ctg_to_targets_map = generate_ctg_to_targets_map(read_to_cid_map, num_ctgs, tot_reads);
  auto new_packed_reads = shuffle_reads_to_targets(packed_reads_list, read_to_cid_map, ctg_to_targets_map, tot_reads);

  // now copy the new packed reads to the old
  for (auto packed_reads : packed_reads_list) {
    delete packed_reads;
  }
  packed_reads_list.clear();
  packed_reads_list.push_back(new PackedReads(qual_offset, *new_packed_reads));
  // FIXME: now the list does not correspond to the original names or files. So any post assembly
  // analysis needs to reload the original reads
  auto num_reads_received = new_packed_reads->size();
  double avg_num_received = (double)reduce_one(num_reads_received, op_fast_add, 0).wait() / rank_n();
  auto max_num_received = reduce_one(num_reads_received, op_fast_max, 0).wait();
  SLOG_VERBOSE("Balance in reads ", fixed, setprecision(3), avg_num_received / max_num_received, "\n");
  auto tot_num_new_reads = reduce_one(new_packed_reads->size(), op_fast_add, 0).wait();
  if (tot_num_new_reads != tot_reads)
    SWARN("Not all reads shuffled, expected ", tot_reads, " but only shuffled ", tot_num_new_reads);
  barrier();
}
