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
#include "upcxx_utils/three_tier_aggr_store.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "utils.hpp"
#include "hash_funcs.h"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

// each read can map to multiple contigs, so first pass is to deterime exactly one cid per read (the highest score)
// generate a map of cids to read vectors from alignments
// determine location of every read using atomics and generate a read_id to location map
// process reads and put in correct locations

intrank_t get_target_rank(int64_t val) {
  return MurmurHash3_x64_64(reinterpret_cast<const void *>(&val), sizeof(int64_t)) % rank_n();
}

using cid_to_reads_map_t = HASH_TABLE<int64_t, vector<int64_t>>;
using read_to_target_map_t = HASH_TABLE<int64_t, int>;

dist_object<cid_to_reads_map_t> process_alns(vector<PackedReads *> &packed_reads_list, Alns &alns, int64_t num_ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  using read_to_cid_map_t = HASH_TABLE<int64_t, pair<int64_t, int>>;
  dist_object<read_to_cid_map_t> read_to_cid_map({});
  ThreeTierAggrStore<tuple<int64_t, int64_t, int>> read_cid_store;
  read_cid_store.set_update_func([&read_to_cid_map](tuple<int64_t, int64_t, int> &&read_cid_info) {
    auto &[read_id, cid, score] = read_cid_info;
    auto it = read_to_cid_map->find(read_id);
    if (it == read_to_cid_map->end()) {
      read_to_cid_map->insert({read_id, {cid, score}});
    } else if (it->second.second < score) {
      it->second.first = cid;
      it->second.second = score;
    }
  });
  int est_update_size = sizeof(tuple<int64_t, int64_t, int>);
  int64_t mem_to_use = 0.1 * get_free_mem() / local_team().rank_n();
  auto max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  read_cid_store.set_size("Read cid store", max_store_bytes);

  for (auto &aln : alns) {
    progress();
    // use abs to ensure both reads in a pair map to the same ctg
    int64_t packed_read_id = abs(PackedRead::to_packed_id(aln.read_id));
    read_cid_store.update(get_target_rank(packed_read_id), {packed_read_id, aln.cid, aln.score1});
  }
  read_cid_store.flush_updates();
  barrier();

  dist_object<cid_to_reads_map_t> cid_to_reads_map({});
  cid_to_reads_map->reserve(num_ctgs);
  ThreeTierAggrStore<pair<int64_t, int64_t>> cid_reads_store;
  cid_reads_store.set_update_func([&cid_to_reads_map](pair<int64_t, int64_t> &&cid_reads_info) {
    auto &[cid, read_id] = cid_reads_info;
    auto it = cid_to_reads_map->find(cid);
    if (it == cid_to_reads_map->end())
      cid_to_reads_map->insert({cid, {read_id}});
    else
      it->second.push_back(read_id);
  });
  est_update_size = sizeof(pair<int64_t, int64_t>);
  mem_to_use = 0.1 * get_free_mem() / local_team().rank_n();
  max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  cid_reads_store.set_size("Read cid store", max_store_bytes);
  for (auto &[read_id, cid_elem] : *read_to_cid_map) {
    progress();
    cid_reads_store.update(get_target_rank(cid_elem.first), {cid_elem.first, read_id});
  }
  cid_reads_store.flush_updates();
  barrier();
  return cid_to_reads_map;
}

dist_object<read_to_target_map_t> compute_read_locations(dist_object<cid_to_reads_map_t> &cid_to_reads_map) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_mapped_reads = 0;
  for (auto &[cid, read_ids] : *cid_to_reads_map) num_mapped_reads += read_ids.size();
  // counted read pairs
  num_mapped_reads *= 2;
  barrier();
  auto all_num_mapped_reads = reduce_all(num_mapped_reads, op_fast_add).wait();
  auto avg_num_mapped_reads = all_num_mapped_reads / rank_n();
  auto max_num_mapped_reads = reduce_one(num_mapped_reads, op_fast_max, 0).wait();
  SLOG("Avg mapped reads per rank ", avg_num_mapped_reads, " max ", max_num_mapped_reads, " balance ",
       (double)avg_num_mapped_reads / max_num_mapped_reads, "\n");
  atomic_domain<int64_t> fetch_add_domain({atomic_op::fetch_add});
  dist_object<global_ptr<int64_t>> read_counter_dobj = (!rank_me() ? new_<int64_t>(0) : nullptr);
  global_ptr<int64_t> read_counter = read_counter_dobj.fetch(0).wait();
  barrier();
  auto read_slot = fetch_add_domain.fetch_add(read_counter, num_mapped_reads, memory_order_relaxed).wait();
  dist_object<read_to_target_map_t> read_to_target_map({});
  read_to_target_map->reserve(avg_num_mapped_reads);
  int block = ceil((double)all_num_mapped_reads / rank_n());
  for (auto &[cid, read_ids] : *cid_to_reads_map) {
    progress();
    for (auto read_id : read_ids) {
      rpc(get_target_rank(read_id),
          [](dist_object<read_to_target_map_t> &read_to_target_map, int64_t read_id, int target) {
            read_to_target_map->insert({read_id, target});
          },
          read_to_target_map, read_id, read_slot / block)
          .wait();
      // each entry is a pair
      read_slot += 2;
    }
  }
  barrier();
  fetch_add_domain.destroy();
  return read_to_target_map;
}

dist_object<vector<PackedRead>> move_reads_to_targets(vector<PackedReads *> &packed_reads_list,
                                                      dist_object<read_to_target_map_t> &read_to_target_map,
                                                      int64_t all_num_reads) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_not_found = 0;
  dist_object<vector<PackedRead>> new_packed_reads({});
  ThreeTierAggrStore<pair<PackedRead, PackedRead>> read_seq_store;
  read_seq_store.set_update_func([&new_packed_reads](pair<PackedRead, PackedRead> &&read_pair_info) {
    new_packed_reads->push_back(read_pair_info.first);
    new_packed_reads->push_back(read_pair_info.second);
  });
  int est_update_size = 600;
  int64_t mem_to_use = 0.1 * get_free_mem() / local_team().rank_n();
  auto max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  read_seq_store.set_size("Read seq store", max_store_bytes);
  future<> fut_chain = make_future();
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    for (int i = 0; i < packed_reads->get_local_num_reads(); i += 2) {
      progress();
      auto &packed_read1 = (*packed_reads)[i];
      auto &packed_read2 = (*packed_reads)[i + 1];
      auto read_id = abs(packed_read1.get_id());
      auto fut = rpc(get_target_rank(read_id),
                     [](dist_object<read_to_target_map_t> &read_to_target_map, int64_t read_id) -> int {
                       const auto it = read_to_target_map->find(read_id);
                       if (it == read_to_target_map->end()) return -1;
                       return it->second;
                     },
                     read_to_target_map, read_id)
                     .then([&num_not_found, &read_seq_store, packed_read1, packed_read2](int target) {
                       if (target == -1) {
                         num_not_found++;
                         target = std::experimental::randint(0, rank_n() - 1);
                       }
                       if (target < 0 || target >= rank_n()) DIE("target out of range ", target);
                       assert(target >= 0 && target < rank_n());
                       read_seq_store.update(target, {packed_read1, packed_read2});
                     });
      fut_chain = when_all(fut_chain, fut);
    }
  }
  fut_chain.wait();
  read_seq_store.flush_updates();
  barrier();
  auto all_num_not_found = reduce_one(num_not_found, op_fast_add, 0).wait();
  SLOG("Didn't find contig targets for ", perc_str(all_num_not_found, all_num_reads / 2), " pairs\n");
  return new_packed_reads;
}

void shuffle_reads(int qual_offset, vector<PackedReads *> &packed_reads_list, Alns &alns, size_t num_ctgs) {
  BarrierTimer timer(__FILEFUNC__);

  int64_t num_reads = 0;
  for (auto packed_reads : packed_reads_list) num_reads += packed_reads->get_local_num_reads();
  auto all_num_reads = reduce_all(num_reads, op_fast_add).wait();

  auto cid_to_reads_map = process_alns(packed_reads_list, alns, num_ctgs);
  auto read_to_target_map = compute_read_locations(cid_to_reads_map);
  auto new_packed_reads = move_reads_to_targets(packed_reads_list, read_to_target_map, all_num_reads);

  // now copy the new packed reads to the old
  for (auto packed_reads : packed_reads_list) delete packed_reads;
  packed_reads_list.clear();
  packed_reads_list.push_back(new PackedReads(qual_offset, *new_packed_reads));
  auto num_reads_received = new_packed_reads->size();
  double avg_num_received = (double)reduce_one(num_reads_received, op_fast_add, 0).wait() / rank_n();
  auto max_reads_received = reduce_one(num_reads_received, op_fast_max, 0).wait();
  SLOG_VERBOSE("Balance in reads ", fixed, setprecision(3), avg_num_received / max_reads_received, "\n");
  auto all_num_new_reads = reduce_one(new_packed_reads->size(), op_fast_add, 0).wait();
  if (all_num_new_reads != all_num_reads)
    SWARN("Not all reads shuffled, expected ", all_num_reads, " but only shuffled ", all_num_new_reads);
  barrier();
}
