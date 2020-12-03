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

void shuffle_reads(int qual_offset, vector<PackedReads *> &packed_reads_list, Alns &alns, size_t num_ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  using read_to_target_map_t = HASH_TABLE<int64_t, pair<int, int64_t>>;
  dist_object<read_to_target_map_t> read_to_target_map({});
  auto num_reads = 0;
  for (auto packed_reads : packed_reads_list) num_reads += packed_reads->get_local_num_reads();
  read_to_target_map->reserve(num_reads);
  using ctg_to_reads_map_t = HASH_TABLE<int64_t, vector<int64_t>>;
  dist_object<ctg_to_reads_map_t> ctg_to_reads_map({});
  auto tot_num_ctgs = reduce_all(num_ctgs, op_fast_add).wait();
  ctg_to_reads_map->reserve((int64_t)ceil((double)tot_num_ctgs / rank_n()));
  // FIXME: this should use an TTAS

  for (auto &aln : alns) {
    progress();
    // using abs ensures that both reads in a pair are mapped to the same location
    int64_t packed_read_id = abs(PackedRead::to_packed_id(aln.read_id));
    // DBG("shuffle insert read id ", aln.read_id, " packed ", packed_read_id, "\n");
    rpc(
        get_target_rank(packed_read_id),
        [](dist_object<read_to_target_map_t> &read_to_target_map, int64_t read_id, int score, int64_t cid) {
          const auto it = read_to_target_map->find(read_id);
          if (it == read_to_target_map->end()) {
            read_to_target_map->insert({read_id, {score, cid}});
          } else {
            if (it->second.first < score) {
              it->second.first = score;
              it->second.second = cid;
            }
          }
        },
        read_to_target_map, packed_read_id, aln.score1, aln.cid)
        .wait();
    rpc(
        get_target_rank(aln.cid),
        [](dist_object<ctg_to_reads_map_t> &ctg_to_reads_map, int64_t cid, int64_t read_id) {
          const auto it = ctg_to_reads_map->find(read_id);
          if (it == ctg_to_reads_map->end())
            ctg_to_reads_map->insert({cid, {read_id}});
          else
            it->second.push_back(read_id);
        },
        ctg_to_reads_map, aln.cid, packed_read_id)
        .wait();
  }
  barrier();
  auto tot_aln_reads = reduce_one(read_to_target_map->size(), op_fast_add, 0).wait();
  auto tot_read_pairs = reduce_one(num_reads, op_fast_add, 0).wait() / 2;
  SLOG_VERBOSE("Number of read pairs with alignments ", perc_str(tot_aln_reads, tot_read_pairs), "\n");
  int64_t num_target_reads = 0;
  for (auto &[cid, reads] : *ctg_to_reads_map) {
    num_target_reads += reads.size();
  }
  auto tot_num_target_reads = reduce_one(num_target_reads, op_fast_add, 0).wait();
  auto avg_num_target_reads = tot_num_target_reads / rank_n();
  auto max_num_target_reads = reduce_one(num_target_reads, op_fast_max, 0).wait();
  SLOG("Avg number of reads per rank ", avg_num_target_reads, " balance ", (double)avg_num_target_reads / max_num_target_reads,
       "\n");
  barrier();
  int64_t num_not_found = 0;
  dist_object<vector<PackedRead>> new_packed_reads({});
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    for (int i = 0; i < packed_reads->get_local_num_reads(); i += 2) {
      auto &packed_read1 = (*packed_reads)[i];
      auto &packed_read2 = (*packed_reads)[i + 1];
      auto read_id = abs(packed_read1.get_id());
      // DBG("shuffle find read id ", packed_read1.get_str_id(), " packed ", packed_read1.get_id(), "\n");
      int64_t cid = rpc(
                        get_target_rank(read_id),
                        [](dist_object<read_to_target_map_t> &read_to_target_map, int64_t read_id) -> int64_t {
                          const auto it = read_to_target_map->find(read_id);
                          if (it == read_to_target_map->end()) return -1;
                          return it->second.second;
                        },
                        read_to_target_map, read_id)
                        .wait();
      if (cid == -1) num_not_found++;
      int target_ctg_rank = (cid == -1 ? std::experimental::randint(0, rank_n() - 1) : get_target_rank(cid));
      if (target_ctg_rank < 0 || target_ctg_rank >= rank_n())
        DIE("target ctg rank out of range [0, ", rank_n(), "] ", target_ctg_rank);
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
  SLOG("Didn't find contig targets for ", perc_str(tot_num_not_found, tot_read_pairs), " pairs\n");
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
  if (tot_num_new_reads != tot_read_pairs * 2)
    SWARN("Not all reads shuffled, expected ", tot_read_pairs * 2, " but only shuffled ", tot_num_new_reads);
  barrier();
}
