#ifndef _LOCALASSM_HPP
#define _LOCALASSM_HPP

#include <iostream>
#include <fstream>
#include <regex>
#include <upcxx/upcxx.hpp>

using std::unordered_map;


struct ReadSeq {
  string read_id;
  string seq;
  string quals;
};


struct CtgWithReads {
  int64_t cid;
  string seq;
  double depth;
  vector<ReadSeq> reads_start;
  vector<ReadSeq> reads_end;
};


class CtgsWithReadsDHT {
  
  using ctgs_map_t = unordered_map<int64_t, CtgWithReads>;
  dist_object<ctgs_map_t> ctgs_map;
  ctgs_map_t::iterator ctgs_map_iter;

  size_t get_target_rank(int64_t cid) {
    return std::hash<int64_t>{}(cid) % rank_n();
  }
  
public:

  CtgsWithReadsDHT(int64_t num_ctgs)
    : ctgs_map({}) {
    // pad the local ctg count a bit for this estimate
    ctgs_map->reserve(num_ctgs * 1.2);
  }

  void add_ctg(Contig &ctg) {
    rpc(get_target_rank(ctg.id),
        [](dist_object<ctgs_map_t> &ctgs_map, int64_t cid, string seq, double depth) {
          const auto it = ctgs_map->find(cid);
          if (it != ctgs_map->end()) DIE("Found duplicate ctg ", cid);
          CtgWithReads ctg_with_reads = { .cid = cid, .seq = seq, .depth = depth, .reads_start = {}, .reads_end = {} };
          ctgs_map->insert({cid, ctg_with_reads });
        }, ctgs_map, ctg.id, ctg.seq, ctg.depth).wait();
  }

  void add_read(int64_t cid, char side, ReadSeq read_seq) {
    rpc(get_target_rank(cid),
        [](dist_object<ctgs_map_t> &ctgs_map, int64_t cid, char side, string read_id, string seq, string quals) {
          const auto it = ctgs_map->find(cid);
          if (it == ctgs_map->end()) DIE("Could not find ctg ", cid);
          if (side == 'S') it->second.reads_start.push_back({read_id, seq, quals});
          else it->second.reads_end.push_back({read_id, seq, quals});
        }, ctgs_map, cid, side, read_seq.read_id, read_seq.seq, read_seq.quals);
  }
  
  int64_t get_num_ctgs() {
    return reduce_one(ctgs_map->size(), op_fast_add, 0).wait();
  }

  int64_t get_local_num_ctgs() {
    return ctgs_map->size();
  }

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


#endif
