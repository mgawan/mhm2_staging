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

#include "alignments.hpp"

#include <fcntl.h>
#include <limits>
#include <sstream>
#include <string>
#include <upcxx/upcxx.hpp>
#include <unordered_set>

#include "contigs.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/thread_pool.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "version.h"
#include "zstr.hpp"

using namespace upcxx_utils;
using std::ostringstream;
using std::string;
using std::to_string;

using upcxx::future;

//
// class Aln
//

Aln::Aln()
    : read_id("")
    , cid(-1)
    , cstart(0)
    , cstop(0)
    , clen(0)
    , rstart(0)
    , rstop(0)
    , rlen(0)
    , score1(0)
    , score2(0)
    , mismatches(0)
    , sam_string({})
    , read_group_id(-1)
    , orient()
    , identity(0) {}

Aln::Aln(const string &read_id, int64_t cid, int rstart, int rstop, int rlen, int cstart, int cstop, int clen, char orient,
         int score1, int score2, int identity, int mismatches, int read_group_id)
    : read_id(read_id)
    , cid(cid)
    , cstart(cstart)
    , cstop(cstop)
    , clen(clen)
    , rstart(rstart)
    , rstop(rstop)
    , rlen(rlen)
    , score1(score1)
    , score2(score2)
    , mismatches(mismatches)
    , sam_string({})
    , read_group_id(read_group_id)
    , orient(orient)
    , identity(identity) {
  // DBG_VERBOSE(read_id, " cid=", cid, " RG=", read_group_id, " mismatches=", mismatches, "\n");
}

// writes out in the format meraligner uses
string Aln::to_string() const {
  ostringstream os;
  os << read_id << "\t" << rstart + 1 << "\t" << rstop << "\t" << rlen << "\t"
     << "Contig" << cid << "\t" << cstart + 1 << "\t" << cstop << "\t" << clen << "\t" << (orient == '+' ? "Plus" : "Minus") << "\t"
     << score1 << "\t" << score2;
  return os.str();
}

bool Aln::is_valid() const {
  assert(rstart >= 0 && "start >= 0");
  assert(rstop <= rlen && "stop <= len");
  assert(cstart >= 0 && "cstart >= 0");
  assert(cstop <= clen && "cstop <= clen");

  return read_group_id >= 0 && (orient == '+' || orient == '-') && mismatches >= 0 && mismatches <= rlen && identity >= 0 &&
         identity <= 100 && cid >= 0 && read_id.size() > 0;
}

//
// class Alns
//

Alns::Alns()
    : num_dups(0) {}

void Alns::clear() {
  alns.clear();
  vector<Aln>().swap(alns);
}

void Alns::add_aln(Aln &aln) {
#ifdef DEBUG
  // check for duplicate alns to this read - do this backwards because only the most recent entries could be for this read
  for (auto it = alns.rbegin(); it != alns.rend(); ++it) {
    // we have no more entries for this read
    if (it->read_id != aln.read_id) break;
    // now check for equality
    if (it->rstart == aln.rstart && it->rstop == aln.rstop && it->cstart == aln.cstart && it->cstop == aln.cstop) {
      num_dups++;
      return;
    }
  }
#endif
  assert(aln.is_valid());
  alns.push_back(aln);
}

void Alns::append(Alns &more_alns) {
  alns.insert(alns.end(), more_alns.alns.begin(), more_alns.alns.end());
  num_dups += more_alns.num_dups;
  more_alns.clear();
}

Aln &Alns::get_aln(int64_t i) { return alns[i]; }

size_t Alns::size() { return alns.size(); }

void Alns::reserve(size_t capacity) { alns.reserve(capacity); }

void Alns::reset() { alns.clear(); }

int64_t Alns::get_num_dups() { return upcxx::reduce_one(num_dups, upcxx::op_fast_add, 0).wait(); }

void Alns::dump_rank_file(string fname) {
  get_rank_path(fname, rank_me());
  zstr::ofstream f(fname);
  dump_all(f, false);
  f.close();
  upcxx::barrier();
}

void Alns::dump_single_file(const string fname) {
  dist_ofstream of(fname);
  dump_all(of, false);
  of.close();
  upcxx::barrier();
}

void Alns::dump_sam_file(const string fname, const vector<string> &read_group_names, const Contigs &ctgs, int min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__);

  string out_str = "";

  dist_ofstream of(fname);
  future<> all_done = make_future();

  // First all ranks dump Sequence tags - @SQ	SN:Contig0	LN:887
  std::unordered_set<int64_t> passed_contigs;
  for (const auto &ctg : ctgs) {
    if (ctg.seq.length() < min_ctg_len) continue;
    assert(ctg.id >= 0);
    passed_contigs.insert({ctg.id});
    of << "@SQ\tSN:Contig" << to_string(ctg.id) << "\tLN:" << to_string(ctg.seq.length()) << "\n";
  }
  // all @SQ headers aggregated to the top of the file
  all_done = of.flush_collective();

  // rank 0 continues with header
  if (!upcxx::rank_me()) {
    // add ReadGroup tags - @RG ID:[0-n] DS:filename
    for (int i = 0; i < read_group_names.size(); i++) {
      string basefilename = upcxx_utils::get_basename(read_group_names[i]);
      of << "@RG\tID:" << to_string(i) << "\tDS:" << basefilename << "\n";
    }
    // add program information
    of << "@PG\tID:MHM2\tPN:MHM2\tVN:" << string(MHM2_VERSION) << "\n";
  }

  // invalidate (temporarily) alignments to short contigs
  for (auto &aln : alns) {
    if (passed_contigs.find(aln.cid) == passed_contigs.end()) {
      assert(aln.cid >= 0);
      aln.cid = std::numeric_limits<int64_t>::min() + aln.cid;  // set the cid negative to invalidate
      assert(aln.cid < 0);
    }
  }
  std::unordered_set<int64_t>().swap(passed_contigs);

  // next alignments.  rank0 will be first with the remaining header fields
  dump_all(of, true);

  // re-validate alignments
  for (auto &aln : alns) {
    if (aln.cid < 0) {
      aln.cid = aln.cid - std::numeric_limits<int64_t>::min();
    }
    assert(aln.cid >= 0);
  }
  
  all_done = when_all(all_done, of.close_async());
  all_done.wait();
  of.report_timings().wait();
}

int Alns::calculate_unmerged_rlen() {
  BarrierTimer timer(__FILEFUNC__);
  // get the unmerged read length - most common read length
  HASH_TABLE<int, int64_t> rlens;
  int64_t sum_rlens = 0;
  for (auto &aln : alns) {
    rlens[aln.rlen]++;
    sum_rlens += aln.rlen;
  }
  auto all_sum_rlens = upcxx::reduce_all(sum_rlens, op_fast_add).wait();
  auto all_nalns = upcxx::reduce_all(alns.size(), op_fast_add).wait();
  auto avg_rlen = all_sum_rlens / all_nalns;
  int most_common_rlen = avg_rlen;
  int64_t max_count = 0;
  for (auto &rlen : rlens) {
    if (rlen.second > max_count) {
      max_count = rlen.second;
      most_common_rlen = rlen.first;
    }
  }
  SLOG_VERBOSE("Computed unmerged read length as ", most_common_rlen, " with a count of ", max_count, " and average of ", avg_rlen,
               "\n");
  return most_common_rlen;
}

future<> Alns::sort_alns() {
  AsyncTimer timer(__FILEFUNC__);
  // execute this in a separate thread so master can continue to communicate freely
  auto fut = execute_in_thread_pool([&alns=this->alns, timer]() {
    timer.start();
    // sort the alns by name and then for the read from best score to worst - this is needed in later stages
    std::sort(alns.begin(), alns.end(), [](const Aln &elem1, const Aln &elem2) {
      if (elem1.read_id == elem2.read_id) {
        // sort by score, then contig len then last by cid to get a deterministic ordering
        if (elem1.score1 == elem2.score1) {
          if (elem1.clen == elem2.clen) return elem1.cid > elem2.cid;
          return elem1.clen > elem2.clen;
        }
        return elem1.score1 > elem2.score1;
      }
      if (elem1.read_id.length() == elem2.read_id.length()) {
        auto rlen = elem1.read_id.length();
        auto cmp = elem1.read_id.compare(0, rlen - 2, elem2.read_id, 0, rlen - 2);
        if (cmp == 0) return (elem1.read_id[rlen - 1] == '1');
        return cmp > 0;
      }
      return elem1.read_id > elem2.read_id;
    });
    timer.stop();
    return timer;
  });
  return fut.then([](AsyncTimer timer){
    // TODO record timer and initate reports after waiting...
  });
}
