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

#include <sstream>
#include <string>
#include <upcxx/upcxx.hpp>

#include "contigs.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "version.h"
#include "zstr.hpp"

using namespace upcxx_utils;
using std::ostringstream;
using std::string;
using std::to_string;

//
// class Aln
//

// writes out in the format meraligner uses
string Aln::to_string() {
  ostringstream os;
  os << read_id << "\t" << rstart + 1 << "\t" << rstop << "\t" << rlen << "\t"
     << "Contig" << cid << "\t" << cstart + 1 << "\t" << cstop << "\t" << clen << "\t" << (orient == '+' ? "Plus" : "Minus") << "\t"
     << score1 << "\t" << score2;
  return os.str();
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
  alns.push_back(aln);
}

Aln &Alns::get_aln(int64_t i) { return alns[i]; }

size_t Alns::size() { return alns.size(); }

int64_t Alns::get_num_dups() { return upcxx::reduce_one(num_dups, upcxx::op_fast_add, 0).wait(); }

void Alns::dump_alns(string fname) {
  get_rank_path(fname, rank_me());
  zstr::ofstream f(fname);
  ostringstream out_buf;
  ProgressBar progbar(alns.size(), "Writing alns");
  size_t bytes_written = 0;
  int64_t i = 0;
  for (auto aln : alns) {
    out_buf << aln.to_string() << std::endl;
    progbar.update();
    i++;
    if (!(i % 1000)) {
      f << out_buf.str();
      bytes_written += out_buf.str().length();
      out_buf = ostringstream();
    }
  }
  if (!out_buf.str().empty()) {
    f << out_buf.str();
    bytes_written += out_buf.str().length();
  }
  f.close();
  progbar.done();
  upcxx::barrier();
}

void Alns::dump_single_file_alns(const string fname, bool as_sam_format, Contigs *ctgs) {
  BarrierTimer timer(__FILEFUNC__);

  string out_str = "";

  // FIXME: first, all ranks must dump contig info to the file, for every contig:
  // @SQ	SN:Contig0	LN:887
  dist_ofstream of(fname);
  future<> all_done;
  for (auto &ctg : *ctgs) {
    of << "@SQ\tSN:Contig" << to_string(ctg.id) << "\tLN:" << to_string(ctg.seq.length()) << "\n";
  }
  all_done = of.flush_collective();

  if (!upcxx::rank_me()) {
    // program information
    of << "@PG\tID:MHM2\tPN:MHM2\tVN:" << string(MHM2_VERSION) << "\n";
  }
  all_done = when_all(all_done, of.flush_collective());

  for (auto aln : alns) {
    if (!as_sam_format)
      of << aln.to_string() << "\n";
    else
      of << aln.sam_string << "\n";
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

void Alns::sort_alns() {
  BarrierTimer timer(__FILEFUNC__);
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
}
