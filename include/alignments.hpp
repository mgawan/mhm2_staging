#pragma once

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
#include <upcxx/upcxx.hpp>
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "zstr.hpp"

using namespace upcxx_utils;

struct Aln {
  string read_id;
  int64_t cid;
  int rstart, rstop, rlen, cstart, cstop, clen;
  char orient;
  int score1, score2;
  string sam_string;

  // writes out in the format meraligner uses
  string to_string() {
    ostringstream os;
    os << read_id << "\t" << rstart + 1 << "\t" << rstop << "\t" << rlen << "\t"
       << "Contig" << cid << "\t" << cstart + 1 << "\t" << cstop << "\t" << clen << "\t"
       << (orient == '+' ? "Plus" : "Minus") << "\t" << score1 << "\t" << score2;
    return os.str();
  }
};


class Alns {

  vector<Aln> alns;
  int64_t num_dups;

public:

  Alns() : num_dups(0) {}

  void clear() {
    alns.clear();
    vector<Aln>().swap(alns);
  }

  void add_aln(Aln &aln) {
    // check for duplicate alns first - do this backwards because only the most recent entries could be for this read
    for (auto it = alns.rbegin(); it != alns.rend(); ++it) {
      // we have no more entries for this read
      if (it->read_id != aln.read_id) break;
      // now check for equality
      if (it->rstart == aln.rstart && it->rstop == aln.rstop && it->cstart == aln.cstart && it->cstop == aln.cstop) {
        num_dups++;
        return;
      }
    }
    alns.push_back(aln);
  }

  Aln &get_aln(int64_t i) {
    return alns[i];
  }

  size_t size() {
    return alns.size();
  }

  int64_t get_num_dups() {
    return upcxx::reduce_one(num_dups, upcxx::op_fast_add, 0).wait();
  }

  auto begin() {
    return alns.begin();
  }

  auto end() {
    return alns.end();
  }

  void dump_alns(string fname) {
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

  void dump_single_file_alns(const string fname, bool as_sam_format=false) {
    BarrierTimer timer(__FILEFUNC__, false, true);

    string out_str = "";
    for (auto aln : alns) {
      if (!as_sam_format) out_str += aln.to_string() + "\n";
      else out_str += aln.sam_string + "\n";
    }
    dump_single_file(fname, out_str);
  }
};


