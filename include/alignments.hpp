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

#include "contigs.hpp"
#include "version.h"

#include "upcxx/upcxx.hpp"

using upcxx::future;

struct Aln {
  Aln();
  Aln(const string &read_id, int64_t cid, int rstart, int rstop, int rlen, int cstart, int cstop, int clen, char orient,
      int score1 = 0, int score2 = 0, int identity = 0, int mismatches = 0, int read_group_id = -1);

  // optimal packing of data fields (does not match constructor exactly)
  string read_id;
  int64_t cid;
  int cstart, cstop, clen;
  int rstart, rstop, rlen;  // TODO can these be int16_t (for short reads only)?
  int score1, score2;       // TODO can this be uint16_t (for short reads only)?
  int mismatches;           // TODO can this be uint16_t (for short reads only)?
  string sam_string;
  int16_t read_group_id;
  char orient;      // TODO can this be bool?
  int8_t identity;  // this can be uint8_t - used as integer 0-100 %  // TODO can it even a member function int identity(int
                    // match_score)

  // writes out in the format meraligner uses
  string to_string() const;
  bool is_valid() const;
};

class Alns {
  vector<Aln> alns;
  int64_t num_dups;

 public:
  Alns();

  void clear();

  bool check_dup(Aln &aln);

  void add_aln(Aln &aln);

  void append(Alns &more_alns);

  const Aln &get_aln(int64_t i) const;

  Aln &get_aln(int64_t i);

  size_t size() const;

  void reserve(size_t capacity);

  void reset();

  int64_t get_num_dups();

  inline auto begin() { return alns.begin(); }
  inline auto end() { return alns.end(); };

  template <typename OSTREAM>
  void dump_all(OSTREAM &os, bool as_sam_format, int min_ctg_len = 0) const {
    // all ranks dump their valid alignments
    for (const auto &aln : alns) {
      if (aln.clen < min_ctg_len) continue;
      if (!as_sam_format)
        os << aln.to_string() << "\n";
      else
        os << aln.sam_string << "\n";
    }
  }

  void dump_single_file(const string fname) const;
  void dump_sam_file(const string fname, const vector<string> &read_group_names, const Contigs &ctgs, int min_contig_len = 0) const;
  void dump_rank_file(const string fname) const;

  int calculate_unmerged_rlen();

  future<> sort_alns();
};
