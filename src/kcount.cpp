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

#include "kcount.hpp"

#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "upcxx_utils/mem_profile.hpp"

#include "utils.hpp"
#include "kmer_dht.hpp"
#include "packed_reads.hpp"
#include "contigs.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;



uint64_t estimate_num_kmers(unsigned kmer_len, vector<PackedReads*> &packed_reads_list) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  int64_t num_reads = 0;
  int64_t tot_num_reads = 0;
  for (auto packed_reads : packed_reads_list) {
    tot_num_reads += packed_reads->get_local_num_reads();
    packed_reads->reset();
    string id, seq, quals;
    ProgressBar progbar(packed_reads->get_local_num_reads(), "Scanning reads to estimate number of kmers");

    for (int i = 0; i < 100000; i++) {
      if (!packed_reads->get_next_read(id, seq, quals)) break;
      progbar.update();
      // do not read the entire data set for just an estimate
      if (seq.length() < kmer_len) continue;
      num_kmers += seq.length() - kmer_len + 1;
      num_reads++;
    }
    progbar.done();
    barrier();
  }
  DBG("This rank processed ", num_reads, " reads, and found ", num_kmers, " kmers\n");
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_tot_num_reads = reduce_one(tot_num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_all(num_kmers, op_fast_add).wait();

  SLOG_VERBOSE("Processed ", perc_str(all_num_reads, all_tot_num_reads), " reads, and estimated a maximum of ",
               all_num_kmers * (all_tot_num_reads / all_num_reads), " kmers\n");
  return num_reads > 0 ? num_kmers * tot_num_reads / num_reads : 0;
}


// Reduce compile time by instantiating templates of common types
// extern template declarations are in kcount.hpp
// template instantiations each happen in src/CMakeLists via kcount-extern-template.in.cpp

/*
 * This is called in kcount-extern-template.in.cpp
 * 
__MACRO_KCOUNT__(32, template);

#if MAX_BUILD_KMER >= 64

__MACRO_KCOUNT__(64,  template);

#endif
#if MAX_BUILD_KMER >= 96

__MACRO_KCOUNT__(96,  template);

#endif
#if MAX_BUILD_KMER >= 128

__MACRO_KCOUNT__(128,  template);

#endif
#if MAX_BUILD_KMER >= 160

__MACRO_KCOUNT__(160,  template);

#endif

 */