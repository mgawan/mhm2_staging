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

#include "alignments.hpp"
#include "contigs.hpp"
#include "packed_reads.hpp"

#ifdef ENABLE_GPUS
#include "adept-sw/driver.hpp"
#endif

struct AlnScoring {
  int match, mismatch, gap_opening, gap_extending, ambiguity;

  string to_string() {
    std::ostringstream oss;
    oss << "match " << match << " mismatch " << mismatch << " gap open " << gap_opening << " gap extend " << gap_extending
        << " ambiguity " << ambiguity;
    return oss.str();
  }
};

template <int MAX_K>
double find_alignments(unsigned kmer_len, std::vector<PackedReads *> &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                       Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit, bool compute_cigar = false, int min_ctg_len = 0,
                       int ranks_per_gpu = 0);

// Reduce compile time by instantiating templates of common types
// extern template declarations are in kmer.hpp
// template instantiations each happen in src/CMakeLists via kmer-extern-template.in.cpp

#define __MACRO_KLIGN__(KMER_LEN, MODIFIER)                                                                                      \
  MODIFIER double find_alignments<KMER_LEN>(unsigned, std::vector<PackedReads *> &, int, int, Contigs &, Alns &, int, int, bool, \
                                            int, int);

__MACRO_KLIGN__(32, extern template);
#if MAX_BUILD_KMER >= 64
__MACRO_KLIGN__(64, extern template);
#endif
#if MAX_BUILD_KMER >= 96
__MACRO_KLIGN__(96, extern template);
#endif
#if MAX_BUILD_KMER >= 128
__MACRO_KLIGN__(128, extern template);
#endif
#if MAX_BUILD_KMER >= 160
__MACRO_KLIGN__(160, extern template);
#endif
