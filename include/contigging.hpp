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
#include "options.hpp"
#include "packed_reads.hpp"

template <int MAX_K>
void contigging(int kmer_len, int prev_kmer_len, int rlen_limit, std::vector<PackedReads *> &packed_reads_list, Contigs &ctgs,
                int &max_expected_ins_size, int &ins_avg, int &ins_stddev, std::shared_ptr<Options> options);

#define __MACRO_CONTIGGING__(KMER_LEN, MODIFIER)                                                                  \
  MODIFIER void contigging<KMER_LEN>(int, int, int, std::vector<PackedReads *> &, Contigs &, int &, int &, int &, \
                                     std::shared_ptr<Options>);

// Reduce compile time by instantiating templates of common types
// extern template declarations are in contigging.hpp
// template instantiations each happen in src/CMakeLists via contigging-extern-template.in.cpp

__MACRO_CONTIGGING__(32, extern template);
#if MAX_BUILD_KMER >= 64
__MACRO_CONTIGGING__(64, extern template);
#endif
#if MAX_BUILD_KMER >= 96
__MACRO_CONTIGGING__(96, extern template);
#endif
#if MAX_BUILD_KMER >= 128
__MACRO_CONTIGGING__(128, extern template);
#endif
#if MAX_BUILD_KMER >= 160
__MACRO_CONTIGGING__(160, extern template);
#endif
