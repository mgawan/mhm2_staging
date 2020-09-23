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

#include "kmer_dht.hpp"

double _dynamic_min_depth = 0;
int _dmin_thres = 2.0;

array<pair<char, int>, 4> ExtCounts::get_sorted() {
  array<pair<char, int>, 4> counts = {make_pair('A', (int)count_A), make_pair('C', (int)count_C), make_pair('G', (int)count_G),
                                      make_pair('T', (int)count_T)};
  sort(std::begin(counts), std::end(counts), [](const auto &elem1, const auto &elem2) {
    if (elem1.second == elem2.second)
      return elem1.first > elem2.first;
    else
      return elem1.second > elem2.second;
  });
  return counts;
}

bool ExtCounts::is_zero() {
  if (count_A + count_C + count_G + count_T == 0) return true;
  return false;
}

void ExtCounts::inc(char ext, int count) {
  switch (ext) {
    case 'A':
      count += count_A;
      count_A = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
      break;
    case 'C':
      count += count_C;
      count_C = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
      break;
    case 'G':
      count += count_G;
      count_G = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
      break;
    case 'T':
      count += count_T;
      count_T = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
      break;
  }
}


char ExtCounts::get_ext(uint16_t count) {
  auto sorted_counts = get_sorted();
  int top_count = sorted_counts[0].second;
  int runner_up_count = sorted_counts[1].second;
  // set dynamic_min_depth to 1.0 for single depth data (non-metagenomes)
  int dmin_dyn = max((int)((1.0 - _dynamic_min_depth) * count), _dmin_thres);
  if (top_count < dmin_dyn) return 'X';
  if (runner_up_count >= dmin_dyn) return 'F';
  return sorted_counts[0].first;
  /*
  // FIXME: this is not very helpful. With qual_cutoff = 20) it increases ctgy & coverage a little bit, but at a cost
  // of increased msa. We really need to try both low q (qual cutoff 10) and hi q, as we do with localassm.
  double dmin_dyn = max(2.0, LASSM_MIN_EXPECTED_DEPTH * count);
  if ((top_count < dmin_dyn && runner_up_count > 0) || (top_count >= dmin_dyn && runner_up_count >= dmin_dyn)) return 'F';
  return sorted_counts[0].first;
  */
}

// Reduce compile time by instantiating templates of common types
// extern template declarations are in kmer_dht.hpp
// template instantiations each happen in src/CMakeLists via kmer_dht-extern-template.in.cpp

/*
 * These are now instantiated in kmer_dht-extern-template.in.cpp
 *

__MACRO_KMER_DHT__(32, template);

#if MAX_BUILD_KMER >= 64

__MACRO_KMER_DHT__(64,  template);

#endif
#if MAX_BUILD_KMER >= 96

__MACRO_KMER_DHT__(96,  template);

#endif
#if MAX_BUILD_KMER >= 128

__MACRO_KMER_DHT__(128,  template);

#endif
#if MAX_BUILD_KMER >= 160

__MACRO_KMER_DHT__(160,  template);

#endif

*/