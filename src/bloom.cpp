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

#include "bloom.hpp"

#include <vector>
#include <cmath>
#include <exception>

#include "hash_funcs.h"

#include "upcxx_utils/log.hpp"

using namespace upcxx_utils;

std::array<uint64_t, 2> bloom_hash(const std::pair<const uint8_t*, std::size_t> &data) {
  std::array<uint64_t, 2> hashValue;
  MurmurHash3_x64_128(data.first, data.second, 0, hashValue.data());
  return hashValue;
}

uint64_t nth_hash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize) {
  return (hashA + n * hashB) % filterSize;
}



  void BloomFilter::init(uint64_t entries, double error) {
    double num = log(error);
    double denom = 0.480453013918201; // ln(2)^2
    double bpe = -(num / denom);

    double dentries = (double)entries;
    uint64_t num_bits = (uint64_t)(dentries * bpe);
    num_hashes = (int)ceil(0.693147180559945 * bpe); // ln(2)
    try {
      m_bits.resize(num_bits, false);
    } catch (std::exception &e) {
      DIE(e.what(), " note: num bits is ", num_bits, " dentries is ", dentries, " bpe is ", bpe);
    }
    SLOG_VERBOSE("Rank 0 created bloom filter with ", num_bits, " bits and ", num_hashes, " hashes (", get_size_str(num_bits/8), ")\n");
  }

  void BloomFilter::clear() {
    std::vector<bool>().swap(m_bits);
  }

  void BloomFilter::add(const std::pair<const uint8_t*, std::size_t> &data) {
    auto hash_values = bloom_hash(data);
    for (int n = 0; n < num_hashes; n++)
      m_bits[nth_hash(n, hash_values[0], hash_values[1], m_bits.size())] = true;
  }

  bool BloomFilter::possibly_contains(const std::pair<const uint8_t*, std::size_t> &data) const {
    auto hash_values = bloom_hash(data);
    for (int n = 0; n < num_hashes; n++)
      if (!m_bits[nth_hash(n, hash_values[0], hash_values[1], m_bits.size())]) return false;
    return true;
  }

  size_t BloomFilter::estimate_num_items() const {
    size_t bits_on = 0, m = m_bits.size(), k = num_hashes;
    for (auto it = m_bits.begin(); it != m_bits.end(); it++) {
        if (*it) bits_on++;
    }
    return (size_t) (- ((double) m/ (double) k) * log( 1.0 - ((double) bits_on / (double) m) ) + 0.5);
  }

  bool BloomFilter::is_initialized() const {
    return m_bits.size() != 0;
  }



