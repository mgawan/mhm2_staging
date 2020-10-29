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

#include <stdint.h>
#include <stdio.h>

#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#ifdef __APPLE__
#include <machine/endian.h>
#define H2BE htonll
#define BE2H ntohll
#else
#include <endian.h>
#define H2BE htobe64
#define BE2H be64toh
#endif

#include <upcxx/upcxx.hpp>

#include "hash_funcs.h"

/* Short description:
 *  - Store kmer strings by using 2 bits per base instead of 8
 *  - Easily return reverse complements of kmers, e.g. TTGG -> CCAA
 *  - Easily compare kmers
 *  - Provide hash of kmers
 *  - Get last and next kmer, e.g. ACGT -> CGTT or ACGT -> AACGT
 *  */
extern const uint64_t TWIN_TABLE[256];

template <int MAX_K>
class Kmer {
  inline static unsigned int k = 0;
  inline static const int N_LONGS = (MAX_K + 31) / 32;
  std::array<uint64_t, N_LONGS> longs;

 public:
  // serialization has to be public
  UPCXX_SERIALIZED_FIELDS(longs);

  Kmer() {
    assert(Kmer::k > 0);
    longs.fill(0);
  }

  Kmer(const Kmer &o) {
    assert(Kmer::k > 0);
    longs = o.longs;
  }

  explicit Kmer(const char *s) {
    assert(Kmer::k > 0);
    set_kmer(s);
  }

  static void set_k(unsigned int k) { Kmer::k = k; }

  static unsigned int get_k() { return Kmer::k; }

  static void get_kmers(unsigned kmer_len, std::string seq, std::vector<Kmer> &kmers) {
    // only need rank 0 to check
    assert(Kmer::k > 0);
    assert(kmer_len == Kmer::k);
    kmers.clear();
    if (seq.size() < Kmer::k) return;
    for (auto &c : seq) c = toupper(c);
    int bufsize = std::max((int)N_LONGS, (int)(seq.size() + 31) / 32) + 2;
    int lastLong = N_LONGS - 1;
    assert(lastLong >= 0 && lastLong < N_LONGS);
    kmers.resize(seq.size() - Kmer::k + 1);
    uint64_t buf[bufsize];
    uint8_t *bufPtr = (uint8_t *)buf;
    memset(buf, 0, bufsize * 8);
    const char *s = seq.c_str();
    // calculate binary along entire sequence
    for (unsigned i = 0; i < seq.size(); ++i) {
      int j = i % 32;
      int l = i / 32;
      assert(*s != '\0');
      size_t x = ((*s) & 4) >> 1;
      buf[l] |= ((x + ((x ^ (*s & 2)) >> 1)) << (2 * (31 - j)));
      s++;
    }
    // fix to big endian
    for (int l = 0; l < bufsize; l++) buf[l] = H2BE(buf[l]);
    const uint64_t mask = ((int64_t)0x3);
    uint64_t endmask = 0;
    if (Kmer::k % 32) {
      endmask = (((uint64_t)2) << (2 * (31 - (k % 32)) + 1)) - 1;
      // k == 0 :                0x0000000000000000
      // k == 1 : 2 << 61  - 1 : 0x3FFFFFFFFFFFFFFF
      // k == 31: 2 << 1   - 1 : 0x0000000000000003
    }
    endmask = ~endmask;
    // make 4 passes for each phase of 2 bits per base
    for (int shift = 0; shift < 4; shift++) {
      if (shift > 0) {
        // shift all bits by 2 in the buffer of longs
        for (int l = 0; l < bufsize - 1; l++) {
          buf[l] = (~mask & (BE2H(buf[l]) << 2)) | (mask & (BE2H(buf[l + 1]) >> 62));
          buf[l] = H2BE(buf[l]);
        }
      }
      // enumerate the kmers in the phase
      for (unsigned i = shift; i < kmers.size(); i += 4) {
        int byteOffset = i / 4;
        assert(byteOffset + N_LONGS * 8 <= bufsize * 8);
        for (int l = 0; l < N_LONGS; l++) {
          kmers[i].longs[l] = BE2H(*((uint64_t *)(bufPtr + byteOffset + l * 8)));
        }
        // set remaining bits to 0
        kmers[i].longs[lastLong] &= endmask;
        //        for (int l = N_LONGS; l < N_LONGS; l++) {
        //          kmers[i].longs[l] = 0;
        //        }
      }
    }
  }

  Kmer &operator=(const Kmer &o) {
    if (this != &o) longs = o.longs;
    return *this;
  }

  bool operator<(const Kmer &o) const {
    for (size_t i = 0; i < N_LONGS; ++i) {
      if (longs[i] < o.longs[i]) return true;
      if (longs[i] > o.longs[i]) return false;
    }
    return false;
  }

  bool operator==(const Kmer &o) const { return longs == o.longs; }

  bool operator!=(const Kmer &o) const { return !(*this == o); }

  void set_kmer(const char *s) {
    size_t i, j, l;
    longs.fill(0);
    for (i = 0; i < Kmer::k; ++i) {
      j = i % 32;
      l = i / 32;
      assert(*s != '\0');
      size_t x = ((*s) & 4) >> 1;
      longs[l] |= ((x + ((x ^ (*s & 2)) >> 1)) << (2 * (31 - j)));
      s++;
    }
  }

  uint64_t hash() const { return MurmurHash3_x64_64(reinterpret_cast<const void *>(longs.data()), N_LONGS * sizeof(uint64_t)); }

  Kmer revcomp() const {
    Kmer km(*this);
    for (size_t i = 0; i < N_LONGS; i++) {
      uint64_t v = longs[i];
      km.longs[N_LONGS - 1 - i] = (TWIN_TABLE[v & 0xFF] << 56) | (TWIN_TABLE[(v >> 8) & 0xFF] << 48) |
                                  (TWIN_TABLE[(v >> 16) & 0xFF] << 40) | (TWIN_TABLE[(v >> 24) & 0xFF] << 32) |
                                  (TWIN_TABLE[(v >> 32) & 0xFF] << 24) | (TWIN_TABLE[(v >> 40) & 0xFF] << 16) |
                                  (TWIN_TABLE[(v >> 48) & 0xFF] << 8) | (TWIN_TABLE[(v >> 56)]);
    }
    size_t shift = (Kmer::k % 32) ? 2 * (32 - (Kmer::k % 32)) : 0;
    uint64_t shiftmask = (Kmer::k % 32) ? (((1ULL << shift) - 1) << (64 - shift)) : 0ULL;
    km.longs[0] = km.longs[0] << shift;
    for (size_t i = 1; i < N_LONGS; i++) {
      km.longs[i - 1] |= (km.longs[i] & shiftmask) >> (64 - shift);
      km.longs[i] = km.longs[i] << shift;
    }
    return km;
  }

  Kmer forward_base(const char b) const {
    Kmer km(*this);
    km.longs[0] = km.longs[0] << 2;
    for (size_t i = 1; i < N_LONGS; i++) {
      km.longs[i - 1] |= (km.longs[i] & (3ULL << 62)) >> 62;
      km.longs[i] = km.longs[i] << 2;
    }
    uint64_t x = (b & 4) >> 1;
    km.longs[N_LONGS - 1] |= (x + ((x ^ (b & 2)) >> 1)) << (2 * (31 - ((k - 1) % 32)));
    return km;
  }

  Kmer backward_base(const char b) const {
    Kmer km(*this);
    km.longs[N_LONGS - 1] = km.longs[N_LONGS - 1] >> 2;
    km.longs[N_LONGS - 1] &= (k % 32) ? (((1ULL << (2 * (k % 32))) - 1) << 2 * (32 - (k % 32))) : ~0ULL;
    for (size_t i = 1; i < N_LONGS; i++) {
      km.longs[N_LONGS - i] |= (km.longs[N_LONGS - i - 1] & 3ULL) << 62;
      km.longs[N_LONGS - i - 1] = km.longs[N_LONGS - i - 1] >> 2;
    }
    uint64_t x = (b & 4) >> 1;
    km.longs[0] |= (x + ((x ^ (b & 2)) >> 1)) << 62;
    return km;
  }

  void to_string(char *s) const {
    size_t i, j, l;
    for (i = 0; i < Kmer::k; i++) {
      j = i % 32;
      l = i / 32;
      switch (((longs[l]) >> (2 * (31 - j))) & 0x03) {
        case 0x00:
          *s = 'A';
          ++s;
          break;
        case 0x01:
          *s = 'C';
          ++s;
          break;
        case 0x02:
          *s = 'G';
          ++s;
          break;
        case 0x03:
          *s = 'T';
          ++s;
          break;
      }
    }
    *s = '\0';
  }

  std::string to_string() const {
    char buf[Kmer::k + 1];
    to_string(buf);
    return std::string(buf);
  }

  std::pair<const uint8_t *, int> get_bytes() const {
    return {reinterpret_cast<const uint8_t *>(longs.data()), N_LONGS * sizeof(uint64_t)};
  }
};

template <int MAX_K>
struct KmerHash {
  size_t operator()(const Kmer<MAX_K> &km) const { return km.hash(); }
};

template <int MAX_K>
struct KmerEqual {
  size_t operator()(const Kmer<MAX_K> &k1, const Kmer<MAX_K> &k2) const { return k1 == k2; }
};

// specialization of std::Hash

namespace std {

template <int MAX_K>
struct hash<Kmer<MAX_K>> {
  typedef std::size_t result_type;
  result_type operator()(Kmer<MAX_K> const &km) const { return km.hash(); }
};
}  // namespace std

template <int MAX_K>
inline std::ostream &operator<<(std::ostream &out, const Kmer<MAX_K> &k) {
  return out << k.to_string();
};

#define __MACRO_KMER__(KMER_LEN, MODIFIER) \
  MODIFIER class Kmer<KMER_LEN>;           \
  MODIFIER struct KmerHash<KMER_LEN>;      \
  MODIFIER struct KmerEqual<KMER_LEN>;

// Reduce compile time by instantiating templates of common types
// extern template declarations are in kmer.hpp
// template instantiations each happen in src/CMakeLists via kmer-extern-template.in.cpp

__MACRO_KMER__(32, extern template);

#if MAX_BUILD_KMER >= 64

__MACRO_KMER__(64, extern template);

#endif
#if MAX_BUILD_KMER >= 96

__MACRO_KMER__(96, extern template);

#endif
#if MAX_BUILD_KMER >= 128

__MACRO_KMER__(128, extern template);

#endif
#if MAX_BUILD_KMER >= 160

__MACRO_KMER__(160, extern template);

#endif
