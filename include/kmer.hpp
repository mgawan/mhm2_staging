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
using longs_t = uint64_t;

inline const longs_t TWIN_TABLE[256] = {
    0xFF, 0xBF, 0x7F, 0x3F, 0xEF, 0xAF, 0x6F, 0x2F, 0xDF, 0x9F, 0x5F, 0x1F, 0xCF, 0x8F, 0x4F, 0x0F, 0xFB, 0xBB, 0x7B, 0x3B,
    0xEB, 0xAB, 0x6B, 0x2B, 0xDB, 0x9B, 0x5B, 0x1B, 0xCB, 0x8B, 0x4B, 0x0B, 0xF7, 0xB7, 0x77, 0x37, 0xE7, 0xA7, 0x67, 0x27,
    0xD7, 0x97, 0x57, 0x17, 0xC7, 0x87, 0x47, 0x07, 0xF3, 0xB3, 0x73, 0x33, 0xE3, 0xA3, 0x63, 0x23, 0xD3, 0x93, 0x53, 0x13,
    0xC3, 0x83, 0x43, 0x03, 0xFE, 0xBE, 0x7E, 0x3E, 0xEE, 0xAE, 0x6E, 0x2E, 0xDE, 0x9E, 0x5E, 0x1E, 0xCE, 0x8E, 0x4E, 0x0E,
    0xFA, 0xBA, 0x7A, 0x3A, 0xEA, 0xAA, 0x6A, 0x2A, 0xDA, 0x9A, 0x5A, 0x1A, 0xCA, 0x8A, 0x4A, 0x0A, 0xF6, 0xB6, 0x76, 0x36,
    0xE6, 0xA6, 0x66, 0x26, 0xD6, 0x96, 0x56, 0x16, 0xC6, 0x86, 0x46, 0x06, 0xF2, 0xB2, 0x72, 0x32, 0xE2, 0xA2, 0x62, 0x22,
    0xD2, 0x92, 0x52, 0x12, 0xC2, 0x82, 0x42, 0x02, 0xFD, 0xBD, 0x7D, 0x3D, 0xED, 0xAD, 0x6D, 0x2D, 0xDD, 0x9D, 0x5D, 0x1D,
    0xCD, 0x8D, 0x4D, 0x0D, 0xF9, 0xB9, 0x79, 0x39, 0xE9, 0xA9, 0x69, 0x29, 0xD9, 0x99, 0x59, 0x19, 0xC9, 0x89, 0x49, 0x09,
    0xF5, 0xB5, 0x75, 0x35, 0xE5, 0xA5, 0x65, 0x25, 0xD5, 0x95, 0x55, 0x15, 0xC5, 0x85, 0x45, 0x05, 0xF1, 0xB1, 0x71, 0x31,
    0xE1, 0xA1, 0x61, 0x21, 0xD1, 0x91, 0x51, 0x11, 0xC1, 0x81, 0x41, 0x01, 0xFC, 0xBC, 0x7C, 0x3C, 0xEC, 0xAC, 0x6C, 0x2C,
    0xDC, 0x9C, 0x5C, 0x1C, 0xCC, 0x8C, 0x4C, 0x0C, 0xF8, 0xB8, 0x78, 0x38, 0xE8, 0xA8, 0x68, 0x28, 0xD8, 0x98, 0x58, 0x18,
    0xC8, 0x88, 0x48, 0x08, 0xF4, 0xB4, 0x74, 0x34, 0xE4, 0xA4, 0x64, 0x24, 0xD4, 0x94, 0x54, 0x14, 0xC4, 0x84, 0x44, 0x04,
    0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20, 0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00};

inline const longs_t ZERO_MASK[32] = {
    0x0000000000000000, 0xC000000000000000, 0xF000000000000000, 0xFC00000000000000, 0xFF00000000000000, 0xFFC0000000000000,
    0xFFF0000000000000, 0xFFFC000000000000, 0xFFFF000000000000, 0xFFFFC00000000000, 0xFFFFF00000000000, 0xFFFFFC0000000000,
    0xFFFFFF0000000000, 0xFFFFFFC000000000, 0xFFFFFFF000000000, 0xFFFFFFFC00000000, 0xFFFFFFFF00000000, 0xFFFFFFFFC0000000,
    0xFFFFFFFFF0000000, 0xFFFFFFFFFC000000, 0xFFFFFFFFFF000000, 0xFFFFFFFFFFC00000, 0xFFFFFFFFFFF00000, 0xFFFFFFFFFFFC0000,
    0xFFFFFFFFFFFF0000, 0xFFFFFFFFFFFFC000, 0xFFFFFFFFFFFFF000, 0xFFFFFFFFFFFFFC00, 0xFFFFFFFFFFFFFF00, 0xFFFFFFFFFFFFFFC0,
    0xFFFFFFFFFFFFFFF0, 0xFFFFFFFFFFFFFFFC};

template <int MAX_K>
class Kmer {
  inline static unsigned int k = 0;
  inline static const int N_LONGS = (MAX_K + 31) / 32;
  std::array<longs_t, N_LONGS> longs;

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

  static unsigned int get_N_LONGS() { return Kmer::N_LONGS; }

  static unsigned int get_MAX_K() { return MAX_K; }

  static void get_kmers(unsigned kmer_len, std::string seq, std::vector<Kmer> &kmers) {
    // only need rank 0 to check
    assert(Kmer::k > 0);
    assert(kmer_len == Kmer::k);
    kmers.clear();
    if (seq.size() < Kmer::k) return;
    for (auto &c : seq) c = toupper(c);
    int bufsize = std::max((int)N_LONGS, (int)(seq.size() + 31) / 32) + N_LONGS;
    int lastLong = N_LONGS - 1;
    assert(lastLong >= 0 && lastLong < N_LONGS);
    kmers.resize(seq.size() - Kmer::k + 1);
    longs_t buf[bufsize];
    uint8_t *bufPtr = (uint8_t *)buf;
    memset(buf, 0, bufsize * 8);
    const char *s = seq.c_str();
    // calculate binary along entire sequence
    for (unsigned i = 0; i < seq.size(); ++i) {
      int j = i % 32;
      int l = i / 32;
      assert(*s != '\0');
      longs_t x = ((*s) & 4) >> 1;
      buf[l] |= ((x + ((x ^ (*s & 2)) >> 1)) << (2 * (31 - j)));
      s++;
    }
    // fix to big endian
    for (int l = 0; l < bufsize; l++) buf[l] = H2BE(buf[l]);
    const longs_t mask = ((int64_t)0x3);
    longs_t endmask = 0;
    if (Kmer::k % 32) {
      endmask = (((longs_t)2) << (2 * (31 - (k % 32)) + 1)) - 1;
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
          kmers[i].longs[l] = BE2H(*((longs_t *)(bufPtr + byteOffset + l * 8)));
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
      assert(((*s >= 'A' && *s <= 'Z') || (*s >= 'a' && *s <= 'z')) && "bases are letters");
      longs_t x;
#if 1
      x = ((*s) & 4) >> 1;
      longs[l] |= ((x + ((x ^ (*s & 2)) >> 1)) << (2 * (31 - j)));
#else
      // This is the same, but broken down...
      x = ((*s) & 4) >> 1;  // i.e. Gg/Tt will set bit 1, so x == 2 | 0
      assert(x == 2 || x == 0);
      x |= (x ^ ((*s) & 2)) >> 1;  // i.e. Cc/Gg will not set bit 0, Aa/Tt will set bit 0
      assert(x >= 0 && x <= 3);
      longs[l] |= (x) << (2 * (31 - j));
#endif
      s++;
    }
  }

  std::string mer_to_string(longs_t mmer, int m) {
    char buf[33] = "";
    char *s = buf;
    for (int j = 0; j < m; j++) {
      switch ((mmer >> (2 * (31 - j))) & 0x03) {
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
    return std::string(buf);
  }
  
  std::string get_minimizer(int m) {
    char s[200];
    to_string(s);
    char *min_s = s;
    for (int i = 1; i <= Kmer::k - m; i++) {
      if (strncmp(s + i, min_s, m) > 0) min_s = s + i;
    }
    min_s[m] = '\0';
    return std::string(min_s);
  }

  std::string get_minimizer_opt(int m) {
    uint64_t minimizer = 0;
//    std::cout << std::endl;
    for (int i = 0; i <= Kmer::k - m; i++) {
      int j = i % 32;
      int l = i / 32;
      longs_t mmer = ((longs[l]) << (2 * j)) & ZERO_MASK[m];
//      std::cout << std::string(i, ' ') << mer_to_string(mmer, m) << "\n";
      if (mmer > minimizer) minimizer = mmer;
    }
//    std::cout << minimizer << " ";
    return mer_to_string(minimizer, m);
  }

  uint64_t minimizer_hash(int m) const {
    char s[200];
    to_string(s);
    char *min_s = s;
    for (int i = 1; i <= Kmer::k - m; i++) {
      if (strncmp(s + i, min_s, m) > 0) min_s = s + i;
    }
    min_s[m] = '\0';
    return MurmurHash3_x64_64(reinterpret_cast<const void *>(min_s), m);
  }

  uint64_t minimizer_hash_opt(int m) const {
    uint64_t minimizer = 0;
    for (int i = 0; i <= Kmer::k - m; i++) {
      int j = i % 32;
      int l = i / 32;
      longs_t mmer = ((longs[l]) << (2 * j)) & ZERO_MASK[m];
      if (mmer > minimizer) minimizer = mmer;
    }
    return MurmurHash3_x64_64(reinterpret_cast<const void *>(&minimizer), sizeof(minimizer));
  }

  uint64_t hash() const { return MurmurHash3_x64_64(reinterpret_cast<const void *>(longs.data()), N_LONGS * sizeof(longs_t)); }

  void set_zeros() {
    // set trailing bits in longs
    auto mod = k % 32;
    auto last_long = k / 32;
    if (mod != 0) {
      longs[last_long] &= ZERO_MASK[mod];
    }

    // set remaining longs, if any, to 0
    for (int l = 1 + last_long; l < N_LONGS; l++) {
      longs[l] = 0;
    }
  }

  Kmer revcomp() const {
    Kmer km;
    auto last_long = (k + 31) / 32;
    assert(last_long <= N_LONGS);
    for (size_t i = 0; i < last_long; i++) {
      longs_t v = longs[i];
      km.longs[last_long - 1 - i] = (TWIN_TABLE[v & 0xFF] << 56) | (TWIN_TABLE[(v >> 8) & 0xFF] << 48) |
                                    (TWIN_TABLE[(v >> 16) & 0xFF] << 40) | (TWIN_TABLE[(v >> 24) & 0xFF] << 32) |
                                    (TWIN_TABLE[(v >> 32) & 0xFF] << 24) | (TWIN_TABLE[(v >> 40) & 0xFF] << 16) |
                                    (TWIN_TABLE[(v >> 48) & 0xFF] << 8) | (TWIN_TABLE[(v >> 56)]);
    }
    longs_t shift = (Kmer::k % 32) ? 2 * (32 - (Kmer::k % 32)) : 0;
    longs_t shiftmask = (Kmer::k % 32) ? (((((longs_t)1) << shift) - 1) << (64 - shift)) : ((longs_t)0);
    km.longs[0] = km.longs[0] << shift;
    for (size_t i = 1; i < last_long; i++) {
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
    longs_t x = (b & 4) >> 1;
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
    longs_t x = (b & 4) >> 1;
    km.longs[0] |= (x + ((x ^ (b & 2)) >> 1)) << 62;
    return km;
  }

  char front() {
    switch (((longs[0]) >> (2 * (31))) & 0x03) {
      case 0x00: return 'A';
      case 0x01: return 'C';
      case 0x02: return 'G';
      case 0x03: return 'T';
      default: return 'N';
    }
  }

  char back() {
    size_t i = Kmer::k - 1;
    int j = i % 32;
    int l = i / 32;
    switch (((longs[l]) >> (2 * (31 - j))) & 0x03) {
      case 0x00: return 'A';
      case 0x01: return 'C';
      case 0x02: return 'G';
      case 0x03: return 'T';
      default: return 'N';
    }
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

  std::string to_hex() const {
    std::ostringstream os;
    assert(N_LONGS == longs.size());
    os << "longs=" << N_LONGS << ",k=" << get_k() << ":" << to_string() << ":";
    for (auto l : longs) os << l << ",";
    os << ":";
    os << std::hex;
    for (auto l : longs) {
      os << l << ",";
    }
    return os.str();
  }

  std::pair<const uint8_t *, int> get_bytes() const {
    return {reinterpret_cast<const uint8_t *>(longs.data()), N_LONGS * sizeof(longs_t)};
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
