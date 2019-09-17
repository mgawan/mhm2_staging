#ifndef __BLOOM
#define __BLOOM

#include <vector>
#include <cmath>
#include <exception>

#include "hash_funcs.h"


inline std::array<uint64_t, 2> bloom_hash(const uint8_t *data, std::size_t len) {
  std::array<uint64_t, 2> hashValue;
  MurmurHash3_x64_128(data, len, 0, hashValue.data());
  return hashValue;
}

inline uint64_t nth_hash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize) {
  return (hashA + n * hashB) % filterSize;
}
  

struct BloomFilter {

  int num_hashes;
  std::vector<bool> m_bits;

public:
  
  void init(uint64_t entries, double error) {
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
    SLOG_VERBOSE("Rank 0 created bloom filter with ", num_bits, " bits and ", num_hashes, " hashes\n");
  }

  void clear() {
    std::vector<bool>().swap(m_bits);
  }
  
  void add(const uint8_t *data, std::size_t len) {
    auto hash_values = bloom_hash(data, len);
    for (int n = 0; n < num_hashes; n++) 
      m_bits[nth_hash(n, hash_values[0], hash_values[1], m_bits.size())] = true;
  }

  bool possibly_contains(const uint8_t *data, std::size_t len) const {
    auto hash_values = bloom_hash(data, len);
    for (int n = 0; n < num_hashes; n++) 
      if (!m_bits[nth_hash(n, hash_values[0], hash_values[1], m_bits.size())]) return false;
    return true;
  }

  size_t estimate_num_items() const {
    size_t bits_on = 0, m = m_bits.size(), k = num_hashes;
    for (auto it = m_bits.begin(); it != m_bits.end(); it++) {
        if (*it) bits_on++;
    }
    return (size_t) (- ((double) m/ (double) k) * log( 1.0 - ((double) bits_on / (double) m) ) + 0.5);
  }

  bool is_initialized() const {
    return m_bits.size() != 0;
  }

};


#endif
