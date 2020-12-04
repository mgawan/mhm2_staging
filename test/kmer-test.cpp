#include "kmer.hpp"
#include "gtest/gtest.h"

#include <map>
#include <string>
#include <vector>
using std::map;
using std::string;
using std::vector;

const char *As = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
const char *Cs = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                 "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
const char *Gs = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                 "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
const char *Ts = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                 "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
const char *ACGTs = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
                    "GTACGTACGTACGTACGTACGTACGTACGT";
const char *TCGAs = "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
                    "GATCGATCGATCGATCGATCGATCGATCGA";
const char *CAGTs = "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCA"
                    "GTCAGTCAGTCAGTCAGTCAGTCAGTCAGT";
const char *RandomRead = "CGCTGTTCCAGATGACGAACCAGGAATTCCGCCAGGTATTCGACTTTATTCGCGAAGTCAAGAAGTTGAACGTCATCAGTGTGAACTACGGTTGCGAAGGCTTCC"
                         "TCGGCAGCTACGAGAAGGATGCACGCATCTGCCCGTTCTTCTGCCGTGCCGGCGTGAACGTGTCCTCGGTGCTTTGCGATGGCAGCATTTCGGCATGCCCGAGC"
                         "T";

string slowrevcomp(const string &seq) {
  string revcomp;
  revcomp.reserve(seq.size());
  for (auto it = seq.rbegin(); it != seq.rend(); it++) {
    char rc = 'N';
    switch (*it) {
      case 'A': rc = 'T'; break;
      case 'C': rc = 'G'; break;
      case 'G': rc = 'C'; break;
      case 'T': rc = 'A'; break;
    }
    revcomp.push_back(rc);
  }
  return revcomp;
}

template <int MAX_K>
void test_get_kmers(int klen, string seq) {
  using K = Kmer<MAX_K>;
  vector<K> vec;
  K::get_kmers(klen, seq, vec);
  EXPECT_EQ(vec.size(), seq.size() - klen + 1) << "Correct num of kmers";
  for (int i = 0; i < vec.size(); i++) {
    auto kstr = vec[i].to_string();
    auto sstr = seq.substr(i, klen);
    EXPECT_STREQ(kstr.c_str(), sstr.c_str()) << "kmers eq string " << kstr << " v " << sstr << " at " << i;
  }
}

template <int MAX_K>
void test_get_kmers(int klen) {
  test_get_kmers<MAX_K>(klen, As);
  test_get_kmers<MAX_K>(klen, Cs);
  test_get_kmers<MAX_K>(klen, Gs);
  test_get_kmers<MAX_K>(klen, Ts);
  test_get_kmers<MAX_K>(klen, ACGTs);
  test_get_kmers<MAX_K>(klen, TCGAs);
  test_get_kmers<MAX_K>(klen, CAGTs);
  test_get_kmers<MAX_K>(klen, RandomRead);
}

template <typename K>
void check_kmer(const K &kmer, string seq) {
  EXPECT_EQ(sizeof(longs_t), 8) << "longs_t is 8 bytes";

  EXPECT_EQ(K::get_N_LONGS(), (K::get_MAX_K() + 31) / 32)
      << "NLONGS for k=" << kmer.get_k() << " is " << K::get_N_LONGS() << " v " << (K::get_MAX_K() + 31) / 32;
  EXPECT_EQ(K::get_N_LONGS(), kmer.get_N_LONGS());

  EXPECT_EQ(seq.size(), kmer.get_k());
  EXPECT_STREQ(seq.c_str(), kmer.to_string().c_str());

  K rev = kmer.revcomp();
  auto kbytes = kmer.get_bytes();
  auto nbytes = ((K::get_MAX_K() + 31) / 32) * 8;
  EXPECT_EQ(kbytes.second, nbytes) << "bytes size for k=" << kmer.get_k() << " is " << kbytes.second << " v " << nbytes;
  string revseq = slowrevcomp(kmer.to_string());
  EXPECT_STREQ(rev.to_string().c_str(), revseq.c_str())
      << "Revcomp to " << kmer.to_string() << " is " << revseq << " but got " << rev.to_string();

  // rev.revcomp() == kmer
  K revrev = rev.revcomp();
  EXPECT_EQ(revrev.hash(), kmer.hash()) << "revrev hash equals kmer " << revrev.hash() << " v " << kmer.hash();
  EXPECT_STREQ(revrev.to_string().c_str(), seq.c_str()) << "revrev seq equals kmer " << revrev.to_string() << " v " << seq;
  EXPECT_EQ(memcmp(revrev.get_bytes().first, kbytes.first, nbytes), 0);
  revseq = slowrevcomp(rev.to_string());
  EXPECT_STREQ(revrev.to_string().c_str(), revseq.c_str())
      << "Revcomp to " << rev.to_string() << " is " << revseq << " but got " << revrev.to_string();

  revrev.set_zeros();
  EXPECT_EQ(revrev.hash(), kmer.hash()) << "revrev hash equals kmer " << revrev.hash() << " v " << kmer.hash();
  EXPECT_STREQ(revrev.to_string().c_str(), seq.c_str()) << "revrev seq equals kmer " << revrev.to_string() << " v " << seq;
  EXPECT_EQ(memcmp(revrev.get_bytes().first, kbytes.first, nbytes), 0);
  EXPECT_STREQ(revrev.to_string().c_str(), revseq.c_str())
      << "Revcomp to " << rev.to_string() << " is " << revseq << " but got " << revrev.to_string();

  // seq construct
  {
    // make a dirty one in the stack
    K tmp(rev);
    EXPECT_STREQ(tmp.to_string().c_str(), rev.to_string().c_str());
    EXPECT_EQ(tmp.hash(), rev.hash()) << "rev hashes equal";
    tmp.set_zeros();
    EXPECT_STREQ(tmp.to_string().c_str(), rev.to_string().c_str());
    EXPECT_EQ(tmp.hash(), rev.hash()) << "rev hashes equal";
  }
  // copy construct
  K cpy(kmer);
  EXPECT_EQ(cpy.hash(), kmer.hash()) << "Copy hashes equal";
  EXPECT_STREQ(cpy.to_string().c_str(), seq.c_str());
  auto cpybytes = cpy.get_bytes();
  ASSERT_EQ(kbytes.second, cpybytes.second);
  EXPECT_EQ(memcmp(kbytes.first, cpybytes.first, nbytes), 0);
  cpy.set_zeros();
  EXPECT_EQ(cpy.hash(), kmer.hash()) << "Copy hashes equal";
  EXPECT_STREQ(cpy.to_string().c_str(), seq.c_str());
  EXPECT_EQ(cpybytes, cpy.get_bytes());
  EXPECT_EQ(memcmp(kbytes.first, cpybytes.first, nbytes), 0);

  K revcpy(rev);
  EXPECT_EQ(revcpy.hash(), rev.hash()) << "Copy hashes equal";
  EXPECT_STREQ(revcpy.to_string().c_str(), rev.to_string().c_str());
  auto revcpybytes = revcpy.get_bytes();
  ASSERT_EQ(memcmp(revcpybytes.first, rev.get_bytes().first, nbytes), 0);

  {
    // make a dirty one in the stack
    K tmp(rev);
    EXPECT_STREQ(tmp.to_string().c_str(), rev.to_string().c_str());
    EXPECT_EQ(tmp.hash(), rev.hash()) << "rev hashes equal";
    tmp.set_zeros();
    EXPECT_STREQ(tmp.to_string().c_str(), rev.to_string().c_str());
    EXPECT_EQ(tmp.hash(), rev.hash()) << "rev hashes equal";
  }
  // use the same stack
  K cpy2(kmer.to_string().c_str());
  EXPECT_EQ(kmer.get_k(), cpy2.get_k());
  EXPECT_EQ(kmer.get_N_LONGS(), cpy2.get_N_LONGS());
  EXPECT_EQ(cpy2.hash(), kmer.hash()) << "Copy hashes equal " << cpy2.to_string() << " v " << kmer.to_string() << " "
                                      << cpy2.to_hex() << " v " << kmer.to_hex();
  EXPECT_EQ(memcmp(kmer.get_bytes().first, cpy2.get_bytes().first, nbytes), 0);
  EXPECT_STREQ(cpy2.to_string().c_str(), seq.c_str());
  auto cpy2bytes = cpy2.get_bytes();
  ASSERT_EQ(kbytes.second, cpy2bytes.second);
  EXPECT_EQ(memcmp(kbytes.first, cpy2bytes.first, nbytes), 0);

  auto hash = cpy2.hash();
  seq = cpy2.to_string();
  EXPECT_EQ(hash, cpy2.hash());
  EXPECT_STREQ(seq.c_str(), cpy2.to_string().c_str());
  cpy2.set_zeros();
  EXPECT_EQ(hash, cpy2.hash());
  EXPECT_STREQ(seq.c_str(), cpy2.to_string().c_str());

  hash = rev.hash();
  seq = rev.to_string();
  EXPECT_EQ(hash, rev.hash());
  EXPECT_STREQ(seq.c_str(), rev.to_string().c_str());
  rev.set_zeros();
  EXPECT_EQ(hash, rev.hash());
  EXPECT_STREQ(seq.c_str(), rev.to_string().c_str());
}

template <int MAX_K>
void test_kmer(int klen) {
  typedef Kmer<MAX_K> K;
  K::set_k(klen);
  vector<string> temps;

  string rread(RandomRead);
  for (int i = 0; i < rread.size() - klen + 1; i++) {
    temps.push_back(rread.substr(i, klen));
  }
  for (int i = 0; i < 10; i++) {
    temps.push_back(string(ACGTs).substr(i, klen));
    temps.push_back(string(TCGAs).substr(i, klen));
    temps.push_back(string(CAGTs).substr(i, klen));
  }

  temps.push_back(string(As).substr(0, klen));
  temps.push_back(string(Cs).substr(0, klen));
  temps.push_back(string(Gs).substr(0, klen));
  temps.push_back(string(Ts).substr(0, klen));
  temps.push_back(string(ACGTs).substr(0, klen));
  temps.push_back(string(TCGAs).substr(0, klen));
  temps.push_back(string(CAGTs).substr(0, klen));

  map<longs_t, string> hashes;
  map<string, longs_t> kmers;
  for (auto &seq : temps) {
    K kmer(seq.c_str());
    check_kmer(kmer, seq);

    EXPECT_STREQ(seq.c_str(), kmer.to_string().c_str()) << seq << " == " << kmer.to_string();
    auto hash = kmer.hash();
    kmers.insert({seq, hash});
    auto find = hashes.find(hash);
    if (find == hashes.end()) {
      EXPECT_TRUE(hashes.find(hash) == hashes.end()) << "All are unique hashes hash=" << hash << " kmer=" << kmer.to_string();
      auto it = hashes.insert({hash, kmer.to_string()});
      EXPECT_TRUE(it.second) << "hash insert succeded";
      EXPECT_TRUE(it.first->first == hash) << "hashes equal";
      EXPECT_STREQ(it.first->second.c_str(), seq.c_str()) << "sequences equal";
    } else {
      EXPECT_FALSE(hashes.find(hash) == hashes.end()) << "Duplicates are not unique hashes";
      EXPECT_STREQ(hashes.find(hash)->second.c_str(), seq.c_str())
          << "Dupliate hash sequences match " << hashes.find(hash)->second << " vs " << kmer.to_string();
    }

    // revcomp
    K rev = kmer.revcomp();
    check_kmer(rev, rev.to_string());

    EXPECT_STREQ(seq.c_str(), kmer.to_string().c_str()) << "revcomp does not modify kmer";
    hash = rev.hash();
    if (rev.to_string().compare(kmer.to_string()) == 0) {
      // a panindrome
      EXPECT_EQ(hash, kmer.hash()) << "Hash for equal palindromes are equal " << rev.to_string() << " " << hash << " v "
                                   << kmer.to_string() << " " << kmer.hash();
    } else {
      EXPECT_NE(hash, kmer.hash());
    }
    kmers.insert({rev.to_string(), hash});
    find = hashes.find(hash);
    if (find == hashes.end()) {
      EXPECT_TRUE(hashes.find(hash) == hashes.end()) << "All are unique hashes hash=" << hash << " rev=" << rev.to_string();
      auto it = hashes.insert({hash, rev.to_string()});
      EXPECT_TRUE(it.second) << "hash insert succeded";
      EXPECT_TRUE(it.first->first == hash) << "hashes equal";
      EXPECT_STREQ(it.first->second.c_str(), rev.to_string().c_str()) << "sequences equal";
    } else {
      EXPECT_FALSE(hashes.find(hash) == hashes.end()) << "Duplicates are not unique hashes";
      EXPECT_STREQ(hashes.find(hash)->second.c_str(), rev.to_string().c_str())
          << "Dupliate hash sequences match " << hashes.find(hash)->second << " vs " << rev.to_string();
    }
  }
};

template <int MAX_K>
void test_kmer_minimizers(int kmer_len) {
  int m_len = 15;
  if (kmer_len < 17) return;
  string seq("AACTGACCAGACGGGGAGGATGCCATGCTGTTGAATTCTCCCCTTTATTAAGTAAGGAAGTCCGGTGATCCAGAATATTCTGCGGAGTTTTCAAATTTATGTTTTTAATTGATCC"
             "CCTGACTTGTAAAGGGAATAGTTCCCTAAAATTAC");
  Kmer<MAX_K>::set_k(kmer_len);
  vector<Kmer<MAX_K>> kmers;
  Kmer<MAX_K>::get_kmers(kmer_len, seq, kmers);
  string prev_minimizer = "";
  string prev_minz_opt = "";
  int i = 0;
  for (auto &kmer : kmers) {
    auto minimizer = kmer.get_minimizer_slow(m_len);
    if (prev_minimizer == "") {
      prev_minimizer = minimizer;
    } else {
      if (prev_minimizer != minimizer) prev_minimizer = minimizer;
    }
    auto opt_minimizer = Kmer<MAX_K>::mer_to_string(kmer.get_minimizer(m_len), m_len);
    if (prev_minz_opt == "") {
      prev_minz_opt = opt_minimizer;
    } else {
      if (prev_minz_opt != opt_minimizer) prev_minz_opt = opt_minimizer;
    }
    auto minz = kmer.get_minimizer(m_len);
    auto minz_rc = Kmer<MAX_K>::revcomp_minimizer(minz, m_len);
    auto minz_rc_back = Kmer<MAX_K>::revcomp_minimizer(minz_rc, m_len);
    EXPECT_EQ(minz, minz_rc_back) << "Revcomp of minimizers should be equal " << minz << " " << minz_rc_back;
    // if (!i && kmer_len == 55) {
    //  std::cout << Kmer<MAX_K>::mer_to_string(minz, m_len) << std::endl << Kmer<MAX_K>::mer_to_string(minz_rc, m_len) <<
    //  std::endl;
    //}
    EXPECT_EQ(minimizer, opt_minimizer) << "Minimizers are equal for slow and opt " << minimizer << " " << opt_minimizer;
    i++;
  }
}

/*
template <int MAX_K>
void test_minimizer_performance(int kmer_len) {
  int m_len = 15;
  if (kmer_len < 17) return;
  string seq("AACTGACCAGACGGGGAGGATGCCATGCTGTTGAATTCTCCCCTTTATTAAGTAAGGAAGTCCGGTGATCCAGAATATTCTGCGGAGTTTTCAAATTTATGTTTTTAATTGATCC"
             "CCTGACTTGTAAAGGGAATAGTTCCCTAAAATTAC");
  Kmer<MAX_K>::set_k(kmer_len);
  vector<Kmer<MAX_K>> kmers;
  Kmer<MAX_K>::get_kmers(kmer_len, seq, kmers);
  for (int i = 0; i < 50000; i++) {
    for (auto &kmer : kmers) {
      kmer.minimizer_hash(m_len);
    }
  }
}

TEST(MHMTest, minimizer_performance) {
  test_minimizer_performance<32>(21);
  test_minimizer_performance<64>(55);
  test_minimizer_performance<96>(77);
}
*/

TEST(MHMTest, kmer) {
  // arrange
  // act
  // assert
  for (int i = 1; i < 96; i++) {
    if (i <= 32) {
      test_kmer<32>(i);
      test_get_kmers<32>(i);
      test_kmer_minimizers<32>(i);
    }
    if (i <= 64) {
      test_kmer<64>(i);
      test_get_kmers<64>(i);
      test_kmer_minimizers<64>(i);
    }
    if (i <= 96) {
      test_kmer<96>(i);
      test_get_kmers<96>(i);
      test_kmer_minimizers<96>(i);
    }
  }
}