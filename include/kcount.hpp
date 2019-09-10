#ifndef _KCOUNT_HPP
#define _KCOUNT_HPP

#include "options.hpp"
#include "kmer_dht.hpp"
#include "contigs.hpp"

uint64_t estimate_cardinality(unsigned kmer_len, vector<string> reads_fname_list);
void count_kmers(unsigned kmer_len, int qual_offset, vector<string> reads_fname_list, dist_object<KmerDHT> &kmer_dht, PASS_TYPE pass_type);
void count_ctg_kmers(unsigned kmer_len, Contigs &ctgs, dist_object<KmerDHT> &kmer_dht);
void add_ctg_kmers(unsigned kmer_len, Contigs &ctgs, dist_object<KmerDHT> &kmer_dht, bool use_bloom);

#endif
