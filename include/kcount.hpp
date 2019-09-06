#ifndef _KCOUNT_HPP
#define _KCOUNT_HPP

#include "options.hpp"
#include "kmer_dht.hpp"

uint64_t estimate_cardinality(shared_ptr<Options> options);
void count_kmers(shared_ptr<Options> options, dist_object<KmerDHT> &kmer_dht, PASS_TYPE pass_type);
void add_ctg_kmers(shared_ptr<Options> options, dist_object<KmerDHT> &kmer_dht, bool use_bloom, int pass_num_mask);

#endif
