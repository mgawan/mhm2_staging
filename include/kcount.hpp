#ifndef _KCOUNT_HPP
#define _KCOUNT_HPP

#include "options.hpp"
#include "kmer_dht.hpp"
#include "contigs.hpp"

uint64_t estimate_num_kmers(unsigned kmer_len, vector<string> reads_fname_list);
void analyze_kmers(unsigned int kmer_len, int qual_offset, vector<string> reads_fname_list, bool use_bloom, int min_depth_cutoff,
                   double dynamic_min_depth, Contigs &ctgs, dist_object<KmerDHT> &kmer_dht);

#endif
