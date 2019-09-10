#ifndef _DBJG_TRAVERSAL_H
#define _DBJG_TRAVERSAL_H

#include "options.hpp"

std::string traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht);

#endif
