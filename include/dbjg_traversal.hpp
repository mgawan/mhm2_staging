#ifndef _DBJG_TRAVERSAL_H
#define _DBJG_TRAVERSAL_H

#include "options.hpp"
#include "contigs.hpp"

void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT> &kmer_dht, Contigs &my_uutigs);

#endif
