#ifndef _DBJG_TRAVERSAL_H
#define _DBJG_TRAVERSAL_H

#include "options.hpp"

void traverse_debruijn_graph(shared_ptr<Options> options, dist_object<KmerDHT> &kmer_dht);

#endif
