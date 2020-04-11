// cgraph 
// Steven Hofmeyr, LBNL, Aug 2018

#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"

#include "ctg_graph.hpp"
#include "utils.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "fastq.hpp"


using namespace std;
using namespace upcxx;
using namespace upcxx_utils;


void build_ctg_graph(CtgGraph *graph, int insert_avg, int insert_stddev, int kmer_len, int read_len,
                     vector<FastqReader*> &fqr_list, Contigs &ctgs, Alns &alns);
void walk_graph(CtgGraph *graph, int max_kmer_len, int kmer_len, int break_scaff_Ns, QualityLevel quality_level, 
                Contigs &ctgs);


void traverse_ctg_graph(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, int read_len, int min_ctg_print_len,
                        vector<FastqReader*> &fqr_list, int break_scaff_Ns, QualityLevel quality_level, 
                        Contigs &ctgs, Alns &alns) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  
  CtgGraph graph;
  build_ctg_graph(&graph, insert_avg, insert_stddev, kmer_len, read_len, fqr_list, ctgs, alns);
  barrier();
  barrier();
  ctgs.clear();
  string graph_fname;
  graph.print_stats(min_ctg_print_len);
  barrier();
  walk_graph(&graph, max_kmer_len, kmer_len, break_scaff_Ns, quality_level, ctgs);
  barrier();
}

