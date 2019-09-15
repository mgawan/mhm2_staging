// cgraph 
// Steven Hofmeyr, LBNL, Aug 2018

#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "ctg_graph.hpp"
#include "utils.hpp"
#include "progressbar.hpp"
#include "contigs.hpp"
#include "alignments.hpp"


using namespace std;
using namespace upcxx;


void build_ctg_graph(CtgGraph *graph, int max_kmer_len, int kmer_len, vector<string> &reads_fname_list, Contigs *ctgs, Alns &alns);
void walk_graph(CtgGraph *graph, int max_kmer_len, int kmer_len, int min_ctg_len, bool break_scaffolds, Contigs *ctgs);


void run_scaffolding(int max_kmer_len, int kmer_len, vector<string> &reads_fname_list, Contigs *ctgs, Alns &alns) {
  Timer timer(__func__);
  CtgGraph graph;
  build_ctg_graph(&graph, max_kmer_len, kmer_len, reads_fname_list, ctgs, alns);
  barrier();
  ctgs->clear();
  string graph_fname;
  graph.print_stats();
  barrier();
  int min_ctg_len = 300;
  bool break_scaffolds = false;
  walk_graph(&graph, max_kmer_len, kmer_len, min_ctg_len, break_scaffolds, ctgs);
  barrier();
}

