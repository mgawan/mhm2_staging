#pragma once

#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <omp.h>

#define NSTREAMS 2

#define NOW std::chrono::high_resolution_clock::now()

namespace gpu_bsw_driver{

// for storing the alignment results
struct alignment_results{
  short* ref_begin;
  short* query_begin;
  short* ref_end;
  short* query_end;
  short* top_scores;
};

void free_alignments(gpu_bsw_driver::alignment_results *alignments);

void
kernel_driver_dna(std::vector<std::string> reads, std::vector<std::string> contigs, unsigned maxReadSize, unsigned maxContigSize, alignment_results *alignments, short scores[4], long long int maxMemAvail, unsigned my_upc_rank, unsigned totRanks);


void
kernel_driver_aa(std::vector<std::string> reads, std::vector<std::string> contigs, alignment_results *alignments, short scoring_matrix[], short openGap, short extendGap);

void
verificationTest(std::string rstFile, short* g_alAbeg, short* g_alBbeg, short* g_alAend,
                 short* g_alBend);

size_t get_tot_gpu_mem();
size_t get_avail_gpu_mem_per_rank(int totRanks);
int get_num_node_gpus();

}
