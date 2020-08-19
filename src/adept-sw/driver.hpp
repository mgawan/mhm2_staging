#pragma once

#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <thread>

#define NSTREAMS 2

#ifndef KLIGN_GPU_BLOCK_SIZE
#define KLIGN_GPU_BLOCK_SIZE 20000
#endif

namespace adept_sw {

// for storing the alignment results
struct AlignmentResults {
  short *ref_begin = nullptr;
  short *query_begin = nullptr;
  short *ref_end = nullptr;
  short *query_end = nullptr;
  short *top_scores = nullptr;
};

size_t get_tot_gpu_mem();
size_t get_avail_gpu_mem_per_rank(int totRanks, int numDevices = 0);
std::string get_device_name(int device_id);
int get_num_node_gpus();

// The first call to cudaMallocHost can take several seconds of real time but no cpu time
// so start it asap in a new thread
std::thread *initialize_gpu();
std::thread *initialize_gpu(double &time_to_initialize);

struct DriverState;

class GPUDriver {
  DriverState *driver_state = nullptr;
  AlignmentResults alignments;

 public:
  ~GPUDriver();
  
  // returns the time to execute
  double init(int upcxx_rank_me, int upcxx_rank_n, short match_score, short mismatch_score, short gap_opening_score,
            short gap_extending_score, int rlen_limit);
  void run_kernel_forwards(std::vector<std::string> &reads, std::vector<std::string> &contigs, unsigned maxReadSize,
                           unsigned maxContigSize);
  void run_kernel_backwards(std::vector<std::string> &reads, std::vector<std::string> &contigs, unsigned maxReadSize,
                            unsigned maxContigSize);
  bool kernel_is_done();
  void kernel_block();

  AlignmentResults &get_aln_results() { return alignments; }
};

}  // namespace adept_sw
