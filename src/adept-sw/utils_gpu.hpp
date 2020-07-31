#pragma once

//#include <cuda_runtime_api.h>

#include "driver.hpp"
#include "gpu_alns.hpp"
#define cudaErrchk(ans)                                                                  \
{                                                                                    \
    utils_gpu::gpuAssert((ans), __FILE__, __LINE__);                                            \
}

namespace utils_gpu {
    
void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true);

unsigned getMaxLength(std::vector<std::string> v);
void asynch_mem_copies_htd(gpu_alignments* gpu_data, unsigned* offsetA_h, unsigned* offsetB_h, char* strA, char* strA_d, char* strB, char* strB_d, unsigned half_length_A, unsigned half_length_B, unsigned totalLengthA, unsigned totalLengthB, int sequences_per_stream, int sequences_stream_leftover, cudaStream_t* streams_cuda);
int get_new_min_length(short* alAend, short* alBend, int blocksLaunched);
void asynch_mem_copies_dth_mid(gpu_alignments* gpu_data, short* alAend, short* alBend, int sequences_per_stream, int sequences_stream_leftover, cudaStream_t* streams_cuda);
void asynch_mem_copies_dth(gpu_alignments* gpu_data, short* alAbeg, short* alBbeg, short* top_scores_cpu, int sequences_per_stream, int sequences_stream_leftover, cudaStream_t* streams_cuda);

};
