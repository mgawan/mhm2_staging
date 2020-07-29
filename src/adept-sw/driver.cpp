#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include "driver.hpp"
#include "utils_gpu.hpp"
#include "gpu_alns.hpp"
#include "kernel.hpp"

static int get_device_count(int totRanks) {
  int deviceCount = 0;
  cudaErrchk(cudaGetDeviceCount(&deviceCount));
  if (deviceCount > totRanks) return totRanks;
  return deviceCount;
}

size_t gpu_bsw_driver::get_avail_gpu_mem_per_rank(int totRanks) {
  int ranksPerDevice = totRanks / get_device_count(totRanks);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  return (prop.totalGlobalMem * 0.8) / ranksPerDevice;
}

size_t gpu_bsw_driver::get_tot_gpu_mem() {
  cudaDeviceProp prop;
  cudaErrchk(cudaGetDeviceProperties(&prop, 0));
  return prop.totalGlobalMem;
}

int gpu_bsw_driver::get_num_node_gpus() {
  int deviceCount = 0;
  cudaErrchk(cudaGetDeviceCount(&deviceCount));
  return deviceCount;
}

struct gpu_bsw_driver::DriverState {
  int device_count;
  int my_gpu_id;
  cudaStream_t streams_cuda[NSTREAMS];
  unsigned* offsetA_h;
  unsigned* offsetB_h;
  char *strA_d, *strB_d;
  char* strA;
  char* strB;
  cudaEvent_t event;
  short matchScore, misMatchScore, startGap, extendGap;
  gpu_alignments *gpu_data;
  unsigned half_length_A = 0;
  unsigned half_length_B = 0;
};

void gpu_bsw_driver::GPUDriver::init(int upcxx_rank_me, int upcxx_rank_n, short match_score, short mismatch_score,
                                     short gap_opening_score, short gap_extending_score, int max_rlen) {
  driver_state = new DriverState();
  driver_state->matchScore = match_score;
  driver_state->misMatchScore = mismatch_score;
  driver_state->startGap = gap_opening_score;
  driver_state->extendGap = gap_extending_score;
  cudaMallocHost(&(alignments.ref_begin), sizeof(short) * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&(alignments.ref_end), sizeof(short) * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&(alignments.query_begin), sizeof(short) * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&(alignments.query_end), sizeof(short) * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&(alignments.top_scores), sizeof(short) * KLIGN_GPU_BLOCK_SIZE);
  driver_state->device_count = get_device_count(upcxx_rank_n);
  driver_state->my_gpu_id = upcxx_rank_me % driver_state->device_count;  
  cudaSetDevice(driver_state->my_gpu_id);
  for (int stm = 0; stm < NSTREAMS; stm++) {
    cudaStreamCreate(&driver_state->streams_cuda[stm]);
  }
  cudaMallocHost(&driver_state->offsetA_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&driver_state->offsetB_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE);

  // FIXME: hack for max contig and read size
  cudaErrchk(cudaMalloc(&driver_state->strA_d, max_rlen * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));
  cudaErrchk(cudaMalloc(&driver_state->strB_d, max_rlen * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));

  cudaMallocHost(&driver_state->strA, sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&driver_state->strB, sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE);
  driver_state->gpu_data = new gpu_alignments(KLIGN_GPU_BLOCK_SIZE);  // gpu mallocs
}
  
gpu_bsw_driver::GPUDriver::~GPUDriver() {
  cudaErrchk(cudaFreeHost(alignments.ref_begin));
  cudaErrchk(cudaFreeHost(alignments.ref_end));
  cudaErrchk(cudaFreeHost(alignments.query_begin));
  cudaErrchk(cudaFreeHost(alignments.query_end));
  cudaErrchk(cudaFreeHost(alignments.top_scores));

  cudaErrchk(cudaFree(driver_state->strA_d));
  cudaErrchk(cudaFree(driver_state->strB_d));
  cudaFreeHost(driver_state->offsetA_h);
  cudaFreeHost(driver_state->offsetB_h);
  cudaFreeHost(driver_state->strA);
  cudaFreeHost(driver_state->strB);
  for (int i = 0; i < NSTREAMS; i++) cudaStreamDestroy(driver_state->streams_cuda[i]);
  delete driver_state->gpu_data;
  delete driver_state;
}

bool gpu_bsw_driver::GPUDriver::kernel_is_done() {
  if (cudaEventQuery(driver_state->event) != cudaSuccess) return false;
  cudaEventDestroy(driver_state->event);
  return true;
}

void gpu_bsw_driver::GPUDriver::run_kernel_forwards(std::vector<std::string> reads, std::vector<std::string> contigs,
                                                    unsigned maxReadSize, unsigned maxContigSize) {
  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length

  short* alAend = alignments.ref_end;
  short* alBend = alignments.query_end;;  // memory on CPU for copying the results

  int blocksLaunched = totalAlignments;
  std::vector<std::string>::const_iterator beginAVec;
  std::vector<std::string>::const_iterator endAVec;
  std::vector<std::string>::const_iterator beginBVec;
  std::vector<std::string>::const_iterator endBVec;
  beginAVec = contigs.begin();
  endAVec = contigs.begin() + totalAlignments;
  beginBVec = reads.begin();
  endBVec = reads.begin() + totalAlignments;

  std::vector<std::string> sequencesA(beginAVec, endAVec);
  std::vector<std::string> sequencesB(beginBVec, endBVec);
  unsigned running_sum = 0;
  int sequences_per_stream = (blocksLaunched) / NSTREAMS;
  int sequences_stream_leftover = (blocksLaunched) % NSTREAMS;
  driver_state->half_length_A = 0;
  driver_state->half_length_B = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    running_sum += sequencesA[i].size();
    driver_state->offsetA_h[i] = running_sum;  // sequencesA[i].size();
    if (i == sequences_per_stream - 1) {
      driver_state->half_length_A = running_sum;
      running_sum = 0;
    }
  }
  unsigned totalLengthA = driver_state->half_length_A + driver_state->offsetA_h[sequencesA.size() - 1];

  running_sum = 0;
  for (int i = 0; i < (int)sequencesB.size(); i++) {
    running_sum += sequencesB[i].size();
    driver_state->offsetB_h[i] = running_sum;  // sequencesB[i].size();
    if (i == sequences_per_stream - 1) {
      driver_state->half_length_B = running_sum;
      running_sum = 0;
    }
  }
  unsigned totalLengthB = driver_state->half_length_B + driver_state->offsetB_h[sequencesB.size() - 1];

  unsigned offsetSumA = 0;
  unsigned offsetSumB = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    char* seqptrA = driver_state->strA + offsetSumA;
    memcpy(seqptrA, sequencesA[i].c_str(), sequencesA[i].size());
    char* seqptrB = driver_state->strB + offsetSumB;
    memcpy(seqptrB, sequencesB[i].c_str(), sequencesB[i].size());
    offsetSumA += sequencesA[i].size();
    offsetSumB += sequencesB[i].size();
  }

  cudaEventCreateWithFlags(&driver_state->event, cudaEventDisableTiming);
  
  asynch_mem_copies_htd(driver_state->gpu_data, driver_state->offsetA_h, driver_state->offsetB_h, driver_state->strA,
                        driver_state->strA_d, driver_state->strB, driver_state->strB_d, driver_state->half_length_A,
                        driver_state->half_length_B, totalLengthA, totalLengthB, sequences_per_stream, sequences_stream_leftover,
                        driver_state->streams_cuda);
  unsigned minSize = (maxReadSize < maxContigSize) ? maxReadSize : maxContigSize;
  unsigned totShmem = 3 * (minSize + 1) * sizeof(short);
  unsigned alignmentPad = 4 + (4 - totShmem % 4);
  size_t ShmemBytes = totShmem + alignmentPad;
  if (ShmemBytes > 48000)
    cudaFuncSetAttribute(gpu_bsw::sequence_dna_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, ShmemBytes);

  gpu_bsw::sequence_dna_kernel<<<sequences_per_stream, minSize, ShmemBytes, driver_state->streams_cuda[0]>>>(
    driver_state->strA_d, driver_state->strB_d, driver_state->gpu_data->offset_ref_gpu, driver_state->gpu_data->offset_query_gpu,
    driver_state->gpu_data->ref_start_gpu, driver_state->gpu_data->ref_end_gpu, driver_state->gpu_data->query_start_gpu,
    driver_state->gpu_data->query_end_gpu, driver_state->gpu_data->scores_gpu, driver_state->matchScore,
    driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  gpu_bsw::sequence_dna_kernel<<<sequences_per_stream + sequences_stream_leftover, minSize, ShmemBytes, driver_state->streams_cuda[1]>>>(
    driver_state->strA_d + driver_state->half_length_A, driver_state->strB_d + driver_state->half_length_B,
    driver_state->gpu_data->offset_ref_gpu + sequences_per_stream,
    driver_state->gpu_data->offset_query_gpu + sequences_per_stream, driver_state->gpu_data->ref_start_gpu + sequences_per_stream,
    driver_state->gpu_data->ref_end_gpu + sequences_per_stream, driver_state->gpu_data->query_start_gpu + sequences_per_stream,
    driver_state->gpu_data->query_end_gpu + sequences_per_stream, driver_state->gpu_data->scores_gpu + sequences_per_stream,
    driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  // copyin back end index so that we can find new min
  asynch_mem_copies_dth_mid(driver_state->gpu_data, alAend, alBend, sequences_per_stream, sequences_stream_leftover,
                            driver_state->streams_cuda);

  cudaEventRecord(driver_state->event);
//  cudaStreamSynchronize(driver_state->streams_cuda[0]);
//  cudaStreamSynchronize(driver_state->streams_cuda[1]);

}

void gpu_bsw_driver::GPUDriver::run_kernel_backwards(std::vector<std::string> reads, std::vector<std::string> contigs,
                                                     unsigned maxReadSize, unsigned maxContigSize) {
  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length
  
  short* alAbeg = alignments.ref_begin;
  short* alBbeg = alignments.query_begin;
  short* alAend = alignments.ref_end;
  short* alBend = alignments.query_end;;  // memory on CPU for copying the results
  short* top_scores_cpu = alignments.top_scores;
  int blocksLaunched = totalAlignments;
  int sequences_per_stream = (blocksLaunched) / NSTREAMS;
  int sequences_stream_leftover = (blocksLaunched) % NSTREAMS;
  unsigned minSize = (maxReadSize < maxContigSize) ? maxReadSize : maxContigSize;
  unsigned totShmem = 3 * (minSize + 1) * sizeof(short);
  unsigned alignmentPad = 4 + (4 - totShmem % 4);
  size_t ShmemBytes = totShmem + alignmentPad;

  int newMin = get_new_min_length(alAend, alBend, blocksLaunched);  // find the new largest of smaller lengths
  
  cudaEventCreateWithFlags(&driver_state->event, cudaEventDisableTiming);
  gpu_bsw::sequence_dna_reverse<<<sequences_per_stream, newMin, ShmemBytes, driver_state->streams_cuda[0]>>>(
    driver_state->strA_d, driver_state->strB_d, driver_state->gpu_data->offset_ref_gpu, driver_state->gpu_data->offset_query_gpu,
    driver_state->gpu_data->ref_start_gpu, driver_state->gpu_data->ref_end_gpu, driver_state->gpu_data->query_start_gpu,
    driver_state->gpu_data->query_end_gpu, driver_state->gpu_data->scores_gpu,
    driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);
    
  gpu_bsw::sequence_dna_reverse<<<sequences_per_stream + sequences_stream_leftover, newMin, ShmemBytes, driver_state->streams_cuda[1]>>>(
    driver_state->strA_d + driver_state->half_length_A, driver_state->strB_d + driver_state->half_length_B,
    driver_state->gpu_data->offset_ref_gpu + sequences_per_stream,
    driver_state->gpu_data->offset_query_gpu + sequences_per_stream, driver_state->gpu_data->ref_start_gpu + sequences_per_stream,
    driver_state->gpu_data->ref_end_gpu + sequences_per_stream, driver_state->gpu_data->query_start_gpu + sequences_per_stream,
    driver_state->gpu_data->query_end_gpu + sequences_per_stream, driver_state->gpu_data->scores_gpu + sequences_per_stream,
    driver_state->matchScore, driver_state->misMatchScore,  driver_state->startGap, driver_state->extendGap);

  asynch_mem_copies_dth(driver_state->gpu_data, alAbeg, alBbeg, top_scores_cpu, sequences_per_stream, sequences_stream_leftover,
                        driver_state->streams_cuda);
  cudaEventRecord(driver_state->event);
  
  alAbeg += totalAlignments;
  alBbeg += totalAlignments;
  alAend += totalAlignments;
  alBend += totalAlignments;
  top_scores_cpu += totalAlignments;
}

