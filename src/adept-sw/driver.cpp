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

static int device_count;
static int my_gpu_id;
static cudaStream_t streams_cuda[NSTREAMS];
static unsigned* offsetA_h;
static unsigned* offsetB_h;
static char *strA_d, *strB_d;
static char* strA;
static char* strB;
static cudaEvent_t event;

void gpu_bsw_driver::init(gpu_bsw_driver::alignment_results *alignments, int max_alignments, int my_upcxx_rank, int totRanks) {
  cudaMallocHost(&(alignments->ref_begin), sizeof(short)*max_alignments);
  cudaMallocHost(&(alignments->ref_end), sizeof(short)*max_alignments);
  cudaMallocHost(&(alignments->query_begin), sizeof(short)*max_alignments);
  cudaMallocHost(&(alignments->query_end), sizeof(short)*max_alignments);
  cudaMallocHost(&(alignments->top_scores), sizeof(short)*max_alignments);
  device_count = get_device_count(totRanks);
  my_gpu_id = my_upcxx_rank % device_count;  
  cudaSetDevice(my_gpu_id);
  for (int stm = 0; stm < NSTREAMS; stm++) {
    cudaStreamCreate(&streams_cuda[stm]);
  }
  cudaMallocHost(&offsetA_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&offsetB_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE);

  // FIXME: hack for max contig and read size
  cudaErrchk(cudaMalloc(&strA_d, 300 * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));
  cudaErrchk(cudaMalloc(&strB_d, 300 * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));

  cudaMallocHost(&strA, sizeof(char) * 300 * KLIGN_GPU_BLOCK_SIZE);
  cudaMallocHost(&strB, sizeof(char) * 300 * KLIGN_GPU_BLOCK_SIZE);
}

void gpu_bsw_driver::fini(gpu_bsw_driver::alignment_results *alignments) {
  cudaErrchk(cudaFreeHost(alignments->ref_begin));
  cudaErrchk(cudaFreeHost(alignments->ref_end));
  cudaErrchk(cudaFreeHost(alignments->query_begin));
  cudaErrchk(cudaFreeHost(alignments->query_end));
  cudaErrchk(cudaFreeHost(alignments->top_scores));

  cudaErrchk(cudaFree(strA_d));
  cudaErrchk(cudaFree(strB_d));
  cudaFreeHost(offsetA_h);
  cudaFreeHost(offsetB_h);
  cudaFreeHost(strA);
  cudaFreeHost(strB);
  for (int i = 0; i < NSTREAMS; i++) cudaStreamDestroy(streams_cuda[i]);
}

bool gpu_bsw_driver::kernel_is_done() {
//  return (cudaStreamQuery(streams_cuda[0]) == cudaSuccess && cudaStreamQuery(streams_cuda[1]) == cudaSuccess);
  return (cudaEventQuery(event) == cudaSuccess);
}

void gpu_bsw_driver::kernel_driver_dna(std::vector<std::string> reads, std::vector<std::string> contigs, unsigned maxReadSize,
                                       unsigned maxContigSize, gpu_bsw_driver::alignment_results* alignments, short scores[4],
                                       long long int maxMemAvail) {
  short matchScore = scores[0], misMatchScore = scores[1], startGap = scores[2], extendGap = scores[3];
  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length

  static gpu_alignments gpu_data(KLIGN_GPU_BLOCK_SIZE);  // gpu mallocs

  short* alAbeg = alignments->ref_begin;
  short* alBbeg = alignments->query_begin;
  short* alAend = alignments->ref_end;
  short* alBend = alignments->query_end;;  // memory on CPU for copying the results
  short* top_scores_cpu = alignments->top_scores;

  int blocksLaunched = 0;
  std::vector<std::string>::const_iterator beginAVec;
  std::vector<std::string>::const_iterator endAVec;
  std::vector<std::string>::const_iterator beginBVec;
  std::vector<std::string>::const_iterator endBVec;
  beginAVec = contigs.begin();
  endAVec = contigs.begin() + totalAlignments;
  beginBVec = reads.begin();
  endBVec = reads.begin() + totalAlignments;
  blocksLaunched = totalAlignments;

  std::vector<std::string> sequencesA(beginAVec, endAVec);
  std::vector<std::string> sequencesB(beginBVec, endBVec);
  unsigned running_sum = 0;
  int sequences_per_stream = (blocksLaunched) / NSTREAMS;
  int sequences_stream_leftover = (blocksLaunched) % NSTREAMS;
  unsigned half_length_A = 0;
  unsigned half_length_B = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    running_sum += sequencesA[i].size();
    offsetA_h[i] = running_sum;  // sequencesA[i].size();
    if (i == sequences_per_stream - 1) {
      half_length_A = running_sum;
      running_sum = 0;
    }
  }
  unsigned totalLengthA = half_length_A + offsetA_h[sequencesA.size() - 1];

  running_sum = 0;
  for (int i = 0; i < (int)sequencesB.size(); i++) {
    running_sum += sequencesB[i].size();
    offsetB_h[i] = running_sum;  // sequencesB[i].size();
    if (i == sequences_per_stream - 1) {
      half_length_B = running_sum;
      running_sum = 0;
    }
  }
  unsigned totalLengthB = half_length_B + offsetB_h[sequencesB.size() - 1];

  unsigned offsetSumA = 0;
  unsigned offsetSumB = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    char* seqptrA = strA + offsetSumA;
    memcpy(seqptrA, sequencesA[i].c_str(), sequencesA[i].size());
    char* seqptrB = strB + offsetSumB;
    memcpy(seqptrB, sequencesB[i].c_str(), sequencesB[i].size());
    offsetSumA += sequencesA[i].size();
    offsetSumB += sequencesB[i].size();
  }

  asynch_mem_copies_htd(&gpu_data, offsetA_h, offsetB_h, strA, strA_d, strB, strB_d, half_length_A, half_length_B, totalLengthA,
                        totalLengthB, sequences_per_stream, sequences_stream_leftover, streams_cuda);
  unsigned minSize = (maxReadSize < maxContigSize) ? maxReadSize : maxContigSize;
  unsigned totShmem = 3 * (minSize + 1) * sizeof(short);
  unsigned alignmentPad = 4 + (4 - totShmem % 4);
  size_t ShmemBytes = totShmem + alignmentPad;
  if (ShmemBytes > 48000)
    cudaFuncSetAttribute(gpu_bsw::sequence_dna_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, ShmemBytes);

  gpu_bsw::sequence_dna_kernel<<<sequences_per_stream, minSize, ShmemBytes, streams_cuda[0]>>>(
    strA_d, strB_d, gpu_data.offset_ref_gpu, gpu_data.offset_query_gpu, gpu_data.ref_start_gpu, gpu_data.ref_end_gpu,
    gpu_data.query_start_gpu, gpu_data.query_end_gpu, gpu_data.scores_gpu, matchScore, misMatchScore, startGap, extendGap);

  gpu_bsw::sequence_dna_kernel<<<sequences_per_stream + sequences_stream_leftover, minSize, ShmemBytes, streams_cuda[1]>>>(
    strA_d + half_length_A, strB_d + half_length_B, gpu_data.offset_ref_gpu + sequences_per_stream,
    gpu_data.offset_query_gpu + sequences_per_stream, gpu_data.ref_start_gpu + sequences_per_stream,
    gpu_data.ref_end_gpu + sequences_per_stream, gpu_data.query_start_gpu + sequences_per_stream,
    gpu_data.query_end_gpu + sequences_per_stream, gpu_data.scores_gpu + sequences_per_stream, matchScore, misMatchScore,
    startGap, extendGap);

  // copyin back end index so that we can find new min
  asynch_mem_copies_dth_mid(&gpu_data, alAend, alBend, sequences_per_stream, sequences_stream_leftover, streams_cuda);

  cudaStreamSynchronize(streams_cuda[0]);
  cudaStreamSynchronize(streams_cuda[1]);

  int newMin = get_new_min_length(alAend, alBend, blocksLaunched);  // find the new largest of smaller lengths

  cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
  gpu_bsw::sequence_dna_reverse<<<sequences_per_stream, newMin, ShmemBytes, streams_cuda[0]>>>(
    strA_d, strB_d, gpu_data.offset_ref_gpu, gpu_data.offset_query_gpu, gpu_data.ref_start_gpu, gpu_data.ref_end_gpu,
    gpu_data.query_start_gpu, gpu_data.query_end_gpu, gpu_data.scores_gpu, matchScore, misMatchScore, startGap, extendGap);

  gpu_bsw::sequence_dna_reverse<<<sequences_per_stream + sequences_stream_leftover, newMin, ShmemBytes, streams_cuda[1]>>>(
    strA_d + half_length_A, strB_d + half_length_B, gpu_data.offset_ref_gpu + sequences_per_stream,
    gpu_data.offset_query_gpu + sequences_per_stream, gpu_data.ref_start_gpu + sequences_per_stream,
    gpu_data.ref_end_gpu + sequences_per_stream, gpu_data.query_start_gpu + sequences_per_stream,
    gpu_data.query_end_gpu + sequences_per_stream, gpu_data.scores_gpu + sequences_per_stream, matchScore, misMatchScore,
    startGap, extendGap);

  asynch_mem_copies_dth(&gpu_data, alAbeg, alBbeg, top_scores_cpu, sequences_per_stream, sequences_stream_leftover,
                        streams_cuda);
  cudaEventRecord(event);
  
  alAbeg += totalAlignments;
  alBbeg += totalAlignments;
  alAend += totalAlignments;
  alBend += totalAlignments;
  top_scores_cpu += totalAlignments;

//  cudaStreamSynchronize(streams_cuda[0]);
//  cudaStreamSynchronize(streams_cuda[1]);
}

