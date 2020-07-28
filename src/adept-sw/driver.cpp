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


void gpu_bsw_driver::free_alignments(gpu_bsw_driver::alignment_results *alignments){
  cudaErrchk(cudaFreeHost(alignments->ref_begin));
  cudaErrchk(cudaFreeHost(alignments->ref_end));
  cudaErrchk(cudaFreeHost(alignments->query_begin));
  cudaErrchk(cudaFreeHost(alignments->query_end));
  cudaErrchk(cudaFreeHost(alignments->top_scores));
}

void gpu_bsw_driver::kernel_driver_dna(std::vector<std::string> reads, std::vector<std::string> contigs, unsigned maxReadSize,
                                       unsigned maxContigSize, gpu_bsw_driver::alignment_results* alignments, short scores[4],
                                       long long int maxMemAvail, unsigned my_upcxx_rank, unsigned totRanks) {
  short matchScore = scores[0], misMatchScore = scores[1], startGap = scores[2], extendGap = scores[3];
  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length

  initialize_alignments(alignments, totalAlignments);  // pinned memory allocation
  auto start = NOW;

  float total_time_cpu = 0;
  int device_count = get_device_count(totRanks);
  int my_gpu_id = my_upcxx_rank % device_count;  
  cudaSetDevice(my_gpu_id);

  cudaStream_t streams_cuda[NSTREAMS];
  for (int stm = 0; stm < NSTREAMS; stm++) {
    cudaStreamCreate(&streams_cuda[stm]);
  }

  gpu_alignments gpu_data(totalAlignments);  // gpu mallocs

  short* alAbeg = alignments->ref_begin;
  short* alBbeg = alignments->query_begin;
  short* alAend = alignments->ref_end;
  short* alBend = alignments->query_end;;  // memory on CPU for copying the results
  short* top_scores_cpu = alignments->top_scores;
  unsigned* offsetA_h;
  cudaMallocHost(&offsetA_h, sizeof(int) * totalAlignments);
  unsigned* offsetB_h;
  cudaMallocHost(&offsetB_h, sizeof(int) * totalAlignments);

  char *strA_d, *strB_d;
  cudaErrchk(cudaMalloc(&strA_d, maxContigSize * totalAlignments * sizeof(char)));
  cudaErrchk(cudaMalloc(&strB_d, maxReadSize * totalAlignments * sizeof(char)));

  char* strA;
  cudaMallocHost(&strA, sizeof(char) * maxContigSize * totalAlignments);
  char* strB;
  cudaMallocHost(&strB, sizeof(char) * maxReadSize * totalAlignments);

  float total_packing = 0;

  auto start2 = NOW;
  auto packing_start = NOW;
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

  auto start_cpu = NOW;

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

  auto end_cpu = NOW;
  std::chrono::duration<double> cpu_dur = end_cpu - start_cpu;

  total_time_cpu += cpu_dur.count();
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

  auto packing_end = NOW;
  std::chrono::duration<double> packing_dur = packing_end - packing_start;

  total_packing += packing_dur.count();

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

  auto sec_cpu_start = NOW;
  int newMin = get_new_min_length(alAend, alBend, blocksLaunched);  // find the new largest of smaller lengths
  auto sec_cpu_end = NOW;
  std::chrono::duration<double> dur_sec_cpu = sec_cpu_end - sec_cpu_start;
  total_time_cpu += dur_sec_cpu.count();

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

  alAbeg += totalAlignments;
  alBbeg += totalAlignments;
  alAend += totalAlignments;
  alBend += totalAlignments;
  top_scores_cpu += totalAlignments;

  auto end1 = NOW;
  std::chrono::duration<double> diff2 = end1 - start2;
  cudaErrchk(cudaFree(strA_d));
  cudaErrchk(cudaFree(strB_d));
  cudaFreeHost(offsetA_h);
  cudaFreeHost(offsetB_h);
  cudaFreeHost(strA);
  cudaFreeHost(strB);

  for (int i = 0; i < NSTREAMS; i++) cudaStreamDestroy(streams_cuda[i]);

  auto end = NOW;
  std::chrono::duration<double> diff = end - start;
}  // end of DNA kernel

