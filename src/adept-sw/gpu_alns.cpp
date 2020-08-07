#include "gpu_alns.hpp"
#include"utils_gpu.hpp"
gpu_alignments::gpu_alignments(int max_alignments){
    cudaErrchk(cudaMalloc(&offset_query_gpu, (max_alignments) * sizeof(int)));
    cudaErrchk(cudaMalloc(&offset_ref_gpu, (max_alignments) * sizeof(int)));
    cudaErrchk(cudaMalloc(&ref_start_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&ref_end_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&query_start_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&query_end_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&scores_gpu, (max_alignments) * sizeof(short)));
}

gpu_alignments::~gpu_alignments(){
    cudaErrchk(cudaFree(offset_query_gpu)); offset_query_gpu = nullptr;
    cudaErrchk(cudaFree(offset_ref_gpu));   offset_ref_gpu = nullptr;
    cudaErrchk(cudaFree(ref_start_gpu));    ref_start_gpu = nullptr;
    cudaErrchk(cudaFree(ref_end_gpu));      ref_end_gpu = nullptr;
    cudaErrchk(cudaFree(query_start_gpu));  query_start_gpu = nullptr;
    cudaErrchk(cudaFree(query_end_gpu));    query_end_gpu = nullptr;
    cudaErrchk(cudaFree(scores_gpu));       scores_gpu = nullptr;
}