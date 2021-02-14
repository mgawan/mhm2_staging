#include <unordered_map>
//#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>
#include <cstring>
#include <fstream>
// #include "helper.hpp"
// #include "kernel.hpp"
#include "driver.hpp"

// size_t get_device_mem(){
//     int gpus;
//     size_t free_mem, total_mem;
//     cudaGetDeviceCount(&gpus);
//     for ( int id = 0; id < gpus; id++ ) {
//         cudaSetDevice(id);
//         cudaMemGetInfo(&free_mem, &total_mem);
//         std::cout << "GPU: " << id << " has free memory (Mbytes):=" << (double)free_mem/(1024*1024) << ", out of total (Mbytes):=" << (double)total_mem/(1024*1024) << std::endl;
//     }
//     return free_mem;
// }

// struct accum_data{
//     std::vector<uint32_t> ht_sizes;
//     std::vector<uint32_t> l_reads_count;
//     std::vector<uint32_t> r_reads_count;
//     std::vector<uint32_t> ctg_sizes;
// };

// std::ofstream ofile("/global/cscratch1/sd/mgawan/local_assem_large/haswell_large/merged/test-results/test-out.dat");
// // void call_kernel(std::vector<CtgWithReads>& data_in, uint32_t max_ctg_size, uint32_t max_read_size, uint32_t max_r_count, uint32_t max_l_count, int mer_len,int max_reads_count, accum_data& sizes_outliers);

// sample cmd line: ./build/ht_loc ../locassm_data/localassm_extend_7-21.dat ./out_file <kmer_size>
int main (int argc, char* argv[]){

    std::string in_file = argv[1];
   // std::string out_file = argv[2];
    int max_mer_len = std::stoi(argv[2]);
    std::vector<CtgWithReads> data_in;
    uint32_t max_ctg_size, total_r_reads, total_l_reads, max_read_size, max_r_count, max_l_count;
    read_locassm_data(&data_in, in_file, max_ctg_size, total_r_reads, total_l_reads, max_read_size, max_r_count, max_l_count);
    timer overall_time;
    overall_time.timer_start();

    std::vector<CtgWithReads> zero_slice, mid_slice, midsup_slice, outlier_slice;
    uint32_t mid_l_max = 0, mid_r_max = 0, outlier_l_max = 0, outlier_r_max = 0, mid_max_contig_sz = 0;
    uint32_t outliers_max_contig_sz = 0, mids_tot_r_reads = 0, mids_tot_l_reads = 0, outliers_tot_r_reads = 0;
    uint32_t outliers_tot_l_reads = 0;
   // std::vector<uint32_t> ht_size_mid, ht_size_midsup, ht_size_outliers;
    accum_data sizes_mid, sizes_outliers;

    for(int i = 0; i < data_in.size(); i++){
        CtgWithReads temp_in = data_in[i];
        if(temp_in.max_reads == 0){
            zero_slice.push_back(temp_in);
        }else if(temp_in.max_reads > 0 && temp_in.max_reads < 10){
            mid_slice.push_back(temp_in);
            uint32_t temp_ht_size = temp_in.max_reads * max_read_size;
            sizes_mid.ht_sizes.push_back(temp_ht_size);
            sizes_mid.ctg_sizes.push_back(temp_in.seq.size());
            sizes_mid.l_reads_count.push_back(temp_in.reads_left.size());
            sizes_mid.r_reads_count.push_back(temp_in.reads_right.size());
            mids_tot_r_reads += temp_in.reads_right.size();
            mids_tot_l_reads += temp_in.reads_left.size();
            if(mid_l_max < temp_in.reads_left.size())
                mid_l_max = temp_in.reads_left.size();
            if(mid_r_max < temp_in.reads_right.size())
                mid_r_max = temp_in.reads_right.size();
            if(mid_max_contig_sz < temp_in.seq.size())
                mid_max_contig_sz = temp_in.seq.size();
        }
        else{
            outlier_slice.push_back(temp_in);
            uint32_t temp_ht_size = temp_in.max_reads * max_read_size;
            sizes_outliers.ht_sizes.push_back(temp_ht_size);
            sizes_outliers.ctg_sizes.push_back(temp_in.seq.size());
            sizes_outliers.l_reads_count.push_back(temp_in.reads_left.size());
            sizes_outliers.r_reads_count.push_back(temp_in.reads_right.size());
            outliers_tot_r_reads += temp_in.reads_right.size();
            outliers_tot_l_reads += temp_in.reads_left.size();
            if(outlier_l_max < temp_in.reads_left.size())
                outlier_l_max = temp_in.reads_left.size();
            if(outlier_r_max < temp_in.reads_right.size())
                outlier_r_max = temp_in.reads_right.size();
            if(outliers_max_contig_sz < temp_in.seq.size())
                outliers_max_contig_sz = temp_in.seq.size();
        }
    }

    print_vals("zeroes, count:", zero_slice.size());
    timer file_time;
    file_time.timer_start();
    // for(int i = 0; i < zero_slice.size(); i++){
    //     ofile << zero_slice[i].cid<<" "<<zero_slice[i].seq<<std::endl;

    // }
    file_time.timer_end();
    print_vals("zeroes file write time:",file_time.get_total_time());
    zero_slice = std::vector<CtgWithReads>();
    data_in = std::vector<CtgWithReads>();
     print_vals("mids calling",  "mids count:", mid_slice.size());
    
    int max_reads_count = 10;
    local_assem_driver(mid_slice, mid_max_contig_sz, max_read_size, mid_r_max, mid_l_max, max_mer_len,max_reads_count, sizes_mid);
    print_vals("midsup calling",  "mids count:", midsup_slice.size());
    // max_reads_count = 100;
    // call_kernel(midsup_slice, midsup_max_contig_sz, midsup_tot_r_reads, midsup_tot_l_reads, max_read_size, midsup_r_max, midsup_l_max, max_mer_len,max_reads_count, ht_size_midsup);

    print_vals("outliers calling", "outliers count:", outlier_slice.size());
    // max_reads_count = 239;
  //  overall_time.timer_start();
    local_assem_driver(outlier_slice, outliers_max_contig_sz, max_read_size, outlier_r_max, outlier_l_max, max_mer_len, max_reads_count, sizes_outliers);
    //overall_time.timer_end();
    overall_time.timer_end();
    
    print_vals("Total Time including file write:", overall_time.get_total_time());

   // ofile.flush();
    //ofile.close();

     return 0;
}