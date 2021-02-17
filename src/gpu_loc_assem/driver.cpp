#include "driver.hpp"
#include "kernel.hpp"
#include <numeric>

#define CUDA_CHECK(ans)                                                                  \
    {                                                                                    \
        gpuAssert((ans), __FILE__, __LINE__);                                            \
    }
inline void
gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
    if(code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if(abort)
            exit(code);
    }
}


void locassm_driver::ctg_bucket::clear(){
    sizes_vec.ht_sizes.clear();
    sizes_vec.ht_sizes.shrink_to_fit();

    sizes_vec.l_reads_count.clear();
    sizes_vec.l_reads_count.shrink_to_fit();

    sizes_vec.r_reads_count.clear();
    sizes_vec.r_reads_count.shrink_to_fit();

    sizes_vec.ctg_sizes.clear();
    sizes_vec.ctg_sizes.shrink_to_fit();

    ctg_vec.clear();
    ctg_vec.shrink_to_fit();
}

size_t locassm_driver::get_device_mem(int ranks_per_gpu, int gpu_id){
    size_t free_mem, total_mem;
    cudaSetDevice(gpu_id);
    cudaMemGetInfo(&free_mem, &total_mem);

    return free_mem/ranks_per_gpu;
}

int locassm_driver::get_gpu_per_node() {
  int deviceCount = 0;
  auto res = cudaGetDeviceCount(&deviceCount);
  if (res != cudaSuccess) return 0;
  return deviceCount;
}

void locassm_driver::local_assem_driver(std::vector<loc_assem_helper::CtgWithReads>& data_in, uint32_t max_ctg_size, uint32_t max_read_size, uint32_t max_r_count, uint32_t max_l_count, int mer_len, int max_kmer_len, accum_data& sizes_vecs, int walk_len_limit, int qual_offset, int ranks, int my_rank)
{

    int total_gpus_avail = get_gpu_per_node();
    int my_gpu_id = my_rank % total_gpus_avail;
   // print_vals("rank:",my_rank, "gpu:", my_gpu_id);
    CUDA_CHECK(cudaSetDevice(my_gpu_id));
    int max_mer_len = max_kmer_len;//mer_len;// max_mer_len needs to come from macro (121) and mer_len is the mer_len for current go

    unsigned tot_extensions = data_in.size();
    uint32_t max_read_count = max_r_count>max_l_count ? max_r_count : max_l_count;
    int max_walk_len = walk_len_limit;
    uint32_t ht_tot_size = std::accumulate(sizes_vecs.ht_sizes.begin(), sizes_vecs.ht_sizes.end(), 0);
    uint32_t total_r_reads = std::accumulate(sizes_vecs.r_reads_count.begin(), sizes_vecs.r_reads_count.end(), 0);
    uint32_t total_l_reads = std::accumulate(sizes_vecs.l_reads_count.begin(), sizes_vecs.l_reads_count.end(), 0);
    uint32_t total_ctg_len = std::accumulate(sizes_vecs.ctg_sizes.begin(), sizes_vecs.ctg_sizes.end(), 0);
    
    size_t gpu_mem_req = sizeof(int32_t) * tot_extensions * 6 + sizeof(int32_t) * total_l_reads
                           + sizeof(int32_t) * total_r_reads + sizeof(char) * total_ctg_len
                           + sizeof(char) * total_l_reads * max_read_size*2 + sizeof(char) * total_r_reads * max_read_size*2 // for quals included
                           + sizeof(double) * tot_extensions + sizeof(char) * total_r_reads * max_read_size 
                           + sizeof(char) * total_l_reads * max_read_size + sizeof(int64_t)*3
                           + sizeof(loc_ht)*ht_tot_size // changed to try the new method
                           + sizeof(char)*tot_extensions * max_walk_len
                           + (max_mer_len + max_walk_len) * sizeof(char) * tot_extensions
                           + sizeof(loc_ht_bool) * tot_extensions * max_walk_len;


  //  print_vals("Total GPU mem required (GBs):", (double)gpu_mem_req/(1024*1024*1024), "my_rank:",my_rank);                     
    size_t gpu_mem_avail = get_device_mem((ranks/total_gpus_avail), my_gpu_id);
    float factor = 0.90;
   // print_vals("GPU Mem using (MB):",((double)gpu_mem_avail*factor)/(1024*1024), "my_rank:",my_rank); 
    int iterations = ceil(((double)gpu_mem_req)/((double)gpu_mem_avail*factor)); // 0.8 is to buffer for the extra mem that is used when allocating once and using again
   // print_vals("Iterations:", iterations, "my_rank:",my_rank);
    unsigned slice_size = tot_extensions/iterations;
    unsigned remaining = tot_extensions % iterations;
    std::vector<uint32_t> max_ht_sizes;
//to get the largest ht size for any iteration and allocate GPU memory for that (once)
    uint32_t max_ht = 0, max_r_rds_its = 0, max_l_rds_its = 0, max_ctg_len_it = 0, test_sum = 0;
    for(int i = 0; i < iterations; i++){
        uint32_t temp_max_ht = 0, temp_max_r_rds = 0, temp_max_l_rds = 0, temp_max_ctg_len = 0;
        if(i < iterations -1 ){
            temp_max_ht = std::accumulate(sizes_vecs.ht_sizes.begin() + i*slice_size, sizes_vecs.ht_sizes.begin()+(i+1)*slice_size, 0 );
            temp_max_r_rds = std::accumulate(sizes_vecs.r_reads_count.begin() + i*slice_size, sizes_vecs.r_reads_count.begin()+(i+1)*slice_size, 0 );
            temp_max_l_rds = std::accumulate(sizes_vecs.l_reads_count.begin() + i*slice_size, sizes_vecs.l_reads_count.begin()+(i+1)*slice_size, 0 );
            temp_max_ctg_len = std::accumulate(sizes_vecs.ctg_sizes.begin() + i*slice_size, sizes_vecs.ctg_sizes.begin()+(i+1)*slice_size, 0 );
        }
        else{
            temp_max_ht = std::accumulate(sizes_vecs.ht_sizes.begin() + i*slice_size, sizes_vecs.ht_sizes.begin()+((i+1)*slice_size) + remaining, 0 );
            temp_max_r_rds = std::accumulate(sizes_vecs.r_reads_count.begin() + i*slice_size, sizes_vecs.r_reads_count.begin()+((i+1)*slice_size) + remaining, 0 );
            temp_max_l_rds = std::accumulate(sizes_vecs.l_reads_count.begin() + i*slice_size, sizes_vecs.l_reads_count.begin()+((i+1)*slice_size) + remaining, 0 );
            temp_max_ctg_len = std::accumulate(sizes_vecs.ctg_sizes.begin() + i*slice_size, sizes_vecs.ctg_sizes.begin()+((i+1)*slice_size) + remaining, 0 );
        }
        if(temp_max_ht > max_ht)
            max_ht = temp_max_ht;
        if(temp_max_r_rds > max_r_rds_its)
            max_r_rds_its = temp_max_r_rds;
        if(temp_max_l_rds > max_l_rds_its)
            max_l_rds_its = temp_max_l_rds; 
        if(temp_max_ctg_len > max_ctg_len_it)
            max_ctg_len_it = temp_max_ctg_len;
        test_sum += temp_max_ht;
    }

    slice_size = slice_size + remaining; // this is the largest slice size, mostly the last iteration handles the leftovers
    //allocating maximum possible memory for a single iteration
   // print_vals("slice size maximum:", slice_size, "my_rank:",my_rank);
    timer mem_timer;
    double cpu_mem_aloc_time = 0, gpu_mem_aloc_time = 0, cpu_mem_dealoc_time = 0, gpu_mem_dealoc_time = 0;

    mem_timer.timer_start();
    char *ctg_seqs_h = new char[max_ctg_size * slice_size];
    uint32_t *cid_h = new uint32_t[slice_size];
    char * ctgs_seqs_rc_h = new char[max_ctg_size * slice_size];// revcomps not requried on GPU, ctg space will be re-used on GPU, but if we want to do right left extensions in parallel, then we need separate space on GPU
    uint32_t *ctg_seq_offsets_h = new uint32_t[slice_size];
    double *depth_h = new double[slice_size];
    char *reads_left_h = new char[max_l_count * max_read_size * slice_size]; 
    char *reads_right_h = new char[max_r_count * max_read_size * slice_size];
    char *quals_right_h = new char[max_r_count * max_read_size * slice_size];
    char *quals_left_h = new char[max_l_count * max_read_size * slice_size];
    uint32_t *reads_l_offset_h = new uint32_t[max_l_count* slice_size];
    uint32_t *reads_r_offset_h = new uint32_t[max_r_count * slice_size];
    uint32_t *rds_l_cnt_offset_h = new uint32_t[slice_size];
    uint32_t *rds_r_cnt_offset_h = new uint32_t[slice_size];
    uint32_t *term_counts_h = new uint32_t[3];
    char* longest_walks_r_h = new char[slice_size * max_walk_len * iterations];// reserve memory for all the walks
    char* longest_walks_l_h = new char[slice_size * max_walk_len * iterations]; // not needed on device, will re-use right walk memory
    uint32_t* final_walk_lens_r_h = new uint32_t[slice_size * iterations]; // reserve memory for all the walks.
    uint32_t* final_walk_lens_l_h = new uint32_t[slice_size * iterations]; // not needed on device, will re use right walk memory
    uint32_t* prefix_ht_size_h = new uint32_t[slice_size];
    mem_timer.timer_end();
    cpu_mem_aloc_time += mem_timer.get_total_time();
    gpu_mem_req = sizeof(int32_t) * slice_size * 6 + sizeof(int32_t) * 3
                           + sizeof(int32_t) * max_l_rds_its
                           + sizeof(int32_t) * max_r_rds_its
                           + sizeof(char) * max_ctg_len_it
                           + sizeof(char) * max_l_rds_its * max_read_size * 2
                           + sizeof(char) * max_r_rds_its * max_read_size  * 2
                           + sizeof(double) * slice_size 
                           + sizeof(loc_ht) * max_ht
                           + sizeof(char) * slice_size * max_walk_len
                           + (max_mer_len + max_walk_len) * sizeof(char) * slice_size
                           + sizeof(loc_ht_bool) * slice_size * max_walk_len;

  //  print_vals("Device Mem requesting per slice (MB):", (double)gpu_mem_req/ (1024*1024), "my_rank:",my_rank);

    uint32_t *cid_d, *ctg_seq_offsets_d, *reads_l_offset_d, *reads_r_offset_d; 
    uint32_t *rds_l_cnt_offset_d, *rds_r_cnt_offset_d, *prefix_ht_size_d;
    char *ctg_seqs_d, *reads_left_d, *reads_right_d, *quals_left_d, *quals_right_d;
    char *longest_walks_d, *mer_walk_temp_d;
    double *depth_d;
    uint32_t *term_counts_d;
    loc_ht *d_ht;
    loc_ht_bool *d_ht_bool;
    uint32_t* final_walk_lens_d;
    mem_timer.timer_start();
    //allocate GPU  memory
    CUDA_CHECK(cudaMalloc(&prefix_ht_size_d, sizeof(uint32_t) * slice_size));
    CUDA_CHECK(cudaMalloc(&cid_d, sizeof(uint32_t) * slice_size));
    CUDA_CHECK(cudaMalloc(&ctg_seq_offsets_d, sizeof(uint32_t) * slice_size));
    CUDA_CHECK(cudaMalloc(&reads_l_offset_d, sizeof(uint32_t) * max_l_rds_its));// changed this with new max
    CUDA_CHECK(cudaMalloc(&reads_r_offset_d, sizeof(uint32_t) * max_r_rds_its)); // changed this with new max
    CUDA_CHECK(cudaMalloc(&rds_l_cnt_offset_d, sizeof(uint32_t) * slice_size));
    CUDA_CHECK(cudaMalloc(&rds_r_cnt_offset_d, sizeof(uint32_t) * slice_size));
    CUDA_CHECK(cudaMalloc(&ctg_seqs_d, sizeof(char) * max_ctg_len_it)); // changed this with new max
    CUDA_CHECK(cudaMalloc(&reads_left_d, sizeof(char) * max_read_size * max_l_rds_its)); // changed
    CUDA_CHECK(cudaMalloc(&reads_right_d, sizeof(char) * max_read_size * max_r_rds_its));//changed
    CUDA_CHECK(cudaMalloc(&depth_d, sizeof(double) * slice_size));
    CUDA_CHECK(cudaMalloc(&quals_right_d, sizeof(char) *max_read_size * max_r_rds_its));//changed this
    CUDA_CHECK(cudaMalloc(&quals_left_d, sizeof(char) * max_read_size * max_l_rds_its));//changed this with new
    CUDA_CHECK(cudaMalloc(&term_counts_d, sizeof(uint32_t)*3));
    // if we separate out kernels for right and left walks then we can use r_count/l_count separately but for now use the max of two
    // also subtract the appropriate kmer length from max_read_size to reduce memory footprint of global ht_loc.
    // one local hashtable for each thread, so total hash_tables equal to vec_size i.e. total contigs
    CUDA_CHECK(cudaMalloc(&d_ht, sizeof(loc_ht)*max_ht)); //**changed for new modifications
    CUDA_CHECK(cudaMalloc(&longest_walks_d, sizeof(char)*slice_size * max_walk_len));
    CUDA_CHECK(cudaMalloc(&mer_walk_temp_d, (max_mer_len + max_walk_len) * sizeof(char) * slice_size));
    CUDA_CHECK(cudaMalloc(&d_ht_bool, sizeof(loc_ht_bool) * slice_size * max_walk_len));
    CUDA_CHECK(cudaMalloc(&final_walk_lens_d, sizeof(uint32_t) * slice_size));
    mem_timer.timer_end();
    gpu_mem_aloc_time += mem_timer.get_total_time();

    timer loop_time;
    loop_time.timer_start();
    double data_mv_tim = 0;
    double packing_tim = 0;
    slice_size = tot_extensions/iterations;
    for(int slice = 0; slice < iterations; slice++){
        uint32_t left_over;
        if(iterations - 1 == slice)
            left_over = tot_extensions % iterations;
        else
            left_over = 0;
        std::vector<loc_assem_helper::CtgWithReads> slice_data (&data_in[slice*slice_size], &data_in[(slice + 1)*slice_size + left_over]);
        uint32_t vec_size = slice_data.size();
        uint32_t ctgs_offset_sum = 0;
        uint32_t prefix_ht_sum = 0;
        uint32_t reads_r_offset_sum = 0;
        uint32_t reads_l_offset_sum = 0;
        int read_l_index = 0, read_r_index = 0;
        timer tim_temp;
        tim_temp.timer_start();
        for(unsigned i = 0; i < slice_data.size(); i++){
            loc_assem_helper::CtgWithReads temp_data = slice_data[i];
            cid_h[i] = temp_data.cid;
            depth_h[i] = temp_data.depth;
            //convert string to c-string
            char *ctgs_ptr = ctg_seqs_h + ctgs_offset_sum;
            memcpy(ctgs_ptr, temp_data.seq.c_str(), temp_data.seq.size());
            ctgs_offset_sum += temp_data.seq.size();
            ctg_seq_offsets_h[i] = ctgs_offset_sum;
            prefix_ht_sum += temp_data.max_reads * max_read_size;
            prefix_ht_size_h[i] = prefix_ht_sum;

            for(int j = 0; j < temp_data.reads_left.size(); j++){
                char *reads_l_ptr = reads_left_h + reads_l_offset_sum;
                char *quals_l_ptr = quals_left_h + reads_l_offset_sum;
                memcpy(reads_l_ptr, temp_data.reads_left[j].seq.c_str(), temp_data.reads_left[j].seq.size());
                //quals offsets will be same as reads offset because quals and reads have same length
                memcpy(quals_l_ptr, temp_data.reads_left[j].quals.c_str(), temp_data.reads_left[j].quals.size());
                reads_l_offset_sum += temp_data.reads_left[j].seq.size();
                reads_l_offset_h[read_l_index] = reads_l_offset_sum;
                read_l_index++;
            }
            rds_l_cnt_offset_h[i] = read_l_index; // running sum of left reads count

            for(int j = 0; j < temp_data.reads_right.size(); j++){
                char *reads_r_ptr = reads_right_h + reads_r_offset_sum;
                char *quals_r_ptr = quals_right_h + reads_r_offset_sum;
                memcpy(reads_r_ptr, temp_data.reads_right[j].seq.c_str(), temp_data.reads_right[j].seq.size());
                //quals offsets will be same as reads offset because quals and reads have same length
                memcpy(quals_r_ptr, temp_data.reads_right[j].quals.c_str(), temp_data.reads_right[j].quals.size());
                reads_r_offset_sum += temp_data.reads_right[j].seq.size();
                reads_r_offset_h[read_r_index] = reads_r_offset_sum;
                read_r_index++;
            }
            rds_r_cnt_offset_h[i] = read_r_index; // running sum of right reads count
        }// data packing for loop ends
        tim_temp.timer_end();
        packing_tim += tim_temp.get_total_time();

        int total_r_reads_slice = read_r_index;
        int total_l_reads_slice = read_l_index;

        for(int i = 0; i < 3; i++){
            term_counts_h[i] = 0;
        }

        tim_temp.timer_start();

        CUDA_CHECK(cudaMemcpy(prefix_ht_size_d, prefix_ht_size_h, sizeof(uint32_t) * vec_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(cid_d, cid_h, sizeof(uint32_t) * vec_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ctg_seq_offsets_d, ctg_seq_offsets_h, sizeof(uint32_t) * vec_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(reads_l_offset_d, reads_l_offset_h, sizeof(uint32_t) * total_l_reads_slice, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(reads_r_offset_d, reads_r_offset_h, sizeof(uint32_t) * total_r_reads_slice, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(rds_l_cnt_offset_d, rds_l_cnt_offset_h, sizeof(uint32_t) * vec_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(rds_r_cnt_offset_d, rds_r_cnt_offset_h, sizeof(uint32_t) * vec_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ctg_seqs_d, ctg_seqs_h, sizeof(char) * ctgs_offset_sum, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(reads_left_d, reads_left_h, sizeof(char) * reads_l_offset_sum, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(reads_right_d, reads_right_h, sizeof(char) * reads_r_offset_sum, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(depth_d, depth_h, sizeof(double) * vec_size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(quals_right_d, quals_right_h, sizeof(char) * reads_r_offset_sum, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(quals_left_d, quals_left_h, sizeof(char) * reads_l_offset_sum, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(term_counts_d, term_counts_h, sizeof(uint32_t)*3, cudaMemcpyHostToDevice));
        tim_temp.timer_end();
        data_mv_tim += tim_temp.get_total_time();
            //call kernel here, one thread per contig
        unsigned total_threads = vec_size*32;// we need one warp (32 threads) per extension, vec_size = extensions
        unsigned thread_per_blk = 512;
        unsigned blocks = (total_threads + thread_per_blk)/thread_per_blk;
        
        int64_t sum_ext=0, num_walks=0;
        uint32_t qual_offset_ = qual_offset;
        iterative_walks_kernel<<<blocks,thread_per_blk>>>(cid_d, ctg_seq_offsets_d, ctg_seqs_d, reads_right_d, quals_right_d, reads_r_offset_d, rds_r_cnt_offset_d, 
        depth_d, d_ht, prefix_ht_size_d, d_ht_bool, mer_len, max_mer_len, term_counts_d, num_walks, max_walk_len, sum_ext, max_read_size, max_read_count, qual_offset_, longest_walks_d, mer_walk_temp_d, final_walk_lens_d, vec_size);

       // CUDA_CHECK(cudaDeviceSynchronize());

        //perform revcomp of contig sequences and launch kernel with left reads, 

        for(unsigned j = 0; j < vec_size; j++){
            int size_lst;
            char* curr_seq;
            char* curr_seq_rc;
            if(j == 0){
                size_lst = ctg_seq_offsets_h[j];
                curr_seq = ctg_seqs_h;
                curr_seq_rc = ctgs_seqs_rc_h;
            }
            else{
                size_lst = ctg_seq_offsets_h[j] - ctg_seq_offsets_h[j-1];
                curr_seq = ctg_seqs_h + ctg_seq_offsets_h[j - 1];
                curr_seq_rc = ctgs_seqs_rc_h + ctg_seq_offsets_h[j - 1];
            }
            loc_assem_helper::revcomp(curr_seq, curr_seq_rc, size_lst);
            #ifdef DEBUG_PRINT_CPU   
            print_vals("orig seq:");
            for(int h = 0; h < size_lst; h++)
                std::cout<<curr_seq[h];
            std::cout << std::endl;
            print_vals("recvomp seq:");
            for(int h = 0; h < size_lst; h++)
                std::cout<<curr_seq_rc[h];
            std::cout << std::endl;    
            #endif   

        }
        tim_temp.timer_start();
        CUDA_CHECK(cudaMemcpy(longest_walks_r_h + slice * max_walk_len * slice_size, longest_walks_d, sizeof(char) * vec_size * max_walk_len, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(final_walk_lens_r_h + slice * slice_size, final_walk_lens_d, sizeof(int32_t) * vec_size, cudaMemcpyDeviceToHost)); 

        //cpying rev comped ctgs to device on same memory as previous ctgs
        CUDA_CHECK(cudaMemcpy(ctg_seqs_d, ctgs_seqs_rc_h, sizeof(char) * ctgs_offset_sum, cudaMemcpyHostToDevice));
        tim_temp.timer_end();
        data_mv_tim += tim_temp.get_total_time();
        iterative_walks_kernel<<<blocks,thread_per_blk>>>(cid_d, ctg_seq_offsets_d, ctg_seqs_d, reads_left_d, quals_left_d, reads_l_offset_d, rds_l_cnt_offset_d, 
        depth_d, d_ht, prefix_ht_size_d, d_ht_bool, mer_len, max_mer_len, term_counts_d, num_walks, max_walk_len, sum_ext, max_read_size, max_read_count, qual_offset_, longest_walks_d, mer_walk_temp_d, final_walk_lens_d, vec_size);
   
        tim_temp.timer_start();
        CUDA_CHECK(cudaMemcpy(longest_walks_l_h + slice * max_walk_len * slice_size , longest_walks_d, sizeof(char) * vec_size * max_walk_len, cudaMemcpyDeviceToHost)); // copy back left walks
        CUDA_CHECK(cudaMemcpy(final_walk_lens_l_h + slice * slice_size , final_walk_lens_d, sizeof(int32_t) * vec_size, cudaMemcpyDeviceToHost)); 
        tim_temp.timer_end();
        data_mv_tim += tim_temp.get_total_time();
    }// the for loop over all slices ends here

    loop_time.timer_end();
   // print_vals("Total Loop Time:", loop_time.get_total_time(), "my_rank:",my_rank);

    //once all the alignments are on cpu, then go through them and stitch them with contigs in front and back.
    int loc_left_over = tot_extensions % iterations;
    for(int j = 0; j < iterations; j++){
        int loc_size = (j == iterations - 1) ? slice_size + loc_left_over : slice_size;

        //TODO: a lot of multiplications in below loop can be optimized (within indices)
        for(int i = 0; i< loc_size; i++){
            if(final_walk_lens_l_h[j*slice_size + i] > 0){
                std::string left(&longest_walks_l_h[j*slice_size*max_walk_len + max_walk_len*i],final_walk_lens_l_h[j*slice_size + i]);
                std::string left_rc = loc_assem_helper::revcomp(left);
                data_in[j*slice_size + i].seq.insert(0,left_rc);  
            }
            if(final_walk_lens_r_h[j*slice_size + i] > 0){
                std::string right(&longest_walks_r_h[j*slice_size*max_walk_len + max_walk_len*i],final_walk_lens_r_h[j*slice_size + i]);
                data_in[j*slice_size + i].seq += right;
            }
        }
    }


    mem_timer.timer_start();  
    CUDA_CHECK(cudaFree(prefix_ht_size_d));
    CUDA_CHECK(cudaFree(term_counts_d));
    CUDA_CHECK(cudaFree(cid_d));
    CUDA_CHECK(cudaFree(ctg_seq_offsets_d));
    CUDA_CHECK(cudaFree(reads_l_offset_d));
    CUDA_CHECK(cudaFree(reads_r_offset_d));
    CUDA_CHECK(cudaFree(rds_l_cnt_offset_d));
    CUDA_CHECK(cudaFree(rds_r_cnt_offset_d));
    CUDA_CHECK(cudaFree(ctg_seqs_d));
    CUDA_CHECK(cudaFree(reads_left_d));
    CUDA_CHECK(cudaFree(reads_right_d));
    CUDA_CHECK(cudaFree(depth_d));
    CUDA_CHECK(cudaFree(quals_right_d));
    CUDA_CHECK(cudaFree(quals_left_d));
    CUDA_CHECK(cudaFree(d_ht)); 
    CUDA_CHECK(cudaFree(longest_walks_d));
    CUDA_CHECK(cudaFree(mer_walk_temp_d));
    CUDA_CHECK(cudaFree(d_ht_bool));
    CUDA_CHECK(cudaFree(final_walk_lens_d));

    mem_timer.timer_end();
    gpu_mem_dealoc_time += mem_timer.get_total_time();
   // print_vals("mem freed, from rank",my_rank);

    mem_timer.timer_start();
    delete[] ctg_seqs_h;
    delete[] cid_h;
    delete[] ctgs_seqs_rc_h;// revcomps not requried on GPU, ctg space will be re-used on GPU, but if we want to do right left extensions in parallel, then we need separate space on GPU
    delete[] ctg_seq_offsets_h;
    delete[] depth_h;
    delete[] reads_left_h; 
    delete[] reads_right_h;
    delete[] quals_right_h;
    delete[] quals_left_h;
    delete[] reads_l_offset_h;
    delete[] reads_r_offset_h;
    delete[] rds_l_cnt_offset_h;
    delete[] rds_r_cnt_offset_h;
    delete[] term_counts_h;
    delete[] longest_walks_r_h;// reserve memory for all the walks
    delete[] longest_walks_l_h; // not needed on device, will re-use right walk memory
    delete[] final_walk_lens_r_h; // reserve memory for all the walks.
    delete[] final_walk_lens_l_h; // not needed on device, will re use right walk memory
    mem_timer.timer_end();
    cpu_mem_dealoc_time += mem_timer.get_total_time();

}
