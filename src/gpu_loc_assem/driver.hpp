#include "locassem_struct.hpp"
namespace locassm_driver{
    struct accum_data{
        std::vector<uint32_t> ht_sizes;
        std::vector<uint32_t> l_reads_count;
        std::vector<uint32_t> r_reads_count;
        std::vector<uint32_t> ctg_sizes;
    };

    struct ctg_bucket{
        std::vector<CtgWithReads> ctg_vec;
        accum_data sizes_vec;
        uint32_t l_max, r_max, max_contig_sz;
	ctg_bucket():l_max{0}, r_max{0}, max_contig_sz{0}{}
        void clear();
    };

    size_t get_device_mem(int ranks_per_gpu, int gpu_id);
    void local_assem_driver(std::vector<CtgWithReads>& data_in, uint32_t max_ctg_size, uint32_t max_read_size, uint32_t max_r_count, uint32_t max_l_count, int mer_len, int max_kmer_len, accum_data& sizes_outliers, int walk_len_limit, int qual_offset, int ranks, int my_rank, int g_rank_me);
    int get_gpu_per_node();

}
