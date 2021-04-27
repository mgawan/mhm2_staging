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

   inline void revcomp(char* str, char* str_rc, int size) {
      int size_rc = 0;
      for (int i = size - 1; i >= 0; i--) {
      	switch (str[i]) {
     	   case 'A': str_rc[size_rc]= 'T'; break;
      	   case 'C': str_rc[size_rc]= 'G'; break;
           case 'G': str_rc[size_rc]= 'C'; break;
           case 'T': str_rc[size_rc]= 'A'; break;
           case 'N': str_rc[size_rc]= 'N'; break;
           case 'U': case 'R': case 'Y': case 'K': case 'M': case 'S': case 'W': case 'B': case 'D': case 'H': case 'V':
             str_rc[size_rc]= 'N';
             break;
           default:
             std::cout<<"Illegal char:"<< str[i]<< "\n";
	     break;
         }
       size_rc++;
      }
   }

    inline std::string revcomp(std::string instr) {
  	std::string str_rc;
  	for (int i = instr.size() - 1; i >= 0; i--) {
    	    switch (instr[i]) {
      		case 'A': str_rc += 'T'; break;
      		case 'C': str_rc += 'G'; break;
	        case 'G': str_rc += 'C'; break;
      		case 'T': str_rc += 'A'; break;
      		case 'N': str_rc += 'N'; break;
      		case 'U': case 'R': case 'Y': case 'K': case 'M': case 'S': case 'W': case 'B': case 'D': case 'H': case 'V':
                  str_rc += 'N';
                  break;
      		default:
                  std::cout<<"Illegal char:"<<instr[i]<<"\n";
		  break;
   	    }
 	 }

  	return str_rc;
	}

    size_t get_device_mem(int ranks_per_gpu, int gpu_id);
    void local_assem_driver(std::vector<CtgWithReads>& data_in, uint32_t max_ctg_size, uint32_t max_read_size, uint32_t max_r_count, uint32_t max_l_count, int mer_len, int max_kmer_len, accum_data& sizes_outliers, int walk_len_limit, int qual_offset, int ranks, int my_rank, int g_rank_me);
    int get_gpu_per_node();

}
