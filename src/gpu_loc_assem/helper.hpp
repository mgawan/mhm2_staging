
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>


//templated functions needs to be in the same translation unit
template<typename T>
void print_log(std::ofstream &file_n, T _log){
    file_n<<_log<<std::endl;
}

template<typename T>
void print_vals(std::ofstream &file_n, T val){
    print_log(file_n, val);
}

template<typename T, typename... Types>
void print_vals(std::ofstream &file_n, T val, Types... val_){
    if(sizeof...(val_) == 0){
        print_vals(file_n, val);
    }else{
        print_vals(file_n, val);
        print_vals(file_n, val_...);
        }
}


struct timer{
double total_time = 0;
std::chrono::time_point<std::chrono::high_resolution_clock> time_begin;
std::chrono::time_point<std::chrono::high_resolution_clock> time_end;
std::chrono::duration<double> diff;

void timer_start(){
    time_begin = std::chrono::high_resolution_clock::now();
}

void timer_end(){
    time_end = std::chrono::high_resolution_clock::now();
}

double get_total_time(){
    diff = time_end - time_begin;
    return diff.count();
}


};

namespace loc_assem_helper{

inline void revcomp(char* str, char* str_rc, int size, int rank) {
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
        std::cout<<"Illegal char:"<< str[i]<< " in revcomp from rank:"<<rank<<", printing string: \n";
	for(auto k = 0; k < size; k++) std::cout<<str[k];
	std::cout<<std::endl;
	break;
    }
    size_rc++;
  }
}

inline std::string revcomp(std::string instr, int rank) {
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
        std::cout<<"Illegal char:"<<instr[i]<< "in revcomp of(string) from rank:"<<rank<<" \n";
	std::cout<<"string:"<<instr<<std::endl;
	break;
    }
  }

  return str_rc;
}

}
