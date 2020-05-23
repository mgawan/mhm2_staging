#pragma once
#include <algorithm>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unistd.h>

#include "upcxx_utils/log.hpp"
using namespace upcxx_utils;

using std::string;
using std::string_view;
using std::to_string;
using std::min;

#ifdef USE_BYTELL
#include "bytell_hash_map.hpp"
#define HASH_TABLE ska::bytell_hash_map
#else
#include <unordered_map>
#define HASH_TABLE std::unordered_map
#endif

// this shouldn't really be defined here, but I didn't want yet another header file
enum class QualityLevel {
  SINGLE_PATH_ONLY,
  DEPTH_RESLN_ONLY,
  ALL
};


inline string revcomp(const string &seq) {
  string seq_rc = "";
  seq_rc.reserve(seq.size());
  for (int i = seq.size() - 1; i >= 0; i--) {
    switch (seq[i]) {
      case 'A': seq_rc += 'T'; break;
      case 'C': seq_rc += 'G'; break;
      case 'G': seq_rc += 'C'; break;
      case 'T': seq_rc += 'A'; break;
      case 'N': seq_rc += 'N'; break;
      default:
        DIE("Illegal char in revcomp of '", seq, "'\n");
    }
  }
  return seq_rc;
}

inline char comp_nucleotide(char ch) {
  switch (ch) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      case 'N': return 'N';
      case '0': return '0';
      default: DIE("Illegal char in comp nucleotide of '", ch, "'\n");
  }
  return 0;
}

inline int hamming_dist(string_view s1, string_view s2, bool require_equal_len=true) {
  if (require_equal_len && s2.size() != s1.size())//abs((int)(s2.size() - s1.size())) > 1)
    DIE("Hamming distance substring lengths don't match, ", s1.size(), ", ", s2.size(), "\n");
  int d = 0;
  int min_size = min(s1.size(), s2.size());
  for (int i = 0; i < min_size; i++)
    d += (s1[i] != s2[i]);
  return d;
}

static string get_merged_reads_fname(const string &reads_fname) {
  // always relative to the current working directory
  return upcxx_utils::remove_file_ext(get_basename(reads_fname)) + "-merged.fastq";
}

inline void switch_orient(int &start, int &stop, int &len) {
  int tmp = start;
  start = len - stop;
  stop = len - tmp;
}

inline int pin_thread(pid_t pid, int cid) {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(cid, &cpu_set);
  if (sched_setaffinity(pid, sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), pid);
    return -1;
  }
  return 0;
}

inline int pin_socket() {
    ifstream cpuinfo("/proc/cpuinfo");
    if(!cpuinfo) {
        return -1;
    }
    std::vector<uint8_t> cpu2socket;
    int maxsocket = -1;
    const string id = "physical id";
    cpu2socket.reserve(256);
    for( std::string line; getline( cpuinfo, line ); ) {
        if (line.find(id) != string::npos) {
            int socket = atoi(line.c_str() + line.find_last_of(' '));
            cpu2socket.push_back(socket);
            if(socket > maxsocket) maxsocket = socket;
        }
    }
    if (maxsocket >= 0) {
        // group local_team ranks by socket (not round robin)
        int num_cpus = cpu2socket.size();
        int my_socket = upcxx::local_team().rank_me() * (maxsocket+1) / num_cpus;
        DBG("Binding to socket ", my_socket, " of ", maxsocket+1, "\n");
        size_t size = CPU_ALLOC_SIZE(num_cpus);
        cpu_set_t *cpu_set_p = CPU_ALLOC(num_cpus);
        if (cpu_set_p == NULL) return -1;
        CPU_ZERO_S(size, cpu_set_p);
        for(int i = 0 ; i < cpu2socket.size(); i++) {
            if (cpu2socket[i] == my_socket) {
                CPU_SET_S(i, size, cpu_set_p);
            }
        }
        if (sched_setaffinity(getpid(), size, cpu_set_p) == -1) {
            if (errno == 3) SWARN("%s, pid: %d", strerror(errno), getpid());
            CPU_FREE(cpu_set_p);
            return -1;
        }
        CPU_FREE(cpu_set_p);
        return 0;
    }
    SWARN("Did not pin to socket\n");
    return -1;
}

