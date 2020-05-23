#pragma once
#include <algorithm>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <utility>
#include <algorithm>

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

using cpu_set_size_t = std::pair<cpu_set_t *, size_t>;
inline cpu_set_size_t get_cpu_mask(bool bySocket = true) {
    cpu_set_size_t ret = {NULL, 0};
    ifstream cpuinfo("/proc/cpuinfo");
    if(!cpuinfo) {
        return ret;
    }
    std::vector<size_t> cpu2socket;
    std::vector<size_t> cpu2core;
    std::vector<size_t> sockets;
    std::vector<size_t> cores;
    cpu2socket.reserve(256);
    cpu2core.reserve(256);
    int socket = 0;
    for( std::string line; getline( cpuinfo, line ); ) {
        if (line.find("physical id") != string::npos) {
            int val = atoi(line.c_str() + line.find_last_of(' '));
            cpu2socket.push_back(val);
            socket = val;
        }
        if (line.find("core id") != string::npos) {
            int val = atoi(line.c_str() + line.find_last_of(' '));
            cpu2core.push_back(val + socket * 16384);
        }
    }
    for(auto id : cpu2core) {
      auto p = std::find(cores.begin(), cores.end(), id);
      if (p == cores.end()) cores.push_back(id);
    }
    for(auto id : cpu2socket) {
      auto p = std::find(sockets.begin(), sockets.end(), id);
      if (p == sockets.end()) sockets.push_back(id);
    }
      
    int num_cpus = cpu2core.size();
    int num_cores = cores.size();
    int num_sockets = sockets.size();
    int num_ids = bySocket ? num_sockets : num_cores;
    int my_id = upcxx::local_team().rank_me() % num_ids;
    DBG("Binding to ", bySocket ? "socket" : "core", " ", my_id, " of ", num_ids, " (num_cores=", num_cores, ", num_sockets=", num_sockets, ")\n");
    size_t size = CPU_ALLOC_SIZE(num_cpus);
    cpu_set_t *cpu_set_p = CPU_ALLOC(num_cpus);
    if (cpu_set_p == NULL) return ret;
    CPU_ZERO_S(size, cpu_set_p);
    for(int i = 0 ; i < cpu2socket.size(); i++) {
        if ((bySocket ? cpu2socket[i] : cpu2core[i]) == (bySocket ? sockets[my_id] : cores[my_id])) {
            CPU_SET_S(i, size, cpu_set_p);
        }
    }
    ret = {cpu_set_p, size};
    return ret;
}

inline int pin_mask(cpu_set_size_t cpu_set_size) {
    if (cpu_set_size.first) {
        if (sched_setaffinity(getpid(), cpu_set_size.second, cpu_set_size.first) == -1) {
            if (errno == 3) SWARN("%s, pid: %d", strerror(errno), getpid());
            CPU_FREE(cpu_set_size.first);
            return -1;
        }
        CPU_FREE(cpu_set_size.first);
        return 0;
    }
    SWARN("Did not pin to process\n");
    return -1;
}

inline int pin_socket() {
    return pin_mask(get_cpu_mask(true));
}

inline int pin_core() {
    return pin_mask(get_cpu_mask(false));
}

