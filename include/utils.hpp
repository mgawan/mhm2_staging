#pragma once

/*
 HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from the U.S. Dept. of Energy).  All rights reserved."

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 (1) Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 (2) Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 (3) Neither the name of the University of California, Lawrence Berkeley National
 Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 endorse or promote products derived from this software without specific prior
 written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 to the features, functionality or performance of the source code ("Enhancements") to
 anyone; however, if you choose to make your Enhancements available either publicly,
 or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 written license agreement for such Enhancements, then you hereby grant the following
 license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 prepare derivative works, incorporate into other computer software, distribute, and
 sublicense such enhancements or derivative works thereof, in binary and source code
 form.
*/


#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <utility>
#include <algorithm>
#include <sstream>
#include <dirent.h>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/timers.hpp"
#include "upcxx_utils/ofstream.hpp"

using namespace upcxx_utils;

using std::string;
using std::string_view;
using std::to_string;
using std::min;
using std::vector;
using std::pair;

#ifdef USE_BYTELL
#include "bytell_hash_map.hpp"
#define HASH_TABLE ska::bytell_hash_map
#else
#include <unordered_map>
#define HASH_TABLE std::unordered_map
#endif

inline size_t estimate_hashtable_memory(size_t num_elements, size_t element_size) {
    // get the hashtable load factor
    HASH_TABLE<char,char> tmp;
    double max_load_factor = tmp.max_load_factor();
    
    // apply the load factor
    size_t expanded_num_elements = num_elements / max_load_factor + 1;
    
    // get the next power of two
    --expanded_num_elements;
    expanded_num_elements |= expanded_num_elements >> 1;
    expanded_num_elements |= expanded_num_elements >> 2;
    expanded_num_elements |= expanded_num_elements >> 4;
    expanded_num_elements |= expanded_num_elements >> 8;
    expanded_num_elements |= expanded_num_elements >> 16;
    expanded_num_elements |= expanded_num_elements >> 32;
    ++expanded_num_elements;

    size_t num_buckets = expanded_num_elements * max_load_factor; 
    
    return expanded_num_elements * element_size + num_buckets * 8;
}

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
      case 'U': case 'R': case 'Y': case 'K': case 'M': case 'S': case 'W': case 'B': case 'D': case 'H': case 'V':
        seq_rc += 'N';
        break;
      default:
        DIE("Illegal char '", seq[i], "' in revcomp of '", seq, "'");
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
      case 'U': case 'R': case 'Y': case 'K': case 'M': case 'S': case 'W': case 'B': case 'D': case 'H': case 'V':
        return 'N';
      default:
        DIE("Illegal char '", ch, "' in comp nucleotide");
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

static string get_merged_reads_fname(string reads_fname) {
  // always relative to the current working directory
  if (reads_fname.find(':') != string::npos) {
      // remove the first pair, if it exists
      reads_fname = reads_fname.substr(reads_fname.find(':'));
  }
  return upcxx_utils::remove_file_ext(get_basename(reads_fname)) + "-merged.fastq";
}

inline void switch_orient(int &start, int &stop, int &len) {
  int tmp = start;
  start = len - stop;
  stop = len - tmp;
}

inline void dump_single_file(const string &fname, const string &out_str, bool append=false) {
  BarrierTimer timer(__FILEFUNC__);
  SLOG_VERBOSE("Writing ", fname, "\n");
  auto fut_tot_bytes_written = upcxx::reduce_one(out_str.size(), upcxx::op_fast_add, 0);
  upcxx_utils::dist_ofstream of(fname, append);
  of << out_str;
  of.close();
  SLOG_VERBOSE("Successfully wrote ", get_size_str(fut_tot_bytes_written.wait()), " bytes to ", fname, "\n");
  assert(rank_me() || of.get_last_known_tellp() == fut_tot_bytes_written.wait());
}

inline vector<string> get_dir_entries(const string &dname, const string &prefix) {
  vector<string> dir_entries;
  DIR *dir = opendir(dname.c_str());
  if (dir) {
    struct dirent *en;
    while ((en = readdir(dir)) != NULL) {
      string dname(en->d_name);
      if (dname.substr(0, prefix.length()) == prefix) dir_entries.push_back(dname);
    }
    closedir(dir);  // close all directory
  } else {
    SWARN("Could not open ", dname);
  }
  return dir_entries;
}

inline std::string &left_trim(std::string &str) {
  auto it = std::find_if(str.begin(), str.end(), [](char ch) { return !std::isspace<char>(ch, std::locale()); });
  str.erase(str.begin(), it);
  return str;
}

inline int pin_clear() {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  for (unsigned i = 0; i < sizeof(cpu_set_t) * 8; i++) {
    CPU_SET(i, &cpu_set);
  }
  if (sched_setaffinity(getpid(), sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), getpid());
    return -1;
  }
  return 0;
}

inline string get_proc_pin() {
  ifstream f("/proc/self/status");
  string line;
  string prefix = "Cpus_allowed_list:";
  while (getline(f, line)) {
    if (line.substr(0, prefix.length()) == prefix) {
      DBG(line, "\n");
      line = line.substr(prefix.length(), line.length() - prefix.length());
      return left_trim(line);
      break;
    }
  }
  return "";
}

inline vector<int> get_pinned_cpus() {
  vector<int> cpus;
  stringstream ss(get_proc_pin());
  while (ss.good()) {
    string s;
    getline(ss, s, ',');
    s = left_trim(s);
    auto dash_pos = s.find('-');
    if (dash_pos != string::npos) {
      int start = std::stoi(s.substr(0, dash_pos));
      int stop = std::stoi(s.substr(dash_pos + 1)) + 1;
      for (int i = start; i < stop; i++) cpus.push_back(i);
    } else {
      cpus.push_back(std::stoi(s));
    }
  }
  return cpus;
}

inline void pin_proc(vector<int> cpus) {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  for (auto cpu : cpus) {
    CPU_SET(cpu, &cpu_set);
  }
  if (sched_setaffinity(getpid(), sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), getpid());
  }
}

inline void pin_cpu() {
  auto pinned_cpus = get_pinned_cpus();
  pin_proc({pinned_cpus[upcxx::rank_me() % pinned_cpus.size()]});
  SLOG("Pinning to logical cpus: process 0 on node 0 pinned to cpu ", get_proc_pin(), "\n");
}

inline void pin_core() {
  string numa_node_dir = "/sys/devices/system/node";
  auto numa_node_entries = get_dir_entries(numa_node_dir, "node");
  if (numa_node_entries.empty()) return;
  vector<int> my_thread_siblings;
  for (auto &entry : numa_node_entries) {
    ifstream f(numa_node_dir + "/" + entry + "/cpulist");
    string buf;
    getline(f, buf);
    f.close();
    int numa_node_i = std::stoi(entry.substr(4));
    auto cpu_entries = get_dir_entries(numa_node_dir + "/" + entry, "cpu");
    for (auto &cpu_entry : cpu_entries) {
      if (cpu_entry != "cpu" + to_string(upcxx::rank_me())) continue;
      if (cpu_entry == "cpulist" || cpu_entry == "cpumap") continue;
      f.open(numa_node_dir + "/" + entry + "/" + cpu_entry + "/topology/thread_siblings_list");
      getline(f, buf);
      stringstream ss(buf);
      while (ss.good()) {
        string s;
        getline(ss, s, ',');
        my_thread_siblings.push_back(std::stoi(s));
      }
      break;
    }
  }
  if (!my_thread_siblings.empty()) {
    pin_proc(my_thread_siblings);
    SLOG("Pinning to cores: process 0 on node 0 pinned to cpus ", get_proc_pin(), "\n");
  }
}

inline void pin_numa() {
  string numa_node_dir = "/sys/devices/system/node";
  auto numa_node_entries = get_dir_entries(numa_node_dir, "node");
  if (numa_node_entries.empty()) return;
  vector<pair<string, vector<int>>> numa_node_list(numa_node_entries.size(), {"", {}});
  int num_cpus = 0;
  int hdw_threads_per_core = 0;
  for (auto &entry : numa_node_entries) {
    ifstream f(numa_node_dir + "/" + entry + "/cpulist");
    string buf;
    getline(f, buf);
    f.close();
    int numa_node_i = std::stoi(entry.substr(4));
    numa_node_list[numa_node_i].first = buf;
    auto cpu_entries = get_dir_entries(numa_node_dir + "/" + entry, "cpu");
    for (auto &cpu_entry : cpu_entries) {
      if (cpu_entry == "cpulist" || cpu_entry == "cpumap") continue;
      if (!hdw_threads_per_core) {
        f.open(numa_node_dir + "/" + entry + "/" + cpu_entry + "/topology/thread_siblings_list");
        getline(f, buf);
        // assume that the threads are separated by commans - is this always true?
        hdw_threads_per_core = std::count(buf.begin(), buf.end(), ',') + 1;
      }
      numa_node_list[numa_node_i].second.push_back(std::stoi(cpu_entry.substr(3)));
      num_cpus++;
    }
  }
  SLOG("On node 0, found a total of ", num_cpus, " hardware threads with ", hdw_threads_per_core,
       " threads per core on ", numa_node_list.size(), " NUMA domains\n");
  // pack onto numa nodes
  int hdw_threads_per_numa_node = num_cpus / numa_node_list.size();
  int cores_per_numa_node = hdw_threads_per_numa_node / hdw_threads_per_core;
  int numa_nodes_to_use = upcxx::local_team().rank_n() / cores_per_numa_node;
  if (numa_nodes_to_use > numa_node_list.size()) numa_nodes_to_use = numa_node_list.size();
  if (numa_nodes_to_use == 0) numa_nodes_to_use = 1;
  int my_numa_node = upcxx::local_team().rank_me() % numa_nodes_to_use;
  vector<int> my_cpu_list = numa_node_list[my_numa_node].second;
  sort(my_cpu_list.begin(), my_cpu_list.end());
  pin_proc(my_cpu_list);
  SLOG("Pinning to ", numa_nodes_to_use, " NUMA domains each with ", cores_per_numa_node, " cores, ", hdw_threads_per_numa_node, " cpus: process 0 on node 0 is pinned to cpus ", get_proc_pin(), "\n");
}


// temporary until it is properly within upcxx_utils
#include <thread>
#include <memory>
#include <functional>
#include <deque>

namespace upcxx_utils {

    // Func no argument returned or given lambda - void()

    template<typename Func>
    future<> execute_in_new_thread(upcxx::persona &persona, Func func) {
        assert(persona.active_with_caller());
        shared_ptr< promise<> > sh_prom = make_shared<promise <> > ();

        shared_ptr<std::thread> sh_run = make_shared<std::thread>(
                [&persona, func, sh_prom] {
                    func();

                    // fulfill only in calling persona
                    persona.lpc_ff([&persona, sh_prom]() {
                        assert(persona.active_with_caller());
                        sh_prom->fulfill_anonymous(1);
                    });
                });
        auto fut_finished = sh_prom->get_future().then( // keep sh_prom and sh_run alive until complete
                [&persona, sh_prom, sh_run] () {
                    assert(persona.active_with_caller());
                    sh_run->join();
                });
        return fut_finished;
    }

    template<typename Func>
    future<> execute_in_new_thread(Func func) {
        return execute_in_new_thread(upcxx::current_persona(), func);
    }
    
    
    
    using LimitedFutureQueue = std::deque< upcxx::future<> >;
    
    inline LimitedFutureQueue &_get_outstanding_queue() {
        static LimitedFutureQueue outstanding_queue;
        return outstanding_queue;
    }
    
    inline upcxx::future<> limit_outstanding_futures(upcxx::future<> fut, int limit = 0, LimitedFutureQueue & outstanding_queue = upcxx_utils::_get_outstanding_queue()) {
        if (limit == 0) limit = upcxx::local_team().rank_n() * 2;
        if (limit == -1) limit = 0;
        outstanding_queue.push_back(fut);
        while (outstanding_queue.size() > limit) {
            auto fut = outstanding_queue.front();
            outstanding_queue.pop_front();
            outstanding_queue.front() = when_all(fut, outstanding_queue.front());
        }
        assert(outstanding_queue.size() <= limit);
        if (outstanding_queue.size() == limit) {
            return outstanding_queue.front();
        } else {
            return make_future();
        }
    }
    
    
    inline future<> flush_outstanding_futures_async(LimitedFutureQueue & outstanding_queue = upcxx_utils::_get_outstanding_queue()) {
        auto all_fut = make_future();
        while (!outstanding_queue.empty()) {
            auto fut = outstanding_queue.front();
            outstanding_queue.pop_front();
            all_fut = when_all(all_fut, fut);
        }
        assert(outstanding_queue.empty());
        return all_fut;
    }
    
    inline void flush_outstanding_futures(LimitedFutureQueue & outstanding_queue = upcxx_utils::_get_outstanding_queue()) {
        flush_outstanding_futures_async(outstanding_queue).wait();
    }
    
    template<typename Result, typename Future>
    upcxx::future<> assign_oustanding_future_result(Result &res, Future fut, int limit = 0, LimitedFutureQueue & outstanding_queue = upcxx_utils::_get_outstanding_queue()) {
        upcxx::future<> res_fut = fut.then(
                [&res](Result val) 
                {
                    res = val;
                });
        return limit_outstanding_futures(res_fut, limit);
    }
    

};
