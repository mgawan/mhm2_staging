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

// FIXME: this is copied from upcxx-utils, DistributedIO branch. It's temporary until Rob finishes his distributed IO implementation
template <typename T, typename BinaryOp>
future<> tmp_reduce_prefix_ring(const T *src, T *dst, size_t count, BinaryOp &op, const upcxx::team &team = upcxx::world(),
                                bool return_final_to_first = false) {
  if (team.from_world(rank_me(), upcxx::rank_n()) == upcxx::rank_n())
    throw std::runtime_error("reduce_prefix called outside of given team");  // not of this team
  DBG("reduce_prefix(src=", src, ", dst=", dst, ", count=", count, ", team.rank_n()=", team.rank_n(),
      ", return_final_to_first=", return_final_to_first, "\n");
  using ShPromise = std::shared_ptr<upcxx::promise<>>;
  using Data = std::tuple<const T *, T *, size_t, ShPromise>;
  using DistData = upcxx::dist_object<Data>;
  using ShDistData = std::shared_ptr<DistData>;

  ShPromise my_prom = std::make_shared<upcxx::promise<>>();
  upcxx::future<> ret_fut = my_prom->get_future();

  if (team.rank_me() == 0) {
    // first is special case of copy without op
    for (size_t i = 0; i < count; i++) dst[i] = src[i];

    if (return_final_to_first) {
      // make ready for the rpc but do not fulfill my_prom yet
      ret_fut = upcxx::make_future();
    } else {
      // make ready
      my_prom->fulfill_anonymous(1);
    }

    // special case
    if (team.rank_n() == 1) {
      if (return_final_to_first) my_prom->fulfill_anonymous(1);
      return ret_fut;
    }
  }

  // create a distributed object holding the data and promise
  ShDistData sh_dist_data = std::make_shared<DistData>(std::make_tuple(src, dst, count, my_prom), team);
  upcxx::intrank_t next_rank = team.rank_me() + 1;
  if (return_final_to_first || next_rank != team.rank_n()) {
    // send rpc to the next rank when my prefix is ready
    ret_fut = ret_fut.then([next_rank, dst, count, sh_dist_data, &op, &team]() -> future<> {
      //DBG("mydst is ready:", dst, ", next_rank=", next_rank, "\n");
      upcxx::rpc_ff(
          team, next_rank % team.rank_n(),
          [&op](DistData &dist_data, upcxx::view<T> prev_prefix, bool just_copy) {
            DBG("Got here: count=", prev_prefix.size(), ", just_copy=", just_copy, "\n");
            const T *mysrc;
            T *mydst;
            size_t mycount;
            std::shared_ptr<upcxx::promise<>> myprom;
            std::tie(mysrc, mydst, mycount, myprom) = *dist_data;
            assert(mycount == prev_prefix.size());
            if (just_copy) {
              for (size_t i = 0; i < mycount; i++) {
                mydst[i] = prev_prefix[i];
              }
            } else {
              for (size_t i = 0; i < mycount; i++) {
                mydst[i] = op(prev_prefix[i], mysrc[i]);
              }
            }
            myprom->fulfill_anonymous(1);  // mydst is ready
          },
          *sh_dist_data, upcxx::make_view(dst, dst + count), next_rank == team.rank_n());
      return upcxx::make_future();
    });
  } else {
    // keep the scope of the DistData object until ready
    ret_fut = ret_fut.then([sh_dist_data]() {
      /* noop */
    });
  }

  if (team.rank_me() == 0 && return_final_to_first) {
    // wait on my_prom to be fulfilled by an rpc from the last rank of the team
    // and keep the scope of the DistData object until ready
    ret_fut = my_prom->get_future().then([sh_dist_data]() {
      /* noop */
    });
  }

  return ret_fut;
};

inline void dump_single_file(const string &fname, const string &out_str, bool append=false) {
  BarrierTimer timer(__FILEFUNC__);
  // write to a temporary file and rename it on completion to ensure that there are no corrupted files should
  // there be a crash
  auto tmp_fname = fname + ".tmp";
  if (append && !upcxx::rank_me()) {
    // rename file to be appended to to temporary
    if (rename(fname.c_str(), tmp_fname.c_str()) != 0) DIE("Could not rename ", fname, " to temporary: ", strerror(errno));
  }
  upcxx::barrier();
  size_t offsets[2];
  auto sz = out_str.length();
  size_t append_file_size = (append && !upcxx::rank_me() ? get_file_size(tmp_fname) : 0);
  offsets[0] = (sz + append_file_size);
  tmp_reduce_prefix_ring(offsets, offsets + 1, 1, upcxx::op_fast_add).wait();
  // the offset is actually the start of this rank's block
  auto my_fpos = offsets[1] - sz;
  upcxx::barrier();
  //WARN("rank ", upcxx::rank_me(), " writing ", sz, " bytes to file ", tmp_fname, " at position ", my_fpos);
  int fileno = -1;
  size_t fsize = 0;
  // last rank creates the file and truncates it to the correct length (only last rank knows the correct total size)
  if (upcxx::rank_me() == upcxx::rank_n() - 1) {
    // rename previous file to tmp file
    fileno = open(tmp_fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) DIE("Could not open file ", tmp_fname, ": ", strerror(errno));
    if (ftruncate(fileno, my_fpos + sz) == -1) WARN("Could not truncate ", tmp_fname, " to ", my_fpos + sz, " bytes");
  }
  upcxx::barrier();
  // wait until rank n-1 has finished setting up the file
  if (rank_me() != rank_n() - 1) fileno = open(tmp_fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) DIE("Could not open file ", tmp_fname, ": ", strerror(errno));
  auto bytes_written = pwrite(fileno, out_str.c_str(), sz, my_fpos);
  close(fileno);
  if (bytes_written != sz) DIE("Could not write all ", sz, " bytes; only wrote ", bytes_written);
  upcxx::barrier();
  if (!upcxx::rank_me() && rename(tmp_fname.c_str(), fname.c_str()) != 0)
	DIE("Could not rename temporary file ", tmp_fname, " to ", fname, ", error: ", strerror(errno));
  auto tot_bytes_written = upcxx::reduce_one(bytes_written, upcxx::op_fast_add, 0).wait();
  SLOG_VERBOSE("Successfully wrote ", get_size_str(tot_bytes_written), " bytes to ", fname, "\n");
  upcxx::barrier();
}

template <typename T>
inline string vec_to_str(const vector<T> &vec, const string &delimiter = ",") {
  std::ostringstream oss;
  for (auto elem : vec) {
    oss << elem;
    if (elem != vec.back()) oss << delimiter;
  }
  return oss.str();
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

inline string get_proc_pin() {
  ifstream f("/proc/" + to_string(getpid()) + "/status");
  string line;
  string prefix = "Cpus_allowed_list";
  while (getline(f, line)) {
    if (line.substr(0, prefix.length()) == prefix) {
      DBG(line);
      return line.substr(prefix.length(), line.length() - prefix.length());
      break;
    }
  }
  return "";
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
  pin_proc({upcxx::rank_me()});
  SLOG("Pinning to logical cpus: process 0 on node 0 pinned to cpu", get_proc_pin(), "\n");
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
    int numa_node_i = atoi(entry.c_str() + 4);
    auto cpu_entries = get_dir_entries(numa_node_dir + "/" + entry, "cpu");
    for (auto &cpu_entry : cpu_entries) {
      if (cpu_entry != "cpu" + to_string(upcxx::rank_me())) continue;
      if (cpu_entry == "cpulist" || cpu_entry == "cpumap") continue;
      f.open(numa_node_dir + "/" + entry + "/" + cpu_entry + "/topology/thread_siblings_list");
      getline(f, buf);
      stringstream iss(buf);
      while (iss.good()) {
        string substr;
        getline(iss, substr, ',');
        my_thread_siblings.push_back(atoi(substr.c_str()));
      }
      break;
    }
  }
  if (!my_thread_siblings.empty()) {
    pin_proc(my_thread_siblings);
    SLOG("Pinning to cores: process 0 on node 0 pinned to cpus", get_proc_pin(), "\n");
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
    int numa_node_i = atoi(entry.c_str() + 4);
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
      numa_node_list[numa_node_i].second.push_back(atoi(cpu_entry.c_str() + 3));
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
  int my_numa_node = upcxx::local_team().rank_me() % numa_nodes_to_use;
  vector<int> my_cpu_list = numa_node_list[my_numa_node].second;
  sort(my_cpu_list.begin(), my_cpu_list.end());
  pin_proc(my_cpu_list);
  SLOG("Pinning to NUMA nodes: process 0 on node 0 is pinned to cpus", get_proc_pin(), "\n");
}
