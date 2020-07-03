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

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/timers.hpp"

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

inline int pin_clear() {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  for(int i = 0; i < sizeof(cpu_set_t) * 8; i++) {
    CPU_SET(i, &cpu_set);
  }
  if (sched_setaffinity(getpid(), sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), getpid());
    return -1;
  }
  return 0;
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
  SLOG_VERBOSE("Successfully wrote ", get_size_str(tot_bytes_written), " bytes to ", fname);
  upcxx::barrier();
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
    // FIXME for summit / power9 with no physical or core id field... use lscpu maybe?
    if (cpu2core.empty()) return ret;
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

