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

#include <map>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/timers.hpp"

using std::shared_ptr;
using std::string;
using std::to_string;

using upcxx::dist_object;
using upcxx::promise;
using upcxx::rank_me;
using upcxx::rank_n;

using upcxx_utils::IntermittentTimer;

#define INT_CEIL(numerator, denominator) (((numerator)-1) / (denominator) + 1)
#define BUF_SIZE 2047

class FastqReader {
  string fname;
  FILE *f;
  off_t file_size;
  int64_t start_read;
  int64_t end_read;
  unsigned max_read_len;
  char buf[BUF_SIZE + 1];
  int qual_offset;
  shared_ptr<FastqReader> fqr2;
  bool first_file;
  bool _is_paired;
  IntermittentTimer io_t;
  struct PromStartStop {
    promise<int64_t> start_prom, stop_prom;
    upcxx::future<> set(FastqReader &fqr) {
      DBG("Setting fqr ", fqr.fname, " start=", fqr.start_read, " end=", fqr.end_read, "\n");
      auto set_start = start_prom.get_future().then([&fqr](int64_t start) { fqr.start_read = start; });
      auto set_end = stop_prom.get_future().then([&fqr](int64_t stop) { fqr.end_read = stop; });
      return when_all(set_start, set_end);
    }
  };
  dist_object<PromStartStop> dist_prom;
  upcxx::future<> open_fut;
  void seek();

  inline static double overall_io_t = 0;

  static void rtrim(string &s);

  bool get_fq_name(string &header);

  int64_t get_fptr_for_next_record(int64_t offset);

 public:
  FastqReader() = delete;  // no default constructor
  FastqReader(const string &_fname, bool wait = false, upcxx::future<> first_wait = make_future());

  // this happens within a separate thread
  upcxx::future<> continue_open(int fd = -1);

  ~FastqReader();

  string get_fname();

  size_t my_file_size();

  void advise(bool will_need);

  size_t get_next_fq_record(string &id, string &seq, string &quals, bool wait_open = true);
  int get_max_read_len();

  double static get_io_time();

  void reset();

  upcxx::future<> get_open_fut() const { return open_fut; }

  bool is_paired() const { return _is_paired; }

  static upcxx::future<> set_matching_pair(FastqReader &fqr1, FastqReader &fqr2, dist_object<PromStartStop> &dist_start_stop1,
                                    dist_object<PromStartStop> &dist_start_stop2);
};

class FastqReaders {
  // singleton class to hold as set of fastq readers open and re-usable
  using ShFastqReader = shared_ptr<FastqReader>;
  std::unordered_map<string, ShFastqReader> readers;

  FastqReaders();
  ~FastqReaders();

 public:
  static FastqReaders &getInstance();

  static FastqReader &open(const string fname);

  template <typename Container>
  static void open_all(Container &fnames) {
    for (string &fname : fnames) {
      open(fname);
    }
  }

  static FastqReader &get(const string fname);

  static void close(const string fname);

  static void close_all();
};
