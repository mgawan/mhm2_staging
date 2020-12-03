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
#include "packed_reads.hpp"

#include <iostream>
// Not available in gcc <= 7
//#include <charconv>
#include <fcntl.h>
#include <unistd.h>

#include <upcxx/upcxx.hpp>

#include "fastq.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "zstr.hpp"

using std::string;
using std::string_view;
using std::to_string;

using upcxx::rank_me;
using upcxx::rank_n;

using namespace upcxx_utils;

using std::max;

PackedRead::PackedRead()
    : read_id(0)
    , read_len(0)
    , packed_read(nullptr) {}

PackedRead::PackedRead(const string &id_str, string_view seq, string_view quals, int qual_offset) {
  read_id = strtol(id_str.c_str() + 1, nullptr, 10) + 1;
  // this uses from_chars because it's the fastest option out there
  // auto res = std::from_chars(id_str.data() + 2, id_str.data() + id_str.size() - 2, read_id);
  // if (res.ec != std::errc()) DIE("Failed to convert string to int64_t, ", res.ec);
  // negative if first of the pair
  if (id_str[id_str.length() - 1] == '1') read_id *= -1;
  // packed is same length as sequence. Set first 3 bits to represent A,C,G,T,N
  // set next five bits to represent quality (from 0 to 32). This doesn't cover the full quality range (only up to 32)
  // but it's all we need since once quality is greater than the qual_thres (20), we treat the base as high quality
  packed_read = new unsigned char[seq.length()];
  for (unsigned i = 0; i < seq.length(); i++) {
    switch (seq[i]) {
      case 'A': packed_read[i] = 0; break;
      case 'C': packed_read[i] = 1; break;
      case 'G': packed_read[i] = 2; break;
      case 'T': packed_read[i] = 3; break;
      case 'N': packed_read[i] = 4; break;
      case 'U':
      case 'R':
      case 'Y':
      case 'K':
      case 'M':
      case 'S':
      case 'W':
      case 'B':
      case 'D':
      case 'H':
      case 'V': packed_read[i] = 4; break;
      default: DIE("Illegal char in comp nucleotide of '", seq[i], "'\n");
    }
    packed_read[i] |= ((unsigned char)std::min(quals[i] - qual_offset, 31) << 3);
  }
  read_len = (uint16_t)seq.length();
}

PackedRead::PackedRead(const PackedRead &copy)
    : read_id(copy.read_id)
    , packed_read(new unsigned char[read_len])
    , read_len(copy.read_len) {
  memcpy(packed_read, copy.packed_read, read_len);
}

PackedRead::PackedRead(PackedRead &&move)
    : read_id(move.read_id)
    , packed_read(move.packed_read)
    , read_len(move.read_len) {
  move.packed_read = nullptr;
  move.clear();
}

PackedRead &PackedRead::operator=(const PackedRead &copy) {
  PackedRead pr(copy);
  std::swap(*this, pr);
  return *this;
}

PackedRead &PackedRead::operator=(PackedRead &&move) {
  PackedRead pr(std::move(move));
  std::swap(*this, pr);
  return *this;
}

PackedRead::~PackedRead() { clear(); }

void PackedRead::clear() {
  if (packed_read) delete[] packed_read;
  packed_read = nullptr;
  read_len = 0;
  read_id = 0;
}

void PackedRead::unpack(string &read_id_str, string &seq, string &quals, int qual_offset) const {
  char pair_id = (read_id < 0 ? '1' : '2');
  read_id_str = "@r" + to_string(labs(read_id)) + '/' + pair_id;
  seq.resize(read_len);
  quals.resize(read_len);
  for (int i = 0; i < read_len; i++) {
    seq[i] = nucleotide_map[packed_read[i] & 7];
    quals[i] = qual_offset + (packed_read[i] >> 3);
  }
  assert(seq.length() == read_len);
  assert(quals.length() == read_len);
}

PackedReads::PackedReads(int qual_offset, const string &fname, bool str_ids)
    : qual_offset(qual_offset)
    , fname(fname)
    , str_ids(str_ids) {}

PackedReads::~PackedReads() { clear(); }

bool PackedReads::get_next_read(string &id, string &seq, string &quals) {
  if (index == packed_reads.size()) return false;
  packed_reads[index].unpack(id, seq, quals, qual_offset);
  if (str_ids) id = read_id_idx_to_str[index];
  index++;
  return true;
}

uint64_t PackedReads::get_read_index() const { return index; }

void PackedReads::get_read(uint64_t index, string &id, string &seq, string &quals) const {
  if (index >= packed_reads.size()) DIE("Invalid get_read(", index, ") - size=", packed_reads.size());
  packed_reads[index].unpack(id, seq, quals, qual_offset);
  if (str_ids) id = read_id_idx_to_str[index];
}

void PackedReads::reset() { index = 0; }

void PackedReads::clear() {
  LOG_MEM("Clearing Packed Reads");
  index = 0;
  fname.clear();
  vector<PackedRead>().swap(packed_reads);
  if (str_ids) vector<string>().swap(read_id_idx_to_str);
  LOG_MEM("Cleared Packed Reads");
}

string PackedReads::get_fname() { return fname; }

unsigned PackedReads::get_max_read_len() { return max_read_len; }

int64_t PackedReads::get_local_num_reads() { return packed_reads.size(); }

void PackedReads::add_read(const string &read_id, const string &seq, const string &quals) {
  packed_reads.emplace_back(read_id, seq, quals, qual_offset);
  if (str_ids) {
    read_id_idx_to_str.push_back(read_id);
    name_bytes += sizeof(string) + read_id.size();
  }
  max_read_len = max(max_read_len, (unsigned)seq.length());
  bases += seq.length();
}

void PackedReads::load_reads(PackedReadsList &packed_reads_list) {
  BarrierTimer timer(__FILEFUNC__);
  upcxx::future<> all_done = upcxx::make_future();
  for (auto pr : packed_reads_list) {
    FastqReaders::open(pr->fname);
  }
  for (auto pr : packed_reads_list) {
    upcxx::discharge();
    upcxx::progress();
    auto fut = pr->load_reads_nb();
    all_done = when_all(all_done, fut);
  }
  FastqReaders::close_all();
  all_done.wait();
}

upcxx::future<> PackedReads::load_reads_nb() {
  // first estimate the number of records
  size_t tot_bytes_read = 0;
  int64_t num_records = 0;
  FastqReader &fqr = FastqReaders::get(fname);
  fqr.advise(true);
  string id, seq, quals;
  for (num_records = 0; num_records < 20000; num_records++) {
    size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
    if (!bytes_read) break;
    tot_bytes_read += bytes_read;
  }
  int64_t bytes_per_record = tot_bytes_read / num_records;
  int64_t estimated_records = fqr.my_file_size() / bytes_per_record;
  int64_t reserve_records = estimated_records * 1.10 + 10000;  // reserve more so there is not a big reallocation if it is under
  packed_reads.reserve(reserve_records);
  fqr.reset();
  ProgressBar progbar(fqr.my_file_size(), "Loading reads from " + fname + " " + get_size_str(fqr.my_file_size()));
  tot_bytes_read = 0;
  while (true) {
    size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
    if (!bytes_read) break;
    tot_bytes_read += bytes_read;
    progbar.update(tot_bytes_read);
    add_read(id, seq, quals);
  }
  fqr.advise(false);
  FastqReaders::close(fname);
  auto fut = progbar.set_done();
  int64_t underestimate = estimated_records - packed_reads.size();
  if (underestimate < 0 && reserve_records < packed_reads.size())
    LOG("NOTICE Underestimated by ", -underestimate, " estimated ", estimated_records, " found ", packed_reads.size(), "\n");
  auto all_under_estimated_fut = upcxx::reduce_one(underestimate < 0 ? 1 : 0, upcxx::op_fast_add, 0);
  auto all_estimated_records_fut = upcxx::reduce_one(estimated_records, upcxx::op_fast_add, 0);
  auto all_num_records_fut = upcxx::reduce_one(packed_reads.size(), upcxx::op_fast_add, 0);
  auto all_num_bases_fut = upcxx::reduce_one(bases, upcxx::op_fast_add, 0);
  return when_all(fut, all_under_estimated_fut, all_estimated_records_fut, all_num_records_fut, all_num_bases_fut)
      .then([max_read_len = this->max_read_len](int64_t all_under_estimated, int64_t all_estimated_records, int64_t all_num_records,
                                                int64_t all_num_bases) {
        SLOG_VERBOSE("Loaded ", all_num_records, " reads (estimated ", all_estimated_records, " with ", all_under_estimated,
                     " ranks underestimated) max_read=", max_read_len, " tot_bases=", all_num_bases, "\n");
      });
}

void PackedReads::load_reads() {
  BarrierTimer timer(__FILEFUNC__);
  load_reads_nb().wait();
  upcxx::barrier();
}

void PackedReads::report_size() {
  auto all_num_records = upcxx::reduce_one(packed_reads.size(), upcxx::op_fast_add, 0).wait();
  auto all_num_bases = upcxx::reduce_one(bases, upcxx::op_fast_add, 0).wait();
  auto all_num_names = upcxx::reduce_one(name_bytes, upcxx::op_fast_add, 0).wait();
  SLOG_VERBOSE("Loaded ", all_num_records, " tot_bases=", all_num_bases, " names=", get_size_str(all_num_names), "\n");
  LOG_MEM("Loaded Packed Reads");
  SLOG_VERBOSE("Estimated memory for PackedReads: ",
               get_size_str(all_num_records * sizeof(PackedRead) + all_num_bases + all_num_names), "\n");
}