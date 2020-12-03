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

#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <upcxx/upcxx.hpp>

using std::max;
using std::string;
using std::string_view;
using std::to_string;
using std::unique_ptr;
using std::vector;

class PackedRead {
  static inline const std::array<char, 5> nucleotide_map = {'A', 'C', 'G', 'T', 'N'};
  // read_id is not packed as it is already reduced to an index number
  // the pair number is indicated in the read id - negative means pair 1, positive means pair 2
  // int64_t read_id : 48;
  int64_t read_id;
  // the read is not going to be larger than 65536 in length, but possibly larger than 256
  // uint64_t read_len : 16;
  uint16_t read_len;
  // each cached read packs the nucleotide into 3 bits (ACGTN), and the quality score into 5 bits
  unsigned char *bytes;

  // overall, we expect the compression to be around 50%. E.g. a read of 150bp would be
  // 16+150=166 vs 13+300=313

 public:
  // default, move and (deep) copy constructors to support containers without unique_ptr overhead
  PackedRead();
  PackedRead(const string &id_str, string_view seq, string_view quals, int qual_offset);
  PackedRead(const PackedRead &copy);
  PackedRead(PackedRead &&Move);
  PackedRead &operator=(const PackedRead &copy);
  PackedRead &operator=(PackedRead &&move);

  ~PackedRead();

  void unpack(string &read_id_str, string &seq, string &quals, int qual_offset);
  void clear();

  int64_t get_id();
  string get_str_id();
  static int64_t to_packed_id(const string &id_str);

  uint16_t get_read_len();

  struct upcxx_serialization {
    template <typename Writer>
    static void serialize(Writer &writer, PackedRead const &packed_read) {
      writer.write(packed_read.read_id);
      writer.write(packed_read.read_len);
      for (int i = 0; i < packed_read.read_len; i++) {
        writer.write(packed_read.bytes[i]);
      }
    }

    template <typename Reader>
    static PackedRead *deserialize(Reader &reader, void *storage) {
      PackedRead *packed_read = new (storage) PackedRead();
      packed_read->read_id = reader.template read<int64_t>();
      packed_read->read_len = reader.template read<uint16_t>();
      packed_read->bytes = new unsigned char[packed_read->read_len];
      for (int i = 0; i < packed_read->read_len; i++) {
        packed_read->bytes[i] = reader.template read<unsigned char>();
      }
      return packed_read;
    }
  };
};

class PackedReads {
  vector<PackedRead> packed_reads;
  // this is only used when we need to know the actual name of the original reads
  vector<string> read_id_idx_to_str;
  unsigned max_read_len = 0;
  uint64_t index = 0;
  uint64_t bases = 0;
  uint64_t name_bytes = 0;
  int qual_offset;
  string fname;
  bool str_ids;

 public:
  using PackedReadsList = vector<PackedReads *>;
  PackedReads(int qual_offset, const string &fname, bool str_ids = false);
  PackedReads(int qual_offset, vector<PackedRead> &packed_reads);
  ~PackedReads();

  bool get_next_read(string &id, string &seq, string &quals);

  PackedRead &operator[](int index);

  void reset();

  void clear();

  string get_fname();

  unsigned get_max_read_len();

  int64_t get_local_num_reads();

  void add_read(const string &read_id, const string &seq, const string &quals);

  void load_reads();

  upcxx::future<> load_reads_nb();

  static void load_reads(PackedReadsList &);

  void report_size();
};
