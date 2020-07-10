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



#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <memory>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"

using std::string;
using std::to_string;
using std::shared_ptr;

using upcxx::rank_me;
using upcxx::rank_n;

using namespace upcxx_utils;

#define INT_CEIL(numerator, denominator) (((numerator) - 1) / (denominator) + 1)
#define BUF_SIZE 2047


class FastqReader {
  FILE *f;
  off_t file_size;
  int64_t start_read;
  int64_t end_read;
  string fname;
  int max_read_len;
  char buf[BUF_SIZE + 1];
  int qual_offset;
  shared_ptr<FastqReader> fqr2;
  bool first_file;

  IntermittentTimer io_t;
  inline static double overall_io_t = 0;

  static void rtrim(string &s) {
    auto pos = s.length() - 1;
    while (pos >= 0 && std::isspace(static_cast<unsigned char>(s[pos]))) pos--;
    s.resize(pos + 1);
  }

  bool get_fq_name(string &header) {
    if (header[0] != '@') {
      //WARN("unknown format ", header, " missing @\n");
      return false;
    }
    // trim off '@'
    header.erase(0, 1);
    // strip trailing spaces
    //header = std::regex_replace(header, std::regex("\\s+$"), std::string(""));
    rtrim(header);
    // convert if new illumina 2 format  or HudsonAlpha format
    int len = header.length();
    if (header[len - 2] != '/') {
      if (header[len - 2] == 'R') {
        // HudsonAlpha format  (@pair-R1,  @pair-R2)
        // replace with @pair/1 or @pair/2 */
        char rNum = header[len - 1];
        header[len - 3] = '/';
        header[len - 2] = rNum;
        header.resize(len - 1);
        return true;
      } else {
        // Latest Illumina header format
        auto end_pos = header.find('\t');
        if (end_pos == string::npos) {
          end_pos = header.find(' ');
          // no comment, just return without modification
          if (end_pos == string::npos) return true;
        }
        // check for @pair/1 @/pair2 vs @pair 1:Y:... @pair 2:Y:....
        if (end_pos > 3 && header[end_pos - 2] == '/' && (header[end_pos - 1] == '1' || header[end_pos - 1] == '2')) {
          // @pair/1 or @pair/2
          // truncate name at comment
          header.resize(end_pos);
          return true;
        }
        if ((len < end_pos + 7) || header[end_pos + 2] != ':' || header[end_pos + 4] != ':' || header[end_pos + 6] != ':' ||
            (header[end_pos + 1] != '1' && header[end_pos + 1] != '2')) {
          // unknown pairing format
          SWARN("unknown format ", header, " end pos ", end_pos, "\n");
          return false;
        }
        // @pair 1:Y:.:.: or @pair 2:Y:.:.:
        // replace with @pair/1 or @pair/2
        header[end_pos] = '/';
        header.resize(end_pos + 2);
      }
    }
    return true;
  }

  int64_t get_fptr_for_next_record(int64_t offset) {
    io_t.start();
    // first record is the first record, include it.  Every other partition will be at least 1 full record after offset.
    if (offset == 0) return 0;
    if (offset >= file_size) return file_size;
    if (fseek(f, offset, SEEK_SET) != 0) DIE("Could not fseek in ", fname, " to ", offset, ": ", strerror(errno));
    // skip first (likely partial) line after this offset to ensure we start at the beginning of a line
    if (!fgets(buf, BUF_SIZE, f)) {
      io_t.stop();
      return ftell(f);
    }

    char last_pair = '\0', this_pair = '\1';
    int64_t last_tell = ftell(f);
    for (int i = 0; ; i++) {
      int64_t this_tell = ftell(f);
      if (!fgets(buf, BUF_SIZE, f)) {
        io_t.stop();
        return ftell(f);
      }

      if (buf[0] == '@') {
        string header(buf);
        rtrim(header);
        if (!get_fq_name(header)) continue;

        // now read another three lines, but check that the third line is a + separator
        int64_t test_tell = ftell(f);
        bool record_found = true;
        for (int j = 0; j < 3; j++) {
          if (!fgets(buf, BUF_SIZE, f)) DIE("Missing record info at pos ", ftell(f));
          if (j == 0 && buf[0] == '@') {
              // previous was quality line, this next line is the header line.
              record_found = false;
              break;
          }
          if (j == 1 && buf[0] != '+') {
            //DIE("muddled fastq, header ", header, " buf '", buf, "' len ", strlen(buf));
            // this is not a correct boundary for a fastq record - move on
            record_found = false;
            break;
          }
        }
        if (record_found) {
            // good
            i+=3;
        } else {
            // rewind and test next line as a potential header
            DBG("Did not find proper pair, rewinding\n");
            if (fseek(f, test_tell, SEEK_SET) != 0) DIE("Could not fseek in ", fname, " to ", offset, ": ", strerror(errno));
        }

        // need this to be the second of the pair
        this_pair = header[header.length() - 1];
        DBG("Found possible pair ", (char) this_pair, ", header: ", header, "\n");
        if (last_pair == this_pair) {
            // one of two split files or a merged file without pairs
            DBG("Same pair, so keep the first\n");
            break;
        } else if (last_pair == '1' && this_pair == '2') {
            // proper pair
            DBG("Found proper pair 1&2\n");
            break;
        }

        // did not find valid new start of (possibly paired) record
        last_tell = this_tell;
        last_pair = this_pair;
      }
      if (i > 20) DIE("Could not find a valid line in the fastq file ", fname, ", last line: ", buf);
    }
    io_t.stop();
    DBG("Found record at ", last_tell, " after offset ", offset, "\n");
    return last_tell;
  }

public:

  FastqReader(const string &_fname)
    : fname(_fname)
    , f(nullptr)
    , max_read_len(0)
    , fqr2(nullptr)
    , first_file(true)
    , io_t("fastq IO for " + fname) {
    string fname2;
    size_t pos;
    if ((pos = fname.find(':')) != string::npos) {
        // colon separating a pair into two files
        fname2 = fname.substr(pos+1);
        fname = fname.substr(0,pos);
    }
    // only one rank gets the file size, to prevent many hits on metadata
    if (!rank_me()) file_size = get_file_size(fname);
    file_size = upcxx::broadcast(file_size, 0).wait();
    f = fopen(fname.c_str(), "r");
    if (!f) {
      SDIE("Could not open file ", fname, ": ", strerror(errno));
    }

    // just a part of the file is read by this thread
    int64_t read_block = INT_CEIL(file_size, rank_n());
    start_read = read_block * rank_me();
    end_read = read_block * (rank_me() + 1);
    if (rank_me()) start_read = get_fptr_for_next_record(start_read);
    if (rank_me() == rank_n() - 1) end_read = file_size;
    else end_read = get_fptr_for_next_record(end_read);
    if (fseek(f, start_read, SEEK_SET) != 0) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
    // tell the OS this file will be accessed sequentially
    posix_fadvise(fileno(f), start_read, end_read - start_read, POSIX_FADV_SEQUENTIAL);
    SLOG_VERBOSE("Reading FASTQ file ", fname, "\n");
    DBG("Reading fastq file ", fname, " at pos ", start_read, " ", ftell(f), "\n");

    if (!fname2.empty()) {
        // this second reader is generally hidden from the user
        LOG("Opening second FastqReader with ", fname2, "\n");
        fqr2 = make_shared<FastqReader>(fname2);
    }
  }

  ~FastqReader() {
    if (f) fclose(f);
    io_t.done();
    FastqReader::overall_io_t += io_t.get_elapsed();
  }

  string get_fname() {
    return fname;
  }

  size_t my_file_size() {
    return end_read - start_read + (fqr2 ? fqr2->my_file_size() : 0);
  }

  size_t get_next_fq_record(string &id, string &seq, string &quals) {
    if (fqr2) {
      // return a single interleaved file
      if (first_file) {
        first_file = false;
      } else {
        first_file = true;
        return fqr2->get_next_fq_record(id, seq, quals);
      }
    }
    if (feof(f) || ftell(f) >= end_read) return 0;
    io_t.start();
    size_t bytes_read = 0;
    id = "";
    for (int i = 0; i < 4; i++) {
      char *bytes = fgets(buf, BUF_SIZE, f);
      if (!bytes) DIE("Read record terminated on file ", fname, " before full record at position ", ftell(f));
      if (i == 0) id.assign(buf);
      else if (i == 1) seq.assign(buf);
      else if (i == 3) quals.assign(buf);
      bytes_read += strlen(buf);
    }
    rtrim(id);
    rtrim(seq);
    rtrim(quals);
    if (id[0] != '@')
      DIE("Invalid FASTQ in ", fname, ": expected read name (@), got: ", id);
    // construct universally formatted name (illumina 1 format)
    if (!get_fq_name(id))
      DIE("Invalid FASTQ in ", fname, ": incorrect name format '", id, "'");
    // get rid of spaces
    replace_spaces(id);
    if (seq.length() != quals.length())
      DIE("Invalid FASTQ in ", fname, ": sequence length ", seq.length(), " != ", quals.length(), " quals length\n",
          "id:   ", id, "\nseq:  ", seq, "\nquals: ", quals);
    if (seq.length() > max_read_len) max_read_len = seq.length();
    io_t.stop();
    return bytes_read;
  }

  int get_max_read_len() {
    return max(max_read_len, fqr2 ? fqr2->get_max_read_len() : 0);
  }

  double static get_io_time() {
    return overall_io_t;
  }

  void reset() {
    if (fseek(f, start_read, SEEK_SET) != 0) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
    if (fqr2) fqr2->reset();
  }

};


