#ifndef _FASTQ_H
#define _FASTQ_H

#include <iostream>
// Not available in gcc <= 7
//#include <charconv>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"

using std::string;
using std::string_view;
using std::to_string;

using upcxx::rank_me;
using upcxx::rank_n;

using namespace upcxx_utils;

#define INT_CEIL(numerator, denominator) (((numerator) - 1) / (denominator) + 1)
#define BUF_SIZE 2047

#define PER_RANK_FILE true

struct CachedRead {
  static inline const std::array<char, 5> nucleotide_map = {'A', 'C', 'G', 'T', 'N'};
  // read_id is not packed as it is already reduced to an index number
  // the pair number is indicated in the read id - negative means pair 1, positive means pair 2
  int64_t read_id;
  // each cached read packs the nucleotide into 3 bits (ACGTN), and the quality score into 5 bits
  unsigned char *packed_read;
  // the read is not going to be larger than 65536 in length, but possibly larger than 256
  uint16_t read_len;
  // overall, we expect the compression to be around 50%. E.g. a read of 150bp would be
  // 8+150+2=160 vs 13+300=313

  CachedRead(const string &id_str, string_view seq, string_view quals, int qual_offset) {
    read_id = strtol(id_str.c_str() + 1, nullptr, 10);
//    auto res = std::from_chars(id_str.data() + 2, id_str.data() + id_str.size() - 2, read_id);
//    if (res.ec != std::errc()) DIE("Failed to convert string to int64_t, ", res.ec);
    // this uses from_chars because it's the fastest option out there
    //std::from_chars(id_str.data() + 2, id_str.data() + id_str.size() - 2, read_id);
    // negative if first of the pair
    if (id_str[id_str.length() - 1] == '1') read_id *= -1;
    // packed is same length as sequence. Set first 3 bits to represent A,C,G,T,N
    // set next five bits to represent quality (from 0 to 32). This doesn't cover the full quality range (only up to 32)
    // but it's all we need since once quality is greater than the qual_thres (20), we treat the base as high quality
    packed_read = new unsigned char [seq.length()];
    for (int i = 0; i < seq.length(); i++) {
      switch (seq[i]) {
        case 'A': packed_read[i] = 0; break;
        case 'C': packed_read[i] = 1; break;
        case 'G': packed_read[i] = 2; break;
        case 'T': packed_read[i] = 3; break;
        case 'N': packed_read[i] = 4; break;
      }
      packed_read[i] |= ((unsigned char)std::min(quals[i] - qual_offset, 31) << 3);
    }
    read_len = (uint16_t)seq.length();
  }

  ~CachedRead() {
    delete[] packed_read;
  }

  void unpack(string &read_id_str, string &seq, string &quals, int qual_offset) {
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
};

class FastqReader {
  FILE *f;
  off_t file_size;
  int64_t start_read;
  int64_t end_read;
  string fname;
  int max_read_len;
  char buf[BUF_SIZE + 1];
  int qual_offset;
  vector<std::unique_ptr<CachedRead>> cached_reads;

  int64_t cache_index;
  bool cached;

  IntermittentTimer io_t;
  inline static double overall_io_t = 0;

  void rtrim(string &s) {
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
          //WARN("unknown format ", header, " end pos ", end_pos, "\n");
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

    for (int i = 0; ; i++) {
      if (!fgets(buf, BUF_SIZE, f)) {
        io_t.stop();
        return ftell(f);
      }
      // keep reading new lines until we find a valid header ending in 2 for the second of the pair
      if (buf[0] == '@') {
        string header(buf);
        rtrim(header);
        if (!get_fq_name(header)) continue;
        // need this to be the second of the pair
        if (header[header.length() - 1] == '2') {
          // now read another three lines, but check that the third line is a + separator
          bool record_found = true;
          for (int j = 0; j < 3; j++) {
            if (!fgets(buf, BUF_SIZE, f)) DIE("Missing record info at pos ", ftell(f));
            if (j == 1 && strncmp(buf, "+", strlen(buf) - 1) != 0) {
              //DIE("muddled fastq, header ", header, " buf '", buf, "' len ", strlen(buf));
              // this is not a correct boundary for a fastq record - move on
              record_found = false;
              break;
            }
          }
          if (record_found) break;
        }
      }
      if (i > 20) DIE("Could not find a valid line in the fastq file ", fname, ", last line: ", buf);
    }
    io_t.stop();
    return ftell(f);
  }

public:

  FastqReader(const string &fname)
    : fname(fname)
    , f(nullptr)
    , max_read_len(0)
    , io_t("fastq IO for " + fname)
    , cached(false)
    , cache_index(0) {
    // only one rank gets the file size, to prevent many hits on metadata
    if (!rank_me()) file_size = get_file_size(fname);
    file_size = upcxx::broadcast(file_size, 0).wait();
    f = fopen(fname.c_str(), "r");
    if (!f) {
      if (!upcxx::rank_me()) DIE("Could not open file ", fname, ": ", strerror(errno));
      upcxx::barrier();
    }
    // just a part of the file is read by this thread
    int64_t read_block = INT_CEIL(file_size, rank_n());
    start_read = read_block * rank_me();
    end_read = read_block * (rank_me() + 1);
    if (rank_me()) start_read = get_fptr_for_next_record(start_read);
    if (rank_me() == rank_n() - 1) end_read = file_size;
    else end_read = get_fptr_for_next_record(end_read);
    if (fseek(f, start_read, SEEK_SET) != 0) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
    posix_fadvise(fileno(f), start_read, end_read - start_read, POSIX_FADV_SEQUENTIAL);
    SLOG_VERBOSE("Reading FASTQ file ", fname, "\n");
    DBG("Reading fastq file ", fname, " at pos ", start_read, " ", ftell(f), "\n");
  }

  ~FastqReader() {
    if (f) fclose(f);
    io_t.done();
    FastqReader::overall_io_t += io_t.get_elapsed();
  }

  size_t my_file_size() {
    return end_read - start_read;
  }

  size_t get_next_fq_record(string &id, string &seq, string &quals) {
    if (cached) {
      if (cache_index == cached_reads.size()) return 0;
      cached_reads[cache_index]->unpack(id, seq, quals, qual_offset);
      cache_index++;
      return id.length() + seq.length() + quals.length();
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
    // FIXME: this should be hexify
    replace_spaces(id);
    if (seq.length() != quals.length())
      DIE("Invalid FASTQ in ", fname, ": sequence length ", seq.length(), " != ", quals.length(), " quals length\n",
          "id:   ", id, "\nseq:  ", seq, "\nquals: ", quals);
    if (seq.length() > max_read_len) max_read_len = seq.length();
    io_t.stop();
    return bytes_read;
  }

  int get_max_read_len() {
    return max_read_len;
  }

  double static get_io_time() {
    return overall_io_t;
  }

  void reset() {
    if (cached) {
      cache_index = 0;
    } else {
      if (fseek(f, start_read, SEEK_SET) != 0) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
    }
  }

  void load_cache(int qual_offset_param) {
    Timer timer(__FILEFUNC__);
    qual_offset = qual_offset_param;
    io_t.start();
    // first estimate the number of records
    size_t tot_bytes_read = 0;
    int64_t num_records = 0;
    string id, seq, quals;
    for (num_records = 0; num_records < 10000; num_records++) {
      size_t bytes_read = get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      tot_bytes_read += bytes_read;
    }
    int64_t bytes_per_record = tot_bytes_read / num_records;
    int64_t estimated_records = my_file_size() / bytes_per_record;
    cached_reads.reserve(estimated_records);
    reset();
    ProgressBar progbar(my_file_size(), "Loading reads into cache");
    tot_bytes_read = 0;
    while (true) {
      size_t bytes_read = get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      cached_reads.push_back(std::make_unique<CachedRead>(id, seq, quals, qual_offset));
    }
    progbar.done();
    upcxx::barrier();
    auto all_estimated_records = upcxx::reduce_one(estimated_records, upcxx::op_fast_add, 0).wait();
    auto all_num_records = upcxx::reduce_one(cached_reads.size(), upcxx::op_fast_add, 0).wait();
    SLOG_VERBOSE("Loaded ", all_num_records, " reads (estimated ", all_estimated_records, ")\n");
    reset();
    cached = true;
    io_t.stop();
    fclose(f);
    f = nullptr;
    // FIXME: this gives a messed up elapsed time
    io_t.done();
    FastqReader::overall_io_t += io_t.get_elapsed();
  }

};


#endif
