#ifndef _FASTQ_H
#define _FASTQ_H


#include <upcxx/upcxx.hpp>

using std::string;

using upcxx::rank_me;
using upcxx::rank_n;

#define INT_CEIL(numerator, denominator) (((numerator) - 1) / (denominator) + 1)
#define BUF_SIZE 2047


class FastqReader {

private:
  FILE *f;
  off_t file_size;
  int64_t start_read;
  int64_t end_read;
  string fname;
  int max_read_len;
  char buf[BUF_SIZE + 1];

  void rtrim(string &s) {
    auto pos = s.length() - 1;
    while (pos >= 0 && std::isspace(static_cast<unsigned char>(s[pos]))) pos--;
    s.resize(pos + 1);
  }

  bool get_fq_name(string header) {
    if (header[0] != '@') return false;
    // convert if new illumina 2 format  or HudsonAlpha format
    int len = header.length();
    if (header[len - 2] != '/') {
      if (header[len - 2] == 'R') {
        // HudsonAlpha format  (@pair-R1,  @pair-R2)
        char *end = &header[len - 3];
        // replace with @pair/1 or @pair/2 */
        char rNum = header[len - 1];
        header[len - 3] = '/';
        header[len - 2] = rNum;
        header.resize(len - 1);
        return true;
      } else {
        // Latest Illumina header format
        auto end_pos = header.find('\t');
        if (end_pos != string::npos) {
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
          return false;
        }
        // @pair 1:Y:.:.: or @pair 2:Y:.:.:
        // replace with @pair/1 or @pair/2
        header[end_pos] = '/';
        header.resize(end_pos + 2);
      }
    }
    assert(header[0] == '@');
    return true;
  }
  
  int64_t get_fptr_for_next_record(int64_t offset) {
    // first record is the first record, include it.  Every other partition will be at least 1 full record after offset.
    if (offset == 0) return 0;
    if (offset >= file_size) return file_size;
    if (fseek(f, offset, SEEK_SET) != 0) DIE("Could not fseek in ", fname, " to ", offset, ": ", strerror(errno));
    // skip first (likely partial) line after this offset to ensure we start at the beginning of a line
    if (!fgets(buf, BUF_SIZE, f)) return ftell(f);

    for (int i = 0; ; i++) {
      if (!fgets(buf, BUF_SIZE, f)) return ftell(f);
      // keep reading new lines until we find a valid header ending in 2 for the second of the pair
      if (buf[0] == '@') {
        string header(buf);
        rtrim(header);
        if (!get_fq_name(header)) continue;
        // need this to be the second of the pair
        if (header[header.length() - 1] == '2') {
          // now read another three lines and be done
          for (int j = 0; j < 3; j++) {
            if (!fgets(buf, BUF_SIZE, f)) DIE("Missing record info at pos ", ftell(f));
          }
          break;
        }
      }
      if (i > 13) DIE("Could not find a valid line in the fastq file ", fname, ", last line: ", buf);
    }
    return ftell(f);
  }
  
public:
  
  FastqReader(const string &fname) : fname(fname), max_read_len(0) {
    // only one rank gets the file size, to prevent many hits on metadata
    if (!rank_me()) {
      struct stat st;
      if (stat(fname.c_str(), &st) != 0) DIE("Could not stat ", fname, ": ", strerror(errno));
      file_size = st.st_size;
    }
    file_size = upcxx::broadcast(file_size, 0).wait();
    f = fopen(fname.c_str(), "r");
    if (!f) DIE("Could not open file ", fname, ": ", strerror(errno));
    // just a part of the file is read by this thread 
    int64_t read_block = INT_CEIL(file_size, rank_n());
    start_read = read_block * rank_me();
    end_read = read_block * (rank_me() + 1);
    if (rank_me()) start_read = get_fptr_for_next_record(start_read);
    if (rank_me() == rank_n() - 1) end_read = file_size;
    else end_read = get_fptr_for_next_record(end_read);
    if (fseek(f, start_read, SEEK_SET) != 0) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
    posix_fadvise(fileno(f), start_read, end_read - start_read, POSIX_FADV_SEQUENTIAL);
    SOUT("Reading FASTQ file ", fname, "\n");
    DBG("Reading fastq file ", fname, " at pos ", start_read, " ", ftell(f), "\n");
  }

  ~FastqReader() {
    fclose(f);
  }
  
  bool get_next_fq_record(string &id, string &seq, string &quals) {
    if (ftell(f) >= end_read) return false;
    for (int i = 0; i < 4; i++) {
      if (!fgets(buf, BUF_SIZE, f)) {
        WARN("Read record terminated before full record at position ", ftell(f));
        return false;
      }
      if (i == 0) id.assign(buf);
      else if (i == 1) seq.assign(buf);
      else if (i == 3) quals.assign(buf);
    }
    rtrim(id);
    rtrim(seq);
    rtrim(quals);
    if (id[0] != '@') DIE("Invalid FASTQ at position ", ftell(f), " expected read name (@), got: ", id);
    // construct universally formatted name (illumina 1 format)
    if (!get_fq_name(id)) DIE("Invalid FASTQ name format: '", id, "'");
    if (seq.length() != quals.length()) DIE("Sequence length ", seq.length(), " != ", quals.length(), " quals length at pos ", ftell(f));
    if (seq.length() > max_read_len) max_read_len = seq.length();
    //DBG(id, "\n", seq, "\n+\n", quals, "\n");
    return true;
  }
  
};

#endif
