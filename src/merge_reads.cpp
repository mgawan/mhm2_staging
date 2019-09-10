#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <chrono>
#ifdef __x86_64__
#include <emmintrin.h>
#include <immintrin.h>
#include <x86intrin.h>
#endif
#include <upcxx/upcxx.hpp>

using namespace std;
using namespace upcxx;

#include "zstr.hpp"
#include "utils.hpp"
#include "fastq.hpp"
#include "progressbar.hpp"

static const double Q2Perror[] = { 1.0,
                                   0.7943,	0.6309,	0.5012,	0.3981,	0.3162,
                                   0.2512,	0.1995,	0.1585,	0.1259,	0.1,
                                   0.07943,	0.06310,	0.05012,	0.03981,	0.03162,
                                   0.02512,	0.01995,	0.01585,	0.01259,	0.01,
                                   0.007943,	0.006310,	0.005012,	0.003981,	0.003162,
                                   0.002512,	0.001995,	0.001585,	0.001259,	0.001,
                                   0.0007943,	0.0006310,	0.0005012,	0.0003981,	0.0003162,
                                   0.0002512,	0.0001995,	0.0001585,	0.0001259,	0.0001,
                                   7.943e-05,	6.310e-05,	5.012e-05,	3.981e-05,	3.162e-05,
                                   2.512e-05,	1.995e-05,	1.585e-05,	1.259e-05,	1e-05,
                                   7.943e-06,	6.310e-06,	5.012e-06,	3.981e-06,	3.162e-06,
                                   2.512e-06,	1.995e-06,	1.585e-06,	1.259e-06,	1e-06,
                                   7.943e-07,	6.310e-07,	5.012e-07,	3.981e-07,	3.1622e-07,
                                   2.512e-07,	1.995e-07,	1.585e-07,	1.259e-07,	1e-07,
                                   7.943e-08,	6.310e-08,	5.012e-08,	3.981e-08,	3.1622e-08,
                                   2.512e-08,	1.995e-08,	1.585e-08,	1.259e-08,	1e-08};


// returns the number of mismatches if it is <= max or a number greater than max (but no the actual count)
int16_t fast_count_mismatches(const char *a, const char *b, int len, int16_t max) {
  assert(len < 32768);
  int16_t mismatches = 0;
  int16_t jumpSize, jumpLen;

#ifdef  __x86_64__
  // 128-bit SIMD
  if (len >= 16) {
    jumpSize = sizeof(__m128i);
    jumpLen = len/jumpSize;
    for(int16_t i = 0; i < jumpLen ; i++) {
      __m128i aa = _mm_loadu_si128((const __m128i *)a); // load 16 bytes from a
      __m128i bb = _mm_loadu_si128((const __m128i *)b); // load 16 bytes from b
      __m128i matched = _mm_cmpeq_epi8(aa, bb); // bytes that are equal are now 0xFF, not equal are 0x00
      uint32_t myMaskMatched = _mm_movemask_epi8(matched); // mask of most significant bit for each byte
      // count mismatches
      mismatches += _popcnt32( (~myMaskMatched) & 0xffff ); // over 16 bits
      if (mismatches > max) break;
      a += jumpSize;
      b += jumpSize;
    }
    len -= jumpLen * jumpSize;
  }
#endif
  // CPU version and fall through 8 bytes at a time
  if (mismatches <= max) {
    assert(len >= 0);
    jumpSize = sizeof(int64_t);
    jumpLen = len/jumpSize;
    for(int16_t i = 0; i < jumpLen ; i++) {
      int64_t *aa = (int64_t*) a, *bb = (int64_t*) b;
      if (*aa != *bb) { // likely
        for (int j =0; j<jumpSize; j++) {
          if (a[j] != b[j]) mismatches++;
        }
        if (mismatches > max) break;
      } // else it matched
      a += jumpSize;
      b += jumpSize;
    }
    len -= jumpLen * jumpSize;
  }
  // do the remaining bytes, if needed
  if (mismatches <= max) {
    assert(len >= 0);
    for(int j = 0; j < len; j++) {
      mismatches += ((a[j] == b[j]) ? 0 : 1);
    }
  }
  return mismatches;
}


void merge_reads(vector<string> reads_fname_list, int qual_offset)
{
  Timer timer(__func__);

  int64_t num_ambiguous = 0;
  int64_t num_merged = 0;
  for (auto const &reads_fname : reads_fname_list) {
    FastqReader fqr(reads_fname);
    string out_fname = get_merged_reads_fname(reads_fname);
    ProgressBar progbar(fqr.my_file_size(), "Merging reads " + reads_fname + " " + get_size_str(fqr.my_file_size()));
    zstr::ofstream out_file(out_fname);
    ostringstream out_buf;
    
    int64_t bytes_written = 0;
    int max_read_len = 0;
    int overlap_len = 0;
    int merged_len = 0;
    
    const int16_t MIN_OVERLAP = 12;
    const int16_t EXTRA_TEST_OVERLAP = 2;
    const int16_t MAX_MISMATCHES = 3; // allow up to 3 mismatches, with MAX_PERROR
    const int Q2PerrorSize = sizeof(Q2Perror) / sizeof(*Q2Perror);
    assert(qual_offset == 33 || qual_offset == 64);

    // illumina reads generally accumulate errors at the end, so allow more mismatches in the overlap as long as differential quality indicates a clear winner
    const double MAX_PERROR = 0.025; // allow up to 2.5% accumulated mismatch probablity of error within the overlap by differential quality score
    const int16_t EXTRA_MISMATCHES_PER_1000 = (int) 150; // allow additional mismatches per 1000 bases of overlap before aborting an overlap test
    const uint8_t MAX_MATCH_QUAL = 41 + qual_offset;
    
    string id1, seq1, quals1, id2, seq2, quals2;
    int64_t num_pairs = 0;
    size_t tot_bytes_read = 0;
    for (; ; num_pairs++) {
      size_t bytes_read1 = fqr.get_next_fq_record(id1, seq1, quals1);
      if (!bytes_read1) break;
      size_t bytes_read2 = fqr.get_next_fq_record(id2, seq2, quals2);
      if (!bytes_read2) break;
      tot_bytes_read += bytes_read1 + bytes_read2;
      progbar.update(tot_bytes_read);
      if (id1.substr(0, id1.length() - 2) != id2.substr(0, id2.length() - 2)) DIE("Mismatched pairs ", id1, " ", id2);
      if (id1[id1.length() - 1] != '1' || id2[id2.length() - 1] != '2') DIE("Mismatched pair numbers ", id1, " ", id2);

      bool is_merged = 0;
      int8_t abort_merge = 0;

      // revcomp the second mate pair and reverse the second quals
      string rc_seq2 = revcomp(seq2);
      string rev_quals2 = quals2;
      reverse(rev_quals2.begin(), rev_quals2.end());

      // use start_i to offset inequal lengths which can be very different but still overlap near the end.  250 vs 178..
      int16_t len = (rc_seq2.length() < seq1.length()) ? rc_seq2.length() : seq1.length();
      int16_t start_i = ((len == seq1.length()) ? 0 : seq1.length() - len);
      int16_t found_i = -1;
      int16_t best_i = -1;
      int16_t best_mm = len;
      double best_perror = -1.0;

      // slide along seq1
      for (int16_t i = 0; i < len - MIN_OVERLAP + EXTRA_TEST_OVERLAP; i++) { // test less overlap than MIN_OVERLAP
        if (abort_merge) break;
        int16_t overlap = len-i;
        int16_t this_max_mismatch = MAX_MISMATCHES + (EXTRA_MISMATCHES_PER_1000 * overlap / 1000);
        int16_t error_max_mismatch = this_max_mismatch * 4 / 3 + 1; // 33% higher
        if (fast_count_mismatches(seq1.c_str() + start_i + i, rc_seq2.c_str(), overlap, error_max_mismatch) > error_max_mismatch) continue;
        int16_t matches = 0, mismatches = 0, bothNs = 0, Ncount = 0;
        int16_t overlapChecked = 0;
        double perror = 0.0;
        for (int16_t j = 0; j < overlap; j++) {
          overlapChecked++;
          char ps = seq1[start_i + i + j];
          char rs = rc_seq2[j];
          if (ps == rs) {
            matches++;
            if (ps == 'N') {
              Ncount += 2;
              if (bothNs++) {
                abort_merge++;
                num_ambiguous++;
                break; // do not match multiple Ns in the same position -- 1 is okay
              }
            }
          } else {
            mismatches++;
            if (ps == 'N') {
              mismatches++; // N still counts as a mismatch
              Ncount++;
              quals1[start_i + i + j] = qual_offset;
              assert(rev_quals2[j] - qual_offset < Q2PerrorSize); 
              assert(rev_quals2[j] - qual_offset >= 0); 
              perror += Q2Perror[rev_quals2[j] - qual_offset];
            } else if (rs == 'N') {
              Ncount++;
              mismatches++; // N still counts as a mismatch
              rev_quals2[j] = qual_offset;
              assert(quals1[start_i + i + j] - qual_offset < Q2PerrorSize); 
              assert(quals1[start_i + i + j] - qual_offset >= 0);
              perror += Q2Perror[quals1[start_i + i + j] - qual_offset];
            }
            if (MAX_PERROR > 0.0) {
              assert(quals1[start_i + i + j] >= qual_offset);
              assert(rev_quals2[j] >= qual_offset);
              uint8_t q1 = quals1[start_i + i + j] - qual_offset;
              uint8_t q2 = rev_quals2[j] - qual_offset;
              if (q1 < 0 || q2 < 0 || q1 >= Q2PerrorSize || q2 >= Q2PerrorSize) 
                DIE("Invalid quality score for read ", id1, " '", quals1[start_i + i + j], "' ", id2, " '", rev_quals2[j],
                    "' assuming common qual_offset of ", qual_offset,
                    ". Check the data and make sure it follows a single consistent quality scoring model (phred+64 vs. phred+33)");

              // sum perror as the difference in q score perrors
              uint8_t diffq = (q1 > q2) ? q1-q2 : q2-q1;
              if (diffq <= 2) {
                perror += 0.5; // cap at flipping a coin when both quality scores are close
              } else {
                assert(diffq < Q2PerrorSize);
                perror += Q2Perror[diffq];
              }
            }
          }
          if (Ncount > 3) {
            abort_merge++;
            num_ambiguous++;
            break; // do not match reads with many Ns
          }
          if (mismatches > error_max_mismatch) break;
        }
        int16_t match_thres = overlap - this_max_mismatch;
        if (match_thres < MIN_OVERLAP) match_thres = MIN_OVERLAP;
        if (matches >= match_thres && overlapChecked == overlap && mismatches <= this_max_mismatch && perror/overlap <= MAX_PERROR) {
          if (best_i < 0 && found_i < 0) {
            best_i = i;
            best_mm = mismatches;
            best_perror = perror;
          } else {
            // another good or ambiguous overlap detected
            num_ambiguous++;
            best_i = -1;
            best_mm = len;
            best_perror = -1.0;
            break;
          }
        } else if (overlapChecked == overlap && mismatches <= error_max_mismatch && perror/overlap <= MAX_PERROR * 4 / 3) {
          // lower threshold for detection of an ambigious overlap
          found_i = i;
          if (best_i >= 0) {
            // ambiguous mapping found after a good one was
            num_ambiguous++;
            best_i = -1;
            best_mm = len;
            best_perror = -1.0;
            break;
          }
        }
      }

      if (best_i >= 0 && !abort_merge) {
        int16_t i = best_i;
        int16_t overlap = len-i;
        // pick the base with the highest quality score for the overlapped region
        for(int16_t j = 0 ; j < overlap ; j++) {
          if (seq1[start_i+i+j] == rc_seq2[j]) {
            // match boost quality up to the limit
            uint16_t newQual = quals1[start_i+i+j] + rev_quals2[j] - qual_offset;
            quals1[start_i+i+j] = ((newQual > MAX_MATCH_QUAL) ? MAX_MATCH_QUAL : newQual);
            assert( quals1[start_i+i+j] >= quals1[start_i+i+j] );
            assert( quals1[start_i+i+j] >= rev_quals2[j] );
          } else {
            uint8_t newQual;
            if (quals1[start_i+i+j] < rev_quals2[j]) {
              // use rev base and discount quality
              newQual = rev_quals2[j] - quals1[start_i+i+j] + qual_offset;
              seq1[start_i+i+j] = rc_seq2[j];
            } else {
              // keep prev base, but still discount quality
              newQual = quals1[start_i+i+j] - rev_quals2[j] + qual_offset;
            }
            quals1[start_i+i+j] = ((newQual > (2+qual_offset)) ? newQual : (2+qual_offset)); // a bit better than random chance here
          }
          assert(quals1[start_i+i+j] >= qual_offset);
        }

        // include the remainder of the rc_seq2 and quals
        seq1 = seq1.substr(0, start_i + i + overlap) + rc_seq2.substr(overlap);
        quals1 = quals1.substr(0, start_i + i + overlap) + rev_quals2.substr(overlap);
                       
        is_merged = true;
        num_merged++;
        out_buf << id1 << "\n" << seq1 << "\n+\n" << quals1 << "\n";
        out_buf << id2 << "\nN\n+\n" << (char)qual_offset << "\n";
        int read_len = seq1.length(); // caculate new merged length
        if (max_read_len < read_len) max_read_len = read_len;
        merged_len += read_len;
        overlap_len += overlap;
      }
      if (!is_merged) {
        // undo the revcomp
        out_buf << id1 << "\n" << seq1 << "\n+\n" << quals1 << "\n";
        out_buf << id2 << "\n" << seq2 << "\n+\n" << quals2 << "\n";
      }
      if (!(num_pairs % 1000)) {
        out_file << out_buf.str();
        bytes_written += out_buf.str().length();
        out_buf = ostringstream();
      }
    }
    if (!out_buf.str().empty()) {
      out_file << out_buf.str();
      bytes_written += out_buf.str().length();
    }
    out_file.close();
    progbar.done();
    SOUT("tot bytes read ", tot_bytes_read, "\n");
    // store the uncompressed size in a secondary file
    write_uncompressed_file_size(get_merged_reads_fname(reads_fname) + ".uncompressedSize", bytes_written);
    
    auto all_num_pairs = upcxx::reduce_one(num_pairs, op_fast_add, 0).wait();
    auto all_num_merged = upcxx::reduce_one(num_merged, op_fast_add, 0).wait();
    auto all_num_ambiguous = upcxx::reduce_one(num_ambiguous, op_fast_add, 0).wait();
    auto all_merged_len = upcxx::reduce_one(merged_len, op_fast_add, 0).wait();
    auto all_overlap_len = upcxx::reduce_one(overlap_len, op_fast_add, 0).wait();
    auto all_max_read_len = upcxx::reduce_one(max_read_len, op_fast_max, 0).wait();
    SOUT("Merged reads in file ", reads_fname, ":\n");
    SOUT("  merged ", perc_str(all_num_merged, all_num_pairs), " pairs\n");
    SOUT("  ambiguous ", perc_str(all_num_ambiguous, all_num_pairs), " ambiguous pairs\n");
    SOUT("  average merged length ", (double)all_merged_len / all_num_merged, "\n");
    SOUT("  average overlap length ", (double)all_overlap_len / all_num_merged, "\n");
    SOUT("  max read length ", all_max_read_len, "\n");

    barrier();
  }
}


