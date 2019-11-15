#include <iostream>
#include <fstream>
#include <regex>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "progressbar.hpp"
#include "zstr.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "fastq.hpp"
#include "localassm.hpp"

using namespace std;
using namespace upcxx;


struct MerFreqs {
  string mer;
  // freqs for all reasonable quality bases
  int base_freqs[4];
  // freqs for high quality bases
  int base_freqs_hi_q[4];
};


// try to get substring without N's that's at least max_len long
static bool get_valid_base_offset(string &s, int subseq_pos, int max_len, int &start, int &len) {
  start = 0;
  int first_N = 0;
  for (int j = subseq_pos; j < len; j++) {
    if (s[j] == 'N') {
      first_N = j;
      if (j - start >= max_len) {
        len = j - start;
        break;
      } else {
        while (s[j] == 'N') j++;
        if (j < len - 1) start = j;
        else len = first_N - start;
      }
    }
    if (j == len - 1) len -= start;
  }
  if (len >= max_len) return true;
  else return false;
}


// split reads into kmers and count frequency of medium and high quality extensions
static void compute_mer_freqs(int mer_len, vector<ReadSeq> &reads, unordered_map<string, MerFreqs> &mer_freqs_ht, int qual_offset) {
  for (auto &read_seq : reads) {
    int start = 0;
    if (mer_len >= read_seq.seq.length()) continue;
    int subseq_len = read_seq.seq.length();
    int subseq_pos = 0;
    // in the while loop, we split the read into substrings without Ns that
    // are at least mer_len long
    while (get_valid_base_offset(read_seq.seq, subseq_pos, mer_len + 1, start, subseq_len)) {
      string subseq_nts = read_seq.seq.substr(subseq_pos + start);
      string subseq_quals = read_seq.quals.substr(subseq_pos + start);
      // iterate through all mers in subsequence, finding quality of extension for each mer
      for (int j = 0; j < subseq_len - mer_len; j++) {
        string mer = subseq_nts.substr(j, mer_len);
        auto it = mer_freqs_ht.find(mer);
        if (it == mer_freqs_ht.end()) {
          MerFreqs mer_freqs = {.mer = mer };
          mer_freqs_ht.insert({mer, mer_freqs});
        }
        it = mer_freqs_ht.find(mer);
        int offs = j + mer_len;
        char extension = subseq_nts[offs];
        int bi = 0;
        switch (extension) {
          case 'A': bi = 0; break;
          case 'C': bi = 1; break;
          case 'G': bi = 2; break;
          case 'T': bi = 3; break;
          default:
            DIE("j ", j, " k ", mer_len, " sublen ", subseq_len, " Invalid base at ", offs, " : ", extension,
                "\nread: ", read_seq.seq, "\nsubseq: ", subseq_nts, "\n");
        }
        int q = subseq_quals[offs] - qual_offset;
        // ignore low quality bases
        if (q >= LASSM_MIN_QUAL) it->second.base_freqs[bi]++;
        // high quality
        if (q > LASSM_MIN_HI_QUAL) it->second.base_freqs_hi_q[bi]++;
        subseq_pos += (start + subseq_len);
        if (subseq_pos >= read_seq.seq.length()) break;
        subseq_len = read_seq.seq.length() - subseq_pos;
      }
    }
  }
}


static string iterative_walks(string &seq, int seq_depth, vector<ReadSeq> &reads, string dirn,
                              int max_mer_len, int kmer_len, int qual_offset) {
  Timer timer(__func__);
  int min_mer_len = LASSM_MIN_KMER_LEN;
  max_mer_len = min(max_mer_len, (int)seq.length());
  // iteratively walk starting from kmer_size, increasing mer size on a fork (F) or repeat (R),
  // and decreasing on an end of path (X)
  // look for the longest walk. Note we have to restart from the beginning for each walk to ensure
  // that all loops will be detected
  int longest_walk_len = 0;
  int longest_walk_depth = 0;
  string longest_walk = "";
  char longest_walk_result = 'X';
  int shift = 0;
  unordered_map<string, MerFreqs> mer_freqs_ht;
  //unordered_map<> mers_ht;
  for (int mer_len = kmer_len; mer_len >= min_mer_len && mer_len <= max_mer_len; mer_len += shift) {
    compute_mer_freqs(mer_len, reads, mer_freqs_ht, qual_offset);
    //analyze_mer_freqs(mer_freqs_ht, mers_ht, seq_depth);
    mer_freqs_ht.clear();
    // make sure the walk starts with the longest possible kmer
    string walk = seq.substr(seq.length() - mer_len);
    int depth = 0;
    //char walk_result = walk_mers(mers_ht, walk, &depth, mer_len);
    char walk_result = 'X';
    //mers_ht.clear();
    int walk_len = walk.length() - mer_len;
    if (walk_len > longest_walk_len) {
      longest_walk = walk.substr(mer_len);
      longest_walk_len = walk_len;
      longest_walk_depth = depth;
      longest_walk_result = walk_result;
    }
    if (walk_result == 'X') {
      // walk reaches a dead-end, downshift, unless we were upshifting
      if (shift == LASSM_SHIFT_SIZE) break;
      shift = -LASSM_SHIFT_SIZE;
    } else {
      // otherwise walk must end with a fork or repeat, so upshift
      if (shift == -LASSM_SHIFT_SIZE) break;
      if (mer_len > seq.length()) break;
      shift = LASSM_SHIFT_SIZE;
    }
  }
  return longest_walk;
}


void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, int qual_offset) {
  Timer timer(__func__, true);
  //_num_forks = 0;
  // _num_terms = 0;
  //_num_repeats = 0;
  // walk should never be more than this. Note we use the maximum insert size from all libraries
  int walk_len_limit = insert_avg + 2 * insert_stddev;
  int64_t num_walks = 0, sum_clen = 0, sum_ext = 0, max_walk_len = 0, num_reads = 0;
  ProgressBar progbar(ctgs_dht.get_local_num_ctgs(), "Extending contigs");
  for (auto ctg = ctgs_dht.get_first_local_ctg(); ctg != nullptr; ctg = ctgs_dht.get_next_local_ctg()) {
    progbar.update();
    sum_clen += ctg->seq.length();
    num_reads += ctg->reads_start.size() + ctg->reads_end.size();
    // have to do right first because the contig needs to be revcomped for the left
    string right_walk = iterative_walks(ctg->seq, ctg->depth, ctg->reads_start, "left", max_kmer_len, kmer_len, qual_offset);
    string seq_rc = revcomp(ctg->seq);
    string left_walk = iterative_walks(seq_rc, ctg->depth, ctg->reads_start, "left", max_kmer_len, kmer_len, qual_offset);
    num_walks += !right_walk.empty() + !left_walk.empty();
    if (max_walk_len < right_walk.length()) max_walk_len = right_walk.length();
    if (max_walk_len < left_walk.length()) max_walk_len = left_walk.length();
    sum_ext += right_walk.length() + left_walk.length();
    // extend the seq
    ctg->seq.insert(0, left_walk);
    ctg->seq += right_walk;
  }
  progbar.done();
  barrier();
  auto tot_num_walks = reduce_one(num_walks, op_fast_add, 0).wait();
  auto tot_sum_ext = reduce_one(sum_ext, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", tot_num_walks, " walks, total extension ",
               perc_str(tot_sum_ext, reduce_one(sum_clen, op_fast_add, 0).wait()), "\n");
  if (tot_num_walks) 
    SLOG_VERBOSE("Average walk length ", tot_sum_ext / tot_num_walks, ", max walk length ",
                 reduce_one(max_walk_len, op_fast_max, 0).wait(), "\n");
}

