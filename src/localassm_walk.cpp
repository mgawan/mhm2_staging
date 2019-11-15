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


struct BaseQual {
  char base;
  int nvotes;
  int nvotes_hi_q;
  int rating;
};


struct MerBase {
  string mer;
  char base;
  int nvotes;
  char fork_bases[4];
  int fork_nvotes[4];
};


struct MerFreqs {
  string mer;
  // freqs for all reasonable quality bases
  int base_freqs[4];
  // freqs for high quality bases
  int base_freqs_hi_q[4];
};


bool get_next_mer_pos(int &start, const string &seq, int mer_len) {
  int end = (start == 0 ? 0 : start + mer_len - 1);
  for (; (end < seq.length()) && (end - start < mer_len); end++) {
    if (seq[end] == 'N') {
      while (end < seq.length() && seq[end] == 'N') end++;
      start = end;
    }
  }
  if (end - start != mer_len) return false;
  return true;
}



static int rate_base_ext(BaseQual &base_qual, int depth) {
  // 0 = No votes
  // 1 = One vote
  // 2 = nVotes < minViable
  // 3 = minDepth > nVotes >= minViable; nHiQ < minViable
  // 4 = minDepth > nVotes >= minViable ; nHiQ >= minViable
  // 5 = nVotes >= minDepth ; nHiQ < minViable
  // 6 = nVotes >= minDepth ; minViable < nHiQ < minDepth
  // 7 = nHiQ >= minDepth
  // make thresholds dependent on contig depth
  double min_viable = max(LASSM_MIN_VIABLE_DEPTH * depth, 2.0);
  double min_expected_depth = max(LASSM_MIN_EXPECTED_DEPTH * depth, 2.0);
  if (base_qual.nvotes == 0) return 0;
  if (base_qual.nvotes == 1) return 1;
  if (base_qual.nvotes < min_viable) return 2;
  if (base_qual.nvotes < min_expected_depth || base_qual.nvotes == min_viable) {
    if (base_qual.nvotes_hi_q < min_viable) return 3;
    else return 4;
  }
  if (base_qual.nvotes_hi_q < min_viable) return 5;
  if (base_qual.nvotes_hi_q < min_expected_depth) return 6;
  return 7;
}


static void compute_mer_freqs(vector<ReadSeq> &reads, unordered_map<string, MerBase> &mers_ht, int seq_depth,
                              int mer_len, int qual_offset) {
  unordered_map<string, MerFreqs> mer_freqs_ht;
  // split reads into kmers and count frequency of medium and high quality extensions
  for (auto &read_seq : reads) {
    int start = 0;
    if (mer_len >= read_seq.seq.length()) continue;
    // iterate through mers that don't contain Ns
    while (get_next_mer_pos(start, read_seq.seq, mer_len)) {
      string mer = read_seq.seq.substr(start, mer_len);
      auto it = mer_freqs_ht.find(mer);
      if (it == mer_freqs_ht.end()) {
        MerFreqs mer_freqs = {.mer = mer, .base_freqs = {0}, .base_freqs_hi_q = {0}};
        mer_freqs_ht.insert({mer, mer_freqs});
      }
      it = mer_freqs_ht.find(mer);
      int extension_pos = start + mer_len;
      if (extension_pos == read_seq.seq.length()) break;
      char extension = read_seq.seq[extension_pos];
      if (extension == 'N') {
        start++;
        continue;
      }
      int bi = 0;
      switch (extension) {
        case 'A': bi = 0; break;
        case 'C': bi = 1; break;
        case 'G': bi = 2; break;
        case 'T': bi = 3; break;
        default:
          DIE("start ", start, " k ", mer_len, " Invalid base at ", extension_pos, " : ", extension, "(", (int)extension, ")",
              "\nread: ", read_seq.seq, "\n");
      }
      int q = read_seq.quals[extension_pos] - qual_offset;
      // ignore low quality bases
      if (q >= LASSM_MIN_QUAL) it->second.base_freqs[bi]++;
      // high quality
      if (q > LASSM_MIN_HI_QUAL) it->second.base_freqs_hi_q[bi]++;
      start++;
    }
  }
  // Full analysis of quality/frequency profile
  const char BASES[4] = { 'A', 'C', 'G', 'T' };
  for (auto &elem : mer_freqs_ht) {
    auto mer_freqs = &elem.second;
    array<BaseQual, 4> base_quals;
    for (int i = 0; i < 4; i++) {
      base_quals[i] = {.base = BASES[i], .nvotes = mer_freqs->base_freqs[i], .nvotes_hi_q = mer_freqs->base_freqs_hi_q[i]};
      base_quals[i].rating = rate_base_ext(base_quals[i], seq_depth);
    }
    sort(base_quals.begin(), base_quals.end(),
         [](const auto &elem1, const auto &elem2) {
           if (elem1.rating > elem2.rating) return 1;
           if (elem1.rating < elem2.rating) return -1;
           if (elem1.nvotes_hi_q > elem2.nvotes_hi_q) return 1;
           if (elem1.nvotes_hi_q < elem2.nvotes_hi_q) return -1;
           if (elem1.nvotes > elem2.nvotes) return 1;
           if (elem1.nvotes < elem2.nvotes) return -1;
           return 0;
         });
    int top_rating = base_quals[0].rating;
    int runner_up = base_quals[1].rating;
    int top_rated_base = base_quals[0].base;
    MerBase mer_base = { .mer = "", .base = ' ', .nvotes = 0, .fork_bases = {0}, .fork_nvotes = {0}};
    // no extension (base = 0) if the runner up is close to the top rating
    // except, if rating is 7 (best quality), then all bases of rating 7 are forks
    if (top_rating > LASSM_RATING_THRES) {         // must have at least minViable bases
      if (top_rating <= 3) {    // must be uncontested
        if (runner_up == 0) mer_base.base = top_rated_base;
      } else if (top_rating < 6) {
        if (runner_up < 3) mer_base.base = top_rated_base;
      } else if (top_rating == 6) {  // viable and fair hiQ support
        if (runner_up < 4) mer_base.base = top_rated_base;
      } else {                     // strongest rating trumps
        if (runner_up < 7) {
          mer_base.base = top_rated_base;
        } else {
          // runner up is also 7, at least 2-way fork
          int k = 0;
          for (int b = 0; b < 4; b++) {
            if (base_quals[b].rating == 7) {
              mer_base.fork_bases[k] = base_quals[b].base;
              mer_base.fork_nvotes[k] = base_quals[b].nvotes;
              k++;
            } else {
              break;
            }
          }
        }
      }
    }
    if (mer_base.base == top_rated_base) mer_base.nvotes = base_quals[0].nvotes;
    // only put in the hash table if it has a valid extension or fork
    if (mer_base.base || mer_base.fork_bases[0]) {
      mer_base.mer = mer_freqs->mer;
      if (mers_ht.find(mer_base.mer) != mers_ht.end()) DIE("Found previous kmer in mers_ht\n");
      mers_ht.insert({mer_base.mer, mer_base});
    }
  }
}


static string iterative_walks(string &seq, int seq_depth, vector<ReadSeq> &reads, string dirn,
                              int max_mer_len, int kmer_len, int qual_offset) {
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
  for (int mer_len = kmer_len; mer_len >= min_mer_len && mer_len <= max_mer_len; mer_len += shift) {
    unordered_map<string, MerBase> mers_ht;
    compute_mer_freqs(reads, mers_ht, seq_depth, mer_len, qual_offset);
    // make sure the walk starts with the longest possible kmer
    string walk = seq.substr(seq.length() - mer_len);
    int depth = 0;
    char walk_result = walk_mers(mers_ht, walk, &depth, mer_len);
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

