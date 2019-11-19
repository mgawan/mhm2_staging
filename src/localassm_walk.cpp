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
#include "kmer_dht.hpp"


using namespace std;
using namespace upcxx;

/*
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
  int max_possible_mers = 0;
  for (auto &read_seq : reads) {
    //int start = 0;
    if (mer_len >= read_seq.seq.length()) continue;
    int num_mers = read_seq.seq.length() - mer_len;
    max_possible_mers += num_mers;
    for (int start = 0; start < num_mers; start++) {
      // skip mers that contain Ns
      if (read_seq.seq.find("N", start) != string::npos) continue;
      int extension_pos = start + mer_len;
      if (extension_pos >= read_seq.seq.length()) DIE("extension_pos ", extension_pos, " out of range ", mer_len);
      char extension = read_seq.seq[extension_pos];
      if (extension == 'N') continue;
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
      
      string mer = read_seq.seq.substr(start, mer_len);

      
      auto it = mer_freqs_ht.find(mer);
      if (it == mer_freqs_ht.end()) {
        mer_freqs_ht.insert({mer, {.mer = mer, .base_freqs = {0}, .base_freqs_hi_q = {0}}});
        it = mer_freqs_ht.find(mer);
      }
    }
  }
  DBG("Compute mer freqs for mer len ", mer_len, ", put ", mer_freqs_ht.size(),
      " mers into the table, max possible ", max_possible_mers, "\n");
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
    MerBase mer_base = { .mer = "", .base = 0, .nvotes = 0, .fork_bases = {0, 0, 0, 0}, .fork_nvotes = {0, 0, 0, 0}};
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
  DBG("Found ", mers_ht.size(), " valid kmers out of ", mer_freqs_ht.size(), " possibilities\n");
}

// return the result of the walk (f, r or x)
static char walk_mers(MerMap &mers_ht, string &mer, string &walk, int &depth, int mer_len, int walk_len_limit) {
  bool have_forked = false;
  int nsteps = 0;
  int tot_ext_votes = 0;
  unordered_map<string, bool> loop_check_ht;
  char walk_result = 'X';
  depth = 0;
  while (1) {
    // check for a cycle in the graph
    if (loop_check_ht.find(mer) != loop_check_ht.end()) {
      walk_result = 'R';
      // if repeated, remove the final mer, since that is what repeated
      // we actually walked one less step
      nsteps--;
      break;
    } else {
      loop_check_ht.insert({mer, true});
    }
    auto it = mers_ht.find(mer);
    if (it == mers_ht.end()) {
      walk_result = 'X';
      // the final base is not part of the walk, since the mer couldn't be found
      //if (nsteps) walk.resize(walk.length() - 1);
      break;
    }
    auto next = &it->second;
    char ext = next->base;
    int nvotes = next->nvotes;
    if (!ext) {
      if (!have_forked && next->fork_bases[0]) {
        // Biallelic positions only (for now) maximum vote path is taken
        if (next->fork_bases[2]) {
          walk_result = 'F';
          break;
        }
        if (next->fork_nvotes[0] > next->fork_nvotes[1]) {
          ext = next->fork_bases[0];
          nvotes = next->fork_nvotes[0];
        } else {
          ext = next->fork_bases[1];
          nvotes = next->fork_nvotes[1];
        }
        have_forked = true;;
      } else {
        if (next->fork_bases[0]) walk_result = 'F';
        else walk_result = 'X';
        break;
      }
    }
    mer.erase(0, 1);
    mer += ext;
    assert(mer.length() == mer_len);;
    walk += ext;
    tot_ext_votes += nvotes;
    nsteps++;
    if (nsteps >= walk_len_limit) break;
  }
  if (nsteps > 0) depth = tot_ext_votes / nsteps;
  if (nsteps < 0) DIE("nsteps is negative ", nsteps, ", walk result ", walk_result, "\n");
  if (depth > 0 && walk.empty()) DIE("walk depth is ", depth, " but there is no walk, nsteps ", nsteps, ", result ", walk_result);
  return walk_result;
}

*/


struct MerFreqs {
  // how many times this kmer has occurred: don't need to count beyond 65536
  // count of high quality extensions
  ExtCounts exts;
  // the final extensions chosen - A,C,G,T, or F,X
  char ext;
  // the count of the final extension
  int count;

  void set_ext(double dynamic_min_depth, int count) {
    auto sorted_exts = exts.get_sorted();
    ext = sorted_exts[0].first;
    int ext_max = sorted_exts[0].second;
    int ext_min = sorted_exts[1].second;
    int dmin_dyn = (1.0 - dynamic_min_depth) * count;      // integer floor
    if (dmin_dyn < 2) dmin_dyn = 2;
    count = ext_max;
    if (ext_max < dmin_dyn) {
      ext = 'X';
      count = 0;
    } else if (ext_min >= dmin_dyn) {
      ext = 'F';
      count = 0;
    }
  }
  
};


using MerMap = unordered_map<string, MerFreqs>;


static void count_mers(vector<ReadSeq> &reads, MerMap &mers_ht, int seq_depth, int mer_len, int qual_offset,
                       double dynamic_min_depth) {
  // split reads into kmers and count frequency of high quality extensions
  for (auto &read_seq : reads) {
    if (mer_len >= read_seq.seq.length()) continue;
    int num_mers = read_seq.seq.length() - mer_len;
    for (int start = 0; start < num_mers; start++) {
      // skip mers that contain Ns
      if (read_seq.seq.find("N", start) != string::npos) continue;
      int ext_pos = start + mer_len;
      if (ext_pos >= read_seq.seq.length()) DIE("extension_pos ", ext_pos, " out of range ", mer_len);
      char ext = read_seq.seq[ext_pos];
      if (read_seq.quals[ext_pos] < qual_offset + QUAL_CUTOFF || ext == 'N') ext = 0;
      string mer = read_seq.seq.substr(start, mer_len);
      auto it = mers_ht.find(mer);
      if (it == mers_ht.end()) {
        mers_ht.insert({mer, {.exts = {0}, .ext = 0, .count = 0}});
        it = mers_ht.find(mer);
      }
      it->second.exts.inc(ext, 1);
    }
  }
  // now set extension choices
  for (auto &elem : mers_ht) {
    elem.second.set_ext(dynamic_min_depth, seq_depth);
  }    
}

// return the result of the walk (f, r or x)
static char walk_mers(MerMap &mers_ht, string &mer, string &walk, int mer_len, int walk_len_limit) {
  bool have_forked = false;
  int nsteps = 0;
  int tot_ext_votes = 0;
  unordered_map<string, bool> loop_check_ht;
  char walk_result = 'X';
  for (int nsteps = 0; nsteps < walk_len_limit; nsteps++) {
    // check for a cycle in the graph
    if (loop_check_ht.find(mer) != loop_check_ht.end()) {
      walk_result = 'R';
      break;
    } else {
      loop_check_ht.insert({mer, true});
    }
    auto it = mers_ht.find(mer);
    if (it == mers_ht.end()) {
      walk_result = 'X';
      break;
    }
    char ext = it->second.ext;
    if (ext == 'F' || ext == 'X') {
      walk_result = ext;
      break;
    }
    mer.erase(0, 1);
    mer += ext;
    walk += ext;
  }
  return walk_result;
}


static string iterative_walks(string &seq, int seq_depth, vector<ReadSeq> &reads, int max_mer_len, int kmer_len,
                              int qual_offset, double dynamic_min_depth, int walk_len_limit, array<int64_t, 3> &term_counts,
                              int64_t &num_walks, int64_t &max_walk_len, int64_t &sum_ext) {
  int min_mer_len = LASSM_MIN_KMER_LEN;
  max_mer_len = min(max_mer_len, (int)seq.length());
  // iteratively walk starting from kmer_size, increasing mer size on a fork (F) or repeat (R),
  // and decreasing on an end of path (X)
  // look for the longest walk. Note we have to restart from the beginning for each walk to ensure
  // that all loops will be detected
  string longest_walk = "";
  int shift = 0;
  DBG("  reads:\n");
#ifdef DEBUG
  for (auto &read_seq : reads) {
    DBG("    ", read_seq.read_id, "\n", read_seq.seq, "\n");
  }
#endif
  for (int mer_len = kmer_len; mer_len >= min_mer_len && mer_len <= max_mer_len; mer_len += shift) {
    //unordered_map<string, MerBase> mers_ht;
    MerMap mers_ht;
    //compute_mer_freqs(reads, mers_ht, seq_depth, mer_len, qual_offset);
    count_mers(reads, mers_ht, seq_depth, mer_len, qual_offset, dynamic_min_depth);
    string mer = seq.substr(seq.length() - mer_len);
    string walk = "";
    char walk_result = walk_mers(mers_ht, mer, walk, mer_len, walk_len_limit);
    int walk_len = walk.length();
    if (walk_len > longest_walk.length()) longest_walk = walk;
    if (walk_result == 'X') {
      term_counts[0]++;
      // walk reaches a dead-end, downshift, unless we were upshifting
      if (shift == LASSM_SHIFT_SIZE) break;
      shift = -LASSM_SHIFT_SIZE;
    } else {
      if (walk_result == 'F') term_counts[1]++;
      else term_counts[2]++;
      // otherwise walk must end with a fork or repeat, so upshift
      if (shift == -LASSM_SHIFT_SIZE) break;
      if (mer_len > seq.length()) break;
      shift = LASSM_SHIFT_SIZE;
    }
  }
  if (!longest_walk.empty()) {
    num_walks++;
    max_walk_len = max(max_walk_len, (int64_t)longest_walk.length());
    sum_ext += longest_walk.length();
  }
  return longest_walk;
}


void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len,
                 int qual_offset, double dynamic_min_depth) {
  Timer timer(__func__, true);
  //_num_forks = 0;
  // _num_terms = 0;
  //_num_repeats = 0;
  // walk should never be more than this. Note we use the maximum insert size from all libraries
  int walk_len_limit = insert_avg + 2 * insert_stddev;
  int64_t num_walks = 0, sum_clen = 0, sum_ext = 0, max_walk_len = 0, num_reads = 0, num_sides = 0;
  array<int64_t, 3> term_counts = {0};
  ProgressBar progbar(ctgs_dht.get_local_num_ctgs(), "Extending contigs");
  for (auto ctg = ctgs_dht.get_first_local_ctg(); ctg != nullptr; ctg = ctgs_dht.get_next_local_ctg()) {
    progbar.update();
    sum_clen += ctg->seq.length();
    if (ctg->reads_right.size()) {
      num_sides++;
      num_reads += ctg->reads_right.size();
      DBG("walk right ctg ", ctg->cid, " ", ctg->depth, "\n", ctg->seq, "\n");
      // have to do right first because the contig needs to be revcomped for the left
      string right_walk = iterative_walks(ctg->seq, ctg->depth, ctg->reads_right, max_kmer_len, kmer_len, qual_offset,
                                          dynamic_min_depth, walk_len_limit, term_counts, num_walks, max_walk_len, sum_ext);
      if (!right_walk.empty()) ctg->seq += right_walk;
    }
    if (ctg->reads_left.size()) {
      num_sides++;
      num_reads += ctg->reads_left.size();
      string seq_rc = revcomp(ctg->seq);
      DBG("walk left ctg ", ctg->cid, " ", ctg->depth, "\n", seq_rc, "\n");
      string left_walk = iterative_walks(seq_rc, ctg->depth, ctg->reads_left, max_kmer_len, kmer_len, qual_offset,
                                         dynamic_min_depth, walk_len_limit, term_counts, num_walks, max_walk_len, sum_ext);
      if (!left_walk.empty()) ctg->seq.insert(0, left_walk);
    }
  }
  progbar.done();
  barrier();
  SLOG_VERBOSE("Walk terminations: ",
               reduce_one(term_counts[0], op_fast_add, 0).wait(), " X, ",
               reduce_one(term_counts[1], op_fast_add, 0).wait(), " F, ",
               reduce_one(term_counts[2], op_fast_add, 0).wait(), " R\n");
              
  auto tot_num_walks = reduce_one(num_walks, op_fast_add, 0).wait();
  auto tot_sum_ext = reduce_one(sum_ext, op_fast_add, 0).wait();
  auto tot_sum_clen = reduce_one(sum_clen, op_fast_add, 0).wait();
  auto tot_max_walk_len = reduce_one(max_walk_len, op_fast_max, 0).wait();
  SLOG_VERBOSE("Could walk ", perc_str(reduce_one(num_sides, op_fast_add, 0).wait(), ctgs_dht.get_num_ctgs() * 2),
               " contig sides\n");
  SLOG_VERBOSE("Found ", tot_num_walks, " walks, total extension length ", tot_sum_ext, " extended ", 
               (double)(tot_sum_ext + tot_sum_clen) / tot_sum_clen, "\n");
  if (tot_num_walks) 
    SLOG_VERBOSE("Average walk length ", tot_sum_ext / tot_num_walks, ", max walk length ", tot_max_walk_len, "\n");
}

