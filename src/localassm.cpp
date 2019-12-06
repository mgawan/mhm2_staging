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
#include "kmer_dht.hpp"


using namespace std;
using namespace upcxx;


enum class AlnStatus { NO_ALN, OVERLAPS_CONTIG, EXTENDS_CONTIG };

class ReadsToCtgsDHT {
  struct CtgInfo {
    int64_t cid;
    char orient;
    char side;
  };

  using reads_to_ctgs_map_t = HASH_TABLE<string, vector<CtgInfo> >;
  dist_object<reads_to_ctgs_map_t> reads_to_ctgs_map;

  size_t get_target_rank(const string &read_id) {
    return std::hash<string>{}(read_id) % rank_n();
  }
  
public:
  ReadsToCtgsDHT(int64_t initial_size) 
    : reads_to_ctgs_map({}) {
    reads_to_ctgs_map->reserve(initial_size);
  }
  
  void add(const string &read_id, int64_t cid, char orient, char side) {
    CtgInfo ctg_info = { .cid = cid, .orient = orient, .side = side };
    rpc(get_target_rank(read_id),
        [](dist_object<reads_to_ctgs_map_t> &reads_to_ctgs_map, string read_id, CtgInfo ctg_info) {
          const auto it = reads_to_ctgs_map->find(read_id);
          if (it == reads_to_ctgs_map->end()) reads_to_ctgs_map->insert({read_id, {ctg_info}});
          else it->second.push_back(ctg_info);
        }, reads_to_ctgs_map, read_id, ctg_info).wait();
  }

  int64_t get_num_mappings() {
    return reduce_one(reads_to_ctgs_map->size(), op_fast_add, 0).wait();
  }

  vector<CtgInfo> get_ctgs(string &read_id) {
    return upcxx::rpc(get_target_rank(read_id),
                      [](upcxx::dist_object<reads_to_ctgs_map_t> &reads_to_ctgs_map, string read_id) -> vector<CtgInfo> {
                        const auto it = reads_to_ctgs_map->find(read_id);
                        if (it == reads_to_ctgs_map->end()) return {};
                        return it->second;
                      }, reads_to_ctgs_map, read_id).wait();
  }

};


struct ReadSeq {
  string read_id;
  string seq;
  string quals;
};


struct CtgWithReads {
  int64_t cid;
  string seq;
  double depth;
  vector<ReadSeq> reads_left;
  vector<ReadSeq> reads_right;
};


class CtgsWithReadsDHT {
  
  using ctgs_map_t = HASH_TABLE<int64_t, CtgWithReads>;
  dist_object<ctgs_map_t> ctgs_map;
  ctgs_map_t::iterator ctgs_map_iter;

  size_t get_target_rank(int64_t cid) {
    return std::hash<int64_t>{}(cid) % rank_n();
  }
  
public:

  CtgsWithReadsDHT(int64_t num_ctgs)
    : ctgs_map({}) {
    // pad the local ctg count a bit for this estimate
    ctgs_map->reserve(num_ctgs * 1.2);
  }

  void add_ctg(Contig &ctg) {
    rpc(get_target_rank(ctg.id),
        [](dist_object<ctgs_map_t> &ctgs_map, int64_t cid, string seq, double depth) {
          const auto it = ctgs_map->find(cid);
          if (it != ctgs_map->end()) DIE("Found duplicate ctg ", cid);
          CtgWithReads ctg_with_reads = { .cid = cid, .seq = seq, .depth = depth, .reads_left = {}, .reads_right = {} };
          ctgs_map->insert({cid, ctg_with_reads });
        }, ctgs_map, ctg.id, ctg.seq, ctg.depth).wait();
  }

  void add_read(int64_t cid, char side, ReadSeq read_seq) {
    rpc(get_target_rank(cid),
        [](dist_object<ctgs_map_t> &ctgs_map, int64_t cid, char side, string read_id, string seq, string quals) {
          const auto it = ctgs_map->find(cid);
          if (it == ctgs_map->end()) DIE("Could not find ctg ", cid);
          if (side == 'L') it->second.reads_left.push_back({read_id, seq, quals});
          else it->second.reads_right.push_back({read_id, seq, quals});
        }, ctgs_map, cid, side, read_seq.read_id, read_seq.seq, read_seq.quals);
  }
  
  int64_t get_num_ctgs() {
    return reduce_one(ctgs_map->size(), op_fast_add, 0).wait();
  }

  int64_t get_local_num_ctgs() {
    return ctgs_map->size();
  }

  CtgWithReads *get_first_local_ctg() {
    ctgs_map_iter = ctgs_map->begin();
    if (ctgs_map_iter == ctgs_map->end()) return nullptr;
    auto ctg = &ctgs_map_iter->second;
    ctgs_map_iter++;
    return ctg;
  }
  
  CtgWithReads *get_next_local_ctg() {
    if (ctgs_map_iter == ctgs_map->end()) return nullptr;
    auto ctg = &ctgs_map_iter->second;
    ctgs_map_iter++;
    return ctg;
  }
};


struct MerFreqs {
  // how many times this kmer has occurred: don't need to count beyond 65536
  // count of high quality extensions and low quality extensions - structure comes from kmer_dht.hpp
  ExtCounts hi_q_exts, low_q_exts;
  // the final extensions chosen - A,C,G,T, or F,X
  char ext;
  // the count of the final extension
  int count;

  struct MerBase {
    char base;
    uint16_t nvotes_hi_q, nvotes, rating;

    uint16_t get_base_rating(int depth) {
      double min_viable = max(LASSM_MIN_VIABLE_DEPTH * depth, 2.0);
      double min_expected_depth = max(LASSM_MIN_EXPECTED_DEPTH * depth, 2.0);
      if (nvotes == 0) return 0;
      if (nvotes == 1) return 1;
      if (nvotes < min_viable) return 2;
      if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q < min_viable) return 3;
      if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q >= min_viable) return 4;
      if (nvotes >= min_expected_depth && nvotes_hi_q < min_viable) return 5;
      if (nvotes >= min_expected_depth && min_viable < nvotes_hi_q && nvotes_hi_q < min_expected_depth) return 6;
      return 7;
    }
  };
  
  void set_ext(double dynamic_min_depth, int seq_depth) {
    // set extension similarly to how it is done with localassm in mhm
    MerBase mer_bases[4] = {{.base = 'A', .nvotes_hi_q = hi_q_exts.count_A, .nvotes = low_q_exts.count_A},
                            {.base = 'C', .nvotes_hi_q = hi_q_exts.count_C, .nvotes = low_q_exts.count_C},
                            {.base = 'G', .nvotes_hi_q = hi_q_exts.count_G, .nvotes = low_q_exts.count_G},
                            {.base = 'T', .nvotes_hi_q = hi_q_exts.count_T, .nvotes = low_q_exts.count_T}};
    for (int i = 0; i < 4; i++) {
      mer_bases[i].rating = mer_bases[i].get_base_rating(seq_depth);
    }
    // sort bases in descending order of quality
    sort(mer_bases, mer_bases + sizeof(mer_bases) / sizeof(mer_bases[0]),
         [](const auto &elem1, const auto &elem2) -> bool {
           if (elem1.rating != elem2.rating) return elem1.rating > elem2.rating;
           if (elem1.nvotes_hi_q != elem2.nvotes_hi_q) return elem1.nvotes_hi_q > elem2.nvotes_hi_q;
           if (elem1.nvotes != elem2.nvotes) return elem1.nvotes > elem2.nvotes;
           return true;
         });
    int top_rating = mer_bases[0].rating;
    int runner_up_rating = mer_bases[1].rating;
    if (top_rating < runner_up_rating) DIE("top_rating ", top_rating, " < ", runner_up_rating, "\n");
    assert(top_rating >= runner_up_rating);
    int top_rated_base = mer_bases[0].base;
    ext = 'X';
    count = 0;
    // no extension (base = 0) if the runner up is close to the top rating
    // except, if rating is 7 (best quality), then all bases of rating 7 are forks
    if (top_rating > LASSM_RATING_THRES) {         // must have at least minViable bases
      if (top_rating <= 3) {    // must be uncontested
        if (runner_up_rating == 0) ext = top_rated_base;
      } else if (top_rating < 6) {
        if (runner_up_rating < 3) ext = top_rated_base;
      } else if (top_rating == 6) {  // viable and fair hiQ support
        if (runner_up_rating < 4) ext = top_rated_base;
      } else {                     // strongest rating trumps
        if (runner_up_rating < 7) {
          ext = top_rated_base;
        } else {
          if (mer_bases[2].rating == 7 || mer_bases[0].nvotes == mer_bases[1].nvotes) ext = 'F'; 
          else if (mer_bases[0].nvotes > mer_bases[1].nvotes) ext = mer_bases[0].base;
          else if (mer_bases[1].nvotes > mer_bases[0].nvotes) ext = mer_bases[1].base;
        }
      }
    }
    for (int i = 0; i < 4; i++) {
      if (mer_bases[i].base == ext) {
        count = mer_bases[i].nvotes;
        break;
      }
    }
  }
  
};


using MerMap = HASH_TABLE<string, MerFreqs>;

static void process_reads(int kmer_len, vector<string> &reads_fname_list, ReadsToCtgsDHT &reads_to_ctgs, CtgsWithReadsDHT &ctgs_dht) {
  Timer timer(__FILEFUNC__);
  int64_t num_reads = 0;
  int64_t num_read_maps_found = 0;
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    FastqReader fqr(merged_reads_fname);
    string id, seq, quals;
    ProgressBar progbar(fqr.my_file_size(), "Processing reads");
    size_t tot_bytes_read = 0;
    while (true) {
      progress();
      size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      // this happens when we have a placeholder entry because reads merged
      if (kmer_len > seq.length()) continue;
      num_reads++;
      string seq_rc = revcomp(seq);
      string quals_rc = quals;
      reverse(quals_rc.begin(), quals_rc.end());
      auto ctgs = reads_to_ctgs.get_ctgs(id);
      if (ctgs.size()) {
        num_read_maps_found++;
        for (auto &ctg : ctgs) {
          if ((ctg.orient == '-' && ctg.side == 'R') || (ctg.orient == '+' && ctg.side == 'L'))
            ctgs_dht.add_read(ctg.cid, ctg.side, {id, seq_rc, quals_rc});
          else
            ctgs_dht.add_read(ctg.cid, ctg.side, {id, seq, quals});
        }
      }
    }
    progbar.done();
  }
  barrier();
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(reduce_one(num_read_maps_found, op_fast_add, 0).wait(), tot_num_reads),
               " reads that map to contigs\n");
}


static void get_best_aln_for_read(Alns &alns, int64_t &i, Aln &best_aln, AlnStatus &best_start_status, AlnStatus &best_end_status,
                                  int64_t &num_alns_found, int64_t &num_alns_invalid) {

  auto classify_aln = [](int runaligned, int cunaligned) ->AlnStatus {
    if (runaligned > cunaligned && cunaligned < ALN_WIGGLE) return AlnStatus::EXTENDS_CONTIG;
    if (runaligned <= cunaligned && runaligned < ALN_WIGGLE) return AlnStatus::OVERLAPS_CONTIG;
    return AlnStatus::NO_ALN;
  };
  
  // choose the highest scoring aln for this read that is useful
  best_start_status = AlnStatus::NO_ALN;
  best_end_status = AlnStatus::NO_ALN;
  string start_read_id = "";
  int best_aln_score = 0;
  best_aln.read_id = "";
  for (; i < alns.size(); i++) {
    Aln aln = alns.get_aln(i);
    // alns for a new read
    if (start_read_id != "" && aln.read_id != start_read_id) return;
    num_alns_found++;
    if (aln.score1 < best_aln_score) continue;
    AlnStatus start_status, end_status;
    if (aln.orient == '+') {
      start_status = classify_aln(aln.rstart - 1, aln.cstart - 1);
      end_status = classify_aln(aln.rlen - aln.rstop, aln.clen - aln.cstop);
    } else {
      // for '-' strand, aln is between read and revcomp of contig
      start_status = classify_aln(aln.rstart - 1, aln.clen - aln.cstop);
      end_status = classify_aln(aln.rlen - aln.rstop, aln.cstart - 1);
    }
    if (start_status == AlnStatus::NO_ALN || end_status == AlnStatus::NO_ALN) {
      num_alns_invalid++;
      continue;
    }
    best_aln = aln;
    best_aln_score = aln.score1;
    best_start_status = start_status;
    best_end_status = end_status;
    start_read_id = aln.read_id;
  }
}


void process_alns(Alns &alns, ReadsToCtgsDHT &reads_to_ctgs, int insert_avg, int insert_stddev) {
  auto pair_overlap = [](Aln &aln, int min_pair_len) -> bool {
    // make sure that the mate won't overlap the same contig
    if (aln.orient == '+') {
      if (min_pair_len - aln.rlen - aln.rstart + 1 <= aln.clen - aln.cstart) return true;
    } else {
      if (min_pair_len - 2 * aln.rlen + aln.rstart - 1 <= aln.cstart) return true;
    }
    return false;
  };

  Timer timer(__FILEFUNC__);
  int64_t num_alns_found = 0, num_alns_invalid = 0, num_direct = 0, num_proj = 0;
  int min_pair_len = insert_avg + 3 * insert_stddev;
  int64_t max_alns = 0;
  IntermittentTimer t_get_alns(__FILENAME__ + string(":") + "get alns reads to contigs");
  int64_t aln_i = 0;
  AlnStatus start_status, end_status;
  ProgressBar progbar(alns.size(), "Getting read-to-contig mappings from alignments");
  while (aln_i < alns.size()) {
    progress();
    Aln aln;
    t_get_alns.start();
    get_best_aln_for_read(alns, aln_i, aln, start_status, end_status, num_alns_found, num_alns_invalid);
    t_get_alns.stop();
    progbar.update(aln_i);
    if (aln.read_id.empty()) continue;
    // add a direct extension to the contig, start or end
    if (start_status == AlnStatus::EXTENDS_CONTIG) {
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient, aln.orient == '+' ? 'L' : 'R');
      num_direct++;
    } else if (end_status == AlnStatus::EXTENDS_CONTIG) {
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient, aln.orient == '+' ? 'R' : 'L');
      num_direct++;
    }
    // FIXME: if read is longer than one pair read length, don't look for mate since it is a merged read
    // add mate pair if feasible
    if (!pair_overlap(aln, min_pair_len)) {
      // indicate the other pair number
      int len = aln.read_id.length();
      assert(len > 1);
      if (aln.read_id[len - 1] == '1') aln.read_id[len - 1] = '2';
      else if (aln.read_id[len - 1] == '2') aln.read_id[len - 1] = '1';
      else DIE("Bad pair number ", (int)aln.read_id[len - 1], " in read: ", aln.read_id);
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient == '+' ? '-' : '+', aln.orient == '+' ? 'R' : 'L');
      num_proj++;
    }
  }
  progbar.done();
  barrier();
  t_get_alns.done_barrier();
  auto tot_alns_found = reduce_one(num_alns_found, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed ", tot_alns_found, " alignments:\n");
  SLOG_VERBOSE("  invalid:   ", perc_str(reduce_one(num_alns_invalid, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("  direct:    ", perc_str(reduce_one(num_direct, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("  projected: ", perc_str(reduce_one(num_proj, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("Added ", reads_to_ctgs.get_num_mappings(), " mappings\n");
}


static void add_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs) {
  Timer timer(__FILEFUNC__);
  // process the local ctgs and insert into the distributed hash table
  ProgressBar progbar(ctgs.size(), "Adding contigs to distributed hash table");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    progbar.update();
    ctgs_dht.add_ctg(*it);
    progress();
  }
  progbar.done();
  barrier();
  SLOG_VERBOSE("Added ", ctgs_dht.get_num_ctgs(), " contigs\n");
}


static void count_mers(vector<ReadSeq> &reads, MerMap &mers_ht, int seq_depth, int mer_len, int qual_offset,
                       double dynamic_min_depth, int64_t &excess_reads) {
  int num_reads = 0;
  // split reads into kmers and count frequency of high quality extensions
  for (auto &read_seq : reads) {
    num_reads++;
    if (num_reads > LASSM_MAX_COUNT_MERS_READS) {
      excess_reads += reads.size() - LASSM_MAX_COUNT_MERS_READS;
      break;
    }
    progress();
    if (mer_len >= read_seq.seq.length()) continue;
    int num_mers = read_seq.seq.length() - mer_len;
    for (int start = 0; start < num_mers; start++) {
      // skip mers that contain Ns
      if (read_seq.seq.find("N", start) != string::npos) continue;
      string mer = read_seq.seq.substr(start, mer_len);
      auto it = mers_ht.find(mer);
      if (it == mers_ht.end()) {
        mers_ht.insert({mer, {.hi_q_exts = {0}, .low_q_exts = {0}, .ext = 0, .count = 0}});
        it = mers_ht.find(mer);
      }
      int ext_pos = start + mer_len;
      assert(ext_pos < read_seq.seq.length());
      char ext = read_seq.seq[ext_pos];
      if (ext == 'N') continue;
      int qual = read_seq.quals[ext_pos] - qual_offset;
      if (qual >= LASSM_MIN_QUAL) it->second.low_q_exts.inc(ext, 1);
      if (qual >= LASSM_MIN_HI_QUAL) it->second.hi_q_exts.inc(ext, 1);
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
  HASH_TABLE<string, bool> loop_check_ht;
  char walk_result = 'X';
  for (int nsteps = 0; nsteps < walk_len_limit; nsteps++) {
    if (!(nsteps % 10)) progress();
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
                              int64_t &num_walks, int64_t &max_walk_len, int64_t &sum_ext, IntermittentTimer &count_mers_timer,
                              IntermittentTimer &walk_mers_timer, int64_t &excess_reads) {
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
    count_mers_timer.start();
    MerMap mers_ht;
    count_mers(reads, mers_ht, seq_depth, mer_len, qual_offset, dynamic_min_depth, excess_reads);
    count_mers_timer.stop();
    string mer = seq.substr(seq.length() - mer_len);
    string walk = "";
    walk_mers_timer.start();
    char walk_result = walk_mers(mers_ht, mer, walk, mer_len, walk_len_limit);
    walk_mers_timer.stop();
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


static void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs, int insert_avg, int insert_stddev, int max_kmer_len,
                        int kmer_len, int qual_offset, double dynamic_min_depth) {
  Timer timer(__FILEFUNC__);
  // walk should never be more than this. Note we use the maximum insert size from all libraries
  int walk_len_limit = insert_avg + 2 * insert_stddev;
  int64_t num_walks = 0, sum_clen = 0, sum_ext = 0, max_walk_len = 0, num_reads = 0, num_sides = 0, max_num_reads = 0,
    excess_reads = 0;
  array<int64_t, 3> term_counts = {0};
  IntermittentTimer count_mers_timer(__FILENAME__ + string(":") + "count_mers"),
    walk_mers_timer(__FILENAME__ + string(":") + "walk_mers");
  ProgressBar progbar(ctgs_dht.get_local_num_ctgs(), "Extending contigs");
  for (auto ctg = ctgs_dht.get_first_local_ctg(); ctg != nullptr; ctg = ctgs_dht.get_next_local_ctg()) {
    progress();
    progbar.update();
    sum_clen += ctg->seq.length();
    if (ctg->reads_right.size()) {
      num_sides++;
      num_reads += ctg->reads_right.size();
      max_num_reads = max(max_num_reads, (int64_t)ctg->reads_right.size());
      DBG("walk right ctg ", ctg->cid, " ", ctg->depth, "\n", ctg->seq, "\n");
      // have to do right first because the contig needs to be revcomped for the left
      string right_walk = iterative_walks(ctg->seq, ctg->depth, ctg->reads_right, max_kmer_len, kmer_len, qual_offset,
                                          dynamic_min_depth, walk_len_limit, term_counts, num_walks, max_walk_len, sum_ext,
                                          count_mers_timer, walk_mers_timer, excess_reads);
      if (!right_walk.empty()) ctg->seq += right_walk;
    }
    if (ctg->reads_left.size()) {
      num_sides++;
      num_reads += ctg->reads_left.size();
      max_num_reads = max(max_num_reads, (int64_t)ctg->reads_left.size());
      string seq_rc = revcomp(ctg->seq);
      DBG("walk left ctg ", ctg->cid, " ", ctg->depth, "\n", seq_rc, "\n");
      string left_walk = iterative_walks(seq_rc, ctg->depth, ctg->reads_left, max_kmer_len, kmer_len, qual_offset,
                                         dynamic_min_depth, walk_len_limit, term_counts, num_walks, max_walk_len, sum_ext,
                                         count_mers_timer, walk_mers_timer, excess_reads);
      if (!left_walk.empty()) {
        left_walk = revcomp(left_walk);
        ctg->seq.insert(0, left_walk);
      }
    }
    ctgs.add_contig({.id = ctg->cid, .seq = ctg->seq, .depth = ctg->depth});
  }
  progbar.done();
  count_mers_timer.done_barrier();
  walk_mers_timer.done_barrier();
  barrier();
  SLOG_VERBOSE("Walk terminations: ",
               reduce_one(term_counts[0], op_fast_add, 0).wait(), " X, ",
               reduce_one(term_counts[1], op_fast_add, 0).wait(), " F, ",
               reduce_one(term_counts[2], op_fast_add, 0).wait(), " R\n");
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto tot_max_num_reads = reduce_one(max_num_reads, op_fast_max, 0).wait();
  auto tot_num_walks = reduce_one(num_walks, op_fast_add, 0).wait();
  auto tot_sum_ext = reduce_one(sum_ext, op_fast_add, 0).wait();
  auto tot_sum_clen = reduce_one(sum_clen, op_fast_add, 0).wait();
  auto tot_max_walk_len = reduce_one(max_walk_len, op_fast_max, 0).wait();
  auto tot_excess_reads = reduce_one(excess_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Used a total of ", tot_num_reads, " reads, max per ctg ", max_num_reads, " avg per ctg ",
               (tot_num_reads / ctgs_dht.get_num_ctgs()), ", dropped ", perc_str(tot_excess_reads, tot_num_reads),
               " excess reads\n");
  SLOG_VERBOSE("Could walk ", perc_str(reduce_one(num_sides, op_fast_add, 0).wait(), ctgs_dht.get_num_ctgs() * 2),
               " contig sides\n");
  if (tot_sum_clen) 
    SLOG_VERBOSE("Found ", tot_num_walks, " walks, total extension length ", tot_sum_ext, " extended ", 
                 (double)(tot_sum_ext + tot_sum_clen) / tot_sum_clen, "\n");
  if (tot_num_walks) 
    SLOG_VERBOSE("Average walk length ", tot_sum_ext / tot_num_walks, ", max walk length ", tot_max_walk_len, "\n");
}


void localassm(int max_kmer_len, int kmer_len, vector<string> &reads_fname_list, int insert_avg, int insert_stddev,
               int qual_offset, double dynamic_min_depth, Contigs &ctgs, Alns &alns) {
  Timer timer(__FILEFUNC__);
  CtgsWithReadsDHT ctgs_dht(ctgs.size());
  add_ctgs(ctgs_dht, ctgs);
  ReadsToCtgsDHT reads_to_ctgs(100);
  // extract read id to ctg id mappings from alignments
  process_alns(alns, reads_to_ctgs, insert_avg, insert_stddev);
  // extract read seqs and add to ctgs
  process_reads(max_kmer_len, reads_fname_list, reads_to_ctgs, ctgs_dht);
  // clear out the local contigs
  ctgs.clear();
  ctgs.set_capacity(ctgs_dht.get_local_num_ctgs());
  // extend contigs using locally mapped reads
  extend_ctgs(ctgs_dht, ctgs, insert_avg, insert_stddev, max_kmer_len, kmer_len, qual_offset, dynamic_min_depth);
}

