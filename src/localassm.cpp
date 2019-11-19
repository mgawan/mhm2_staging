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


void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, int insert_avg, int insert_stddev, int max_mer_len, int kmer_len, int qual_offset,
                 double dynamic_min_depth);


enum class AlnStatus { NO_ALN, OVERLAPS_CONTIG, EXTENDS_CONTIG };

struct CtgInfo {
  int64_t cid;
  char orient;
  char side;
};


class ReadsToCtgsDHT {
  
  using reads_to_ctgs_map_t = unordered_map<string, vector<CtgInfo> >;
  dist_object<reads_to_ctgs_map_t> reads_to_ctgs_map;

  size_t get_target_rank(string read_id) {
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


static void process_reads(int kmer_len, vector<string> &reads_fname_list, ReadsToCtgsDHT &reads_to_ctgs, CtgsWithReadsDHT &ctgs_dht) {
  Timer timer(__func__, true);
  int64_t num_reads = 0;
  int64_t num_read_maps_found = 0;
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    FastqReader fqr(merged_reads_fname, PER_RANK_FILE);
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

  Timer timer(__func__, true);
  int64_t num_alns_found = 0, num_alns_invalid = 0, num_direct = 0, num_proj = 0;
  int min_pair_len = insert_avg + 3 * insert_stddev;
  int64_t max_alns = 0;
  IntermittentTimer t_get_alns("get alns reads to contigs");
  int64_t aln_i = 0;
  AlnStatus start_status, end_status;
  ProgressBar progbar(alns.size(), "Getting read-to-contig mappings from alignments");
  while (aln_i < alns.size()) {
    Aln aln;
    t_get_alns.start();
    get_best_aln_for_read(alns, aln_i, aln, start_status, end_status, num_alns_found, num_alns_invalid);
    t_get_alns.stop();
    progbar.update(aln_i);
    int pair_num = aln.read_id.back() - '0';
    // add a direct extension to the contig, start or end
    if (start_status == AlnStatus::EXTENDS_CONTIG) {
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient, aln.orient == '+' ? 'L' : 'R');
      num_direct++;
    } else if (end_status == AlnStatus::EXTENDS_CONTIG) {
      reads_to_ctgs.add(aln.read_id, aln.cid, aln.orient, aln.orient == '+' ? 'R' : 'L');
      num_direct++;
    }
    // add mate pair if feasible
    if (!pair_overlap(aln, min_pair_len)) {
      // indicate the other pair_num
      string pair_read_id = aln.read_id;
      pair_read_id[aln.read_id.length() - 1] = (pair_num == 1 ? '2' : '1');
      reads_to_ctgs.add(pair_read_id, aln.cid, aln.orient == '+' ? '-' : '+', aln.orient == '+' ? 'R' : 'L');
      num_proj++;
    }
    progress();
  }
  progbar.done();
  barrier();
  auto tot_alns_found = reduce_one(num_alns_found, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed ", tot_alns_found, " alignments:\n");
  SLOG_VERBOSE("  invalid:   ", perc_str(reduce_one(num_alns_invalid, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("  direct:    ", perc_str(reduce_one(num_direct, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("  projected: ", perc_str(reduce_one(num_proj, op_fast_add, 0).wait(), tot_alns_found), "\n");
  SLOG_VERBOSE("Added ", reads_to_ctgs.get_num_mappings(), " mappings\n");
}


void add_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs) {
  Timer timer(__func__, true);
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


void localassm(int max_kmer_len, int kmer_len, vector<string> &reads_fname_list, int insert_avg, int insert_stddev,
               int qual_offset, double dynamic_min_depth, Contigs &ctgs, Alns &alns) {
  Timer timer(__func__, true);
  CtgsWithReadsDHT ctgs_dht(ctgs.size());
  add_ctgs(ctgs_dht, ctgs);
  ReadsToCtgsDHT reads_to_ctgs(100);
  // extract read id to ctg id mappings from alignments
  process_alns(alns, reads_to_ctgs, insert_avg, insert_stddev);
  // extract read seqs and add to ctgs
  process_reads(max_kmer_len, reads_fname_list, reads_to_ctgs, ctgs_dht);
  // extend contigs using locally mapped reads
  extend_ctgs(ctgs_dht, insert_avg, insert_stddev, max_kmer_len, kmer_len, qual_offset, dynamic_min_depth);
  // add extended ctg sequences back to the original ctgs
  
}
