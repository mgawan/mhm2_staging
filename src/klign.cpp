/*
  Sequence alignment similar to meraligner.
  c shofmeyr@lbl.gov
  Nov 2018
*/


#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "progressbar.hpp"
#include "ssw.hpp"
#include "contigs.hpp"

using namespace std;
using namespace upcxx;


using cid_t = int64_t;

struct CtgLoc {
  cid_t cid;
  global_ptr<char> seq_gptr;
  int clen;
  int pos_in_ctg;
  bool is_rc;
};

int64_t _num_dropped = 0;

class KmerHashTable {
private:
  // no need for a global ptr - just need to know where to find the sequence with an rpc to the destination
  // so the string ptr is only valid on the destination rank
  using kmer_map_t = unordered_map<string, vector<CtgLoc> >;
  dist_object<kmer_map_t> kmer_map;
  int64_t num_perfect_alns;
  ofstream alns_file;
  IntermittentTimer t_ssw;

  // the key is to_string(cid) + rname
  unordered_map<string, bool> cid_reads_map;
  
  // default aligner and filter
  StripedSmithWaterman::Aligner ssw_aligner;
  StripedSmithWaterman::Filter ssw_filter;
  
  size_t get_target_rank(const string &kmer) {
    return std::hash<string>{}(kmer) % rank_n();
  }

public:

  int kmer_len;
  int64_t num_alns;
  
  
  // aligner construction: match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty - no entry for ambiquity, which should be 2
  // note SSW internal defaults are 2 2 3 1
  KmerHashTable(int kmer_len) : kmer_map({}), num_perfect_alns(0), ssw_aligner(1, 3, 5, 2, 2), num_alns(0),
                                t_ssw("SSW"), kmer_len(kmer_len) {
    ssw_filter.report_cigar = false;
  }

  void clear() {
    for (auto it = kmer_map->begin(); it != kmer_map->end(); ) {
      it = kmer_map->erase(it);
    }
  }
  
  ~KmerHashTable() {
    clear();
  }

  int64_t get_num_kmers(bool all = false) {
    if (!all) return reduce_one(kmer_map->size(), op_fast_add, 0).wait();
    else return reduce_all(kmer_map->size(), op_fast_add).wait();
  }
  
  int64_t get_num_perfect_alns(bool all = false) {
    if (!all) return reduce_one(num_perfect_alns, op_fast_add, 0).wait();
    else return reduce_all(num_perfect_alns, op_fast_add).wait();
  }

  int64_t get_num_dropped(bool all = false) {
    if (!all) return reduce_one(_num_dropped, op_fast_add, 0).wait();
    else return reduce_all(_num_dropped, op_fast_add).wait();
  }

  void open_alns_file(const string fname) {
    alns_file.open(fname);
  }

  void close_alns_file() {
    alns_file.close();
  }

  /*
  double count_entries() {
    int64_t tot = 0;
    int max_num = 0;
    for (auto entry : *kmer_map) {
      tot += entry.second.size();
      max_num = max((int)entry.second.size(), max_num);
    }
    SOUT("max num ", max_num, "\n");
    return (double)tot / kmer_map->size();
  }
  */

  void add_kmer(const string &kmer, const CtgLoc &ctg_loc, promise<> &prom) {
    rpc(get_target_rank(kmer), operation_cx::as_promise(prom),
        [](dist_object<kmer_map_t> &kmer_map, string kmer, CtgLoc ctg_loc) {
          const auto it = kmer_map->find(kmer);
          if (it == kmer_map->end()) {
            kmer_map->insert({kmer, {ctg_loc}});
          } else {
            vector<CtgLoc> *ctg_locs = &it->second;
            // limit the number of matching contigs to any given kmer - this is an explosion in the graph anyway
            if (ctg_locs->size() >= MAX_KMER_MAPPINGS) {
              _num_dropped++;
              return;
            }
            // only add it if we don't already have a mapping from this kmer to this ctg
            for (auto cl : it->second) 
              if (cl.cid == ctg_loc.cid) return;
            it->second.push_back(ctg_loc);
          }
        }, kmer_map, kmer, ctg_loc);
  }

  future<vector<CtgLoc> > get_ctgs_with_kmer(const string &kmer) {
    return rpc(get_target_rank(kmer),
               [](dist_object<kmer_map_t> &kmer_map, string kmer) -> vector<CtgLoc> {
                 const auto it = kmer_map->find(kmer);
                 if (it == kmer_map->end()) return {};
                 return it->second;
               }, kmer_map, kmer);
  }

  void align_read(const string &rname, const string &rseq, char *seq_buf, int rstart,
                  int start_pos, int end_pos, const CtgLoc &ctg_loc, char orient, int subseq_len) {
    bool perfect_match = true;
    if (strncmp(seq_buf, (orient == '-' ? rseq.c_str() : rseq.c_str()) + rstart, subseq_len) == 0) num_perfect_alns++;
    else perfect_match = false;
    int rend = rstart + subseq_len;
    progress();
    int offset = rstart;
    if (orient == '-') {
      rstart = rseq.length() - rend;
      rend = rseq.length() - offset;
    }

#ifdef TRACE_SEQS
    if (!rank_me()) {
      cerr << "\n" << rname << "\t" << ctg_loc.cid << "\t" << (orient == '+' ? "Plus " : "Minus ") << ctg_loc.pos_in_ctg << endl;
      if (offset) cerr << string(offset, ' ');
      cerr << seq_buf << endl;
      cerr << rseq << endl;
      if (perfect_match) cerr << " perfect\n";
      else cerr << " IMPERFECT\n";
    }
#endif

    if (perfect_match) {
      alns_file << "MERALIGNER-0\t" << rname << "\t" << rstart + 1 << "\t" << rend << "\t" << rseq.length() << "\t"
                << "Contig" + to_string(ctg_loc.cid) << "\t" << start_pos + 1 << "\t" << end_pos << "\t" << ctg_loc.clen << "\t"
                << (orient == '+' ? "Plus" : "Minus") << "\t0\t0\t0\t0\t" << (rend - rstart) << "\t0\n";
    } else {
      // contig is the ref, read is the query
      StripedSmithWaterman::Alignment aln;
      const char *query_str = rseq.c_str();
      t_ssw.start();
      ssw_aligner.Align(query_str, seq_buf, strlen(seq_buf), ssw_filter, &aln, max((int)(rseq.length() / 2), 15));
      t_ssw.stop();

      rstart = aln.query_begin;
      rend = aln.query_end + 1;
      if (orient == '-') {
        rstart = rseq.length() - aln.query_end - 1;
        rend = rseq.length() - aln.query_begin;
      }
      alns_file << "MERALIGNER-1\t" << rname << "\t" << rstart + 1 << "\t" << rend << "\t" << rseq.length() << "\t"
                << "Contig" + to_string(ctg_loc.cid) << "\t"
                << start_pos + aln.ref_begin + 1 << "\t" << start_pos + aln.ref_end + 1 << "\t" << ctg_loc.clen << "\t"
                << (orient == '+' ? "Plus" : "Minus") << "\t0\t0\t0\t0\t" << aln.sw_score << "\t" << aln.sw_score_next_best << endl;
    }
  }
};


static void build_alignment_index(KmerHashTable &kmer_ht, unsigned kmer_len, Contigs &ctgs) 
{
  Timer timer(__func__);
  int64_t tot_num_kmers = 0;
  promise<> prom;
  ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    global_ptr<char> seq_gptr = allocate<char>(ctg->seq.length() + 1);
    strcpy(seq_gptr.local(), ctg->seq.c_str());
    CtgLoc ctg_loc = { .cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length() };
    auto num_kmers = ctg->seq.length() - kmer_len + 1;
    tot_num_kmers += num_kmers;
    for (int i = 0; i < num_kmers; i++) {
      string&& kmer = ctg->seq.substr(i, kmer_len);
      string kmer_rc = revcomp(kmer);
      ctg_loc.pos_in_ctg = i;
      if (kmer_rc <= kmer) {
        ctg_loc.is_rc = false;
        kmer_ht.add_kmer(kmer, ctg_loc, prom);
      } else {
        ctg_loc.is_rc = true;
        kmer_ht.add_kmer(kmer_rc, ctg_loc, prom);
      }
      progress();
    }
  }
  prom.finalize().wait();
  progbar.done();
  barrier();
  auto num_kmers_added = reduce_one(tot_num_kmers, op_fast_add, 0).wait();
  auto num_kmers_in_ht = kmer_ht.get_num_kmers();
  SOUT("Processed ", num_kmers_added, " kmers\n");
  auto num_dropped = kmer_ht.get_num_dropped();
  if (num_dropped) {
    SOUT("Dropped ", num_dropped, " kmer-to-contig mappings (", 
         setprecision(2), fixed, (100.0 * num_dropped / num_kmers_added), "%)\n");
  }
}

/*
static void compute_alns_for_kmer(KmerHashTable *kmer_ht, vector<CtgLoc> &aligned_ctgs, unordered_map<cid_t, int> &cids_mapped,
                                  int pos_in_read, bool is_rc, const string &rname, const string &rseq, promise<> &prom) {
  string rseq_rc = "";
  int rlen = rseq.length();
  for (auto ctg_loc : aligned_ctgs) {
    progress();
    const auto it = cids_mapped.find(ctg_loc.cid);
    if (it == cids_mapped.end()) {
      cids_mapped.insert({ctg_loc.cid, pos_in_read});
      kmer_ht->num_alns++;
      char orient = '+';
      if (ctg_loc.is_rc != is_rc) {
        orient = '-';
        pos_in_read = rlen - (kmer_ht->kmer_len + pos_in_read);
        revcomp(rseq, &rseq_rc);
      }
      int pos_in_ctg_primer = max(ctg_loc.pos_in_ctg - pos_in_read - rlen, 0);
      int start_pos = pos_in_ctg_primer;
      int end_pos = min(ctg_loc.pos_in_ctg + rlen * 2, ctg_loc.clen);
      int subseq_len = end_pos - start_pos;
      if (subseq_len < 0) DIE("end_pos <= start_pos, ", end_pos, " ", start_pos);
      char *seq_buf = new char[subseq_len + 1];

      rget(ctg_loc.seq_gptr + start_pos, seq_buf, subseq_len,
           operation_cx::as_promise(prom)|operation_cx::as_future()).then(
             [=]() {
               seq_buf[subseq_len] = '\0';
               int rstart = 0;
               if (pos_in_read && subseq_len < rlen) {
                 if (orient == '-') rstart = max(pos_in_read - ctg_loc.pos_in_ctg, 0);
                 else rstart = pos_in_read;
               }
               kmer_ht->align_read(rname, orient == '+' ? rseq : rseq_rc, seq_buf, rstart, start_pos, end_pos, ctg_loc, orient, subseq_len);
               delete[] seq_buf;
             });

    }
  }
}
*/

static void do_alignments(KmerHashTable &kmer_ht, unsigned kmer_len)
{
  /*
  Timer timer(__func__);
  int64_t tot_num_kmers = 0;
  int64_t num_reads = 0;
  int64_t num_reads_aligned = 0;
  string rname, rseq, dummy;
  chrono::duration<double> get_ctgs_dt(0);
  {
    kmer_ht->open_alns_file(options.out_fname);
    barrier();
    ifstream reads_file(options.reads_fname);
    ProgressBar progbar(&reads_file, "Parsing reads file", options.one_file_per_rank);
    int64_t start_rank = options.one_file_per_rank ? 0 : rank_me();
    int64_t stop_rank = options.one_file_per_rank ? rank_n() : rank_me() + 1;
    auto start_offset = get_file_offset_for_rank(reads_file, start_rank);
    auto stop_offset =  get_file_offset_for_rank(reads_file, stop_rank);
    reads_file.seekg(start_offset);
    // only one alignment needed per read-ctg
    unordered_map<cid_t, int> cids_mapped;
    cids_mapped.reserve(1000);
    while (!reads_file.eof()) {
      for (auto buf : {&rname, &rseq, &dummy, &dummy}) 
        getline(reads_file, *buf);
      if (reads_file.tellg() > stop_offset) break;
      progbar.update();
      //adjust_read_name(rname);
      if (options.kmer_len > rseq.length()) continue;
      auto num_kmers = rseq.length() - options.kmer_len + 1;
      tot_num_kmers += num_kmers;
      // only one alignment needed per read-ctg
      cids_mapped.clear();
      promise<> compute_prom;
      // split read into kmers
      for (int i = 0; i < num_kmers; i += options.seed_space) {
        string&& kmer = rseq.substr(i, options.kmer_len); 
        string kmer_rc;
        bool is_rc = false;
        revcomp(kmer, &kmer_rc);
        if (cond_revcomp(kmer, &kmer_rc)) {
          kmer = kmer_rc;
          is_rc = true;
        }
        auto t = NOW();

        compute_prom.require_anonymous(1);
        kmer_ht->get_ctgs_with_kmer(kmer).then(
          [=, &cids_mapped, &compute_prom, &get_ctgs_dt](vector<CtgLoc> aligned_ctgs) {
            get_ctgs_dt += (NOW() - t);
            compute_alns_for_kmer(kmer_ht, aligned_ctgs, cids_mapped, i, is_rc, rname, rseq, compute_prom);
            compute_prom.fulfill_anonymous(1);
          });
      }
      compute_prom.finalize().wait();
      if (cids_mapped.size()) num_reads_aligned++;
      num_reads++;
    }
    progbar.done();
    barrier();
    kmer_ht->close_alns_file();
  }
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto tot_num_reads_aligned = reduce_one(num_reads_aligned, op_fast_add, 0).wait();
  SOUT("Parsed ", tot_num_reads, " reads, with ",
       reduce_one(tot_num_kmers, op_fast_add, 0).wait(), " kmers, and found ",
       reduce_one(kmer_ht->num_alns, op_fast_add, 0).wait(), " alignments, of which ",
       kmer_ht->get_num_perfect_alns(), " are perfect\n");
  SOUT("Reads mapped ", setprecision(2), fixed, 100.0 * tot_num_reads_aligned / tot_num_reads, "%\n");
  double av_ssw_secs = reduce_one(kmer_ht->ssw_dt.count(), op_fast_add, 0).wait() / rank_n();
  double max_ssw_secs = reduce_one(kmer_ht->ssw_dt.count(), op_fast_max, 0).wait();
  SOUT("Average SSW time ", setprecision(2), fixed, av_ssw_secs, " s, ",
       "max ", setprecision(2), fixed, max_ssw_secs, " s, ",
       "balance ", setprecision(2), fixed, av_ssw_secs / max_ssw_secs, "\n");
  double av_get_ctgs_secs = reduce_one(get_ctgs_dt.count(), op_fast_add, 0).wait() / rank_n();
  double max_get_ctgs_secs = reduce_one(get_ctgs_dt.count(), op_fast_max, 0).wait();
  SOUT("Average time to get ctgs ", setprecision(2), fixed, av_get_ctgs_secs, " s, ",
       "max ", max_get_ctgs_secs, " s, ",
       "balance ", av_get_ctgs_secs / max_get_ctgs_secs, "\n");
  */
}


void find_alignments(unsigned kmer_len, Contigs &ctgs)
{
  Timer(__func__);
  KmerHashTable kmer_ht(kmer_len);
  build_alignment_index(kmer_ht, kmer_len, ctgs);
  do_alignments(kmer_ht, kmer_len);
  barrier();
}

