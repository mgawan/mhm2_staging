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
#include "kmer.hpp"

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

class KmerCtgDHT {
private:
  // no need for a global ptr - just need to know where to find the sequence with an rpc to the destination
  // so the string ptr is only valid on the destination rank
  using kmer_map_t = unordered_map<Kmer, vector<CtgLoc> >;
  dist_object<kmer_map_t> kmer_map;
  int64_t num_perfect_alns;
  ofstream alns_file;
  IntermittentTimer t_ssw;

  // the key is to_string(cid) + rname
  unordered_map<string, bool> cid_reads_map;
  
  // default aligner and filter
  StripedSmithWaterman::Aligner ssw_aligner;
  StripedSmithWaterman::Filter ssw_filter;
  
  size_t get_target_rank(Kmer &kmer) {
    return std::hash<Kmer>{}(kmer) % rank_n();
  }

public:

  int kmer_len;
  int64_t num_alns;
  
  
  // aligner construction: match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty - no entry for ambiquity, which should be 2
  // note SSW internal defaults are 2 2 3 1
  KmerCtgDHT(int kmer_len) : kmer_map({}), num_perfect_alns(0), ssw_aligner(1, 3, 5, 2, 2), num_alns(0),
                             t_ssw("SSW"), kmer_len(kmer_len) {
    ssw_filter.report_cigar = false;
  }

  void clear() {
    for (auto it = kmer_map->begin(); it != kmer_map->end(); ) {
      it = kmer_map->erase(it);
    }
  }
  
  ~KmerCtgDHT() {
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

  void add_kmer(Kmer &kmer, const CtgLoc &ctg_loc, promise<> &prom) {
    rpc(get_target_rank(kmer), operation_cx::as_promise(prom),
        [](dist_object<kmer_map_t> &kmer_map, Kmer::MerArray merarr, CtgLoc ctg_loc) {
          Kmer kmer(merarr);
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
        }, kmer_map, kmer.get_array(), ctg_loc);
  }

  future<vector<CtgLoc> > get_ctgs_with_kmer(Kmer &kmer) {
    return rpc(get_target_rank(kmer),
               [](dist_object<kmer_map_t> &kmer_map, Kmer::MerArray merarr) -> vector<CtgLoc> {
                 Kmer kmer(merarr);
                 const auto it = kmer_map->find(kmer);
                 if (it == kmer_map->end()) return {};
                 return it->second;
               }, kmer_map, kmer.get_array());
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


static void build_alignment_index(KmerCtgDHT &kmer_ctg_dht, unsigned kmer_len, Contigs &ctgs) 
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
    
    auto kmers = Kmer::get_kmers(ctg->seq);
    tot_num_kmers += kmers.size();
    for (int i = 0; i < kmers.size(); i++) {
      Kmer kmer_rc = kmers[i].twin();
      if (kmer_rc < kmers[i]) {
        ctg_loc.is_rc = false;
        kmer_ctg_dht.add_kmer(kmers[i], ctg_loc, prom);
      } else {
        ctg_loc.is_rc = true;
        kmer_ctg_dht.add_kmer(kmer_rc, ctg_loc, prom);
      }
      ctg_loc.pos_in_ctg = i;
      progress();
    }
  }
  prom.finalize().wait();
  progbar.done();
  barrier();
  auto num_kmers_added = reduce_one(tot_num_kmers, op_fast_add, 0).wait();
  auto num_kmers_in_ht = kmer_ctg_dht.get_num_kmers();
  SOUT("Found ", num_kmers_added, " kmers in contigs\n");
  auto num_dropped = kmer_ctg_dht.get_num_dropped();
  if (num_dropped) {
    SOUT("Dropped ", num_dropped, " kmer-to-contig mappings (", 
         setprecision(2), fixed, (100.0 * num_dropped / num_kmers_added), "%)\n");
  }
}

/*
static void compute_alns_for_kmer(KmerCtgDHT *kmer_ctg_dht, vector<CtgLoc> &aligned_ctgs, unordered_map<cid_t, int> &cids_mapped,
                                  int pos_in_read, bool is_rc, const string &rname, const string &rseq, promise<> &prom) {
  string rseq_rc = "";
  int rlen = rseq.length();
  for (auto ctg_loc : aligned_ctgs) {
    progress();
    const auto it = cids_mapped.find(ctg_loc.cid);
    if (it == cids_mapped.end()) {
      cids_mapped.insert({ctg_loc.cid, pos_in_read});
      kmer_ctg_dht->num_alns++;
      char orient = '+';
      if (ctg_loc.is_rc != is_rc) {
        orient = '-';
        pos_in_read = rlen - (kmer_ctg_dht->kmer_len + pos_in_read);
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
               kmer_ctg_dht->align_read(rname, orient == '+' ? rseq : rseq_rc, seq_buf, rstart, start_pos, end_pos, ctg_loc, orient, subseq_len);
               delete[] seq_buf;
             });

    }
  }
}
*/

static void do_alignments(KmerCtgDHT &kmer_ctg_dht, unsigned kmer_len, unsigned seed_space)
{
  /*
  Timer timer(__func__);
  int64_t tot_num_kmers = 0;
  int64_t num_reads = 0;
  int64_t num_reads_aligned = 0;
  string rname, rseq, dummy;
  IntermittentTimer t_get_ctgs("get_ctgs");
  kmer_ctg_dht->open_alns_file(options.out_fname);
  barrier();
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    FastqReader fqr(merged_reads_fname, PER_RANK_FILE);
    string id, seq, quals;
    ProgressBar progbar(fqr.my_file_size(), "Aligning reads to contigs");
    size_t tot_bytes_read = 0;
    // only one alignment needed per read-ctg
    unordered_map<cid_t, int> cids_mapped;
    cids_mapped.reserve(1000);
    while (true) {
      size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      num_reads++;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      if (kmer_len > seq.length()) continue;
      auto num_kmers = seq.length() - kmer_len + 1;
      tot_num_kmers += num_kmers;
      // only one alignment needed per read-ctg
      cids_mapped.clear();
      promise<> compute_prom;
      // split read into kmers
      for (int i = 0; i < num_kmers; i += seed_space) {
        string&& kmer = seq.substr(i, kmer_len); 
        string kmer_rc;
        bool is_rc = false;
        revcomp(kmer, &kmer_rc);
        if (cond_revcomp(kmer, &kmer_rc)) {
          kmer = kmer_rc;
          is_rc = true;
        }
        auto t = NOW();

        compute_prom.require_anonymous(1);
        kmer_ctg_dht->get_ctgs_with_kmer(kmer).then(
          [=, &cids_mapped, &compute_prom, &get_ctgs_dt](vector<CtgLoc> aligned_ctgs) {
            get_ctgs_dt += (NOW() - t);
            compute_alns_for_kmer(kmer_ctg_dht, aligned_ctgs, cids_mapped, i, is_rc, rname, rseq, compute_prom);
            compute_prom.fulfill_anonymous(1);
          });
      }
      compute_prom.finalize().wait();
      if (cids_mapped.size()) num_reads_aligned++;
      num_reads++;
    }
    progbar.done();
    barrier();
    kmer_ctg_dht->close_alns_file();
  }
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto tot_num_reads_aligned = reduce_one(num_reads_aligned, op_fast_add, 0).wait();
  SOUT("Parsed ", tot_num_reads, " reads, with ",
       reduce_one(tot_num_kmers, op_fast_add, 0).wait(), " kmers, and found ",
       reduce_one(kmer_ctg_dht->num_alns, op_fast_add, 0).wait(), " alignments, of which ",
       kmer_ctg_dht->get_num_perfect_alns(), " are perfect\n");
  SOUT("Reads mapped ", setprecision(2), fixed, 100.0 * tot_num_reads_aligned / tot_num_reads, "%\n");
  double av_ssw_secs = reduce_one(kmer_ctg_dht->ssw_dt.count(), op_fast_add, 0).wait() / rank_n();
  double max_ssw_secs = reduce_one(kmer_ctg_dht->ssw_dt.count(), op_fast_max, 0).wait();
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


void find_alignments(unsigned kmer_len, unsigned seed_space, Contigs &ctgs)
{
  Timer timer(__func__);
  KmerCtgDHT kmer_ctg_dht(kmer_len);
  build_alignment_index(kmer_ctg_dht, kmer_len, ctgs);
  do_alignments(kmer_ctg_dht, kmer_len, seed_space);
  barrier();
}

