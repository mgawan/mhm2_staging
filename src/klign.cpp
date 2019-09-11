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
#include "aggr_store.hpp"
#include "zstr.hpp"
#include "fastq.hpp"


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

// global variables to avoid passing dist objs to rpcs
static int64_t _num_dropped = 0;
static int64_t _num_ctg_kmers_matched = 0;

class KmerCtgDHT {
private:

  struct MerarrAndCtgLoc {
    Kmer::MerArray merarr;
    CtgLoc ctg_loc;
  };

  using MerarrAndRead = tuple<Kmer::MerArray, string, string, int, bool>;
  
  // no need for a global ptr - just need to know where to find the sequence with an rpc to the destination
  // so the string ptr is only valid on the destination rank
  using kmer_map_t = unordered_map<Kmer, vector<CtgLoc> >;
  dist_object<kmer_map_t> kmer_map;
    
  AggrStore<MerarrAndCtgLoc> kmer_store;
  AggrStore<MerarrAndRead> try_align_store;
  
  int64_t num_perfect_alns;
  zstr::ofstream *alns_file;
  IntermittentTimer t_ssw;

  // the key is to_string(cid) + rname
  unordered_map<string, bool> cid_reads_map;
  
  // default aligner and filter
  StripedSmithWaterman::Aligner ssw_aligner;
  StripedSmithWaterman::Filter ssw_filter;
  
  size_t get_target_rank(Kmer &kmer) {
    return std::hash<Kmer>{}(kmer) % rank_n();
  }

  struct InsertKmer {
    void operator()(MerarrAndCtgLoc &merarr_and_ctg_loc, dist_object<kmer_map_t> &kmer_map) {
      Kmer kmer(merarr_and_ctg_loc.merarr);
      CtgLoc ctg_loc = merarr_and_ctg_loc.ctg_loc;
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
    }
  };
  dist_object<InsertKmer> insert_kmer;

  struct FindSeedAndAlign {
    void operator()(MerarrAndRead &merarr_and_read, dist_object<kmer_map_t> &kmer_map) {
      Kmer::MerArray merarr;
      string id, seq;
      int pos;
      bool is_rc;
      tie(merarr, id, seq, pos, is_rc) = merarr_and_read;
      Kmer kmer(merarr);
      const auto it = kmer_map->find(kmer);
      if (it != kmer_map->end()) _num_ctg_kmers_matched++;
    }
  };
  dist_object<FindSeedAndAlign> find_seed_and_align;
  
public:

  int kmer_len;
  int64_t num_alns;
  
  // aligner construction: match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty - no entry for ambiquity, which should be 2
  // note SSW internal defaults are 2 2 3 1
  KmerCtgDHT(int kmer_len, int max_store_size) : kmer_map({}), kmer_store({}), insert_kmer({}), try_align_store({}),
                                                 find_seed_and_align({}), num_perfect_alns(0), ssw_aligner(1, 3, 5, 2, 2),
                                                 num_alns(0), t_ssw("SSW"), kmer_len(kmer_len) {
    ssw_filter.report_cigar = false;
    kmer_store.set_size("insert ctg seeds", max_store_size);
    try_align_store.set_size("check read alns", max_store_size);
  }

  void clear() {
    for (auto it = kmer_map->begin(); it != kmer_map->end(); ) {
      it = kmer_map->erase(it);
    }
    kmer_store.clear();
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

  int64_t get_num_alns(bool all = false) {
    if (!all) return reduce_one(num_alns, op_fast_add, 0).wait();
    else return reduce_all(num_alns, op_fast_add).wait();
  }

  int64_t get_num_dropped(bool all = false) {
    if (!all) return reduce_one(_num_dropped, op_fast_add, 0).wait();
    else return reduce_all(_num_dropped, op_fast_add).wait();
  }

  int64_t get_num_ctg_kmers_matched(bool all = false) {
    if (!all) return reduce_one(_num_ctg_kmers_matched, op_fast_add, 0).wait();
    else return reduce_all(_num_ctg_kmers_matched, op_fast_add).wait();
  }
  
  void open_alns_file(zstr::ofstream *f) {
    alns_file = f;
  }

  void close_alns_file() {
    alns_file->close();
  }

  void add_kmer(Kmer kmer, CtgLoc &ctg_loc) {
    Kmer kmer_rc = kmer.revcomp();
    ctg_loc.is_rc = false;
    if (kmer_rc < kmer) {
      kmer = kmer_rc;
      ctg_loc.is_rc = true;
    }
    MerarrAndCtgLoc merarr_and_ctg_loc = { kmer.get_array(), ctg_loc };
    kmer_store.update(get_target_rank(kmer), merarr_and_ctg_loc, insert_kmer, kmer_map);
  }

  void flush_add_kmers() {
    Timer timer(__func__);
    barrier();
    kmer_store.flush_updates(insert_kmer, kmer_map);
    barrier();
  }

  void try_align_read(Kmer kmer, string &id, string &seq, int pos) {
    Kmer kmer_rc = kmer.revcomp();
    bool is_rc = false;
    if (kmer_rc < kmer) {
      kmer = kmer_rc;
      is_rc = true;
    }
    MerarrAndRead merarr_and_read = { kmer.get_array(), id, seq, pos, is_rc };
    try_align_store.update(get_target_rank(kmer), merarr_and_read, find_seed_and_align, kmer_map);
  }

  void flush_try_align_read() {
    Timer timer(__func__);
    barrier();
    try_align_store.flush_updates(find_seed_and_align, kmer_map);
    barrier();
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
      *alns_file << "MERALIGNER-0\t" << rname << "\t" << rstart + 1 << "\t" << rend << "\t" << rseq.length() << "\t"
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
      *alns_file << "MERALIGNER-1\t" << rname << "\t" << rstart + 1 << "\t" << rend << "\t" << rseq.length() << "\t"
                 << "Contig" + to_string(ctg_loc.cid) << "\t"
                 << start_pos + aln.ref_begin + 1 << "\t" << start_pos + aln.ref_end + 1 << "\t" << ctg_loc.clen << "\t"
                 << (orient == '+' ? "Plus" : "Minus") << "\t0\t0\t0\t0\t" << aln.sw_score << "\t" << aln.sw_score_next_best << endl;
    }
  }

  // this is really only for debugging
  void dump_ctg_kmers(int kmer_len) {
    Timer timer(__func__);
    string dump_fname = "ctg_kmers-" + to_string(kmer_len) + ".txt.gz";
    get_rank_path(dump_fname, rank_me());
    zstr::ofstream dump_file(dump_fname);
    ostringstream out_buf;
    ProgressBar progbar(kmer_map->size(), "Dumping kmers to " + dump_fname);
    int64_t i = 0;
    for (auto &elem : *kmer_map) {
      auto ctg_loc = &elem.second[0];
      out_buf << elem.first << " " << elem.second.size() << " " << ctg_loc->clen << " " << ctg_loc->pos_in_ctg << " "
              << ctg_loc->is_rc << endl;
      i++;
      if (!(i % 1000)) {
        dump_file << out_buf.str();
        out_buf = ostringstream();
      }
      progbar.update();
    }
    if (!out_buf.str().empty()) dump_file << out_buf.str();
    dump_file.close();
    progbar.done();
    SOUT("Dumped ", this->get_num_kmers(), " kmers\n");
  }

};


static void build_alignment_index(KmerCtgDHT &kmer_ctg_dht, unsigned kmer_len, Contigs &ctgs) 
{
  Timer timer(__func__);
  int64_t num_kmers = 0;
  ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    global_ptr<char> seq_gptr = allocate<char>(ctg->seq.length() + 1);
    strcpy(seq_gptr.local(), ctg->seq.c_str());
    CtgLoc ctg_loc = { .cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length() };
    auto kmers = Kmer::get_kmers(ctg->seq);
    num_kmers += kmers.size();
    for (int i = 0; i < kmers.size(); i++) {
      ctg_loc.pos_in_ctg = i;
      kmer_ctg_dht.add_kmer(kmers[i], ctg_loc);
      progress();
    }
  }
  progbar.done();
  kmer_ctg_dht.flush_add_kmers();
  barrier();
  auto tot_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto num_kmers_in_ht = kmer_ctg_dht.get_num_kmers();
  SOUT("Processed ", tot_num_kmers, " seeds from contigs, added ", num_kmers_in_ht, "\n");
  auto num_dropped = kmer_ctg_dht.get_num_dropped();
  if (num_dropped) {
    SOUT("Dropped ", num_dropped, " seed-to-contig mappings (", 
         setprecision(2), fixed, (100.0 * num_dropped / tot_num_kmers), "%)\n");
  }
}

/*
static void compute_alns_for_kmer(KmerCtgDHT *kmer_ctg_dht, vector<CtgLoc> &aligned_ctgs, unordered_map<cid_t, int> &cids_mapped,
                                  int pos_in_read, bool is_rc, const string &rname, const string &rseq) {
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

static void do_alignments(KmerCtgDHT &kmer_ctg_dht, unsigned kmer_len, unsigned seed_space, vector<string> &reads_fname_list)
{
  Timer timer(__func__);
  int64_t tot_num_kmers = 0;
  int64_t num_reads = 0;
  int64_t num_reads_aligned = 0;
  IntermittentTimer t_get_ctgs("get_ctgs");
  string dump_fname = "klign-" + to_string(kmer_len) + ".alns.gz";
  get_rank_path(dump_fname, rank_me());
  zstr::ofstream dump_file(dump_fname);
  kmer_ctg_dht.open_alns_file(&dump_file);
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
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      if (kmer_len > seq.length()) continue;
      auto kmers = Kmer::get_kmers(seq);
      // only one alignment needed per read-ctg
      cids_mapped.clear();
      for (int i = 0; i < kmers.size(); i += seed_space) {
        tot_num_kmers++;
        kmer_ctg_dht.try_align_read(kmers[i], id, seq, i);
        progress();
      }
      if (cids_mapped.size()) num_reads_aligned++;
      num_reads++;
    }
    progbar.done();
    kmer_ctg_dht.flush_try_align_read();
    barrier();
    kmer_ctg_dht.close_alns_file();
  }
  SOUT("Found ", perc_str(kmer_ctg_dht.get_num_ctg_kmers_matched(), kmer_ctg_dht.get_num_kmers()), " seeds in contigs\n");
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto num_alns = kmer_ctg_dht.get_num_alns();
  SOUT("Parsed ", tot_num_reads, " reads, with ", reduce_one(tot_num_kmers, op_fast_add, 0).wait(), " seeds, and found ",
       num_alns, " alignments, of which ", perc_str(kmer_ctg_dht.get_num_perfect_alns(), num_alns), " are perfect\n");
  SOUT("Mapped ", perc_str(reduce_one(num_reads_aligned, op_fast_add, 0).wait(), tot_num_reads), " reads to contigs\n");
}


void find_alignments(unsigned kmer_len, unsigned seed_space, vector<string> &reads_fname_list, int max_store_size, Contigs &ctgs)
{
  Timer timer(__func__);
  _num_ctg_kmers_matched = 0;
  _num_dropped = 0;
  KmerCtgDHT kmer_ctg_dht(kmer_len, max_store_size);
  build_alignment_index(kmer_ctg_dht, kmer_len, ctgs);
#ifdef DEBUG
  kmer_ctg_dht.dump_ctg_kmers(kmer_len);
#endif
  do_alignments(kmer_ctg_dht, kmer_len, seed_space, reads_fname_list);
  barrier();
}

