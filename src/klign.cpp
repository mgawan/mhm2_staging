/*
  Sequence alignment similar to meraligner.
  c shofmeyr@lbl.gov
  Nov 2018
*/


#include <iostream>
#include <algorithm>
#include <chrono>
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
#include "alignments.hpp"

using namespace std;
using namespace upcxx;

#define NOW std::chrono::high_resolution_clock::now

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

static const int ALN_EXTRA = 10;

using aligned_ctgs_map_t = unordered_map<cid_t, tuple<int, bool, CtgLoc>>;


class KmerCtgDHT {

  struct MerarrAndCtgLoc {
    MerArray merarr;
    CtgLoc ctg_loc;
  };

  AggrStore<MerarrAndCtgLoc> kmer_store;
  
  using kmer_map_t = unordered_map<Kmer, vector<CtgLoc> >;
  dist_object<kmer_map_t> kmer_map;
#ifdef DEBUG  
  zstr::ofstream *alns_file;
#endif  
  int64_t num_alns;
  int64_t num_perfect_alns;
  chrono::duration<double> ssw_dt;
  int max_ctg_seq_cache_size;
  unordered_map<cid_t, string> ctg_seq_cache;
  int num_ctg_seq_cache_hits;
  int64_t ctg_seq_bytes_fetched;
  int64_t sum_hamming_dist;
  
  Alns *alns;
  
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
          assert(ctg_locs->size() <= MAX_KMER_MAPPINGS);
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

  void align_read(const string &rname, int64_t cid, const string &rseq, const string &cseq, int rstart, int rend, int rlen,
                  int cstart, int cend, int clen, char orient, int start_overlap, int end_overlap, 
                  const string &full_rseq, const string &full_cseq, int pr, int pc, bool ctg_rc) {
    num_alns++;
    int hdist = 0;
    Aln aln;
    StripedSmithWaterman::Alignment ssw_aln;

    if (rseq == cseq) {
      num_perfect_alns++;
      int aln_len = rseq.length();
      ssw_aln.ref_begin = 0;
      ssw_aln.ref_end = aln_len;
      ssw_aln.query_begin = 0;
      ssw_aln.query_end = aln_len;
      ssw_aln.sw_score = aln_len;
      ssw_aln.sw_score_next_best = 0;
    } else {
      // make sure upcxx progress is done before starting alignment
      discharge();
#ifdef DEBUG
      // sanity check that kmer matches
      for (int i = 0; i < rseq.length(); i++) {
        // allow matches to any Ns
        if (rseq[i] == 'N') continue;
        if (rseq[i] != cseq[i]) {
          // we have a mismatch
          if (i >= start_overlap && i < start_overlap + kmer_len) {
            // this is in the seed region - there should be no mismatches
            //if (full_cseq.length() < 500) 
              DIE("Mismatch found within seed at position ", i, ": ", rseq[i], " != ", cseq[i],
                  " start overlap ", start_overlap, " end overlap ", end_overlap, " orient ", orient,
                  " rlen ", rlen, " clen ", clen, " pr ", pr, " pc ", pc, " ctg rc ", ctg_rc, " rid ", rname, " cid ", cid,
                  "\n", rseq, "\n", cseq, "\n", full_rseq, "\n", full_cseq);
          }
        }
      }
#endif
      // contig is the ref, read is the query - done this way so that we can do multiple alns to each read
      auto t = NOW();
      ssw_aligner.Align(cseq.c_str(), rseq.c_str(), rseq.length(), ssw_filter, &ssw_aln, max((int)(rseq.length() / 2), 15));
      ssw_dt += (NOW() - t);
    }
    rend = rstart + ssw_aln.ref_end + 1;
    rstart = rstart + ssw_aln.ref_begin;
    if (orient == '-') {
      int tmp = rstart;
      rstart = rlen - rend;
      rend = rlen - tmp;
    }
    aln = { .read_id = rname, .cid = cid,
            .rstart = rstart, .rstop = rend, .rlen = rlen,
            .cstart = cstart + ssw_aln.query_begin, .cstop = cstart + ssw_aln.query_end + 1, .clen = clen,
            .orient = orient, .score1 = ssw_aln.sw_score, .score2 = ssw_aln.sw_score_next_best };
    int aln_len = rend - rstart;
    // check against hamming distance
    hdist = hamming_dist(rseq.substr(ssw_aln.ref_begin, aln_len), cseq.substr(ssw_aln.query_begin, aln_len), false);
    sum_hamming_dist += hdist;
#ifdef DEBUG
    *alns_file << "MERALIGNER\t" << aln.to_string() << endl;
#endif
    // adjust for use with cgraph algos
    if (aln.orient == '-') {
      int tmp = aln.cstart;
      aln.cstart = aln.clen - aln.cstop;
      aln.cstop = aln.clen - tmp;
    }
    alns->add_aln(aln);
  }
  
public:

  int kmer_len;  
  
  // aligner construction: match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty - no entry for ambiquity, which should be 2
  // note SSW internal defaults are 2 2 3 1
  KmerCtgDHT(int kmer_len, int max_store_size, int max_ctg_cache, Alns *alns)
    : kmer_map({})
    , kmer_store({})
    , insert_kmer({})
    , ssw_aligner(1, 3, 5, 2, 2)
    , num_alns(0)
    , num_perfect_alns(0)
    , ssw_dt(0)
    , kmer_len(kmer_len)
    , num_ctg_seq_cache_hits(0)
    , ctg_seq_bytes_fetched(0)
    , max_ctg_seq_cache_size(max_ctg_cache)
    , alns(alns)
    , sum_hamming_dist(0) {
    
    ssw_filter.report_cigar = false;
    kmer_store.set_size("insert ctg seeds", max_store_size);
    if (max_ctg_seq_cache_size) ctg_seq_cache.reserve(max_ctg_cache);

#ifdef DEBUG
    string dump_fname = "klign-" + to_string(kmer_len) + ".alns.gz";
    get_rank_path(dump_fname, rank_me());
    alns_file = new zstr::ofstream(dump_fname);
#endif
  }

  void clear() {
    for (auto it = kmer_map->begin(); it != kmer_map->end(); ) {
      it = kmer_map->erase(it);
    }
    kmer_store.clear();
  }
  
  ~KmerCtgDHT() {
    clear();
#ifdef DEBUG
    alns_file->close();
    delete alns_file;
#endif
  }

  int64_t get_num_kmers(bool all = false) {
    if (!all) return reduce_one(kmer_map->size(), op_fast_add, 0).wait();
    else return reduce_all(kmer_map->size(), op_fast_add).wait();
  }
  
  int64_t get_num_perfect_alns(bool all = false) {
    if (!all) return reduce_one(num_perfect_alns, op_fast_add, 0).wait();
    else return reduce_all(num_perfect_alns, op_fast_add).wait();
  }

  int64_t get_my_num_alns() {
    return num_alns;
  }
  
  int64_t get_num_alns(bool all = false) {
    if (!all) return reduce_one(num_alns, op_fast_add, 0).wait();
    else return reduce_all(num_alns, op_fast_add).wait();
  }

  int64_t get_num_dropped(bool all = false) {
    if (!all) return reduce_one(_num_dropped, op_fast_add, 0).wait();
    else return reduce_all(_num_dropped, op_fast_add).wait();
  }

  double get_av_ssw_secs() {
    return reduce_one(ssw_dt.count(), op_fast_add, 0).wait() / rank_n();
  }
  
  double get_max_ssw_secs() {
    return reduce_one(ssw_dt.count(), op_fast_max, 0).wait();
  }

  int64_t get_ctg_seq_cache_hits() {
    return reduce_one(num_ctg_seq_cache_hits, op_fast_add, 0).wait();
  }

  int64_t get_ctg_seq_bytes_fetched() {
    return reduce_one(ctg_seq_bytes_fetched, op_fast_add, 0).wait();
  }    

  int64_t get_sum_hamming_dist() {
    return reduce_one(sum_hamming_dist, op_fast_add, 0).wait();
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

  future<vector<CtgLoc> > get_ctgs_with_kmer(Kmer &kmer) {
    return rpc(get_target_rank(kmer),
               [](MerArray merarr, dist_object<kmer_map_t> &kmer_map) -> vector<CtgLoc> {
                 Kmer kmer(merarr);
                 const auto it = kmer_map->find(kmer);
                 if (it == kmer_map->end()) return {};
                 return it->second;
               }, kmer.get_array(), kmer_map);
  }

  // this is really only for debugging
  void dump_ctg_kmers() {
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
    SLOG_VERBOSE("Dumped ", this->get_num_kmers(), " kmers\n");
  }

  future<> compute_alns_for_read(aligned_ctgs_map_t *aligned_ctgs_map, const string &rname, string rseq) {
    int rlen = rseq.length();
    string rseq_rc = revcomp(rseq);
    future<> all_fut = make_future();
    for (auto &elem : *aligned_ctgs_map) {
      progress();
      int pos_in_read;
      bool read_kmer_is_rc;
      CtgLoc ctg_loc;
      tie(pos_in_read, read_kmer_is_rc, ctg_loc) = elem.second;
      char orient = '+';
      string *rseq_ptr = &rseq;
      if (ctg_loc.is_rc != read_kmer_is_rc) {
        // it's revcomp in either contig or read, but not in both or neither
        orient = '-';
        pos_in_read = rlen - (kmer_len + pos_in_read);
        rseq_ptr = &rseq_rc;
      }
      // assume the alignment is anchored to the seed and no indels (only for illumina) 
      // so look for max possible overlap
      int start_overlap = min(pos_in_read, ctg_loc.pos_in_ctg);
      int end_overlap = min(rlen - pos_in_read - kmer_len, ctg_loc.clen - ctg_loc.pos_in_ctg - kmer_len);
      int rstart = pos_in_read - start_overlap;
      int rend = pos_in_read + kmer_len + end_overlap;
      int cstart = ctg_loc.pos_in_ctg - start_overlap;
      int cend = ctg_loc.pos_in_ctg + kmer_len + end_overlap;
      int overlap_len = cend - cstart;
      assert(rstart >= 0);
      assert(rend <= rlen);
      assert(cstart >= 0);
      assert(cend <= ctg_loc.clen);
      assert(overlap_len == rend - rstart);
      string read_subseq = rseq_ptr->substr(rstart, overlap_len);
      bool ctg_in_cache = false;
      // if max_ctg_seq_cache_size is > 0, then we are using a cache
      if (max_ctg_seq_cache_size) {
        auto it = ctg_seq_cache.find(ctg_loc.cid);
        if (it != ctg_seq_cache.end()) {
          num_ctg_seq_cache_hits++;
          // only get the subsequence needed
          string ctg_subseq = it->second.substr(cstart, overlap_len);
          // align only the subseqs
          align_read(rname, ctg_loc.cid, read_subseq, ctg_subseq, rstart, rend, rlen, cstart, cend, ctg_loc.clen, orient,
                     start_overlap, end_overlap, *rseq_ptr, it->second, pos_in_read, ctg_loc.pos_in_ctg, ctg_loc.is_rc);
          ctg_in_cache = true;
        }
      }
      if (!ctg_in_cache) {
        int fetch_start_pos, fetch_end_pos;
        if (!max_ctg_seq_cache_size) {
          // not caching - only fetch subseq
          fetch_start_pos = cstart;
          fetch_end_pos = cend;
        } else {
          // caching - fetch whole contig
          fetch_start_pos = 0;
          fetch_end_pos = ctg_loc.clen;
        }
        int fetch_seq_len = fetch_end_pos - fetch_start_pos;
        // space for whatever is getting fetched plus a null terminator
        char *seq_buf = new char[fetch_seq_len + 1];
        // get the full or subseq
        auto fut = rget(ctg_loc.seq_gptr + fetch_start_pos, seq_buf, fetch_seq_len).then(
          [=]() {
            string ctg_seq(seq_buf, fetch_seq_len);
            delete[] seq_buf;
            ctg_seq_bytes_fetched += fetch_seq_len;
            if (ctg_seq_cache.size() < max_ctg_seq_cache_size) {
              // if cache is in use and there is space, stash the whole contig
              ctg_seq_cache[ctg_loc.cid] = ctg_seq;
              assert(ctg_seq_cache[ctg_loc.cid].length() == ctg_loc.clen);
            }
            // get the subseq for alignment
            string ctg_subseq;
            if (!max_ctg_seq_cache_size) {
              // if we are not caching, we fetched just the substring anyway
              ctg_subseq = ctg_seq;
            } else {
              // we had to fetch the whole seq, so get just the substring
              ctg_subseq = ctg_seq.substr(cstart, overlap_len);
            }
            align_read(rname, ctg_loc.cid, read_subseq, ctg_subseq, rstart, rend, rlen, cstart, cend, ctg_loc.clen, orient,
                       start_overlap, end_overlap, *rseq_ptr, ctg_seq, pos_in_read, ctg_loc.pos_in_ctg, ctg_loc.is_rc);
          });
        all_fut = when_all(all_fut, fut);
      }
    }
    return all_fut;
  }
  
};


static void build_alignment_index(KmerCtgDHT &kmer_ctg_dht, Contigs &ctgs) {
  Timer timer(__func__);
  int64_t num_kmers = 0;
  ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    global_ptr<char> seq_gptr = allocate<char>(ctg->seq.length() + 1);
    strcpy(seq_gptr.local(), ctg->seq.c_str());
    CtgLoc ctg_loc = { .cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length() };
    auto kmers = Kmer::get_kmers(kmer_ctg_dht.kmer_len, ctg->seq);
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
  SLOG_VERBOSE("Processed ", tot_num_kmers, " seeds from contigs, added ", num_kmers_in_ht, "\n");
  auto num_dropped = kmer_ctg_dht.get_num_dropped();
  if (num_dropped) {
    SLOG_VERBOSE("Dropped ", num_dropped, " seed-to-contig mappings (", 
                 setprecision(2), fixed, (100.0 * num_dropped / tot_num_kmers), "%)\n");
  }
}

static void do_alignments(KmerCtgDHT *kmer_ctg_dht, unsigned seed_space, vector<string> &reads_fname_list) {
  Timer timer(__func__);
  int64_t tot_num_kmers = 0;
  int64_t num_reads = 0;
  int64_t num_reads_aligned = 0;
  barrier();
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    FastqReader fqr(merged_reads_fname, PER_RANK_FILE);
    string read_id, read_seq, quals;
    ProgressBar progbar(fqr.my_file_size(), "Aligning reads to contigs");
    size_t tot_bytes_read = 0;
    // keep track of all futures to enable asynchrony at the file level
    vector<future<> > reads_futures;
    while (true) {
      size_t bytes_read = fqr.get_next_fq_record(read_id, read_seq, quals);
      if (!bytes_read) break;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      if (kmer_ctg_dht->kmer_len > read_seq.length()) continue;
      // accumulate all the ctgs that align to this read at any position in a hash table to filter out duplicates
      // this is dynamically allocated and deleted in the final .then callback when the computation is over
      // a mapping of cid to {pos in read, is_rc, ctg location)
      auto aligned_ctgs_map = new aligned_ctgs_map_t();
      auto kmers = Kmer::get_kmers(kmer_ctg_dht->kmer_len, read_seq);
      tot_num_kmers += kmers.size();
      // get all the seeds/kmers for a read, and add all the potential ctgs for aln to the aligned_ctgs_map
      // when the read future chain is completed, all the ctg info will be collected and the alignment can happen
      // FIXME: if using a seed > 1 and no alns are found, repeat with a seed of 1
      future<> read_fut_chain = make_future();
      for (int i = 0; i < kmers.size(); i += seed_space) {
        Kmer kmer = kmers[i];
        Kmer kmer_rc = kmer.revcomp();
        bool is_rc = false;
        if (kmer_rc < kmer) {
          kmer = kmer_rc;
          is_rc = true;
        }
        auto t = NOW();
        // add the fetched ctg into the unordered map
        auto fut = kmer_ctg_dht->get_ctgs_with_kmer(kmer).then(
          [=](vector<CtgLoc> aligned_ctgs) {
            for (auto ctg_loc : aligned_ctgs) {
              // ensures only the first kmer to cid mapping is retained - copies will not be inserted
              aligned_ctgs_map->insert({ctg_loc.cid, {i, is_rc, ctg_loc}});
            }
          });
        read_fut_chain = when_all(read_fut_chain, fut);
        progress();
      }
      // when all the ctgs are fetched, do the alignments
      auto fut = read_fut_chain.then(
        [=, &num_reads_aligned]() {
          if (aligned_ctgs_map->size()) num_reads_aligned++;
          return kmer_ctg_dht->compute_alns_for_read(aligned_ctgs_map, read_id, read_seq).then(
            [aligned_ctgs_map]() {
              delete aligned_ctgs_map;
            });
        });
      reads_futures.push_back(fut);
      //fut.wait();
      num_reads++;
    }
    //for (auto fut : reads_futures) fut.wait();
    progbar.done();
    barrier();
  }
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Parsed ", tot_num_reads, " reads, with ", reduce_one(tot_num_kmers, op_fast_add, 0).wait(), " seeds\n");
  auto tot_num_alns = kmer_ctg_dht->get_num_alns();
  SLOG_VERBOSE("Found ", tot_num_alns, " alignments, of which ", perc_str(kmer_ctg_dht->get_num_perfect_alns(), tot_num_alns),
               " are perfect\n");
  auto tot_num_reads_aligned = reduce_one(num_reads_aligned, op_fast_add, 0).wait();
  SLOG("Mapped ", perc_str(tot_num_reads_aligned, tot_num_reads), " reads to contigs\n");
  SLOG_VERBOSE("Average mappings per read ", (double)tot_num_alns / tot_num_reads_aligned, "\n");

  SLOG_VERBOSE("Ctg cache hits ", perc_str(kmer_ctg_dht->get_ctg_seq_cache_hits(), tot_num_alns), "\n");
  SLOG_VERBOSE("Fetched ", get_size_str(kmer_ctg_dht->get_ctg_seq_bytes_fetched()), " of contig sequences\n");
  
  double av_ssw_secs = kmer_ctg_dht->get_av_ssw_secs();
  double max_ssw_secs = kmer_ctg_dht->get_max_ssw_secs();
  SLOG_VERBOSE("Average SSW time ", setprecision(2), fixed, av_ssw_secs, " s, ",
               "max ", setprecision(2), fixed, max_ssw_secs, " s, ",
               "balance ", setprecision(2), fixed, av_ssw_secs / max_ssw_secs, "\n");
}

void find_alignments(unsigned kmer_len, unsigned seed_space, vector<string> &reads_fname_list, int max_store_size,
                     int max_ctg_cache, Contigs &ctgs, Alns *alns) {
  Timer timer(__func__, true);
  _num_dropped = 0;
  SLOG("Aligning with seed length ", kmer_len, " and seed space ", seed_space, "\n");
  //_get_ctgs_dt = std::chrono::duration<double>(0);
  KmerCtgDHT kmer_ctg_dht(kmer_len, max_store_size, max_ctg_cache, alns);
  barrier();
  build_alignment_index(kmer_ctg_dht, ctgs);
#ifdef DEBUG
  kmer_ctg_dht.dump_ctg_kmers();
#endif
  do_alignments(&kmer_ctg_dht, seed_space, reads_fname_list);
  assert(kmer_ctg_dht.get_my_num_alns() == alns->size());
  barrier();
  auto sum_hamming_dist = kmer_ctg_dht.get_sum_hamming_dist();
  auto num_alns = kmer_ctg_dht.get_num_alns();
  SOUT("Average hamming distance in alignments ", sum_hamming_dist / num_alns, "\n");
  SOUT("Number of duplicate alignments ", perc_str(alns->get_num_dups(), num_alns), "\n");
  barrier();
}

