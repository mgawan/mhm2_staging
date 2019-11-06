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

#define DUMP_ALNS

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

using aligned_ctgs_map_t = unordered_map<cid_t, tuple<int, bool, CtgLoc>>;


class KmerCtgDHT {

  struct MerarrAndCtgLoc {
    MerArray merarr;
    CtgLoc ctg_loc;
  };

  AggrStore<MerarrAndCtgLoc> kmer_store;
  
  using kmer_map_t = unordered_map<Kmer, vector<CtgLoc> >;
  dist_object<kmer_map_t> kmer_map;
#ifdef DUMP_ALNS
  zstr::ofstream *alns_file;
#endif  
  int64_t num_alns;
  int64_t num_perfect_alns;
  int64_t num_excess_alns_reads;
  int64_t num_overlaps;
  
  chrono::duration<double> ssw_dt;
  int max_ctg_seq_cache_size;
  unordered_map<cid_t, string> ctg_seq_cache;
  int num_ctg_seq_cache_hits;
  int64_t ctg_seq_bytes_fetched;
  
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

  void align_read(const string &rname, int64_t cid, const string &rseq, const string &cseq, int rstart, int rlen,
                  int cstart, int clen, char orient, Alns *read_alns, int extra_offset) {
    Aln aln;
    StripedSmithWaterman::Alignment ssw_aln;

    // use hamming distance for checking match first - works for > 90% of matches and reduces time spent doing SSW
    //int hdist = hamming_dist(rseq, cseq, true);
    //if (hdist < 3) {
    //if (rseq == cseq) {
    if (cseq.compare(extra_offset, rseq.length(), rseq) == 0) {
      num_perfect_alns++;
      int aln_len = rseq.length();
      ssw_aln.ref_begin = 0;
      ssw_aln.ref_end = aln_len - 1;
      ssw_aln.query_begin = extra_offset;
      ssw_aln.query_end = extra_offset + aln_len - 1;
      // the mismatch penalty is 3 
      ssw_aln.sw_score = aln_len;// - 3 * hdist;
      ssw_aln.sw_score_next_best = 0;
    } else {
      // make sure upcxx progress is done before starting alignment
      discharge();
      // contig is the ref, read is the query - done this way so that we can potentially do multiple alns to each read
      auto t = NOW();
      ssw_aligner.Align(cseq.c_str(), rseq.c_str(), rseq.length(), ssw_filter, &ssw_aln, max((int)(rseq.length() / 2), 15));
      ssw_dt += (NOW() - t);
    }

    int rstop = rstart + ssw_aln.ref_end + 1;
    rstart = rstart + ssw_aln.ref_begin;
    int cstop = cstart + ssw_aln.query_end + 1;
    cstart = cstart + ssw_aln.query_begin;

    if (orient == '-') switch_orient(rstart, rstop, rlen);
    aln = { .read_id = rname, .cid = cid,
            .rstart = rstart, .rstop = rstop, .rlen = rlen,
            .cstart = cstart, .cstop = cstop, .clen = clen,
            .orient = orient, .score1 = ssw_aln.sw_score, .score2 = ssw_aln.sw_score_next_best };
    int aln_len = rstop - rstart;
#ifdef DUMP_ALNS
    *alns_file << "MERALIGNER\t" << aln.to_string() << endl;
#endif
          
    // FIXME: check to see if this alignment is a duplicate of an existing contig alignment. It's a duplicate
    // if another alignment aligns to exactly the same position in the contig.
    
    read_alns->add_aln(aln);
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
    , num_excess_alns_reads(0)
    , num_overlaps(0)
    , ssw_dt(0)
    , kmer_len(kmer_len)
    , num_ctg_seq_cache_hits(0)
    , ctg_seq_bytes_fetched(0)
    , max_ctg_seq_cache_size(max_ctg_cache)
    , alns(alns) {
    
    ssw_filter.report_cigar = false;
    kmer_store.set_size("insert ctg seeds", max_store_size);
    if (max_ctg_seq_cache_size) ctg_seq_cache.reserve(max_ctg_cache);

#ifdef DUMP_ALNS
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
#ifdef DUMP_ALNS
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

  int64_t get_num_alns(bool all = false) {
    if (!all) return reduce_one(num_alns, op_fast_add, 0).wait();
    else return reduce_all(num_alns, op_fast_add).wait();
  }

  int64_t get_num_excess_alns_reads(bool all = false) {
    if (!all) return reduce_one(num_excess_alns_reads, op_fast_add, 0).wait();
    else return reduce_all(num_excess_alns_reads, op_fast_add).wait();
  }

  int64_t get_num_overlaps(bool all = false) {
    if (!all) return reduce_one(num_overlaps, op_fast_add, 0).wait();
    else return reduce_all(num_overlaps, op_fast_add).wait();
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
  
  void compute_alns_for_read(aligned_ctgs_map_t *aligned_ctgs_map, const string &rname, string rseq) {
    int rlen = rseq.length();
    string rseq_rc = revcomp(rseq);
    future<> all_fut = make_future();
    auto read_alns = new Alns();
    for (auto &elem : *aligned_ctgs_map) {
      progress();
      // drop alns for reads with too many alns
      if (read_alns->size() >= MAX_ALNS_PER_READ) continue;
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

      // FIXME: keep track of all alignments to each contig, and for each seed, look to see if there's already an alignment to
      // the contig that covers that seed. If so, don't align it.

      // calculate available bases before and after the seeded kmer
      int ctg_bases_left_of_kmer = ctg_loc.pos_in_ctg;
      int ctg_bases_right_of_kmer = ctg_loc.clen - ctg_bases_left_of_kmer - kmer_len;
      int read_bases_left_of_kmer = pos_in_read;
      int read_bases_right_of_kmer = rlen - kmer_len - pos_in_read;
      int left_of_kmer = min(read_bases_left_of_kmer, ctg_bases_left_of_kmer);
      int right_of_kmer = min(read_bases_right_of_kmer, ctg_bases_right_of_kmer);
      
      int cstart = ctg_loc.pos_in_ctg - left_of_kmer;
      int rstart = pos_in_read - left_of_kmer;
      int overlap_len = left_of_kmer + kmer_len + right_of_kmer;

      // add a few extra on either end if possible
      int ctg_aln_len = overlap_len + min(ctg_loc.clen - (cstart + overlap_len), ALIGN_EXPAND_BASES);
      int extra_offset = min(cstart, ALIGN_EXPAND_BASES);
      cstart -= extra_offset;
      ctg_aln_len += extra_offset;
      
      string read_subseq = rseq_ptr->substr(rstart, overlap_len);
      string ctg_subseq;
      bool ctg_in_cache = false;
      // if max_ctg_seq_cache_size is > 0, then we are using a cache
      if (max_ctg_seq_cache_size) {
        auto it = ctg_seq_cache.find(ctg_loc.cid);
        // found it in the cache, so fetch the subsequence needed
        if (it != ctg_seq_cache.end()) {
          num_ctg_seq_cache_hits++;
          // only get the subsequence needed
          ctg_subseq = it->second.substr(cstart, ctg_aln_len);
        } else {
          // not in cache, fetch the whole sequence
          // space for whatever is getting fetched plus a null terminator
          char *seq_buf = new char[ctg_loc.clen + 1];
          // get the full seq
          rget(ctg_loc.seq_gptr, seq_buf, ctg_loc.clen).wait();
          ctg_seq_bytes_fetched += ctg_loc.clen;
          string ctg_seq(seq_buf, ctg_loc.clen);
          delete[] seq_buf;
          // if space in cache, add the contig
          if (ctg_seq_cache.size() < max_ctg_seq_cache_size) ctg_seq_cache[ctg_loc.cid] = ctg_seq;
          // get the subseq for alignment
          ctg_subseq = ctg_seq.substr(cstart, ctg_aln_len);
        }
      } else {
        // fetch only the substring
        char *seq_buf = new char[ctg_aln_len + 1];
        rget(ctg_loc.seq_gptr + cstart, seq_buf, ctg_aln_len).wait();
        ctg_subseq = string(seq_buf, ctg_aln_len);
        delete[] seq_buf;
        ctg_seq_bytes_fetched += ctg_aln_len;
      }
      align_read(rname, ctg_loc.cid, read_subseq, ctg_subseq, rstart, rlen, cstart, ctg_loc.clen, orient, read_alns, extra_offset);
      num_alns++;
    }

    delete aligned_ctgs_map;
    // only add alns if there are not too many for this read
    if (read_alns->size() < MAX_ALNS_PER_READ) {
      for (int i = 0; i < read_alns->size(); i++) {
        alns->add_aln(read_alns->get_aln(i));
      }
    } else {
      num_excess_alns_reads += 1;
    }
    delete read_alns;
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
    SLOG_VERBOSE("Dropped ", num_dropped, " excessive seed-to-contig mappings (", 
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
      read_fut_chain.wait();
      if (aligned_ctgs_map->size()) num_reads_aligned++;
      kmer_ctg_dht->compute_alns_for_read(aligned_ctgs_map, read_id, read_seq);
      num_reads++;
    }
    progbar.done();
    barrier();
  }
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Parsed ", tot_num_reads, " reads, with ", reduce_one(tot_num_kmers, op_fast_add, 0).wait(), " seeds\n");
  auto tot_num_alns = kmer_ctg_dht->get_num_alns();
  SLOG_VERBOSE("Found ", tot_num_alns, " alignments");
  if (!ALN_EXTRA) SLOG_VERBOSE(" of which ", perc_str(kmer_ctg_dht->get_num_perfect_alns(), tot_num_alns), " are perfect\n");
  else SLOG_VERBOSE("\n");
  auto num_excess_alns_reads = kmer_ctg_dht->get_num_excess_alns_reads();
  if (num_excess_alns_reads)
    SLOG_VERBOSE("Dropped ", num_excess_alns_reads, " reads because of alignments in excess of ", MAX_ALNS_PER_READ, "\n");
  auto num_overlaps = kmer_ctg_dht->get_num_overlaps();
  if (num_overlaps)
    SLOG_VERBOSE("Dropped ", perc_str(num_overlaps, tot_num_alns), " alignments becasue of overlaps\n");
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
  barrier();
  auto num_alns = kmer_ctg_dht.get_num_alns();
  SLOG_VERBOSE("Number of duplicate alignments ", perc_str(alns->get_num_dups(), num_alns), "\n");
  barrier();
}

