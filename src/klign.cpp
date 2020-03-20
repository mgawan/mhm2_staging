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
#include "alignments.hpp"

using namespace std;
using namespace upcxx;

//#define DUMP_ALNS

using cid_t = int64_t;

struct CtgLoc {
  cid_t cid;
  global_ptr<char> seq_gptr;
  int clen;
  int pos_in_ctg;
  bool is_rc;
};


struct ReadAndCtgLoc {
  int pos_in_read;
  bool read_is_rc;
  CtgLoc ctg_loc;
};

template<int MAX_K>
struct KmerCtgLoc {
  Kmer<MAX_K> kmer;
  CtgLoc ctg_loc;
  UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
};
  

// global variables to avoid passing dist objs to rpcs
static int64_t _num_dropped_seed_to_ctgs = 0;

template<int MAX_K>
class KmerCtgDHT {

  struct KmerAndCtgLoc {
    Kmer<MAX_K> kmer;
    CtgLoc ctg_loc;
    UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
  };

  AggrStore<KmerAndCtgLoc> kmer_store;
  // maps a kmer to a contig - the first part of the pair is set to true if this is a conflict,
  // with a kmer mapping to multiple contigs
  using kmer_map_t = HASH_TABLE<Kmer<MAX_K>, pair<bool, CtgLoc>>;
  dist_object<kmer_map_t> kmer_map;
#ifdef DUMP_ALNS
  zstr::ofstream *alns_file;
#endif  
  int64_t num_alns;
  int64_t num_perfect_alns;
  int64_t num_overlaps;
  
  int64_t max_ctg_seq_cache_size;
  HASH_TABLE<cid_t, string> ctg_seq_cache;
  int64_t num_ctg_seq_cache_hits;
  int64_t ctg_seq_bytes_fetched;
  
  Alns *alns;
  
  // default aligner and filter
  StripedSmithWaterman::Aligner ssw_aligner;
  StripedSmithWaterman::Filter ssw_filter;
  
  struct InsertKmer {
    void operator()(KmerAndCtgLoc &kmer_and_ctg_loc, dist_object<kmer_map_t> &kmer_map) {
      CtgLoc ctg_loc = kmer_and_ctg_loc.ctg_loc;
      const auto it = kmer_map->find(kmer_and_ctg_loc.kmer);
      if (it == kmer_map->end()) {
        kmer_map->insert({kmer_and_ctg_loc.kmer, {false, ctg_loc}});
      } else {
        // in this case, we have a conflict, i.e. the kmer maps to multiple contigs
        it->second.first = true;
        _num_dropped_seed_to_ctgs++;
      }
    }
  };
  dist_object<InsertKmer> insert_kmer;

  void align_read(const string &rname, int64_t cid, const string &rseq, const string &cseq, int rstart, int rlen,
                  int cstart, int clen, char orient, Alns &read_alns, int ctg_extra_offset, int read_extra_offset,
                  int overlap_len) {
    Aln aln;
    StripedSmithWaterman::Alignment ssw_aln;
    // use hamming distance for checking match first - works for > 90% of matches and reduces time spent doing SSW
    //int hdist = hamming_dist(rseq.substr(read_extra_offset, overlap_len), cseq.substr(ctg_extra_offset, overlap_len) , true);
    // allow for ~1% read error and at least 1
    // int max_hdist = max(overlap_len / 100, 1);
    //if (hdist < max_hdist) {
    if (cseq.compare(ctg_extra_offset, overlap_len, rseq, read_extra_offset, overlap_len) == 0) {
      num_perfect_alns++;
      ssw_aln.ref_begin = read_extra_offset;
      ssw_aln.ref_end = read_extra_offset + overlap_len - 1;
      ssw_aln.query_begin = ctg_extra_offset;
      ssw_aln.query_end = ctg_extra_offset + overlap_len - 1;
      // the mismatch penalty is 3 
      ssw_aln.sw_score = overlap_len;
      //ssw_aln.sw_score = overlap_len - 3 * hdist;
      ssw_aln.sw_score_next_best = 0;
    } else {
      // make sure upcxx progress is done before starting alignment
      discharge();
      // contig is the ref, read is the query - done this way so that we can potentially do multiple alns to each read
      // this is also the way it's done in meraligner
      ssw_aligner.Align(cseq.c_str(), rseq.c_str(), rseq.length(), ssw_filter, &ssw_aln, max((int)(rseq.length() / 2), 15));
    }

    int rstop = rstart + ssw_aln.ref_end + 1;
    rstart = rstart + ssw_aln.ref_begin;
    int cstop = cstart + ssw_aln.query_end + 1;
    cstart = cstart + ssw_aln.query_begin;

    if (orient == '-') switch_orient(rstart, rstop, rlen);

    // for some reason, on Cori icc this causes an internal compiler error:
    // internal error: assertion failed at: "shared/cfe/edgcpfe/overload.c", line 9538
    aln = { .read_id = rname, .cid = cid,
            .rstart = rstart, .rstop = rstop, .rlen = rlen,
            .cstart = cstart, .cstop = cstop, .clen = clen,
            .orient = orient, .score1 = ssw_aln.sw_score,
            .score2 = ssw_aln.sw_score_next_best };
    /*
    aln.read_id = rname; aln.cid = cid;
    aln.rstart = rstart; aln.rstop = rstop; aln.rlen = rlen;
    aln.cstart = cstart; aln.cstop = cstop; aln.clen = clen;
    aln.orient = orient; aln.score1 = ssw_aln.sw_score; aln.score2 = ssw_aln.sw_score_next_best;
    */
#ifdef DUMP_ALNS
    *alns_file << "MERALIGNER\t" << aln.to_string() << endl;
#endif
    read_alns.add_aln(aln);
  }
  
public:

  int kmer_len;  
  
  // aligner construction: SSW internal defaults are 2 2 3 1
  KmerCtgDHT(int kmer_len, int max_store_size, int max_ctg_cache, Alns &alns)
    : kmer_map({})
    , kmer_store({})
    , insert_kmer({})
    , ssw_aligner(SSW_MATCH_SCORE, SSW_MISMATCH_COST, SSW_GAP_OPENING_COST, SSW_GAP_EXTENDING_COST, SSW_AMBIGUITY_COST)
    , num_alns(0)
    , num_perfect_alns(0)
    , num_overlaps(0)
    , kmer_len(kmer_len)
    , num_ctg_seq_cache_hits(0)
    , ctg_seq_bytes_fetched(0)
    , max_ctg_seq_cache_size(max_ctg_cache)
    , alns(&alns) {
    
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

  size_t get_target_rank(Kmer<MAX_K> &kmer) {
    return std::hash<Kmer<MAX_K>>{}(kmer) % rank_n();
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

  int64_t get_num_overlaps(bool all = false) {
    if (!all) return reduce_one(num_overlaps, op_fast_add, 0).wait();
    else return reduce_all(num_overlaps, op_fast_add).wait();
  }

  int64_t get_num_dropped_seed_to_ctgs(bool all = false) {
    if (!all) return reduce_one(_num_dropped_seed_to_ctgs, op_fast_add, 0).wait();
    else return reduce_all(_num_dropped_seed_to_ctgs, op_fast_add).wait();
  }

  double get_av_ssw_secs() {
    return reduce_one(ssw_aligner.get_ssw_secs(), op_fast_add, 0).wait() / rank_n();
  }
  
  double get_max_ssw_secs() {
    return reduce_one(ssw_aligner.get_ssw_secs(), op_fast_max, 0).wait();
  }

  int64_t get_ctg_seq_cache_hits() {
    return reduce_one(num_ctg_seq_cache_hits, op_fast_add, 0).wait();
  }

  int64_t get_ctg_seq_bytes_fetched() {
    return reduce_one(ctg_seq_bytes_fetched, op_fast_add, 0).wait();
  }    

  void add_kmer(Kmer<MAX_K> kmer, CtgLoc &ctg_loc) {
    Kmer<MAX_K> kmer_rc = kmer.revcomp();
    ctg_loc.is_rc = false;
    if (kmer_rc < kmer) {
      kmer = kmer_rc;
      ctg_loc.is_rc = true;
    }
    KmerAndCtgLoc kmer_and_ctg_loc = { kmer, ctg_loc };
    kmer_store.update(get_target_rank(kmer), kmer_and_ctg_loc, insert_kmer, kmer_map);
  }

  void flush_add_kmers() {
    Timer timer(__FILEFUNC__);
    barrier();
    kmer_store.flush_updates(insert_kmer, kmer_map);
    barrier();
  }

  future<vector<KmerCtgLoc<MAX_K>>> get_ctgs_with_kmers(int target_rank, vector<Kmer<MAX_K>> &kmers) {
    return rpc(target_rank,
               [](vector<Kmer<MAX_K>> kmers, dist_object<kmer_map_t> &kmer_map) {
                 vector<KmerCtgLoc<MAX_K>> kmer_ctg_locs;
                 for (auto &kmer : kmers) {
                   const auto it = kmer_map->find(kmer);
                   if (it == kmer_map->end()) continue;
                   // skip conflicts
                   if (it->second.first) continue;
                   // now add it
                   kmer_ctg_locs.push_back({kmer, it->second.second});
                 }
                 return kmer_ctg_locs;
               }, kmers, kmer_map);
  }

#ifdef DEBUG
  void dump_ctg_kmers() {
    Timer timer(__FILEFUNC__);
    string dump_fname = "ctg_kmers-" + to_string(kmer_len) + ".txt.gz";
    get_rank_path(dump_fname, rank_me());
    zstr::ofstream dump_file(dump_fname);
    ostringstream out_buf;
    ProgressBar progbar(kmer_map->size(), "Dumping kmers to " + dump_fname);
    int64_t i = 0;
    for (auto &elem : *kmer_map) {
      auto ctg_loc = &elem.second;
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
#endif
  
  void compute_alns_for_read(HASH_TABLE<cid_t, ReadAndCtgLoc> *aligned_ctgs_map, const string &rname, string rseq) {
    int rlen = rseq.length();
    string rseq_rc = revcomp(rseq);
    Alns read_alns;
    for (auto &elem : *aligned_ctgs_map) {
      progress();
      int pos_in_read = elem.second.pos_in_read;
      bool read_kmer_is_rc = elem.second.read_is_rc;
      CtgLoc ctg_loc = elem.second.ctg_loc;
      char orient = '+';
      string *rseq_ptr = &rseq;
      if (ctg_loc.is_rc != read_kmer_is_rc) {
        // it's revcomp in either contig or read, but not in both or neither
        orient = '-';
        pos_in_read = rlen - (kmer_len + pos_in_read);
        rseq_ptr = &rseq_rc;
      }
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
      int ctg_aln_len = overlap_len + min(ctg_loc.clen - (cstart + overlap_len), KLIGN_EXPAND_BASES);
      int ctg_extra_offset = min(cstart, KLIGN_EXPAND_BASES);
      cstart -= ctg_extra_offset;
      ctg_aln_len += ctg_extra_offset;

      // use the whole read, to account for possible indels
      int read_aln_len = rlen;
      int read_extra_offset = rstart;
      rstart = 0;
      
      string read_subseq = rseq_ptr->substr(rstart, read_aln_len);
      string ctg_subseq;
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
      align_read(rname, ctg_loc.cid, read_subseq, ctg_subseq, rstart, rlen, cstart, ctg_loc.clen, orient, read_alns,
                 ctg_extra_offset, read_extra_offset, overlap_len);
      num_alns++;
    }
    // sort the alns from best score to worst - this could be used in spanner later
    sort(read_alns.begin(), read_alns.end(), 
         [](const auto &elem1, const auto &elem2) {
           return elem1.score1 > elem2.score1;
         });
    for (int i = 0; i < read_alns.size(); i++) {
      alns->add_aln(read_alns.get_aln(i));
    }
  }
};


template<int MAX_K>
static void build_alignment_index(KmerCtgDHT<MAX_K> &kmer_ctg_dht, Contigs &ctgs) {
  Timer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
  vector<Kmer<MAX_K>> kmers;
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    global_ptr<char> seq_gptr = allocate<char>(ctg->seq.length() + 1);
    strcpy(seq_gptr.local(), ctg->seq.c_str());
    CtgLoc ctg_loc = { .cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length() };
    Kmer<MAX_K>::get_kmers(kmer_ctg_dht.kmer_len, ctg->seq, kmers);
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
  auto num_dropped_seed_to_ctgs = kmer_ctg_dht.get_num_dropped_seed_to_ctgs(); 
  if (num_dropped_seed_to_ctgs) 
    SLOG_VERBOSE("Dropped ", num_dropped_seed_to_ctgs, " excessive seed-to-contig mappings (", 
                 setprecision(2), fixed, (100.0 * num_dropped_seed_to_ctgs / tot_num_kmers), "%)\n");
}


struct ReadRecord {
  string id;
  string seq;
  string quals;
  
  HASH_TABLE<cid_t, ReadAndCtgLoc> aligned_ctgs_map;
  
  ReadRecord(const string &id, const string &seq, const string &quals) : id(id), seq(seq), quals(quals) {}
};


struct KmerToRead {
  ReadRecord *read_record;
  int pos_in_read;
  bool is_rc;
};


template<int MAX_K>
static int align_kmers(KmerCtgDHT<MAX_K> &kmer_ctg_dht, HASH_TABLE<Kmer<MAX_K>, vector<KmerToRead>> &kmer_read_map,
                       vector<ReadRecord*> &read_records, IntermittentTimer &compute_alns_timer,
                       IntermittentTimer &get_ctgs_timer, int64_t &num_excess_alns_reads) {
  // extract a list of kmers for each target rank
  auto kmer_lists = new vector<Kmer<MAX_K>>[rank_n()];
  for (auto &elem : kmer_read_map) {
    auto kmer = elem.first;
    kmer_lists[kmer_ctg_dht.get_target_rank(kmer)].push_back(kmer);
  }
  //size_t min_kmers = 10000000, max_kmers = 0, num_kmers = 0;
  get_ctgs_timer.start();
  future<> fut_chain = make_future();
  // fetch ctgs for each set of kmers from target ranks
  for (int i = 0; i < rank_n(); i++) {
    progress();
    //min_kmers = min(min_kmers, kmer_lists[i].size());
    //max_kmers = max(max_kmers, kmer_lists[i].size());
    //num_kmers += kmer_lists[i].size();
    auto fut = kmer_ctg_dht.get_ctgs_with_kmers(i, kmer_lists[i]).then(
      [&](vector<KmerCtgLoc<MAX_K>> kmer_ctg_locs) {
        // iterate through the kmers, each one has an associated ctg location
        for (auto &kmer_ctg_loc : kmer_ctg_locs) {
          // get the reads that this kmer mapped to
          auto kmer_read_map_it = kmer_read_map.find(kmer_ctg_loc.kmer);
          if (kmer_read_map_it == kmer_read_map.end()) DIE("Could not find kmer ", kmer_ctg_loc.kmer);
          // this is a list of the reads
          auto &kmer_to_reads = kmer_read_map_it->second;
          // now add the ctg loc to all the reads
          for (auto &kmer_to_read : kmer_to_reads) {
            auto read_record = kmer_to_read.read_record;
            int pos_in_read = kmer_to_read.pos_in_read;
            bool read_is_rc = kmer_to_read.is_rc;
            if (read_record->aligned_ctgs_map.size() >= KLIGN_MAX_ALNS_PER_READ) {
              // too many mappings for this read, stop adding to it
              num_excess_alns_reads++;
              continue;
            }
            // this here ensures that we don't insert duplicate mappings
            read_record->aligned_ctgs_map.insert({kmer_ctg_loc.ctg_loc.cid, {pos_in_read, read_is_rc, kmer_ctg_loc.ctg_loc}});
          }
        }
      });
    fut_chain = when_all(fut_chain, fut);
  }
  fut_chain.wait();
  get_ctgs_timer.stop();
  delete[] kmer_lists;
  kmer_read_map.clear();
  compute_alns_timer.start();
  int num_reads_aligned = 0;
  // create a new list of records with all reads having < KLIGN_MAX_ALNS_PER_READ, i.e. those without excessive mappings
  vector<ReadRecord*> good_read_records;
  good_read_records.reserve(read_records.size());
  for (auto read_record : read_records) {
    if (read_record->aligned_ctgs_map.size() < KLIGN_MAX_ALNS_PER_READ)
      good_read_records.push_back(read_record);
  }
  // fetch contigs for each read
  for (auto read_record : good_read_records) {
    progress();
    // compute alignments
    if (read_record->aligned_ctgs_map.size()) {
      num_reads_aligned++;
      // when all the ctgs are fetched, do the alignments
      //num_reads_aligned++;
      //compute_alns_timer.start();
      kmer_ctg_dht.compute_alns_for_read(&read_record->aligned_ctgs_map, read_record->id, read_record->seq);
      //compute_alns_timer.stop();
    }
    delete read_record;
  }
  read_records.clear();
  compute_alns_timer.stop();
  return num_reads_aligned;
}

template<int MAX_K>
static void do_alignments(KmerCtgDHT<MAX_K> &kmer_ctg_dht, vector<FastqReader*> &fqr_list) {
  Timer timer(__FILEFUNC__);
  int64_t tot_num_kmers = 0;
  int64_t num_reads = 0;
  int64_t num_reads_aligned = 0, num_excess_alns_reads = 0;
  IntermittentTimer compute_alns_timer(__FILENAME__ + string(":") + "Compute alns");
  IntermittentTimer get_reads_timer(__FILENAME__ + string(":") + "Get reads");
  IntermittentTimer get_ctgs_timer(__FILENAME__ + string(":") + "Get ctgs with kmer");
  barrier();
  for (auto fqr : fqr_list) {
    fqr->reset();
    string read_id, read_seq, quals;
    ProgressBar progbar(fqr->my_file_size(), "Aligning reads to contigs");
    size_t tot_bytes_read = 0;
    vector<ReadRecord*> read_records;
    HASH_TABLE<Kmer<MAX_K>, vector<KmerToRead>> kmer_read_map;
    vector<Kmer<MAX_K>> kmers;
    while (true) {
      progress();
      get_reads_timer.start();
      size_t bytes_read = fqr->get_next_fq_record(read_id, read_seq, quals);
      get_reads_timer.stop();
      if (!bytes_read) break;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      // this happens when a placeholder read with just a single N character is added after merging reads
      if (kmer_ctg_dht.kmer_len > read_seq.length()) continue;
      Kmer<MAX_K>::get_kmers(kmer_ctg_dht.kmer_len, read_seq, kmers);
      tot_num_kmers += kmers.size();
      ReadRecord *read_record = new ReadRecord(read_id, read_seq, quals);
      read_records.push_back(read_record);
      bool filled = false;
      for (int i = 0; i < kmers.size(); i += KLIGN_SEED_SPACE) {
        Kmer<MAX_K> kmer = kmers[i];
        Kmer<MAX_K> kmer_rc = kmer.revcomp();
        bool is_rc = false;
        if (kmer_rc < kmer) {
          kmer = kmer_rc;
          is_rc = true;
        }
        auto it = kmer_read_map.find(kmer);
        if (it == kmer_read_map.end()) it = kmer_read_map.insert({kmer, {}}).first;
        it->second.push_back({read_record, i, is_rc});
        if (kmer_read_map.size() >= KLIGN_CTG_FETCH_BUF_SIZE) filled = true;
      }
      if (filled) 
        num_reads_aligned += align_kmers(kmer_ctg_dht, kmer_read_map, read_records, compute_alns_timer, get_ctgs_timer,
                                         num_excess_alns_reads);
      num_reads++;
    }
    if (read_records.size())
      num_reads_aligned += align_kmers(kmer_ctg_dht, kmer_read_map, read_records, compute_alns_timer, get_ctgs_timer,
                                       num_excess_alns_reads);
    progbar.done();
    barrier();
  }
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Parsed ", tot_num_reads, " reads, with ", reduce_one(tot_num_kmers, op_fast_add, 0).wait(), " seeds\n");
  auto tot_num_alns = kmer_ctg_dht.get_num_alns();
  SLOG_VERBOSE("Found ", tot_num_alns, " alignments of which ", perc_str(kmer_ctg_dht.get_num_perfect_alns(), tot_num_alns), 
               " are perfect\n");
  auto tot_excess_alns_reads = reduce_one(num_excess_alns_reads, op_fast_add, 0).wait();
  if (num_excess_alns_reads)
    SLOG_VERBOSE("Dropped ", tot_excess_alns_reads, " reads because of alignments in excess of ", KLIGN_MAX_ALNS_PER_READ, "\n");
  auto num_overlaps = kmer_ctg_dht.get_num_overlaps();
  if (num_overlaps)
    SLOG_VERBOSE("Dropped ", perc_str(num_overlaps, tot_num_alns), " alignments becasue of overlaps\n");
  auto tot_num_reads_aligned = reduce_one(num_reads_aligned, op_fast_add, 0).wait();
  SLOG_VERBOSE("Mapped ", perc_str(tot_num_reads_aligned, tot_num_reads), " reads to contigs\n");
  SLOG_VERBOSE("Average mappings per read ", (double)tot_num_alns / tot_num_reads_aligned, "\n");

  SLOG_VERBOSE("Ctg cache hits ", perc_str(kmer_ctg_dht.get_ctg_seq_cache_hits(), tot_num_alns), "\n");
  SLOG_VERBOSE("Fetched ", get_size_str(kmer_ctg_dht.get_ctg_seq_bytes_fetched()), " of contig sequences\n");
  
  double av_ssw_secs = kmer_ctg_dht.get_av_ssw_secs();
  double max_ssw_secs = kmer_ctg_dht.get_max_ssw_secs();
  SLOG_VERBOSE("Average SSW time ", setprecision(2), fixed, av_ssw_secs, " s, ",
               "max ", setprecision(2), fixed, max_ssw_secs, " s, ",
               "balance ", setprecision(2), fixed, av_ssw_secs / max_ssw_secs, "\n");
  compute_alns_timer.done_barrier();
  get_reads_timer.done_barrier();
  get_ctgs_timer.done_barrier();
}

template<int MAX_K>
void find_alignments(unsigned kmer_len, vector<FastqReader*> &fqr_list, int max_store_size, int max_ctg_cache,
                     Contigs &ctgs, Alns &alns) {
  Timer timer(__FILEFUNC__);
  _num_dropped_seed_to_ctgs = 0;
  Kmer<MAX_K>::set_k(kmer_len);
  KmerCtgDHT<MAX_K> kmer_ctg_dht(kmer_len, max_store_size, max_ctg_cache, alns);
  barrier();
  build_alignment_index(kmer_ctg_dht, ctgs);
#ifdef DEBUG
  //kmer_ctg_dht.dump_ctg_kmers();
#endif
  do_alignments(kmer_ctg_dht, fqr_list);
  barrier();
  auto num_alns = kmer_ctg_dht.get_num_alns();
  SLOG_VERBOSE("Number of duplicate alignments ", perc_str(alns.get_num_dups(), num_alns), "\n");
  barrier();
}

template 
void find_alignments<32>(unsigned kmer_len, vector<FastqReader*> &fqr_list, int max_store_size, int max_ctg_cache, 
                         Contigs &ctgs, Alns &alns);
template 
void find_alignments<64>(unsigned kmer_len, vector<FastqReader*> &fqr_list, int max_store_size, int max_ctg_cache, 
                         Contigs &ctgs, Alns &alns);
template 
void find_alignments<96>(unsigned kmer_len, vector<FastqReader*> &fqr_list, int max_store_size, int max_ctg_cache, 
                         Contigs &ctgs, Alns &alns);
template 
void find_alignments<128>(unsigned kmer_len, vector<FastqReader*> &fqr_list, int max_store_size, int max_ctg_cache, 
                          Contigs &ctgs, Alns &alns);
template 
void find_alignments<160>(unsigned kmer_len, vector<FastqReader*> &fqr_list, int max_store_size, int max_ctg_cache, 
                          Contigs &ctgs, Alns &alns);
