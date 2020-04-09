// kcount - kmer counting
// Steven Hofmeyr, LBNL, June 2019

#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "upcxx_utils/flat_aggr_store.hpp"
#include "upcxx_utils/mem_profile.hpp"

#include "utils.hpp"
#include "kmer_dht.hpp"
#include "fastq.hpp"
#include "contigs.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

//#define DBG_DUMP_KMERS

//#define DBG_ADD_KMER DBG
#define DBG_ADD_KMER(...)

extern ofstream _dbgstream;
extern ofstream _logstream;

uint64_t estimate_num_kmers(unsigned kmer_len, vector<FastqReader*> &fqr_list) {
  Timer timer(__FILEFUNC__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t num_kmers = 0;
  int64_t total_records_processed = 0;
  int64_t estimated_total_records = 0;
  for (auto fqr : fqr_list) {
    fqr->reset();
    string id, seq, quals;
    ProgressBar progbar(fqr->my_file_size(), "Scanning reads file to estimate number of kmers");
    size_t tot_bytes_read = 0;
    int64_t records_processed = 0; 
    while (true) {
      size_t bytes_read = fqr->get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      num_lines += 4;
      num_reads++;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      records_processed++;
      // do not read the entire data set for just an estimate
      if (records_processed > 100000) break; 
      if (seq.length() < kmer_len) continue;
      num_kmers += seq.length() - kmer_len + 1;
    }
    total_records_processed += records_processed;
    int64_t bytes_per_record = tot_bytes_read / records_processed;
    estimated_total_records += fqr->my_file_size() / bytes_per_record;
    progbar.done();
    barrier();
  }
  double fraction = (double) total_records_processed / (double) estimated_total_records;
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads), and found ", num_kmers, " kmers\n");
  auto all_num_lines = reduce_one(num_lines / fraction, op_fast_add, 0).wait();
  auto all_num_reads = reduce_one(num_reads / fraction, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_all(num_kmers / fraction, op_fast_add).wait();
  int percent = 100.0 * fraction;
  SLOG_VERBOSE("Processed ", percent, "% of the estimated total of ", (int64_t)all_num_lines,
               " lines (", (int64_t)all_num_reads, " reads), and found an estimated maximum of ", 
               (int64_t)all_num_kmers, " kmers\n");
  return num_kmers / fraction;
}

template<int MAX_K>
static void count_kmers(unsigned kmer_len, int qual_offset, vector<FastqReader*> &fqr_list, 
                        dist_object<KmerDHT<MAX_K>> &kmer_dht, PASS_TYPE pass_type) {
  Timer timer(__FILEFUNC__);
  // probability of an error is P = 10^(-Q/10) where Q is the quality cutoff
  // so we want P = 0.5*1/k (i.e. 50% chance of 1 error)
  // and Q = -10 log10(P)
  // eg qual_cutoff for k=21 is 16, for k=99 is 22.
  //int qual_cutoff = -10 * log10(0.5 / kmer_len); 
  //SLOG_VERBOSE("Using quality cutoff ", qual_cutoff, "\n");
  int qual_cutoff = KCOUNT_QUAL_CUTOFF;
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t num_kmers = 0;
  int64_t num_bad_quals = 0;
  int64_t num_Ns = 0;
  string progbar_prefix = "";
  switch (pass_type) {
    case BLOOM_SET_PASS: progbar_prefix = "Pass 1: Parsing reads file to setup bloom filter"; break;
    case BLOOM_COUNT_PASS: progbar_prefix = "Pass 2: Parsing reads file to count kmers"; break;
    case NO_BLOOM_PASS: progbar_prefix = "Parsing reads file to count kmers"; break;
    default: DIE("Should never get here");
  };
  kmer_dht->set_pass(pass_type);
  barrier();
  for (auto fqr : fqr_list) {
    fqr->reset();
    string id, seq, quals;
    ProgressBar progbar(fqr->my_file_size(), progbar_prefix);
    size_t tot_bytes_read = 0;
    vector<Kmer<MAX_K>> kmers;
    while (true) {
      size_t bytes_read = fqr->get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      num_lines += 4;
      num_reads++;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      if (seq.length() < kmer_len) continue;
      // split into kmers
      Kmer<MAX_K>::get_kmers(kmer_len, seq, kmers);
#ifdef KCOUNT_FILTER_BAD_QUAL_IN_READ
      size_t found_bad_qual_pos = seq.length();
      if (pass_type != BLOOM_SET_PASS) {
        // disable kmer counting of kmers after a bad quality score (of 2) in the read
        // ... but allow extension counting (if an extension q score still passes the QUAL_CUTOFF)
        found_bad_qual_pos = quals.find_first_of(qual_offset + 2);
        if (found_bad_qual_pos == string::npos) found_bad_qual_pos = seq.length(); 
        else num_bad_quals++;
      }
#endif
      // skip kmers that contain an N
      size_t found_N_pos = seq.find_first_of('N');
      if (found_N_pos == string::npos) found_N_pos = seq.length();
      else num_Ns++;
      for (int i = 1; i < kmers.size() - 1; i++) {
        // skip kmers that contain an N
        if (i + kmer_len > found_N_pos) {
          i = found_N_pos; // skip
          // find the next N
          found_N_pos = seq.find_first_of('N', found_N_pos + 1);
          if (found_N_pos == string::npos) found_N_pos = seq.length();
          else num_Ns++;
          continue;
        }
        int count = 1;
#ifdef KCOUNT_FILTER_BAD_QUAL_IN_READ
        if (pass_type != BLOOM_SET_PASS && i + kmer_len > found_bad_qual_pos) count = 0;
#endif
        char left_base = seq[i - 1];
        if (quals[i - 1] < qual_offset + qual_cutoff) left_base = '0';
        char right_base = seq[i + kmer_len];
        if (quals[i + kmer_len] < qual_offset + qual_cutoff) right_base = '0';
        kmer_dht->add_kmer(kmers[i], left_base, right_base, count);
        DBG_ADD_KMER("kcount add_kmer ", kmers[i].to_string(), " count ", count, "\n");
        num_kmers++;
      }
      progress();
    }
    progbar.done();
    kmer_dht->flush_updates();
  }
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads)\n");
  auto all_num_lines = reduce_one(num_lines, op_fast_add, 0).wait();
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto all_num_bad_quals = reduce_one(num_bad_quals, op_fast_add, 0).wait();
  auto all_num_Ns = reduce_one(num_Ns, op_fast_add, 0).wait();
  auto all_distinct_kmers = kmer_dht->get_num_kmers();
  SLOG_VERBOSE("Processed a total of ", all_num_lines, " lines (", all_num_reads, " reads)\n");
  if (pass_type != BLOOM_SET_PASS) {
    SLOG_VERBOSE("Found ", perc_str(all_distinct_kmers, all_num_kmers), " unique kmers\n");
    if (all_num_bad_quals) SLOG_VERBOSE("Found ", perc_str(all_num_bad_quals, all_num_kmers), " bad quality positions\n");
    if (all_num_Ns) SLOG_VERBOSE("Found ", perc_str(all_num_Ns, all_num_kmers), " kmers with Ns\n");
  }
}

// count ctg kmers if using bloom
template<int MAX_K>
static void count_ctg_kmers(unsigned kmer_len, Contigs &ctgs, dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  Timer timer(__FILEFUNC__);
  ProgressBar progbar(ctgs.size(), "Counting kmers in contigs");
  int64_t num_kmers = 0;
  vector<Kmer<MAX_K>> kmers;
  kmer_dht->set_pass(CTG_BLOOM_SET_PASS);
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() >= kmer_len) {
      Kmer<MAX_K>::get_kmers(kmer_len, ctg->seq, kmers);
      if (kmers.size() != ctg->seq.length() - kmer_len + 1)
        DIE("kmers size mismatch ", kmers.size(), " != ", (ctg->seq.length() - kmer_len + 1), " '", ctg->seq, "'");
      for (int i = 1; i < ctg->seq.length() - kmer_len; i++) {
        kmer_dht->add_kmer(kmers[i], ctg->seq[i - 1], ctg->seq[i + kmer_len], 1);
      }
      num_kmers += kmers.size();
    }
    progress();
  }
  progbar.done();
  kmer_dht->flush_updates();
  DBG("This rank processed ", ctgs.size(), " contigs and ", num_kmers , " kmers\n");
  auto all_num_ctgs = reduce_one(ctgs.size(), op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed a total of ", all_num_ctgs, " contigs and ", all_num_kmers, " kmers\n");
  barrier();
}

template<int MAX_K>
static void add_ctg_kmers(unsigned kmer_len, unsigned prev_kmer_len, Contigs &ctgs, dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  Timer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  int64_t num_prev_kmers = kmer_dht->get_num_kmers();
#ifdef USE_KMER_DEPTH
  double tot_depth_diff = 0;
  double max_depth_diff = 0;
#endif
  ProgressBar progbar(ctgs.size(), "Adding extra contig kmers from kmer length " + to_string(prev_kmer_len));
  vector<Kmer<MAX_K>> kmers;
  kmer_dht->set_pass(CTG_KMERS_PASS);
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() >= kmer_len + 2) {
      Kmer<MAX_K>::get_kmers(kmer_len, ctg->seq, kmers);
      if (kmers.size() != ctg->seq.length() - kmer_len + 1)
        DIE("kmers size mismatch ", kmers.size(), " != ", (ctg->seq.length() - kmer_len + 1), " '", ctg->seq, "'");
      for (int i = 1; i < ctg->seq.length() - kmer_len; i++) {
        uint16_t depth = ctg->depth;
#ifdef USE_KMER_DEPTHS
        uint16_t kmer_depth = ctg->get_kmer_depth(i, kmer_len, prev_kmer_len);
        tot_depth_diff += (double)(kmer_depth - depth) / depth;
        max_depth_diff = max(max_depth_diff, abs(kmer_depth - depth));
        depth = kmer_depth;
#endif
        kmer_dht->add_kmer(kmers[i], ctg->seq[i - 1], ctg->seq[i + kmer_len], depth);
        num_kmers++;
      }
    }
    progress();
  }
  progbar.done();
  kmer_dht->flush_updates();
  DBG("This rank processed ", ctgs.size(), " contigs and ", num_kmers , " kmers\n");
  auto all_num_ctgs = reduce_one(ctgs.size(), op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed a total of ", all_num_ctgs, " contigs and ", all_num_kmers, " kmers\n");
  SLOG_VERBOSE("Found ", perc_str(kmer_dht->get_num_kmers() - num_prev_kmers, all_num_kmers), " additional unique kmers\n");
#ifdef USE_KMER_DEPTH
  auto all_tot_depth_diff = reduce_one(tot_depth_diff, op_fast_add, 0).wait();
  SLOG_VERBOSE(KLRED, "Average depth diff ", all_tot_depth_diff / all_num_kmers, " max depth diff ", 
               reduce_one(max_depth_diff, op_fast_max, 0).wait(), KNORM, "\n");
#endif
}

template<int MAX_K>
void analyze_kmers(unsigned kmer_len, unsigned prev_kmer_len, int qual_offset, vector<FastqReader*> &fqr_list,
                   double dynamic_min_depth, int dmin_thres, Contigs &ctgs, dist_object<KmerDHT<MAX_K>> &kmer_dht,
                   double &num_kmers_factor) {
  Timer timer(__FILEFUNC__);
  
  _dynamic_min_depth = dynamic_min_depth;
  _dmin_thres = dmin_thres;
    
  if (kmer_dht->get_use_bloom()) {
    count_kmers(kmer_len, qual_offset, fqr_list, kmer_dht, BLOOM_SET_PASS);
    if (ctgs.size()) count_ctg_kmers(kmer_len, ctgs, kmer_dht);
    kmer_dht->reserve_space_and_clear_bloom1();
    count_kmers(kmer_len, qual_offset, fqr_list, kmer_dht, BLOOM_COUNT_PASS);
  } else {
    count_kmers(kmer_len, qual_offset, fqr_list, kmer_dht, NO_BLOOM_PASS);
  }
  barrier();
  kmer_dht->print_load_factor();
  barrier();
  num_kmers_factor = kmer_dht->get_num_kmers_factor();
  kmer_dht->purge_kmers(2);
  int64_t new_count = kmer_dht->get_num_kmers();
  SLOG_VERBOSE("After purge of kmers < 2, there are ", new_count, " unique kmers\n");
  barrier();
  if (ctgs.size()) {
    add_ctg_kmers(kmer_len, prev_kmer_len, ctgs, kmer_dht);
    kmer_dht->purge_kmers(1);
  }
  barrier();
  kmer_dht->compute_kmer_exts();
#ifdef DEBUG
  // FIXME: dump if an option specifies
#ifdef DBG_DUMP_KMERS
  kmer_dht->dump_kmers(kmer_len);
#endif
#endif
  barrier();
  //kmer_dht->purge_fx_kmers();
}

#define AK(KMER_LEN) \
  template \
  void analyze_kmers<KMER_LEN>(unsigned, unsigned, int, vector<FastqReader*>&, double, int, Contigs&, dist_object< KmerDHT<KMER_LEN> >&, double&)

AK(32);
#if MAX_BUILD_KMER >= 64
AK(64);
#endif
#if MAX_BUILD_KMER >= 96
AK(96);
#endif
#if MAX_BUILD_KMER >= 128
AK(128);
#endif
#if MAX_BUILD_KMER >= 160
AK(160);
#endif

#undef AK

