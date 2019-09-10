// kcount - kmer counting
// Steven Hofmeyr, LBNL, June 2019

#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "kcount.hpp"
#include "progressbar.hpp"
#include "kmer_dht.hpp"
#include "fastq.hpp"


using namespace std;
using namespace upcxx;


extern ofstream _dbgstream;

uint64_t estimate_cardinality(unsigned kmer_len, vector<string> reads_fname_list)
{
  Timer timer(__func__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t num_kmers = 0;
  int64_t estimated_total_records = 0;
  int64_t total_records_processed = 0;
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    FastqReader fqr(merged_reads_fname, PER_RANK_FILE);
    string id, seq, quals;
    ProgressBar progbar(fqr.my_file_size(), "Scanning reads file to estimate cardinality");
    size_t tot_bytes_read = 0;
    int64_t records_processed = 0; 
    while (true) {
      size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
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
    estimated_total_records += records_processed * fqr.my_file_size() / fqr.tell();
    progbar.done();
    barrier();
  }
  //auto all_records_processed = reduce_one(total_records_processed, op_fast_add, 0).wait();
  //SOUT("all records processed ", all_records_processed, "\n");
  double fraction = (double) total_records_processed / (double) estimated_total_records;
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads), and found ", num_kmers, " kmers\n");
  auto all_num_lines = reduce_one(num_lines / fraction, op_fast_add, 0).wait();
  auto all_num_reads = reduce_one(num_reads / fraction, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_all(num_kmers / fraction, op_fast_add).wait();
  int percent = 100.0 * fraction;
  SOUT("Processed ", percent, " % of the estimated total of ", all_num_lines,
       " lines (", all_num_reads, " reads), and found a maximum of ", all_num_kmers, " kmers\n");
  int my_cardinality = all_num_kmers / rank_n();
  SOUT("Cardinality estimated as ", my_cardinality, "\n");
  return my_cardinality;
}

void count_kmers(unsigned kmer_len, int qual_offset, vector<string> reads_fname_list, dist_object<KmerDHT> &kmer_dht, PASS_TYPE pass_type)
{
  Timer timer(__func__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t num_kmers = 0;
  string progbar_prefix = "";
  switch (pass_type) {
    case BLOOM_SET_PASS: progbar_prefix = "Pass 1: Parsing reads file to setup bloom filter"; break;
    case BLOOM_COUNT_PASS: progbar_prefix = "Pass 2: Parsing reads file to count kmers"; break;
    case NO_BLOOM_PASS: progbar_prefix = "Parsing reads file to count kmers"; break;
  };
  char special = qual_offset + 2;
  IntermittentTimer t_io("reads IO");
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    FastqReader fqr(merged_reads_fname, PER_RANK_FILE);
    string id, seq, quals;
    ProgressBar progbar(fqr.my_file_size(), progbar_prefix);
    size_t tot_bytes_read = 0;
    while (true) {
      size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      num_lines += 4;
      num_reads++;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      if (seq.length() < kmer_len) continue;
      
      // split into kmers
      auto kmers = Kmer::get_kmers(seq);
      // disable kmer counting of kmers after a bad quality score (of 2) in the read
      // ... but allow extension counting (if an extention q score still passes the QUAL_CUTOFF)
      size_t found_bad_qual = quals.find_first_of(special);
      if (found_bad_qual == string::npos) found_bad_qual = seq.length();   // remember that the last valid position is length()-1
      int found_bad_qual_kmer = found_bad_qual - kmer_len + 1;
      assert( (int) kmers.size() >= found_bad_qual_kmer );
      // skip kmers that contain an N
      size_t found_N = seq.find_first_of('N');
      if (found_N == string::npos) found_N = seq.length();
      for (int i = 0; i < kmers.size(); i++) {
        // skip kmers that contain an N
        if (i + kmer_len > found_N) {
          i = found_N; // skip
          // find the next N
          found_N = seq.find_first_of('N', found_N+1);
          if (found_N == string::npos) found_N = seq.length();
          continue;
        }
        char left_base = '0';
        if (i > 0 && quals[i - 1] >= qual_offset + QUAL_CUTOFF) {
          left_base = seq[i - 1];
        }
        char right_base = '0';
        if (i + kmer_len < seq.length() && quals[i + kmer_len] >= qual_offset + QUAL_CUTOFF) {
          right_base = seq[i + kmer_len];
        }
        int count = (i < found_bad_qual_kmer) ? 1 : 0;
        kmer_dht->add_kmer(kmers[i], left_base, right_base, count, pass_type);
        DBG("kcount add_kmer ", kmers[i].to_string(), " count ", count, "\n");
        num_kmers++;
      }
      progress();
    }
    progbar.done();
    kmer_dht->flush_updates(pass_type);
  }
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads), max read length ", max_read_len, "\n");
  auto all_num_lines = reduce_one(num_lines, op_fast_add, 0).wait();
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto all_distinct_kmers = kmer_dht->get_num_kmers();
  SOUT("Processed a total of ", all_num_lines, " lines (", all_num_reads, " reads)\n");
  if (pass_type != BLOOM_SET_PASS) SOUT("Found ", perc_str(all_distinct_kmers, all_num_kmers), " unique kmers\n");
}


// count ctg kmers if using bloom
void count_ctg_kmers(unsigned kmer_len, Contigs &ctgs, dist_object<KmerDHT> &kmer_dht)
{
  Timer timer(__func__);
  ProgressBar progbar(ctgs.size(), "Counting kmers in contigs");
  int64_t num_kmers = 0;
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() >= kmer_len) {
      auto kmers = Kmer::get_kmers(ctg->seq);
      if (kmers.size() != ctg->seq.length() - kmer_len + 1)
        DIE("kmers size mismatch ", kmers.size(), " != ", (ctg->seq.length() - kmer_len + 1), " '", ctg->seq, "'");
      for (int i = 1; i < ctg->seq.length() - kmer_len; i++) {
        kmer_dht->add_kmer(kmers[i], ctg->seq[i - 1], ctg->seq[i + kmer_len], 1, CTG_BLOOM_SET_PASS);
      }
      num_kmers += kmers.size();
    }
    progress();
  }
  kmer_dht->flush_updates(CTG_BLOOM_SET_PASS);
  progbar.done();
  DBG("This rank processed ", ctgs.size(), " contigs and ", num_kmers , " kmers\n");
  auto all_num_ctgs = reduce_one(ctgs.size(), op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  SOUT("Processed a total of ", all_num_ctgs, " contigs and ", all_num_kmers, " kmers\n");
  barrier();
}


void add_ctg_kmers(unsigned kmer_len, Contigs &ctgs, dist_object<KmerDHT> &kmer_dht, bool use_bloom)
{
  Timer timer(__func__);
  int64_t num_kmers = 0;
  int64_t num_prev_kmers = kmer_dht->get_num_kmers();
  ProgressBar progbar(ctgs.size(), "Adding extra contig kmers");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() >= kmer_len + 2) {
      auto kmers = Kmer::get_kmers(ctg->seq);
      if (kmers.size() != ctg->seq.length() - kmer_len + 1)
        DIE("kmers size mismatch ", kmers.size(), " != ", (ctg->seq.length() - kmer_len + 1), " '", ctg->seq, "'");
      for (int i = 1; i < ctg->seq.length() - kmer_len; i++) {
        kmer_dht->add_kmer(kmers[i], ctg->seq[i - 1], ctg->seq[i + kmer_len], ctg->depth, CTG_KMERS_PASS);
        num_kmers++;
      }
    }
    progress();
  }
  progbar.done();
  kmer_dht->flush_updates(CTG_KMERS_PASS);
  DBG("This rank processed ", ctgs.size(), " contigs and ", num_kmers , " kmers\n");
  auto all_num_ctgs = reduce_one(ctgs.size(), op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  SOUT("Processed a total of ", all_num_ctgs, " contigs and ", all_num_kmers, " kmers\n");
  SOUT("Found ", perc_str(kmer_dht->get_num_kmers() - num_prev_kmers, all_num_kmers), " additional unique kmers\n");
}

