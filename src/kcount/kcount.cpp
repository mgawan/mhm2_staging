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
#include "options_kcount.hpp"
#include "progressbar.hpp"
#include "kmer_dht.hpp"
#include "kcount.h"
#include "dbjg_traversal.hpp"

using namespace std;
using namespace upcxx;

static const int QUAL_CUTOFF = 20;

extern ofstream _dbgstream;

uint64_t estimate_cardinality(shared_ptr<OptionsKcount> options)
{
  Timer timer(__func__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t num_kmers = 0;
  string read_record[4];
  int64_t estimated_total_records = 0;
  int64_t total_records_processed = 0;
  for (auto const &reads_fname : options->reads_fname_list) {
    int64_t records_processed = 0;
    bool isCompressed = hasEnding(reads_fname, ".gz");
    int64_t fileSize = get_file_size(reads_fname);
    zstr::ifstream reads_file(reads_fname);
    int64_t bytes_read = 0;
    ProgressBar progbar(reads_fname, &reads_file, "Scanning reads file to estimate cardinality");
    while (!reads_file.eof()) {
      bool done = false;
      for (int i = 0; i < 4; i++) {
        getline(reads_file, read_record[i]);
        if (read_record[i] == "") {
          done = true;
          break;
        }
        bytes_read += read_record[i].length();
        progbar.update(bytes_read);
        num_lines++;
      }
      if (records_processed++ > 100000) break; // do not read the entire data set for just an estimate
      if (done) break;
      num_reads++;
      if (read_record[1].length() < options->kmer_len) continue;
      num_kmers += read_record[1].length() - options->kmer_len + 1;
    }
    total_records_processed += records_processed;
    estimated_total_records += records_processed * fileSize / reads_file.zstr_tellg();
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

void count_kmers(shared_ptr<OptionsKcount> options, dist_object<KmerDHT> &kmer_dht, PASS_TYPE pass_type)
{
  Timer timer(__func__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  string read_record[4];
  int max_read_len = 0;
  int64_t num_kmers = 0;
  string progbar_prefix = "";
  switch (pass_type) {
    case BLOOM_SET_PASS: progbar_prefix = "Pass 1: Parsing reads file to setup bloom filter"; break;
    case BLOOM_COUNT_PASS: progbar_prefix = "Pass 2: Parsing reads file to count kmers"; break;
    case NO_BLOOM_PASS: progbar_prefix = "Parsing reads file to count kmers"; break;
  };
  char special = options->qual_offset + 2;
  vector<int> maxReadLengths;
  IntermittentTimer t_io("reads IO");
  for (auto const &reads_fname : options->reads_fname_list) {
    int64_t bytes_read = 0;
    int maxReadLength = 0;
    zstr::ifstream reads_file(reads_fname);
    /*
    stringstream reads_file_buf;
    {
      t_io.start();
      zstr::ifstream reads_file(reads_fname);
      reads_file_buf << reads_file.rdbuf();
      t_io.stop();
    }
    */
    ProgressBar progbar(reads_fname, &reads_file, progbar_prefix);
    while (!reads_file.eof()) {
      bool done = false;
      for (int i = 0; i < 4; i++) {
        getline(reads_file, read_record[i]);
        if (read_record[i] == "") {
          done = true;
          break;
        }
        bytes_read += read_record[i].length();
        progbar.update(bytes_read);
        num_lines++;
        if (i == 1 && maxReadLength < read_record[i].length()) maxReadLength = read_record[i].length();
      }
      if (done) break;
      string seq = move(read_record[1]);
      string quals = move(read_record[3]);
      num_reads++;
      if (seq.length() < options->kmer_len) continue;
      if (seq.length() > max_read_len) max_read_len = seq.length();

      // split into kmers
      auto kmers = Kmer::getKmers(seq);

      // disable kmer counting of kmers after a bad quality score (of 2) in the read
      // ... but allow extension counting (if an extention q score still passes the QUAL_CUTOFF)
      size_t foundBadQual = quals.find_first_of(special);
      if (foundBadQual == string::npos) foundBadQual = seq.length();   // remember that the last valid position is length()-1
      int foundBadQualKmer = foundBadQual - options->kmer_len + 1;
      assert( (int) kmers.size() >= foundBadQualKmer );

      // skip kmers that contain an N
      size_t foundN = seq.find_first_of('N');
      if (foundN == string::npos) foundN = seq.length();

      for (int i = 0; i < kmers.size(); i++) {
        // skip kmers that contain an N
        if (i + options->kmer_len > foundN) {
          i = foundN; // skip
          // find the next N
          foundN = seq.find_first_of('N', foundN+1);
          if (foundN == string::npos) foundN = seq.length();
          continue;
        }
        char left_base = '0';
        if (i > 0 && quals[i - 1] >= options->qual_offset + QUAL_CUTOFF) {
          left_base = seq[i - 1];
        }
        char right_base = '0';
        if (i + options->kmer_len < seq.length() && quals[i + options->kmer_len] >= options->qual_offset + QUAL_CUTOFF) {
          right_base = seq[i + options->kmer_len];
        }
        int count = (i < foundBadQualKmer) ? 1 : 0;
        kmer_dht->add_kmer(kmers[i], left_base, right_base, count, pass_type);
        num_kmers++;
      }
      progress();
    }
    progbar.done();
    kmer_dht->flush_updates(pass_type);
    maxReadLength = reduce_one(maxReadLength, op_max, 0).wait();
    if (!rank_me()) {
      SOUT("Detected maxReadLength of ", maxReadLength, " for ", reads_fname, "\n");
      string maxFile = reads_fname;
      size_t pos = maxFile.find_last_of('/');
      if (pos != string::npos) maxFile = maxFile.substr(pos+1);
      maxFile = "per_thread/" + maxFile + ".maxReadLen.txt";
      if (! does_file_exist(maxFile) ) {
        write_num_file(maxFile, maxReadLength);
      }
    }
  }
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads), max read length ", max_read_len, "\n");
  auto all_num_lines = reduce_one(num_lines, op_fast_add, 0).wait();
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto all_distinct_kmers = kmer_dht->get_num_kmers();
  auto all_max_read_len = reduce_one(max_read_len, op_fast_max, 0).wait();
  SOUT("Processed a total of ", all_num_lines, " lines (", all_num_reads, " reads), max read length ", all_max_read_len, "\n");
  if (pass_type != BLOOM_SET_PASS) SOUT("Found ", perc_str(all_distinct_kmers, all_num_kmers), " unique kmers\n");
}

void add_ctg_kmers(shared_ptr<OptionsKcount> options, dist_object<KmerDHT> &kmer_dht, bool use_bloom, int pass_num_mask)
{
  Timer timer(__func__);
  string ctgs_fname = options->cached_io ? "/dev/shm/" : "./";
  ctgs_fname += options->ctgs_fname;
  get_rank_path(ctgs_fname, rank_me());
  // first pass over ctgs file - count the number of kmers
  if ((pass_num_mask&1) == 1) {
    kmer_dht->num_prev_mers_from_ctgs = 0;
    zstr::ifstream ctgs_file(ctgs_fname);
    ProgressBar progbar(ctgs_fname, &ctgs_file, "Parsing contigs for extra kmers: pass 1");
    string cname, seq;
    int64_t bytes_read = 0;
    while (!ctgs_file.eof()) {
      getline(ctgs_file, cname);
      if (cname == "") break;
      getline(ctgs_file, seq);
      if (seq == "") break;
      bytes_read += cname.length() + seq.length();
      int64_t num_mers = seq.length() - options->kmer_len - 1;
      if (num_mers < 0) num_mers = 0;
      kmer_dht->num_prev_mers_from_ctgs += seq.length() - options->prev_kmer_len + 1;
      if (use_bloom && seq.length() >= options->kmer_len) {
        auto kmers = Kmer::getKmers(seq);
        if (kmers.size() != seq.length() - options->kmer_len + 1)
          DIE("kmers size mismatch ", kmers.size(), " != ", (seq.length() - options->kmer_len + 1), " '", seq, "'");
        for (int i = 1; i < seq.length() - options->kmer_len; i++) {
          kmer_dht->add_kmer(kmers[i], seq[i - 1], seq[i + options->kmer_len], 1, CTG_BLOOM_SET_PASS);
        }
      }
      progbar.update(bytes_read);
      progress();
    }
    if (use_bloom) kmer_dht->flush_updates(CTG_BLOOM_SET_PASS);
    progbar.done();
    DBG("This rank found ", kmer_dht->num_prev_mers_from_ctgs, " prev mers from contigs\n");
    auto all_num_prev_mers_from_ctgs = reduce_one(kmer_dht->num_prev_mers_from_ctgs, op_fast_add, 0).wait();
    SOUT("Found ", all_num_prev_mers_from_ctgs, " prev mers from contigs\n");
  }
  // second pass over contigs file - insert the new kmers into the dht
  if ((pass_num_mask&2) == 2) {
    // read all the kmer depths from the binary depths file
    string depths_fname(options->cached_io ? "/dev/shm/" : "./");
    depths_fname += options->ctg_depths_fname;
    get_rank_path(depths_fname, rank_me());
    zstr::ifstream depths_file(depths_fname, std::ios::binary);
    int64_t prefix_buf_sz = kmer_dht->num_prev_mers_from_ctgs * sizeof(int64_t);
    int64_t *prefix_buf = (int64_t*)malloc(kmer_dht->num_prev_mers_from_ctgs * sizeof(int64_t));
    depths_file.read((char*)prefix_buf, prefix_buf_sz);
    if (!depths_file)
      DIE("could not read all ", kmer_dht->num_prev_mers_from_ctgs, " kmers from depths file ", depths_fname,
          " only read ", depths_file.gcount(), " bytes");
    if (depths_file.gcount() != prefix_buf_sz) DIE("read ", depths_file.gcount(), " but expected ", prefix_buf_sz, "\n");
    string cname, seq;
    int64_t num_kmers = 0;
    int64_t num_ctgs = 0;
    int64_t num_prev_kmers = kmer_dht->get_num_kmers();
    int64_t bytes_read = 0;
    int64_t offset_in_buf = 0;
    zstr::ifstream ctgs_file(ctgs_fname);
    ProgressBar progbar(ctgs_fname, &ctgs_file, "Parsing contigs for extra kmers: pass 2");
    while (!ctgs_file.eof()) {
      getline(ctgs_file, cname);
      if (cname == "") break;
      getline(ctgs_file, seq);
      if (seq == "") break;
      bytes_read += cname.length() + seq.length();
      num_ctgs++;
      int64_t *cur_prefix_buf = (int64_t *)&prefix_buf[offset_in_buf];
      offset_in_buf += seq.length() - options->prev_kmer_len + 1;
      if (offset_in_buf > kmer_dht->num_prev_mers_from_ctgs) DIE("Offset in buf out of range ", offset_in_buf, " max ", kmer_dht->num_prev_mers_from_ctgs);
      if (seq.length() >= options->kmer_len + 2) {
        auto kmers = Kmer::getKmers(seq);
        if (kmers.size() != seq.length() - options->kmer_len + 1)
          DIE("kmers size mismatch ", kmers.size(), " != ", (seq.length() - options->kmer_len + 1), " '", seq, "'");
        for (int i = 1; i < seq.length() - options->kmer_len; i++) {
          int64_t depth_sum = cur_prefix_buf[i + options->kmer_len - options->prev_kmer_len] - cur_prefix_buf[i - 1];
          depth_sum = (int)((depth_sum + options->kmer_len - options->prev_kmer_len) / (options->kmer_len - options->prev_kmer_len + 1));
          if (depth_sum > numeric_limits<uint16_t>::max()) depth_sum = numeric_limits<uint16_t>::max();
          kmer_dht->add_kmer(kmers[i], seq[i - 1], seq[i + options->kmer_len], depth_sum, CTG_KMERS_PASS);
          num_kmers++;
        }
      }
      progbar.update(bytes_read);
      progress();
    }
    progbar.done();
    kmer_dht->flush_updates(CTG_KMERS_PASS);
    DBG("This rank processed ", num_ctgs, " contigs and ", num_kmers , " kmers\n");
    auto all_num_ctgs = reduce_one(num_ctgs, op_fast_add, 0).wait();
    auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
    SOUT("Processed a total of ", all_num_ctgs, " contigs and ", all_num_kmers, " kmers\n");
    SOUT("Found ", perc_str(kmer_dht->get_num_kmers() - num_prev_kmers, all_num_kmers), " additional unique kmers\n");
    free(prefix_buf);
  }
}

int kcount_main(int argc, char **argv)
{
  SOUT("max kmer size ", MAX_KMER_SIZE, "\n");
  auto start_t = chrono::high_resolution_clock::now();

  SOUT("Running on ", rank_n(), " processes\n");
  
#ifdef DBG_ON
  time_t curr_t = std::time(nullptr);
  string dbg_fname = "kcount" + to_string(curr_t) + ".dbg"; // never in cached_io
  get_rank_path(dbg_fname, rank_me());
  _dbgstream.open(dbg_fname);
  SOUT(KRED, "DEBUG mode - expect low performance\n", KNORM);
#endif

  double start_mem_free = get_free_mem_gb();
  SOUT("Initial free memory on node 0: ", start_mem_free, "GB\n");
    
  {
    Timer timer(__func__);
    shared_ptr<OptionsKcount> options = make_shared<OptionsKcount>();
    options->load(argc, argv);
    auto my_cardinality = estimate_cardinality(options);
    Kmer::set_k(options->kmer_len);
    dist_object<KmerDHT> kmer_dht(world(), my_cardinality, options->max_kmer_store, options->min_depth_cutoff,
                                  options->dynamic_min_depth, options->use_bloom);
    barrier();
    if (options->use_bloom) {
      count_kmers(options, kmer_dht, BLOOM_SET_PASS);
      if (options->ctgs_fname != "") {
        SOUT("Scanning contigs file to populate bloom2\n");
        add_ctg_kmers(options, kmer_dht, true, 1);
      }
      kmer_dht->reserve_space_and_clear_bloom1();
      count_kmers(options, kmer_dht, BLOOM_COUNT_PASS);
    } else {
      count_kmers(options, kmer_dht, NO_BLOOM_PASS);
    }
    barrier();
    SOUT("kmer DHT load factor: ", kmer_dht->load_factor(), "\n");
    barrier();
    kmer_dht->write_histogram();
    barrier();
    kmer_dht->purge_kmers(options->min_depth_cutoff);
    int64_t newCount = kmer_dht->get_num_kmers();
    SOUT("After purge of kmers <", options->min_depth_cutoff, " there are ", newCount, " unique kmers\n");
    barrier();
    if (options->ctgs_fname != "") {
      add_ctg_kmers(options, kmer_dht, options->use_bloom, options->use_bloom ? 2 : 3);
      kmer_dht->purge_kmers(1);
    }
    barrier();
    kmer_dht->compute_kmer_exts();
    kmer_dht->dump_kmers(options->kmer_len, options->cached_io);
    barrier();
    kmer_dht->purge_fx_kmers();
    traverse_debruijn_graph(options, kmer_dht);
    double end_mem_free = get_free_mem_gb();
    SOUT("Final free memory on node 0: ", end_mem_free, "GB, used ", (start_mem_free - end_mem_free), "GB\n");
    barrier();
  }
  Timer lastly("Reductions");
  auto tot_upc_mem_leak = reduce_one(upc_mem_alloced - upc_mem_freed, op_fast_add, 0).wait();
  auto tot_upc_mem_peak = reduce_one(upc_mem_peak, op_fast_add, 0).wait();
  SOUT("Peak UPC memory ", (double)tot_upc_mem_peak / ONE_GB / rank_n(), "GB per rank\n");
  if (tot_upc_mem_leak) SOUT("Apparent memory leak of ", tot_upc_mem_leak, " across all ranks\n");

  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SOUT("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(), "\n"); 
  barrier();

#ifdef DBG_ON
  _dbgstream.flush();
  _dbgstream.close();
#endif

  upcxx::barrier();
  
  return 0;
}

