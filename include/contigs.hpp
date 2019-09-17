#ifndef _CONTIGS_HPP
#define _CONTIGS_HPP

#include <iostream>
#include <vector>
#include <upcxx/upcxx.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "zstr.hpp"
#include "progressbar.hpp"

using std::pair;
using std::vector;
using std::string;
using std::endl;
using std::max;
using std::memory_order_relaxed;

using upcxx::rank_me;
using upcxx::rank_n;
using upcxx::reduce_one;
using upcxx::op_fast_add;
using upcxx::op_fast_max;
using upcxx::barrier;
using upcxx::atomic_domain;
using upcxx::atomic_op;
using upcxx::global_ptr;
using upcxx::new_;

struct Contig {
  int64_t id;
  string seq;
  double depth;
};

class Contigs {
  
  vector<Contig> contigs;
  
public:

  Contigs() {}

  void clear() {
    contigs.clear();
    vector<Contig>().swap(contigs);
  }
    
  void add_contig(Contig &contig) {
    contigs.push_back(contig);
  }

  size_t size() {
    return contigs.size();
  }
  
  auto begin() {
    return contigs.begin();
  }

  auto end() {
    return contigs.end();
  }

  /*
  string dump_contigs(string fname, unsigned kmer_len) {
    fname += "-" + to_string(kmer_len) + ".fasta.gz";
    get_rank_path(fname, rank_me());
    {
      zstr::ofstream f(fname);
      ostringstream out_buf;
      ProgressBar progbar(contigs.size(), "Writing contigs");
      size_t bytes_written = 0;
      int64_t i = 0;
      for (auto it = contigs.begin(); it != contigs.end(); ++it) {
        auto uutig = it;
        out_buf << ">Contig" << uutig->id << " " << << uutig->depth << endl;
        string rc_uutig = revcomp(uutig->seq);
        if (rc_uutig < uutig->seq) uutig->seq = rc_uutig;
        // fold output
        for (int p = 0; p < uutig->seq.length(); p += 50)
          out_buf << uutig->seq.substr(p, 50) << endl;
        progbar.update();
        i++;
        if (!(i % 1000)) {
          f << out_buf.str();
          bytes_written += out_buf.str().length();
          out_buf = ostringstream();
        }
      }
      if (!out_buf.str().empty()) {
        f << out_buf.str();
        bytes_written += out_buf.str().length();
      }
      f.close();
      progbar.done();
      SLOG("Wrote ", reduce_one(contigs.size(), op_fast_add, 0).wait(), " contigs to ", fname, "\n");
      write_uncompressed_file_size(fname + ".uncompressedSize", bytes_written);
    }
    barrier();
    return fname;
  }
  */
  
  void print_stats(int min_ctg_len) {
    int64_t tot_len = 0, max_len = 0;
    double tot_depth = 0;
    vector<pair<int, int64_t>> length_sums({ {1, 0}, {5, 0}, {10, 0}, {25, 0}, {50, 0}});
    int64_t num_ctgs = 0;
    int64_t num_ns = 0;
    vector<int> lens;
    lens.reserve(contigs.size());
    for (auto ctg : contigs) {
      auto len = ctg.seq.length();
      if (len < min_ctg_len) continue;
      num_ctgs++;
      tot_len += len;
      tot_depth += ctg.depth;
      max_len = max(max_len, static_cast<int64_t>(len));
      for (auto &length_sum : length_sums) {
        if (len >= length_sum.first * 1000) length_sum.second += len;
      }
      num_ns += count(ctg.seq.begin(), ctg.seq.end(), 'N');
      lens.push_back(len);
    }
    // Compute local N50 and then take the median across all of them. This gives a very good approx of the exact N50 and is much
    // cheaper to compute
    sort(lens.rbegin(), lens.rend());
    int64_t sum_lens = 0;
    int64_t n50 = 0;
    for (auto l : lens) {
      sum_lens += l;
      if (sum_lens >= tot_len / 2) {
        n50 = l;
        break;
      }
    }
    barrier();
    int64_t all_num_ctgs = reduce_one(num_ctgs, op_fast_add, 0).wait();
    int64_t all_tot_len = reduce_one(tot_len, op_fast_add, 0).wait();
    int64_t all_max_len = reduce_one(max_len, op_fast_max, 0).wait();
    double all_tot_depth = reduce_one(tot_depth, op_fast_add, 0).wait();
    int64_t all_num_ns = reduce_one(num_ns, op_fast_add, 0).wait();
    int64_t all_n50s = reduce_one(n50, op_fast_add, 0).wait();

    SLOG("Assembly statistics (contig lengths >= ", min_ctg_len, ")\n");
    SLOG("    Number of contigs:       ", all_num_ctgs, "\n");
    SLOG("    Total assembled length:  ", all_tot_len, "\n");
    SLOG("    Average contig depth:    ", all_tot_depth / all_num_ctgs, "\n");
    SLOG("    Number of Ns/100kbp:     ", (double)all_num_ns * 100000.0 / all_tot_len, " (", all_num_ns, ")", KNORM, "\n");
    // FIXME: need to this to be the median, not average - median is much more reliable
    SLOG("    Approx N50:              ", all_n50s / rank_n(), " (rank 0 only ", n50, ")\n");
    SLOG("    Max. contig length:      ", all_max_len, "\n");
    SLOG("    Contig lengths:\n");
    for (auto &length_sum : length_sums) {
      SLOG("        > ", std::left, std::setw(19), to_string(length_sum.first) + "kbp:", 
           perc_str(reduce_one(length_sum.second, op_fast_add, 0).wait(), all_tot_len), "\n");
    }
  }

  void dump_contigs(const string &fname, int min_ctg_len) {
    Timer timer(__func__);
    string tmpfname = fname + ".tmp"; // make a .tmp file and rename on success
    string fasta = "";
    for (auto it = contigs.begin(); it != contigs.end(); ++it) {
      auto ctg = it;
      if (ctg->seq.length() < min_ctg_len) continue;
      fasta += ">Contig" + to_string(ctg->id) + " " + to_string(ctg->depth) + "\n";
      string rc_uutig = revcomp(ctg->seq);
      string seq = (rc_uutig < ctg->seq ? rc_uutig : ctg->seq);
      //for (int64_t i = 0; i < ctg->seq.length(); i += 50) fasta += ctg->seq.substr(i, 50) + "\n";
      fasta += ctg->seq + "\n";
    }
    auto sz = fasta.size();
    atomic_domain<size_t> ad({atomic_op::fetch_add, atomic_op::load});
    global_ptr<size_t> fpos = nullptr;
    if (!rank_me()) fpos = new_<size_t>(0);
    fpos = broadcast(fpos, 0).wait();
    size_t my_fpos = ad.fetch_add(fpos, sz, memory_order_relaxed).wait();
    // wait until all ranks have updated the global counter
    barrier();
    int bytes_written = 0;
    int fileno = -1;
    size_t fsize = 0;
    if (!rank_me()) {
      fsize = ad.load(fpos, memory_order_relaxed).wait();
      // rank 0 creates the file and truncates it to the correct length
      fileno = open(tmpfname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
      if (fileno == -1) WARN("Error trying to create file ", tmpfname, ": ", strerror(errno), "\n");
      if (ftruncate(fileno, fsize) == -1) WARN("Could not truncate ", tmpfname, " to ", fsize, " bytes\n");
    }
    barrier();
    ad.destroy();
    // wait until rank 0 has finished setting up the file
    if (rank_me()) fileno = open(tmpfname.c_str(), O_WRONLY, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
    if (fileno == -1) WARN("Error trying to open file ", tmpfname, ": ", strerror(errno), "\n");
    bytes_written = pwrite(fileno, fasta.c_str(), sz, my_fpos);
    close(fileno);

    if (bytes_written != sz)
      DIE("Could not write all ", sz, " bytes; only wrote ", bytes_written, "\n");
    barrier();
    if (rank_me() == 0) {
      string new_fname = fname + ".fasta";
      if (rename(tmpfname.c_str(), new_fname.c_str()) != 0)
        SDIE("Could not rename ", tmpfname, " to ", new_fname);
      SLOG_VERBOSE("Successfully wrote ", fsize, " bytes to ", new_fname, "\n");
    }
  }

  /*
  void load_contigs(string ctgs_fname) {
    string cname, seq, buf;
    size_t bytes_read = 0;
    zstr::ifstream ctgs_file(ctgs_fname);
    ProgressBar progbar(ctgs_fname, &ctgs_file, "Parsing contigs");
    while (!ctgs_file.eof()) {
      getline(ctgs_file, cname);
      if (cname == "") break;
      // seq is folded, so read until we encounter '>' or eof
      seq = "";
      while (true) {
        getline(ctgs_file, buf);
        if (buf == "") break;
        seq += buf;
        auto next_ch = ctgs_file.peek();
        if (next_ch == (int)'>' || next_ch == EOF) break;
      }
      if (seq == "") break;
      bytes_read += cname.length() + seq.length();
      // extract the id
      size_t firstspace = cname.find_first_of(" ");
      if (firstspace == std::string::npos) DIE("Ctgs file ", ctgs_fname, " is incorrect format on line: '", cname, "'");
      int64_t id = stol(cname.substr(7, firstspace));
      // depth is the last field in the cname
      size_t lastspace = cname.find_last_of(" ");
      if (lastspace == std::string::npos) DIE("Depth is missing from ctgs file ", ctgs_fname, " on line: '", cname, "'");
      double depth = stod(cname.substr(lastspace));
      Contig contig = { .id = id, .seq = seq, .depth = depth };
      add_contig(contig);
      progbar.update(bytes_read);
    }
    progbar.done();
    DBG("This rank processed ", contigs.size(), " contigs\n");
    auto all_num_ctgs = reduce_one(contigs.size(), op_fast_add, 0).wait();
    SLOG("Processed a total of ", all_num_ctgs, " contigs\n");
    barrier();
  }
  */
/*
  void load_contigs(const string &ctgs_fname) {
  auto get_file_offset_for_rank = [](ifstream &f, int rank, string &ctg_prefix) -> size_t {
    f.seekg (0, f.end);
    auto sz = f.tellg();
    if (rank == 0) return 0;
    if (rank == rank_n()) return sz;
    size_t offset = sz / rank_n() * rank;
    f.seekg(offset);
    string line;
    while (getline(f, line)) {
      if (line.substr(0, ctg_prefix.size()) == ctg_prefix) {
        getline(f, line);
        break;
      }
    }
    return f.tellg();
  };

  Timer timer(__func__);
  string line;
  string ctg_prefix = ">Contig_";
  int64_t num_ctgs_found = 0;
  cid_t max_cid = 0;
  int64_t tot_len = 0;
  int64_t tot_num_kmers = 0;
  KmerHashTable *kmer_ht = new KmerHashTable(options.kmer_len);
  {
    ifstream ctgs_file(options.ctgs_fname);
    if (!ctgs_file.is_open()) DIE("Could not open ctgs file '", options.ctgs_fname, "': ", strerror(errno));
    ProgressBar progbar(&ctgs_file, "Parsing ctgs file", options.one_file_per_rank);
    int64_t start_rank = options.one_file_per_rank ? 0 : rank_me();
    int64_t stop_rank = options.one_file_per_rank ? rank_n() : rank_me() + 1;
    auto start_offset = get_file_offset_for_rank(ctgs_file, start_rank, ctg_prefix);
    auto stop_offset = get_file_offset_for_rank(ctgs_file, stop_rank, ctg_prefix);
    ctgs_file.seekg(start_offset);
#ifdef ASYNC
    promise<> prom;
#endif
    int64_t num_ctgs = 0;
    string seq;
    while (getline(ctgs_file, line)) {
      if (line.substr(0, ctg_prefix.size()) == ctg_prefix) {
        cid_t cid = stol(line.substr(ctg_prefix.size()));
        max_cid = max(max_cid, cid);
        num_ctgs_found++;
        // now read in the sequence - all on a single line
        getline(ctgs_file, seq);
        tot_len += seq.length();
        // skip any that are shorter than the seed
        if (options.kmer_len > seq.length()) continue;
        global_ptr<char> seq_gptr = allocate<char>(seq.length() + 1);
        strcpy(seq_gptr.local(), seq.c_str());
        CtgLoc ctg_loc = { .cid = cid, .seq_gptr = seq_gptr, .clen = (int)seq.length() };
        auto num_kmers = seq.length() - options.kmer_len + 1;
        tot_num_kmers += num_kmers;
        for (int i = 0; i < num_kmers; i++) {
          string&& kmer = seq.substr(i, options.kmer_len);
          string kmer_rc;
          ctg_loc.pos_in_ctg = i;
#ifndef ASYNC
          promise<> prom;
#endif
          if (!cond_revcomp(kmer, &kmer_rc)) {
            ctg_loc.is_rc = false;
            kmer_ht->add_kmer(kmer, ctg_loc, prom);
          } else {
            ctg_loc.is_rc = true;
            kmer_ht->add_kmer(kmer_rc, ctg_loc, prom);
          }
#ifndef ASYNC
          prom.finalize().wait();
#endif
          progress();
        }
        num_ctgs++;
      }
      progbar.update();
      if (ctgs_file.tellg() >= stop_offset) break;
    }
#ifdef ASYNC
    prom.finalize().wait();
#endif
    progbar.done();
    barrier();
  }
  SOUT("Found ", reduce_one(num_ctgs_found, op_fast_add, 0).wait(), " contigs, max id ",
       reduce_one(max_cid, op_fast_max, 0).wait(), "\n");
  SOUT("Total length ", reduce_one(tot_len, op_fast_add, 0).wait(), "\n");

  auto num_kmers_added = reduce_one(tot_num_kmers, op_fast_add, 0).wait();
  auto num_kmers_in_ht = kmer_ht->get_num_kmers();
  SOUT("Processed ", num_kmers_added, " kmers\n");
  auto num_dropped = kmer_ht->get_num_dropped();
  if (num_dropped) 
    SOUT("Dropped ", num_dropped, " kmer-to-contig mappings (", 
         setprecision(2), fixed, (100.0 * num_dropped / num_kmers_added), "%)\n");
  return kmer_ht;
*/
};
 
#endif
