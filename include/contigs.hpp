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
using upcxx::dist_object;


struct Contig {
  int64_t id;
  string seq;
  double depth;
  //vector<uint16_t> kmer_depths;

  /*
  uint16_t get_kmer_depth(int start_pos, int kmer_len, int prev_kmer_len) {
    int len_diff = kmer_len - prev_kmer_len;
    int d = 0;
    for (int i = start_pos; i < start_pos + len_diff; i++) {
      assert(i < kmer_depths.size());
      d += kmer_depths[i];
    }
    d /= len_diff;
    return d;
  }
  */
};

class Contigs {
  
  vector<Contig> contigs;
  
public:

  Contigs() {}

  void clear() {
    contigs.clear();
    vector<Contig>().swap(contigs);
  }

  void set_capacity(int64_t sz) {
    contigs.reserve(sz);
  }
  
  void add_contig(Contig contig) {
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
    dist_object<int64_t> n50_val(n50);
    int64_t median_n50 = 0;
    if (!rank_me()) {
      vector<int64_t> n50_vals;
      n50_vals.push_back(n50);
      for (int i = 1; i < rank_n(); i++) {
        n50_vals.push_back(n50_val.fetch(i).wait());
        //SOUT(i, " n50 fetched ", n50_vals.back(), "\n");
      }
      sort(n50_vals.begin(), n50_vals.end());
      median_n50 = n50_vals[rank_n() / 2];
    }
    // barrier to ensure the other ranks dist_objects don't go out of scope before rank 0 is done
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
    //SLOG("    Approx N50 (average):    ", all_n50s / rank_n(), " (rank 0 only ", n50, ")\n");
    SLOG("    Approx N50:              ", median_n50, "\n");
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
    contigs.clear();
    string line;
    string ctg_prefix = ">Contig";
    string cname, seq, buf;
    size_t bytes_read = 0;
    ifstream ctgs_file(ctgs_fname);
    if (!ctgs_file.is_open()) DIE("Could not open ctgs file '", ctgs_fname, "': ", strerror(errno));
    int64_t start_rank = rank_me();
    int64_t stop_rank = rank_me() + 1;
    auto start_offset = get_file_offset_for_rank(ctgs_file, start_rank, ctg_prefix);
    auto stop_offset = get_file_offset_for_rank(ctgs_file, stop_rank, ctg_prefix);
    ProgressBar progbar(stop_offset - start_offset, "Parsing contigs");
    ctgs_file.seekg(start_offset);
    int64_t tot_len = 0;
    while (!ctgs_file.eof()) {
      getline(ctgs_file, cname);
      if (cname == "") break;
      getline(ctgs_file, seq);
      if (seq == "") break;
      tot_len += seq.length();
      bytes_read += cname.length() + seq.length();
      progbar.update(bytes_read);
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
      if (ctgs_file.tellg() >= stop_offset) break;
    }
    progbar.done();
    barrier();
    SOUT("Found ", reduce_one(contigs.size(), op_fast_add, 0).wait(), " contigs\n");
    SOUT("Total length ", reduce_one(tot_len, op_fast_add, 0).wait(), "\n");
  }
  
};
 
#endif
