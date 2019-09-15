#ifndef _CONTIGS_HPP
#define _CONTIGS_HPP

#include <iostream>
#include <vector>
#include <upcxx/upcxx.hpp>

#include "zstr.hpp"
#include "progressbar.hpp"

using std::vector;
using std::string;
using std::endl;

using upcxx::rank_me;
using upcxx::reduce_one;
using upcxx::op_fast_add;
using upcxx::barrier;


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
        out_buf << ">Contig" << uutig->id << " " << (uutig->seq.length() - kmer_len + 1) << " " << uutig->depth << endl;
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
      SOUT("Wrote ", reduce_one(contigs.size(), op_fast_add, 0).wait(), " contigs to ", fname, "\n");
      write_uncompressed_file_size(fname + ".uncompressedSize", bytes_written);
    }
    barrier();
    return fname;
  }

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
    SOUT("Processed a total of ", all_num_ctgs, " contigs\n");
    barrier();
  }
  
};
 
#endif
