#ifndef _ALIGNMENT_HPP
#define _ALIGNMENT_HPP

#include <fcntl.h>
#include <upcxx/upcxx.hpp>
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "zstr.hpp"

using namespace upcxx_utils;

struct Aln {
  string read_id;
  int64_t cid;
  int rstart, rstop, rlen, cstart, cstop, clen;
  char orient;
  int score1, score2;
  string sam_string;

  // writes out in the format meraligner uses
  string to_string() {
    ostringstream os;
    os << read_id << "\t" << rstart + 1 << "\t" << rstop << "\t" << rlen << "\t"
       << "Contig" << cid << "\t" << cstart + 1 << "\t" << cstop << "\t" << clen << "\t"
       << (orient == '+' ? "Plus" : "Minus") << "\t" << score1 << "\t" << score2;
    return os.str();
  }
};


class Alns {

  vector<Aln> alns;
  int64_t num_dups;

public:

  Alns() : num_dups(0) {}

  void clear() {
    alns.clear();
    vector<Aln>().swap(alns);
  }

  void add_aln(Aln &aln) {
    // check for duplicate alns first - do this backwards because only the most recent entries could be for this read
    for (auto it = alns.rbegin(); it != alns.rend(); ++it) {
      // we have no more entries for this read
      if (it->read_id != aln.read_id) break;
      // now check for equality
      if (it->rstart == aln.rstart && it->rstop == aln.rstop && it->cstart == aln.cstart && it->cstop == aln.cstop) {
        num_dups++;
        return;
      }
    }
    alns.push_back(aln);
  }

  Aln &get_aln(int64_t i) {
    return alns[i];
  }

  size_t size() {
    return alns.size();
  }

  int64_t get_num_dups() {
    return upcxx::reduce_one(num_dups, upcxx::op_fast_add, 0).wait();
  }

  auto begin() {
    return alns.begin();
  }

  auto end() {
    return alns.end();
  }

  void dump_alns(string fname) {
    get_rank_path(fname, rank_me());
    zstr::ofstream f(fname);
    ostringstream out_buf;
    ProgressBar progbar(alns.size(), "Writing alns");
    size_t bytes_written = 0;
    int64_t i = 0;
    for (auto aln : alns) {
      out_buf << aln.to_string() << std::endl;
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
    upcxx::barrier();
  }

  void dump_single_file_alns(const string fname, bool as_sam_format=false) {
    BarrierTimer timer(__FILEFUNC__, false, true);

    string out_str = "";
    for (auto aln : alns) {
      if (!as_sam_format) out_str += aln.to_string() + "\n";
      else out_str += aln.sam_string + "\n";
    }
    upcxx::atomic_domain<size_t> ad({upcxx::atomic_op::fetch_add, upcxx::atomic_op::load});
    upcxx::global_ptr<size_t> fpos = nullptr;
    if (!upcxx::rank_me()) fpos = upcxx::new_<size_t>(0);
    fpos = upcxx::broadcast(fpos, 0).wait();
    auto sz = out_str.length();
    size_t my_fpos = ad.fetch_add(fpos, sz, std::memory_order_relaxed).wait();
    // wait until all ranks have updated the global counter
    upcxx::barrier();
    int fileno = -1;
    size_t fsize = 0;
    if (!upcxx::rank_me()) {
      fsize = ad.load(fpos, std::memory_order_relaxed).wait();
      // rank 0 creates the file and truncates it to the correct length
      fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
      if (fileno == -1) WARN("Error trying to create file ", fname, ": ", strerror(errno), "\n");
      if (ftruncate(fileno, fsize) == -1) WARN("Could not truncate ", fname, " to ", fsize, " bytes\n");
    }
    upcxx::barrier();
    ad.destroy();
    // wait until rank 0 has finished setting up the file
    if (rank_me()) fileno = open(fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) WARN("Error trying to open file ", fname, ": ", strerror(errno), "\n");
    auto bytes_written = pwrite(fileno, out_str.c_str(), sz, my_fpos);
    close(fileno);
    if (bytes_written != sz) DIE("Could not write all ", sz, " bytes; only wrote ", bytes_written, "\n");
    upcxx::barrier();
    auto tot_bytes_written = upcxx::reduce_one(bytes_written, upcxx::op_fast_add, 0).wait();
    upcxx::barrier();
    SLOG_VERBOSE("Successfully wrote ", get_size_str(tot_bytes_written), " bytes to ", fname, "\n");
  }
};


#endif

