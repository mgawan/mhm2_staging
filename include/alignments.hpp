#ifndef _ALIGNMENT_HPP
#define _ALIGNMENT_HPP


struct Aln {
  string read_id;
  int64_t cid;
  int rstart, rstop, rlen, cstart, cstop, clen;
  char orient;
  int score1, score2;
  
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
  
public:

  Alns() {}

  void clear() {
    alns.clear();
    vector<Aln>().swap(alns);
  }
    
  void add_aln(Aln &aln) {
    alns.push_back(aln);
  }

  Aln &get_aln(int64_t i) {
    return alns[i];
  }
  
  size_t size() {
    return alns.size();
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
    for (auto &aln : alns) {
      out_buf << aln.to_string() << endl;
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
    barrier();
  }
};

  
#endif

