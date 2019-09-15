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
};

  
#endif

