#include <iostream>
#include <fstream>
#include <regex>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "progressbar.hpp"
#include "zstr.hpp"
#include "ctg_graph.hpp"
#include "contigs.hpp"
#include "alignments.hpp"
#include "fastq.hpp"

using namespace std;
using namespace upcxx;


// contig is the ref (subject), read is the query - done this way so that we can potentially do multiple alns to each read
/*
struct Aln {
  string read_id;
  int64_t cid; 
  // query
  int rstart, rstop, rlen;
  // ref
  int cstart, cstop, clen;
  char orient;
  int score1, score2;

  // writes out in the format meraligner uses
  string to_string() {
    ostringstream os;
    os << read_id << "\t" << rstart + 1 << "\t" << rstop << "\t" << rlen << "\t"
       << "Contig" << cid << "\t" << cstart + 1 << "\t" << cstop << "\t" << clen << "\t"
       << (orient == '+' ? "Plus" : "Minus") << "\t" << score1 << "\t" << score2;
    return os.str();
  }
}
*/


static bool get_best_span_aln(vector<Aln> &alns, Aln &best_aln, string &readStatus, string &typeStatus,
                              int64_t *reject_5_trunc, int64_t *reject_3_trunc, int64_t *reject_uninf) {
  if (alns.size() == 0) return false;
  int min_rstart = -1;
  readStatus = "";
  // sort in order of highest to lowest aln scores
  sort(alns.begin(), alns.end(),
       [](auto &aln1, auto &aln2) {
         return aln1.score1 > aln2.score1;
       });
  for (auto aln : alns) {
    //if (aln.orient == '-') switch_orient(aln.rstart, aln.rstop, aln.rlen);
    // Assess alignment for completeness (do this before scaffold coordinate conversion!)
    // Important: Truncations are performed before reverse complementation
    // and apply to the end of the actual read
    int rstart = aln.rstart + 1;
    int cstart = aln.cstart + 1;
    int unalignedStart = rstart - 1;
    int projectedStart = (aln.orient == '+' ? cstart - unalignedStart : aln.cstop + unalignedStart);
    string startStatus = "";
    int projectedOff = 0;
    if (projectedStart < 1) projectedOff = 1 - projectedStart;
    else if (projectedStart > aln.clen) projectedOff = projectedStart - aln.clen;
    int missingStartBases = unalignedStart - projectedOff;
    if (unalignedStart == 0) {
      startStatus = "FULL";
    } else if (projectedOff > 0 && missingStartBases < FIVE_PRIME_WIGGLE_ROOM) {
      startStatus = "GAP";
    } else if (unalignedStart < FIVE_PRIME_WIGGLE_ROOM) {
      startStatus = "INC";
    } else {
      (*reject_5_trunc)++;
      continue;
    }

    int unalignedEnd = aln.rlen - aln.rstop;
    int projectedEnd = (aln.orient == '+' ? aln.cstop + unalignedEnd : cstart - unalignedEnd);
    string endStatus = "";
    projectedOff = 0;
    if (projectedEnd < 1) projectedOff = 1 - projectedEnd;
    else if (projectedEnd > aln.clen) projectedOff = projectedEnd - aln.clen;
    int missingEndBases = unalignedEnd - projectedOff;
    if (unalignedEnd == 0) {
      endStatus = "FULL";
    } else if (projectedOff > 0 && missingEndBases < THREE_PRIME_WIGGLE_ROOM) {
      endStatus = "GAP";
    } else if (unalignedEnd < THREE_PRIME_WIGGLE_ROOM) {
      endStatus = "INC";
    } else {
      (*reject_3_trunc)++;
      continue;
    }

    int endDistance = INSERT_AVG + 3 * INSERT_STDDEV;
    // Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle)
    if ((aln.orient == '+' && projectedStart > aln.clen - endDistance) ||
        (aln.orient == '-' && projectedStart < endDistance)) {
      if (min_rstart == -1 || rstart < min_rstart) {
        readStatus = startStatus + "." + endStatus + ".OUT";
        if (startStatus == "GAP") typeStatus = "OUTGAP";
        else if (endStatus == "GAP") typeStatus = "INTGAP";
        else typeStatus = "ANCHOR";
        best_aln = aln;
        min_rstart = rstart;
      } else {
        (*reject_uninf)++;
      }
    } else {
      (*reject_uninf)++;
    }
  }
  if (min_rstart == -1) return false;
  return true;
}


static string get_ctg_aln_str(Aln &aln, const string &read_status) {
  ostringstream os;
  os << "Contig" << aln.cid << "." << (aln.orient == '+' ? 3 : 5) << "\t["
     << read_status << " " << aln.read_id << " " << (aln.rstart + 1) << " " << aln.rstop << " " << aln.rlen
     << " Contig" << aln.cid << " " << (aln.cstart + 1) << " " << aln.cstop << " " << aln.clen << " "
     << (aln.orient == '+' ? "Plus" : "Minus") << "]";
  return os.str();
}


// gets all the alns for a single read, and returns true if there are more alns
static bool get_all_alns_for_read(Alns &alns, int64_t &i, vector<Aln> &alns_for_read) {
  string start_read_id = "";
  for (; i < alns.size(); i++) {
    Aln aln = alns.get_aln(i);
    // alns for a new read
    if (start_read_id != "" && aln.read_id != start_read_id) return true;
    alns_for_read.push_back(aln);
    start_read_id = aln.read_id;
  }
  return false;
}


// so the way meraculous spanner seems to work is that it only makes a pair out of the first (best) alignments
// for a read. It relies on the alignments for a read being sorted by score, best to worst. It builds a list of
// all viable alns for each pair separately, and then takes the best from each list and uses them as the final pair

// WRONG: the first aln is the one with the lowest rstart!

void run_spanner(int max_kmer_len, int kmer_len, Alns &alns) {
  Timer timer(__func__);
  IntermittentTimer t_get_alns("get alns spanner");
  ProgressBar progbar(alns.size(), "Running spanner");
  int64_t reject_5_trunc = 0;
  int64_t reject_3_trunc = 0;
  int64_t reject_uninf = 0;
  int64_t num_pairs = 0;
  int64_t aln_i = 0;
  Aln prev_best_aln = { .read_id = "" };
  string read_status = "", type_status = "", prev_read_status = "", prev_type_status = "";
  string dump_fname = "spanner-" + to_string(kmer_len) + ".spans.gz";
  get_rank_path(dump_fname, rank_me());
  zstr::ofstream spans_file(dump_fname);
  while (true) {
    vector<Aln> alns_for_read;
    t_get_alns.start();
    if (!get_all_alns_for_read(alns, aln_i, alns_for_read)) break;
    t_get_alns.stop();
    progbar.update(aln_i);
    Aln best_aln = { .read_id = "" };
    if (get_best_span_aln(alns_for_read, best_aln, read_status, type_status, &reject_5_trunc, &reject_3_trunc, &reject_uninf)) {
      if (!prev_best_aln.read_id.empty()) {
        auto read_id_len = best_aln.read_id.length();
        if (best_aln.read_id.compare(0, read_id_len - 2, prev_best_aln.read_id, 0, read_id_len - 2) == 0) {
          if (best_aln.cid == prev_best_aln.cid) {
            reject_uninf++;
          } else {
            num_pairs++;
            //add_span(max_kmer_len, kmer_len, best_aln_pre, best_aln);
            spans_file << "PAIR\t" << prev_type_status << "." << type_status << "\t"
                       << get_ctg_aln_str(prev_best_aln, prev_read_status) << "\t"
                       << get_ctg_aln_str(best_aln, read_status) << endl;
            // there will be no previous one next time 
            prev_best_aln.read_id = "";
            prev_read_status = "";
            prev_type_status = "";
            continue;
          }
        }
      }
    }
    prev_best_aln = best_aln;
    prev_read_status = read_status;
    prev_type_status = type_status;
  }
  progbar.done();
  spans_file.close();
  barrier();
  exit(0);
}
