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

void get_spans_from_alns(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, Alns &alns, CtgGraph *graph);


struct AlnStats {
  int64_t nalns, unaligned, short_alns, containments, circular, mismatched, conflicts, bad_overlaps, empty_spans;
  
  string to_string() const {
    int64_t tot_nalns = reduce_one(nalns, op_fast_add, 0).wait();
    int64_t tot_unaligned = reduce_one(unaligned, op_fast_add, 0).wait();
    int64_t tot_containments = reduce_one(containments, op_fast_add, 0).wait();
    int64_t tot_short_alns = reduce_one(short_alns, op_fast_add, 0).wait();
    int64_t tot_circular = reduce_one(circular, op_fast_add, 0).wait();
    int64_t tot_mismatched = reduce_one(mismatched, op_fast_add, 0).wait();
    int64_t tot_conflicts = reduce_one(conflicts, op_fast_add, 0).wait();
    int64_t tot_bad_overlaps = reduce_one(bad_overlaps, op_fast_add, 0).wait();
    int64_t tot_empty_spans = reduce_one(empty_spans, op_fast_add, 0).wait();
    int64_t nbad = tot_unaligned + tot_containments + tot_short_alns + tot_circular +
      tot_mismatched + tot_conflicts + tot_bad_overlaps;
    
    ostringstream os;
    os.precision(2);
    os << fixed;
    os << "Edge stats:\n";
    //os << "    used:               " << perc_str(tot_nalns - nbad, tot_nalns) << endl;
    os << "    unaligned:          " << perc_str(tot_unaligned, tot_nalns) << endl;
    os << "    containments:       " << perc_str(tot_containments, tot_nalns) << endl;
    os << "    short alns:         " << perc_str(tot_short_alns, tot_nalns) << endl;
    os << "    circular:           " << perc_str(tot_circular, tot_nalns) << endl;
    os << "    mismatched:         " << perc_str(tot_mismatched, tot_nalns) << endl;
    os << "    conflicts:          " << perc_str(tot_conflicts, tot_nalns) << endl;
    os << "    bad overlaps:       " << perc_str(tot_bad_overlaps, tot_nalns) << endl;
    os << "    empty spans:        " << perc_str(tot_empty_spans, tot_nalns) << endl;
    return os.str();
  }
};


static CtgGraph *_graph = nullptr;


static void add_vertices_from_ctgs(Contigs *ctgs) {
  Timer timer(__func__);
  ProgressBar progbar(ctgs->size(), "Adding contig vertices to graph");
  for (auto ctg : *ctgs) {
    Vertex v = { .cid = ctg.id, .clen = (int)ctg.seq.length(), .depth = ctg.depth };
    _graph->add_vertex(v, ctg.seq);
    progbar.update();
  }
  progbar.done();
  barrier();
  auto num_vertices = _graph->get_num_vertices();
  SLOG_VERBOSE("Added ", num_vertices, " vertices\n");
}


// gets all the alns for a single read, and returns true if there are more alns
static bool get_alns_for_read(Alns &alns, int64_t &i, vector<Aln> &alns_for_read, int64_t *nalns, int64_t *unaligned) {
  string start_read_id = "";
  for (; i < alns.size(); i++) {
    Aln aln = alns.get_aln(i);
    // alns for a new read
    if (start_read_id != "" && aln.read_id != start_read_id) return true;
    (*nalns)++;
    // convert to coords for use here
    if (aln.orient == '-') {
      int tmp = aln.cstart;
      aln.cstart = aln.clen - aln.cstop;
      aln.cstop = aln.clen - tmp;
    }
    int unaligned_left = min(aln.rstart, aln.cstart);
    int unaligned_right = min(aln.rlen - aln.rstop, aln.clen - aln.cstop);
    if (unaligned_left > ALN_SLOP || unaligned_right > ALN_SLOP) {
      (*unaligned)++;
      DBG("unaligned ", aln.read_id, " ", aln.rstart, " ", aln.rstop, " ", aln.rlen, " ",
          aln.cid, " ", aln.cstart, " ", aln.cstop, " ", aln.clen, " ", aln.orient, " ", aln.score1, " ", aln.score2, "\n");
    } else {
      alns_for_read.push_back(aln);
    }
    start_read_id = aln.read_id;
    progress();
  }
  return false;
}


static bool add_splint(const Aln *aln1, const Aln *aln2, AlnStats &stats) {
  struct AlnCoords {
    int start, stop;
  };

  auto get_aln_coords = [](const Aln *aln) -> AlnCoords {
    // get contig start and end positions in read coordinate set
    return { .start = aln->rstart - aln->cstart, .stop = aln->rstop + (aln->clen - aln->cstop) };
  };

  AlnCoords ctg1 = get_aln_coords(aln1);
  AlnCoords ctg2 = get_aln_coords(aln2);
  
  if (aln1->cid == aln2->cid) {
    stats.circular++;
    return false;
  }
  auto is_contained = [](const AlnCoords &ctg1, const AlnCoords &ctg2) -> bool {
    if (ctg1.start >= ctg2.start && ctg1.stop <= ctg2.stop) return true;
    else return false;
  };
  
  if (is_contained(ctg1, ctg2) || is_contained(ctg2, ctg1)) {
    stats.containments++;
    return false;
  }
  int end1, end2;
  int gap;
  int gap_start;
  if (ctg1.start <= ctg2.start) {
    end1 = aln1->orient == '+' ? 3 : 5;
    end2 = aln2->orient == '+' ? 5 : 3;
    gap = ctg2.start - ctg1.stop;
    gap_start = ctg1.stop;
  } else {
    end1 = aln1->orient == '+' ? 5 : 3;
    end2 = aln2->orient == '+' ? 3 : 5;
    gap = ctg1.start - ctg2.stop;
    gap_start = ctg2.stop;
  }
  int min_aln_len = min(aln1->rstop - aln1->rstart, aln2->rstop - aln2->rstart);
  int min_aln_score = min(aln1->score1, aln2->score1);
  if (gap < -min_aln_len) {
    stats.short_alns++;
    return false;
  }
  if (gap < -aln1->rlen) DIE("Gap is too small: ", gap, ", read length ", aln1->rlen, "\n");
  // check for bad overlaps
  if (gap < 0 && (aln1->clen < -gap || aln2->clen < -gap)) {
    stats.bad_overlaps++;
    return false;
  }
  char orient1 = aln1->orient;
  CidPair cids = { .cid1 = aln1->cid, .cid2 = aln2->cid };
  if (aln1->cid < aln2->cid) {
    swap(end1, end2);
    swap(cids.cid1, cids.cid2);
    orient1 = aln2->orient;
  }
  Edge edge = { .cids = cids, .end1 = end1, .end2 = end2, .gap = gap, .support = 1, .aln_len = min_aln_len,
                .aln_score = min_aln_score, .edge_type = EdgeType::SPLINT, .seq = "",
                .mismatch_error = false, .conflict_error = false, .excess_error = false, .short_aln = false, .gap_reads = {}};
  if (edge.gap > 0) {
    edge.gap_reads = vector<GapRead>{GapRead(aln1->read_id, gap_start, -1, -1, orient1, cids.cid1)};
    _graph->add_pos_gap_read(aln1->read_id);
  }
  _graph->add_or_update_edge(edge);
  return true;
}


static void set_nbs(AlnStats &stats)
{
  Timer timer(__func__);
  {
    int64_t clen_excess = 0;
    int64_t num_excess_ctgs = 0;
    int max_excess_degree = 0;
    ProgressBar progbar(_graph->get_local_num_edges(), "Updating edges");
    // add edges to vertices
    for (auto edge = _graph->get_first_local_edge(); edge != nullptr; edge = _graph->get_next_local_edge()) {
      progbar.update();
      for (cid_t cid : {edge->cids.cid1, edge->cids.cid2}) {
        auto v = _graph->get_vertex(cid);
        assert(!edge->excess_error);
        if (v->end5.size() + v->end3.size() > MAX_CTG_GRAPH_DEGREE) {
          edge->excess_error = true;
          clen_excess += v->clen;
          num_excess_ctgs++;
          max_excess_degree = std::max(max_excess_degree, v->clen);
          break;
        }
      }
      if (edge->excess_error) continue;
      // minor race condition here but it shouldn't lead to very high degrees
      _graph->add_vertex_nb(edge->cids.cid1, edge->cids.cid2, edge->end1);
      _graph->add_vertex_nb(edge->cids.cid2, edge->cids.cid1, edge->end2);
    }
    progbar.done();
    barrier();
    auto tot_clen_excess = reduce_one(clen_excess, op_fast_add, 0).wait();
    auto tot_num_excess_ctgs = reduce_one(num_excess_ctgs, op_fast_add, 0).wait();
    int all_max_excess_clen = reduce_one(max_excess_degree, op_fast_max, 0).wait();
    if (tot_num_excess_ctgs) 
      SLOG_VERBOSE("Average excess clen ", ((double)tot_clen_excess / tot_num_excess_ctgs), " max ", all_max_excess_clen, "\n");
    int64_t num_excess_edges = reduce_one(_graph->purge_excess_edges(), op_fast_add, 0).wait();
    if (num_excess_edges) {
      if (rank_me() == 0 && tot_num_excess_ctgs == 0) SDIE("We have no excess ctgs, but ", num_excess_edges, " excess edges");
      SLOG_VERBOSE("Purged ", num_excess_edges, " excess degree edges\n");
    }
  }
  barrier();
  SLOG_VERBOSE(stats.to_string());
}


void add_pos_gap_read(Edge* edge, Aln &aln) {
  int end = (aln.cid == edge->cids.cid1 ? edge->end1 : edge->end2);
  int gap_start;
  if (aln.cid == edge->cids.cid1) {
    if ((edge->end1 == 3 && aln.orient == '+') || (edge->end1 == 5 && aln.orient == '-')) gap_start = aln.rstop;
    else gap_start = aln.rlen - aln.rstart;
    if (gap_start == aln.rlen) return;
  } else if (aln.cid == edge->cids.cid2) {
    // head
    if ((edge->end2 == 5 && aln.orient == '+') || (edge->end2 == 3 && aln.orient == '-')) gap_start = aln.rstart;
    else gap_start = aln.rlen - aln.rstop;
    if (gap_start == 0) return;
  } else {
    DIE("cid doesn't match in pos edge\n");
  }
  edge->gap_reads.push_back(GapRead(aln.read_id, gap_start, aln.rstart, aln.rstop, aln.orient, aln.cid));
  _graph->add_pos_gap_read(aln.read_id);
}


static AlnStats get_splints_from_alns(Alns &alns) {
  Timer timer(__func__);
  AlnStats stats = {0};
  IntermittentTimer t_get_alns("get alns splints");
  ProgressBar progbar(alns.size(), "Adding edges to graph from splints");
  int64_t aln_i = 0;
  int64_t num_splints = 0;
  while (true) {
    vector<Aln> alns_for_read;
    t_get_alns.start();
    if (!get_alns_for_read(alns, aln_i, alns_for_read, &stats.nalns, &stats.unaligned)) break;
    t_get_alns.stop();
    progbar.update(aln_i);
    for (int i = 0; i < alns_for_read.size(); i++) {
      auto aln = &alns_for_read[i];
      for (int j = i + 1; j < alns_for_read.size(); j++) {
        progress();
        auto other_aln = &alns_for_read[j];
        if (other_aln->read_id != aln->read_id) DIE("Mismatched read ids: ", other_aln->read_id, " != ", aln->read_id, "\n");
        if (add_splint(other_aln, aln, stats)) num_splints++;
      }
    }
  }
  progbar.done();
  barrier();
  auto tot_nalns = reduce_one(stats.nalns, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed ", tot_nalns, " alignments and found ",
               perc_str(reduce_one(num_splints, op_fast_add, 0).wait(), tot_nalns), " splints\n");
  return stats;
}


string get_consensus_seq(const vector<string> &seqs, int max_len)
{
  static char bases[5] = {'A', 'C', 'G', 'T', 'N'};
  auto base_freqs = new int[max_len][4]();
  for (auto seq : seqs) {
    for (int i = 0; i < seq.size(); i++) {
      switch (seq[i]) {
        case 'A': base_freqs[i][0]++; break;
        case 'C': base_freqs[i][1]++; break;
        case 'G': base_freqs[i][2]++; break;
        case 'T': base_freqs[i][3]++; break;
        case 'N': break;
        default:
          WARN("Unknown base at pos ", i, " (", seq.size(), "): ", seq[i], "\n", seq);
      }
    }
  }
  string consensus_seq = "";
  for (int i = 0; i < max_len; i++) {
    int max_freq = 0;
    // start off with N
    int highest_idx = 4;
    for (int j = 0; j < 4; j++) {
      if (base_freqs[i][j] > max_freq) {
        max_freq = base_freqs[i][j];
        highest_idx = j;
      }
    }
    consensus_seq += bases[highest_idx];
  }
  delete[] base_freqs;
  // trim any Ns off the front
  consensus_seq.erase(0, consensus_seq.find_first_not_of('N'));
  return consensus_seq;
}

  
static string get_splint_edge_seq(int kmer_len, Edge *edge)
{
  vector<string> seqs;
  // tail and end for checking primer matches
  int gap_size = edge->gap + 2 * (kmer_len - 1);
  for (auto gap_read : edge->gap_reads) {
    auto seq = _graph->get_read_seq(gap_read.read_name);
    if (seq == "") DIE("Could not find read seq for read ", gap_read.read_name, "\n");
    if (gap_read.gap_start < kmer_len) {
      //WARN("Positive gap overlap is less than kmer length, ", gap_read.gap_start, " < ", kmer_len, "\n");
      continue;
    }
    int rstart = gap_read.gap_start - (kmer_len - 1);
    string gap_seq = seq.substr(rstart, gap_size);
    if (edge->gap_reads.size() > 1) {
      string gap_seq_rc = revcomp(gap_seq);
      // these sequences should all be similar, so choosing the lexicographically least should ensure they have the same orientation
      if (gap_seq > gap_seq_rc) gap_seq = gap_seq_rc;
    }
    seqs.push_back(gap_seq);
  }
  return get_consensus_seq(seqs, gap_size);
}


static string get_span_edge_seq(int kmer_len, Edge *edge, bool tail)
{
  string ctg_seq = "";
  cid_t cid = (tail ? edge->cids.cid1 : edge->cids.cid2);
  vector<string> seqs;
  char buf[100];
  int max_len = 0;
  for (auto gap_read : edge->gap_reads) {
    if (gap_read.cid != cid) continue;
    sprintf(buf, "gs %3d rs %3d rp %3d %c ", gap_read.gap_start, gap_read.rstart, gap_read.rstop, gap_read.orient);
    auto gap_seq = _graph->get_read_seq(gap_read.read_name);
    if (gap_seq == "") DIE("Could not find read seq for read ", gap_read.read_name, "\n");
    if (tail) {
      if ((edge->end1 == 5 && gap_read.orient == '+') || (edge->end1 == 3 && gap_read.orient == '-')) gap_seq = revcomp(gap_seq);
      if (gap_read.gap_start > kmer_len) {
        gap_seq.erase(0, gap_read.gap_start - kmer_len);
        gap_read.gap_start = kmer_len;
      }
      DBG_SPANS(buf, gap_seq, "\n");
    } else {
      if ((edge->end2 == 3 && gap_read.orient == '+') || (edge->end2 == 5 && gap_read.orient == '-')) gap_seq = revcomp(gap_seq);
      if (gap_read.gap_start + kmer_len < gap_seq.size()) gap_seq.erase(gap_read.gap_start + kmer_len);
      // pad the front of the gap sequence with Ns to make them all the same length
      string offset_padding(1 + _graph->max_read_len - kmer_len - gap_read.gap_start, 'N');
      gap_seq = offset_padding + gap_seq;
      DBG_SPANS(buf, gap_seq, "\n");
    }
    seqs.push_back(gap_seq);
    if (gap_seq.size() > max_len) max_len = gap_seq.size();
  }
  if (seqs.empty()) {
    if (tail) {
      auto vertex = _graph->get_vertex(edge->cids.cid1);
      ctg_seq = _graph->get_vertex_seq(vertex->seq_gptr, vertex->clen);
      if (edge->end1 == 5) ctg_seq = revcomp(ctg_seq);
      int tail_len = vertex->clen - kmer_len;
      if (tail_len < 0) tail_len = 0;
      ctg_seq = ctg_seq.substr(tail_len);
      DBG_SPANS("TAIL contig", vertex->cid, ".", edge->end1, "\t", ctg_seq, "\n");
    } else {
      auto vertex = _graph->get_vertex(edge->cids.cid2);
      ctg_seq = _graph->get_vertex_seq(vertex->seq_gptr, vertex->clen);
      if (edge->end2 == 3) ctg_seq = revcomp(ctg_seq);
      ctg_seq = ctg_seq.substr(0, kmer_len);
      DBG_SPANS("HEAD contig", vertex->cid, ".", edge->end2, "\t", ctg_seq, "\n");
    }
    return ctg_seq;
  }
  return get_consensus_seq(seqs, max_len);
}


static void parse_reads(int kmer_len, const vector<string> &reads_fname_list) {
  Timer timer(__func__);

  int64_t num_seqs_added = 0;
  int64_t num_reads = 0;
  int max_read_len = 0;
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    FastqReader fqr(merged_reads_fname, PER_RANK_FILE);
    string id, seq, quals;
    ProgressBar progbar(fqr.my_file_size(), "Parsing reads for gap sequences");
    size_t tot_bytes_read = 0;
    while (true) {
      progress();
      size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
      if (!bytes_read) break;
      tot_bytes_read += bytes_read;
      progbar.update(tot_bytes_read);
      // this happens when we have a placeholder entry because reads merged
      if (kmer_len > seq.length()) continue;
      if (_graph->update_read_seq(id, seq)) num_seqs_added++;
      if (seq.length() > max_read_len) max_read_len = seq.length();
      num_reads++;
    }
    progbar.done();
    barrier();
  }
  _graph->max_read_len = reduce_all(max_read_len, op_fast_max).wait();  
  auto tot_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed a total of ", tot_num_reads, " reads, found max read length ", _graph->max_read_len, "\n");
  SLOG_VERBOSE("Extracted ", perc_str(reduce_one(num_seqs_added, op_fast_add, 0).wait(), tot_num_reads),
               " read sequences for pos gaps\n");

  int64_t num_pos_gaps = 0;
  int64_t num_pos_spans = 0;
  int64_t num_pos_spans_closed = 0;
  int64_t num_pos_spans_w_ns = 0;
  {
    ProgressBar progbar(_graph->get_local_num_edges(), "Fill pos gaps");
    for (auto edge = _graph->get_first_local_edge(); edge != nullptr; edge = _graph->get_next_local_edge()) {
      progbar.update();
      if (edge->gap <= 0) continue;
      num_pos_gaps++;
      if (edge->edge_type == EdgeType::SPAN) {
        num_pos_spans++;
        DBG_SPANS("SPAN pos gap ", edge->gap, " with ", edge->gap_reads.size(), " reads\n");
        string tail_seq = get_span_edge_seq(kmer_len, edge, true);
        if (tail_seq == "") continue;
        string head_seq = get_span_edge_seq(kmer_len, edge, false);
        if (head_seq == "") continue;
        DBG_SPANS("tail consensus         ", tail_seq, "\n");
        int offset = kmer_len + edge->gap - head_seq.size() + kmer_len;
        if (offset < 1) offset = 1;
        DBG_SPANS("head consensus         ", string(offset, ' '), head_seq, "\n");
        // now merge tail_seq and head_seq using the best (lowest hamming dist) overlap
        int min_len = min(tail_seq.size(), head_seq.size());
        int max_len = max(tail_seq.size(), head_seq.size());
        auto min_dist = min_hamming_dist(tail_seq, head_seq, min_len);
        if (is_overlap_mismatch(min_dist.first, min_dist.second)) {
          min_dist.first = min_len;
          min_dist.second = -1;
          for (int i = 0; i < max_len - min_len; i++) {
            int dist = hamming_dist(tail_seq.substr(0, min_len), head_seq.substr(0, min_len));
            if (dist < min_dist.first) {
              min_dist.first = dist;
              min_dist.second = i;
            }
          }
          if (is_overlap_mismatch(min_dist.first, min_dist.second)) {
            DBG_SPANS("overlap mismatch: hdist ", min_dist.first, " best overlap ", min_dist.second, " original gap ", edge->gap, "\n");
            if (tail_seq.size() + head_seq.size() < 2 * kmer_len + edge->gap) {
              // the gap is not closed - fill with Ns
              int num_ns = 2 * kmer_len + edge->gap - tail_seq.size() - head_seq.size();
              edge->seq = tail_seq + string(num_ns, 'N') + head_seq;
              // remove one of either end because the later checking will use kmer_len - 1
              edge->seq = edge->seq.substr(1, edge->seq.size() - 2);
              DBG_SPANS("using orig gap ", edge->gap, " with ", num_ns, " Ns: ", edge->seq, "\n");
              num_pos_spans_w_ns++;
              num_pos_spans_closed++;
            }
            continue;
          }
        }
        num_pos_spans_closed++;
        int gap_size = tail_seq.size() + head_seq.size() - min_dist.second - 2 * kmer_len;
        DBG_SPANS("overlap is ", min_dist.second, " original gap is ", edge->gap, " corrected gap is ", gap_size, "\n");
        edge->gap = gap_size;
        DBG_SPANS(tail_seq, "\n");
        DBG_SPANS(string(tail_seq.size() - min_dist.second, ' '), head_seq, "\n");
        tail_seq.erase(tail_seq.size() - min_dist.second);
        edge->seq = tail_seq + head_seq;
        edge->seq = edge->seq.substr(1, edge->seq.size() - 2);
        DBG_SPANS(edge->seq, "\n");
        if (edge->seq.size() != 2 * (kmer_len - 1) + gap_size)
          WARN("fill mismatch ", edge->seq.size(), " != ", 2 * kmer_len + gap_size);
      } else {
        edge->seq = get_splint_edge_seq(kmer_len, edge);
      }
    }
    progbar.done();
    barrier();
  }
  auto tot_pos_gaps = reduce_one(num_pos_gaps, op_fast_add, 0).wait();
  SLOG_VERBOSE("Filled ", tot_pos_gaps, " positive gaps with ", _graph->get_num_read_seqs(), " reads\n");
  auto tot_pos_spans = reduce_one(num_pos_spans, op_fast_add, 0).wait();
  auto tot_pos_spans_closed = reduce_one(num_pos_spans_closed, op_fast_add, 0).wait();
  auto tot_pos_spans_w_ns = reduce_one(num_pos_spans_w_ns, op_fast_add, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(tot_pos_spans, tot_pos_gaps), " positive spans, of which ",
               perc_str(tot_pos_spans_closed, tot_pos_spans), " were closed - ",
               perc_str(tot_pos_spans_w_ns, tot_pos_spans_closed), " with Ns\n");
}


static bool merge_end(Vertex *curr_v, const vector<cid_t> &nb_cids, vector<vector<cid_t> > &nb_cids_merged, 
               IntermittentTimer &t_merge_get_nbs, IntermittentTimer &t_merge_sort_nbs, IntermittentTimer &t_merge_output_nbs) 
{
  nb_cids_merged.clear();

  // dead-end
  if (!nb_cids.size()) {
    DBG_BUILD("\tdead end\n");
    progress();
    return false;
  }

  struct NbPair {
    shared_ptr<Vertex> vertex;
    shared_ptr<Edge> edge;
  };
  
  t_merge_get_nbs.start();
  vector<NbPair> nbs;
  // select the neighbors that are valid for connections
  for (auto nb_cid : nb_cids) {
    auto nb_edge = _graph->get_edge_cached(curr_v->cid, nb_cid);
    // drop low support for high degree nodes - this reduces computational overhead and load imbalance, which can get severe
    if (nb_cids.size() > 50 && nb_edge->support < 2) continue;
    if (nb_cids.size() > 100 && nb_edge->support < 3) continue;
    if (nb_cids.size() > 150 && nb_edge->support < 4) continue;
    nbs.push_back({_graph->get_vertex_cached(nb_cid), nb_edge});
  }
  t_merge_get_nbs.stop();
  if (nbs.empty()) return false;
  // just one nb found
  if (nbs.size() == 1) {
    nb_cids_merged.push_back({nbs[0].vertex->cid});
    DBG_BUILD("\t", nbs[0].vertex->cid, " gap ", nbs[0].edge->gap, " len ", nbs[0].vertex->clen, " depth ", nbs[0].vertex->depth, "\n");
    return false;
  }
  t_merge_sort_nbs.start();
  // found multiple nbs, check for overlaps that can be merged
  // first, sort nbs by gap size
  sort(nbs.begin(), nbs.end(), 
       [](const auto &elem1, const auto &elem2) {
         return elem1.edge->gap < elem2.edge->gap;
       });
  t_merge_sort_nbs.stop();

  // gather a vector of merged paths (there can be more than one because of forks)
  vector<vector<NbPair*> > all_next_nbs = {};
  // attempt to merge all neighbors as overlaps 
  for (int i = 0; i < nbs.size(); i++) {
    NbPair *nb = &nbs[i];
    DBG_BUILD("\t", nb->vertex->cid, " gap ", nb->edge->gap, " len ", nb->vertex->clen, " depth ", nb->vertex->depth, "\n");
    if (i == 0) {
      all_next_nbs.push_back({nb});
      continue;
    }
    bool found_merge = false;
    for (auto &next_nbs : all_next_nbs) {
      progress();
      auto prev_nb = next_nbs.back();
      int g1 = -prev_nb->edge->gap;
      int g2 = -nb->edge->gap;
      int gdiff = g1 - g2;
      // check gap spacing to see if current nb overlaps previous nb
      if ((prev_nb->edge->gap == nb->edge->gap) || (gdiff >= prev_nb->vertex->clen - 10) || 
          (gdiff + nb->vertex->clen <= prev_nb->vertex->clen)) {
        DBG_BUILD("\tMerge conflict ", prev_nb->vertex->cid, " ", nb->vertex->cid, "\n");
        continue;
      }
      // now check that there is an edge from this vertex to the previous one
      auto intermediate_edge = _graph->get_edge_cached(nb->vertex->cid, prev_nb->vertex->cid);
      if (!intermediate_edge) {
        DBG_BUILD("\tNo edge found between ", prev_nb->vertex->cid, " and ", nb->vertex->cid, "\n");
        continue;
      } 
      // now check the overlaps are correct
      DBG_BUILD("\tMERGE ", prev_nb->vertex->cid, " ", nb->vertex->cid, "\n");
      next_nbs.push_back(nb);
      found_merge = true;
      break;
    }
    if (!found_merge) all_next_nbs.push_back({nb});
  }
  if (all_next_nbs.empty()) DIE("empty all_next_nbs. How?\n");

  t_merge_output_nbs.start();
  // update the merged nbs and write out the merged paths found
  for (auto &next_nbs : all_next_nbs) {
    progress();
    nb_cids_merged.push_back({});
    auto last_nb = next_nbs.back();
    DBG_BUILD("\tmerged path (len ", last_nb->vertex->clen + last_nb->edge->gap, ", depth ", last_nb->vertex->depth, "): ");
    for (auto &next_nb : next_nbs) {
      DBG_BUILD(next_nb->vertex->cid, " ");
      nb_cids_merged.back().push_back(next_nb->vertex->cid);
    }
    DBG_BUILD("\n");
  }    
  t_merge_output_nbs.stop();

  return true;
}

static void merge_nbs()
{
  barrier();
  Timer timer(__func__);
  int64_t num_merges = 0;
  int64_t num_orphans = 0;
  int max_orphan_len = 0;
  int64_t num_nbs = 0, num_merged_nbs = 0, max_nbs = 0;
  double max_orphan_depth = 0;
  {
    IntermittentTimer t_merge_ends("merge ends"), t_merge_get_nbs("merge get nbs"), t_merge_sort_nbs("merge sort nbs"), 
      t_merge_output_nbs("merge output nbs");
    ProgressBar progbar(_graph->get_local_num_vertices(), "Merge nbs");
    // mark all the vertices that have forks and the side of the forks. Note that in many cases what look like forks are
    // actually vertices that should be merged into a single neighbor
    for (auto v = _graph->get_first_local_vertex(); v != nullptr; v = _graph->get_next_local_vertex()) {
      DBG_BUILD(v->cid,  " len ",  v->clen,  " depth ",  v->depth,  ":\n");
      if (v->end5.empty() && v->end3.empty()) {
        num_orphans++;
        if (v->clen > max_orphan_len) max_orphan_len = v->clen;
        if (v->depth > max_orphan_depth) max_orphan_depth = v->depth;
        DBG_BUILD("\torphan: ", v->cid, " len ", v->clen, " depth ", v->depth, "\n");
      }
      t_merge_ends.start();
      DBG_BUILD("  5-end:\n");
      if (merge_end(v, v->end5, v->end5_merged, t_merge_get_nbs, t_merge_sort_nbs, t_merge_output_nbs)) num_merges++;
      DBG_BUILD("  3-end:\n");
      if (merge_end(v, v->end3, v->end3_merged, t_merge_get_nbs, t_merge_sort_nbs, t_merge_output_nbs)) num_merges++;
      t_merge_ends.stop();
      int v_num_nbs = v->end5.size() + v->end3.size();
      num_nbs += v_num_nbs;
      num_merged_nbs += v->end5_merged.size() + v->end3_merged.size();
      if (v_num_nbs > max_nbs) max_nbs = v_num_nbs;
      progbar.update();
    }
    DBG("Number of nbs ", num_nbs, " avg degree ", ((double)num_nbs / _graph->get_local_num_vertices()), 
        " merged ", num_merged_nbs, " avg degree ", ((double)num_merged_nbs / _graph->get_local_num_vertices()), 
        " max degree ", max_nbs, "\n");
    progbar.done();
  }
  barrier();
  auto tot_merges = reduce_one(num_merges, op_fast_add, 0).wait();
  SLOG_VERBOSE("Merged ", perc_str(tot_merges, 2 *_graph->get_num_vertices()), " vertices\n");
  auto tot_orphans = reduce_one(num_orphans, op_fast_add, 0).wait();
  auto all_max_orphan_len = reduce_one(max_orphan_len, op_fast_max, 0).wait();
  auto all_max_orphan_depth = reduce_one(max_orphan_depth, op_fast_max, 0).wait();
  SLOG_VERBOSE("Found ", perc_str(tot_orphans, _graph->get_num_vertices()), " orphaned vertices (no edges), max length ", 
               all_max_orphan_len, ", max depth ", all_max_orphan_depth, "\n");
}


void run_spanner(int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len, Alns &alns, CtgGraph *graph);
  
void mark_short_aln_edges(int max_kmer_len) {
  Timer timer(__func__);
  // make sure we don't use an out-of-date edge
  barrier();
  _graph->clear_caches();
  {
    ProgressBar progbar(_graph->get_local_num_vertices(), "Mark short aln edges");
    for (auto v = _graph->get_first_local_vertex(); v != nullptr; v = _graph->get_next_local_vertex()) {
      for (auto cid_list : { v->end5, v->end3 }) {
        vector<CidPair> drop_edges = {};
        bool long_aln_found = false;
        for (auto i = 0; i < cid_list.size(); i++) {
          auto edge = _graph->get_edge(v->cid, cid_list[i]);
          if (edge->aln_len >= max_kmer_len) long_aln_found = true;
          else drop_edges.push_back(edge->cids);
        }
        if (long_aln_found) {
          for (auto cids : drop_edges) _graph->mark_edge_short_aln(cids);
        }
      }
      progbar.update();
    }
    progbar.done();
  }
  barrier();
  int64_t num_edges = _graph->get_num_edges();
  int64_t num_short = _graph->purge_short_aln_edges();
  SLOG_VERBOSE("Purged ", perc_str(reduce_one(num_short, op_fast_add, 0).wait(), num_edges), " short aln edges\n");
}

void build_ctg_graph(CtgGraph *graph, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len,
                     vector<string> &reads_fname_list, Contigs *ctgs, Alns &alns) {
  Timer timer(__func__);
  _graph = graph;
  add_vertices_from_ctgs(ctgs);
  auto aln_stats = get_splints_from_alns(alns);
  get_spans_from_alns(insert_avg, insert_stddev, max_kmer_len, kmer_len, alns, graph);
  _graph->purge_error_edges(&aln_stats.mismatched, &aln_stats.conflicts, &aln_stats.empty_spans);
  barrier();
  set_nbs(aln_stats);
  //mark_short_aln_edges(max_kmer_len);
  parse_reads(kmer_len, reads_fname_list);
  merge_nbs();
}
