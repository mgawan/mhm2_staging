#include "ssw.hpp"
#include "gtest/gtest.h"
#include "klign.hpp"

#include <sstream>
#include <string>
#include <vector>

using std::vector;
using std::string;

using namespace StripedSmithWaterman;

string aln2string(Alignment &aln) {
  std::stringstream ss;
  ss << "score=" << aln.sw_score << " score2=" << aln.sw_score_next_best;
  ss << " rbegin=" << aln.ref_begin << " rend=" << aln.ref_end;
  ss << " qbegin=" << aln.query_begin << " qend=" << aln.query_end;
  ss << " rend2=" << aln.ref_end_next_best << " mismatches=" << aln.mismatches;
  ss << " cigarstr=" << aln.cigar_string;
  return ss.str();
}

AlnScoring aln_scoring = {.match = ALN_MATCH_SCORE,
                          .mismatch = ALN_MISMATCH_COST,
                          .gap_opening = ALN_GAP_OPENING_COST,
                          .gap_extending = ALN_GAP_EXTENDING_COST,
                          .ambiguity = ALN_AMBIGUITY_COST};
AlnScoring cigar_aln_scoring = {.match = 2, .mismatch = 4, .gap_opening = 4, .gap_extending = 2, .ambiguity = 1};

Aligner ssw_aligner;
Aligner ssw_aligner_mhm2(aln_scoring.match, aln_scoring.mismatch, aln_scoring.gap_opening, aln_scoring.gap_extending,
                         aln_scoring.ambiguity);
Aligner ssw_aligner_cigar(cigar_aln_scoring.match, cigar_aln_scoring.mismatch, cigar_aln_scoring.gap_opening,
                          cigar_aln_scoring.gap_extending, cigar_aln_scoring.ambiguity);

Filter ssw_filter(true, false, 0, 32767), ssw_filter_cigar(true, true, 0, 32767);

void test_aligns(vector<Alignment> &alns, string query, string ref) {
  alns.resize(6);
  auto reflen = ref.size();
  auto qlen = query.size();
  auto masklen = max((int)min(reflen, qlen) / 2, 15);
  ssw_aligner.Align(query.c_str(), ref.c_str(), reflen, ssw_filter, &alns[0], masklen);
  ssw_aligner.Align(query.c_str(), ref.c_str(), reflen, ssw_filter_cigar, &alns[1], masklen);

  ssw_aligner_mhm2.Align(query.c_str(), ref.c_str(), reflen, ssw_filter, &alns[2], masklen);
  ssw_aligner_mhm2.Align(query.c_str(), ref.c_str(), reflen, ssw_filter_cigar, &alns[3], masklen);

  ssw_aligner_cigar.Align(query.c_str(), ref.c_str(), reflen, ssw_filter, &alns[4], masklen);
  ssw_aligner_cigar.Align(query.c_str(), ref.c_str(), reflen, ssw_filter_cigar, &alns[5], masklen);
}

void check_alns(vector<Alignment> &alns, int qstart, int qend, int rstart, int rend, int mismatches, string query = "",
                string ref = "", string cigar = "") {
  for (Alignment &aln : alns) {
    EXPECT_EQ(aln.ref_begin, rstart) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.ref_end, rend) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.query_begin, qstart) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.query_end, qend) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    if (!aln.cigar_string.empty()) {  // mismatches should be recorded...
      EXPECT_EQ(aln.mismatches, mismatches) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
      if (!cigar.empty()) EXPECT_STREQ(aln.cigar_string.c_str(), cigar.c_str());
    }
  }
}
void check_not_alns(vector<Alignment> &alns, string query = "", string ref = "") {
  for (Alignment &aln : alns) {
    EXPECT_EQ(aln.ref_begin, 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.ref_end, 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.query_begin, 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.query_end, 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.mismatches, 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.sw_score, 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.sw_score_next_best, 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
  }
}

TEST(MHMTest, ssw) {
  // arrange
  // act
  // assert

  EXPECT_EQ(ssw_filter.report_cigar, false);
  EXPECT_EQ(ssw_filter_cigar.report_cigar, true);

  vector<Alignment> alns;
  string ref = "ACGT";
  string query = ref;
  test_aligns(alns, query, ref);
  check_alns(alns, 0, 3, 0, 3, 0, query, ref);
  ref = "AACGT";
  test_aligns(alns, query, ref);
  check_alns(alns, 0, 3, 1, 4, 0, query, ref);
  query = "TACGT";
  test_aligns(alns, query, ref);
  check_alns(alns, 1, 4, 1, 4, 0, query, ref);
  ref = "ACGT";
  query = "TACGT";
  test_aligns(alns, query, ref);
  check_alns(alns, 1, 4, 0, 3, 0, query, ref);

  string r = "AAAATTTTCCCCGGGG";
  string q = "AAAATTTTCCCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 15, 0, 15, 0, q, r);

  // 1 subst
  q = "AAAATTTTACCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 15, 0, 15, 1, q, r);

  // 1 insert
  q = "AAAATTTTACCCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 16, 0, 15, 1, q, r);

  // 1 del
  q = "AAAATTTCCCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 14, 0, 15, 1, q, r);

  // no match
  q = "ACGTCGTAGTACTACGT";
  test_aligns(alns, q, r);
  check_not_alns(alns, q, r);
}