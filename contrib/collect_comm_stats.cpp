#include <iostream>
#include <fstream>
#include <regex>

#include "upcxx/upcxx.hpp"

using namespace upcxx;
using namespace std;

vector<string> split_string(const string &content, string in_pattern) {
  vector<string> split_content;
  std::regex pattern(in_pattern);
  copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1), std::sregex_token_iterator(),
       back_inserter(split_content));
  return split_content;
}


// 0> Total gets:                  5758  avg/min/max/total sz = 420426.398/584/678504/2420815200
// 0>  AMREQUEST_MEDIUM:         492466  avg/min/max/total sz = 111.411/0/4032/54866232
// 0>  AMREQUEST_LONG:              192  avg/min/max/total sz = 17407.500/8/131072/3342240

void get_sizes(const string &s, double &avg_size, long &min_size, long &max_size, long &tot_size) {
  auto sizes = split_string(s, "/");
  avg_size = stod(sizes[0]);
  min_size = stol(sizes[1]);
  max_size = stol(sizes[2]);
  tot_size = stol(sizes[3]);
}

tuple<long, long, long, long, double> reduce_sizes(long num, long min_size, long max_size, long tot) {
  auto all_num = reduce_one(num, op_fast_add, 0).wait() / rank_n();
  auto all_min_size = reduce_one(min_size, op_fast_min, 0).wait();
  auto all_max_size = reduce_one(max_size, op_fast_max, 0).wait();
  auto all_tot_size = reduce_one(tot, op_fast_add, 0).wait() / rank_n();
  auto max_tot_size = reduce_one(tot, op_fast_max, 0).wait();
  return {all_num, all_min_size, all_max_size, all_tot_size, (double)max_tot_size / (all_tot_size * max_tot_size)};
}

int main(int argc, char **argv) {
  upcxx::init();
  string my_fname = argv[1] + to_string(rank_me());
  ifstream stats_file(my_fname);
  string line;
  long num_gets = 0;
  double avg_gets_size = 0;
  long min_gets_size = 0, max_gets_size = 0, tot_gets_size = 0;
  long num_am_medium = 0;
  double avg_am_medium_size = 0;
  long min_am_medium_size = 0, max_am_medium_size = 0, tot_am_medium_size = 0;
  long num_am_long = 0;
  double avg_am_long_size = 0;
  long min_am_long_size = 0, max_am_long_size = 0, tot_am_long_size = 0;
  
  while (getline(stats_file, line)) {
    auto tokens = split_string(line, "\\s+");
    if (tokens.size() < 3) continue;
    if (tokens[1] == "Total" && tokens[2] == "gets:") {
      num_gets = stol(tokens[3]);
      if (num_gets > 0) get_sizes(tokens[7], avg_gets_size, min_gets_size, max_gets_size, tot_gets_size);
    } else if (tokens[1] == "AMREQUEST_MEDIUM:") {
      num_am_medium = stol(tokens[2]);
      if (num_am_medium > 0) get_sizes(tokens[6], avg_am_medium_size, min_am_medium_size, max_am_medium_size, tot_am_medium_size);
    } else if (tokens[1] == "AMREQUEST_LONG:") {
      num_am_long = stol(tokens[2]);
      if (num_am_long > 0) get_sizes(tokens[6], avg_am_long_size, min_am_long_size, max_am_long_size, tot_am_long_size);
    }
  }
  barrier();
  auto [all_num_gets, all_min_gets_size, all_max_gets_size, all_tot_gets_size, gets_balance] =
    reduce_sizes(num_gets, min_gets_size, max_gets_size, tot_gets_size);
  auto [all_num_am_medium, all_min_am_medium_size, all_max_am_medium_size, all_tot_am_medium_size, am_medium_balance] =
    reduce_sizes(num_am_medium, min_am_medium_size, max_am_medium_size, tot_am_medium_size);
  auto [all_num_am_long, all_min_am_long_size, all_max_am_long_size, all_tot_am_long_size, am_long_balance] =
    reduce_sizes(num_am_long, min_am_long_size, max_am_long_size, tot_am_long_size);
  if (!rank_me()) {
    cout << "Stats are (per rank): count   avg_sz   min_sz   max_sz   tot_sz   balance\n";
    cout << "gets " << all_num_gets << " " << (double)all_tot_gets_size / all_num_gets << " " << all_min_gets_size
         << " " << all_max_gets_size  << " " << all_tot_gets_size << " " << gets_balance << endl;
    cout << "am_medium " << all_num_am_medium << " " << (double)all_tot_am_medium_size / all_num_am_medium << " "
         << all_min_am_medium_size << " " << all_max_am_medium_size  << " " << all_tot_am_medium_size
         << " " << am_medium_balance << endl;
    cout << "am_long " << all_num_am_long << " " << (double)all_tot_am_long_size / all_num_am_long << " " << all_min_am_long_size
         << " " << all_max_am_long_size  << " " << all_tot_am_long_size
         << " " << am_long_balance << endl;
  }
  stats_file.close();
  upcxx::finalize();
  return 0;
}
