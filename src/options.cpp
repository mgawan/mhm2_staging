/*
 HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from the U.S. Dept. of Energy).  All rights reserved."

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 (1) Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 (2) Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 (3) Neither the name of the University of California, Lawrence Berkeley National
 Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 endorse or promote products derived from this software without specific prior
 written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 to the features, functionality or performance of the source code ("Enhancements") to
 anyone; however, if you choose to make your Enhancements available either publicly,
 or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 written license agreement for such Enhancements, then you hereby grant the following
 license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 prepare derivative works, incorporate into other computer software, distribute, and
 sublicense such enhancements or derivative works thereof, in binary and source code
 form.
*/

#include "options.hpp"

#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <upcxx/upcxx.hpp>

#include "CLI11.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "version.h"

using namespace upcxx_utils;
using namespace std;

#define YES_NO(X) ((X) ? "YES" : "NO")

vector<string> Options::splitter(string in_pattern, string &content) {
  vector<string> split_content;
  std::regex pattern(in_pattern);
  copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1), std::sregex_token_iterator(),
       back_inserter(split_content));
  return split_content;
}

bool Options::extract_previous_lens(vector<unsigned> &lens, unsigned k) {
  for (size_t i = 0; i < lens.size(); i++) {
    if (lens[i] == k) {
      lens.erase(lens.begin(), lens.begin() + i + 1);
      return true;
    }
  }
  return false;
}

bool Options::find_restart(string stage_type, int k) {
  string new_ctgs_fname(stage_type + "-" + to_string(k) + ".fasta");
  if (!file_exists(new_ctgs_fname)) return false;
  if (stage_type == "contigs") {
    if (k == kmer_lens.back() && stage_type == "contigs") {
      max_kmer_len = kmer_lens.back();
      kmer_lens = {};
      stage_type = "scaff-contigs";
      k = scaff_kmer_lens[0];
    } else {
      if (!extract_previous_lens(kmer_lens, k)) {
        SWARN("Cannot find kmer length ", k, " in configuration: ", vec_to_str(kmer_lens));
        return false;
      }
      prev_kmer_len = k;
    }
    ctgs_fname = new_ctgs_fname;
  } else if (stage_type == "uutigs") {
    if (!extract_previous_lens(kmer_lens, k)) {
      SWARN("Cannot find kmer length ", k, " in configuration: ", vec_to_str(kmer_lens));
      return false;
    }
    kmer_lens.insert(kmer_lens.begin(), k);
    prev_kmer_len = k;
    ctgs_fname = new_ctgs_fname;
  } else if (stage_type == "scaff-contigs") {
    max_kmer_len = kmer_lens.back();
    kmer_lens = {};
    if (k == scaff_kmer_lens.front()) {
      ctgs_fname = "contigs-" + to_string(k) + ".fasta";
    } else {
      if (k == scaff_kmer_lens.back()) k = scaff_kmer_lens[scaff_kmer_lens.size() - 2];
      ctgs_fname = "scaff-contigs-" + to_string(k) + ".fasta";
    }
    if (!extract_previous_lens(scaff_kmer_lens, k)) {
      SWARN("Cannot find kmer length ", k, " in configuration: ", vec_to_str(scaff_kmer_lens));
      return false;
    }
  } else {
    SWARN("Invalid previous stage '", stage_type, "' k ", k, ", could not restart");
    return false;
  }
  SLOG(KLBLUE, "Restart options:\n", "  kmer-lens =              ", vec_to_str(kmer_lens), "\n",
       "  scaff-kmer-lens =        ", vec_to_str(scaff_kmer_lens), "\n", "  prev-kmer-len =          ", prev_kmer_len, "\n",
       "  max-kmer-len =           ", max_kmer_len, "\n", "  contigs =                ", ctgs_fname, KNORM, "\n");
  if (!upcxx::rank_me() && !file_exists(ctgs_fname)) {
    SWARN("File ", ctgs_fname, " not found. Did the previous run have --checkpoint enabled?");
    return false;
  }
  return true;
}

void Options::get_restart_options() {
  // check directory for most recent contigs file dump
  bool found = false;
  for (auto it = scaff_kmer_lens.rbegin(); it != scaff_kmer_lens.rend(); ++it) {
    if ((found = find_restart("scaff-contigs", *it)) == true) break;
  }
  if (!found) {
    for (auto it = kmer_lens.rbegin(); it != kmer_lens.rend(); ++it) {
      if ((found = find_restart("contigs", *it)) == true) break;
      if ((found = find_restart("uutigs", *it)) == true) break;
    }
  }
  if (!found) {
    SWARN("No previously completed stage found, restarting from the beginning\n");
  }
}

void Options::setup_output_dir() {
  if (!upcxx::rank_me()) {
    // create the output directory and stripe it if not doing a restart
    if (restart) {
      if (access(output_dir.c_str(), F_OK) == -1) {
        SDIE("Output directory ", output_dir, " for restart does not exist");
      }
    } else {
      if (mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */) == -1) {
        // could not create the directory
        if (errno == EEXIST) {
          cerr << KLRED << "WARNING: " << KNORM << "Output directory " << output_dir
               << " already exists. May overwrite existing files\n";
        } else {
          ostringstream oss;
          oss << KLRED << "ERROR: " << KNORM << " Could not create output directory " << output_dir << ": " << strerror(errno)
              << endl;
          throw std::runtime_error(oss.str());
        }
      } else {
        // created the directory - now stripe it if possible
        auto status = std::system("which lfs 2>&1 > /dev/null");
        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
          string cmd = "lfs setstripe -c -1 " + output_dir;
          auto status = std::system(cmd.c_str());
          if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
            cout << "Set Lustre striping on the output directory\n";
          else
            cout << "Failed to set Lustre striping on output directory: " << WEXITSTATUS(status) << endl;

          // ensure per_thread dir exists and has stripe 1
          string per_thread = output_dir + "/per_thread";
          mkdir(per_thread.c_str(), S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);  // ignore any errors
          cmd = "lfs setstripe -c 1 " + per_thread;
          status = std::system(cmd.c_str());
          if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
            cout << "Set Lustre striping on the per_thread output directory\n";
          else
            cout << "Failed to set Lustre striping on per_thread output directory: " << WEXITSTATUS(status) << endl;
          // this should avoid contention on the filesystem when ranks start racing to creating these top levels
          for (int i = 0; i < rank_n(); i += 1000) {
            char basepath[256];
            sprintf(basepath, "%s/%08d", per_thread.c_str(), i);
            auto ret = mkdir(basepath, S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);
            if (ret != 0) break;  // ignore any errors, just stop
          }
        }
      }
    }
  }

  upcxx::barrier();
  // after we change to the output directory, relative paths could be incorrect, so make sure we have the correct path of the
  // reads files
  char cwd_str[FILENAME_MAX];
  if (!getcwd(cwd_str, FILENAME_MAX)) SDIE("Cannot get current working directory: ", strerror(errno));
  for (auto &fname : reads_fnames) {
    if (fname[0] != '/') fname = string(cwd_str) + "/" + fname;
  }
  // all change to the output directory
  if (chdir(output_dir.c_str()) == -1 && !upcxx::rank_me()) {
    DIE("Cannot change to output directory ", output_dir, ": ", strerror(errno));
  }
  upcxx::barrier();
}

void Options::setup_log_file() {
  if (!upcxx::rank_me()) {
    // check to see if mhm2.log exists. If so, and not restarting, rename it
    if (file_exists("mhm2.log") && !restart) {
      string new_log_fname = "mhm2-" + get_current_time(true) + ".log";
      cerr << KLRED << "WARNING: " << KNORM << output_dir << "/mhm2.log exists. Renaming to " << output_dir << "/" << new_log_fname
           << endl;
      if (rename("mhm2.log", new_log_fname.c_str()) == -1) DIE("Could not rename mhm2.log: ", strerror(errno));
    } else if (!file_exists("mhm2.log") && restart) {
      DIE("Could not restart - missing mhm2.log in tis directory");
    }
  }
  upcxx::barrier();
}

Options::Options() {
  output_dir = string("mhm2-run-<reads_fname[0]>-n") + to_string(upcxx::rank_n()) + "-N" +
               to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + get_current_time(true);
}
Options::~Options() {
  flush_logger();
  cleanup();
}

void Options::cleanup() {
  // cleanup and close loggers that Options opened in load
  close_logger();
#ifdef DEBUG
  close_dbg();
#endif
}

bool Options::load(int argc, char **argv) {
  // MHM2 version v0.1-a0decc6-master (Release) built on 2020-04-08T22:15:40 with g++
  string full_version_str = "MHM2 version " + string(MHM2_VERSION) + "-" + string(MHM2_BRANCH) + " with upcxx-utils " +
                            string(UPCXX_UTILS_VERSION) + " built on " + string(MHM2_BUILD_DATE);
  CLI::App app(full_version_str);
  // basic options - see user guide
  app.add_option("-r, --reads", reads_fnames, "Files containing merged and unmerged reads in FASTQ format (comma separated).")
      ->delimiter(',')
      ->check(CLI::ExistingFile);
  app.add_option("-p, --paired-reads", paired_fnames,
                 "Pairs of files containing paired reads for the same insert in FASTQ format (comma separated).")
      ->delimiter(',')
      ->check(CLI::ExistingFile);
  app.add_option("-i, --insert", insert_size, "Insert size (average:stddev) (autodetected by default).")
      ->delimiter(':')
      ->expected(2)
      ->check(CLI::Range(1, 10000));
  auto *kmer_lens_opt = app.add_option("-k, --kmer-lens", kmer_lens, "kmer lengths (comma separated) for contigging.")
                            ->delimiter(',')
                            ->capture_default_str();
  auto *scaff_kmer_lens_opt = app.add_option("-s, --scaff-kmer-lens", scaff_kmer_lens,
                                             "kmer lengths (comma separated) for scaffolding (set to 0 to disable scaffolding).")
                                  ->delimiter(',')
                                  ->capture_default_str();
  app.add_option("--min-ctg-print-len", min_ctg_print_len, "Minimum length required for printing a contig in the final assembly.")
      ->capture_default_str()
      ->check(CLI::Range(0, 100000));
  auto *output_dir_opt = app.add_option("-o,--output", output_dir, "Output directory.")->capture_default_str();
  app.add_flag("--checkpoint", checkpoint, "Enable checkpointing.")
      ->default_val(checkpoint ? "true" : "false")
      ->capture_default_str()
      ->multi_option_policy();
  app.add_flag("--restart", restart,
               "Restart in previous directory where a run failed (must specify the previous directory with -o).")
      ->capture_default_str();
  app.add_flag("--post-asm-align", post_assm_aln, "Align reads to final assembly")->capture_default_str();
  app.add_flag("--post-asm-abd", post_assm_abundances, "Compute and output abundances for final assembly (used by MetaBAT).")
      ->capture_default_str();
  app.add_flag("--post-asm-only", post_assm_only, "Only run post assembly (alignment and/or abundances).")->capture_default_str();
  app.add_flag("--write-gfa", dump_gfa, "Write scaffolding contig graphs in GFA2 format.")->capture_default_str();
  app.add_option("-Q, --quality-offset", qual_offset, "Phred encoding offset (auto-detected by default).")
      ->check(CLI::IsMember({0, 33, 64}));
  app.add_flag("--progress", show_progress, "Show progress bars for operations.");
  app.add_flag("-v, --verbose", verbose, "Verbose output: lots of detailed information (always available in the log).");
  auto *cfg_opt = app.set_config("--config", "", "Load options from a configuration file.");

  // advanced options
  // restarts
  app.add_option("-c, --contigs", ctgs_fname, "FASTA file containing contigs used for restart.");
  app.add_option("--max-kmer-len", max_kmer_len, "Maximum contigging kmer length for restart (only needed if not contigging).")
      ->check(CLI::Range(0, 159));
  app.add_option("--prev-kmer-len", prev_kmer_len,
                 "Previous contigging kmer length for restart (needed if contigging and contig file is specified).")
      ->check(CLI::Range(0, 159));
  // quality tuning
  app.add_option("--break-scaff-Ns", break_scaff_Ns, "Number of Ns allowed before a scaffold is broken.")
      ->check(CLI::Range(0, 1000));
  app.add_option("--min-depth-thres", dmin_thres, "Absolute mininimum depth threshold for DeBruijn graph traversal")
      ->check(CLI::Range(1, 100));
  // performance trade-offs
  app.add_option("--max-kmer-store", max_kmer_store_mb, "Maximum size for kmer store in MB per rank (set to 0 for auto 1% memory).")
      ->check(CLI::Range(0, 1000));
  app.add_option("--max-rpcs-in-flight", max_rpcs_in_flight,
                 "Maximum number of RPCs in flight, per process (set to 0 for unlimited).")
      ->check(CLI::Range(0, 10000));
  app.add_flag("--use-heavy-hitters", use_heavy_hitters, "Enable the Heavy Hitter Streaming Store (experimental).");
  app.add_option("--ranks-per-gpu", ranks_per_gpu, "Number of processes multiplexed to each GPU (default depends on hardware).")
      ->check(CLI::Range(0, (int)upcxx::local_team().rank_n() * 8));
  auto *bloom_opt = app.add_flag("--force-bloom", force_bloom, "Always use bloom filters.")
//      ->default_val(force_bloom ? "true" : "false")
      ->multi_option_policy();
  app.add_flag("--pin", pin_by,
               "Restrict processes according to logical CPUs, cores (groups of hardware threads), "
               "or NUMA domains (cpu, core, numa, none).")
      ->check(CLI::IsMember({"cpu", "core", "numa", "none"}));
  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    if (upcxx::rank_me() == 0) {
      if (e.get_exit_code() != 0) cerr << "\nError (" << e.get_exit_code() << ") in command line:\n";
      app.exit(e);
    }
    return false;
  }

  if (!paired_fnames.empty()) {
    // convert pairs to colon ':' separated single files for FastqReader to process
    if (paired_fnames.size() % 2 != 0) {
      if (!rank_me()) cerr << "Did not get pairs of files in -p: " << paired_fnames.size() << endl;
      return false;
    }
    while (paired_fnames.size() >= 2) {
      reads_fnames.push_back(paired_fnames[0] + ":" + paired_fnames[1]);
      paired_fnames.erase(paired_fnames.begin());
      paired_fnames.erase(paired_fnames.begin());
    }
  }

  // we can get restart incorrectly set in the config file from restarted runs
  if (*cfg_opt) {
    restart = false;
    app.get_option("--restart")->default_val("false");
  }

  if (!restart && reads_fnames.empty()) {
    if (!rank_me())
      cerr << "\nError in command line:\nRequire read names if not restarting\nRun with --help for more information\n";
    return false;
  }

  if (post_assm_only && ctgs_fname.empty()) ctgs_fname = "final_assembly.fasta";

  upcxx::barrier();

  if (!*output_dir_opt) {
    if (restart) {
      if (!rank_me())
        cerr << "\nError in command line:\nRequire output directory when restarting run\nRun with --help for more information\n";
      return false;
    }
    string first_read_fname = remove_file_ext(get_basename(reads_fnames[0]));
    output_dir = "mhm2-run-" + first_read_fname + "-n" + to_string(upcxx::rank_n()) + "-N" +
                 to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + get_current_time(true);
    output_dir_opt->default_val(output_dir);
  }

  if (restart) {
    try {
      app.parse_config(output_dir + "/mhm2.config");
    } catch (const CLI::ConfigError &e) {
      if (!upcxx::rank_me()) {
        cerr << "\nError (" << e.get_exit_code() << ") in config file (" << output_dir << "/mhm2.config" << "):\n";
        app.exit(e);
      }
      return false;
    }
  }

  // make sure we only use defaults for kmer lens if none of them were set by the user
  if (!*kmer_lens_opt && !*scaff_kmer_lens_opt) {
    kmer_lens = {21, 33, 55, 77, 99};
    scaff_kmer_lens = {99, 33};
    // set the option default strings so that the correct values are printed and saved to the .config file
    kmer_lens_opt->default_str("21 33 55 77 99");
    scaff_kmer_lens_opt->default_str("99 33");
  } else if (kmer_lens.size() && *kmer_lens_opt && !*scaff_kmer_lens_opt) {
    // use the last and second from kmer_lens for scaffolding
    auto n = kmer_lens.size();
    if (n == 1) {
      scaff_kmer_lens = kmer_lens;
      scaff_kmer_lens_opt->default_str(to_string(scaff_kmer_lens[0]));
    } else {
      scaff_kmer_lens = {kmer_lens[n - 1], kmer_lens[n == 2 ? 0 : 1]};
      scaff_kmer_lens_opt->default_str(to_string(scaff_kmer_lens[0]) + " " + to_string(scaff_kmer_lens[1]));
    }
  } else if (scaff_kmer_lens.size() == 1 && scaff_kmer_lens[0] == 0) {
    // disable scaffolding rounds
    scaff_kmer_lens.clear();
  }

  setup_output_dir();
  setup_log_file();

  if (max_kmer_store_mb == 0) {
    // use 1% of the minimum available memory
    max_kmer_store_mb = get_free_mem() / 1024 / 1024 / 100;
    max_kmer_store_mb = upcxx::reduce_all(max_kmer_store_mb / local_team().rank_n(), upcxx::op_fast_min).wait();
  }

  auto logger_t = chrono::high_resolution_clock::now();
  if (upcxx::local_team().rank_me() == 0) {
    // open 1 log per node
    // rank0 has mhm2.log in rundir, all others have logs in per_thread
    init_logger("mhm2.log", verbose, rank_me());
  }

  barrier();
  chrono::duration<double> logger_t_elapsed = chrono::high_resolution_clock::now() - logger_t;
  SLOG_VERBOSE("init_logger took ", setprecision(2), fixed, logger_t_elapsed.count(), " s at ", get_current_time(), "\n");

#ifdef DEBUG
  open_dbg("debug");
#endif

  SLOG(KLBLUE, full_version_str, KNORM, "\n");

  if (restart) get_restart_options();

  if (upcxx::rank_me() == 0) {
    // print out all compiler definitions
    SLOG_VERBOSE(KLBLUE, "_________________________", KNORM, "\n");
    SLOG_VERBOSE(KLBLUE, "Compiler definitions:", KNORM, "\n");
    std::istringstream all_defs_ss(ALL_DEFNS);
    vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
    for (auto &def : all_defs) SLOG_VERBOSE("  ", def, "\n");
    SLOG_VERBOSE("_________________________\n");
    SLOG(KLBLUE, "Options:", KNORM, "\n");
    auto all_opts_strings = app.config_to_str_vector(true, false);
    for (auto &opt_str : all_opts_strings) SLOG(KLBLUE, opt_str, KNORM, "\n");
    SLOG(KLBLUE, "_________________________", KNORM, "\n");
  }
  auto num_nodes = upcxx::rank_n() / upcxx::local_team().rank_n();
  SLOG("Starting run with ", upcxx::rank_n(), " processes on ", num_nodes, " node", (num_nodes > 1 ? "s" : ""), " at ",
       get_current_time(), "\n");
#ifdef DEBUG
  SWARN("Running low-performance debug mode");
#endif
  if (!upcxx::rank_me()) {
    // write out configuration file for restarts
    ofstream ofs("mhm2.config");
    ofs << app.config_to_str(true, true);
  }
  upcxx::barrier();
  return true;
}
