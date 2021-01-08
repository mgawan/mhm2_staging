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
#include <string_view>
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
  if (output_dir.empty()) DIE("Invalid empty ouput_dir");
  if (!upcxx::rank_me()) {
    // create the output directory (and possibly stripe it)

    if (restart) {
      // it must already exist for a restart
      if (access(output_dir.c_str(), F_OK) == -1) {
        SDIE("Output directory ", output_dir, " for restart does not exist");
      }
    }

    // always try to make the output_dir
    if (mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */) == -1) {
      // could not create the directory
      if (errno == EEXIST) {
        // okay but warn if not restarting
        if (!restart) SWARN("Output directory ", output_dir, " already exists. May overwrite existing files\n");
      } else {
        SDIE("Could not create output directory '", output_dir, "': ", strerror(errno), "\n");
      }
    }

    // always ensure striping is set or reset wide when lustre is available
    auto status = std::system("which lfs 2>&1 > /dev/null");
    bool set_lfs_stripe = WIFEXITED(status) & (WEXITSTATUS(status) == 0);
    int num_osts = 0;
    if (set_lfs_stripe) {
      // detect the number of OSTs with "lfs df -l $output_dir" and count the entries with OST:
      string cmd("lfs df -l ");
      cmd += output_dir;
      FILE *f_osts = popen(cmd.c_str(), "r");
      if (f_osts) {
        char buf[256];
        while (char *ret = fgets(buf, 256, f_osts)) {
          std::string_view s(buf);
          if (s.find("OST:") != string::npos) num_osts++;
        }
        fclose(f_osts);
      }
      // reduce to the minimum of 90% or rank_n()
      num_osts = std::min(9 * num_osts / 10, std::min((int)72, (int)rank_n()));  // see Issue #70
    }
    set_lfs_stripe &= num_osts > 0;
    if (set_lfs_stripe) {
      // stripe with count -1 to use all the OSTs
      string cmd = "lfs setstripe --stripe-count " + std::to_string(num_osts) + " --stripe-size 16M " + output_dir;
      auto status = std::system(cmd.c_str());
      if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
        cout << "Set Lustre striping on the output directory: count=" << num_osts << ", size=16M\n";
      else
        cout << "Failed to set Lustre striping on output directory: " << WEXITSTATUS(status) << endl;
    }

    // ensure per_rank dir exists (and additionally has stripe 1, if lfs exists)
    string per_rank = output_dir + "/per_rank";
    status = mkdir(per_rank.c_str(), S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);
    if (status != 0 && errno != EEXIST) SDIE("Could not create '", per_rank, "'! ", strerror(errno));
    if (set_lfs_stripe) {
      // stripe per_rank directory with count 1
      string cmd = "lfs setstripe --stripe-count 1 " + per_rank;
      status = std::system(cmd.c_str());
      if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
        cout << "Set Lustre striping on the per_rank output directory: " << per_rank << "\n";
      else
        cout << "Failed to set Lustre striping on per_rank output directory: " << WEXITSTATUS(status) << endl;
    }

    // this should avoid contention on the filesystem when ranks start racing to creating these top levels
    char basepath[256];
    for (int i = 0; i < rank_n(); i += 1000) {
      sprintf(basepath, "%s/%08d", per_rank.c_str(), i);
      status = mkdir(basepath, S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);
      if (status != 0 && errno != EEXIST) {
        SWARN("Could not create '", basepath, "'! ", strerror(errno));
        break;  // ignore any errors, just stop
      }
    }
    sprintf(basepath, "%s/00000000/%08d", per_rank.c_str(), 0);
    status = mkdir(basepath, S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);
    if (status != 0 && errno != EEXIST) SDIE("Could not mkdir rank 0 per thread directory '", basepath, "'! ", strerror(errno));

    cout << "Using output dir: " << output_dir << "\n";  // required for mhm2.py to find the output_dir when not verbose!

    cout.flush();
    cerr.flush();
  }

  upcxx::barrier();
  // after we change to the output directory, relative paths could be incorrect, so make sure we have the correct path of the
  // reads files
  char cwd_str[FILENAME_MAX];
  if (!getcwd(cwd_str, FILENAME_MAX)) DIE("Cannot get current working directory: ", strerror(errno));
  for (auto &fname : reads_fnames) {
    if (fname[0] != '/') {
      string dir = string(cwd_str) + "/";
      auto spos = fname.find_first_of(':');
      if (spos == string::npos || spos == fname.size() - 1) {
        // interleaved (no colon) or single read (colon at the end)
        fname = dir + fname;
      } else {
        // paired read (colon in between)
        fname = dir + fname.substr(0, spos) + ":" + dir + fname.substr(spos + 1);
      }
    }
  }
  // all change to the output directory
  auto chdir_attempts = 0;
  while (chdir(output_dir.c_str()) != 0) {
    // failed, retry for 5 more seconds - Issue #69
    if (chdir_attempts++ > 10) DIE("Cannot change to output directory ", output_dir, ": ", strerror(errno));
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }
  upcxx::barrier();
}

void Options::setup_log_file() {
  if (!upcxx::rank_me()) {
    // check to see if mhm2.log exists. If so, and not restarting, rename it
    if (file_exists("mhm2.log") && !restart) {
      string new_log_fname = "mhm2-" + setup_time + ".log";
      cerr << KLRED << "WARNING: " << KNORM << output_dir << "/mhm2.log exists. Renaming to " << output_dir << "/" << new_log_fname
           << endl;
      if (rename("mhm2.log", new_log_fname.c_str()) == -1) DIE("Could not rename mhm2.log: ", strerror(errno));
      // also unlink the rank0 per_rank file (a hard link if it exists)
      unlink("per_rank/00000000/00000000/mhm2.log");  // ignore any errors
    } else if (!file_exists("mhm2.log") && restart) {
      DIE("Could not restart - missing mhm2.log in this directory");
    }
  }
  upcxx::barrier();
}

Options::Options() {
  char buf[32];
  if (!upcxx::rank_me()) {
    setup_time = get_current_time(true);
    strncpy(buf, setup_time.c_str(), sizeof(buf) - 1);
  }
  upcxx::broadcast(buf, sizeof(buf), 0, world()).wait();
  setup_time = string(buf);
  output_dir = string("mhm2-run-<reads_fname[0]>-n") + to_string(upcxx::rank_n()) + "-N" +
               to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + setup_time;
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
  app.add_option("-r, --reads", reads_fnames,
                 "Files containing interleaved paired reads in FASTQ format (comma or space separated).")
      ->delimiter(',')
      ->check(CLI::ExistingFile);
  app.add_option("-p, --paired-reads", paired_fnames,
                 "Alternating read files containing separate paired reads in FASTQ format (comma or space separated).")
      ->delimiter(',')
      ->check(CLI::ExistingFile);
  app.add_option("-u, --unpaired-reads", unpaired_fnames, "Unpaired or single reads in FASTQ format (comma or space separated).")
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
  app.add_flag("--checkpoint-merged", checkpoint_merged,
               "(debugging option) enables checkpointing of merged fastq files in the output directory")
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
    // convert pairs to colon ':' separated unpaired/single files for FastqReader to process
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

  if (!unpaired_fnames.empty()) {
    // append a ':' to the file name, signaling to FastqReader that this is a unpaired/single file
    // a paired file would have another name after the ':'
    for (auto name : unpaired_fnames) {
      reads_fnames.push_back(name + ":");
    }
    unpaired_fnames.clear();
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
    string first_read_fname = reads_fnames[0];
    // strip out the paired or unpaired/single the possible ':' in the name string
    auto colpos = first_read_fname.find_last_of(':');
    if (colpos != string::npos) first_read_fname = first_read_fname.substr(0, colpos);
    first_read_fname = remove_file_ext(get_basename(first_read_fname));
    auto spos = first_read_fname.find_first_of(':');
    if (spos != string::npos) first_read_fname = first_read_fname.substr(0, spos);
    output_dir = "mhm2-run-" + first_read_fname + "-n" + to_string(upcxx::rank_n()) + "-N" +
                 to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + setup_time;
    output_dir_opt->default_val(output_dir);
  }

  setup_output_dir();
  setup_log_file();

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

  // save to per_rank, but hardlink to output_dir
  string config_file = "per_rank/mhm2.config";
  string linked_config_file = "mhm2.config";
  if (restart) {
    // use per_rank to read/write this small file, hardlink to top run level
    try {
      app.parse_config(linked_config_file);
    } catch (const CLI::ConfigError &e) {
      if (!upcxx::rank_me()) {
        cerr << "\nError (" << e.get_exit_code() << ") in config file (" << config_file << "):\n";
        app.exit(e);
      }
      return false;
    }
  }

  if (max_kmer_store_mb == 0) {
    // use 1% of the minimum available memory
    max_kmer_store_mb = get_free_mem() / 1024 / 1024 / 100;
    max_kmer_store_mb = upcxx::reduce_all(max_kmer_store_mb / upcxx::local_team().rank_n(), upcxx::op_fast_min).wait();
  }

  auto logger_t = chrono::high_resolution_clock::now();
  if (upcxx::local_team().rank_me() == 0) {
    // open 1 log per node
    // all have logs in per_rank
    if (rank_me() == 0 && restart) {
      auto ret = rename("mhm2.log", "per_rank/00000000/00000000/mhm2.log");
      if (ret != 0)
        SWARN("For this restart, could not rename mhm2.log to per_rank/00000000/00000000/mhm2.log: ", strerror(errno), "\n");
    }
    init_logger("mhm2.log", verbose, true);
    // if not restarting, hardlink just the rank0 log to the output dir
    // this ensures a stripe count of 1 even when the output dir is striped wide
    if (rank_me() == 0) {
      auto ret = link("per_rank/00000000/00000000/mhm2.log", "mhm2.log");
      if (ret != 0) SWARN("Could not hard link mhm2.log from per_rank/00000000/00000000/mhm2.log: ", strerror(errno), "\n");
    }
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
    ofstream ofs(config_file);
    ofs << app.config_to_str(true, true);
    ofs.close();
    unlink(linked_config_file.c_str());  // ignore errors
    auto ret = link(config_file.c_str(), linked_config_file.c_str());
    if (ret != 0 && !restart) LOG("Could not hard link config file, continuing\n");
  }
  upcxx::barrier();
  return true;
}
