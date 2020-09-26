#pragma once

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

#include <string>
#include <string_view>
#include <utility>
#include <vector>

using std::pair;
using std::string;
using std::string_view;
using std::vector;

#ifdef USE_BYTELL
#include "bytell_hash_map.hpp"
#define HASH_TABLE ska::bytell_hash_map
#else
#include <unordered_map>
#define HASH_TABLE std::unordered_map
#endif

size_t estimate_hashtable_memory(size_t num_elements, size_t element_size);

string revcomp(const string &seq);

char comp_nucleotide(char ch);

int hamming_dist(string_view s1, string_view s2, bool require_equal_len = true);

string get_merged_reads_fname(string reads_fname);

void switch_orient(int &start, int &stop, int &len);

void dump_single_file(const string &fname, const string &out_str, bool append = false);

vector<string> get_dir_entries(const string &dname, const string &prefix);

string &left_trim(string &str);

int pin_clear();

string get_proc_pin();

vector<int> get_pinned_cpus();

void pin_proc(vector<int> cpus);

void pin_cpu();

void pin_core();

void pin_numa();

// temporary until it is properly within upcxx_utils
#include <deque>
#include <functional>
#include <memory>
#include <thread>

#include <upcxx/upcxx.hpp>

namespace upcxx_utils {

/*
 * TODO #include <future> // for std::async
 * requires removal of using namespace std across the board
 * templates must be in header files
 * and std::future makes calls to upcxx::future ambiguous (same with promise)
 */

/*
#include <future>
template <typename F, typename... Ts>
inline auto reallyAsync(F&& f, Ts&&... params) {
    return std::async(std::launch::async, std::forward<F>(f),
                std::forward<Ts>(params)...);
}

template<typename Func, typename... Args>
auto execute_async(Func&& func, Args&&... args) {
    auto t_start = Timer::now();
    auto sh_prom = make_shared< upcxx::promise <> > ();
    const upcxx::persona &persona = upcxx::current_persona();
    DBG("Starting async sh_prom=", sh_prom.get(), "\n");

    auto returning_func =
    [t_start, sh_prom, &persona] (Func&& f, Args&& ... a) {
        auto ret = f(a...);
        duration_seconds sec = Timer::now() - t_start;
        DBG("Completed running async sh_prom=", sh_prom.get(), " in ", sec.count(), "s\n");

        // fulfill only in calling persona
        persona.lpc_ff([t_start, sh_prom]() {
            duration_seconds sec = Timer::now() - t_start;
            DBG_VERBOSE("Fulfilling promised async sh_prom=", sh_prom.get(), " in ", sec.count(), "s\n");
            sh_prom->fulfill_anonymous(1);
        });

        sec = Timer::now() - t_start;
        DBG_VERBOSE("Returning async result sh_prom=", sh_prom.get(), " in ", sec.count(), "s\n");
        return ret;
    };
    auto async_fut = reallyAsync(returning_func, args...);

    return sh_prom->get_future().then(
            [t_start, sh_prom, &persona, async_fut]() {
                assert(persona.active_with_caller());

                duration_seconds sec = Timer::now() - t_start;
                DBG_VERBOSE("Waiting for completed async ", async_fut.valid(), " sh_prom=", sh_prom.get(), " in ", sec.count(), "s\n");

                async_fut.wait(); // should be noop but there could be a short race between lpc and returning the value
                assert(async_fut.valid());

                sec = Timer::now() - t_start;
                DBG("Returning completed async ", async_fut.valid(), " sh_prom=", sh_prom.get(), " in ", sec.count(), "s\n");
                return async_fut.get();
            });
}

template <typename Func>
upcxx::future<> execute_in_new_thread(Func &&func) {
  return execute_async(func);
};

*/

// Func no argument returned or given lambda - void()

template <typename Func>
upcxx::future<> execute_in_new_thread(upcxx::persona &persona, Func func) {
  assert(persona.active_with_caller());
  std::shared_ptr<upcxx::promise<> > sh_prom = std::make_shared<upcxx::promise<> >();

  std::shared_ptr<std::thread> sh_run = std::make_shared<std::thread>([&persona, func, sh_prom] {
    func();

    // fulfill only in calling persona
    persona.lpc_ff([&persona, sh_prom]() {
      assert(persona.active_with_caller());
      sh_prom->fulfill_anonymous(1);
    });
  });
  auto fut_finished = sh_prom->get_future().then(  // keep sh_prom and sh_run alive until complete
      [&persona, sh_prom, sh_run]() {
        assert(persona.active_with_caller());
        sh_run->join();
      });
  return fut_finished;
}

template <typename Func>
upcxx::future<> execute_in_new_thread(Func func) {
  return execute_in_new_thread(upcxx::current_persona(), func);
}

/**/

};  // namespace upcxx_utils
