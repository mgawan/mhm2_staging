#ifndef _AGGR_STORE_HPP
#define _AGGR_STORE_HPP

#include <deque>
#include <algorithm>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"

using std::vector;

using upcxx::intrank_t;
using upcxx::rank_me;
using upcxx::rank_n;
using upcxx::barrier;
using upcxx::dist_object;
using upcxx::rpc;
using upcxx::reduce_one;
using upcxx::reduce_all;
using upcxx::op_fast_add;
using upcxx::op_fast_max;
using upcxx::progress;
using upcxx::make_view;
using upcxx::view;
using upcxx::future;
using upcxx::make_future;
using upcxx::promise;

// this class aggregates updates into local buffers and then periodically does an rpc to dispatch them

template<typename T, typename... Data>
class AggrStore {
  using RankStore = vector<T>;
  using Store = vector<RankStore>;
  using UpdateFunc = std::function<void(T, Data&...)>;

  Store store;
  int64_t max_store_size;
  int max_rpcs_in_flight;
  vector<int> rpcs_sent;
  int64_t tot_rpcs_sent;
  dist_object<int64_t> tot_rpcs_expected;
  dist_object<int64_t> tot_rpcs_processed;

  // save the update function to use in both update and flush
  dist_object<UpdateFunc> update_func = dist_object<UpdateFunc>({});
  // save all associated data structures as a tuple of a variable number of parameters
  std::tuple<Data...> data;

  // operates on a vector of elements in the store
  static void update_remote(AggrStore *astore, intrank_t target_rank, Data &...data) {
    astore->tot_rpcs_sent++;
    astore->rpcs_sent[target_rank]++;
    rpc_ff(target_rank,
           [](dist_object<UpdateFunc> &update_func, view<T> rank_store, dist_object<int64_t> &tot_rpcs_processed,
              Data &...data) {
             (*tot_rpcs_processed)++;
             for (auto elem : rank_store) {
               (*update_func)(elem, data...);
             }
           }, astore->update_func, make_view(astore->store[target_rank].begin(), astore->store[target_rank].end()),
           astore->tot_rpcs_processed, data...);
    astore->store[target_rank].clear();
  }

public:

  AggrStore(Data&... data)
    : store({})
    , max_store_size(0)
    , rpcs_sent({})
    , tot_rpcs_sent(0)
    , tot_rpcs_expected(0)
    , tot_rpcs_processed(0)
    , max_rpcs_in_flight(AGGR_STORE_MAX_RPCS_IN_FLIGHT)
    , data(data...) {
    rpcs_sent.resize(upcxx::rank_n(), 0);
  }

  virtual ~AggrStore() {
    clear();
  }

  void set_size(const string &desc, int64_t max_store_bytes) {
    int64_t tmp_max_rpcs_in_flight = 500L * rank_n() / 100L + 1;
    // hard maximum limit
    max_rpcs_in_flight = tmp_max_rpcs_in_flight > max_rpcs_in_flight ? max_rpcs_in_flight : tmp_max_rpcs_in_flight;
    max_store_size = max_store_bytes / (sizeof(T) * (rank_n() + max_rpcs_in_flight));
    if (max_store_size <= 1) {
      // no reason for delay and storage of 1 entry (i.e. small max mem at large scale), still uses max_rpcs_in_flight
      max_store_size = 0;
      store.clear();
    } else {
      store.resize(rank_n(), {});
    }
    // reduce max in flight if necessary
    int64_t tmp_inflight_bytes = max_store_bytes - (max_store_size * rank_n() * sizeof(T)); // calc remaining memory for rpcs
    // hard minimum limit
    int64_t max_inflight_bytes = tmp_inflight_bytes > AGGR_STORE_MIN_INFLIGHT_BYTES ?
                                 tmp_inflight_bytes : AGGR_STORE_MIN_INFLIGHT_BYTES;
    int64_t per_rpc_bytes = (max_store_size > 0 ? max_store_size : 1) * sizeof(T);
    if (max_rpcs_in_flight * per_rpc_bytes > max_inflight_bytes)
      max_rpcs_in_flight = max_inflight_bytes / per_rpc_bytes + 1;
    size_t max_target_buf = max_store_size * sizeof(T);
    SLOG_VERBOSE(desc, ": using an aggregating store for each rank of max ", get_size_str(max_store_bytes / rank_n()),
                 " per target rank\n");
    SLOG_VERBOSE("  - buffers: max ", max_store_size, " entries of ", get_size_str(sizeof(T)),
         " per target rank (", get_size_str(max_target_buf), ")\n");
    SLOG_VERBOSE("  - buffers: max over all target ranks ", get_size_str(max_target_buf * rank_n()), "\n");
    SLOG_VERBOSE("  - RPCs in flight: max ", max_rpcs_in_flight, " RPCs of ", get_size_str(per_rpc_bytes),
         " max per RPC (", get_size_str(max_rpcs_in_flight * per_rpc_bytes), ")\n");
    SLOG_VERBOSE("  - max possible memory: ", get_size_str(max_target_buf * rank_n() + per_rpc_bytes * max_rpcs_in_flight),
                 "\n");
    barrier();
  }

  void set_update_func(UpdateFunc update_func) {
    *(this->update_func) = update_func;
  }

  void clear() {
    for (auto s : store) {
      if (!s.empty()) throw string("rank store is not empty!");
    }
    for (int i = 0; i < rpcs_sent.size(); i++) {
      rpcs_sent[i] = 0;
    }
    tot_rpcs_sent = 0;
    *tot_rpcs_expected = 0;
    *tot_rpcs_processed = 0;
  }

  void update(intrank_t target_rank, T &elem) {
    assert(max_store_size > 0);
    store[target_rank].push_back(elem);
    if (store[target_rank].size() < max_store_size) return;
    std::apply(update_remote, std::tuple_cat(std::make_tuple(this, target_rank), data));
  }

  void flush_updates() {
    for (int target_rank = 0; target_rank < rank_n(); target_rank++) {
      if (max_store_size > 0 && store[target_rank].size()) {
        std::apply(update_remote, std::tuple_cat(std::make_tuple(this, target_rank), data));
      }
      // tell the target how many rpcs we sent to it
      rpc(target_rank,
          [](dist_object<int64_t> &tot_rpcs_expected, int64_t rpcs_sent) {
            (*tot_rpcs_expected) += rpcs_sent;
          }, tot_rpcs_expected, rpcs_sent[target_rank]).wait();
    }
    barrier();
    // now wait for all of our rpcs
    while (*tot_rpcs_expected != *tot_rpcs_processed) {
      progress();
    }
    barrier();
    clear();
    barrier();
  }

};




#endif
