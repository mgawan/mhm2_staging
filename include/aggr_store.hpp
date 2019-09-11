#ifndef _AGGR_STORE_HPP
#define _AGGR_STORE_HPP

#include <deque>
#include <algorithm>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"

using std::vector;
using std::deque;

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

template<typename T>
class AggrStore {
private:
  using RankStore = vector<T>;
  using Store = vector<RankStore>;
  using RpcFutures = deque< future<> >;

  Store store;
  RpcFutures rpc_futures;
  int64_t max_store_size;
  int max_rpcs_in_flight; // Limit for the number of rpcs in flight. This limit exists to prevent the dispatch buffers from growing indefinitely

  void wait_max_rpcs() {
    progress();
    while (!rpc_futures.empty() && rpc_futures.size() >= max_rpcs_in_flight) {
      rpc_futures.front().wait();
      rpc_futures.pop_front();
    }
  }

  // operation on 1 element (i.e. no store)
  template<typename FuncDistObj, typename ...Args>  
  void update_remote1(intrank_t target_rank, T &elem, FuncDistObj &func, Args &...args) {
    wait_max_rpcs();
    auto fut = rpc(target_rank, 
                   [](FuncDistObj &func, T elem, Args &...args) {
                     (*func)(elem, args...);
                   }, func, elem, args...);

    rpc_futures.push_back(fut);
  }

  // operate on a vector of elements in the store
  template<typename FuncDistObj, typename ...Args>  
  void update_remote(intrank_t target_rank, FuncDistObj &func, Args &...args) {
    wait_max_rpcs();
    auto fut = rpc(target_rank, 
                   [](FuncDistObj &func, view<T> rank_store, Args &...args) {
                     for (auto elem : rank_store) (*func)(elem, args...);
                   }, func, make_view(store[target_rank].begin(), store[target_rank].end()), args...);

    rpc_futures.push_back(fut);
    store[target_rank].clear();
  }
  
public:
  
  AggrStore() : store({}), rpc_futures(), max_store_size(0), max_rpcs_in_flight(MAX_RPCS_IN_FLIGHT) {}
  virtual ~AggrStore() {
    clear();
  }

  void set_size(int64_t max_store_bytes) {
    int64_t tmp_max_rpcs_in_flight = 300L * rank_n() / 100L + 1; // assume 150% of rpcs as there are ranks
    max_rpcs_in_flight = tmp_max_rpcs_in_flight > max_rpcs_in_flight ? max_rpcs_in_flight : tmp_max_rpcs_in_flight; // hard maximum limit
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
    int64_t max_inflight_bytes = tmp_inflight_bytes > MIN_INFLIGHT_BYTES ? tmp_inflight_bytes : MIN_INFLIGHT_BYTES; // hard minimum limit
    int64_t per_rpc_bytes = (max_store_size > 0 ? max_store_size : 1) * sizeof(T);
    if (max_rpcs_in_flight * per_rpc_bytes > max_inflight_bytes) {
        max_rpcs_in_flight = max_inflight_bytes / per_rpc_bytes + 1;
    }
    SOUT("Using a store of max ", get_size_str(max_store_bytes / rank_n()), " per target rank, giving max ",
         max_store_size, " of ", get_size_str(sizeof(T)), " entries per target rank (",get_size_str(max_store_size*sizeof(T)),", ",
         get_size_str(max_store_size * sizeof(T) * rank_n()), ") and ", max_rpcs_in_flight,
         " rpcs in flight (", get_size_str(max_rpcs_in_flight * per_rpc_bytes ), ")\n");
    SOUT("Max possible memory for store for rank 0: ", get_size_str(max_store_size * sizeof(T) * rank_n() + per_rpc_bytes * max_rpcs_in_flight), "\n");
    barrier();
  }

  void clear() {
    for (auto s : store) {
      if (!s.empty()) throw string("rank store is not empty!");
    }
    if (!rpc_futures.empty()) throw string("rpc_futures are not empty!");
    Store().swap(store);
    RpcFutures().swap(rpc_futures);
  }

  template<typename FuncDistObj, typename ...Args>  
  void update(intrank_t target_rank, T &elem, FuncDistObj &func, Args &...args) {
    if (max_store_size <= 0) {
        update_remote1(target_rank, elem, func, args...);
        return;
    }
    store[target_rank].push_back(elem);
    if (store[target_rank].size() < max_store_size) return;
    update_remote(target_rank, func, args...);
  }

  template<typename FuncDistObj, typename ...Args>  
  void flush_updates(FuncDistObj &func, Args &...args) {
    for (int target_rank = 0; target_rank < rank_n(); target_rank++) {
      if (max_store_size > 0 && store[target_rank].size()) update_remote(target_rank, func, args...);
    }
    for (auto fut : rpc_futures) fut.wait();
    rpc_futures.clear();
    barrier();
  }

};




#endif
