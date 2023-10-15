#pragma once

#include <pthread.h>
#include <iostream>
#include <atomic>
#include <thread>
#include "../misc/dimensions.h"
#include "../misc/coordinate.h"
#include "../common.h"
#include "../misc/range.h"
#include "../misc/random.h"

NAMESPACE_PMT

class ThreadPool
{
public:
  using range_t = Range<size_t>;

  ThreadPool(size_t max_threads);
  ~ThreadPool();
  
  void parallel(
    void (*user_f)(void *, thread_nr_t thread_nr),
    void *user_data);

  template <typename functor_t>
  void for_all_blocks(size_t n_blocks, functor_t const &user_block_f);

  template <typename primitives_t, typename functor_t>
  void for_all_blocks(
    Dimensions<primitives_t::n_dimensions> const &dims,
    functor_t const &user_block_f);

  template <typename functor_t>
  void for_all(
    size_t n,
    functor_t const &user_f,
    size_t n_items_per_block = default_n_items_per_block);

  constexpr size_t max_threads() { return max_threads_; }

  constexpr size_t n_active_threads() const { return n_active_threads_; }

  constexpr void set_n_active_threads(size_t n)
  {
    if (n == 0 || n > max_threads_)
    {
      n_active_threads_ = max_threads_;
      return;
    }

    n_active_threads_ = n;
  }  

private:
  enum state {running, terminating};

  struct ThreadDataUnpadded
  {
    std::atomic<Range<size_t>> range_;
    ThreadPool* pool_;
    pthread_t id_;
#ifdef PMT_DEBUG
    size_t concurr_reads_;
    size_t concurr_writes_;
#endif
    bool flag_;
    thread_nr_t thread_nr_;
    
  };

  struct ThreadData
  {
    union
    {
      ThreadDataUnpadded data_;
      // avoid false sharing
      uint8_t data_plus_padding_[
        ((sizeof(ThreadDataUnpadded) + cacheline_len - 1U) / cacheline_len) *
        cacheline_len];
    };
  };

  void wake_threads();
  void wait_ready();
  void create_threads();
  void join_threads();
  static void* thread_pool_loop(void* data);
  static bool find_work(ThreadData* td, thread_nr_t thread_nr, size_t n_threads, size_t* block_nr);

  template <typename functor_t>
  void iterate_range(
    size_t begin,
    size_t end,
    functor_t const &user_elem_f,
    thread_nr_t thread_nr);

  size_t max_threads_;  
  size_t n_active_threads_;
  ThreadData* thread_data_;
  pthread_mutex_t mutex_ = PTHREAD_MUTEX_INITIALIZER;
  pthread_cond_t cond_ = PTHREAD_COND_INITIALIZER;
  bool main_flag_{false};
  std::atomic<size_t> n_ready_{0U};
  state state_; // only changed by the main thread and if the mutex is locked
  void (*user_f_)(void*, size_t){nullptr};
  void* user_data_{nullptr};
  cpu_set_t cpu_set_;
  pthread_attr_t pthread_attr_;  
};

template <typename functor_t>
void ThreadPool::for_all_blocks(size_t n_blocks, functor_t const &user_block_f)
{
  using prim = primitives<size_t>;
  using dim_t = Dimensions<prim::n_dimensions>;
  using vec_t = Coordinate<prim>;

  for_all_blocks<prim>(dim_t({n_blocks}),
    [=](vec_t const& b, thread_nr_t thread_nr) ALWAYS_INLINE
    {
      user_block_f(b[0], thread_nr);
    });
}

template <typename prim, typename functor_t>
void ThreadPool::for_all_blocks(
  Dimensions<prim::n_dimensions> const &grid_dims,
  functor_t const &user_block_f)
{   
  using vec_t = Coordinate<prim>;
  using index_t = typename prim::index_t;

  struct SharedData
  {
    functor_t const &user_block_f_;
    Dimensions<prim::n_dimensions> const grid_dims_;
    ThreadData* thread_data_;
    size_t n_threads_;
  };

  size_t n_blocks = grid_dims.length();

  if (n_blocks == 0U)
  {
    return;
  }

  if (n_blocks == 1U)
  {
    vec_t v = vec_t::from_index(index_t(0), grid_dims);
    user_block_f(v, 0U);
    return;
  }

  if (n_active_threads_ == 1U)
  {
    vec_t v = vec_t::from_index(index_t(0), grid_dims);

    size_t begin = 0;
    size_t end = grid_dims.length();

    for (; begin < end; ++begin)
    {
      user_block_f(v, 0U);
      v.inc_index(grid_dims);
    }

    return;
  }

  size_t n_threads = std::min(n_active_threads_, n_blocks);

  SharedData data{user_block_f, grid_dims, thread_data_, n_threads};

#ifndef PMT_DEBUG
  size_t per_thread = n_blocks / n_threads;
  size_t remainder = n_blocks % n_threads;
  size_t begin = 0;
#endif  
  Range<size_t> range;
  
  for (size_t i = 0; i < n_threads; ++i)
  {
#ifdef PMT_DEBUG
    range = {0U, 0U};
#else  
    range.begin_ = begin;
    begin += per_thread + size_t(i < remainder);
    range.end_ = begin;
    thread_data_[i].data_.range_.store(range, std::memory_order_relaxed);
#endif
  }

#ifdef PMT_DEBUG
  thread_data_[0].data_.range_.store({0U, n_blocks}, std::memory_order_relaxed);
#endif

  user_f_ = [](void *data_p, thread_nr_t thread_nr)
  {      
    SharedData& data = *reinterpret_cast<SharedData *>(data_p);

    if (thread_nr >= data.n_threads_)
    {
      return;
    }

    size_t block_nr = 0U;

    while (true)
    {
      if (!find_work(data.thread_data_, thread_nr, data.n_threads_, &block_nr))
      {
        return;
      }

      vec_t v = vec_t::from_index(block_nr, data.grid_dims_);
      data.user_block_f_(v, thread_nr);
    }
  };

  parallel(user_f_, &data);

#ifdef PMT_DEBUG
  size_t total_cr = 0;
  size_t total_cw = 0;

  for (size_t i = 0; i < max_threads_; ++i)
  {
    printf("thread %ld; cr=%ld, cw=%ld\n",
      i,
      thread_data_[i].data_.concurr_reads_,
      thread_data_[i].data_.concurr_writes_);

    total_cr += thread_data_[i].data_.concurr_reads_;
    total_cw += thread_data_[i].data_.concurr_writes_;
  }

  printf("\ntotal cr=%ld, cw=%ld\n", total_cr, total_cw);
#endif
}

template <typename functor_t>
NO_INLINE void ThreadPool::iterate_range(
  size_t begin,
  size_t end,
  functor_t const &user_elem_f,
  thread_nr_t thread_nr)
{    
  for (; begin < end; ++begin)
  {
    user_elem_f(begin, thread_nr);
  }
}

template <typename functor_t>
void ThreadPool::for_all(
  size_t n,
  functor_t const& user_f,
  size_t n_items_per_block)
{
  size_t n_blocks = div_roundup(n, n_items_per_block);

  for_all_blocks(n_blocks, [=](size_t block_nr, thread_nr_t thread_nr) {
    size_t begin = block_nr * n_items_per_block;
    size_t end = begin + n_items_per_block;

    if (block_nr == n_blocks - 1U)
    {
      end = n;
    }

    iterate_range(begin, end, user_f, thread_nr);
  });
}

NAMESPACE_PMT_END