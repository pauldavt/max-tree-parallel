#include "../../include/parallel/thread_pool.h"
#include "../../include/misc/logger.h"

NAMESPACE_PMT

size_t hardware_concurrency = std::thread::hardware_concurrency();
//size_t hardware_concurrency = 128U;

ThreadPool thread_pool(std::min(
  std::max(hardware_concurrency, size_t(1)),
  default_max_threads_limit));

ThreadPool::ThreadPool(size_t max_threads) :
  max_threads_(max_threads),
  n_active_threads_(max_threads),
  state_(running)
{
  thread_data_ = reinterpret_cast<ThreadData*>(
    aligned_alloc(cacheline_len, sizeof(ThreadData) * max_threads));

  check(thread_data_ != nullptr);

  pthread_attr_init(&pthread_attr_); 

  create_threads();
  wait_ready();
  // mutex locked
}

ThreadPool::~ThreadPool()
{
  // mutex locked
  state_ = terminating;
  wake_threads();
  join_threads();

  free(thread_data_);
  pthread_attr_destroy(&pthread_attr_); 
}

void ThreadPool::parallel(void (*user_f)(void *, size_t), void *user_data)
{
  // mutex locked 

  cpu_set_t saved;
  pthread_getaffinity_np(thread_data_[0].data_.id_, sizeof(cpu_set_t), &saved);
  CPU_ZERO(&cpu_set_);
  CPU_SET(0U, &cpu_set_);
  pthread_setaffinity_np(thread_data_[0].data_.id_, sizeof(cpu_set_t), &cpu_set_);
  memcpy(&cpu_set_, &saved, sizeof(cpu_set_t));

  main_flag_ = !main_flag_;
  user_data_ = user_data;
  user_f_ = user_f;
  
  wake_threads();    
  // mutex unlocked

  user_f_(user_data, 0U);

  wait_ready();
  // mutex locked

  pthread_setaffinity_np(thread_data_[0].data_.id_, sizeof(cpu_set_t), &cpu_set_);
}  

void ThreadPool::wake_threads()
{
  // mutex locked
  pthread_mutex_unlock(&mutex_);
  pthread_cond_broadcast(&cond_);
  // mutex unlocked
}

void ThreadPool::wait_ready()
{
  // mutex unlocked
  while (n_ready_.load(std::memory_order_consume) != max_threads_ - 1U)
  {
    // busy wait
  }

  pthread_mutex_lock(&mutex_);  
  // all other threads are waiting on the condition var
  n_ready_.store(0, std::memory_order_relaxed);
  // mutex locked
}

void ThreadPool::create_threads()
{
  ThreadDataUnpadded &td = thread_data_[0].data_;
  td.thread_nr_ = 0;
  td.pool_ = this;
  td.id_ = pthread_self();
  td.flag_ = false;
#ifdef PMT_DEBUG
  td.concurr_reads_ = 0;
  td.concurr_writes_ = 0;
#endif

  // outside a parallel operation, thread 0 can use all cores.
  // in a parallel operation, thread 0 has an affinity for core 0.

  for (size_t i = 1; i < max_threads_; ++i)
  {   
    ThreadDataUnpadded &td = thread_data_[i].data_;
    td.thread_nr_ = i;
    td.pool_ = this;
    td.flag_ = false;
#ifdef PMT_DEBUG
  td.concurr_reads_ = 0;
  td.concurr_writes_ = 0;
#endif

    CPU_ZERO(&cpu_set_);
    CPU_SET(i, &cpu_set_);
    pthread_attr_setaffinity_np(&pthread_attr_, sizeof(cpu_set_t), &cpu_set_);
    
    pthread_create(
      &td.id_,
      &pthread_attr_,
      thread_pool_loop,
      &td);
  }  
}

void ThreadPool::join_threads()
{
  for (size_t i = 1; i < max_threads_; ++i)
  {   
    ThreadDataUnpadded &td = thread_data_[i].data_;
    pthread_join(td.id_, nullptr);
  }       
}

void* ThreadPool::thread_pool_loop(void* data)
{
  ThreadDataUnpadded &td = *reinterpret_cast<ThreadDataUnpadded*>(data);
  ThreadPool &tp = *td.pool_;

  while (true)
  {
    pthread_mutex_lock(&tp.mutex_);

    tp.n_ready_.fetch_add(1U, std::memory_order_release);

    do
    {      
      if (tp.state_ == ThreadPool::terminating)
      {
        pthread_mutex_unlock(&tp.mutex_);
        return nullptr;
      }

      pthread_cond_wait(&tp.cond_, &tp.mutex_);
    } while (tp.main_flag_ == td.flag_); // avoid spurious wakeup

    pthread_mutex_unlock(&tp.mutex_);

    td.flag_ = !td.flag_;

    if (td.thread_nr_ >= tp.n_active_threads_)
    {
      continue;
    }

    tp.user_f_(tp.user_data_, td.thread_nr_);
  }
}

bool ThreadPool::find_work(
  ThreadData* tds,
  thread_nr_t thread_nr,
  size_t n_threads,
  size_t* block_nr)
{
  ThreadDataUnpadded& td = tds[thread_nr].data_;
  range_t current = td.range_.load(std::memory_order_relaxed);
#ifdef PMT_DEBUG
  ++td.concurr_reads_;
#endif

  range_t desired;

  do 
  {
    if (current.begin_ >= current.end_)
    {
      thread_nr_t k = thread_nr + 1U;

      if (k == n_threads)
      {
        k = 0U;
      }

      size_t i = 0;
      for (; i < n_threads - 1U; ++i)
      {
        if (k == thread_nr)
        {
          continue;
        }

        ThreadDataUnpadded& other = tds[k].data_;

        range_t other_range = other.range_.load(std::memory_order_relaxed);
#ifdef PMT_DEBUG
        ++td.concurr_reads_;
#endif
        range_t other_desired;
        size_t length;

        do
        {
          length = other_range.end_ - other_range.begin_;

          if (length == 0U)
          {
            break;
          }

          size_t half = (length + 1U) / 2U;
          other_desired = {other_range.begin_, other_range.end_ - half};
#ifdef PMT_DEBUG
          ++td.concurr_reads_;
          ++td.concurr_writes_;
#endif
        } while (!other.range_.compare_exchange_weak(
            other_range,
            other_desired,
            std::memory_order_relaxed));

        if (length > 0U)
        {
          *block_nr = other_desired.end_;
          range_t new_range{other_desired.end_ + 1U, other_range.end_};
          td.range_.store(new_range, std::memory_order_relaxed);
#ifdef PMT_DEBUG
          ++td.concurr_writes_;
#endif

          return true;          
        }

        ++k;
        if (k == n_threads)
        {
          k = 0U;
        }
      }

      if (i == n_threads - 1U)
      {
        return false;
      }
    }

    desired = {current.begin_ + 1U, current.end_};
#ifdef PMT_DEBUG
    ++td.concurr_reads_;
    ++td.concurr_writes_;
#endif
  } while (!td.range_.compare_exchange_weak(current, desired, std::memory_order_relaxed));     

  *block_nr = current.begin_;
  return true;
}

NAMESPACE_PMT_END