#include "../include/parallel/thread_pool.h"
#include "../include/misc/dimensions.h"
#include "../include/misc/coordinate.h"
#include "../include/misc/logger.h"

constexpr unsigned N = 1024U * 1024U * 64U;

NO_INLINE size_t sum_array(uint32_t* values, uint32_t n)
{
  size_t sum = 0;

  for (unsigned i = 0; i < n; ++i)
  {
    sum += values[i];
  }

  return sum;
}

int main()
{
  uint32_t* vals = new uint32_t[N];

  size_t sum = 0;

  for (unsigned i = 0; i < N; ++i)
  {
    vals[i] = i;
    sum += i;
  }

  size_t max_threads = pmt::thread_pool.max_threads();
  out("max threads = " << max_threads);

  size_t partial_sums[max_threads];
  size_t* partial_sums_ptr = partial_sums;

  for (unsigned i = 0; i < max_threads; ++i)
  {
    partial_sums[i] = 0;
  }

  pmt::thread_pool.for_all(N, [=](size_t i, pmt::thread_nr_t thread_nr) ALWAYS_INLINE {
    partial_sums_ptr[thread_nr] += i;
  });

  size_t sum2 = 0;
  for (unsigned i = 0; i < max_threads; ++i)
  {
    sum2 += partial_sums[i];
  }

  check(sum == sum2);

  out("Success.");

  delete[] vals;
}
