#include "../include/image/image_blocks.h"
#include "../include/misc/random.h"
#include "../include/parallel/thread_pool.h"

using index_t = uint32_t;

int main()
{
  using prim = pmt::primitives<index_t, index_t, 3>;
  using image_t = pmt::Image<prim>;
  using image_blocks_t = pmt::ImageBlocks<prim>;

  index_t h = 100U;
  index_t w = 100U;
  index_t d = 100U;
  index_t n = h * w * d;

  index_t* data = new index_t[n];

  for (size_t i = 0; i < n; ++i)
  {
    data[i] = i;
  }

  typename pmt::rng<index_t>::type r;

  pmt::random_shuffle(data, n, r);


  image_t image(data, {w, h, d});
  image_blocks_t image_blocks(image);


  size_t max_threads = pmt::thread_pool.max_threads();
  size_t* sums = new size_t[max_threads];

  for (size_t i = 0; i < max_threads; ++i)
  {
    sums[i] = 0;
  }

  using vec_t = pmt::Coordinate<prim>;
  pmt::thread_pool.for_all_blocks<prim>(image_blocks.dimensions(), [=, &image_blocks](vec_t loc, pmt::thread_nr_t t) {
    pmt::ImageBlock<prim> block(image_blocks, loc);
    
    block.apply([=](index_t global_index, index_t local_index) {
      sums[t] += data[global_index];
    });
  });

  size_t sum = 0;
  for (size_t i = 0; i < max_threads; ++i)
  {
    sum += sums[i];
  }

  check(sum == (n - 1) * size_t(n) / 2);

  info("Success.");

  delete[] sums;
}