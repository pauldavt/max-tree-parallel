
#include <cstdint>
#include <iostream>

#include "../include/image/image_blocks.h"
#include "../include/common.h"
#include "../include/misc/timer.h"
#include "../include/maxtree/maxtree.h"
#include "../include/misc/edge.h"
#include "../include/maxtree/tree_scan.h"
#include "../include/maxtree/tree_scan_seq.h"

using index_t = uint32_t;

template <typename value_t>
void construct()
{
  index_t W = 1024U;
  index_t H = W;
  index_t D = 1U;
  index_t N = W * H * D;

  value_t* vals = new value_t[N];

  using image_t = typename pmt::image<index_t, value_t, 2, 4>::type;  
  image_t img(vals, {W, H});
  check(img.dimensions().length() == N);

  using rng = typename pmt::rng<index_t>::type;

  rng* rand = new rng[pmt::thread_pool.max_threads()];

  pmt::thread_pool.for_all(N, [=](index_t i, pmt::thread_nr_t t) {
    if (std::is_floating_point<value_t>::value)
    {
      vals[i] = pmt::random_fp(rand[t]);
    }
    else
    {
      vals[i] = rand[t]();
    }
  });

  index_t* parents = new index_t[N];

  {
    pmt::maxtree(img, parents);
  }

  index_t* area = new index_t[N];
  index_t* area2 = new index_t[N];

  {
    auto const &weight = [](index_t i) ALWAYS_INL_L(index_t) {
      return 1U;
    };

    auto const &plus = [](index_t a, index_t b) ALWAYS_INL_L(index_t) {
      return a + b;
    };

    pmt::tree_scan(parents, N, area, weight, plus);
  }

  {
    auto const &weight = [](index_t i) ALWAYS_INL_L(index_t) {
      return 1U;
    };

    auto const &plus = [](index_t a, index_t b) ALWAYS_INL_L(index_t) {
      return a + b;
    };

    pmt::tree_scan_seq(parents, N, area2, weight, plus);
  }

  for (size_t i = 0; i < N; ++i)
  {
    check(area[i] == area2[i]);
  }

  info("Areas match.");

  delete[] area;
  delete[] area2;
  delete[] parents;
  delete[] rand;
  delete[] vals;
}

int main(int argc, char** argv)
{
  construct<uint32_t>();

  return 0;
}