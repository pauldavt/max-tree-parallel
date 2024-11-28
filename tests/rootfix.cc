
#include <cstdint>
#include <iostream>

#include "../include/image/image_blocks.h"
#include "../include/common.h"
#include "../include/misc/timer.h"
#include "../include/maxtree/maxtree.h"
#include "../include/misc/edge.h"
#include "../include/maxtree/tree_scan.h"
#include "../include/maxtree/rootfix.h"
#include "../include/maxtree/rootfix_seq.h"
#include "../include/maxtree/euler_tour_scan.h"

using index_t = uint32_t;
using attribute_t = uint32_t;

template <typename value_t>
void construct(index_t width, index_t height, index_t depth)
{
  index_t n = width * height * depth;
  size_t max_threads = pmt::thread_pool.max_threads();

  value_t* vals = new value_t[n];
  value_t* vals_out = new value_t[n];
  value_t* vals_out2 = new value_t[n];

  using image_t = typename pmt::image<index_t, value_t, 2, 4>::type;  
  image_t img(vals, {width, height});
  check(img.dimensions().length() == n);

  using rng = typename pmt::rng<index_t>::type;
  rng* rand = new rng[max_threads];

  pmt::thread_pool.for_all(n, [=](index_t i, pmt::thread_nr_t t) {
    if (std::is_floating_point<value_t>::value)
    {
      vals[i] = pmt::random_fp(rand[t]);
    }
    else
    {
      vals[i] = rand[t]();
    }
  });

  delete[] rand;

  index_t* parents = new index_t[n];

  pmt::maxtree(img, parents);

  // binary associative function
  auto const &plus =
    [](attribute_t a, attribute_t b) ALWAYS_INL_L(attribute_t)
    {
      return a + b;
    };
  {
    auto const &w = [=](index_t i) ALWAYS_INL_L(attribute_t)
    {
      index_t parent = parents[i];

      if (parent == i)
      {
        return vals[i];
      }

      return vals[i] - vals[parent];
    };
  
    {
      pmt::rootfix(parents, n, vals_out, w, plus);
    }

    {
      pmt::rootfix_seq(parents, n, vals_out2, w, plus);
    }
  }

  for (size_t i = 0; i < n; ++i)
  {
    check(vals_out2[i] == vals[i]);
    check(vals_out2[i] == vals_out[i]);
  }

  out("Seems correct.");

  delete[] parents;
  delete[] vals;
  delete[] vals_out;
  delete[] vals_out2;
}

int main(int argc, char** argv)
{
  index_t width = 1024U;
  index_t height = width;
  index_t depth = 1U;
  
  construct<uint32_t>(width, height, depth);

  return 0;
}