
#include <cstdint>
#include <iostream>

#include "../include/image/image_blocks.h"
#include "../include/common.h"
#include "../include/misc/timer.h"
#include "../include/maxtree/maxtree.h"
#include "../include/misc/edge.h"
#include "../include/maxtree/tree_scan.h"
#include "../include/maxtree/reconstruct_image.h"
#include "../include/maxtree/reconstruct_image_seq.h"

using index_t = uint32_t;
using attribute_t = uint32_t;

template <typename value_t>
void construct(index_t width, index_t height, index_t depth)
{
  index_t n = width * height * depth;
  size_t max_threads = pmt::thread_pool.max_threads();

  value_t* vals = new value_t[n];
  value_t* vals_out = new value_t[n];
  attribute_t* areas = new attribute_t[n];

  using image_t = typename pmt::image<index_t, value_t, 2, 4>::type;  
  image_t img(vals, {width, height});
  check(img.dimensions().length() == n);

  using rng = typename pmt::rng<index_t>::type;
  rng* rand = new rng[max_threads];

  index_t* parents = new index_t[n];

  // binary associative function
  auto const &plus =
    [](attribute_t a, attribute_t b) ALWAYS_INL_L(attribute_t)
    {
      return a + b;
    };

  auto const &w = [](index_t i) ALWAYS_INL_L(attribute_t)
    {
      return 1U;
    };


  size_t lambda = 10000;

  {
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

    pmt::maxtree(img, parents);
    pmt::tree_scan(parents, n, areas, w, plus);
  
    {
      auto const &criterion = [=](index_t i) ALWAYS_INL_L(bool)
      {
        return areas[i] > lambda;
      };

      pmt::reconstruct_image(vals, n, vals_out, parents, criterion);
    }
  }

  delete[] rand; 

  {
    auto const &criterion = [=](index_t i) ALWAYS_INL_L(bool)
    {
      return areas[i] > lambda;
    };
  
    {
      pmt::reconstruct_image_seq(vals, n, vals, parents, criterion);
    }
  }

  for (size_t i = 0; i < n; ++i)
  {
    check(vals[i] == vals_out[i]);
  }

  out("Seems correct.");


  delete[] parents;
  delete[] areas;
  delete[] vals;
  delete[] vals_out;
}

int main(int argc, char** argv)
{
  index_t width = 1024U;
  index_t height = width;
  index_t depth = 1U;
  
  construct<uint32_t>(width, height, depth);

  return 0;
}