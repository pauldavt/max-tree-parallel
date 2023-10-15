#pragma once

#include "../common.h"
#include "../misc/bit_array.h"
#include "../parallel/thread_pool.h"
#include "../parallel/iterative_select2_compact1.h"
#include "../misc/random.h"
#include "../misc/integerhash.h"
#include "../misc/edge.h"
#include <vector>

NAMESPACE_PMT

template <typename index_t, typename functor_t>
ALWAYS_INLINE_F index_t compress_path(index_t i, index_t *roots, functor_t const &criterion)
{
  if (criterion(i))
  {
    return i;
  }

  index_t root = i;

  do
  {
    root = roots[root];
  } while (!criterion(root) && roots[root] != root);

  while (roots[i] != root)
  {
    index_t next = roots[i];
    roots[i] = root;
    i = next;
  }  

  return root;
}

template <typename index_t, typename value_t, typename functor_t>
void reconstruct_image_seq(value_t const*values, size_t n, value_t* values_out, index_t* parents, functor_t const &criterion)
{
  index_t* roots = new index_t[n];

  for (size_t i = 0; i < n; ++i)
  {
    roots[i] = parents[i];
  }

  for (size_t i = 0; i < n; ++i)
  {
    index_t root = compress_path(index_t(i), roots, criterion);

    values_out[i] = values[root];
  }

  delete[] roots;
}

NAMESPACE_PMT_END