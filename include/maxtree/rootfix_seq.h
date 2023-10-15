#pragma once

#include "../common.h"
#include "../misc/bit_array.h"
#include "../parallel/thread_pool.h"
#include "../parallel/iterative_select2_compact1.h"
#include "../misc/random.h"
#include "../misc/integerhash.h"
#include "../misc/edge.h"
#include "../misc/dynamic_stack.h"
#include <vector>

NAMESPACE_PMT

template <typename index_t, typename attribute_t, typename functor1_t, typename functor2_t>
void rootfix_seq(index_t* parents, size_t n, attribute_t *RESTRICT attributes, functor1_t const& w, functor2_t const& plus)
{
  index_t *RESTRICT roots = new index_t[n];

  for (size_t i = 0; i < n; ++i)
  {
    roots[i] = parents[i];
    attributes[i] = w(i);
  }

  index_t on_stack[4096];
  pmt::DynamicStack<index_t> stack(&on_stack[0], 4096);

  for (size_t i = 0; i < n; ++i)
  {
    index_t root = i;

    while (root != roots[root])
    {
      stack.insert(root);
      root = roots[root];
    }

    attribute_t curr_attr = attributes[root];

    while (stack.length())
    {
      index_t current = stack.remove();
      curr_attr = plus(curr_attr, attributes[current]);
      attributes[current] = curr_attr;
      roots[current] = current;
    }
  }

  delete[] roots;
}

NAMESPACE_PMT_END