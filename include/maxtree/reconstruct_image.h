#pragma once

#include "../common.h"
#include "../misc/bit_array.h"
#include "../parallel/thread_pool.h"
#include "../parallel/iterative_select2_compact1.h"
#include "../misc/random.h"
#include "../misc/integerhash.h"
#include "../misc/edge.h"
//#include "rootfix2.h"
//#include "rootfix_seq.h"
#include <vector>
#include <limits>

NAMESPACE_PMT

template <typename index_t, typename value_t, typename functor_t>
void reconstruct_image(value_t const*values, size_t n, value_t* values_out, index_t* parents, functor_t const &criterion)
{
  index_t* roots = new index_t[n];
  index_t* node_indices = new index_t[n];
  index_t* aux = new index_t[n];
  
  IterativeSelect2Compact1<index_t> select(n);

  select.item_blocks().select([=](index_t i, index_t o) ALWAYS_INL_L(bool) {
    if (criterion(i)) {
      values_out[i] = values[i];
      roots[i] = i;
      return false;
    }

    node_indices[o] = i;
    roots[i] = parents[i];
    return true;
  });

  typename rng<index_t>::type r;
  IntegerHash<index_t> hash;

  std::vector<size_t> update_later;  

  size_t total_compacted = 0;
  while (select.length() > 0)
  {
    hash.generate_vars(r);

    auto const& f = [=](index_t const& x, index_t* out) ALWAYS_INL_L(uint_fast8_t)
    {
      *out = x;

      if (hash(x, 1) == 0)
      {
        return 1;
      }

      index_t root = roots[x];

      //check(root != x);

      if (hash(root, 1) == 1)
      {
        return 1;
      }

      return 2;
    };

    size_t n_compacted = 0;
    select.iterate(node_indices, aux + total_compacted, f, &n_compacted);
    update_later.push_back(n_compacted);
    total_compacted += n_compacted;

    select.item_blocks().apply([=](index_t i)
    {
      index_t x = node_indices[i];
      
      index_t root = roots[x];

      if (hash(root, 1) == 0)
      {
        return;
      }

      index_t root_of_root = roots[root];

      if (hash(root_of_root, 1) == 1)
      {
        return;
      }

      roots[x] = roots[roots[x]];
    });
  }

  index_t* nodes_to_update = aux + total_compacted;

  for (size_t i = update_later.size(); i--;)
  {
    size_t len = update_later[i];
    nodes_to_update -= len;

    thread_pool.for_all(len, [=](index_t i, thread_nr_t t) {
      index_t x = nodes_to_update[i];      
      values_out[x] = values_out[roots[x]];
    });
  }

  delete[] node_indices;
  delete[] aux;
  delete[] roots;
}

NAMESPACE_PMT_END