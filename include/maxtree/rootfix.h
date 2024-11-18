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

template <typename index_t, typename attribute_t, typename functor1_t, typename functor2_t>
void rootfix(index_t* parents, size_t n, attribute_t *attributes, functor1_t const& w, functor2_t const& plus)
{
//  using edge_t = Edge<index_t>;

  IterativeSelect2Compact1<index_t> select(n);

  index_t* RESTRICT roots = new index_t[n];
  index_t* RESTRICT node_indices = new index_t[2U * n];
  index_t* RESTRICT aux = node_indices + n;
  char* RESTRICT merging = new char[n];

  select.item_blocks().apply([=](index_t i) ALWAYS_INLINE
    {
      roots[i] = parents[i];
      node_indices[i] = i;
      attributes[i] = w(i);
    });

  typename rng<index_t>::type r;
  IntegerHash<index_t> hash;
  std::vector<size_t> update_later;  
  size_t total_compacted = 0;

  while (select.length() > 0)
  {
    hash.generate_vars(r);

    auto const& f =
      [=](index_t const& x, index_t* out) ALWAYS_INL_L(uint_fast8_t)
      {
        index_t root = roots[x];

        bool is_merging = hash(x, 1) != 0 && hash(root, 1) == 0;

        merging[x] = is_merging;
        
        if (root == x) return 0;

        *out = x;

        if (is_merging)
        {
          return 2;
        }

        return 1;
      };

    size_t n_compacted = 0;
    select.iterate(node_indices, aux + total_compacted, f, &n_compacted);
    update_later.push_back(n_compacted);
    total_compacted += n_compacted;

    select.item_blocks().apply([=](index_t i)
    {
      index_t x = node_indices[i];
      index_t root = roots[x];

      if (!merging[root])
      {
        return;
      }

      index_t root_root = roots[root];
      
      attributes[x] = plus(attributes[root], attributes[x]);
      roots[x] = root_root;
    });
  }

  index_t* nodes_to_update = aux + total_compacted;

  for (size_t i = update_later.size(); i--;)
  {
    size_t len = update_later[i];
    nodes_to_update -= len;

    thread_pool.for_all(len, [=](index_t i, thread_nr_t t) {
      index_t x = nodes_to_update[i];
      index_t root = roots[x];
      attributes[x] = plus(attributes[root], attributes[x]);
    });
  }

  delete[] merging;
  delete[] node_indices;
  delete[] roots;
}

NAMESPACE_PMT_END