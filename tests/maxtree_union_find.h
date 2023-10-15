#pragma once

#include "../include/common.h"
#include "../include/misc/unsigned_conversion.h"
#include "../include/sort/radix_sort_seq.h"

NAMESPACE_PMT

#define UNION_FIND_CHECK_NEIGHBOR(offset) \
{ \
  neighbor = find_root<index_t>(roots, current + (offset)); \
  if (neighbor != current && values[neighbor] >= values[current]) \
  { \
    roots[neighbor] = current; \
    parents[neighbor] = current; \
  } \
}

template <typename index_t>
index_t find_root(index_t *roots, index_t i)
{
  index_t root = i;

  while (roots[root] != root)
  {
    root = roots[root];
  }

  while (roots[i] != root)
  {
    index_t tmp = roots[i];
    roots[i] = root;
    i = tmp;    
  }

  return root;
}

template <typename value_t, typename index_t>
void maxtree_union_find_2d_4n(
  value_t *values,
  index_t *rank_to_index,
  index_t *parents,
  size_t width,
  size_t height)
{
  size_t n = width * height;
  index_t *roots = new index_t[n];

  for (size_t i = 0; i < n; ++i)
  {
    roots[i] = i;
  }

  for (size_t i = n; i--;)
  {
    index_t current = find_root(roots, rank_to_index[i]); 
    index_t neighbor;
    index_t x = current % width;
    index_t y = current / width;

    if (x > 0)           UNION_FIND_CHECK_NEIGHBOR(-1U)
    if (x < width - 1U)  UNION_FIND_CHECK_NEIGHBOR(1U)
    if (y > 0)           UNION_FIND_CHECK_NEIGHBOR(-width)
    if (y < height - 1U) UNION_FIND_CHECK_NEIGHBOR(width)        
  }

  index_t root = find_root(roots, rank_to_index[0]);

  parents[root] = root;

  delete[] roots;
}

template <typename value_t, typename index_t>
void maxtree_union_find_3d_6n(
  value_t *values,
  index_t *rank_to_index,
  index_t *parents,
  size_t width,
  size_t height,
  size_t depth)
{
  size_t n = width * height * depth;
  index_t *roots = new index_t[n];

  for (size_t i = 0; i < n; ++i)
  {
    roots[i] = i;
  }

  for (size_t i = n; i--;)
  {
    index_t current = find_root(roots, rank_to_index[i]); 
    index_t neighbor;
    index_t x = current % width;
    index_t y = (current % (width * height)) / width;
    index_t z = current / (width * height);

    if (x > 0)           UNION_FIND_CHECK_NEIGHBOR(-1U)
    if (x < width - 1U)  UNION_FIND_CHECK_NEIGHBOR(1U)
    if (y > 0)           UNION_FIND_CHECK_NEIGHBOR(-width)
    if (y < height - 1U) UNION_FIND_CHECK_NEIGHBOR(width)        
    if (z > 0)           UNION_FIND_CHECK_NEIGHBOR(-(width * height))
    if (z < depth - 1U)  UNION_FIND_CHECK_NEIGHBOR(width * height) 
  }

  index_t root = find_root(roots, rank_to_index[0]);

  parents[root] = root;

  delete[] roots;
}


template <typename value_t, typename index_t>
void maxtree_union_find_2d_8n(
  value_t *values,
  index_t *rank_to_index,
  index_t *parents,
  size_t width,
  size_t height)
{
  size_t n = width * height;
  index_t *roots = new index_t[n];

  for (size_t i = 0; i < n; ++i)
  {
    roots[i] = i;
  }

  for (size_t i = n; i--;)
  {
    index_t current = find_root(roots, rank_to_index[i]); 
    index_t neighbor;
    index_t x = current % width;
    index_t y = current / width;

    if (x > 0)                             UNION_FIND_CHECK_NEIGHBOR(-1U)
    if (x < width - 1U)                    UNION_FIND_CHECK_NEIGHBOR(1U)

    if (x > 0 && y > 0)                    UNION_FIND_CHECK_NEIGHBOR(-width - 1U)
    if (y > 0)                             UNION_FIND_CHECK_NEIGHBOR(-width)
    if (x < width - 1U && y > 0)           UNION_FIND_CHECK_NEIGHBOR(-width + 1U)
    
    if (x > 0 && y < height - 1U)          UNION_FIND_CHECK_NEIGHBOR(width - 1U)        
    if (y < height - 1U)                   UNION_FIND_CHECK_NEIGHBOR(width)        
    if (x < width - 1U && y < height - 1U) UNION_FIND_CHECK_NEIGHBOR(width + 1U)        
  }

  index_t root = find_root(roots, rank_to_index[0]);

  parents[root] = root;

  delete[] roots;
}

template <typename value_t, typename index_t>
void maxtree_union_find(
  value_t *values,
  index_t *parents,
  size_t width,
  size_t height,
  size_t depth,
  size_t n_neighbors = 0)
{
  size_t n = width * height * depth;
  check(n != 0);

  using uvalue_t = decltype(pmt::unsigned_conversion(value_t(0)));
  using sort_pair_t = pmt::SortPair<uvalue_t, index_t>;

  sort_pair_t *aux1 = nullptr;
  sort_pair_t *aux2 = nullptr;

  size_t n_digits = pmt::radix_sort_n_digits<uvalue_t>();

  if (n_digits > 1)
  {
    aux1 = new sort_pair_t[n];

    if (n_digits > 2)
    {
      aux2 = new sort_pair_t[n];
    }
  }

  auto const &f_initial_item = [=](index_t i) ALWAYS_INL_L(sort_pair_t)
  {
    return {pmt::unsigned_conversion(values[i]), i};
  };

  auto const &f_out = [=](index_t &out, sort_pair_t const &item) ALWAYS_INLINE
  {
    out = item.data();
  };

  index_t *rank_to_index = new index_t[n];

  pmt::radix_sort_seq(rank_to_index, n, aux1, aux2, f_initial_item, f_out);    

  delete[] aux2;
  delete[] aux1;

  if (depth == 1U)
  {
    if (n_neighbors == 8U)
    {
      maxtree_union_find_2d_8n(values, rank_to_index, parents, width, height);
    }
    else
    {
      maxtree_union_find_2d_4n(values, rank_to_index, parents, width, height);
    }
  }
  else
  {
    maxtree_union_find_3d_6n(values, rank_to_index, parents, width, height, depth);
  }

  delete[] rank_to_index;
}


NAMESPACE_PMT_END