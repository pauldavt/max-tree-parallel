#pragma once

#include "../common.h"

NAMESPACE_PMT

template <typename Index>
struct RankSet
{
  static constexpr size_t n_digits = div_roundup(sizeof(Index) * CHAR_BIT, 8U);

  ALWAYS_INLINE_F void reset(Index k)
  {
    store(k);
    is_set_root_ = true;
    rank_ = 0;

    debug(load() == k);
  }

  ALWAYS_INLINE_F void store(Index k)
  {
    data_ = k;
  }

  ALWAYS_INLINE_F Index load() const
  {
    return data_;
  }

  Index data_;

  // For 64-bit indices: rank is log2(64) worst case when same sized sets merge.
  uint8_t rank_:7;
  // True if set root, and cc_root_ refers to the connected component root.
  // False if not a set root, and set_root_ refers to the next set on the
  // path to the set root.
  uint8_t is_set_root_:1;
};

template <typename Index>
Index compress_path(RankSet<Index>* sets, Index i)
{
  if (sets[i].is_set_root_)
  {
    return i;
  }

  Index root = sets[i].load();
  while (!sets[root].is_set_root_)
  {
    root = sets[root].load();
  }

  while (!sets[sets[i].load()].is_set_root_)
  {
    Index tmp = sets[i].load();
    sets[i].store(root);
    i = tmp;
  }

  return root;
}

template <typename Index>
void merge_sets(RankSet<Index>* sets, Index set_a, Index set_b, Index cc_root)
{
  uint_fast8_t rank_a = sets[set_a].rank_;
  uint_fast8_t rank_b = sets[set_b].rank_;
 
  if (rank_b > rank_a)
  {
    sets[set_a].store(set_b);
    sets[set_a].is_set_root_ = false;
    sets[set_b].store(cc_root);
    return;
  }

  sets[set_b].store(set_a);
  sets[set_b].is_set_root_ = false;
  sets[set_a].store(cc_root);

  if (rank_b < rank_a)
  {
    return;
  }

  // rank_a == rank_b  
  //chk(rank_a < 64U);

  ++sets[set_a].rank_;
}

NAMESPACE_PMT_END