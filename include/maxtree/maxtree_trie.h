#pragma once

#include "../common.h"
#include "../image/image.h"
#include "../misc/bit_array.h"
#include "../sort/sort_item.h"
#include "../sort/radix_sort_seq.h"
#include "../sort/sort.h"
#include "../misc/trie_queue.h"
#include "../misc/coordinate.h"

NAMESPACE_PMT

template <typename Primitives>
struct MaxtreeTrie;

template<typename Primitives>
void maxtree_trie(
  Image<Primitives> const& img,
  typename Primitives::index_t* parents,
  typename Primitives::index_t* rank_to_index,
  BitArray* visited,
  TrieQueue<typename Primitives::index_t>* queue)
{
  MaxtreeTrie<Primitives>(img, parents, rank_to_index, parents, visited, queue);
}

template<typename Primitives>
void maxtree_trie(Image<Primitives> const& img, typename Primitives::index_t* parents)
{
  using index_t = typename Primitives::index_t;
  using visited_t = BitArray;
  using queue_t = TrieQueue<index_t>;
  using value_t = typename Primitives::value_t;
  using uvalue_t = decltype(unsigned_conversion(value_t(0)));
  using dim_t = Dimensions<Primitives::n_dimensions>;
  using sort_pair_t = SortPair<uvalue_t, index_t>;

  size_t n = img.dimensions().length();
  const dim_t& dims = img.dimensions();
  value_t const* vals = img.values();

  size_t n_digits = radix_sort_n_digits<uvalue_t>();
  
  sort_pair_t* alloc = nullptr;
  sort_pair_t* aux1 = nullptr;
  sort_pair_t* aux2 = nullptr;
  if (n_digits == 2)
  {
    alloc = new sort_pair_t[dims.length()];
    aux1 = alloc;
  }
  else if (n_digits > 2)
  {
    alloc = new sort_pair_t[2U * dims.length()];
    aux1 = alloc;
    aux2 = alloc + dims.length();
  }

  auto const& f_initial_item = [=](index_t i) ALWAYS_INL_L(sort_pair_t)
  {
    return {unsigned_conversion(vals[i]), i};
  };

  auto const& f_out = [=](index_t& out, sort_pair_t const& item) ALWAYS_INLINE
  {
    out = item.data();
  };

  index_t* rank_to_index = new index_t[n];
  index_t* index_to_rank = parents;

  radix_sort_seq<index_t>(rank_to_index, n, aux1, aux2, f_initial_item, f_out);  

  delete[] alloc;

  for (size_t rank = 0; rank != n; ++rank)
  {
    index_to_rank[rank_to_index[rank]] = rank;
  }

  visited_t visited(n);
  queue_t queue(n - size_t(1), false);

  MaxtreeTrie<Primitives>(img, parents, rank_to_index, index_to_rank, &visited, &queue);

  delete[] rank_to_index;
}

template <typename Primitives>
struct MaxtreeTrie
{
private:
  using index_t = typename Primitives::index_t;
  using image_t = Image<Primitives>;
  using dim_t = Dimensions<Primitives::n_dimensions>;
  using vec_t = Coordinate<Primitives>;
  using visited_t = BitArray;
  using queue_t = TrieQueue<index_t>;
  using value_t = typename Primitives::value_t;
  using uvalue_t = decltype(unsigned_conversion(value_t(0)));
  using sort_pair_t = SortPair<uvalue_t, index_t>;

  friend void maxtree_trie<Primitives>(
    Image<Primitives> const& img,
    typename Primitives::index_t* parents,
    typename Primitives::index_t* rank_to_index,
    BitArray* visited,
    TrieQueue<typename Primitives::index_t>* queue);

  friend void maxtree_trie<Primitives>(
    Image<Primitives> const& img,
    typename Primitives::index_t* parents);

  MaxtreeTrie(
    image_t const& img,
    index_t* parents,
    index_t* rank_to_index,
    index_t* index_to_rank,
    visited_t* visited,
    queue_t* queue);
    
  bool check_neighbor(
    index_t current,
    index_t current_rank,
    index_t next,
    index_t &next_rank);
  bool next_neighbor(
    index_t current,
    index_t current_rank,
    index_t const *skip,
    index_t &next,
    index_t &next_rank);
  void main_loop();

  image_t const& image_;
  index_t* parents_;
  index_t* RESTRICT rank_to_index_;
  index_t* index_to_rank_;
  visited_t& visited_;
  queue_t& queue_;
};

template <typename Primitives>
MaxtreeTrie<Primitives>::MaxtreeTrie(
  image_t const& img,
  index_t* parents,
  index_t* rank_to_index,
  index_t* index_to_rank,
  visited_t* visited,
  queue_t* queue) : 
  image_(img),
  parents_(parents),
  rank_to_index_(rank_to_index),
  index_to_rank_(index_to_rank),
  visited_(*visited),
  queue_(*queue)
{
  static_assert(Primitives::n_dimensions > 0, "");
  
  static_assert(
    Primitives::n_neighbors == 2U * Primitives::n_dimensions ||
    (Primitives::n_neighbors == 8 && Primitives::n_dimensions == 2), "");

  visited_.clear();
  queue_.clear();

  main_loop();
}

template <typename prim>
ALWAYS_INLINE_F bool
MaxtreeTrie<prim>::check_neighbor(
  index_t current,
  index_t current_rank,
  index_t next,
  index_t &next_rank)
{
  if (visited_.is_set(next))
  {
    return false;
  }

  visited_.set(next);

  index_t l_next_rank = index_to_rank_[next];

  if (l_next_rank <= current_rank)
  {
    queue_.insert(l_next_rank);
    return false;
  }

  next_rank = l_next_rank;
  return true;
}

#define CHECK_NEIGHBOR(cond, offset) \
{ \
  index_t l_next = current + offset; \
  if ((cond) && \
    check_neighbor(current, current_rank, l_next, next_rank)) \
  { \
    next = l_next; \
    return true; \
  } \
}

template <typename prim>
ALWAYS_INLINE_F bool
MaxtreeTrie<prim>::next_neighbor(
  index_t current,
  index_t current_rank,
  index_t const* skip,
  index_t &next,
  index_t &next_rank)
{
  const dim_t& dims = image_.dimensions();
  vec_t pos = vec_t::from_index(current, dims);

  CHECK_NEIGHBOR(pos[0] > 0U, -index_t(1))
  CHECK_NEIGHBOR(pos[0] < dims[0] - size_t(1), index_t(1))

  if (prim::n_dimensions == 2 && prim::n_neighbors == 8)
  {
    index_t width = skip[0];

    if (pos[1] > 0)
    {
      CHECK_NEIGHBOR(pos[0] > 0U, -width -index_t(1))
      CHECK_NEIGHBOR(true, -width)
      CHECK_NEIGHBOR(pos[0] < width - size_t(1), -width + index_t(1))
    }

    if (pos[1] < dims[1] - size_t(1))
    {
      CHECK_NEIGHBOR(pos[0] > 0U, width -index_t(1))
      CHECK_NEIGHBOR(true, width)
      CHECK_NEIGHBOR(pos[0] < width - size_t(1), width + index_t(1))
    }

    return false;
  }

  for (dim_idx_t d = 1; d < prim::n_dimensions; ++d)
  {
    index_t distance = skip[d - dim_idx_t(1)];

    CHECK_NEIGHBOR(pos[d] > 0, -distance)
    CHECK_NEIGHBOR(pos[d] < dims[d] - size_t(1), distance)
  }

  return false;
}

template <typename Primitives>
void MaxtreeTrie<Primitives>::main_loop()
{
  index_t current_rank = 0;
  index_t current = rank_to_index_[current_rank];
  parents_[current] = current;
  visited_.set(current);

  const dim_t& dims = image_.dimensions();
  index_t skip[Primitives::n_dimensions - dim_idx_t(1)];

  if (Primitives::n_dimensions > 1)
  {
    skip[0] = dims[0];
  }

  for (dim_idx_t d = 1; d < Primitives::n_dimensions - dim_idx_t(1); ++d)
  {
    skip[d] = dims[d] * skip[d - dim_idx_t(1)];
  }

  while (true)
  {
    index_t next;
    index_t next_rank;

    bool ret = next_neighbor(current, current_rank, skip, next, next_rank);

    if (ret)
    {
      // moving to a higher level
      queue_.insert(current_rank);
      current = next;
      current_rank = next_rank;
      continue;
    }

    if (queue_.empty()) break;

    index_t parent_rank = queue_.top();
    index_t parent = rank_to_index_[parent_rank];
    queue_.remove();
    parents_[current] = parent;      
    current = parent;
    current_rank = parent_rank;
  }
}

NAMESPACE_PMT_END