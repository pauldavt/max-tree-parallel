#pragma once
#include "../common.h"
#include "../misc/edge.h"

NAMESPACE_PMT

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
class EulerTourScan;

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
void euler_tour_scan(
  index_t *parents,
  size_t n,
  attribute_t *attributes,
  functor1_t const &weight,
  functor2_t const &plus,
  functor3_t const &inverse,
  attribute_t const &identity)
{
  EulerTourScan<
    index_t,
    attribute_t,
    functor1_t,
    functor2_t,
    functor3_t> ets(
    parents,
    n,
    attributes,
    weight,
    plus,
    inverse,
    identity);
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
class EulerTourScan
{  
  using edge_t = pmt::SortableEdgeByStart<index_t>;

  //static constexpr size_t items_per_block = 8192U;

  static constexpr unsigned n_hash_bits = 8U;

  friend void euler_tour_scan<
    index_t,
    attribute_t,
    functor1_t,
    functor2_t,
    functor3_t>(
    index_t *parents,
    size_t n,
    attribute_t *attributes,
    functor1_t const &weight,
    functor2_t const &plus,
    functor3_t const &inverse,
    attribute_t const &identity);

  EulerTourScan(
    index_t *parents,
    size_t n,
    attribute_t *attributes,
    functor1_t const &weight,
    functor2_t const &plus,
    functor3_t const &inverse,
    attribute_t identity);

  ~EulerTourScan();

  void create_forward_edges();
  void euler_tour();
  void linked_list_scan();

  index_t *parents_;
  size_t n_;
  size_t ll_n_;
  attribute_t *attributes_;    
  attribute_t *ll_attributes_;  

  edge_t *forward_;  
  union
  {
    edge_t *aux_;
    index_t *next_;
    index_t *roots_;
  };
  index_t *first_;
  index_t *nodes_;
  index_t *update_later_;
  functor1_t const &weight_;
  functor2_t const &plus_;
  functor3_t const &inverse_;
  attribute_t identity_;
  index_t root_;
  index_t begin_;
  index_t end_;
};

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
EulerTourScan<
  index_t,
  attribute_t,
  functor1_t,
  functor2_t,
  functor3_t>::
EulerTourScan(
  index_t *parents,
  size_t n,
  attribute_t *attributes,
  functor1_t const &weight,
  functor2_t const &plus,
  functor3_t const &inverse,
  attribute_t identity) :
  parents_(parents),
  n_(n),
  attributes_(attributes),
  weight_(weight),
  plus_(plus),
  inverse_(inverse),
  identity_(identity)
{
  forward_ = new edge_t[n];
  aux_ = new edge_t[n];
  
  create_forward_edges();

  next_ = reinterpret_cast<index_t *>(aux_);
  first_ = new index_t[n];

  euler_tour();

  delete[] first_;

  ll_n_ = 2U * n - 1U;
  ll_attributes_ = new attribute_t[ll_n_];
  nodes_ = new index_t[ll_n_];
  update_later_ = new index_t[ll_n_];

  linked_list_scan();

  // compute branch differences
  thread_pool.for_all(n - 1U, [=](index_t i, thread_nr_t t) {    
    edge_t const& current = forward_[i];
    index_t y = current.b_;

    attributes_[y] = plus(ll_attributes_[y], inverse(ll_attributes_[n + i]));    
  });

  attributes_[root_] = ll_attributes_[root_];

  delete[] update_later_;
  delete[] nodes_;
  delete[] ll_attributes_;
  delete[] aux_;
  delete[] forward_;
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
EulerTourScan<
  index_t,
  attribute_t,
  functor1_t,
  functor2_t,
  functor3_t>::
~EulerTourScan()
{  
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
void EulerTourScan<
  index_t,
  attribute_t,
  functor1_t,
  functor2_t,
  functor3_t>::
euler_tour()
{
  thread_pool.for_all(n_, [=](index_t i, thread_nr_t thread_nr) {
    next_[i] = parents_[i];
    first_[i] = i;
  });

  thread_pool.for_all(n_ - 1U, [=](index_t i, thread_nr_t thread_nr) {
    edge_t const& current = forward_[i];
    edge_t const& right = forward_[i + 1];

    if (right.a_ != current.a_)
    {
      first_[current.a_] = n_ + i;
    }
    else
    {
      next_[right.b_] = n_ + i;
    }
  });

  thread_pool.for_all(n_ - 1U, [=](index_t i, thread_nr_t thread_nr) {
    next_[n_ + i] = first_[forward_[i].b_];
  });

  end_ = root_;
  begin_ = first_[end_];
  next_[end_] = begin_;
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
void EulerTourScan<
  index_t,
  attribute_t,
  functor1_t,
  functor2_t,
  functor3_t>::
create_forward_edges()
{
  auto const& f_initial_item = [=](index_t i) ALWAYS_INL_L(edge_t)
  {
    if (parents_[i] == i)
    {
      return {~index_t(0), i};
    }

    return {parents_[i], i};
  };

  edge_t* sorted = pmt::radix_sort_parallel(forward_, aux_, n_, f_initial_item);  

  if (sorted != forward_)
  {
    std::swap(forward_, aux_);
  }

  root_ = forward_[n_ - 1U].b_;
}

#define PMT_MERGEABLE(x) (x != begin_ && x != end_ && hash(x, n_hash_bits) != 0)

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t,
  typename functor3_t>
void EulerTourScan<
  index_t,
  attribute_t,
  functor1_t,
  functor2_t,
  functor3_t>::
linked_list_scan()
{ 
  size_t ll_n = 2U * n_ - 1U;

  thread_pool.for_all(
    ll_n,
    [=](index_t i, thread_nr_t thread_nr) ALWAYS_INLINE
    {
      ll_attributes_[i] = i < n_ ? weight_(i) : identity_;
      nodes_[i] = i;
    });

  size_t total_compacted = 0;
  std::vector<size_t> to_update;
  typename pmt::rng<index_t>::type rng;
  IntegerHash<index_t> hash;
  IterativeSelect2Compact1<index_t> select(ll_n);

  while (select.length() > 2)
  {
    hash.generate_vars(rng);

    auto const& f =
      [=](index_t x, index_t* out) ALWAYS_INL_L(uint_fast8_t)
      {
        bool mergeable = PMT_MERGEABLE(x);

        *out = x;

        if (mergeable) return 2U;

        index_t current = roots_[x];
        while (PMT_MERGEABLE(current))
        {
          index_t root = roots_[current];
          ll_attributes_[root] = plus_(ll_attributes_[root], ll_attributes_[current]);
          roots_[current] = x;
          current = root;
          root = roots_[current];
        }

        roots_[x] = current;

        return 1U;
      };

    size_t n_compacted = 0;
    select.iterate(nodes_, update_later_ + total_compacted, f, &n_compacted);
    to_update.push_back(n_compacted);
    total_compacted += n_compacted;
  }

  ll_attributes_[end_] = plus_(ll_attributes_[end_], ll_attributes_[begin_]);

  index_t* nodes_to_update = update_later_ + total_compacted;

  for (size_t i = to_update.size(); i--;)
  {
    size_t len = to_update[i];
    nodes_to_update -= len;

    thread_pool.for_all(len, [=](index_t k, thread_nr_t t) 
      {
        index_t x = nodes_to_update[k];
        ll_attributes_[x] = plus_(ll_attributes_[x], ll_attributes_[roots_[x]]);
      });
  }
}

NAMESPACE_PMT_END