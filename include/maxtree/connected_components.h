#pragma once

#include "../common.h"
#include "../misc/edge.h"
#include "../misc/integerhash.h"
#include "../misc/random.h"
#include "../parallel/iterative_select2_compact1.h"
#include <vector>

NAMESPACE_PMT

template <typename index_t, typename value_t>
class ConnectedComponents;

/* 
 * Used to determine connected components of a node weighted graph.
 * The edges and aux arrays have at least length n_edges.
 * Weight function is given by values[x].
 * A component root is the node with minimal tuple (values[x], x).
 * Component root references will be stored in the roots array.
 */
template <typename index_t, typename value_t>
void connected_components(
  Edge<index_t>* edges,
  size_t n_edges,
  Edge<index_t>* aux,
  value_t const* values,
  index_t* roots)
{
  ConnectedComponents<index_t, value_t> cc(edges, n_edges, aux, values, roots);
}

template <typename index_t, typename value_t>
class ConnectedComponents
{
private:
  friend void connected_components<index_t, value_t>(
    Edge<index_t>* edges,
    size_t n_edges,
    Edge<index_t>* aux,
    value_t const* values,
    index_t* roots);

  using edge_t = Edge<index_t>;
  using select_t = IterativeSelect2Compact1<edge_t>;
  static constexpr uint_fast8_t remove = select_t::remove;
  static constexpr uint_fast8_t keep = select_t::copy_to_first_array;
  static constexpr uint_fast8_t remove_update_later = select_t::copy_to_second_array;
  static constexpr size_t block_length = 256U * 1024U / sizeof(edge_t);

  ConnectedComponents(
    edge_t* edges,
    size_t n_edges,
    edge_t* aux,
    value_t const* values,
    index_t* roots);

  void update_roots();
  void change_roots_to_minima();
  void contract();
  void update_edges();

  edge_t* edges_;
  edge_t* aux_;
  value_t const* RESTRICT values_;
  index_t* RESTRICT roots_;
  IterativeSelect2Compact1<edge_t> select_;
  IntegerHash<index_t> hash_;
  size_t total_compacted_;
  std::vector<size_t> update_later_;
};

template <typename index_t, typename value_t>
ConnectedComponents<index_t, value_t>::ConnectedComponents(
  edge_t* edges,
  size_t n_edges,
  edge_t* aux,
  value_t const* values,
  index_t* roots) :
  edges_(edges),
  aux_(aux),
  values_(values),
  roots_(roots),
  select_(n_edges, block_length)
{
  typename rng<index_t>::type r;    

  total_compacted_ = 0;

  while (select_.length() > 0)
  {
    hash_.generate_vars(r);

    contract();
    update_edges();
  }

  update_roots();
  change_roots_to_minima();
}

template <typename index_t, typename value_t>
void ConnectedComponents<index_t, value_t>::update_roots()
{
  edge_t* to_update = aux_ + total_compacted_;
  edge_t* to_copy = edges_ + total_compacted_;

  for (size_t i = update_later_.size(); i--;)
  {
    size_t n = update_later_[i];
    to_update -= n;
    to_copy -= n;

    thread_pool.for_all(n, [=](size_t k, thread_nr_t t) ALWAYS_INLINE {
      edge_t& edge = to_update[k];

      debug(roots_[roots_[edge.b_]] == roots_[edge.b_]);
      
      roots_[edge.a_] = roots_[edge.b_];
      edge = {edge.a_, roots_[edge.b_]};
      to_copy[k] = {edge.a_, roots_[edge.b_]};
    });
  }

  check(to_update == aux_);

#ifdef PMT_DEBUG
  parallel_for_all(total_compacted_,
    [=](size_t i, thread_nr_t t) ALWAYS_INLINE
  {
    edge_t const& edge = aux_[i];
    check(roots_[edge.b_] == edge.b_);
    check(roots_[edge.a_] == edge.b_);
    check(aux_[i].a_ == edges_[i].a_ && aux_[i].b_ == edges_[i].b_);
  }); 
#endif
}

template <typename index_t, typename value_t>
void ConnectedComponents<index_t, value_t>::change_roots_to_minima()
{
  ItemBlocks& ib = select_.item_blocks();    
  ib.reset(total_compacted_);   

  while (ib.length() > 0)
  {      
    ib.select([=](size_t i, size_t o) ALWAYS_INL_L(bool)
    {
      edge_t edge = aux_[i];
      
      index_t min_root = roots_[edge.b_];

      index_t candidate = min_root;
      value_t val_a = values_[edge.a_];
      value_t val_candidate = values_[candidate];
      
      if (val_a > val_candidate || (val_a == val_candidate && edge.a_ >= candidate))
      {
        return false;
      }

      roots_[edge.b_] = edge.a_;
      aux_[o] = edge; 
      return true;
    });
  }

#ifdef PMT_DEBUG
  parallel_for_all(total_compacted_,
    [=](size_t i, thread_nr_t t) ALWAYS_INLINE
  {
    edge_t const& edge = edges_[i];
    check(roots_[edge.a_] == edge.b_);
    check(values_[roots_[edge.b_]] <= values_[edge.a_]);
  }); 
#endif      
  
  thread_pool.for_all(total_compacted_,
    [=](size_t i, thread_nr_t t) ALWAYS_INLINE
  {
    edge_t& edge = edges_[i];      

    roots_[edge.a_] = roots_[edge.b_];
  });

#ifdef PMT_DEBUG
  parallel_for_all(total_compacted_,
    [=](size_t i, thread_nr_t t) ALWAYS_INLINE
  {
    edge_t const& edge = edges_[i];

    check(roots_[roots_[edge.a_]] == roots_[edge.a_]);
    check(values_[roots_[edge.a_]] <= values_[edge.a_])
  });    
#endif          

}

template <typename index_t, typename value_t>
void ConnectedComponents<index_t, value_t>::contract()
{
  select_.item_blocks().apply([=](size_t i) ALWAYS_INLINE
  {
    edge_t const& edge = edges_[i];

    bool hash_a = hash_(edge.a_, 1);
    bool hash_b = hash_(edge.b_, 1);

    if (!(hash_a ^ hash_b)) return;

    if (hash_a)
    {
      roots_[edge.b_] = edge.a_;
      return;
    }

    roots_[edge.a_] = edge.b_;
  });
}

template <typename index_t, typename value_t>
void ConnectedComponents<index_t, value_t>::update_edges()
{
  auto const& f = [=](
    edge_t const& edge,
    edge_t* out) ALWAYS_INL_L(uint_fast8_t)
  {
    bool hash_a = hash_(edge.a_, 1U);
    bool hash_b = hash_(edge.b_, 1U);

    if (!(hash_a ^ hash_b))
    {
      if (!hash_a)
      {
        if (roots_[edge.a_] == roots_[edge.b_])
        {
          return remove;
        }

        *out = {roots_[edge.a_], roots_[edge.b_]};
      }
      else
      {
        *out = edge;
      }

      return keep;
    }

    if (hash_a)
    {
      if (roots_[edge.b_] == edge.a_)
      {
        *out = {edge.b_, edge.a_}; // To update the root later.
        return remove_update_later;
      }

      if (!hash_b)
      {
        *out = {edge.a_, roots_[edge.b_]};
      }
      else
      {
        *out = edge;
      }
      return keep;
    }

    if (roots_[edge.a_] == edge.b_)
    {
      *out = edge; // To update the root later.
      return remove_update_later;
    }

    *out = {roots_[edge.a_], edge.b_};
    return keep;
  };

  size_t n_compacted = 0;
  select_.iterate(edges_, aux_ + total_compacted_, f, &n_compacted);
  update_later_.push_back(n_compacted);
  total_compacted_ += n_compacted;
}

NAMESPACE_PMT_END