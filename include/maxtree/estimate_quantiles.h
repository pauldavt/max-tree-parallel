#pragma once

#include "../common.h"
#include "graph.h"

NAMESPACE_PMT

template <typename index_t, typename value_t>
class EstimateQuantiles;

template <typename index_t, typename value_t>
void estimate_quantiles(
  Graph<index_t> const& graph,
  value_t const* values,
  size_t n_partitions,
  Quantile<value_t, index_t>* quantiles,
  void* aux1,
  void* aux2) 
{
  EstimateQuantiles<index_t, value_t> eq(graph, values, n_partitions, quantiles, aux1, aux2);
}

template <typename index_t, typename value_t>
class EstimateQuantiles
{
private:
  friend void estimate_quantiles<index_t, value_t>(
  Graph<index_t> const& graph,
  value_t const* values,
  size_t n_partitions,
  Quantile<value_t, index_t>* quantiles,
  void* aux1,
  void* aux2);

  using graph_t = Graph<index_t>;
  using quantile_t = Quantile<value_t, index_t>;
  using uvalue_t = decltype(unsigned_conversion(value_t(0)));
  using sort_index_t = SortValue<index_t>;
  using sort_pair_t = SortPair<uvalue_t, index_t>;
  using edge_t = Edge<index_t>;

  // n_samples = n_samples_factor * n_partitions^2
  static constexpr size_t n_samples_factor = 384U;

  EstimateQuantiles(
    graph_t const& graph,
    value_t const* values,
    size_t n_partitions,
    quantile_t* quantiles,
    void* aux1,
    void* aux2);

  ~EstimateQuantiles();

  sort_pair_t* sort_everything();
  sort_pair_t* create_sorted_sample();

  void determine_sample_n_per_subgraph(size_t total_sample_n_approx);  
  void determine_quantiles(sort_pair_t* uvalue_sorted);

  graph_t const& graph_;
  value_t const* values_;
  size_t n_partitions_;
  quantile_t* quantiles_;
  size_t* n_selected_ = nullptr;
  size_t* offsets_ = nullptr;
  size_t sample_n_;
  
  union
  {
    void* aux1_;
    sort_index_t* index_aux1_;
    sort_pair_t* pair_aux1_;
  };

  union
  {
    void* aux2_;
    sort_index_t* index_aux2_;
    sort_pair_t* pair_aux2_;
  };
};

template <typename index_t, typename value_t>
EstimateQuantiles<index_t, value_t>::EstimateQuantiles(
  graph_t const& graph,
  value_t const* values,
  size_t n_partitions,
  quantile_t* quantiles,
  void* aux1,
  void* aux2) :
  graph_(graph),
  values_(values),
  n_partitions_(n_partitions),
  quantiles_(quantiles),
  aux1_(aux1),
  aux2_(aux2)
{
  size_t n_edges = graph.n_edges();
  check(n_edges > 0);
  size_t n_subgraphs = graph.n_subgraphs();
  n_selected_ = new size_t[2U * n_subgraphs];
  offsets_ = n_selected_ + n_subgraphs;
  size_t total_sample_n_approx = 384U * n_partitions_ * n_partitions_;

  sort_pair_t* sorted;
  if (total_sample_n_approx >= n_edges / 2)
  {
    warn("sorting all edges");
    sorted = sort_everything();
  }
  else
  {
    determine_sample_n_per_subgraph(total_sample_n_approx);      
    sorted = create_sorted_sample();      
  }
  
  determine_quantiles(sorted);
}

template <typename index_t, typename value_t>
EstimateQuantiles<index_t, value_t>::~EstimateQuantiles()
{
  delete[] n_selected_;
}

template <typename index_t, typename value_t>
typename EstimateQuantiles<index_t, value_t>::sort_pair_t*
EstimateQuantiles<index_t, value_t>::sort_everything()
{
  size_t n_edges = graph_.n_edges();
  size_t n_subgraphs = graph_.n_subgraphs();

  sample_n_ = 0;
  for (size_t i = 0; i != n_subgraphs; ++i)
  {
    offsets_[i] = sample_n_;
    size_t n_edges_in_subgraph =
      graph_.global_edge_count(i) + graph_.local_edge_count(i);
    sample_n_ += n_edges_in_subgraph;
  }

  check(sample_n_ == n_edges);

  thread_pool.for_all_blocks(
    n_subgraphs,
    [=](index_t subgraph_nr, thread_nr_t t)
    {
      sort_pair_t* to_sort = pair_aux1_ + offsets_[subgraph_nr];
      edge_t* edges = graph_.subgraph(subgraph_nr);
      size_t n_edges_in_subgraph =
        graph_.global_edge_count(subgraph_nr) +
        graph_.local_edge_count(subgraph_nr);

      for (size_t i = 0; i != n_edges_in_subgraph; ++i)
      {
        to_sort[i] = {pmt::unsigned_conversion(values_[edges[i].a_]), edges[i].a_};
      }
    });

  return radix_sort_parallel(pair_aux2_, pair_aux1_, sample_n_);
}

template <typename index_t, typename value_t>
void EstimateQuantiles<index_t, value_t>::determine_sample_n_per_subgraph(
  size_t total_sample_n_approx)
{
  size_t n_edges = graph_.n_edges();
  size_t n_subgraphs = graph_.n_subgraphs();

  typename rng<size_t>::type r;
  
  sample_n_ = 0;
  for (size_t i = 0; i != n_subgraphs; ++i)
  {
    offsets_[i] = sample_n_;

    size_t n_edges_in_subgraph =
      graph_.global_edge_count(i) + graph_.local_edge_count(i);

    // multiplication may overflow if n_edges_in_subgraph or n_samples is too large
    size_t n_selected = n_edges_in_subgraph * total_sample_n_approx / n_edges;
    size_t remainder = (n_edges_in_subgraph * total_sample_n_approx) % n_edges;

    if (remainder > 0 && random_uint(r, n_edges - 1U) < remainder)
    {
      ++n_selected;
    }

    n_selected_[i] = n_selected;
    sample_n_ += n_selected;
  }
}

template <typename index_t, typename value_t>
typename EstimateQuantiles<index_t, value_t>::sort_pair_t*
EstimateQuantiles<index_t, value_t>::create_sorted_sample()
{
  using rng_t = typename rng<size_t>::type;
  using sort_index_t = SortValue<index_t>;

  size_t n_subgraphs = graph_.n_subgraphs();
  rng_t* rs = new rng_t[thread_pool.max_threads()];
  
  thread_pool.for_all_blocks(
    n_subgraphs,
    [=](index_t subgraph_nr, thread_nr_t t)
    {
      rng_t& r = rs[t];
      sort_index_t* locally_selected = index_aux1_ + offsets_[subgraph_nr];
      size_t n_edges_in_subgraph =
        graph_.global_edge_count(subgraph_nr) +
        graph_.local_edge_count(subgraph_nr);
      edge_t* edges = graph_.subgraph(subgraph_nr);
      size_t max_idx = n_edges_in_subgraph - 1U;

      for (size_t i = 0; i < n_selected_[subgraph_nr]; ++i)
      {
        locally_selected[i].set_value(edges[pmt::random_uint(r, max_idx)].a_);
      }
    });

  // items placed in aux1
  sort_index_t* index_sorted =
    pmt::radix_sort_parallel(index_aux2_, index_aux1_, sample_n_);
  sort_pair_t* to_sort = pair_aux1_;
  sort_pair_t* to_sort_aux = pair_aux2_;

  // items placed in index_sorted
  if (index_sorted == index_aux2_)
  {
    to_sort = pair_aux2_;
    to_sort_aux = pair_aux1_;
  }

  auto const& f_initial = [=](size_t i) ALWAYS_INL_L(sort_pair_t)
    {
      return {
        unsigned_conversion(values_[index_sorted[i].unsigned_value()]),
        index_sorted[i].unsigned_value()};
    };


  sort_pair_t* sorted =
    radix_sort_parallel(to_sort_aux, to_sort, sample_n_, f_initial);
  
  delete[] rs;

  return sorted;
}

template <typename index_t, typename value_t>
void EstimateQuantiles<index_t, value_t>::determine_quantiles(
  sort_pair_t* uvalue_sorted)
{
  quantiles_[0] = {std::numeric_limits<value_t>::min(), index_t(0)};

  for (size_t i = 1; i < n_partitions_; ++i)
  {
    size_t offset = i * sample_n_ / n_partitions_;

    value_t v;
    pmt::undo_unsigned_conversion(uvalue_sorted[offset].unsigned_value(), v);

    quantiles_[i] = {v, uvalue_sorted[offset].data()};
  }    
}

NAMESPACE_PMT_END