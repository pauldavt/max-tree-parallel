#pragma once

#include "../common.h"
#include "../image/image.h"
#include "reduce_edges.h"
#include "graph.h"
#include "../image/image_blocks.h"
#include "../sort/sort_item.h"
#include "rank_set.h"
#include "union_by_rank.h"
#include "estimate_quantiles.h"
#include "graph_partitioning.h"
#include "maxtree_trie.h"

NAMESPACE_PMT

template <typename Primitives>
class Maxtree;

template <typename prim>
void maxtree(Image<prim> const& image, typename prim::index_t* parents)
{
  Maxtree<prim> mp(image, parents);
}

template <typename Primitives>
class Maxtree
{
private:
  using prim = Primitives;
  using index_t = typename Primitives::index_t;
  using image_t = Image<Primitives>;
  using image_blocks_t = ImageBlocks<Primitives>;
  using graph_t = Graph<index_t>;
  using value_t = typename Primitives::value_t;
  using uvalue_t = decltype(pmt::unsigned_conversion(value_t(0)));
  using edge_t = Edge<index_t>;
  using edge_sortpair_t = SortPair<uvalue_t, edge_t>;
  using dim_t = typename image_t::dim_t;
  using rank_set_t = RankSet<index_t>;
  using quantile_t = Quantile<value_t, index_t>;

  friend void maxtree<prim>(
    Image<prim> const& image,
    typename prim::index_t* parents);

  constexpr static size_t n_dimensions = prim::n_dimensions;
  constexpr static size_t n_neighbors = prim::n_neighbors;

  Maxtree(image_t const& image, index_t* parents);
  ~Maxtree();
  void determine_partition_offsets(graph_t* graph);    
  void create_partition_image(graph_t* graph);
  void export_edges(graph_t* graph);
  edge_t* sort_exported_edges(size_t n_edges);
  void union_by_rank_partitions(edge_t* sorted_edges);
  size_t determine_max_edges();

  image_t const& image_;
  index_t* parents_;

  union 
  {
    void* aux1_ = nullptr;
    edge_t* edges_aux1_;
    edge_sortpair_t* sort_edges_aux1_;
    rank_set_t* rank_sets_aux1_;
  };

  union 
  {
    void* aux2_ = nullptr;
    edge_t* edges_aux2_;
    edge_sortpair_t* sort_edges_aux2_;
    rank_set_t* rank_sets_aux2_;
  };
  
  size_t n_;
  image_blocks_t ib_;
  size_t max_partitions_;
  quantile_t* quantiles_ = nullptr;
  uint8_t* partition_img_ = nullptr;
  size_t* partition_offsets_ = nullptr;
  size_t* partition_offsets_per_subgraph_ = nullptr;
};

template <typename prim>
Maxtree<prim>::Maxtree(image_t const& image, index_t* parents) :
  image_(image), parents_(parents), n_(image.dimensions().length()), ib_(image)
{
  if (n_ == 0) return;
  if (n_ == 1)
  {
    parents[0] = 0;
    return;
  }

  size_t n_edges = 0;

  {
    graph_t graph(
      ib_.dimensions().length(),
      image.dimensions().length(),
      determine_max_edges());

    reduce_edges(ib_, parents, &graph);

    n_edges = graph.n_edges();

    if (n_edges == 0)
    {
      return;
    }

    size_t alloc_sz = std::max(n_edges * sizeof(edge_sortpair_t), sizeof(rank_set_t) * n_);

    aux1_ = malloc(alloc_sz);
    aux2_ = malloc(alloc_sz);
    
    check(aux1_ != nullptr && aux2_ != nullptr);

    max_partitions_ = 1U << pmt::log2(hardware_concurrency);
    quantiles_ = new quantile_t[max_partitions_];
    partition_offsets_ = new size_t[max_partitions_ + 1U];
    partition_offsets_per_subgraph_ = new size_t[graph.n_subgraphs() * max_partitions_];
    
    if (max_partitions_ > 1)
    {
      estimate_quantiles(graph, ib_.image().values(), max_partitions_, quantiles_, aux1_, aux2_);       
      create_partition_image(&graph);
      partition_graph(
        ib_,
        &graph,
        parents,
        partition_img_,
        max_partitions_,
        partition_offsets_per_subgraph_,
        edges_aux1_,
        edges_aux2_);
    }
    else
    {
      for (size_t i = 0; i < graph.n_subgraphs(); ++i)
      {
        partition_offsets_per_subgraph_[i] = graph.edge_count(i);
      }
    }
    
    determine_partition_offsets(&graph);

    check(graph.n_edges() == partition_offsets_[max_partitions_]);

    export_edges(&graph);

    n_edges = graph.n_edges();
    // memory used by graph is freed
  }

  edge_t* sorted_edges = sort_exported_edges(n_edges);

#ifdef PMT_DEBUG
  parallel_for_all_blocks(max_partitions_, [=](size_t p, thread_nr_t t)
    {
      index_t begin = partition_offsets_[p];
      index_t end = partition_offsets_[p + 1U];

      for (; begin != end; ++begin)
      {
        edge_t edge = sorted_edges[begin];

        check(partition_img_[edge.a_] == p);
        check(partition_img_[edge.b_] == p);
        check(edge.a_ != edge.b_);
      }
    });
#endif        

  union_by_rank_partitions(sorted_edges);
}

template <typename prim>
Maxtree<prim>::~Maxtree()
{
  delete[] partition_offsets_per_subgraph_;
  delete[] partition_offsets_;
  delete[] partition_img_;
  delete[] quantiles_;
  free(aux1_);
  free(aux2_);
}

template <typename prim>
void Maxtree<prim>::determine_partition_offsets(graph_t* graph)
{
  for (size_t k = 0; k < max_partitions_; ++k)
  {
    partition_offsets_[k] = 0;
  }

  for (size_t i = 0; i < graph->n_subgraphs(); ++i)
  {
    size_t* histo = partition_offsets_per_subgraph_ + i * max_partitions_;

    for (size_t k = 0; k < max_partitions_; ++k)
    {
      size_t tmp = histo[k];
      histo[k] = partition_offsets_[k];
      partition_offsets_[k] += tmp;
    }
  }

  exclusive_sum(partition_offsets_, partition_offsets_ + max_partitions_ + 1U);    

  for (size_t i = 0; i < graph->n_subgraphs(); ++i)
  {
    size_t* histo = partition_offsets_per_subgraph_ + i * max_partitions_;

    for (size_t k = 0; k < max_partitions_; ++k)
    {
      histo[k] += partition_offsets_[k];
    }
  }    
}
  

template <typename prim>
void Maxtree<prim>::create_partition_image(graph_t* graph)
{
  {
    partition_img_ = new partition_t[n_];    
    value_t const* values = ib_.image().values();

    thread_pool.for_all(n_, [=](index_t i, thread_nr_t t) ALWAYS_INLINE {     
      partition_img_[i] = quantile_t::determine_partition(values[i], i, quantiles_, max_partitions_);
    });
  }
}

template <typename prim>
void Maxtree<prim>::export_edges(graph_t* graph)
{
  value_t const* values = image_.values();
  size_t n_subgraphs = graph->n_subgraphs();

  if (max_partitions_ == 1)
  {
    thread_pool.for_all_blocks(n_subgraphs, [=](size_t subgraph_nr, thread_nr_t t) {
      edge_t* edges_begin = graph->subgraph(subgraph_nr);
      edge_t* edges_end = edges_begin + graph->edge_count(subgraph_nr);
      size_t& out_offsets = partition_offsets_per_subgraph_[subgraph_nr];

      for (; edges_begin != edges_end; ++edges_begin)      
      {
        edge_t const& edge = *edges_begin;          
        sort_edges_aux2_[out_offsets++] = {unsigned_conversion(values[edge.a_]), edge};
      }
    });      

    return;
  }

  thread_pool.for_all_blocks(n_subgraphs, [=](size_t subgraph_nr, thread_nr_t t) {
    edge_t* edges_begin = graph->subgraph(subgraph_nr);
    edge_t* edges_end = edges_begin + graph->edge_count(subgraph_nr);
    size_t* out_offsets = partition_offsets_per_subgraph_ + subgraph_nr * max_partitions_;

    for (; edges_begin != edges_end; ++edges_begin)      
    {
      edge_t const& edge = *edges_begin;
      
      partition_t p = partition_img_[edge.a_];

      sort_edges_aux2_[out_offsets[p]++] = {unsigned_conversion(values[edge.a_]), edge};
    }
  });
}

template <typename prim>
typename Maxtree<prim>::edge_t *
Maxtree<prim>::sort_exported_edges(size_t n_edges)
{
  edge_t* sorted_edges =
    radix_sort_n_digits<uvalue_t>() & 1 ? edges_aux1_ : edges_aux2_;

  auto const& f_initial =
    [=](size_t i) ALWAYS_INL_L(edge_sortpair_t)
    {
      return sort_edges_aux2_[i];
    };      

  auto const& f_out =
    [=](edge_t& out, edge_sortpair_t const& item) ALWAYS_INLINE
    {
      out = item.data();
    };      

  radix_sort_parallel(
    sorted_edges,
    sort_edges_aux1_,
    sort_edges_aux2_,
    n_edges,
    0,
    sizeof(uvalue_t) * CHAR_BIT,
    f_initial,
    f_out);

#ifdef PMT_DEBUG
  value_t const* values = image_.values();

  for (size_t i = 1; i < n_edges; ++i)
  {
    check(values[sorted_edges[i].a_] <= values[sorted_edges[i].b_]);
    check(values[sorted_edges[i].a_] >= values[sorted_edges[i - 1].a_]);
  }
#endif        

  return sorted_edges;
}

template <typename prim>
void Maxtree<prim>::union_by_rank_partitions(edge_t* sorted_edges)
{
  rank_set_t* rank_sets =
    sorted_edges == edges_aux1_ ? rank_sets_aux2_ : rank_sets_aux1_;

  thread_pool.for_all(
    n_,
    [=](index_t i, thread_nr_t thread_nr) ALWAYS_INL_L(void)
    {
      rank_sets[i].reset(i);
    });    

  thread_pool.for_all_blocks(max_partitions_, [=](size_t p, thread_nr_t t)
    {
      index_t begin = partition_offsets_[p];
      index_t end = partition_offsets_[p + 1U];

      maxtree_union_by_rank(sorted_edges, begin, end, rank_sets, parents_);
    });  
}

template <typename prim>
size_t Maxtree<prim>::determine_max_edges()
{
  dim_t const& dims = ib_.image().dimensions();
  dim_t const& grid_dims = ib_.dimensions();

  size_t size = dims.length();
  size_t max_edges = size;

  for (dim_idx_t d = 0; d != n_dimensions; ++d)
  {
    size_t n_edges = size / dims[d] * (grid_dims[d] - size_t(1));

    max_edges += n_edges;
  }

  if (n_dimensions == 2 && n_neighbors == 8)
  {
    // NW & NE horizontal edges
    max_edges += 2U * (grid_dims[1] - size_t(1)) * (dims[0] - size_t(1));

    // NW & SW vertical edges
    max_edges += 2U * (grid_dims[0] - size_t(1)) * (dims[1] - size_t(1));

    // double edges in corners NW and NE
    max_edges -= 2U * (grid_dims[0] - size_t(1)) * (grid_dims[1] - size_t(1));
  }

  return max_edges;
}

NAMESPACE_PMT_END