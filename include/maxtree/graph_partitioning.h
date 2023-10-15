#pragma once

#include "../common.h"
#include "graph.h"
#include "../misc/bits.h"
#include "../misc/bit_array.h"
#include "../image/image_blocks.h"
#include "connected_components.h"

NAMESPACE_PMT

template <typename Primitives>
class GraphPartitioning;

template <typename prim>
void partition_graph(
  ImageBlocks<prim> const& ib,
  Graph<typename prim::index_t>* graph,
  typename prim::index_t* parents,
  partition_t* img,
  size_t max_partitions,
  size_t* total_partition_counts,
  Edge<typename prim::index_t>* aux1,
  Edge<typename prim::index_t>* aux2)
{
  GraphPartitioning<prim> gp(
    ib,
    graph,
    parents,
    img,
    max_partitions,
    total_partition_counts,
    aux1,
    aux2);
}

template <typename Primitives>
class GraphPartitioning
{
private:
  using prim = Primitives;
  using index_t = typename Primitives::index_t;
  using value_t = typename Primitives::value_t;
  using graph_t = Graph<index_t>;
  using edge_t = Edge<index_t>;
  using image_block_t = ImageBlock<Primitives>;
  using image_blocks_t = ImageBlocks<Primitives>;
  using vec_t = Coordinate<Primitives>;

  friend void partition_graph<prim>(
  ImageBlocks<prim> const& ib,
  Graph<typename prim::index_t>* graph,
  typename prim::index_t* parents,
  partition_t* img,
  size_t max_partitions,
  size_t* total_partition_counts,
  Edge<typename prim::index_t>* aux1,
  Edge<typename prim::index_t>* aux2);


  GraphPartitioning(
    image_blocks_t const& ib,
    graph_t* graph,
    index_t* parents,
    partition_t* img,
    size_t max_partitions,
    size_t* partition_counts_per_subgraph,
    edge_t* aux1,
    edge_t* aux2);

  ~GraphPartitioning();

  bool completed() const { return completed_; }
  void partition();
  // from aux1_ -> aux2_
  void compact_edges01_and_edges11();
  void export_local_edges(thread_nr_t t, size_t subgraph_nr, edge_t** local_edges_end_ptr, edge_t** out_edges01_begin_ptr);
  void export_global_edges(thread_nr_t t, size_t subgraph_nr, edge_t* local_edges_end, edge_t* out_edges01_begin);
  // copy from edges -> aux1_
  void export_edges01_and_edges11();

  ALWAYS_INLINE_F bool bit(index_t i) const
  {
    return (img_[i] >> partition_t(msb_)) & partition_t(1);
  }

  image_blocks_t const& ib_;
  graph_t& graph_;
  index_t* parents_;
  partition_t* img_;
  edge_t* aux1_;
  edge_t* aux2_;
  size_t max_partitions_;
  size_t msb_;
  size_t* aux_subgraph_offsets_;
  size_t* edges01_counts_;
  size_t* edges11_counts_;
  size_t* edges01_offsets_;
  size_t* edges11_offsets_;
  index_t* roots_ = nullptr;
  size_t* partition_counts_per_subgraph_;
  bool completed_;

#ifdef PMT_DEBUG
  size_t* checksums_;
#endif  
};

template <typename prim>
GraphPartitioning<prim>::GraphPartitioning(
  image_blocks_t const& ib,
  graph_t* graph,
  index_t* parents,
  partition_t* img,
  size_t max_partitions,
  size_t* partition_counts_per_subgraph,
  edge_t* aux1,
  edge_t* aux2) :
  ib_(ib),
  graph_(*graph),
  parents_(parents),
  img_(img),
  aux1_(aux1),
  aux2_(aux2),
  max_partitions_(max_partitions),
  msb_(pmt::log2(max_partitions) - 1U),
  partition_counts_per_subgraph_(partition_counts_per_subgraph)
{
  completed_ = max_partitions <= 1;    

  if (completed_) return;

  size_t n_subgraphs = graph_.n_subgraphs();

  aux_subgraph_offsets_ = new size_t[5U * (n_subgraphs + 1U)];
  edges01_counts_ = aux_subgraph_offsets_ + n_subgraphs + 1U;
  edges11_counts_ = edges01_counts_ + n_subgraphs + 1U;
  edges01_offsets_ = edges11_counts_ + n_subgraphs + 1U;
  edges11_offsets_ = edges01_offsets_ + n_subgraphs + 1U;
  roots_ = new index_t[graph_.max_nodes()];

#ifdef PMT_DEBUG
  checksums_ = new size_t[n_subgraphs + 1U];
#endif

  size_t offset = 0;
  for (size_t i = 0; i < n_subgraphs; ++i)
  {
    aux_subgraph_offsets_[i] = offset;
    offset += graph_.edge_count(i);
  }

  aux_subgraph_offsets_[n_subgraphs] = offset;
  debug(offset == graph_.n_edges());

  while (!completed())
  {
    partition();
  }

  // determine_partition_counts(total_partition_counts);
}

template <typename prim>
GraphPartitioning<prim>::~GraphPartitioning()
{
  delete[] roots_;
  delete[] aux_subgraph_offsets_;

#ifdef PMT_DEBUG
  delete[] checksums_;
#endif    
}

template <typename prim>
void GraphPartitioning<prim>::partition()
{
  check(!completed_);

  size_t n_subgraphs = graph_.n_subgraphs();


#ifdef PMT_DEBUG
  for (size_t i = 0; i < n_subgraphs; ++i)
  {
    checksums_[i] = 0;
  }
#endif

    // from edges -> aux1_ -> (compact) aux2_
  export_edges01_and_edges11();

  edge_t* edges11 = aux2_ + edges01_offsets_[n_subgraphs];
  size_t n_edges11 = edges11_offsets_[n_subgraphs];
  {
    if (n_edges11 > 0)
    {
      connected_components(edges11, n_edges11, aux1_ + edges01_offsets_[n_subgraphs], ib_.image().values(), roots_);
    }
  }

  {
    thread_pool.for_all(edges01_offsets_[n_subgraphs], [=](size_t i, thread_nr_t) {
      aux1_[i] = aux2_[i];
    });
  }
  
  size_t n_edges01 = edges01_offsets_[n_subgraphs];
  value_t const* values =  ib_.image().values();
  {
    ItemBlocks select(n_edges01);

    while (select.length() > 0)
    {
      select.select([=](size_t i, size_t o) ALWAYS_INL_L(bool) {
        edge_t const& edge = aux2_[i];

        index_t cc_root = roots_[edge.b_];
        index_t parent_cc = parents_[cc_root];

        if (parent_cc == cc_root || values[edge.a_] > values[parent_cc] ||
          (values[edge.a_] == values[parent_cc] && edge.a_ > parent_cc))
        {
          parents_[cc_root] = edge.a_;
          aux2_[o] = edge;
          return true;
        }

        return false;
      });
    }
  }

  {
    thread_pool.for_all_blocks(graph_.n_subgraphs(), [=](index_t subgraph_nr, thread_nr_t t)
    {
      size_t n_new = edges01_counts_[subgraph_nr];
      edge_t* edges = aux1_ + edges01_offsets_[subgraph_nr];
      edge_t* out = graph_.subgraph(subgraph_nr) + graph_.edge_count(subgraph_nr);
      size_t* part_counts = partition_counts_per_subgraph_ + subgraph_nr * max_partitions_;

      size_t ctr = 0;
      for (size_t i = 0; i < n_new; ++i)
      {
        edge_t edge = edges[i];
        index_t parent = parents_[roots_[edge.b_]];
        
        debug(roots_[roots_[edge.b_]] == roots_[edge.b_]);
        debug(parent != roots_[edge.b_]);
        debug(!bit(edge.a_));
        debug(!bit(parent));
        debug(values[parent] <= values[edge.b_]);
        debug(values[edge.a_] <= values[parent]);

        if (edge.a_ != parent)
        {
          if (msb_ == 0)
          {
            ++part_counts[img_[edge.a_]];
          }

          out[ctr++] = {edge.a_, parent};
        }
      }

      graph_.set_global_edge_count(subgraph_nr, graph_.global_edge_count(subgraph_nr) + ctr);      
    });
  }

  if (msb_ == 0)
  {
    completed_ = true;

    graph_.determine_n_edges();

    return;
  }

  --msb_;
}

// from aux1_ -> aux2_
template <typename prim>
void GraphPartitioning<prim>::compact_edges01_and_edges11()
{
  //pmt::Timer t("compact");
  size_t n_subgraphs = graph_.n_subgraphs();

  edge_t* out_edges01 = aux2_;
  edge_t* out_edges11 = aux2_ + edges01_offsets_[n_subgraphs];

  thread_pool.for_all_blocks(n_subgraphs, [=](index_t subgraph_nr, thread_nr_t t) {
    edge_t* begin = aux1_ + aux_subgraph_offsets_[subgraph_nr];
    edge_t* end = begin + edges01_counts_[subgraph_nr];
    edge_t* out = out_edges01 + edges01_offsets_[subgraph_nr];

#ifdef PMT_DEBUG
    size_t checksum = 0;
#endif

    for (; begin != end; ++begin)
    {
      edge_t edge = *begin;
#ifdef PMT_DEBUG
      check(img_[edge.a_] <= img_[edge.b_]);
      check(edge.a_ < graph_.max_nodes() && edge.b_ < graph_.max_nodes());
      check(roots_[roots_[edge.b_]] == roots_[edge.b_]);
      check(!bit(edge.a_) && bit(edge.b_));
      checksum += edge.a_ + edge.b_;
#endif
      

      *out++ = {edge.a_, roots_[edge.b_]};
    }

    end = aux1_ + aux_subgraph_offsets_[subgraph_nr + 1U];
    begin = end - edges11_counts_[subgraph_nr];
    out = out_edges11 + edges11_offsets_[subgraph_nr];
    
    for (; end-- != begin;)
    {
      edge_t edge = *end;
#ifdef PMT_DEBUG
      check(img_[edge.a_] <= img_[edge.b_]);
      check(edge.a_ < graph_.max_nodes() && edge.b_ < graph_.max_nodes());
      checksum += edge.a_ + edge.b_;
      check(bit(edge.a_) && bit(edge.b_));
      check(roots_[roots_[edge.a_]] == roots_[edge.a_]);
      check(roots_[roots_[edge.b_]] == roots_[edge.b_]);
#endif

      *out++ = {roots_[edge.a_], roots_[edge.b_]};
    }

#ifdef PMT_DEBUG
    if (checksum != checksums_[subgraph_nr]) {
      out(checksum << " " << checksums_[subgraph_nr]);
    };
    check(checksum == checksums_[subgraph_nr]);
#endif      
  });
}

template <typename prim>
void GraphPartitioning<prim>::export_local_edges(thread_nr_t t, size_t subgraph_nr, edge_t** local_edges_end_ptr, edge_t** out_edges01_begin_ptr)
{
  edge_t* out_edges01_begin = aux1_ + aux_subgraph_offsets_[subgraph_nr];
  edge_t* local_edges = graph_.subgraph(subgraph_nr);
  edge_t* const local_edges_end = local_edges + graph_.local_edge_count(subgraph_nr);
  edge_t* compact = local_edges;
  size_t* part_counts = partition_counts_per_subgraph_ + subgraph_nr * max_partitions_;

  for (; local_edges != local_edges_end; ++local_edges)
  {
    edge_t const edge = *local_edges;

    debug(img_[edge.a_] <= img_[edge.b_]);

    bool bit_a = bit(edge.a_);
    bool bit_b = bit(edge.b_);
    bool is_edge01 = !bit_a && bit_b;
    bool is_edge11 = bit_a;

    debug(!is_edge11 || bit_b);

    if (is_edge01)
    {
#ifdef PMT_DEBUG
      checksums_[subgraph_nr] += edge.a_ + edge.b_;
#endif

      *out_edges01_begin++ = edge;
      continue;
    }
    else if (is_edge11)
    {
      roots_[edge.b_] = roots_[edge.a_];
      debug(roots_[edge.a_] == roots_[roots_[edge.a_]]);
    }

    if (msb_ == 0)
    {
      ++part_counts[img_[edge.a_]];
    }

    *compact++ = edge;
  } 

  graph_.set_local_edge_count(subgraph_nr, compact - graph_.subgraph(subgraph_nr)); 

  *local_edges_end_ptr = local_edges_end;
  *out_edges01_begin_ptr = out_edges01_begin;
}

template <typename prim>
void GraphPartitioning<prim>::export_global_edges(thread_nr_t t, size_t subgraph_nr, edge_t* local_edges_end, edge_t* out_edges01_begin)
{
  edge_t* out_edges11_end = aux1_ + aux_subgraph_offsets_[subgraph_nr + 1U];
  edge_t* global_edges = local_edges_end;
  edge_t* compact = graph_.subgraph(subgraph_nr) + graph_.local_edge_count(subgraph_nr);
  edge_t* const global_edges_end = global_edges + graph_.global_edge_count(subgraph_nr);
  size_t* part_counts = partition_counts_per_subgraph_ + subgraph_nr * max_partitions_;

  for (; global_edges != global_edges_end; ++global_edges)
  {
    edge_t edge = *global_edges;

    debug(img_[edge.a_] <= img_[edge.b_]);

    bool bit_a = bit(edge.a_);
    bool bit_b = bit(edge.b_);
    bool is_edge01 = !bit_a && bit_b;
    bool is_edge11 = bit_a;

    debug(!is_edge11 || bit_b);

#ifdef PMT_DEBUG
    if (is_edge01 || is_edge11)
    {
      checksums_[subgraph_nr] += edge.a_ + edge.b_;
    }
#endif      

    if (is_edge01)
    {
      *out_edges01_begin++ = edge;
      continue;
    }       
    else if (is_edge11)
    {
      *(--out_edges11_end) = edge;
    }

    if (msb_ == 0)
    {
      ++part_counts[img_[edge.a_]];
    }

    *compact++ = edge;
  }

  debug(out_edges01_begin <= out_edges11_end);

  graph_.set_global_edge_count(subgraph_nr, compact - graph_.subgraph(subgraph_nr) - graph_.local_edge_count(subgraph_nr));

  size_t edges01_count = out_edges01_begin - aux1_ - aux_subgraph_offsets_[subgraph_nr];
  edges01_counts_[subgraph_nr] = edges01_count;
  edges01_offsets_[subgraph_nr] = edges01_count;
    
  size_t edges11_count = aux1_ + aux_subgraph_offsets_[subgraph_nr + 1U] - out_edges11_end;
  edges11_counts_[subgraph_nr] = edges11_count;
  edges11_offsets_[subgraph_nr] = edges11_count;
}

// copy from edges -> aux1_
template <typename prim>
void GraphPartitioning<prim>::export_edges01_and_edges11()
{
  size_t n_subgraphs = graph_.n_subgraphs();

  thread_pool.for_all_blocks<prim>(ib_.dimensions(), [=](vec_t const& block_loc, thread_nr_t t) {
    image_block_t block(ib_, block_loc);
    size_t subgraph_nr = block.block_nr();
    
    using block_index_t = typename image_block_t::block_index_t;
    block.apply([=](index_t global_index, block_index_t block_index) ALWAYS_INLINE {
      roots_[global_index] = global_index;
    });

    size_t* part_counts = partition_counts_per_subgraph_ + subgraph_nr * max_partitions_;      

    if (msb_ == 0)
    {
      for (size_t i = 0; i < max_partitions_; ++i)
      {
        part_counts[i] = 0;
      }
    }
    
    edge_t* local_edges_end;
    edge_t* out_edges01_begin;   
    export_local_edges(t, subgraph_nr, &local_edges_end, &out_edges01_begin);
    export_global_edges(t, subgraph_nr, local_edges_end, out_edges01_begin);
  });

  exclusive_sum(edges01_offsets_, edges01_offsets_ + n_subgraphs + 1U);
  exclusive_sum(edges11_offsets_, edges11_offsets_ + n_subgraphs + 1U);

  compact_edges01_and_edges11();
}

NAMESPACE_PMT_END