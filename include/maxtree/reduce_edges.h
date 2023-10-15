#pragma once

#include "../common.h"
#include "../misc/edge.h"
#include "../image/image_blocks.h"
#include "graph.h"
#include "../sort/radix_sort_parallel.h"
#include "../sort/radix_sort_seq.h"
#include "../image/image_blocks.h"
#include "../misc/exclusive_sum.h"
#include "../misc/trie_queue.h"
#include "../maxtree/maxtree_trie.h"

NAMESPACE_PMT

template <typename Primitives>
class ReduceEdges;

template <typename prim>
void reduce_edges(
  ImageBlocks<prim> const& ib,
  typename prim::index_t* parents,
  Graph<typename prim::index_t>* graph)
{
  ReduceEdges<prim> re(ib, parents, graph);
}

template <typename Primitives>
class ReduceEdges
{
private:
  using prim = Primitives;
  using index_t = typename prim::index_t;
  using value_t = typename prim::value_t;
  using uvalue_t = decltype(pmt::unsigned_conversion(value_t(0)));
  using graph_t = Graph<index_t>;
  using dim_t = pmt::Dimensions<prim::n_dimensions>;
  using vec_t = pmt::Coordinate<prim>;
  using edge_t = Edge<index_t>;

  using image_blocks_t = ImageBlocks<prim>;
  using image_block_t = ImageBlock<prim>;
  using block_t = pmt::Block<prim::n_dimensions>;
  using block_index_t = typename block_t::block_index_t;
  using block_vec_t = typename image_block_t::block_vec_t;
  using block_loc_t = typename image_block_t::block_loc_t;
  using sort_pair_t = pmt::SortPair<uvalue_t, block_index_t>;
  using queue_t = pmt::TrieQueue<block_index_t>;
  using visited_t = pmt::BitArray;

  friend void reduce_edges<prim>(
    ImageBlocks<prim> const& ib,
    typename prim::index_t* parents,
    Graph<typename prim::index_t>* graph);

  constexpr static size_t n_dimensions = prim::n_dimensions;
  constexpr static size_t n_neighbors = prim::n_neighbors;

  struct thread_data
  {

    thread_data() :
      sort_space_1(&sort_space_2[max_items_per_block]),
      visited(max_items_per_block),
      queue(max_items_per_block - index_t(1), false)
    {
    }

    static constexpr size_t max_items_per_block = block_t::max_length;
    static constexpr size_t n_aux_sort_arrays = radix_sort_n_digits<uvalue_t>() == 1 ? 1 : 2;

    union
    {
      sort_pair_t sort_space_2[n_aux_sort_arrays * max_items_per_block];      
      index_t local_to_global_index[max_items_per_block];
    };
    
    sort_pair_t* sort_space_1;
    block_index_t parents[max_items_per_block];
    block_index_t rank_to_index[max_items_per_block];
    visited_t visited;
    queue_t queue;
  };

  ReduceEdges(image_blocks_t const& ib, index_t* parents, graph_t* graph);
  void determine_local_edges(image_block_t const& block, vec_t const& block_loc, index_t block_nr, thread_data* data);
  void iterate_blocks_parallel();
  void determine_edge_offsets();
  void determine_edge_offsets_2d_8n();
  void add_edge(edge_t* out, index_t current, index_t neighbor);

  size_t add_global_edges(
    dim_t const& block_dims,
    dim_idx_t d,
    dim_idx_t d_exclude,
    index_t index_offset,
    size_t edge_offset,
    index_t* skip);

  void add_global_edges(image_block_t const& block, vec_t const& block_loc, index_t block_nr);
  void add_global_edges_2d_8n(image_block_t const& block, vec_t const& block_loc, index_t block_nr);

  image_blocks_t const& ib_;
  index_t* parents_;
  graph_t& graph_;
};

template <typename prim>
ReduceEdges<prim>::ReduceEdges(image_blocks_t const& ib, index_t* parents, graph_t* graph) :
  ib_(ib),
  parents_(parents),
  graph_(*graph)
{
  if (n_dimensions == 2 && n_neighbors == 8)
  {
    determine_edge_offsets_2d_8n();
  }
  else
  {
    determine_edge_offsets();
  }
  
  iterate_blocks_parallel();

  graph_.determine_n_edges();
}

template <typename prim>
void ReduceEdges<prim>::determine_local_edges(image_block_t const& block, vec_t const& block_loc, index_t block_nr, thread_data* data)
{
  size_t n_items_in_block = block.dimensions().length();
  value_t const* vals = ib_.image().values();
  sort_pair_t* sort_space_2 = data->sort_space_2;
  block_index_t* rank_to_index = data->rank_to_index;
  block_index_t* index_to_rank = data->parents;
  index_t* local_to_global_index = data->local_to_global_index;

  // init sort pairs, to make index_to_rank and rank_to_index mappings
  block.apply([=](index_t global_index, block_index_t local_index) ALWAYS_INLINE {
    parents_[global_index] = global_index; // init parent
    sort_space_2[local_index] = {unsigned_conversion(vals[global_index]), local_index};
  });          


  auto const& f_initial_item = [=](block_index_t i) ALWAYS_INL_L(sort_pair_t)
  {
    return sort_space_2[i];
  };

  auto const& f_out = [=](block_index_t& out, sort_pair_t const& item) ALWAYS_INLINE
  {
    out = item.data();
  };      

  pmt::radix_sort_seq(
    rank_to_index,
    n_items_in_block,
    data->sort_space_1,
    data->sort_space_2,
    f_initial_item,
    f_out);

  // index_to_rank and rank_to_index mappings
  for (size_t rank = 0; rank != n_items_in_block; ++rank)
  {
    index_to_rank[rank_to_index[rank]] = rank;
  }

  using block_image_t = typename image<
    block_index_t,
    block_index_t,
    n_dimensions,
    prim::n_neighbors>::type;

  block_image_t block_img(index_to_rank, block.dimensions()); 

  // create a max-tree of the block
  pmt::maxtree_trie(block_img, data->parents, rank_to_index, &data->visited, &data->queue);

  // change to a boundary tree

  data->visited.clear();
  block.set_boundaries(&data->visited);

  for (size_t i = n_items_in_block; i--;)
  {
    block_index_t k = rank_to_index[i];

    if (data->visited.is_set(k))
    {
      data->visited.set(data->parents[k]);
    }
  }

  block.apply([=](index_t global_index, block_index_t local_index) ALWAYS_INLINE {
    local_to_global_index[local_index] = global_index;    
  });

  edge_t* out = graph_.subgraph(block_nr);
        
  for (size_t i = 0; i != n_items_in_block; ++i)
  {
    block_index_t k = rank_to_index[i];
    if (data->visited.is_set(k) && data->parents[k] != k)
    {
      // local edges, actual parent of k may be outside of the block
      *out++ = {local_to_global_index[data->parents[k]], local_to_global_index[k]};
    }
    else
    {
      // these global parents can already be set
      parents_[local_to_global_index[k]] = local_to_global_index[data->parents[k]];
    }
  }          

  graph_.set_local_edge_count(block_nr, out - graph_.subgraph(block_nr));
}

template <typename prim>
void ReduceEdges<prim>::iterate_blocks_parallel()
{
  thread_data* ts = new thread_data[thread_pool.max_threads()];

  thread_pool.for_all_blocks<prim>(ib_.dimensions(), [=](vec_t const& block_loc, thread_nr_t thread_nr) {
    //printf("thread %ld doing block %d %d\n", thread_nr, block_loc[1], block_loc[0]);
    image_block_t block(ib_, block_loc);
    index_t block_nr = block.block_nr();
    thread_data& data = ts[thread_nr];

    determine_local_edges(block, block_loc, block_nr, &data);

    if (n_dimensions == 2 && n_neighbors == 8)
    {
      add_global_edges_2d_8n(block, block_loc, block_nr);
    }
    else
    {
      add_global_edges(block, block_loc, block_nr);
    }
  });

  graph_.determine_n_edges();

//    out(graph_.n_edges());

  delete[] ts;
}

template <typename prim>
void ReduceEdges<prim>::determine_edge_offsets()
{
  vec_t block_loc;
  block_loc.init_zeros();

  dim_t grid_dims = ib_.dimensions();
  size_t offset = 0;
  size_t n_subgraphs = grid_dims.length();

  for (size_t i = 0; i < n_subgraphs; ++i)
  {
    dim_t block_dims;
    image_block_t::determine_dimensions(ib_.image(), block_loc, &block_dims);
    size_t block_size = block_dims.length();
    size_t n_edges = block_size; // local edges, connecting nodes in the block

    for (dim_idx_t d = 0; d != n_dimensions; ++d)
    {
      if (block_loc[d] == 0) continue;

      n_edges += block_size / block_dims[d]; // global edges, connecting to previous blocks
    }

    graph_.set_subgraph_offset(i, offset);
    offset += n_edges;       
    block_loc.inc_index(grid_dims);      
  }

  graph_.set_subgraph_offset(n_subgraphs, offset);

  check(offset == graph_.max_edges());
}

template <typename prim>
void ReduceEdges<prim>::determine_edge_offsets_2d_8n()
{
  vec_t block_loc;
  block_loc.init_zeros();

  dim_t grid_dims = ib_.dimensions();
  size_t offset = 0;
  size_t n_subgraphs = grid_dims.length();

  for (size_t i = 0; i < n_subgraphs; ++i)
  {
    dim_t block_dims;
    image_block_t::determine_dimensions(ib_.image(), block_loc, &block_dims);
    size_t block_size = block_dims.length();
    size_t n_edges = block_size; // local edges, connecting nodes in the block
    size_t block_width = block_dims[0];
    size_t block_height = block_dims[1];

    if (block_loc[1] > 0)
    {
      // first pixel
      if (block_loc[0] > 0)
      {
        ++n_edges; // NW
      }      

      ++n_edges; // N

      if (block_width > 1)
      {
        n_edges += 1; // NE first pixel on first col
        n_edges += 3U * (block_width - 2U); // NW, N and NE
        n_edges += 2; // N and NW last pixel on first col
      }      

      if (block_loc[0] < grid_dims[0] - size_t(1))
      {
        ++n_edges; // NE last pixel on first col
      }      
    }

    if (block_loc[0] > 0)
    {
      // first pixel
      /*
        // already added
      if (block_loc[1] > 0)
      {
        ++n_edges; // NW
      } */

      ++n_edges; // W

      if (block_height > 1)
      {
        n_edges += 1; // SW first pixel on first row
        n_edges += 3U * (block_height - 2U); // NW, W and SW
        n_edges += 2; // NW and W first pixel on last row
      }
    }    

    graph_.set_subgraph_offset(i, offset);
    offset += n_edges;       
    block_loc.inc_index(grid_dims);      
  }

  graph_.set_subgraph_offset(n_subgraphs, offset);

  //out(offset << " " << graph_.max_edges());
  check(offset == graph_.max_edges());
}

template <typename prim>
ALWAYS_INLINE_F void
ReduceEdges<prim>::add_edge(edge_t* out, index_t current, index_t neighbor)
{
  value_t const* values = ib_.image().values();

  if (values[neighbor] > values[current])
    *out = {current, neighbor};
  else
    *out = {neighbor, current};    
}

template <typename prim>
size_t ReduceEdges<prim>::add_global_edges(
  dim_t const& block_dims,
  dim_idx_t d,
  dim_idx_t d_exclude,
  index_t index_offset,
  size_t edge_offset,
  index_t* skip)
{
  value_t const* values = ib_.image().values();
  edge_t* edges = graph_.edges();

  if (d == 0)
  {
    for (size_t i = 0; i < block_dims[0]; ++i)
    {
      index_t neighbor = index_offset - skip[d_exclude];        

      //printf("adding %d %d\n", index_offset, neighbor);

      if (values[neighbor] > values[index_offset])
        edges[edge_offset] = {index_offset, neighbor};
      else
        edges[edge_offset] = {neighbor, index_offset};

      ++index_offset;
      ++edge_offset;
    }

    return block_dims[0];
  }
  if (d == d_exclude)
  {      
    return add_global_edges(block_dims, d - dim_idx_t(1), d_exclude, index_offset, edge_offset, skip);
  }
  if (d == 1 && d_exclude == 0)
  {
    for (size_t i = 0; i < block_dims[1]; ++i)
    {
      index_t neighbor = index_offset - index_t(1);

      if (values[neighbor] > values[index_offset])
        edges[edge_offset] = {index_offset, neighbor};
      else
        edges[edge_offset] = {neighbor, index_offset};

      index_offset += skip[1];
      ++edge_offset;
    }

    return block_dims[1];
  }

  size_t ctr = 0;

  for (size_t i = 0; i < block_dims[d]; ++i)
  {
    size_t n_added = add_global_edges(block_dims, d - dim_idx_t(1), d_exclude, index_offset, edge_offset, skip);
    ctr += n_added;
    edge_offset += n_added;
    index_offset += skip[d];
  }

  return ctr;
}

template <typename prim>
void ReduceEdges<prim>::add_global_edges(image_block_t const& block, vec_t const& block_loc, index_t block_nr)
{
  size_t ctr = 0;
  dim_t const& dims = ib_.image().dimensions();
  index_t skip[n_dimensions];

  skip[0] = 1;

  for (dim_idx_t d = 1; d != n_dimensions; ++d)
  {
    skip[d] = skip[d - dim_idx_t(1)] * dims[d - dim_idx_t(1)];
  }

  size_t edge_offset = graph_.subgraph_offset(block_nr) + graph_.local_edge_count(block_nr);

  for (dim_idx_t d = 0; d != n_dimensions; ++d)
  {
    if (block_loc[d] == 0) continue;

    size_t n_added = add_global_edges(block.dimensions(), n_dimensions - dim_idx_t(1), d, block.global_offset(), edge_offset, skip);
    edge_offset += n_added;
    ctr += n_added;
  }

  debug(edge_offset - (graph_.subgraph_offset(block_nr) + graph_.local_edge_count(block_nr)) == ctr);
  debug(graph_.subgraph_offset(block_nr) + ctr <= graph_.subgraph_offset(block_nr + 1U));

  graph_.set_global_edge_count(block_nr, ctr);
}

template <typename prim>
void ReduceEdges<prim>::add_global_edges_2d_8n(image_block_t const& block, vec_t const& block_loc, index_t block_nr)
{
  debug(n_dimensions == 2 && n_neighbors == 8);

  dim_t const& block_dims = block.dimensions();
  dim_t const& dims = ib_.image().dimensions();
  edge_t* out = graph_.subgraph(block_nr) + graph_.local_edge_count(block_nr);

  size_t img_width = dims[0];
  size_t block_width = block_dims[0];
  size_t block_height = block_dims[1];

  // ONLY edges linking to nodes with lesser index

  if (block_loc[1] > 0)
  {
    // north

    index_t current = block.global_offset();

    if (block_loc[0] > 0)
    {
      add_edge(out++, current, current - img_width - index_t(1));        
    }      
    add_edge(out++, current, current - img_width);
    if (block_width > 1)
    {
      add_edge(out++, current, current - img_width + index_t(1));
    }      

    for (size_t w = 1; w < block_width - size_t(1); ++w)
    {      
      ++current;
      add_edge(out++, current, current - img_width - index_t(1));
      add_edge(out++, current, current - img_width);
      add_edge(out++, current, current - img_width + index_t(1));
    }

    if (block_width > 1)
    {
      ++current;
      add_edge(out++, current, current - img_width - index_t(1));        
      add_edge(out++, current, current - img_width);
    }
    if (block_loc[0] < ib_.dimensions()[0] - size_t(1))
    {
      add_edge(out++, current, current - img_width + index_t(1));
    }      
  }

  if (block_loc[0] > 0)
  {
    // west

    index_t current = block.global_offset();

    add_edge(out++, current, current - index_t(1));
    if (block_height > 1)
    {
      add_edge(out++, current, current + img_width - index_t(1));
    }      

    for (size_t h = 1; h < block_height - size_t(1); ++h)
    {      
      current += img_width;
      add_edge(out++, current, current - img_width - index_t(1));
      add_edge(out++, current, current - index_t(1));
      add_edge(out++, current, current + img_width - index_t(1));
    }

    if (block_height > 1)
    {
      current += img_width;
      add_edge(out++, current, current - img_width - index_t(1));        
      add_edge(out++, current, current - index_t(1));
    }
  }

  size_t ctr = out - (graph_.subgraph(block_nr) + graph_.local_edge_count(block_nr));

  graph_.set_global_edge_count(block_nr, ctr);
}  

NAMESPACE_PMT_END