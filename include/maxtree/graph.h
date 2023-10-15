#pragma once

#include "../common.h"
#include "../misc/edge.h"
#include <limits>
#include "../misc/integerhash.h"
#include "../misc/random.h"
#include "../misc/quantile.h"
#include "../sort/radix_sort_parallel.h"
#include "../misc/unsigned_conversion.h"

NAMESPACE_PMT

/*
 * A graph is partitioned into subgraphs.
 * Local edges are edges between nodes in the same subgraph.
 * Global edges are edges between any two nodes in the graph.
 */
template <typename index_t>
class Graph
{
public:  
  using edge_t = Edge<index_t>;

  Graph(size_t n_subgraphs, size_t max_nodes, size_t max_edges);
  ~Graph();

  inline size_t max_edges() const { return max_edges_; }
  inline size_t max_nodes() const { return max_nodes_; }
  inline edge_t* edges() const { return edges_; }

  inline edge_t* subgraph(index_t subgraph_nr) const
  {
    return edges_ + subgraph_offsets_[subgraph_nr];
  }  

  inline void set_subgraph_offset(index_t subgraph_nr, size_t offset)
  {
    subgraph_offsets_[subgraph_nr] = offset;    
  }

  inline void set_local_edge_count(index_t subgraph_nr, size_t count)
  {
    local_edge_counts_[subgraph_nr] = count;
  }

  inline void set_global_edge_count(index_t subgraph_nr, size_t count)
  {
    global_edge_counts_[subgraph_nr] = count;
  }

  inline size_t subgraph_offset(index_t subgraph_nr) const
  {
    return subgraph_offsets_[subgraph_nr];
  }

  inline size_t local_edge_count(index_t subgraph_nr) const
  {
    return local_edge_counts_[subgraph_nr];
  }

  inline size_t global_edge_count(index_t subgraph_nr) const
  {
    return global_edge_counts_[subgraph_nr];
  }

  void determine_n_edges();

  inline size_t n_edges() const
  {
    return n_edges_;
  }

  inline size_t n_subgraphs() const
  {
    return n_subgraphs_;
  }

  inline size_t edge_count(index_t subgraph_nr) const
  {
    return global_edge_counts_[subgraph_nr] + local_edge_counts_[subgraph_nr];
  }


private:
  edge_t* edges_ = nullptr;
  size_t* subgraph_offsets_ = nullptr;
  size_t* local_edge_counts_;
  size_t* global_edge_counts_;
  size_t max_nodes_ = 0;
  size_t max_edges_ = 0;
  size_t n_edges_ = 0;
  size_t n_subgraphs_ = 0;
};

template <typename index_t>
Graph<index_t>::Graph(size_t n_subgraphs, size_t max_nodes, size_t max_edges) :
  max_nodes_(max_nodes), max_edges_(max_edges)
{
  if (n_subgraphs == 0) return;

  n_subgraphs_ = n_subgraphs;
  subgraph_offsets_ = new size_t[3U * (n_subgraphs + 1U)];
  local_edge_counts_ = subgraph_offsets_ + n_subgraphs + 1U;
  global_edge_counts_ = local_edge_counts_ + n_subgraphs + 1U;

  edges_ = new edge_t[max_edges];
}

template <typename index_t>
Graph<index_t>::~Graph()
{
  delete[] subgraph_offsets_;
  delete[] edges_;
}

template <typename index_t>
void Graph<index_t>::determine_n_edges()
{
  size_t total = 0;

  for (size_t i = 0; i < n_subgraphs_; ++i)
  {
    total += local_edge_counts_[i];
    total += global_edge_counts_[i];
  }

  n_edges_ = total;
}


NAMESPACE_PMT_END