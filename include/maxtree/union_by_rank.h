#pragma once

#include "../common.h"
#include "../misc/edge.h"
#include "rank_set.h"

NAMESPACE_PMT
 
template <typename index_t>
void maxtree_union_by_rank(
  Edge<index_t>* sorted_edges,
  size_t idx_begin,
  size_t idx_end,
  RankSet<index_t>* sets,
  index_t* parents)
{
  Edge<index_t>* edges_end = sorted_edges + idx_end;
  Edge<index_t>* edges_begin = sorted_edges + idx_begin;

  while (edges_end-- > edges_begin)
  {
    Edge<index_t>& e = *edges_end;

    if (e.a_ == e.b_) continue;

    index_t set_a = compress_path(sets, e.a_);
    index_t set_b = compress_path(sets, e.b_);

    if (set_a == set_b) continue;

    index_t cc_root_b = sets[set_b].load();
    index_t cc_root_a = sets[set_a].load();
    index_t cc_root;

    if (parents[cc_root_b] != cc_root_b)
    { 
      // Parent of cc_root_b was set earlier. Make sure to set the parent
      // of cc_root_a to cc_root_b, to have a path to the lower partition.
      // This can happen if there are nodes in the connected component 
      // with the same image value (and higher index) as the connected
      // component root, but are processed earlier.
      cc_root = cc_root_b;
      
      parents[cc_root_a] = cc_root_b;
    }
    else
    {
      cc_root = cc_root_a;
      parents[cc_root_b] = cc_root_a;
    }

    merge_sets(sets, set_a, set_b, cc_root);
  }
} 

NAMESPACE_PMT_END