#include "../common.h"
#include "../misc/edge.h"
#include "../sort/radix_sort_parallel.h"
#include "../misc/integerhash.h"
#include "../misc/random.h"
#include <vector>
#include "../parallel/iterative_select2_compact1.h"
#include "../misc/dynamic_stack.h"
#include "../parallel/thread_pool.h"

NAMESPACE_PMT

template <
  typename Index,
  typename Attribute,
  typename Functor1,
  typename Functor2>
class TreeContract;

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
void tree_scan(
  index_t* parents,
  size_t n,
  attribute_t* attributes,
  functor1_t const& w,
  functor2_t const& plus)
{
  TreeContract<
    index_t,
    attribute_t,
    functor1_t,
    functor2_t> tc(parents, n, attributes, w, plus);
}

template <
  typename Index,
  typename Attribute,
  typename Functor1,
  typename Functor2>
class TreeContract
{
public:  
  static constexpr uint_fast8_t n_hash_bits = 8U;

//  static constexpr size_t items_per_block = 8192U;

private:
  using index_t = Index;
  using edge_t = SortableEdgeByStart<index_t>;
  using functor1_t = Functor1;
  using functor2_t = Functor2;
  using attribute_t = Attribute;

  friend void tree_scan<index_t, attribute_t, functor1_t, functor2_t>(
    index_t *parents,
    size_t n,
    attribute_t *attributes,
    functor1_t const &w,
    functor2_t const &plus);

  struct EdgeArray
  {
    edge_t* edges_;
    size_t len_;
  };  

  TreeContract(
    index_t* parents,
    size_t n,
    attribute_t* attributes,
    functor1_t const& w,
    functor2_t const& plus);
  ~TreeContract();
  
  void merge_first_excluded_descendant();
  void select_first(IterativeSelect2Compact1<index_t>* select);
  void init();
  edge_t* sort_edges();

  ALWAYS_INLINE_F bool is_balanced(index_t node) const
  {
    return childs_[node] == n_;
  }

  ALWAYS_INLINE_F bool is_leaf(index_t node) const
  {
    return childs_[node] == node;
  }

  ALWAYS_INLINE_F bool is_linked_list_node(index_t node) const
  {
    return !is_balanced(node) && !is_leaf(node);
  }

  bool try_merge_and_check_if_leaf(index_t i);
  void contract(IterativeSelect2Compact1<index_t>* select, index_t root);

  index_t const* RESTRICT parents_;
  size_t n_;
  attribute_t* RESTRICT attributes_;
  edge_t* RESTRICT forward_ = nullptr;
  edge_t* RESTRICT forward_aux_ = nullptr;
  index_t* RESTRICT edge_indices_;
  index_t* RESTRICT edge_indices_aux_;
  index_t* RESTRICT childs_;
  IntegerHash<index_t> hash_;
  typename pmt::rng<index_t>::type rand_;
  size_t n_forward_;
  size_t total_compacted_;
  functor1_t const& w_;
  functor2_t const& plus_;
  std::vector<size_t> to_update_later_;
};

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
TreeContract<index_t, attribute_t, functor1_t, functor2_t>::TreeContract(
  index_t *parents,
  size_t n,
  attribute_t *attributes,
  functor1_t const& w,
  functor2_t const& plus) :
  parents_(parents), n_(n), attributes_(attributes), w_(w), plus_(plus)
{
  if (n <= 1) return;

  check(n <= ~index_t(0));

  childs_ = new index_t[n];
  forward_ = new edge_t[n];
  forward_aux_ = new edge_t[n];

  edge_t* sorted = sort_edges();

  if (sorted == forward_aux_)
  {
    std::swap(forward_, forward_aux_);
  }

  check(sorted == forward_);


  edge_indices_ = reinterpret_cast<index_t*>(forward_aux_);
  edge_indices_aux_ = edge_indices_ + n;

  n_forward_ = n - 1U; // n - 1 edges, the self loop for the root is moved to the end

  init();

  IterativeSelect2Compact1<index_t> select(n_forward_);


  select_first(&select);

  //check(forward_[n_ - 1U].a_ = n);

  {
    total_compacted_ = 0;
    while (select.length() > 0)
    {
      hash_.generate_vars(rand_);
      
      contract(&select, forward_[n_ - 1U].b_);
    }      
  }

  merge_first_excluded_descendant();
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
TreeContract<index_t, attribute_t, functor1_t, functor2_t>::~TreeContract()
{
  delete[] forward_aux_;
  delete[] forward_;
  delete[] childs_;
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
void TreeContract<index_t, attribute_t, functor1_t, functor2_t>::
merge_first_excluded_descendant()
{
  index_t* indices = edge_indices_aux_ + total_compacted_;
  for (size_t i = to_update_later_.size(); i--;)
  {            
    size_t n_to_update = to_update_later_[i];
    indices -= n_to_update;

    thread_pool.for_all(n_to_update, [=](index_t k, thread_nr_t thread_nr) ALWAYS_INLINE {
      edge_t const& edge = forward_[indices[k]];
      attributes_[edge.a_] = plus_(attributes_[edge.a_], attributes_[edge.b_]);
    });
  }
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
void TreeContract<index_t, attribute_t, functor1_t, functor2_t>::
select_first(IterativeSelect2Compact1<index_t>* select)
{
  auto const& flag = [=](index_t i, index_t o) ALWAYS_INL_L(bool)
  {
    edge_t const& current = forward_[i];
    edge_t const& right = forward_[i + index_t(1)];

    bool right_start_different =
      i + index_t(1) == n_forward_ || current.a_ != right.a_;

    if (right_start_different)
    {
      edge_indices_[o] = i;
      return true;
    }

    return false;
  };

  select->item_blocks().select(flag);
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
void TreeContract<index_t, attribute_t, functor1_t, functor2_t>::init()
{
  thread_pool.for_all(n_, [=](index_t i, thread_nr_t thread_nr) ALWAYS_INLINE {
    attributes_[i] = w_(i);
    childs_[i] = i;
  });

  // determine in-degrees and last child of every node
  thread_pool.for_all(n_forward_, [=](index_t i, pmt::thread_nr_t thread_nr) {
    edge_t const& left = forward_[i - index_t(1)];
    edge_t const& current = forward_[i];
    edge_t const& right = forward_[i + index_t(1)];

    bool left_end_different = i == 0 || current.a_ != left.a_;
    bool right_end_different = i + index_t(1) == n_forward_ || current.a_ != right.a_;

    if (right_end_different)
    {
      if (left_end_different)
      {
        childs_[current.a_] = current.b_;
        return;
      }

      childs_[current.a_] = n_;
    }    
  });
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
typename TreeContract<index_t, attribute_t, functor1_t, functor2_t>::edge_t
*TreeContract<index_t, attribute_t, functor1_t, functor2_t>::sort_edges()
{
  auto const& f_initial_item = [=](index_t i) ALWAYS_INL_L(edge_t)
  {
    if (parents_[i] == i)
    {
      return {index_t(n_), i};
    }

    return {parents_[i], i};
  };

  edge_t* sorted = pmt::radix_sort_parallel(forward_, forward_aux_, n_, f_initial_item);

  return sorted;
}

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
bool TreeContract<index_t, attribute_t, functor1_t, functor2_t>::
try_merge_and_check_if_leaf(index_t i)
{
  edge_t edge = forward_[i];
  //printf("try merge from %d -> %d\n", edge.a_, edge.b_);

  constexpr index_t init_stack_sz = 4096;
  index_t on_stack[init_stack_sz];

  DynamicStack<index_t> stack(&on_stack[0], init_stack_sz);

  //printf("edge.b_ = %d is linked_list_node: %d\n", edge.b_, childs_[edge.b_]);

  index_t current = edge.b_;
  while (is_linked_list_node(current) && hash_(current, n_hash_bits) != 0)
  {
    stack.insert(current);
    current = childs_[current];
  }

  bool future_leaf = is_leaf(current);

  if (!future_leaf && stack.length() == 0)
  {
    return false;
  }

  if (future_leaf)
  {
    index_t tail = current;      

    while (stack.length() > 0)
    {
      index_t next = stack.remove();
      childs_[next] = next;
      attributes_[next] = plus_(attributes_[next], attributes_[tail]);
      //merge_(next, tail);
      tail = next;        
    }

    attributes_[edge.a_] = plus_(attributes_[edge.a_], attributes_[tail]);
    //merge_(edge.a_, tail);            
    return true;
  }

  index_t relink_node = current;
  index_t tail = stack.remove();

  while (stack.length() > 0)
  {
    index_t next = stack.remove();
    childs_[next] = relink_node;
    attributes_[next] = plus_(attributes_[next], attributes_[tail]);
    //merge_(next, tail);
    tail = next;
  }
  
  attributes_[edge.a_] = plus_(attributes_[edge.a_], attributes_[tail]);
  //merge_(edge.a_, tail);
  forward_[i].b_ = relink_node;
  return false;
}


template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
void TreeContract<index_t, attribute_t, functor1_t, functor2_t>::
contract(IterativeSelect2Compact1<index_t>* select, index_t root)
{
  auto const& f = [=](index_t i)
  {                
    index_t edge_idx = edge_indices_[i];
    index_t start_point = forward_[edge_idx].a_;

    bool is_ll_node = edge_idx == 0 || forward_[edge_idx - 1U].a_ != start_point;

    if (is_ll_node && start_point != root && hash_(start_point, n_hash_bits) != 0)
    {
      // start_point is merged by another node
      return;
    }

    if (is_ll_node)
    {
      bool is_leaf = try_merge_and_check_if_leaf(edge_idx);

      if (is_leaf)
      {
        forward_[edge_idx].b_ = start_point; // mark as next iteration leaf
      }

      return;
    }

    index_t compacted_index = edge_idx;
    edge_t* forward_compacting = forward_ + edge_idx;
    size_t n_childs = 0;
    do
    {            
      if (!try_merge_and_check_if_leaf(edge_idx))
      {
        ++n_childs;
        *forward_compacting-- = forward_[edge_idx];
      }
    } while (edge_idx > 0 && forward_[--edge_idx].a_ == start_point);

    if (n_childs == 0)
    {
      forward_[compacted_index].b_ = start_point; // mark as next iteration leaf
      return;
    }

    // make sure removed edges to the left are not used
    if (forward_compacting >= forward_ && forward_compacting->a_ == start_point)
    {
      forward_compacting->a_ = n_;
    }
  };

  auto const& flag = [=](index_t const& in, index_t* out) ALWAYS_INL_L(uint_fast8_t)
  {
    index_t i = in;
    edge_t& edge = forward_[i];
    index_t start_point = forward_[i].a_;

    if (edge.b_ == edge.a_)
    {
      // merged leaf segment
      childs_[start_point] = start_point;
      return 0;
    }

    if (childs_[start_point] == start_point)
    {
      // merged leaf segment
      return 0;
    }

    *out = in;

    if (i > 0 && forward_[i - 1].a_ == start_point)
    {
      // balanced node
      //check(childs_[start_point] == ~index_t(0));
      return 1;
    }

    if (childs_[start_point] == n_ || start_point == root || hash_(start_point, n_hash_bits) == 0)
    {
      // linked list node
      childs_[start_point] = edge.b_;
      return 1;
    }

    // linked list node was contracted
    forward_[i] = {edge.a_, childs_[edge.a_]};
    return 2;
  };

  size_t n_compacted = 0;
  select->item_blocks().apply(f);
  select->iterate(edge_indices_, edge_indices_aux_ + total_compacted_, flag, &n_compacted);      
  to_update_later_.push_back(n_compacted);
  total_compacted_ += n_compacted;    
}



NAMESPACE_PMT_END