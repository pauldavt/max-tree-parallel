#include "../common.h"
#include "../sort/radix_sort_seq.h"
#include "rootfix_seq.h"

NAMESPACE_PMT

template <
  typename index_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
void tree_scan_seq(
  index_t* parents,
  size_t n,
  attribute_t* attributes,
  functor1_t const& w,
  functor2_t const& plus)
{
  index_t* root_distance = new index_t[n];

  auto const &rootfix_plus =
    [](attribute_t a, attribute_t b) ALWAYS_INL_L(attribute_t)
    {
      return a + 1;
    };

  auto const &rootfix_w = [=](index_t i) ALWAYS_INL_L(attribute_t)
    {
      return 0;
    };

  rootfix_seq(parents, n, root_distance, rootfix_w, rootfix_plus);

  using item_t = SortPair<index_t, index_t>;  

  item_t* items = new item_t[n];
  item_t* aux = new item_t[n];

  auto const &f_item =
    [=](index_t i) ALWAYS_INL_L(item_t)
    {
      return {root_distance[i], i};
    };  

  item_t* sorted = pmt::radix_sort_seq(aux, items, n, f_item);

  delete[] root_distance;

  for (size_t i = 0; i < n; ++i)
  {
    attributes[i] = w(i);    
  }

  for (size_t i = n; i-- > 1;)
  {

    index_t node = sorted[i].data();
    index_t parent = parents[node];

    attributes[parent] = plus(attributes[parent], attributes[node]);
  }

  delete[] aux;
  delete[] items;  
}





NAMESPACE_PMT_END