#include "../common.h"
#include "../misc/edge.h"
#include "../sort/radix_sort_seq.h"
#include "../misc/integerhash.h"

NAMESPACE_PMT

template <
  typename index_t,
  typename value_t,
  typename attribute_t,
  typename functor1_t,
  typename functor2_t>
void tree_scan_seq(
  index_t* parents,
  size_t n,
  value_t* values,
  attribute_t* attributes,
  functor1_t const& w,
  functor2_t const& plus)
{
  using uvalue_t = decltype(unsigned_conversion(value_t(0)));
  using item_t = SortPair<uvalue_t, index_t>;  

  item_t* items = new item_t[n];
  item_t* aux = new item_t[n];

  for (size_t i = 0; i < n; ++i)
  {
    attributes[i] = w(i);
  }

  auto const &f_item =
    [=](index_t i) ALWAYS_INL_L(item_t)
    {
      return {unsigned_conversion(values[i]), i};
    };  

  item_t* sorted = pmt::radix_sort_seq(aux, items, n, f_item);

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