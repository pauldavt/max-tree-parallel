#include "../common.h"
#include <algorithm>

NAMESPACE_PMT

template <typename Value, typename Index>
void validate_sort(std::string const& info, Value* values, Index* rank_to_index, size_t n)
{
  Index* counts = new Index[n];
  std::fill(counts, counts + n, Index(0));

  // Check if rank_to_index is a permutation
  for (size_t i = 0; i < n; ++i)
  {
    ++counts[rank_to_index[i]];
  }

  for (size_t i = 0; i < n; ++i)
  {
    check(counts[i] == 1U);
  }

  delete[] counts;

  // Check order
  for (size_t i = 1; i < n; ++i)
  {
    check(values[rank_to_index[i]] > values[rank_to_index[i - 1]] || 
      (values[rank_to_index[i]] == values[rank_to_index[i - 1]] &&
      rank_to_index[i] > rank_to_index[i - 1]));
  }
}

NAMESPACE_PMT_END