#pragma once

#include "../common.h"

NAMESPACE_PMT


template <typename Iter>
void exclusive_sum(Iter iter, Iter const iter_end, Iter out)
{
  using Value = typename std::iterator_traits<Iter>::value_type;

  Value sum = 0;

  while(iter != iter_end)
  {
    Value tmp = *iter;
    *out++ = sum;
    sum += tmp;
    ++iter;
  }
}

template <typename Iter>
void exclusive_sum(Iter iter, Iter const iter_end)
{
  exclusive_sum(iter, iter_end, iter);
}



NAMESPACE_PMT_END