#pragma once

#include "../common.h"
#include <climits>

NAMESPACE_PMT

constexpr uint_fast8_t histo_sz_log2 = 8U;
constexpr uint_fast32_t block_size_sorting = 777U * 1024U;

template <typename uvalue_t>
constexpr size_t radix_sort_n_digits()
{
  return div_roundup(sizeof(uvalue_t) * CHAR_BIT, histo_sz_log2);
}

NAMESPACE_PMT_END
