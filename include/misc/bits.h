#pragma once

#include <climits>
#include "../common.h"

NAMESPACE_PMT

constexpr unsigned log2(unsigned int i)
{
  return CHAR_BIT * sizeof(i) - __builtin_clz(i) - 1U;
}

constexpr unsigned log2(unsigned long i)
{
  return CHAR_BIT * sizeof(i) - __builtin_clzl(i) - 1U;
}

constexpr unsigned log2(unsigned long long i)
{
  return CHAR_BIT * sizeof(i) - __builtin_clzll(i) - 1U;
}

constexpr size_t bits_per_word_log2 =
  pmt::log2(sizeof(size_t) * CHAR_BIT);
constexpr size_t bits_per_word = sizeof(size_t) * CHAR_BIT;
constexpr size_t bits_word_mask = sizeof(size_t) * CHAR_BIT - 1U;

template <typename U>
void set_bit(U& value, unsigned bit_nr)
{
  value |= U(1) << bit_nr;
}

template <typename U>
void clear_bit(U& value, unsigned bit_nr)
{
  value &= ~(U(1) << bit_nr);
}

template <typename U>
bool is_bit_set(U& value, unsigned bit_nr)
{
  return (value >> bit_nr) & 1U;
}

NAMESPACE_PMT_END