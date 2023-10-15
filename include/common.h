#pragma once

#include <cstdlib>
#include <cstdint>
#include <type_traits>
#include <algorithm>

#define ALWAYS_INLINE __attribute__((always_inline))
#define ALWAYS_INLINE_F ALWAYS_INLINE inline
#define INLINE ALWAYS_INLINE_F
#define NO_INLINE __attribute__ ((noinline))
#define RESTRICT __restrict__

#define L_RET(...) -> __VA_ARGS__
#if defined(__clang__) || __GNUC__ < 9 || (__GNUC__ == 9 && __GNUC_MINOR__ <= 1)
#  define ALWAYS_INL_L(...) ALWAYS_INLINE -> __VA_ARGS__
#else
#  define ALWAYS_INL_L(...) -> __VA_ARGS__ ALWAYS_INLINE
#endif

#define PACKED __attribute__((packed))

#define NAMESPACE_PMT namespace pmt {
#define NAMESPACE_PMT_END }

NAMESPACE_PMT

using dim_idx_t = uint_fast8_t;
using thread_nr_t = uint_fast16_t;
using partition_t = uint8_t; // change to uint16_t if more than 256 threads

constexpr size_t cacheline_len = 64U;

constexpr size_t default_max_threads_limit = 256;
extern size_t hardware_concurrency;
class ThreadPool;
extern ThreadPool thread_pool;

constexpr size_t default_n_items_per_block = 64U * 1024U;

template <
  typename Index,
  typename Value = Index,
  size_t NDimensions = 1,
  size_t NNeighbors = 2 * NDimensions,
  std::enable_if_t<
    std::is_integral<Index>::value &&
    std::is_unsigned<Index>::value &&
    (std::is_integral<Value>::value ||
    std::is_floating_point<Value>::value), bool> = true>
struct primitives
{
  using index_t = Index;
  using value_t = Value;

  static constexpr size_t n_dimensions = NDimensions;
  static constexpr size_t n_neighbors = NNeighbors;
};

template <typename A, typename B>
constexpr A div_roundup(A const& a, B const& b)
{
  A result = a / b;

  if (a % b)
  {
    ++result;
  }

  return result;
}

NAMESPACE_PMT_END
