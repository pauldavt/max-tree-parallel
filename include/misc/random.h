#pragma once

#include "../common.h"
#include "../pcg/pcg_random.hpp"
#include "../misc/bits.h"
#include <cmath>

NAMESPACE_PMT

template <typename U> struct rng {};
template <> struct rng<uint8_t> { using type = pcg32_unique; };
template <> struct rng<uint16_t> { using type = pcg32_unique; };
template <> struct rng<uint32_t> { using type = pcg32_unique; };
template <> struct rng<uint64_t> { using type = pcg64_unique; };

template <typename random_t, typename u_t>
ALWAYS_INLINE_F u_t random_uint(random_t& rng, u_t max_val)
{
  u_t const mask = ~u_t(0U) >> (sizeof(u_t) * CHAR_BIT - 1U - pmt::log2(max_val));
  u_t next;

  do
  {
    next = rng() & mask;
  } while (next > max_val);

  return next;
}

ALWAYS_INLINE_F float random_fp(pcg32_unique& rng)
{
  uint32_t u;
  float* f_ = reinterpret_cast<float*>(&u);
  float& f = *f_;

  do
  {
    u = rng();
  } while (std::isnan(f) || std::isinf(f));

  return f;
}

ALWAYS_INLINE_F double random_fp(pcg64_unique& rng)
{
  uint64_t u;
  double* f_ = reinterpret_cast<double*>(&u);
  double& f = *f_;

  do
  {
    u = rng();
  } while (std::isnan(f) || std::isinf(f));

  return f;
}

template <typename item_rai_t, typename index_t, typename random_t>
void random_shuffle(item_rai_t items, index_t len, random_t& rng)
{
  using item_t = typename std::iterator_traits<item_rai_t>::value_type;

  for (index_t i = 0; i < len; ++i)
  {
    index_t k = random_uint(rng, len - 1);

    item_t tmp = items[k];
    items[k] = items[i];
    items[i] = tmp;
  }
}

NAMESPACE_PMT_END
