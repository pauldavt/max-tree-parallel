
#pragma once

#include "../common.h"
#include "bits.h"
#include "../misc/random.h"

NAMESPACE_PMT

template <typename U>
class IntegerHash
{
public:
  // i < 2^(sizeof(U) * CHAR_BIT) - 1
  INLINE U operator()(U i, uint_fast8_t n_bits) const
  {
    debug(i != ~U(0));
    return (a_ * (i + U(1)) + b_) >> uint_fast8_t(sizeof(U) * CHAR_BIT - n_bits);
  }

  template <typename R>
  void generate_vars(R& rng)
  {
    a_ = rng() | U(1);
    b_ = rng() & (~U(0) >> 1);
  }

private:
  U a_;
  U b_;
};

NAMESPACE_PMT_END