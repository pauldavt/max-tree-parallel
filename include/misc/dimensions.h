#pragma once
#include "../common.h"

NAMESPACE_PMT

/* Limited dimensions, for typically 1D to 4D integer vectors */

template <size_t NDimensions>
class Dimensions
{
public:
  static constexpr size_t n_dimensions = NDimensions;

  constexpr size_t length() const
  {
    size_t x = xs_[0];

    for (size_t d = 1; d < n_dimensions; ++d)
    {
      x *= xs_[d];
    }

    return x;
  }

  constexpr size_t& operator[](size_t d)
  {
    return xs_[d];
  }

  constexpr size_t const& operator[](size_t d) const
  {
    return xs_[d];
  }

  size_t xs_[n_dimensions];
};

NAMESPACE_PMT_END