#pragma once
#include "../common.h"
#include "dimensions.h"

NAMESPACE_PMT

/* Typically 1D to 4D integer vectors */

template <typename Primitives>
class Coordinate
{
public:
  using prim = Primitives;
  using index_t = typename prim::index_t;
  using dim_t = Dimensions<prim::n_dimensions>;

  static constexpr size_t n_dimensions = prim::n_dimensions;

  static Coordinate from_index(index_t index, dim_t const& dims);

  constexpr index_t index(dim_t const& dims) const
  {
    index_t result = 0;
    index_t factor = 1;

    for (dim_idx_t d = 0; d != n_dimensions; ++d)
    {
      result += factor * xs_[d];      
      factor *= dims[d];
    }

    return result;    
  }
  
  constexpr index_t& operator[](size_t d)
  {
    return xs_[d];
  }

  constexpr index_t const& operator[](size_t d) const
  {
    return xs_[d];
  }  

  void init_zeros();
  void inc_index(dim_t const& dims);
   
  index_t xs_[n_dimensions];
};

template <typename prim>
Coordinate<prim> Coordinate<prim>::from_index(index_t index, dim_t const& dims)
{
  Coordinate result;

  for (dim_idx_t d = 0; d != n_dimensions; ++d)
  {
    result.xs_[d] = index % dims[d];
    
    index /= dims[d];
  }

  return result;
}

template <typename prim>
void Coordinate<prim>::init_zeros()
{
  for (dim_idx_t d = 0; d < n_dimensions; ++d)
  {
    xs_[d] = 0;
  }
}

template <typename prim>
void Coordinate<prim>::inc_index(dim_t const& dims)
{
  ++xs_[0];
  for (dim_idx_t d = 0; d < n_dimensions - 1U; ++d)
  {
    if (xs_[d] == dims[d])
    {
      xs_[d] = 0;
    }
    else if (xs_[d] != 0)
    {
      break;
    }

    ++xs_[d + dim_idx_t(1)];
  }
}


NAMESPACE_PMT_END