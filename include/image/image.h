#pragma once

#include "../common.h"
#include "../misc/dimensions.h"

NAMESPACE_PMT

/*
 * An image is a hypermatrix with connectivity relations between
 * neighboring elements, which can be represented in a graph.
 */
template <typename Primitives>
class Image
{
public:
  static constexpr size_t n_dimensions = Primitives::n_dimensions;

  using prim = Primitives;
  using index_t = typename prim::index_t;
  using value_t = typename prim::value_t;
  using dim_t = Dimensions<n_dimensions>;  
  
  Image(Image const &img) = delete;
  Image(Image &&img) = delete;

  Image(value_t const* values, dim_t const& dimensions);

  constexpr dim_t const &dimensions() const { return dim_; }
  constexpr value_t const *values() const { return values_; }

private:
  value_t const* values_;
  dim_t const dim_;
};

template <typename prim>
Image<prim>::Image(value_t const* values, dim_t const& dimensions) :
    values_(values),
    dim_(dimensions)
{
}

template <
  typename index_t,
  typename value_t,
  size_t n_dimensions,
  size_t n_neighbors = 2U * n_dimensions>
struct image
{
  using prim = primitives<index_t, value_t, n_dimensions, n_neighbors>;
  using type = Image<prim>;
};

NAMESPACE_PMT_END
