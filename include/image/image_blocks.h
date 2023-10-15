#pragma once

#include "../common.h"
#include "../misc/dimensions.h"
#include "image.h"
#include "image_block.h"

NAMESPACE_PMT

/*
 * An image partitioned into blocks
 */
template <typename Primitives>
class ImageBlocks
{
public:  
  using prim = Primitives;
  using index_t = typename prim::index_t;
  using value_t = typename prim::value_t;
  using image_t = Image<prim>;
  using dim_t = typename image_t::dim_t;
  using block_t = Block<prim::n_dimensions>;

  static constexpr unsigned n_dimensions = prim::n_dimensions;

  ImageBlocks(image_t const& image);

  constexpr image_t const& image() const { return image_; }
  constexpr dim_t const& dimensions() const { return dimensions_; }

private:
  image_t const& image_;
  dim_t dimensions_;  
};

template <typename prim>
ImageBlocks<prim>::ImageBlocks(image_t const& image) : image_(image)
{
  for (dim_idx_t d = 0; d != n_dimensions; ++d)
  {
    dimensions_[d] =
      div_roundup(image.dimensions()[d], block_t::max_dimensions[d]);
  }
} 


NAMESPACE_PMT_END