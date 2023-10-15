#pragma once

#include "../common.h"
#include "../misc/dimensions.h"
#include "../misc/coordinate.h"
#include "../misc/bit_array.h"
#include "image.h"
#include "image_blocks.h"

NAMESPACE_PMT

template <typename Primitives>
class ImageBlocks;

template <size_t n_dims>
struct Block;

/* 
 * An image block represents a piece of an image with constant maximum dimension
 * lengths.
 */
template <typename Primitives>
class ImageBlock
{
public:  
  using prim = Primitives;
  using index_t = typename prim::index_t;
  using value_t = typename prim::value_t;
  using block_t = Block<prim::n_dimensions>;
  using block_index_t = typename block_t::block_index_t;
  using image_blocks_t = ImageBlocks<prim>;
  using prim_block =
    primitives<block_index_t, value_t, prim::n_dimensions, prim::n_neighbors>;
  using dim_t = Dimensions<prim::n_dimensions>;
  using block_vec_t = Coordinate<prim_block>;
  using block_loc_t = Coordinate<prim>;
  using image_t = Image<prim>;

  static constexpr unsigned n_dimensions = prim::n_dimensions;

  ImageBlock(image_blocks_t const& ib, block_loc_t const& block_loc);

  /*
   * Apply a function to every element in the block
   */
  template <typename functor_t>
  void apply(functor_t const& f) const;

  index_t block_nr() const;
  
  /*
   * Assign the true value to boundary elements
   */
  void set_boundaries(BitArray* bit_array) const;

  constexpr dim_t const& dimensions() const { return dim_; }
  /*
   * Location in the block grid
   */
  constexpr block_loc_t const& location() const { return block_loc_; }
  /*
   * Offset of the first block element in the image
   */
  constexpr index_t global_offset() const { return global_offset_; }

  static void determine_dimensions(
    image_t const& img,
    block_loc_t const& block_loc,
    dim_t* dims);

private:
    
  template <typename functor_t>
  void apply_line(
    functor_t const& f,
    index_t global_index,
    block_index_t index_in_block) const;

  template <typename functor_t>
  void apply_dim(
    functor_t const& f,
    index_t global_index,
    block_index_t index_in_block,
    dim_idx_t d,
    index_t* skip_img,
    block_index_t* skip_block) const;

  void set_boundaries(
    BitArray* bit_array,
    index_t offset,
    index_t* skip,
    dim_idx_t d,
    dim_idx_t d_exclude) const;


  void determine_global_offset();

  image_blocks_t ib_;
  index_t global_offset_;
  block_loc_t block_loc_;
  dim_t dim_;
};

template <typename Primitives>
ImageBlock<Primitives>::ImageBlock(
  image_blocks_t const& ib,
  block_loc_t const& block_loc) :
  ib_(ib),
  block_loc_(block_loc)
{    
  determine_dimensions(ib_.image(), block_loc, &dim_);
  determine_global_offset();
}


template <typename prim>
template <typename functor_t>
void ImageBlock<prim>::apply_line(
  functor_t const& f,
  index_t global_index,
  block_index_t index_in_block) const
{
  size_t i = global_index;
  size_t global_line_end = i + dim_[0];

  while (i != global_line_end)
  {
    f(index_t(i), index_in_block);

    ++i;
    ++index_in_block;
  }
}

template <typename prim>
template <typename functor_t>
void ImageBlock<prim>::apply_dim(
  functor_t const& f,
  index_t global_index,
  block_index_t index_in_block,
  dim_idx_t d,
  index_t* skip_img,
  block_index_t* skip_block) const
{
  if (d == 0)
  {
    apply_line(f, global_index, index_in_block);
    return;
  }

  for (size_t i = 0; i != dim_[d]; ++i)
  {
    apply_dim(f, global_index, index_in_block, d - 1U, skip_img, skip_block);
    global_index += skip_img[d - 1U];
    index_in_block += skip_block[d - 1U];
  }
}

/*
  * Apply a function to every element in the block
  */
template <typename prim>
template <typename functor_t>
void ImageBlock<prim>::apply(functor_t const& f) const
{
  index_t global_index = global_offset_;
  block_index_t index_in_block = 0;
  index_t skip_img[n_dimensions - 1U];
  block_index_t skip_block[n_dimensions - 1U];
  dim_t img_dims = ib_.image().dimensions();

  skip_img[0] = img_dims[0];
  skip_block[0] = dim_[0];

  for (size_t d = 1; d < n_dimensions - 1U; ++d)
  {
    skip_img[d] = skip_img[d - 1U] * img_dims[d];
    skip_block[d] = skip_block[d - 1U] * dim_[d];
  }

  apply_dim(
    f,
    global_index,
    index_in_block,
    n_dimensions - 1U,
    skip_img,
    skip_block);
}

template <typename prim>
typename prim::index_t ImageBlock<prim>::block_nr() const
{
  dim_t grid_dims = ib_.dimensions();
  index_t nr = block_loc_[0];
  index_t factor = 1;

  for (dim_idx_t d = 1; d < n_dimensions; ++d)
  {
    factor *= grid_dims[d - dim_idx_t(1)];
    nr += factor * block_loc_[d];
  }

  return nr;
}   

template <typename prim>
void ImageBlock<prim>::set_boundaries(
  BitArray* bit_array,
  index_t offset,
  index_t* skip,
  dim_idx_t d,
  dim_idx_t d_exclude) const
{
  if (d == 0 && d_exclude != 0)
  {      
    bit_array->set_range(offset, offset + dim_[0] - index_t(1));
    return;
  }

  if (d == 1 && d_exclude == 0U)
  {
    for (size_t i = 0; i != dim_[1]; ++i)
    {
      bit_array->set(offset);
      offset += dim_[0];
    }
    return;
  }

  if (d == d_exclude)
  {
    set_boundaries(bit_array, offset, skip, d - dim_idx_t(1), d_exclude);
    return;
  }

  for (index_t i = 0; i < dim_[d]; ++i)
  {
    set_boundaries(bit_array, offset, skip, d - dim_idx_t(1), d_exclude);
    offset += skip[d];
  }
}

/*
  * Given the dimensions of this block, set the boundaries of bit_array to true
  */
template <typename prim>
void ImageBlock<prim>::set_boundaries(BitArray* bit_array) const
{
  check(dim_.length() <= bit_array->length());

  dim_t const& grid_dims = ib_.dimensions();

  if (n_dimensions == 1U)
  {
    if (block_loc_[0] > 0)
    {
      bit_array->set(0);
    }

    if (block_loc_[0] < grid_dims[0] - size_t(1))
    {
      bit_array->set(dim_[0] - size_t(1));
    }
    return;      
  }

  // a node is on the boundary if it's on the boundary in a single dimension

  index_t skip[n_dimensions];
  skip[0] = 1;

  for (dim_idx_t d = 1; d < n_dimensions; ++d)
  {
    skip[d] = skip[d - dim_idx_t(1)] * dim_[d - dim_idx_t(1)];
  }

  for (dim_idx_t d_exclude = 0; d_exclude < n_dimensions; ++d_exclude)
  {
    index_t offset = 0;

    if (block_loc_[d_exclude] > 0)
    {
      set_boundaries(
        bit_array,
        offset,
        skip,
        n_dimensions - dim_idx_t(1),
        d_exclude);
    }

    offset += (dim_[d_exclude] - size_t(1)) * skip[d_exclude];

    if (block_loc_[d_exclude] < grid_dims[d_exclude] - size_t(1))
    {
      set_boundaries(
        bit_array,
        offset,
        skip,
        n_dimensions - dim_idx_t(1),
        d_exclude);
    }
  }
}

template <typename prim>
void ImageBlock<prim>::determine_dimensions(
  image_t const& img,
  block_loc_t const& block_loc,
  dim_t* dims)
{
  dim_t const& img_dims = img.dimensions();

  for (size_t d = 0; d != n_dimensions; ++d)
  {
    size_t len = block_t::max_dimensions[d];
    index_t x = block_loc[d] * block_t::max_dimensions[d];

    if (len > img_dims[d] - x)
    {
      len = img_dims[d] - x;
    }

    (*dims)[d] = len;
  }
}

template <typename prim>
void ImageBlock<prim>::determine_global_offset()
{
  global_offset_ = block_loc_[0] * block_t::max_dimensions[0];

  dim_t const& img_dims = ib_.image().dimensions();

  index_t skip = 1;
  for (size_t d = 1; d < n_dimensions; ++d)
  {
    skip *= img_dims[d - 1U];
    global_offset_ += skip * block_loc_[d] * block_t::max_dimensions[d];
  }
}

template <size_t n_dims>
constexpr size_t block_determine_max_length(size_t const *dims)
{
  size_t ret = dims[0];

  for (size_t d = 1; d != n_dims; ++d)
  {
    ret *= dims[d];
  }

  return ret;
}

template <size_t n_dims>
struct Block
{
  using block_index_t = uint16_t;
  static constexpr size_t max_dimensions[1] = {65536};
  static constexpr size_t max_length = block_determine_max_length<1>(max_dimensions);
};

template <>
struct Block<2U>
{
  using block_index_t = uint16_t;
  static constexpr size_t max_dimensions[2] = {256, 256};
  static constexpr size_t max_length = block_determine_max_length<2>(max_dimensions);
};

template <>
struct Block<3U>
{
  using block_index_t = uint16_t;
  static constexpr size_t max_dimensions[3] = {64, 32, 32};  
  static constexpr size_t max_length = block_determine_max_length<3>(max_dimensions);
};

NAMESPACE_PMT_END