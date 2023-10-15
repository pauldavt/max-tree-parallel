#pragma once

#include "../common.h"
#include "thread_pool.h"

NAMESPACE_PMT

struct ItemBlock
{
  size_t offset_;
  size_t length_;
};

class ItemBlocks
{
public:
  template <typename Item>
  friend class IterativeSelect2Compact1;
  using item_block_t = ItemBlock;

  ItemBlocks(size_t n, size_t max_block_length = default_n_items_per_block) :
    max_n_(n), n_(n), max_block_length_(max_block_length)
  {
    determine_blocks();
  }

  ~ItemBlocks()
  {
    delete[] blocks_;
    delete[] partitions_;
  }

  void reset(size_t n)
  {
    check(n <= max_n_);
    n_ = n;
    determine_blocks();
  }

  template <typename functor_t>
  void apply(functor_t const& f) const
  {
    thread_pool.for_all_blocks(n_partitions_, [=](size_t p, thread_nr_t t) {
      size_t b_begin = partitions_[p];
      size_t b_end = partitions_[p + 1U];

      for (size_t b = b_begin; b != b_end; ++b)
      {
        apply_range(b, f);
      }
    });
  }

  template <typename functor_t>
  void select(functor_t const& f)
  { 
    thread_pool.for_all_blocks(n_partitions_, [=](size_t p, thread_nr_t t) {
      size_t b_begin = partitions_[p];
      size_t b_end = partitions_[p + 1U];

      debug(b_begin < b_end);

      size_t offset = blocks_[b_begin].offset_;

      for (size_t b = b_begin; b != b_end; ++b)
      {
        offset = select_range(f, b, offset);
        blocks_[b].length_ = 0;
      }

      blocks_[b_begin].length_ = offset - blocks_[b_begin].offset_;
    });

    concat_blocks();
  }

  size_t length() const
  {
    return n_;
  }

  size_t n_partitions() const
  {
    return n_partitions_;
  }


private:
  template <typename functor_t>
  NO_INLINE size_t select_range(
    functor_t const& f,
    size_t b,
    size_t offset)
  {
    size_t begin = blocks_[b].offset_;
    size_t end = begin + blocks_[b].length_;

    for (; begin != end; ++begin)
    {
      debug(offset <= begin);
      if (f(begin, offset))
      {
        ++offset;
      }
    }

    return offset;
  }

  void concat_blocks()
  {
    if (n_blocks_ == 0)
    {
      return;
    }

    n_ = 0;
    partitions_[0] = 0;
    n_partitions_ = 0;
    //size_t b_begin = 0;
    size_t block_len = 0;
    size_t new_b = 0;
    for (size_t b = 0; b != n_blocks_; ++b)
    {
      if (blocks_[b].length_ == 0)
      {
        continue;
      }

      n_ += blocks_[b].length_;

      if (block_len + blocks_[b].length_ > max_block_length_)
      {
        partitions_[++n_partitions_] = new_b;
        block_len = 0;
      }

      block_len += blocks_[b].length_;
      blocks_[new_b++] = blocks_[b];
    }

#ifdef PMT_DEBUG
    check(block_len > 0 || n_ == 0);
#endif

    if (n_ == 0)
    {
      return;
    }
    
    n_blocks_ = new_b;
    partitions_[++n_partitions_] = new_b;
  }

  template <typename functor_t>
  NO_INLINE void apply_range(size_t b, functor_t const& f) const
  {
    size_t begin = blocks_[b].offset_;
    size_t end = begin + blocks_[b].length_;

    for (; begin != end; ++begin)
    {
      f(begin);
    }
  }

  void determine_blocks()
  {
    n_blocks_ = div_roundup(n_, max_block_length_);

    if (n_blocks_ == 0)
    {
      n_partitions_ = 0;
      return;
    }

    if (blocks_ == nullptr)
    {
      blocks_ = new item_block_t[n_blocks_];
    }

    size_t offset = 0;
    for (size_t i = 0; i < n_blocks_ - 1U; ++i)
    {
      blocks_[i] = {offset, max_block_length_};
      offset += max_block_length_;
    }

    blocks_[n_blocks_ - 1U] = {offset, n_ - offset};

    if (partitions_ == nullptr)
    {
      partitions_ = new size_t[n_blocks_ + 1U];
    }

    for (size_t i = 0; i < n_blocks_; ++i)
    {
      partitions_[i] = i;
    }

    partitions_[n_blocks_] = n_blocks_;
    n_partitions_ = n_blocks_;
  }

  size_t max_n_ = 0;
  size_t n_ = 0;
  size_t max_block_length_ = 0;
  size_t n_blocks_ = 0;
  item_block_t* blocks_ = nullptr;
  size_t* partitions_ = nullptr;
  size_t n_partitions_ = 0;
};

NAMESPACE_PMT_END
