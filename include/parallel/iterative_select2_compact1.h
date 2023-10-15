#pragma once

#include "../common.h"
#include "item_blocks.h"
#include "../misc/exclusive_sum.h"

NAMESPACE_PMT

template <typename Item>
class IterativeSelect2Compact1
{
public:
  static constexpr uint_fast8_t remove = 0;
  static constexpr uint_fast8_t copy_to_first_array = 1;
  static constexpr uint_fast8_t copy_to_second_array = 2;

  using item_t = Item;

  IterativeSelect2Compact1(size_t n, size_t max_block_length = 8192) :
    item_blocks_(n, max_block_length)
  {
    array2_lengths = new size_t[n + 1U];
    size_t n_threads = std::min(thread_pool.max_threads(), item_blocks_.n_blocks_);
    buffers = new item_t[max_block_length * n_threads];
  }

  ~IterativeSelect2Compact1()
  {
    delete[] buffers;
    delete[] array2_lengths;
  }

  size_t length() const
  {
    return item_blocks_.length();
  }

  ItemBlocks& item_blocks()
  {
    return item_blocks_;
  }

  template <typename functor_t>
  void iterate(
    item_t* items,
    item_t* compact,
    functor_t const& f,
    size_t* n_compacted)
  {
    size_t* partitions = item_blocks_.partitions_;
    ItemBlock* blocks = item_blocks_.blocks_;
    size_t max_block_length = item_blocks_.max_block_length_;
    
    thread_pool.for_all_blocks(item_blocks_.n_partitions_, [=](size_t p, thread_nr_t t) {
      size_t b_begin = partitions[p];
      size_t b_end = partitions[p + 1U];
      item_t* buffer = buffers + t * max_block_length;
      item_t* buffer_begin = buffer;
      item_t* buffer_end = buffer + max_block_length;

      size_t len1 = 0;
      size_t len2 = 0;

      for (size_t b = b_begin; b != b_end; ++b)
      {
        copy_to_buffer(items, b, buffer_begin, buffer_end, f, &len1, &len2);

        buffer_begin += len1;
        buffer_end -= len2;

        blocks[b].length_ = 0;
        array2_lengths[b] = 0;
      }

      // total lengths
      len1 = buffer_begin - buffer;
      len2 = buffer + max_block_length - buffer_end;

      blocks[b_begin].length_ = len1;
      array2_lengths[b_begin] = len2;

      copy_from_buffer(items + blocks[b_begin].offset_, buffer, len1, len2);
    });

    compact_array2(items, compact);

    *n_compacted = array2_lengths[item_blocks_.n_blocks_];

    //out(">>> " << *n_compacted);

    item_blocks_.concat_blocks();    
  }

  NO_INLINE void copy_from_buffer(item_t* items_begin, item_t* buffer, size_t len1, size_t len2)
  {
    item_t* items_end = items_begin + len1;
    item_t* buffer_begin = buffer;
    item_t* buffer_end = buffer + item_blocks_.max_block_length_;

    for (; items_begin != items_end; ++items_begin)
    {
      *items_begin = *buffer_begin++;
    }

    items_end = items_begin + len2;

    for (; items_begin != items_end; ++items_begin)
    {
      *items_begin = *(--buffer_end);
    }
  }

  void compact_array2(item_t* items, item_t* compact)
  {
    size_t* partitions = item_blocks_.partitions_;
    ItemBlock* blocks = item_blocks_.blocks_;
//    size_t max_block_length = item_blocks_.max_block_length_;

    exclusive_sum(array2_lengths, array2_lengths + item_blocks_.n_blocks_ + 1U);

    thread_pool.for_all_blocks(item_blocks_.n_partitions_, [=](size_t p, thread_nr_t t) {
      size_t b = partitions[p];
      size_t length = array2_lengths[b + 1U] - array2_lengths[b];
      item_t* items_begin = items + blocks[b].offset_ + blocks[b].length_;
      item_t* compact_begin = compact + array2_lengths[b];
      item_t* compact_end = compact_begin + length;

      for (; compact_begin != compact_end; ++compact_begin)
      {
        *compact_begin = *items_begin++;
      }
    });    
  }

  template <typename functor_t>
  NO_INLINE void copy_to_buffer(
    item_t* items,
    size_t b,
    item_t* buffer_begin,
    item_t* buffer_end,
    functor_t const& f,
    size_t* len1_ptr,
    size_t* len2_ptr)
  {
    ItemBlock* blocks = item_blocks_.blocks_;
    size_t len1 = 0;
    size_t len2 = 0;
    item_t out;
    item_t* items_begin = items + blocks[b].offset_;
    item_t* items_end = items_begin + blocks[b].length_;

    for (;items_begin != items_end; ++items_begin)
    {
      switch(f(*items_begin, &out))
      {
        case copy_to_first_array:
          *buffer_begin++ = out;
          ++len1;
          break;
        case copy_to_second_array:
          *(--buffer_end) = out;
          ++len2;
        default:;
      }
    }

    *len1_ptr = len1;
    *len2_ptr = len2;
  }

private:  
  ItemBlocks item_blocks_;
  size_t* array2_lengths = nullptr;
  item_t* buffers = nullptr;
};

NAMESPACE_PMT_END