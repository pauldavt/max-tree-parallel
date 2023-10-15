#include "../../include/misc/bit_array.h"

NAMESPACE_PMT

void BitArray::clear()
{
  std::fill_n(data_, n_words(), size_t(0));
}

void BitArray::set_range(size_t x_begin, size_t x_end_inclusive)
{
  size_t first_word = (x_begin >> bits_per_word_log2);
  size_t last_word = (x_end_inclusive >> bits_per_word_log2);

  // begin_mask: 
  //   set bits 0 to (x_begin & bits_word_mask) exclusive
  //   to 0, otherwise 1
  size_t begin_mask = ~size_t(0U) << (x_begin & bits_word_mask);

  // end_mask: 
  //   set bits 0 to ((x_end - 1U) & bits_word_mask) inclusive
  //   to 1, otherwise 0
  //
  // if (x_end - 1U) == 0, 
  //   end_mask should be 1.
  // if (x_end - 1U) == bits_word_mask
  //   end_mask should be all ones.
  size_t n_bits = bits_word_mask - (x_end_inclusive & bits_word_mask);
  size_t end_mask = ~size_t(0U) >> n_bits;

  if (first_word == last_word)
  {
    data_[first_word] |= begin_mask & end_mask;
    return;
  }

  data_[first_word] |= begin_mask;

  for (++first_word; first_word < last_word; ++first_word)
  {
    data_[first_word] = ~size_t(0);
  }

  debug(first_word < n_words());

  data_[first_word] |= end_mask;
}  

NAMESPACE_PMT_END