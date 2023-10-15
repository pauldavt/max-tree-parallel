#pragma once

#include "../common.h"
#include "bits.h"
#include "logger.h"

NAMESPACE_PMT

class BitArray
{
public:
  static constexpr size_t bits_per_word = sizeof(size_t) * CHAR_BIT;

  BitArray(size_t length);
  ~BitArray();
  BitArray(BitArray const&) = delete;
  BitArray& operator=(BitArray const&) = delete;

  size_t length() const;
  void set(size_t i);
  void clear(size_t i);
  bool is_set(size_t i) const;
  void clear();
  void set_range(size_t x_begin, size_t x_end_inclusive);

private:
  constexpr size_t n_words() const
  {
    return div_roundup(length_, bits_per_word);
  }

  size_t* data_;
  size_t length_;
};

inline BitArray::BitArray(size_t length) : length_(length)
{
  data_ = new size_t[n_words()];
}

inline BitArray::~BitArray()
{
  delete[] data_;
}

ALWAYS_INLINE_F size_t BitArray::length() const
{
  return length_;
}

ALWAYS_INLINE_F void BitArray::set(size_t i)
{
  size_t word_nr = i / bits_per_word;
  size_t bit_nr = i % bits_per_word;

  debug(word_nr < n_words());

  pmt::set_bit(data_[word_nr], bit_nr);
}

ALWAYS_INLINE_F void BitArray::clear(size_t i)
{
  size_t word_nr = i / bits_per_word;
  size_t bit_nr = i % bits_per_word;

  debug(word_nr < n_words());

  pmt::clear_bit(data_[word_nr], bit_nr);
}

ALWAYS_INLINE_F bool BitArray::is_set(size_t i) const
{    
  size_t word_nr = i / bits_per_word;
  size_t bit_nr = i % bits_per_word;

  debug(word_nr < n_words());

  return pmt::is_bit_set(data_[word_nr], bit_nr);
}

NAMESPACE_PMT_END