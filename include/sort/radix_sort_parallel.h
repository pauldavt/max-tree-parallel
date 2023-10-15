#pragma once

#include "../common.h"
#include "sort.h"
#include "../misc/exclusive_sum.h"
#include "sort_item.h"
#include "../parallel/thread_pool.h"
#include "../misc/range.h"
#include "../misc/timer.h"

NAMESPACE_PMT

template<
  typename Item,
  typename Sorted,
  typename InitialItemF,
  typename LastItemF>
class RadixSortParallel
{
public:
  using index_t = size_t;
  using uvalue_t = typename Item::uvalue_t;
  using item_t = Item;
  using sorted_t = Sorted;
  using range_t = Range<size_t>;
  using initial_item_f = InitialItemF;
  using last_item_f = LastItemF;

  static constexpr size_t histo_sz = 1U << histo_sz_log2;
  static constexpr size_t histo_mask = histo_sz - 1U;

  // ensure that the number of bytes are divisible by 16, to avoid misaligned primitives,
  // which could cause race conditions near the start/end of item blocks
  static constexpr size_t items_per_block = 16U * ((block_size_sorting / sizeof(item_t)) / 16U);

  using histogram_t = index_t[histo_sz];

  RadixSortParallel(
    sorted_t* sorted,
    item_t* aux1,
    item_t* aux2,
    size_t n,
    unsigned bit_start,
    unsigned bit_end,
    InitialItemF const& f_initial_item,
    LastItemF const& f_last_item) :
    sorted_(sorted),
    aux1_(aux1),
    aux2_(aux2),
    n_(n),
    f_initial_item_(f_initial_item),
    f_last_item_(f_last_item),
    bits_{bit_start, bit_end}    
  {    
    global_offsets_ = new histogram_t[n_blocks()];

    sort_digits();
  }

  ~RadixSortParallel()
  {
    delete[] global_offsets_;
  }

private:
  template <typename item_f>
  void create_histograms(unsigned shift, item_f const& f_item)
  {
    thread_pool.for_all_blocks(n_blocks(), [=](size_t b, thread_nr_t thread_nr)
    {
      range_t range = make_range(b);
      histogram_t& g = global_offsets_[b];

      std::fill(g, g + histo_sz, 0);

      for (size_t i = range.begin_; i < range.end_; ++i)
      {
        uvalue_t u = f_item(i).unsigned_value() >> shift;
        ++g[u & histo_mask];
      }
    });
  }

  void make_offsets()
  {
    index_t* sums = new index_t[histo_sz + 1U];

    std::fill(sums, sums + histo_sz + 1U, index_t(0));

    for (size_t b = 0; b < n_blocks(); ++b)
    {
      histogram_t& g = global_offsets_[b];

      for (size_t i = 0; i < histo_sz; ++i)
      {
        index_t tmp = g[i];
        g[i] = sums[i];
        sums[i] += tmp;
      }
    }

    exclusive_sum(sums, sums + histo_sz + 1U);
    debug(sums[histo_sz] == n_);

    for (size_t b = 0; b < n_blocks(); ++b)
    {
      histogram_t& g = global_offsets_[b];

      for (size_t i = 0; i < histo_sz; ++i)
      {
        g[i] += sums[i];
      }
    }

    delete[] sums;
  }

  template <typename out_t, typename item_f, typename out_f>
  void scatter_digit(
    out_t* out,
    uint8_t shift,
    item_f const& f_item,
    out_f const& f_out)
  {    
    thread_pool.for_all_blocks(n_blocks(), [=](size_t b, thread_nr_t thread_nr) NO_INLINE
    {
      range_t range = make_range(b);          
      histogram_t& g = global_offsets_[b];

      for (size_t i = range.begin_; i < range.end_; ++i)
      {
        item_t const& item = f_item(i); 

        uvalue_t u = item.unsigned_value() >> shift;
        f_out(out[g[u & histo_mask]++], item);
      }
    });
  }

  void sort_digits()
  {
    unsigned n_digits = div_roundup(bits_[1] - bits_[0], histo_sz_log2);

    if (n_digits == 1)
    {
      sort_digit(sorted_, 0, f_initial_item_, f_last_item_);
      return;
    } 

    auto const& f_out =
      [=](item_t& out, item_t const& item) ALWAYS_INLINE
      {
        out = item;
      };

    sort_digit(aux1_, bits_[0], f_initial_item_, f_out);

    auto const& f_item =
      [=](size_t i) ALWAYS_INL_L(item_t const&)
      {
        return aux1_[i];
      };

    for (unsigned d = 1; d < n_digits - 1U; ++d)
    {
      sort_digit(aux2_, bits_[0] + d * histo_sz_log2, f_item, f_out);
      std::swap(aux1_, aux2_);
    }

    sort_digit(sorted_, bits_[0] + (n_digits - 1U) * histo_sz_log2, f_item, f_last_item_);
  }
    
  template <typename out_t, typename item_f, typename out_f>
  void sort_digit(
    out_t* out,
    unsigned shift,
    item_f const& f_item,
    out_f const& f_out)
  {
    create_histograms(shift, f_item);
    make_offsets();
    scatter_digit(out, shift, f_item, f_out);
  }

  constexpr range_t make_range(size_t block_nr) const
  {
    size_t begin = block_nr * items_per_block;
    size_t end = begin + std::min(n_ - begin, items_per_block);

    return {begin, end};
  }

  constexpr size_t n_blocks() const
  {
    return div_roundup(n_, items_per_block);
  }

  sorted_t* RESTRICT sorted_;
  item_t* RESTRICT aux1_;
  item_t* RESTRICT aux2_;
  size_t n_;
  initial_item_f const& f_initial_item_;
  last_item_f const& f_last_item_;
  histogram_t* RESTRICT global_offsets_;
  unsigned bits_[2];
};

template<
  typename sorted_t,
  typename item_t,
  typename initial_item_f,
  typename last_item_f>
void radix_sort_parallel(
  sorted_t* sorted,
  item_t* aux1,
  item_t* aux2,
  size_t n,
  unsigned bit_start,
  unsigned bit_end,
  initial_item_f const& f_initial_item,
  last_item_f const& f_last_item)
{
  if (bit_end <= bit_start || n <= 1) return;

  RadixSortParallel<item_t, sorted_t, initial_item_f, last_item_f> sorter(
      sorted, aux1, aux2, n, bit_start, bit_end, f_initial_item, f_last_item);
}

template<
  typename sorted_t,
  typename item_t,
  typename initial_item_f,
  typename last_item_f>
void radix_sort_parallel(
  sorted_t* sorted,
  size_t n,
  item_t* aux1,
  item_t* aux2,
  initial_item_f const& f_initial_item,
  last_item_f const& f_last_item)
{
  using uvalue_t = typename item_t::uvalue_t;

  RadixSortParallel<item_t, sorted_t, initial_item_f, last_item_f> sorter(
      sorted, aux1, aux2, n, 0, sizeof(uvalue_t) * CHAR_BIT, f_initial_item, f_last_item);
}

template<typename item_t, typename initial_item_f>
item_t* radix_sort_parallel(
  item_t* aux1,
  item_t* aux2,
  size_t n,
  initial_item_f const& f_initial_item)
{
  using uvalue_t = decltype(f_initial_item(0).unsigned_value());

  unsigned n_digits = radix_sort_n_digits<uvalue_t>();

  item_t* sorted = aux2;

  if (n_digits & 1)
  {
    sorted = aux1;
  }

  auto const& f_out =
    [=](item_t& out, item_t const& item) ALWAYS_INLINE
    {
      out = item;
    };

  radix_sort_parallel(sorted, aux1, aux2, n, 0U, sizeof(uvalue_t) * CHAR_BIT, f_initial_item, f_out);

  return sorted;
}

template<typename item_t>
item_t* radix_sort_parallel(
  item_t* aux,
  item_t* to_sort,
  size_t n)
{
  using uvalue_t = typename item_t::uvalue_t;

  unsigned n_digits = radix_sort_n_digits<uvalue_t>();

  item_t* sorted = to_sort;

  if (n_digits & 1)
  {
    sorted = aux;
  }

  auto const& f_item =
    [=](size_t i) ALWAYS_INL_L(item_t const&)
    {
      return to_sort[i];
    };  

  auto const& f_out =
    [=](item_t& out, item_t const& item) ALWAYS_INLINE
    {
      out = item;
    };

  radix_sort_parallel(sorted, aux, to_sort, n, 0, sizeof(uvalue_t) * CHAR_BIT, f_item, f_out);

  return sorted;
}

NAMESPACE_PMT_END