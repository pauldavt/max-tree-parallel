#include "../common.h"
#include "unsigned_conversion.h"

NAMESPACE_PMT


template <typename value_t, typename index_t>
struct Quantile
{
  using uvalue_t = decltype(pmt::unsigned_conversion(value_t(0)));

  Quantile() {}
  Quantile(value_t v, index_t i) : val_(v), index_(i) {}
  
  INLINE uvalue_t unsigned_conversion()
  {
    return pmt::unsigned_conversion(val_);
  }

  INLINE bool less_than_or_equal(value_t v, index_t i)
  {
    return val_ < v || (val_ == v && index_ <= i);
  }

  INLINE static size_t determine_partition(
    value_t v,
    index_t i,
    Quantile* quantiles,
    size_t n_partitions)
  {
    size_t min = 0;
    size_t max = n_partitions;

    while (min + 1U < max)
    {
      size_t mid = min + (max - min) / 2U;

      if (quantiles[mid].less_than_or_equal(v, i))
      {
        min = mid;
      }
      else
      {
        max = mid;
      }
    }

    return min;
  }


  INLINE Quantile const& data()
  {
    return *this;
  }

private:
  value_t val_;
  index_t index_;
};

NAMESPACE_PMT_END