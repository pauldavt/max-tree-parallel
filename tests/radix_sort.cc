#include "../include/sort/radix_sort_parallel.h"
#include "../include/sort/radix_sort_seq.h"
#include "../include/sort/sort_item.h"
#include "../include/sort/validate.h"
#include "../include/misc/timer.h"
#include "../include/misc/random.h"

using index_t = uint32_t;

constexpr index_t N = 33U * 1024U * 1024U;

template <typename value_t>
void check_value_t(char const* s)
{  
  using uvalue_t = decltype(pmt::unsigned_conversion(value_t(0)));

  using SortPair = pmt::SortPair<uvalue_t, index_t>;
  
  value_t* values = new value_t[N];
  index_t* sorted = new index_t[N];
  SortPair* aux_1 = new SortPair[N];
  SortPair* aux_2 = new SortPair[N];

  typename pmt::rng<uvalue_t>::type rnd;

  for (size_t i = 0; i < N; ++i)
  {
    if (std::is_floating_point<value_t>::value)
    {
      values[i] = pmt::random_fp(rnd);;
    }
    else
    {
      values[i] = rnd();
    }
  }

  auto const& f_item =
    [=](index_t i) ALWAYS_INL_L(SortPair) ALWAYS_INLINE
    {
      return {pmt::unsigned_conversion(values[i]), i};
    };  

  auto const& f_out =
    [=](index_t& out, SortPair const& item) ALWAYS_INLINE
    {
      out = item.data();
    };

  for (size_t i = 0; i < N; ++i)
  {
    sorted[i] = 0;
  }

  {
    pmt::Timer t;
    pmt::radix_sort_parallel(sorted, N, aux_1, aux_2, f_item, f_out);

    printf("Sorted %.2e random %s values in %f seconds\n", (double)N, s, t.stop());
  }

  pmt::validate_sort("parallel radix sort", values, sorted, N);

  for (size_t i = 0; i < N; ++i)
  {
    sorted[i] = 0;
  }

  {
    pmt::Timer t;
    pmt::radix_sort_seq(sorted, N, aux_1, aux_2, f_item, f_out);

    printf("Sorted %.2e random %s values in %f seconds (seq)\n", (double)N, s, t.stop());
  }  

  pmt::validate_sort("sequential radix sort", values, sorted, N);

  delete[] sorted;
  delete[] aux_2;
  delete[] aux_1;
  delete[] values;
}

int main()
{
  check_value_t<int8_t>("int8_t");
  check_value_t<uint8_t>("uint8_t");
  check_value_t<int16_t>("int16_t");
  check_value_t<uint16_t>("uint16_t");
  check_value_t<int32_t>("int32_t");
  check_value_t<uint32_t>("uint32_t");
  check_value_t<int64_t>("int64_t");
  check_value_t<uint64_t>("uint64_t");
  check_value_t<float>("float"); 
  check_value_t<double>("double");
}
