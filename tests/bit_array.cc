
#include "../include/common.h"
#include "../include/misc/bit_array.h"

int main()
{
  pmt::BitArray bit_array1(0);
  pmt::BitArray bit_array2(1);
  pmt::BitArray bit_array3(70);

  check(bit_array1.length() == 0);
  bit_array1.clear();

  check(bit_array2.length() == 1);
  bit_array2.clear();
  check(bit_array2.is_set(0) == false);
  bit_array2.set(0);
  check(bit_array2.is_set(0) == true);

  check(bit_array3.length() == 70);
  bit_array3.clear();

  for (size_t i = 0; i < 70; ++i)
  {
    check(!bit_array3.is_set(i));
  }

  bit_array3.set_range(0, 65);

  for (size_t i = 0; i <= 65; ++i)
  {
    check(bit_array3.is_set(i));
  }

  for (size_t i = 66; i < 70; ++i)
  {
    check(!bit_array3.is_set(i));
  }

  info("Success.");
  return 0;
}