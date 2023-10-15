#pragma once

#include "../common.h"

NAMESPACE_PMT

template <typename UValue, typename Data>
struct SortPair
{
  using uvalue_t = UValue;

  SortPair() {}  
  SortPair(UValue val, Data data) : val_(val), data_(data) {}

  INLINE UValue unsigned_value() const
  {
    return val_;
  }  

  INLINE Data data() const
  {
    return data_;
  }

  INLINE void set_value(UValue val)
  {
    val_ = val;
  }

  INLINE void set_data(Data data)
  {
    data_ = data;
  }

private:
  UValue val_;
  Data data_;
};

template <typename UValue>
struct SortValue
{
  using uvalue_t = UValue;

  SortValue() {}
  SortValue(UValue val) : val_(val) {}

  INLINE UValue unsigned_value() const
  {
    return val_;
  }  

  INLINE UValue data() const
  {
    return data_;
  }

  INLINE void set_value(UValue val)
  {
    val_ = val;
  }

  INLINE void set_data(UValue data)
  {
    data_ = data;
  }  
  
private:
  union
  {
    UValue val_;
    UValue data_;
  };
};

NAMESPACE_PMT_END