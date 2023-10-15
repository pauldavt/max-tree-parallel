#pragma once

#include "../common.h"

NAMESPACE_PMT

template <typename SubType>
struct Edge
{
  SubType a_;
  SubType b_;
};


template<typename SubType>
struct SortableEdgeByStart
{
  using uvalue_t = SubType;

  SubType a_;
  SubType b_;

  INLINE SubType unsigned_value() const
  {
    return a_;
  }
};

template<typename SubType>
struct SortableEdgeByEnd
{
  using uvalue_t = SubType;

  SubType a_;
  SubType b_;

  INLINE SubType unsigned_value() const
  {
    return b_;
  }
};


NAMESPACE_PMT_END
