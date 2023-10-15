#pragma once

#include "../common.h"

NAMESPACE_PMT

template <typename Sub>
struct Range
{
  using sub_t = Sub;

  sub_t begin_;
  sub_t end_;
};

NAMESPACE_PMT_END