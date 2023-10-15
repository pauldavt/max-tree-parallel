#pragma once

#include "../common.h"
#include <cstdlib>

NAMESPACE_PMT

template <typename index_t>
class DynamicStack
{
public:
  DynamicStack(index_t* stack, size_t max_n) :
    stack_(stack), n_(0), max_n_(max_n), is_reallocated_(false)
  {    
  }

  ~DynamicStack()
  {
    if (is_reallocated_)
    {
      free(stack_);
    }
  }

  INLINE void insert(index_t x)
  {
    if (n_ >= max_n_)
    {
      reallocate();
    }

    stack_[n_++] = x;
  }

  INLINE index_t top() const
  {
    return stack_[n_ - 1U];
  }

  INLINE index_t remove()
  {
    return stack_[--n_];
  }

  NO_INLINE void reallocate()
  {
    if (!is_reallocated_)
    {    
      index_t* new_stack = (index_t*)malloc(2U * max_n_ * sizeof(index_t));
      check(new_stack != nullptr);

      for (size_t i = 0; i < max_n_; ++i)
      {
        new_stack[i] = stack_[i];
      }

      stack_ = new_stack;
      max_n_ = 2U * max_n_;
      is_reallocated_ = true;
      return;
    }

    stack_ = (index_t*)realloc(stack_, 2U * max_n_ * sizeof(index_t));
    check(stack_ != nullptr);
    max_n_ = 2U * max_n_;
  }

  INLINE size_t length() const
  {
    return n_;
  }

private:
  index_t* stack_;
  size_t n_;
  size_t max_n_;
  bool is_reallocated_;
};

NAMESPACE_PMT_END