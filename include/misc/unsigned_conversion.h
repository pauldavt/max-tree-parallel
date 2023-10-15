#pragma once

#include "../common.h"

#include <climits>

NAMESPACE_PMT

template <typename U>
constexpr U msb()
{
  return U(1U) << ((sizeof(U) * CHAR_BIT) - 1U);
}

template <typename U, typename I>
constexpr U conv_generic_int(I const& i)
{
  U u = i;
  
  return u ^ msb<U>();
}

template <typename I, typename U>
constexpr I undo_conv_generic_int(U u)
{
  return u ^ msb<U>();
}

template <typename U, typename F>
constexpr U conv_generic_float(F const& v)
{
  U const& u = *(U*)(&v);
  
  if ((u & msb<U>()) == msb<U>())
  {
    return ~u;
  }

  return u | msb<U>();
}

template <typename F, typename U>
constexpr F undo_conv_generic_float(U u)
{
  if ((u & msb<U>()) == msb<U>())
  {
    u = u ^ msb<U>();
  }
  else
  {
    u = ~u;
  }

  F* f = ((F*)&u);
  return *f;
}


constexpr uint8_t unsigned_conversion(uint8_t const& u)
{
  return u; 
}

constexpr uint16_t unsigned_conversion(uint16_t const& u)
{
  return u; 
}

constexpr uint32_t unsigned_conversion(uint32_t const& u)
{
  return u; 
}

constexpr uint64_t unsigned_conversion(uint64_t const& u)
{
  return u; 
}

constexpr uint32_t unsigned_conversion(float const& f)
{
  return conv_generic_float<uint32_t>(f);
}

constexpr uint64_t unsigned_conversion(double const& f)
{
  return conv_generic_float<uint64_t>(f);
}

constexpr uint8_t unsigned_conversion(int8_t i)
{
  return conv_generic_int<uint8_t>(i);
}

constexpr uint16_t unsigned_conversion(int16_t i)
{
  return conv_generic_int<uint16_t>(i);
}

constexpr uint32_t unsigned_conversion(int32_t i)
{
  return conv_generic_int<uint32_t>(i);
}

constexpr uint64_t unsigned_conversion(int64_t i)
{
  return conv_generic_int<uint64_t>(i);
}



constexpr void undo_unsigned_conversion(uint32_t u, float& result)
{
  result = undo_conv_generic_float<float>(u);
}

constexpr void undo_unsigned_conversion(uint64_t u, double& result)
{
  result = undo_conv_generic_float<double>(u);
}

constexpr void undo_unsigned_conversion(uint8_t const& u, uint8_t& result)
{
  result = u;
}

constexpr void undo_unsigned_conversion(uint16_t const& u, uint16_t& result)
{
  result = u;
}

constexpr void undo_unsigned_conversion(uint32_t const& u, uint32_t& result)
{
  result = u;
}

constexpr void undo_unsigned_conversion(uint64_t const& u, uint64_t& result)
{
  result = u;
}


constexpr void undo_unsigned_conversion(uint8_t u, int8_t& result)
{
  result = undo_conv_generic_int<int8_t>(u);
}

constexpr void undo_unsigned_conversion(uint16_t u, int16_t& result)
{
  result = undo_conv_generic_int<int16_t>(u);
}

constexpr void undo_unsigned_conversion(uint32_t u, int32_t& result)
{
  result = undo_conv_generic_int<int32_t>(u);
}

constexpr void undo_unsigned_conversion(uint64_t u, int64_t& result)
{
  result = undo_conv_generic_int<int64_t>(u);
}




NAMESPACE_PMT_END