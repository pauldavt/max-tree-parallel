#pragma once

#include <iostream>
#include <sstream>
#include <iomanip>
#include "../common.h"

NAMESPACE_PMT

extern bool show_info;

enum warn_t
{
  INFO,
  WARNING,
  ERROR
};

extern std::string warnings[];

#define PMT_STR_(...) #__VA_ARGS__
#define PMT_STR(...) PMT_STR_(__VA_ARGS__)
#define PMT_FILE_LINE __FILE__ ":" PMT_STR(__LINE__)
#define pmt_logger(level, msg) \
  std::cout << "[" << pmt::warnings[level] << "] " << msg << " [" << PMT_FILE_LINE << "]\n";
#define pmt_logger_no_loc(level, msg) \
  std::cout << "[" << pmt::warnings[level] << "] " << msg << '\n';

#define info(m) if (pmt::show_info) { pmt_logger(pmt::INFO, m) };
#define warn(m) pmt_logger(pmt::WARNING, m);
#define err(m){ pmt_logger(pmt::ERROR, m); exit(1); }

#define check(...) if (!(__VA_ARGS__)) \
{ \
  std::ostringstream ss; \
  ss << "false assertion: " << PMT_STR(__VA_ARGS__); \
  err(ss.str()); \
}

#ifdef PMT_DEBUG
#  define debug(cond) check(cond)
#else
#  define debug(cond)
#endif

#define out(...) \
{ \
  std::ostringstream os; \
  os << __VA_ARGS__; \
  pmt_logger_no_loc(pmt::INFO, os.str()); \
}

NAMESPACE_PMT_END