#ifndef UTIL_CHECK_H_
#define UTIL_CHECK_H_

#include <cstdlib>
#include <iostream>

#define CHECK(condition)                                  \
  do {                                                    \
    if (!(condition)) {                                   \
      std::cerr << __FILE__ << ":" << __LINE__            \
                << ": CHECK failed: " #condition << "\n"; \
      std::exit(1);                                       \
    }                                                     \
  } while (false)

#define CHECK_FAIL()                                            \
  std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK_FAIL\n"; \
  std::exit(1)

#endif  // UTIL_CHECK_H_
