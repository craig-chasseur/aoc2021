#ifndef UTIL_IO_H_
#define UTIL_IO_H_

#include <cstddef>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "absl/strings/string_view.h"
#include "util/check.h"

namespace aoc2021 {

// Returns the whole contents of the text file specified by `filename`.
std::string ReadFile(const char* filename);

// Returns each of the lines in the text file specified by `filename`.
std::vector<std::string> ReadLinesFromFile(const char* filename);

// Returns comma-delimited strings from `filename`.
std::vector<std::string> ReadCommaDelimitedFile(const char* filename);

template <typename IntType = int>
std::vector<IntType> ParseIntegers(const std::vector<std::string>& strings) {
  std::vector<IntType> integers(strings.size(), 0);
  for (std::size_t idx = 0; idx < strings.size(); ++idx) {
    CHECK(absl::SimpleAtoi(strings[idx], &(integers[idx])));
  }
  return integers;
}

template <typename IntType = int>
IntType ParseBinary(const absl::string_view str) {
  IntType parsed = 0;
  for (const char c : str) {
    parsed <<= 1;
    switch (c) {
      case '0':
        break;
      case '1':
        parsed |= 1;
        break;
      default:
        CHECK_FAIL();
    }
  }
  return parsed;
}

template <typename IntType = int>
std::vector<IntType> ParseMultiBinary(const std::vector<std::string>& strings) {
  std::vector<IntType> integers(strings.size(), 0);
  for (std::size_t idx = 0; idx < strings.size(); ++idx) {
    integers[idx] = ParseBinary<IntType>(strings[idx]);
  }
  return integers;
}

std::vector<std::vector<std::string>> SplitByEmptyStrings(
    std::vector<std::string> strs);

}  // namespace aoc2021

#endif  // UTIL_IO_H_
