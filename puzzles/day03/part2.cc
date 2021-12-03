#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "util/check.h"
#include "util/io.h"

namespace {

constexpr int kDigits = 12;

bool CommonBit(const std::vector<int>& numbers, int pos, bool most_common) {
  int ones_count = 0;
  const int mask = 1 << pos;
  for (const int num : numbers) {
    if ((mask & num) != 0) ++ones_count;
  }
  if (most_common) {
    return ones_count >= std::ceil(numbers.size() / 2.0);
  }
  return ones_count < std::ceil(numbers.size() / 2.0);
}

int Filter(std::vector<int> numbers, bool most_common) {
  for (int pos = kDigits - 1; pos >= 0; --pos) {
    const int match_bit = CommonBit(numbers, pos, most_common) * (1 << pos);
    std::vector<int> filtered_numbers;
    for (int num : numbers) {
      if ((num & (1 << pos)) == match_bit) {
        filtered_numbers.push_back(num);
      }
    }
    if (filtered_numbers.size() == 1) return filtered_numbers.front();
    CHECK(!filtered_numbers.empty());
    numbers = std::move(filtered_numbers);
  }
  CHECK_FAIL();
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<int> numbers;
  for (const std::string& line : lines) {
    int value = 0;
    CHECK(absl::numbers_internal::safe_strto32_base(line, &value, 2));
    numbers.push_back(value);
  }

  const int oxygen = Filter(numbers, true);
  const int scrubber = Filter(numbers, false);

  std::cout << (oxygen * scrubber) << "\n";

  return 0;
}
