#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "util/check.h"
#include "util/io.h"

namespace {

bool GammaBit(const std::vector<int>& numbers, int pos) {
  int ones_count = 0;
  const int mask = 1 << pos;
  for (const int num : numbers) {
    if ((mask & num) != 0) ++ones_count;
  }
  return ones_count >= (numbers.size() / 2);
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

  int gamma = 0;
  for (int pos = 0; pos < 12; ++pos) {
    if (GammaBit(numbers, pos)) {
      gamma |= (1 << pos);
    }
  }
  const int epsilon = ~gamma & ((1 << 12) - 1);

  std::cout << (gamma * epsilon) << "\n";

  return 0;
}
