#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

int64_t FuelForTarget(const std::vector<int64_t>& positions,
                      const int64_t target) {
  int64_t total_fuel = 0;
  for (const int64_t pos : positions) {
    const int64_t distance = std::abs(target - pos);
    total_fuel += (distance * (distance + 1)) / 2;
  }
  return total_fuel;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> input = aoc2021::ReadCommaDelimitedFile(argv[1]);
  std::vector<int64_t> positions;
  for (absl::string_view inputnum : input) {
    int64_t parsed = 0;
    CHECK(absl::SimpleAtoi(inputnum, &parsed));
    positions.push_back(parsed);
  }

  const int64_t min_pos = *std::min_element(positions.begin(), positions.end());
  const int64_t max_pos = *std::max_element(positions.begin(), positions.end());
  int64_t min_fuel = std::numeric_limits<int64_t>::max();
  int64_t target = 0;
  for (int64_t candidate_target = min_pos; candidate_target <= max_pos;
       ++candidate_target) {
    const int64_t fuel = FuelForTarget(positions, candidate_target);
    if (fuel < min_fuel) {
      min_fuel = fuel;
      target = candidate_target;
    }
  }

  std::cout << "target = " << target << "\n";
  std::cout << min_fuel << "\n";

  return 0;
}
