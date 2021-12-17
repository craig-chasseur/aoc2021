#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <optional>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "re2/re2.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

namespace {

struct YInfo {
  int64_t min_step = -1;
  int64_t max_step = -1;
  int64_t peak = 0;
};

std::optional<YInfo> StepsForInitialDy(int64_t dy, const int64_t y_min,
                                       const int64_t y_max) {
  YInfo info;
  int64_t y_pos = 0;
  int64_t step = 0;
  while (y_pos >= y_min) {
    if (y_pos <= y_max) {
      if (info.min_step == -1) info.min_step = step;
      info.max_step = step;
    }
    info.peak = std::max(info.peak, y_pos);
    y_pos += dy;
    --dy;
    ++step;
  }

  if (info.min_step >= 0 and info.max_step >= 0) {
    return info;
  }
  return std::nullopt;
}

int64_t XAfterSteps(int64_t dx, int64_t step) {
  if (step >= dx) {
    return (dx + 1) * dx / 2;
  }
  return (dx + (dx - step + 1)) * step / 2;
}

}  // namespace

int main(int argc, char** argv) {
  const std::string target_str = aoc2021::ReadFile(argv[1]);

  int64_t x_min = 0;
  int64_t x_max = 0;
  CHECK(RE2::PartialMatch(target_str, R"re(x=(-?\d+)\.\.(-?\d+))re", &x_min,
                          &x_max));

  int64_t y_min = 0;
  int64_t y_max = 0;
  CHECK(RE2::PartialMatch(target_str, R"re(y=(-?\d+)\.\.(-?\d+))re", &y_min,
                          &y_max));

  int64_t highest_peak = 0;
  absl::flat_hash_set<aoc2021::grid2::Vec> valid_velocities;
  for (int64_t candidate_dy = y_min; candidate_dy < 1000; ++candidate_dy) {
    auto step_range = StepsForInitialDy(candidate_dy, y_min, y_max);
    if (!step_range.has_value()) continue;
    for (int64_t steps = step_range->min_step; steps <= step_range->max_step;
         ++steps) {
      for (int64_t candidate_dx = 0; candidate_dx <= x_max; ++candidate_dx) {
        int64_t x_pos = XAfterSteps(candidate_dx, steps);
        if (x_pos >= x_min && x_pos <= x_max) {
          highest_peak = std::max(highest_peak, step_range->peak);
          valid_velocities.emplace(
              aoc2021::grid2::Vec{candidate_dx, candidate_dy});
        }
      }
    }
  }

  std::cout << "Highest peak (part 1): " << highest_peak << "\n";
  std::cout << "Total valid velocities (part 2): " << valid_velocities.size()
            << "\n";
  return 0;
}
