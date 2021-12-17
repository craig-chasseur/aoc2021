#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <optional>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

namespace {

/*
constexpr int64_t kTargetXMin = 20;
constexpr int64_t kTargetXMax = 30;
constexpr int64_t kTargetYMin = -10;
constexpr int64_t kTargetYMax = -5;
*/

constexpr int64_t kTargetXMin = 57;
constexpr int64_t kTargetXMax = 116;
constexpr int64_t kTargetYMin = -198;
constexpr int64_t kTargetYMax = -148;

struct YInfo {
  int64_t min_step = -1;
  int64_t max_step = -1;
  int64_t peak = 0;
};

std::optional<YInfo> StepsForInitialDy(int64_t dy) {
  YInfo info;
  int64_t y_pos = 0;
  int64_t step = 0;
  while (y_pos >= kTargetYMin) {
    if (y_pos <= kTargetYMax) {
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
  auto test_steps = StepsForInitialDy(3);
  CHECK(test_steps.has_value());

  int64_t highest_peak = 0;
  absl::flat_hash_set<aoc2021::grid2::Vec> valid_velocities;
  for (int64_t candidate_dy = kTargetYMin; candidate_dy < 1000;
       ++candidate_dy) {
    auto step_range = StepsForInitialDy(candidate_dy);
    if (!step_range.has_value()) continue;
    for (int64_t steps = step_range->min_step; steps <= step_range->max_step;
         ++steps) {
      for (int64_t candidate_dx = 0; candidate_dx <= kTargetXMax;
           ++candidate_dx) {
        int64_t x_pos = XAfterSteps(candidate_dx, steps);
        if (x_pos >= kTargetXMin && x_pos <= kTargetXMax) {
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
