#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "util/a_star.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

using namespace aoc2021::grid2;

namespace {

Grid<int> BlowupGrid(const Grid<int>& original_grid) {
  Grid<int> new_grid(original_grid.XSize() * 5, original_grid.YSize() * 5);
  const Vec x_stride{.dx = static_cast<int64_t>(original_grid.XSize())};
  const Vec y_stride{.dy = static_cast<int64_t>(original_grid.YSize())};

  for (const Point& original_point : original_grid.Points()) {
    for (int x_repeat = 0; x_repeat < 5; ++x_repeat) {
      for (int y_repeat = 0; y_repeat < 5; ++y_repeat) {
        Point new_point =
            original_point + x_repeat * x_stride + y_repeat * y_stride;
        int new_risk = original_grid[original_point] + x_repeat + y_repeat;
        if (new_risk > 9) {
          new_risk -= 9;
        }
        new_grid[new_point] = new_risk;
      }
    }
  }

  return new_grid;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  Grid<int> risk_grid = BlowupGrid(Grid<int>::ReadFromDigits(lines));

  const Point goal{static_cast<int64_t>(risk_grid.XSize()) - 1,
                   static_cast<int64_t>(risk_grid.YSize()) - 1};
  auto result = aoc2021::AStar(
      Points::kOrigin, goal,
      [&risk_grid](const Point& p) {
        return risk_grid.FilterInRange(p.AdjacentCardinal());
      },
      [&risk_grid](const Point& from, const Point& to) {
        return risk_grid[to];
      });
  CHECK(result.has_value());

  std::cout << result->cost << "\n";

  return 0;
}
