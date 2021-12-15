#include <iostream>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

using namespace aoc2021::grid2;

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);

  Grid<int> risk_grid = Grid<int>::ReadFromDigits(lines);

  absl::flat_hash_map<Point, int> risk_to_point{{Points::kOrigin, 0}};
  absl::flat_hash_set<Point> new_points{Points::kOrigin};
  while (!new_points.empty()) {
    absl::flat_hash_set<Point> next_points;
    for (const Point& point : new_points) {
      for (const Point& adjacent : point.AdjacentCardinal()) {
        if (!risk_grid.InRange(adjacent)) continue;
        const int risk = risk_to_point[point] + risk_grid[adjacent];
        auto [risk_it, inserted] = risk_to_point.try_emplace(adjacent, risk);
        if (inserted || risk_it->second > risk) {
          risk_it->second = risk;
          next_points.emplace(adjacent);
        }
      }
    }
    new_points = std::move(next_points);
  }

  auto destination_it =
      risk_to_point.find({static_cast<int64_t>(risk_grid.XSize()) - 1,
                          static_cast<int64_t>(risk_grid.YSize()) - 1});
  CHECK(destination_it != risk_to_point.end());
  std::cout << destination_it->second << "\n";

  return 0;
}
