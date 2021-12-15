#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
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

struct PointAndRisk {
  PointAndRisk(Point point_in, int risk_in) : point(point_in), risk(risk_in) {}

  Point point;
  int risk = 0;
};

struct GreaterRisk {
  bool operator()(const PointAndRisk& a, const PointAndRisk& b) {
    return a.risk > b.risk;
  }
};

int BestPath(const Grid<int>& risk_grid) {
  const Point goal{static_cast<int64_t>(risk_grid.XSize()) - 1,
                   static_cast<int64_t>(risk_grid.YSize()) - 1};

  absl::flat_hash_map<Point, int> path_risks{{Points::kOrigin, 0}};
  std::priority_queue<PointAndRisk, std::vector<PointAndRisk>, GreaterRisk>
      point_queue;
  point_queue.emplace(Points::kOrigin, 0);
  while (!point_queue.empty()) {
    PointAndRisk current = point_queue.top();
    point_queue.pop();
    if (current.point == goal) return current.risk;

    // It's possible that in the time since a path to this point was put on the
    // queue, a shorter path to the point was already found and processed first.
    // Skip re-processing this point if that happened.
    if (current.risk > path_risks[current.point]) continue;

    for (const Point& adjacent : current.point.AdjacentCardinal()) {
      if (!risk_grid.InRange(adjacent)) continue;
      const int tentative_risk = current.risk + risk_grid[adjacent];
      auto [risk_it, inserted] =
          path_risks.try_emplace(adjacent, tentative_risk);
      if (inserted || tentative_risk < risk_it->second) {
        risk_it->second = tentative_risk;
        point_queue.emplace(adjacent, tentative_risk);
      }
    }
  }

  CHECK_FAIL();
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  Grid<int> risk_grid = BlowupGrid(Grid<int>::ReadFromDigits(lines));
  std::cout << BestPath(risk_grid) << "\n";

  return 0;
}
