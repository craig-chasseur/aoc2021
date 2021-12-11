#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "re2/re2.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

namespace {

int DangerCells(const aoc2021::grid2::Grid<int>& grid) {
  return std::count_if(grid.begin(), grid.end(),
                       [](int val) { return val >= 2; });
}

class Line {
 public:
  explicit Line(absl::string_view linerep) {
    static re2::LazyRE2 pattern = {R"re((\d+),(\d+) -> (\d+),(\d+))re"};
    CHECK(
        re2::RE2::FullMatch(linerep, *pattern, &p1_.x, &p1_.y, &p2_.x, &p2_.y));
  }
  int64_t MaxX() const { return std::max(p1_.x, p2_.x); }

  int64_t MaxY() const { return std::max(p1_.y, p2_.y); }

  void MarkOnGrid(aoc2021::grid2::Grid<int>& grid) const {
    const aoc2021::grid2::Vec step = (p2_ - p1_).UnitStep();
    for (aoc2021::grid2::Point p = p1_; p != p2_ + step; p += step) {
      ++grid[p];
    }
  }

 private:
  aoc2021::grid2::Point p1_, p2_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> inputlines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<Line> lines;
  for (absl::string_view inputline : inputlines) {
    lines.emplace_back(inputline);
  }

  int64_t max_x = 0;
  int64_t max_y = 0;
  for (const Line& line : lines) {
    max_x = std::max(max_x, line.MaxX());
    max_y = std::max(max_y, line.MaxY());
  }
  aoc2021::grid2::Grid grid(max_x + 1, max_y + 1);

  for (const Line& line : lines) {
    line.MarkOnGrid(grid);
  }

  std::cout << DangerCells(grid) << "\n";

  return 0;
}
