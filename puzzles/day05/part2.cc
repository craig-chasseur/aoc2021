#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "re2/re2.h"
#include "util/check.h"
#include "util/io.h"

namespace {

class Grid {
 public:
  explicit Grid(const int max_x, const int max_y)
      : rep_(max_x + 1, std::vector<int>(max_y + 1, 0)) {}

  int DangerCells() const {
    int count = 0;
    for (const std::vector<int>& row : rep_) {
      count += std::count_if(row.begin(), row.end(),
                             [](int val) { return val >= 2; });
    }
    return count;
  }

 private:
  friend class Line;

  std::vector<std::vector<int>> rep_;
};

class Line {
 public:
  explicit Line(absl::string_view linerep) {
    static re2::LazyRE2 pattern = {R"re((\d+),(\d+) -> (\d+),(\d+))re"};
    CHECK(re2::RE2::FullMatch(linerep, *pattern, &x1_, &y1_, &x2_, &y2_));
  }

  bool HorizontalOrVertical() const { return (x1_ == x2_) || (y1_ == y2_); }

  int MaxX() const { return std::max(x1_, x2_); }

  int MaxY() const { return std::max(y1_, y2_); }

  void MarkOnGrid(Grid& grid) const {
    if (x1_ == x2_) {
      const int min_y = std::min(y1_, y2_);
      const int max_y = std::max(y1_, y2_);
      for (int y = min_y; y <= max_y; ++y) {
        ++grid.rep_[x1_][y];
      }
      return;
    }

    if (y1_ == y2_) {
      const int min_x = std::min(x1_, x2_);
      const int max_x = std::max(x1_, x2_);
      for (int x = min_x; x <= max_x; ++x) {
        ++grid.rep_[x][y1_];
      }
      return;
    }

    const int min_x = std::min(x1_, x2_);
    const int max_x = std::max(x1_, x2_);
    const bool left_1 = min_x == x1_;
    const int start_y = left_1 ? y1_ : y2_;
    const int end_y = left_1 ? y2_ : y1_;
    const bool up = end_y > start_y;
    CHECK(max_x - min_x == (up ? (end_y - start_y) : (start_y - end_y)));
    for (int i = 0; i <= (max_x - min_x); ++i) {
      const int x = min_x + i;
      const int y = start_y + (up ? i : -i);
      ++grid.rep_[x][y];
    }
  }

 private:
  int x1_ = 0;
  int y1_ = 0;
  int x2_ = 0;
  int y2_ = 0;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> inputlines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<Line> lines;
  for (absl::string_view inputline : inputlines) {
    lines.emplace_back(inputline);
  }

  int max_x = 0;
  int max_y = 0;
  for (const Line& line : lines) {
    max_x = std::max(max_x, line.MaxX());
    max_y = std::max(max_y, line.MaxY());
  }
  Grid grid(max_x, max_y);

  for (const Line& line : lines) {
    line.MarkOnGrid(grid);
  }

  std::cout << grid.DangerCells() << "\n";

  return 0;
}
