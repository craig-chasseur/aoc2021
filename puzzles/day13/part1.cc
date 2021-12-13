#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

namespace {

struct Fold {
  bool x = false;
  int64_t value = 0;
};

absl::flat_hash_set<aoc2021::grid2::Point> ApplyFold(
    absl::flat_hash_set<aoc2021::grid2::Point> points, const Fold& fold) {
  absl::flat_hash_set<aoc2021::grid2::Point> new_points;
  for (const aoc2021::grid2::Point& point : points) {
    if (fold.x) {
      CHECK(point.x != fold.value);
      if (point.x < fold.value) {
        new_points.emplace(point);
      } else {
        const int64_t dist = point.x - fold.value;
        aoc2021::grid2::Point folded(point);
        folded.x = fold.value - dist;
        new_points.emplace(folded);
      }
    } else {
      CHECK(point.y != fold.value);
      if (point.y < fold.value) {
        new_points.emplace(point);
      } else {
        const int64_t dist = point.y - fold.value;
        aoc2021::grid2::Point folded(point);
        folded.y = fold.value - dist;
        new_points.emplace(folded);
      }
    }
  }
  return new_points;
}

void RenderPoints(const absl::flat_hash_set<aoc2021::grid2::Point>& points) {
  int64_t max_x = 0;
  int64_t max_y = 0;

  for (const aoc2021::grid2::Point& point : points) {
    max_x = std::max(max_x, point.x);
    max_y = std::max(max_y, point.y);
  }

  std::vector<std::string> view(max_y + 1, std::string(max_x + 1, ' '));
  for (const aoc2021::grid2::Point& point : points) {
    view[point.y][point.x] = '#';
  }

  for (const std::string& line : view) {
    std::cout << line << "\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> input = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<std::vector<std::string>> split =
      aoc2021::SplitByEmptyStrings(std::move(input));
  CHECK(split.size() == 2);

  absl::flat_hash_set<aoc2021::grid2::Point> points;
  for (const std::string& point_str : split.front()) {
    std::vector<absl::string_view> coords = absl::StrSplit(point_str, ',');
    CHECK(coords.size() == 2);
    aoc2021::grid2::Point point;
    CHECK(absl::SimpleAtoi(coords.front(), &point.x));
    CHECK(absl::SimpleAtoi(coords.back(), &point.y));
    points.emplace(point);
  }

  std::vector<Fold> folds;
  for (const std::string& fold_str : split.back()) {
    Fold fold;
    if (fold_str[11] == 'x') {
      fold.x = true;
    } else {
      CHECK(fold_str[11] == 'y');
      fold.x = false;
    }

    CHECK(absl::SimpleAtoi(fold_str.substr(13), &fold.value));
    folds.emplace_back(fold);
  }
  std::cout << "Points after first fold: "
            << ApplyFold(points, folds.front()).size() << "\n\n";

  for (const Fold& fold : folds) {
    points = ApplyFold(std::move(points), fold);
  }
  RenderPoints(points);

  return 0;
}
