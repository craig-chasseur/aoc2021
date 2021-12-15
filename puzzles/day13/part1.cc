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

using namespace aoc2021::grid2;

namespace {

absl::flat_hash_set<Point> ApplyFold(absl::flat_hash_set<Point> points,
                                     const Line& fold) {
  absl::flat_hash_set<Point> new_points;
  for (const Point& point : points) {
    Vec normal = point - fold;
    if (normal.dx > 0 || normal.dy > 0) {
      new_points.emplace(fold.Reflect(point));
    } else {
      new_points.emplace(point);
    }
  }
  return new_points;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> input = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<std::vector<std::string>> split =
      aoc2021::SplitByEmptyStrings(std::move(input));
  CHECK(split.size() == 2);

  absl::flat_hash_set<Point> points;
  for (const absl::string_view point_str : split.front()) {
    points.emplace(Point::ParseCoordPair(point_str));
  }

  std::vector<Line> folds;
  for (const absl::string_view fold_str : split.back()) {
    Line fold;
    if (fold_str[11] == 'x') {
      fold.fixed_coord = Line::FixedCoord::kX;
    } else {
      CHECK(fold_str[11] == 'y');
      fold.fixed_coord = Line::FixedCoord::kY;
    }

    CHECK(absl::SimpleAtoi(fold_str.substr(13), &fold.fixed_coord_value));
    folds.emplace_back(fold);
  }
  std::cout << "Points after first fold: "
            << ApplyFold(points, folds.front()).size() << "\n\n";

  for (const Line& fold : folds) {
    points = ApplyFold(std::move(points), fold);
  }
  std::cout << RenderPoints(points) << "\n";

  return 0;
}
