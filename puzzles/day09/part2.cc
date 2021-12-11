#include <array>
#include <cstdint>
#include <functional>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

namespace {

class HeightMap {
 public:
  explicit HeightMap(const std::vector<std::string>& lines)
      : heights_(aoc2021::grid2::Grid<int>::ReadFromDigits(lines)) {}

  int Top3BasinSizeProduct() const {
    std::array<int, 3> top_3_basins{-1, -1, -1};
    for (const aoc2021::grid2::Point point : heights_.Points()) {
      if (IsLowPoint(point)) {
        const int basin_size = BasinSize(point);
        auto min_basin_it =
            std::min_element(top_3_basins.begin(), top_3_basins.end());
        if ((basin_size) > *min_basin_it) *min_basin_it = basin_size;
      }
    }
    return std::accumulate(top_3_basins.begin(), top_3_basins.end(), 1,
                           std::multiplies<int>());
  }

 private:
  bool IsLowPoint(const aoc2021::grid2::Point point) const {
    const int height = heights_[point];
    for (const aoc2021::grid2::Point adjacent :
         heights_.AdjacentCardinal(point)) {
      if (height >= heights_[adjacent]) return false;
    }
    return true;
  }

  int BasinSize(aoc2021::grid2::Point low_point) const {
    absl::flat_hash_set<aoc2021::grid2::Point> basin;
    absl::flat_hash_set<aoc2021::grid2::Point> frontier{low_point};
    while (!frontier.empty()) {
      basin.insert(frontier.begin(), frontier.end());
      absl::flat_hash_set<aoc2021::grid2::Point> new_frontier;
      for (const aoc2021::grid2::Point cell : frontier) {
        for (const aoc2021::grid2::Point adjacent :
             heights_.AdjacentCardinal(cell)) {
          if (heights_[adjacent] != 9 && !basin.contains(adjacent)) {
            new_frontier.insert(adjacent);
          }
        }
      }
      frontier = std::move(new_frontier);
    }
    return basin.size();
  }

  aoc2021::grid2::Grid<int> heights_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  HeightMap map(lines);
  std::cout << map.Top3BasinSizeProduct() << "\n";
  return 0;
}
