#include <cstdint>
#include <functional>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "util/check.h"
#include "util/io.h"

namespace {

class HeightMap {
 public:
  explicit HeightMap(const std::vector<std::string>& lines) {
    for (const std::string& line : lines) {
      heights_.emplace_back();
      for (const char c : line) {
        const int height = c - '0';
        CHECK(height >= 0);
        CHECK(height <= 9);
        heights_.back().emplace_back(height);
      }
    }
  }

  int TotalRisk() const {
    int64_t total_risk = 0;
    for (int x = 0; x < heights_.size(); ++x) {
      for (int y = 0; y < heights_[x].size(); ++y) {
        if (IsLowPoint({x, y})) {
          total_risk += HeightAt({x, y}) + 1;
        }
      }
    }
    return total_risk;
  }

  int Top3BasinSizeProduct() const {
    std::array<int, 3> top_3_basins{-1, -1, -1};
    for (int x = 0; x < heights_.size(); ++x) {
      for (int y = 0; y < heights_[x].size(); ++y) {
        if (IsLowPoint({x, y})) {
          const int basin_size = BasinSize({x, y});
          auto min_basin_it =
              std::min_element(top_3_basins.begin(), top_3_basins.end());
          if ((basin_size) > *min_basin_it) *min_basin_it = basin_size;
        }
      }
    }
    return std::accumulate(top_3_basins.begin(), top_3_basins.end(), 1,
                           std::multiplies<int>());
  }

 private:
  using Coords = std::pair<int, int>;

  int HeightAt(const Coords coords) const {
    return heights_[coords.first][coords.second];
  }

  absl::InlinedVector<Coords, 4> AdjacentCells(const Coords coords) const {
    absl::InlinedVector<Coords, 4> cells;
    if (coords.first > 0) cells.emplace_back(coords.first - 1, coords.second);
    if (coords.first < heights_.size() - 1)
      cells.emplace_back(coords.first + 1, coords.second);
    if (coords.second > 0) cells.emplace_back(coords.first, coords.second - 1);
    if (coords.second < heights_.size() - 1)
      cells.emplace_back(coords.first, coords.second + 1);
    return cells;
  }

  bool IsLowPoint(const Coords coords) const {
    const int height = HeightAt(coords);
    for (const Coords adjacent : AdjacentCells(coords)) {
      if (height >= HeightAt(adjacent)) return false;
    }
    return true;
  }

  int BasinSize(Coords low_point) const {
    absl::flat_hash_set<std::pair<int, int>> basin;
    absl::flat_hash_set<std::pair<int, int>> frontier{low_point};
    while (!frontier.empty()) {
      basin.insert(frontier.begin(), frontier.end());
      absl::flat_hash_set<std::pair<int, int>> new_frontier;
      for (const Coords cell : frontier) {
        for (const Coords adjacent : AdjacentCells(cell)) {
          if (HeightAt(adjacent) != 9 && !basin.contains(adjacent)) {
            new_frontier.insert(adjacent);
          }
        }
      }
      frontier = std::move(new_frontier);
    }
    return basin.size();
  }

  std::vector<std::vector<int>> heights_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  HeightMap map(lines);
  std::cout << map.Top3BasinSizeProduct() << "\n";
  return 0;
}
