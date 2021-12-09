#include <algorithm>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
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
        if (IsLowPoint(x, y)) {
          total_risk += heights_[x][y] + 1;
        }
      }
    }
    return total_risk;
  }

  int Top3BasinSizeProduct() const {
    std::vector<int> basin_sizes;
    for (int x = 0; x < heights_.size(); ++x) {
      for (int y = 0; y < heights_[x].size(); ++y) {
        if (IsLowPoint(x, y)) {
          basin_sizes.emplace_back(BasinSize(x, y));
        }
      }
    }
    CHECK(basin_sizes.size() >= 3);
    std::partial_sort(basin_sizes.begin(), basin_sizes.begin() + 3,
                      basin_sizes.end(), std::greater<int>());
    return basin_sizes[0] * basin_sizes[1] * basin_sizes[2];
  }

 private:
  bool IsLowPoint(const int x, const int y) const {
    const int height = heights_[x][y];
    if (x > 0 && height >= heights_[x - 1][y]) return false;
    if (x < heights_.size() - 1 && height >= heights_[x + 1][y]) return false;
    if (y > 0 && height >= heights_[x][y - 1]) return false;
    if (y < heights_[x].size() - 1 && height >= heights_[x][y + 1])
      return false;
    return true;
  }

  int BasinSize(const int x, const int y) const {
    absl::flat_hash_set<std::pair<int, int>> basin;
    absl::flat_hash_set<std::pair<int, int>> frontier{{x, y}};
    while (!frontier.empty()) {
      absl::flat_hash_set<std::pair<int, int>> new_frontier;
      for (const std::pair<int, int>& cell : frontier) {
        basin.insert(cell);
        AddAdjacentCellsToSet(cell, new_frontier);
      }

      frontier.clear();
      for (const std::pair<int, int>& cell : new_frontier) {
        if (!basin.contains(cell)) frontier.insert(cell);
      }
    }
    return basin.size();
  }

  void AddAdjacentCellsToSet(
      std::pair<int, int> coords,
      absl::flat_hash_set<std::pair<int, int>>& set) const {
    if (coords.first > 0 && heights_[coords.first - 1][coords.second] != 9) {
      set.emplace(coords.first - 1, coords.second);
    }
    if (coords.first < heights_.size() - 1 &&
        heights_[coords.first + 1][coords.second] != 9) {
      set.emplace(coords.first + 1, coords.second);
    }
    if (coords.second > 0 && heights_[coords.first][coords.second - 1] != 9) {
      set.emplace(coords.first, coords.second - 1);
    }
    if (coords.second < heights_[coords.first].size() - 1 &&
        heights_[coords.first][coords.second + 1] != 9) {
      set.emplace(coords.first, coords.second + 1);
    }
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
