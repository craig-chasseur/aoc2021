#include <cstddef>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "util/check.h"
#include "util/io.h"

namespace {

class FloorMap {
 public:
  explicit FloorMap(std::vector<std::string> grid) : grid_(std::move(grid)) {}

  bool Step() {
    std::vector<std::pair<size_t, size_t>> east_movers;
    for (size_t y = 0; y < grid_.size(); ++y) {
      for (size_t x = 0; x < grid_[y].size(); ++x) {
        auto pt = std::make_pair(y, x);
        if (At(pt) == '>' && At(East(pt)) == '.') {
          east_movers.emplace_back(pt);
        }
      }
    }

    for (const auto& pt : east_movers) {
      At(pt) = '.';
      At(East(pt)) = '>';
    }

    std::vector<std::pair<size_t, size_t>> south_movers;
    for (size_t y = 0; y < grid_.size(); ++y) {
      for (size_t x = 0; x < grid_[y].size(); ++x) {
        auto pt = std::make_pair(y, x);
        if (At(pt) == 'v' && At(South(pt)) == '.') {
          south_movers.emplace_back(pt);
        }
      }
    }

    for (const auto& pt : south_movers) {
      At(pt) = '.';
      At(South(pt)) = 'v';
    }

    return !(east_movers.empty() && south_movers.empty());
  }

  void Print() const {
    for (const std::string& row : grid_) {
      std::cout << row << "\n";
    }
    std::cout << "\n";
  }

 private:
  std::pair<size_t, size_t> East(std::pair<size_t, size_t> pt) const {
    if (++pt.second == grid_[pt.first].size()) pt.second = 0;
    return pt;
  }

  std::pair<size_t, size_t> South(std::pair<size_t, size_t> pt) const {
    if (++pt.first == grid_.size()) pt.first = 0;
    return pt;
  }

  char& At(const std::pair<size_t, size_t>& pt) {
    return grid_[pt.first][pt.second];
  }

  std::vector<std::string> grid_;
};

}  // namespace

int main(int argc, char** argv) {
  FloorMap floor_map(aoc2021::ReadLinesFromFile(argv[1]));
  int steps = 0;
  while (floor_map.Step()) {
    ++steps;
  }
  ++steps;
  std::cout << steps << "\n";
  return 0;
}
