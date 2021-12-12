#include <cstdint>
#include <string>
#include <vector>

#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

namespace {

class HeightMap {
 public:
  explicit HeightMap(const std::vector<std::string>& lines)
      : heights_(aoc2021::grid2::Grid<int>::ReadFromDigits(lines)) {}

  int TotalRisk() const {
    int64_t total_risk = 0;
    for (const aoc2021::grid2::Point point : heights_.Points()) {
      if (IsLowPoint(point)) total_risk += heights_[point] + 1;
    }
    return total_risk;
  }

 private:
  bool IsLowPoint(const aoc2021::grid2::Point point) const {
    const int height = heights_[point];
    for (const aoc2021::grid2::Point adjacent :
         heights_.FilterInRange(point.AdjacentCardinal())) {
      if (height >= heights_[adjacent]) return false;
    }
    return true;
  }

  aoc2021::grid2::Grid<int> heights_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  HeightMap map(lines);
  std::cout << map.TotalRisk() << "\n";
  return 0;
}
