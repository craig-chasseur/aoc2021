#include <cstdint>
#include <string>
#include <vector>

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

 private:
  bool IsLowPoint(const int x, const int y) const {
    const int height = heights_[x][y];
    if (x > 0 && height >= heights_[x-1][y]) return false;
    if (x < heights_.size() - 1 && height >= heights_[x+1][y]) return false;
    if (y > 0 && height >= heights_[x][y-1]) return false;
    if (y < heights_[x].size() - 1 && height >= heights_[x][y+1]) return false;
    return true;
  }

  std::vector<std::vector<int>> heights_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  HeightMap map(lines);
  std::cout << map.TotalRisk() << "\n";
  return 0;
}
