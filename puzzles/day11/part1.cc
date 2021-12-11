#include <iostream>
#include <numeric>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/inlined_vector.h"
#include "util/check.h"
#include "util/io.h"

namespace {

constexpr bool kPart2 = true;

class OctoMap {
 public:
  explicit OctoMap(const std::vector<std::string>& lines) {
    for (const std::string& line : lines) {
      octopi_.emplace_back();
      for (const char c : line) {
        const int energy = c - '0';
        CHECK(energy >= 0);
        CHECK(energy <= 9);
        octopi_.back().emplace_back(energy);
      }
    }
  }

  int Simulate(const int steps) {
    int total_flashes = 0;
    for (int step = 0; step < steps; ++step) {
      total_flashes += Step();
    }
    return total_flashes;
  }

  int SimultaneousFlashStep() {
    const int total_cells = std::accumulate(
        octopi_.begin(), octopi_.end(), 0,
        [](int sum, const std::vector<int>& row) { return sum + row.size(); });
    int step = 1;
    while (Step() != total_cells) ++step;
    return step;
  }

 private:
  using Coords = std::pair<int, int>;

  int Step() {
    std::queue<Coords> pending_flashers;
    for (int x = 0; x < octopi_.size(); ++x) {
      for (int y = 0; y < octopi_[x].size(); ++y) {
        if (++octopi_[x][y] == 10) {
          pending_flashers.emplace(x, y);
        }
      }
    }

    int total_flashes = 0;
    while (!pending_flashers.empty()) {
      ++total_flashes;
      Coords flasher = pending_flashers.front();
      pending_flashers.pop();
      for (Coords adjacent : AdjacentCells(flasher)) {
        if (++octopi_[adjacent.first][adjacent.second] == 10) {
          pending_flashers.emplace(adjacent);
        }
      }
    }

    for (int x = 0; x < octopi_.size(); ++x) {
      for (int y = 0; y < octopi_[x].size(); ++y) {
        if (octopi_[x][y] > 9) {
          octopi_[x][y] = 0;
        }
      }
    }

    return total_flashes;
  }

  absl::InlinedVector<Coords, 8> AdjacentCells(const Coords coords) const {
    absl::InlinedVector<Coords, 8> cells;
    const bool left = coords.first > 0;
    const bool right = coords.first < octopi_.size() - 1;
    const bool up = coords.second > 0;
    const bool down = coords.second < octopi_[coords.first].size() - 1;

    if (left) cells.emplace_back(coords.first - 1, coords.second);
    if (left && up) cells.emplace_back(coords.first - 1, coords.second - 1);
    if (up) cells.emplace_back(coords.first, coords.second - 1);
    if (right && up) cells.emplace_back(coords.first + 1, coords.second - 1);
    if (right) cells.emplace_back(coords.first + 1, coords.second);
    if (right && down) cells.emplace_back(coords.first + 1, coords.second + 1);
    if (down) cells.emplace_back(coords.first, coords.second + 1);
    if (left && down) cells.emplace_back(coords.first - 1, coords.second + 1);

    return cells;
  }

  std::vector<std::vector<int>> octopi_;
};

};  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  OctoMap octo_map(lines);
  const int solution =
      kPart2 ? octo_map.SimultaneousFlashStep() : octo_map.Simulate(100);
  std::cout << solution << "\n";
  return 0;
}
