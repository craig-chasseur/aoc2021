#include <iostream>
#include <numeric>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

namespace {

constexpr bool kPart2 = true;

class OctoMap {
 public:
  explicit OctoMap(const std::vector<std::string>& lines)
      : octopi_(aoc2021::grid2::Grid<int>::ReadFromDigits(lines)) {}

  int Simulate(const int steps) {
    int total_flashes = 0;
    for (int step = 0; step < steps; ++step) {
      total_flashes += Step();
    }
    return total_flashes;
  }

  int SimultaneousFlashStep() {
    const int total_cells = octopi_.XSize() * octopi_.YSize();
    int step = 1;
    while (Step() != total_cells) ++step;
    return step;
  }

 private:
  int Step() {
    std::queue<aoc2021::grid2::Point> pending_flashers;
    for (aoc2021::grid2::Point octopus : octopi_.Points()) {
      if (++octopi_[octopus] == 10) pending_flashers.emplace(octopus);
    }

    int total_flashes = 0;
    while (!pending_flashers.empty()) {
      ++total_flashes;
      aoc2021::grid2::Point flasher = pending_flashers.front();
      pending_flashers.pop();
      for (aoc2021::grid2::Point adjacent :
           octopi_.AdjacentWithDiagonal(flasher)) {
        if (++octopi_[adjacent] == 10) pending_flashers.emplace(adjacent);
      }
    }

    for (int& energy : octopi_) {
      if (energy > 9) energy = 0;
    }

    return total_flashes;
  }

  aoc2021::grid2::Grid<int> octopi_;
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
