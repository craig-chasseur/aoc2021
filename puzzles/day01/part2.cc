#include <iostream>
#include <string>
#include <vector>

#include "util/io.h"

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<int> depths = aoc2021::ParseIntegers(lines);
  int count = 0;
  for (auto iter = depths.begin() + 3; iter != depths.end(); ++iter) {
    const int prev = *(iter - 3) + *(iter - 2) + *(iter - 1);
    const int current = *(iter - 2) + *(iter - 1) + *iter;
    if (current > prev) ++count;
  }
  std::cout << count << "\n";
  return 0;
}
