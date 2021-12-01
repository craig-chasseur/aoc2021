#include <iostream>
#include <string>
#include <vector>

#include "util/io.h"

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<int> depths = aoc2021::ParseIntegers(lines);
  int count = 0;
  for (auto iter = depths.begin() + 1; iter != depths.end(); ++iter) {
    if (*iter > *(iter - 1)) ++count;
  }
  std::cout << count << "\n";
  return 0;
}
