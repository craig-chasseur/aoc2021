#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

int main(int argc, char** argv) {
  int depth = 0;
  int horizontal = 0;
  int aim = 0;

  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  for (const std::string& line : lines) {
    std::vector<absl::string_view> parts = absl::StrSplit(line, ' ');
    CHECK(parts.size() == 2);
    int value = 0;
    CHECK(absl::SimpleAtoi(parts.back(), &value));
    if (parts.front() == "forward") {
      horizontal += value;
      depth += aim * value;
    } else if (parts.front() == "up") {
      aim -= value;
    } else if (parts.front() == "down") {
      aim += value;
    } else {
      CHECK_FAIL();
    }
  }

  std::cout << (depth * horizontal) << "\n";

  return 0;
}
