#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

class Signal {
 public:
  Signal() = default;

  explicit Signal(const absl::string_view str) {
    for (char c : str) {
      const int offset = c - 'a';
      CHECK(offset >= 0);
      CHECK(offset < 7);
      segments_ |= (1 << offset);
    }
  }

  int ActiveSegments() const {
    return std::popcount(segments_);
  }

 private:
  uint8_t segments_ = 0;
};

class Display {
 public:
  explicit Display(const absl::string_view rep) {
    std::vector<absl::string_view> parts = absl::StrSplit(rep, ' ');
    CHECK(parts.size() == 15);
    for (int pattern_idx = 0; pattern_idx < 10; ++pattern_idx) {
      all_patterns_[pattern_idx] = Signal(parts[pattern_idx]);
    }
    CHECK(parts[10] == "|");
    for (int output_idx = 0; output_idx < 4; ++output_idx) {
      outputs_[output_idx] = Signal(parts[11 + output_idx]);
    }
  }

  int UniqueDigitsInOutput() const {
    int count = 0;
    for (const Signal& signal : outputs_) {
      switch (signal.ActiveSegments()) {
        case 2:
        case 3:
        case 4:
        case 7:
          ++count;
          break;
        default:
          break;
      }
    }
    return count;
  }

 private:
  std::array<Signal, 10> all_patterns_;
  std::array<Signal, 4> outputs_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<Display> displays;
  for (const std::string& line : lines) {
    displays.emplace_back(line);
  }

  int total_unique_digits = 0;
  for (const Display& display : displays) {
    total_unique_digits += display.UniqueDigitsInOutput();
  }
  std::cout << total_unique_digits << "\n";

  return 0;
}
