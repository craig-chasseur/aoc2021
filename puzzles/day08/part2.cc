#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
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

  bool operator==(const Signal& other) const {
    return segments_ == other.segments_;
  }

  bool operator!=(const Signal& other) const {
    return segments_ != other.segments_;
  }

  template <typename H>
  friend H AbslHashValue(H h, const Signal& signal) {
    return H::combine(std::move(h), signal.segments_);
  }

  Signal operator+(const Signal& other) const {
    Signal result(*this);
    result.segments_ |= other.segments_;
    return result;
  }

  Signal operator-(const Signal& other) const {
    Signal result(*this);
    result.segments_ &= ~other.segments_;
    return result;
  }

  bool Contains(const Signal& other) const {
    return (segments_ & other.segments_) == other.segments_;
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

  int Solve() const {
    Signal sig_1, sig_4, sig_7, sig_8;
    int found = 0;
    for (const Signal& pattern : all_patterns_) {
      switch (pattern.ActiveSegments()) {
        case 2:
          sig_1 = pattern;
          ++found;
          break;
        case 3:
          sig_7 = pattern;
          ++found;
          break;
        case 4:
          sig_4 = pattern;
          ++found;
          break;
        case 7:
          sig_8 = pattern;
          ++found;
          break;
        default:
          break;
      }
    }
    CHECK(found = 4);

    Signal sig_a = sig_7 - sig_1;
    Signal sig_bd = sig_4 - sig_1;
    Signal sig_abd = sig_a + sig_bd;
    Signal sig_abcdf = sig_abd + sig_1;

    Signal sig_9;
    found = 0;
    for (const Signal& pattern : all_patterns_) {
      if (pattern.Contains(sig_abcdf) && pattern != sig_8) {
        sig_9 = pattern;
        ++found;
        break;
      }
    }
    CHECK(found == 1);
    Signal sig_e = sig_8 - sig_9;

    Signal sig_6, sig_0;
    found = 0;
    for (const Signal& pattern : all_patterns_) {
      if (pattern != sig_9 && pattern.ActiveSegments() == 6) {
        if (pattern.Contains(sig_1)) {
          sig_0 = pattern;
        } else {
          sig_6 = pattern;
        }
        ++found;
      }
    }
    CHECK(found == 2);
    Signal sig_c = sig_8 - sig_6;
    Signal sig_d = sig_8 - sig_0;
    Signal sig_f = sig_1 - sig_c;

    Signal sig_2, sig_3, sig_5;
    found = 0;
    for (const Signal& pattern : all_patterns_) {
      if (pattern.ActiveSegments() == 5) {
        if (pattern.Contains(sig_1)) {
          sig_3 = pattern;
        } else if (pattern.Contains(sig_c)) {
          sig_2 = pattern;
        } else {
          CHECK(pattern.Contains(sig_f));
          sig_5 = pattern;
        }
        ++found;
      }
    }
    CHECK(found == 3);

    const absl::flat_hash_map<Signal, int> lookup_table {
      {sig_0, 0},
      {sig_1, 1},
      {sig_2, 2},
      {sig_3, 3},
      {sig_4, 4},
      {sig_5, 5},
      {sig_6, 6},
      {sig_7, 7},
      {sig_8, 8},
      {sig_9, 9},
    };

    int display_value = 0;
    for (const Signal& output_signal : outputs_) {
      const auto value_it = lookup_table.find(output_signal);
      CHECK(value_it != lookup_table.end());
      display_value = display_value * 10 + value_it->second;
    }
    return display_value;
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

  int total_outputs = 0;
  for (const Display& display : displays) {
    total_outputs += display.Solve();
  }
  std::cout << total_outputs << "\n";

  return 0;
}
