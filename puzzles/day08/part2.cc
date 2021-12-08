#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
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
    std::array<Signal, 10> known_signals;
    std::bitset<10> known;
    std::vector<Signal> seg5, seg6;
    for (const Signal& pattern : all_patterns_) {
      switch (pattern.ActiveSegments()) {
        case 2:
          known_signals[1] = pattern;
          known.set(1);
          break;
        case 3:
          known_signals[7] = pattern;
          known.set(7);
          break;
        case 4:
          known_signals[4] = pattern;
          known.set(4);
          break;
        case 5:
          seg5.push_back(pattern);
          break;
        case 6:
          seg6.push_back(pattern);
          break;
        case 7:
          known_signals[8] = pattern;
          known.set(8);
          break;
        default:
          break;
      }
    }
    CHECK(known.count() == 4);
    CHECK(seg5.size() == 3);
    CHECK(seg6.size() == 3);

    const Signal sig_abcdf = known_signals[4] + known_signals[7];
    for (const Signal& pattern : seg6) {
      if (pattern.Contains(sig_abcdf)) {
        known_signals[9] = pattern;
        known.set(9);
      } else if (pattern.Contains(known_signals[1])) {
        known_signals[0] = pattern;
        known.set(0);
      } else {
        known_signals[6] = pattern;
        known.set(6);
      }
    }
    CHECK(known.count() == 7);

    const Signal sig_c = known_signals[8] - known_signals[6];
    for (const Signal& pattern : seg5) {
      if (pattern.Contains(known_signals[1])) {
        known_signals[3] = pattern;
        known.set(3);
      } else if (pattern.Contains(sig_c)) {
        known_signals[2] = pattern;
        known.set(2);
      } else {
        known_signals[5] = pattern;
        known.set(5);
      }
    }
    CHECK(known.count() == 10);

    absl::flat_hash_map<Signal, int> lookup_table;
    for (int i = 0; i < 10; ++i) {
      lookup_table[known_signals[i]] = i;
    }

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
