#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

class CompressedPolymer {
 public:
  CompressedPolymer() = default;

  explicit CompressedPolymer(absl::string_view uncompressed)
      : first_(uncompressed.front()), last_(uncompressed.back()) {
    for (int pos = 1; pos < uncompressed.size(); ++pos) {
      ++pair_counts_[{uncompressed[pos - 1], uncompressed[pos]}];
    }
  }

  CompressedPolymer ApplyRules(
      const absl::flat_hash_map<std::array<char, 2>, char>& rules) const {
    CompressedPolymer next;
    next.first_ = first_;
    next.last_ = last_;

    for (auto [pair, count] : pair_counts_) {
      auto rule_it = rules.find(pair);
      if (rule_it == rules.end()) {
        next.pair_counts_[pair] += count;
        continue;
      }
      next.pair_counts_[{pair[0], rule_it->second}] += count;
      next.pair_counts_[{rule_it->second, pair[1]}] += count;
    }
    return next;
  }

  int64_t Score() const {
    absl::flat_hash_map<char, int64_t> element_counts;
    for (auto [pair, count] : pair_counts_) {
      element_counts[pair[0]] += count;
      element_counts[pair[1]] += count;
    }
    ++element_counts[first_];
    ++element_counts[last_];

    for (auto& [_, count] : element_counts) {
      CHECK(count % 2 == 0);
      count /= 2;
    }

    auto [min_it, max_it] = std::minmax_element(
        element_counts.begin(), element_counts.end(),
        [](const auto l, const auto r) { return l.second < r.second; });
    return max_it->second - min_it->second;
  }

 private:
  absl::flat_hash_map<std::array<char, 2>, int64_t> pair_counts_;
  char first_ = 0;
  char last_ = 0;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  CompressedPolymer polymer(lines.front());

  absl::flat_hash_map<std::array<char, 2>, char> rules;
  for (auto rule_it = lines.begin() + 2; rule_it != lines.end(); ++rule_it) {
    absl::string_view rule_str = *rule_it;
    rules.try_emplace({rule_str[0], rule_str[1]}, rule_str.back());
  }

  for (int i = 0; i < 40; ++i) {
    polymer = polymer.ApplyRules(rules);
  }
  std::cout << polymer.Score() << "\n";

  return 0;
}
