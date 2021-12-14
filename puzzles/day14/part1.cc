#include <algorithm>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

std::string ApplyRules(
    const std::string& polymer,
    const absl::flat_hash_map<absl::string_view, char>& rules) {
  std::string new_polymer(1, polymer.front());
  for (int pos = 1; pos < polymer.size(); ++pos) {
    auto rule_it = rules.find(polymer.substr(pos - 1, 2));
    if (rule_it != rules.end()) {
      new_polymer.push_back(rule_it->second);
    }
    new_polymer.push_back(polymer[pos]);
  }
  return new_polymer;
}

absl::flat_hash_map<char, int> ElementCounts(const std::string& polymer) {
  absl::flat_hash_map<char, int> counts;
  for (const char element : polymer) {
    ++counts[element];
  }
  return counts;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::string polymer = lines.front();

  absl::flat_hash_map<absl::string_view, char> rules;
  for (auto rule_it = lines.begin() + 2; rule_it != lines.end(); ++rule_it) {
    absl::string_view rule_str = *rule_it;
    rules.try_emplace(rule_str.substr(0, 2), rule_str.back());
  }

  for (int i = 0; i < 10; ++i) {
    polymer = ApplyRules(polymer, rules);
  }

  const auto element_counts = ElementCounts(polymer);
  auto [min_it, max_it] = std::minmax_element(
      element_counts.begin(), element_counts.end(),
      [](const auto l, const auto r) { return l.second < r.second; });
  std::cout << (max_it->second - min_it->second) << "\n";

  return 0;
}
