#include <algorithm>
#include <cstdint>
#include <iostream>
#include <optional>
#include <stack>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

std::optional<int64_t> CompletionScore(absl::string_view line) {
  std::stack<char> parens;
  for (char c : line) {
    switch (c) {
      case '(':
      case '[':
      case '{':
      case '<':
        parens.push(c);
        break;
      case ')':
        if (parens.empty() || parens.top() != '(') {
          return std::nullopt;
        }
        parens.pop();
        break;
      case ']':
        if (parens.empty() || parens.top() != '[') {
          return std::nullopt;
        }
        parens.pop();
        break;
      case '}':
        if (parens.empty() || parens.top() != '{') {
          return std::nullopt;
        }
        parens.pop();
        break;
      case '>':
        if (parens.empty() || parens.top() != '<') {
          return std::nullopt;
        }
        parens.pop();
        break;
      default:
        CHECK_FAIL();
    }
  }

  int64_t completion_score = 0;
  while (!parens.empty()) {
    completion_score *= 5;
    switch (parens.top()) {
      case '(':
        completion_score += 1;
        break;
      case '[':
        completion_score += 2;
        break;
      case '{':
        completion_score += 3;
        break;
      case '<':
        completion_score += 4;
        break;
      default:
        CHECK_FAIL();
    }
    parens.pop();
  }
  return completion_score;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);

  std::vector<int64_t> completion_scores;
  for (absl::string_view line : lines) {
    std::optional<int64_t> completion = CompletionScore(line);
    if (!completion.has_value()) continue;
    completion_scores.push_back(*completion);
  }
  CHECK(completion_scores.size() % 2 == 1);

  std::sort(completion_scores.begin(), completion_scores.end());
  std::cout << completion_scores[completion_scores.size() / 2] << "\n";

  return 0;
}
