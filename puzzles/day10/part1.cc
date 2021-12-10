#include <iostream>
#include <optional>
#include <stack>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

std::optional<int> ErrorScore(absl::string_view line) {
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
          return 3;
        }
        parens.pop();
        break;
      case ']':
        if (parens.empty() || parens.top() != '[') {
          return 57;
        }
        parens.pop();
        break;
      case '}':
        if (parens.empty() || parens.top() != '{') {
          return 1197;
        }
        parens.pop();
        break;
      case '>':
        if (parens.empty() || parens.top() != '<') {
          return 25137;
        }
        parens.pop();
        break;
      default:
        CHECK_FAIL();
    }
  }

  return std::nullopt;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  int error_score = 0;
  for (absl::string_view line : lines) {
    error_score += ErrorScore(line).value_or(0);
  }
  std::cout << error_score << "\n";
  return 0;
}
