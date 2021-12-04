#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

class Board {
 public:
  explicit Board(const std::vector<std::string>& lines) {
    CHECK(lines.size() == 5);
    for (int row = 0; row < 5; ++row) {
      int column = 0;
      for (absl::string_view line_part :
           absl::StrSplit(lines[row], ' ', absl::SkipEmpty())) {
        CHECK(column < 5);
        int number = 0;
        CHECK(absl::SimpleAtoi(line_part, &number));
        numbers_[row][column++] = number;
      }
      CHECK(column == 5);
    }
  }

  bool HasWon() const {
    for (int row = 0; row < 5; ++row) {
      if (WinningRow(row)) {
        return true;
      }
    }
    for (int column = 0; column < 5; ++column) {
      if (WinningColumn(column)) {
        return true;
      }
    }
    return false;
  }

  void Mark(const int number) {
    for (int row = 0; row < 5; ++row) {
      for (int column = 0; column < 5; ++column) {
        if (number == numbers_[row][column]) {
          marked_[row][column] = true;
        }
      }
    }
  }

  int Score(const int called_number) const {
    int unmarked_sum = 0;
    for (int row = 0; row < 5; ++row) {
      for (int column = 0; column < 5; ++column) {
        if (!marked_[row][column]) unmarked_sum += numbers_[row][column];
      }
    }
    return unmarked_sum * called_number;
  }

 private:
  bool WinningRow(int row) const {
    return std::count(marked_[row].begin(), marked_[row].end(), true) == 5;
  }

  bool WinningColumn(int column) const {
    for (int row = 0; row < 5; ++row) {
      if (!marked_[row][column]) return false;
    }
    return true;
  }

  std::array<std::array<int, 5>, 5> numbers_ = {};
  std::array<std::array<bool, 5>, 5> marked_ = {};
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<std::vector<std::string>> segments =
      aoc2021::SplitByEmptyStrings(lines);
  CHECK(segments.size() > 0);

  CHECK(segments.front().size() == 1);
  std::vector<int> calls;
  for (absl::string_view numstr :
       absl::StrSplit(segments.front().front(), ',')) {
    int call = 0;
    CHECK(absl::SimpleAtoi(numstr, &call));
    calls.push_back(call);
  }

  std::vector<Board> boards;
  for (auto segment_iter = segments.begin() + 1; segment_iter != segments.end();
       ++segment_iter) {
    boards.emplace_back(*segment_iter);
  }

  for (const int call : calls) {
    std::vector<Board> remaining_boards;
    for (Board& board : boards) {
      board.Mark(call);
      if (!board.HasWon()) {
        remaining_boards.emplace_back(board);
      } else if (boards.size() == 1) {
        std::cout << board.Score(call) << "\n";
        return 0;
      }
    }
    boards = std::move(remaining_boards);
  }

  CHECK_FAIL();

  return 0;
}
