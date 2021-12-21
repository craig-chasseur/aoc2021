#include <initializer_list>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "re2/re2.h"
#include "util/check.h"
#include "util/io.h"

namespace {

int ParsePosition(absl::string_view line) {
  static re2::LazyRE2 pattern{R"re(starting position: (\d+))re"};
  int position = 0;
  CHECK(re2::RE2::PartialMatch(line, *pattern, &position));
  CHECK(position >= 1 && position <= 10);
  return position - 1;
}

class Dice {
 public:
  Dice() = default;

  int Roll() {
    ++rolled_;
    const int roll = next_roll_;
    if (++next_roll_ > 100) {
      next_roll_ = 1;
    }
    return roll;
  }

  int TimesRolled() const { return rolled_; }

 private:
  int next_roll_ = 1;
  int rolled_ = 0;
};

struct PlayerState {
  int position = 0;
  int score = 0;
};

class GameState {
 public:
  explicit GameState(std::initializer_list<int> player_starting_positions) {
    for (int pos : player_starting_positions) {
      players_.emplace_back(PlayerState{.position = pos});
    }
  }

  int PlayUntilWin() {
    for (;;) {
      for (int player_idx = 0; player_idx < 2; ++player_idx) {
        if (TakeTurn(players_[player_idx])) {
          return players_[player_idx == 0 ? 1 : 0].score * dice_.TimesRolled();
        }
      }
    }
  }

 private:
  bool TakeTurn(PlayerState& player) {
    int next_position =
        (dice_.Roll() + dice_.Roll() + dice_.Roll() + player.position) % 10;
    player.position = next_position;
    player.score += next_position + 1;
    return player.score >= 1000;
  }

  std::vector<PlayerState> players_;
  Dice dice_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  CHECK(lines.size() == 2);
  GameState game({ParsePosition(lines.front()), ParsePosition(lines.back())});
  std::cout << game.PlayUntilWin() << "\n";
  return 0;
}
