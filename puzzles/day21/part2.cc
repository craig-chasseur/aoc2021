#include <cstdint>
#include <initializer_list>
#include <queue>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "util/check.h"
#include "util/io.h"

namespace {

struct PlayerState {
  int position = 0;
  int score = 0;

  bool operator==(const PlayerState& other) const {
    return position == other.position && score == other.score;
  }

  template <typename H>
  friend H AbslHashValue(H h, const PlayerState& p) {
    return H::combine(std::move(h), p.position, p.score);
  }

  bool Win() const { return score >= 21; }

  PlayerState Move(int spaces) const {
    PlayerState next;
    next.position = (position + spaces) % 10;
    next.score = score + next.position + 1;
    return next;
  }
};

struct GameState {
  PlayerState player1;
  PlayerState player2;
  bool player2_next = false;

  int TotalScore() const { return player1.score + player2.score; }

  bool Win() const { return player1.Win() || player2.Win(); }

  absl::flat_hash_map<GameState, int64_t> ExpandUniverses() const {
    absl::flat_hash_map<GameState, int64_t> expanded;
    for (std::pair<int, int> move_universes :
         {std::make_pair(3, 1), std::make_pair(4, 3), std::make_pair(5, 6),
          std::make_pair(6, 7), std::make_pair(7, 6), std::make_pair(8, 3),
          std::make_pair(9, 1)}) {
      GameState next(*this);
      if (player2_next) {
        next.player2 = player2.Move(move_universes.first);
        next.player2_next = false;
        CHECK(next.TotalScore() > TotalScore());
        CHECK(expanded.emplace(next, move_universes.second).second);
      } else {
        next.player1 = player1.Move(move_universes.first);
        next.player2_next = true;
        CHECK(next.TotalScore() > TotalScore());
        CHECK(expanded.emplace(next, move_universes.second).second);
      }
    }
    return expanded;
  }

  bool operator==(const GameState& other) const {
    return player1 == other.player1 && player2 == other.player2 &&
           player2_next == other.player2_next;
  }

  template <typename H>
  friend H AbslHashValue(H h, const GameState& g) {
    return H::combine(std::move(h), g.player1, g.player2, g.player2_next);
  }
};

struct GreaterTotalScore {
  bool operator()(const GameState& a, const GameState& b) const {
    return a.TotalScore() > b.TotalScore();
  }
};

}  // namespace

int main(int argc, char** argv) {
  GameState game;
  game.player1.position = 3;
  game.player2.position = 2;

  absl::flat_hash_map<GameState, int64_t> game_to_universes{{game, 1}};
  absl::flat_hash_set<GameState> enqueued{game};
  std::priority_queue<GameState, std::vector<GameState>, GreaterTotalScore>
      to_process;
  to_process.emplace(game);
  while (!to_process.empty()) {
    GameState current = to_process.top();
    to_process.pop();
    int64_t universes_in = game_to_universes.at(current);
    for (auto [next_state, expansion_factor] : current.ExpandUniverses()) {
      if (!next_state.Win() && enqueued.emplace(next_state).second) {
        to_process.emplace(next_state);
      }
      game_to_universes[next_state] += universes_in * expansion_factor;
    }
  }

  int64_t player1_wins = 0;
  int64_t player2_wins = 0;
  for (const auto& [game, universes] : game_to_universes) {
    if (game.player1.Win()) player1_wins += universes;
    if (game.player2.Win()) player2_wins += universes;
  }
  std::cout << std::max(player1_wins, player2_wins) << "\n";

  return 0;
}
