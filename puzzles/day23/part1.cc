#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "util/a_star.h"
#include "util/check.h"
#include "util/io.h"

namespace {

enum class Cell : char { kEmpty = 0, kA, kB, kC, kD };

Cell ParseCell(const char c) {
  switch (c) {
    case '.':
      return Cell::kEmpty;
    case 'A':
      return Cell::kA;
    case 'B':
      return Cell::kB;
    case 'C':
      return Cell::kC;
    case 'D':
      return Cell::kD;
    default:
      CHECK_FAIL();
  }
}

int RoomPosition(const Cell cell) {
  switch (cell) {
    case Cell::kEmpty:
      CHECK_FAIL();
    case Cell::kA:
      return 2;
    case Cell::kB:
      return 4;
    case Cell::kC:
      return 6;
    case Cell::kD:
      return 8;
  }
}

int64_t Energy(const Cell cell) {
  switch (cell) {
    case Cell::kEmpty:
      CHECK_FAIL();
    case Cell::kA:
      return 1;
    case Cell::kB:
      return 10;
    case Cell::kC:
      return 100;
    case Cell::kD:
      return 1000;
  }
}

struct State {
  std::array<Cell, 11> hallway;
  std::array<Cell, 2> room_a;
  std::array<Cell, 2> room_b;
  std::array<Cell, 2> room_c;
  std::array<Cell, 2> room_d;

  std::vector<std::pair<State, int64_t>> SuccessorStates() const {
    std::vector<std::pair<State, int64_t>> successors;

    // Moving Amphipods -> Final Rooms.
    for (int hallway_pos = 0; hallway_pos < 11; ++hallway_pos) {
      const Cell source_cell = hallway[hallway_pos];
      if (source_cell == Cell::kEmpty) continue;
      int destination = RoomPosition(source_cell);
      bool path_clear = true;
      if (hallway_pos < destination) {
        for (int pos = hallway_pos + 1; pos <= destination; ++pos) {
          if (hallway[pos] != Cell::kEmpty) {
            path_clear = false;
            break;
          }
        }
      } else {
        CHECK(hallway_pos > destination);
        for (int pos = hallway_pos - 1; pos >= destination; --pos) {
          if (hallway[pos] != Cell::kEmpty) {
            path_clear = false;
            break;
          }
        }
      }
      if (!path_clear) continue;
      const int hallway_path_length = std::abs(hallway_pos - destination);
      const std::array<Cell, 2>& dest = DestinationRoom(source_cell);
      if (dest[0] != Cell::kEmpty) continue;
      if (dest[1] == Cell::kEmpty) {
        State successor = *this;
        successor.hallway[hallway_pos] = Cell::kEmpty;
        successor.DestinationRoom(source_cell)[1] = source_cell;
        const int64_t cost = (hallway_path_length + 2) * Energy(source_cell);
        successors.emplace_back(std::move(successor), cost);
      } else if (dest[1] == source_cell) {
        State successor = *this;
        successor.hallway[hallway_pos] = Cell::kEmpty;
        successor.DestinationRoom(source_cell)[0] = source_cell;
        const int64_t cost = (hallway_path_length + 1) * Energy(source_cell);
        successors.emplace_back(std::move(successor), cost);
      }
    }

    // Moving from rooms -> hallway.
    for (const Cell room_type : {Cell::kA, Cell::kB, Cell::kC, Cell::kD}) {
      const std::array<Cell, 2>& source_room = DestinationRoom(room_type);
      if (source_room[0] == room_type) continue;
      int initial_room_pos = 0;
      if (source_room[0] != Cell::kEmpty) {
        initial_room_pos = 0;
      } else {
        if (source_room[1] == Cell::kEmpty || source_room[1] == room_type)
          continue;
        initial_room_pos = 1;
      }

      const int starting_room_position = RoomPosition(room_type);
      CHECK(hallway[starting_room_position] == Cell::kEmpty);
      for (int final_hallway_pos = starting_room_position;
           final_hallway_pos >= 0; --final_hallway_pos) {
        if (hallway[final_hallway_pos] != Cell::kEmpty) break;
        if (final_hallway_pos == 2 || final_hallway_pos == 4 ||
            final_hallway_pos == 6 || final_hallway_pos == 8) {
          continue;
        }
        State successor = *this;
        successor.DestinationRoom(room_type)[initial_room_pos] = Cell::kEmpty;
        successor.hallway[final_hallway_pos] = source_room[initial_room_pos];
        const int64_t cost = (initial_room_pos + 1 + starting_room_position -
                              final_hallway_pos) *
                             Energy(source_room[initial_room_pos]);
        successors.emplace_back(std::move(successor), cost);
      }
      for (int final_hallway_pos = starting_room_position;
           final_hallway_pos < 11; ++final_hallway_pos) {
        if (hallway[final_hallway_pos] != Cell::kEmpty) break;
        if (final_hallway_pos == 2 || final_hallway_pos == 4 ||
            final_hallway_pos == 6 || final_hallway_pos == 8) {
          continue;
        }
        State successor = *this;
        successor.DestinationRoom(room_type)[initial_room_pos] = Cell::kEmpty;
        successor.hallway[final_hallway_pos] = source_room[initial_room_pos];
        const int64_t cost = (initial_room_pos + 1 + final_hallway_pos -
                              starting_room_position) *
                             Energy(source_room[initial_room_pos]);
        successors.emplace_back(std::move(successor), cost);
      }
    }

    return successors;
  }

  bool operator==(const State& other) const {
    return hallway == other.hallway && room_a == other.room_a &&
           room_b == other.room_b && room_c == other.room_c &&
           room_d == other.room_d;
  }

  template <typename H>
  friend H AbslHashValue(H h, const State& state) {
    return H::combine(std::move(h), state.hallway, state.room_a, state.room_b,
                      state.room_c, state.room_d);
  }

 private:
  const std::array<Cell, 2>& DestinationRoom(const Cell cell) const {
    switch (cell) {
      case Cell::kEmpty:
        CHECK_FAIL();
      case Cell::kA:
        return room_a;
      case Cell::kB:
        return room_b;
      case Cell::kC:
        return room_c;
      case Cell::kD:
        return room_d;
    }
  }

  std::array<Cell, 2>& DestinationRoom(const Cell cell) {
    switch (cell) {
      case Cell::kEmpty:
        CHECK_FAIL();
      case Cell::kA:
        return room_a;
      case Cell::kB:
        return room_b;
      case Cell::kC:
        return room_c;
      case Cell::kD:
        return room_d;
    }
  }
};

constexpr State kGoalState{.hallway = {},
                           .room_a = {Cell::kA, Cell::kA},
                           .room_b = {Cell::kB, Cell::kB},
                           .room_c = {Cell::kC, Cell::kC},
                           .room_d = {Cell::kD, Cell::kD}};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  CHECK(lines.size() == 5);
  State initial_state{
      .hallway = {},
      .room_a = {ParseCell(lines[2][3]), ParseCell(lines[3][3])},
      .room_b = {ParseCell(lines[2][5]), ParseCell(lines[3][5])},
      .room_c = {ParseCell(lines[2][7]), ParseCell(lines[3][7])},
      .room_d = {ParseCell(lines[2][9]), ParseCell(lines[3][9])}};

  absl::flat_hash_map<State, std::vector<State>> successors;
  absl::flat_hash_map<std::pair<State, State>, int64_t> edge_costs;

  auto get_adjacent = [&successors,
                       &edge_costs](const State& start) -> std::vector<State> {
    auto known_it = successors.find(start);
    if (known_it != successors.end()) {
      return known_it->second;
    }

    std::vector<std::pair<State, int64_t>> state_and_cost =
        start.SuccessorStates();

    std::vector<State> states;
    for (const auto& [state, _] : state_and_cost) {
      states.emplace_back(state);
    }
    CHECK(successors.emplace(start, states).second);

    for (auto& [state, cost] : state_and_cost) {
      CHECK(edge_costs.emplace(std::make_pair(start, std::move(state)), cost)
                .second);
    }

    return states;
  };

  auto get_cost = [&edge_costs](const State& start,
                                const State& next) -> int64_t {
    auto cost_it = edge_costs.find({start, next});
    CHECK(cost_it != edge_costs.end());
    return cost_it->second;
  };

  auto path = aoc2021::AStar(initial_state, kGoalState, get_adjacent, get_cost);
  CHECK(path.has_value());
  std::cout << path->cost << "\n";

  return 0;
}
