#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_split.h"
#include "util/check.h"
#include "util/io.h"

namespace {

int CountPaths(
    const absl::flat_hash_map<std::string, absl::flat_hash_set<std::string>>&
        edges,
    absl::flat_hash_set<std::string>& visited_small, std::string current_cave) {
  if (current_cave == "end") return 1;

  int paths = 0;
  for (const std::string& next : edges.at(current_cave)) {
    if (!visited_small.contains(next)) {
      if (!absl::ascii_isupper(next.front())) visited_small.insert(next);
      paths += CountPaths(edges, visited_small, next);
      if (!absl::ascii_isupper(next.front())) visited_small.erase(next);
    }
  }
  return paths;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> edge_lines = aoc2021::ReadLinesFromFile(argv[1]);
  absl::flat_hash_map<std::string, absl::flat_hash_set<std::string>> edges;
  for (const std::string& edge_line : edge_lines) {
    std::vector<std::string> caves = absl::StrSplit(edge_line, '-');
    CHECK(caves.size() == 2);
    edges[caves.front()].insert(caves.back());
    edges[caves.back()].insert(caves.front());
  }

  absl::flat_hash_set<std::string> visited_small{"start"};
  std::cout << CountPaths(edges, visited_small, "start") << "\n";

  return 0;
}
