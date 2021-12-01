#include "util/io.h"

#include <fstream>
#include <ios>
#include <string>
#include <utility>
#include <vector>

#include "absl/strings/str_split.h"
#include "util/check.h"

namespace aoc2021 {

std::string ReadFile(const char* filename) {
  std::ifstream stream(filename);
  CHECK(stream);

  std::string buffer;
  stream.seekg(0, std::ios::end);
  buffer.resize(stream.tellg());
  stream.seekg(0, std::ios::beg);
  stream.read(&buffer[0], buffer.size());
  stream.close();

  return buffer;
}

std::vector<std::string> ReadLinesFromFile(const char* filename) {
  std::ifstream stream(filename);
  CHECK(stream);

  std::vector<std::string> all_lines;
  std::string line;
  while (std::getline(stream, line)) {
    all_lines.emplace_back(std::move(line));
  }
  stream.close();
  return all_lines;
}

std::vector<std::string> ReadCommaDelimitedFile(const char* filename) {
  std::string contents = ReadFile(filename);
  return absl::StrSplit(contents, ',');
}

std::vector<std::vector<std::string>> SplitByEmptyStrings(
    std::vector<std::string> strs) {
  std::vector<std::vector<std::string>> splits;
  std::vector<std::string> current_split;
  for (std::string& str : strs) {
    if (!str.empty()) {
      current_split.emplace_back(std::move(str));
      continue;
    }

    if (current_split.empty()) continue;

    splits.emplace_back(std::move(current_split));
    current_split.clear();
  }
  if (!current_split.empty()) splits.emplace_back(std::move(current_split));
  return splits;
}

}  // namespace aoc2021
