#include <cstdint>
#include <deque>
#include <numeric>
#include <string>

#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

static constexpr bool kPart2 = true;
static constexpr int kDays = kPart2 ? 256 : 80;

}  // namespace

int main(int argc, char** argv) {
  std::string input = aoc2021::ReadFile(argv[1]);
  std::deque<int64_t> fish(9, 0);
  for (absl::string_view str : absl::StrSplit(input, ',')) {
    int timer = 0;
    CHECK(absl::SimpleAtoi(str, &timer));
    CHECK(timer >= 0);
    CHECK(timer < 9);
    ++fish[timer];
  }

  for (int day = 0; day < kDays; ++day) {
    const int64_t zeroes = fish.front();
    fish.pop_front();
    fish[6] += zeroes;
    fish.push_back(zeroes);
  }

  std::cout << std::accumulate(fish.begin(), fish.end(), int64_t{0}) << "\n";

  return 0;
}
