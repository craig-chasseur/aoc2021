#include <cstdint>
#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/grid2.h"
#include "util/io.h"

using namespace aoc2021::grid2;

namespace {

std::vector<bool> ParseEnhancementProgram(const absl::string_view rep) {
  std::vector<bool> program;
  for (char c : rep) {
    switch (c) {
      case '#':
        program.push_back(true);
        break;
      case '.':
        program.push_back(false);
        break;
      default:
        CHECK_FAIL();
    }
  }
  CHECK(program.size() == 512);
  return program;
}

class Image {
 public:
  Image() = default;

  explicit Image(const std::vector<std::string>& rep) {
    for (int64_t y = 0; y < rep.size(); ++y) {
      for (int64_t x = 0; x < rep[y].size(); ++x) {
        switch (rep[y][x]) {
          case '#':
            bright_pixels_.emplace(Point{x, y});
            break;
          case '.':
            break;
          default:
            CHECK_FAIL();
        }
      }
    }
  }

  Image ApplyProgram(const std::vector<bool>& program) const {
    Image result;
    if (!inverted_ && program.front()) {
      result.inverted_ = true;
    } else if (inverted_ && !program.back()) {
      result.inverted_ = false;
    }

    Point min = MinDimensions(bright_pixels_) - Vec{.dx = 1, .dy = 1};
    Point max = MaxDimensions(bright_pixels_) + Vec{.dx = 1, .dy = 1};

    for (int64_t y = min.y; y <= max.y; ++y) {
      for (int64_t x = min.x; x <= max.x; ++x) {
        const bool bright = program[ValueForPixel(Point{x, y})];
        if (bright ^ result.inverted_) {
          result.bright_pixels_.emplace(Point{x, y});
        }
      }
    }
    return result;
  }

  int64_t PixelsLit() const {
    CHECK(!inverted_);
    return bright_pixels_.size();
  }

 private:
  int ValueForPixel(Point p) const {
    int value = 0;
    for (int64_t dy : {-1, 0, 1}) {
      for (int64_t dx : {-1, 0, 1}) {
        value =
            (value << 1) | bright_pixels_.contains(p + Vec{.dx = dx, .dy = dy});
      }
    }
    if (inverted_) {
      value = ~value & 511;
    }
    CHECK(value >= 0);
    CHECK(value < 512);
    return value;
  }

  absl::flat_hash_set<Point> bright_pixels_;
  bool inverted_ = false;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<std::vector<std::string>> segments =
      aoc2021::SplitByEmptyStrings(std::move(lines));
  CHECK(segments.size() == 2);

  CHECK(segments.front().size() == 1);
  std::vector<bool> program = ParseEnhancementProgram(segments.front().front());

  Image image(segments.back());
  for (int iteration = 0; iteration < 50; ++iteration) {
    image = image.ApplyProgram(program);
    if (iteration == 1) {
      std::cout << "Part 1: " << image.PixelsLit() << "\n";
    }
  }
  std::cout << "Part 2: " << image.PixelsLit() << "\n";

  return 0;
}
