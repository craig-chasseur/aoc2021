#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "re2/re2.h"
#include "util/check.h"
#include "util/dimension_grid.h"
#include "util/io.h"

using Grid3 = aoc2021::DimensionGrid<3>;

namespace {

struct Procedure {
  Grid3::Orthotope box;
  bool on = false;

  std::optional<Procedure> Clamp() const {
    Procedure clamped;
    clamped.on = on;
    for (size_t d = 0; d < 3; ++d) {
      if (box.min_point.coords[d] > 50) return std::nullopt;
      if (box.max_point.coords[d] < -50) return std::nullopt;
      clamped.box.min_point.coords[d] =
          std::max(box.min_point.coords[d], int64_t{-50});
      clamped.box.max_point.coords[d] =
          std::min(box.max_point.coords[d], int64_t{50});
    }
    return clamped;
  }
};

Procedure ParseProcedure(absl::string_view line) {
  static re2::LazyRE2 pattern = {
      R"re((\w+) x=(-?\d+)\.\.(-?\d+),y=(-?\d+)\.\.(-?\d+),z=(-?\d+)\.\.(-?\d+))re"};
  Procedure parsed;
  std::string on_or_off;
  CHECK(re2::RE2::FullMatch(
      line, *pattern, &on_or_off, &parsed.box.min_point.coords[0],
      &parsed.box.max_point.coords[0], &parsed.box.min_point.coords[1],
      &parsed.box.max_point.coords[1], &parsed.box.min_point.coords[2],
      &parsed.box.max_point.coords[2]));
  if (on_or_off == "on") {
    parsed.on = true;
  } else if (on_or_off == "off") {
    parsed.on = false;
  } else {
    CHECK_FAIL();
  }
  return parsed;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<Procedure> procedures;
  for (absl::string_view line : lines) {
    procedures.emplace_back(ParseProcedure(line));
  }

  Grid3::MultiVolume on;
  for (const Procedure& procedure : procedures) {
    if (procedure.on) {
      on += procedure.box;
    } else {
      on -= procedure.box;
    }
  }
  std::cout << on.HyperVolumeInclusive() << "\n";

  return 0;
}
