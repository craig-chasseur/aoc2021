#include "util/grid2.h"

#include <vector>

#include "util/check.h"

namespace aoc2021::grid2 {

Point Point::Rotate90(const Point& pivot, int turns) const {
  return pivot + (*this - pivot).Rotate90(turns);
}

std::vector<Point> Point::AdjacentCardinal() const {
  return *this + Vecs::kCardinal;
}

std::vector<Point> Point::AdjacentWithDiagonal() const {
  return *this + Vecs::kAdjacent8;
}

Vec Vec::UnitStep() const {
  CHECK(Horizontal() || Vertical() || Diag45());
  Vec step;

  if (dx < 0) {
    step.dx = -1;
  } else if (dx == 0) {
    step.dx = 0;
  } else {
    step.dx = 1;
  }

  if (dy < 0) {
    step.dy = -1;
  } else if (dy == 0) {
    step.dy = 0;
  } else {
    step.dy = 1;
  }

  return step;
}

}  // namespace aoc2021::grid2
