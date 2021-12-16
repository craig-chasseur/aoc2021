#ifndef UTIL_DIMENSION_GRID_H_
#define UTIL_DIMENSION_GRID_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "util/check.h"

namespace aoc2021 {

template <size_t dim>
class DimensionGrid {
 public:
  DimensionGrid() = delete;

  struct Vec;

  struct Point {
    std::array<int64_t, dim> coords = {};

    bool operator==(const Point& other) const { return coords == other.coords; }

    bool operator!=(const Point& other) const { return coords != other.coords; }

    Point& operator+=(const Vec& vec);
    Point& operator-=(const Vec& vec);

    std::vector<Point> AdjacentCardinal() const;
    std::vector<Point> AdjacentWithDiagonal() const;

    template <typename H>
    friend H AbslHashValue(H h, const Point& point) {
      return H::combine(std::move(h), point.coords);
    }
  };

  class Points {
   public:
    Points() = delete;

    inline static constexpr Point kOrigin{};
  };

  struct Vec {
    std::array<int64_t, dim> deltas = {};

    bool operator==(const Vec& other) const { return deltas == other.deltas; }

    bool operator!=(const Vec& other) const { return deltas != other.deltas; }

    Vec operator-() const {
      Vec negated(*this);
      for (int64_t& delta : negated.deltas) {
        delta = -delta;
      }
      return negated;
    }

    Vec& operator+=(const Vec& other) {
      for (size_t d = 0; d < dim; ++d) {
        deltas[d] += other.deltas[d];
      }
      return *this;
    }

    Vec operator+(const Vec& other) const {
      Vec result(*this);
      result += other;
      return result;
    }

    Vec& operator-=(const Vec& other) {
      *this += (-other);
      return *this;
    }

    Vec operator-(const Vec& other) const {
      Vec result(*this);
      result -= other;
      return result;
    }

    // Scaling.
    Vec& operator*=(const int64_t factor) {
      for (int64_t d : deltas) {
        d *= factor;
      }
      return *this;
    }

    Vec operator*(const int64_t factor) const {
      Vec result(*this);
      result *= factor;
      return result;
    }

    int64_t ManhattanDistance() const {
      return std::accumulate(
          deltas.begin(), deltas.end(), int64_t{0},
          [](int64_t sum, int64_t val) { return sum + std::abs(val); });
    }

    double Magnitude() const {
      double sum_of_squares = 0.0;
      for (const int64_t d : deltas) {
        sum_of_squares += d * d;
      }
      return std::sqrt(sum_of_squares);
    }

    bool Cardinal() const {
      return std::count_if(deltas.begin(), deltas.end(),
                           [](int64_t d) { return d != 0; }) == 1;
    }

    bool CardinalOrDiagonal() const {
      auto iter = std::find_if(deltas.begin(), deltas.end(),
                               [](int64_t d) { return d != 0; });
      if (iter == deltas.end()) return false;
      const int64_t mag = std::abs(*iter);
      return std::all_of(deltas.begin(), deltas.end(), [mag](int64_t d) {
        return d == 0 || std::abs(d) == mag;
      });
    }

    Vec UnitStep() const {
      CHECK(CardinalOrDiagonal());
      Vec step;
      for (size_t d = 0; d < dim; ++d) {
        if (deltas[d] < 0) {
          step.deltas[d] = -1;
        } else if (deltas[d] == 0) {
          step.deltas[d] = 0;
        } else {
          step.deltas[d] = 1;
        }
      }
      return step;
    }

    template <typename H>
    friend H AbslHashValue(H h, const Vec& vec) {
      return H::combine(std::move(h), vec.deltas);
    }
  };

  class Vecs {
   public:
    Vecs() = delete;

    static std::vector<Vec> Cardinal() {
      std::vector<Vec> cardinal;
      for (size_t d = 0; d < dim; ++d) {
        Vec minus, plus;
        minus.deltas[d] = -1;
        plus.deltas[d] = 1;
        cardinal.emplace_back(std::move(minus));
        cardinal.emplace_back(std::move(plus));
      }
      return cardinal;
    }

    static std::vector<Vec> CardinalAndDiagonal() {
      std::vector<Vec> result(1);
      for (size_t d = 0; d < dim; ++d) {
        std::vector<Vec> expanded_result;
        for (const Vec& v : result) {
          Vec minus(v);
          minus.deltas[d] = -1;
          expanded_result.emplace_back(std::move(minus));

          Vec plus(v);
          plus.deltas[d] = 1;
          expanded_result.emplace_back(std::move(plus));

          if (d == dim - 1 && v.ManhattanDistance() == 0) continue;
          expanded_result.emplace_back(v);
        }
        result = std::move(expanded_result);
      }
      return result;
    }
  };

  // Represents an infinite hyperplane with dimensionality dim - 1. Currently
  // only handles hyperplanes that are fixed on one axis.
  struct Hyperplane {
    size_t fixed_coord = 0;
    int64_t fixed_coord_value;

    Point Reflect(const Point& p) const;
  };

  template <typename PointContainer>
  static Point MinDimensions(const PointContainer& container) {
    Point min;
    std::fill(min.coords.begin(), min.coords.end(),
              std::numeric_limits<int64_t>::max());
    for (const Point& point : container) {
      for (size_t d = 0; d < dim; ++d) {
        min.coords[d] = std::min(min.coords[d], point.coords[d]);
      }
    }
    return min;
  }

  template <typename PointContainer>
  static Point MaxDimensions(const PointContainer& container) {
    Point max;
    std::fill(max.coords.begin(), max.coords.end(),
              std::numeric_limits<int64_t>::min());
    for (const Point& point : container) {
      for (size_t d = 0; d < dim; ++d) {
        max.coords[d] = std::max(max.coords[d], point.coords[d]);
      }
    }
    return max;
  }
};

// Printing.
template <size_t dim>
std::ostream& operator<<(std::ostream& os,
                         const typename DimensionGrid<dim>::Point& p) {
  os << "(";
  for (size_t d = 0; d < dim; ++d) {
    os << p.coords[d];
    if (d != dim - 1) os << ", ";
  }
  os << ")";
  return os;
}

template <size_t dim>
std::ostream& operator<<(std::ostream& os,
                         const typename DimensionGrid<dim>::Vec& v) {
  os << "Vec{";
  for (size_t d = 0; d < dim; ++d) {
    os << v.deltas[d];
    if (d != dim - 1) os << ", ";
  }
  os << "}";
  return os;
}

template <size_t dim>
std::ostream& operator<<(std::ostream& os,
                         const typename DimensionGrid<dim>::Hyperplane& h) {
  os << "Hyperplane(dim_" << h.fixed_coord << " = " << h.fixed_coord_value
     << ")";
  return os;
}

// Arithmetic between Points and Vecs.
template <size_t dim>
typename DimensionGrid<dim>::Point& DimensionGrid<dim>::Point::operator+=(
    const Vec& vec) {
  for (size_t d = 0; d < dim; ++d) {
    this->coords[d] += vec.deltas[d];
  }
  return *this;
}

template <size_t dim>
typename DimensionGrid<dim>::Point& DimensionGrid<dim>::Point::operator-=(
    const Vec& vec) {
  *this += (-vec);
  return *this;
}

template <size_t dim>
typename DimensionGrid<dim>::Vec operator-(
    const typename DimensionGrid<dim>::Point& a,
    const typename DimensionGrid<dim>::Point& b) {
  typename DimensionGrid<dim>::Vec diff;
  for (size_t d = 0; d < dim; ++d) {
    diff.deltas[d] = a.coords[d] - b.coords[d];
  }
  return diff;
}

template <size_t dim>
typename DimensionGrid<dim>::Point operator+(
    const typename DimensionGrid<dim>::Point& p,
    const typename DimensionGrid<dim>::Vec& v) {
  typename DimensionGrid<dim>::Point result;
  for (size_t d = 0; d < dim; ++d) {
    result.coords[d] = p.coords[d] + v.deltas[d];
  }
}

template <size_t dim>
typename DimensionGrid<dim>::Point operator-(
    const typename DimensionGrid<dim>::Point& p,
    const typename DimensionGrid<dim>::Vec& v) {
  return p + (-v);
}

template <size_t dim, typename VecContainer>
std::enable_if_t<std::is_same_v<typename DimensionGrid<dim>::Vec,
                                typename VecContainer::value_type>,
                 std::vector<typename DimensionGrid<dim>::Point>>
operator+(const typename DimensionGrid<dim>::Point& p,
          const VecContainer& vecs) {
  std::vector<typename DimensionGrid<dim>::Point> result;
  for (const typename DimensionGrid<dim>::Vec& v : vecs) {
    result.emplace_back(p + v);
  }
  return result;
}

// Arithmetic between Points and Hyperplanes.

// Returns the shortest vector from `h` to `p`, which is also guaranteed to be
// normal to `h`.
template <size_t dim>
typename DimensionGrid<dim>::Vec operator-(
    const typename DimensionGrid<dim>::Point& p,
    const typename DimensionGrid<dim>::Hyperplane& h) {
  typename DimensionGrid<dim>::Vec diff;
  diff.deltas[h.fixed_coord] = p.coords[h.fixed_coord] - h.fixed_coord_value;
  return diff;
}

// Out-of-line override for scaling a Vec with the scaling factor first.
template <size_t dim>
typename DimensionGrid<dim>::Vec operator*(
    const int64_t factor, const typename DimensionGrid<dim>::Vec& vec) {
  return vec * factor;
}

template <size_t dim>
std::vector<typename DimensionGrid<dim>::Point>
DimensionGrid<dim>::Point::AdjacentCardinal() const {
  return *this + DimensionGrid<dim>::Vecs::Cardinal();
}

template <size_t dim>
std::vector<typename DimensionGrid<dim>::Point>
DimensionGrid<dim>::Point::AdjacentWithDiagonal() const {
  return *this + DimensionGrid<dim>::Vecs::CardinalAndDiagonal();
}

template <size_t dim>
typename DimensionGrid<dim>::Point DimensionGrid<dim>::Hyperplane::Reflect(
    const typename DimensionGrid<dim>::Point& p) const {
  return p - 2 * (p - *this);
}

};  // namespace aoc2021

#endif  // UTIL_DIMENSION_GRID_H_
