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

#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "util/check.h"

namespace aoc2021 {

template <size_t dim>
class DimensionGrid {
 public:
  DimensionGrid() = delete;

  struct Vec;
  struct Rotation;

  struct Point {
    std::array<int64_t, dim> coords = {};

    static Point ParseCommaSeparated(const absl::string_view rep) {
      std::vector<absl::string_view> coord_strings = absl::StrSplit(rep, ',');
      CHECK(coord_strings.size() == dim);
      Point parsed;
      for (size_t d = 0; d < dim; ++d) {
        CHECK(absl::SimpleAtoi(coord_strings[d], &parsed.coords[d]));
      }
      return parsed;
    }

    bool operator==(const Point& other) const { return coords == other.coords; }

    bool operator!=(const Point& other) const { return coords != other.coords; }

    Point& operator+=(const Vec& vec) {
      for (size_t d = 0; d < dim; ++d) {
        this->coords[d] += vec.deltas[d];
      }
      return *this;
    }

    Point operator+(const Vec& vec) const {
      Point result(*this);
      result += vec;
      return result;
    }

    Point& operator-=(const Vec& vec) {
      *this += (-vec);
      return *this;
    }

    Point operator-(const Vec& vec) const {
      Point result(*this);
      result -= vec;
      return result;
    }

    friend Vec operator-(const typename DimensionGrid<dim>::Point& a,
                         const typename DimensionGrid<dim>::Point& b) {
      typename DimensionGrid<dim>::Vec diff;
      for (size_t d = 0; d < dim; ++d) {
        diff.deltas[d] = a.coords[d] - b.coords[d];
      }
      return diff;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& p) {
      os << "(";
      for (size_t d = 0; d < dim; ++d) {
        os << p.coords[d];
        if (d != dim - 1) os << ", ";
      }
      os << ")";
      return os;
    }

    std::vector<Point> AdjacentCardinal() const {
      return *this + DimensionGrid<dim>::Vecs::Cardinal();
    }

    std::vector<Point> AdjacentWithDiagonal() const {
      return *this + DimensionGrid<dim>::Vecs::CardinalAndDiagonal();
    }

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

    Vec& operator*=(const Rotation& rot) {
      *this = *this * rot;
      return *this;
    }

    Vec operator*(const Rotation& rot) const {
      Vec result;
      for (size_t d = 0; d < dim; ++d) {
        for (size_t k = 0; k < dim; ++k) {
          result.deltas[d] += rot.matrix[d][k] * deltas[k];
        }
      }
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

    friend std::ostream& operator<<(std::ostream& os, const Vec& v) {
      os << "Vec{";
      for (size_t d = 0; d < dim; ++d) {
        os << v.deltas[d];
        if (d != dim - 1) os << ", ";
      }
      os << "}";
      return os;
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

  // Scaling Vec with scale factor first.
  friend Vec operator*(const int64_t factor, Vec& vec) { return vec * factor; }

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

    friend std::ostream& operator<<(std::ostream& os, const Hyperplane& h) {
      os << "Hyperplane(dim_" << h.fixed_coord << " = " << h.fixed_coord_value
         << ")";
      return os;
    }

    Point Reflect(const Point& p) const { return p - 2 * (p - *this); }
  };

  struct Rotation {
    std::array<std::array<int8_t, dim>, dim> matrix = {};

    static Rotation AboutAxis90(size_t axis, int turns) {
      CHECK(dim == 3);
      CHECK(axis < dim);

      int8_t cos = 0;
      int8_t sin = 0;
      switch (turns % 4) {
        case 0:
          cos = 1;
          sin = 0;
          break;
        case 1:
          cos = 0;
          sin = 1;
          break;
        case 2:
          cos = -1;
          sin = 0;
          break;
        case 3:
          cos = 0;
          sin = -1;
      }

      Rotation rot;
      rot.matrix[axis][axis] = 1;
      switch (axis) {
        case 0:
          rot.matrix[1][1] = cos;
          rot.matrix[1][2] = -sin;
          rot.matrix[2][1] = sin;
          rot.matrix[2][2] = cos;
          break;
        case 1:
          rot.matrix[0][0] = cos;
          rot.matrix[0][2] = sin;
          rot.matrix[2][0] = -sin;
          rot.matrix[2][2] = cos;
          break;
        case 2:
          rot.matrix[0][0] = cos;
          rot.matrix[0][1] = -sin;
          rot.matrix[1][0] = sin;
          rot.matrix[1][1] = cos;
          break;
      }
      return rot;
    }

    Rotation operator*(const Rotation& other) const {
      Rotation product;
      for (size_t row = 0; row < dim; ++row) {
        for (size_t col = 0; col < dim; ++col) {
          for (size_t k = 0; k < dim; ++k) {
            product.matrix[row][col] += matrix[row][k] * other.matrix[k][col];
          }
        }
      }
      return product;
    }

    Rotation& operator*=(const Rotation& other) {
      *this = *this * other;
      return *this;
    }
  };

  class Rotations {
   public:
    Rotations() = delete;

    static const std::array<Rotation, 24>& AllOrientations() {
      // Only 3D is currently supported.
      CHECK(dim == 3);

      static const std::array<Rotation, 24>* const orientations =
          ComputeAllOrientations();
      return *orientations;
    }

   private:
    static const std::array<Rotation, 24>* ComputeAllOrientations() {
      // Only 3D is currently supported.
      CHECK(dim == 3);
      auto orientations = new std::array<Rotation, 24>();
      for (int z_turns = 0; z_turns < 4; ++z_turns) {
        Rotation z_rot = Rotation::AboutAxis90(2, z_turns);
        for (int x_turns = 0; x_turns < 4; ++x_turns) {
          (*orientations)[z_turns * 4 + x_turns] =
              z_rot * Rotation::AboutAxis90(0, x_turns);
        }
      }
      for (int y_turns : {1, 3}) {
        Rotation y_rot = Rotation::AboutAxis90(1, y_turns);
        for (int x_turns = 0; x_turns < 4; ++x_turns) {
          (*orientations)[16 + (y_turns == 3 ? 4 : 0) + x_turns] =
              y_rot * Rotation::AboutAxis90(0, x_turns);
        }
      }
      return orientations;
    }
  };

  // Arithmetic between Points and Vecs.
  template <typename VecContainer>
  std::enable_if_t<std::is_same_v<Vec, typename VecContainer::value_type>,
                   std::vector<Point>> friend
  operator+(const Point& p, const VecContainer& vecs) {
    std::vector<Point> result;
    for (const Vec& v : vecs) {
      result.emplace_back(p + v);
    }
    return result;
  }

  // Arithmetic between Vecs and Rotations.
  template <typename VecContainer>
  friend std::enable_if_t<
      std::is_same_v<Vec, typename VecContainer::value_type>, std::vector<Vec>>
  operator*(const VecContainer& vecs, const Rotation& rot) {
    std::vector<Vec> result;
    for (const Vec& v : vecs) {
      result.emplace_back(v * rot);
    }
    return result;
  }

  // Arithmetic between Points and Hyperplanes.

  // Returns the shortest vector from `h` to `p`, which is also guaranteed to be
  // normal to `h`.
  friend Vec operator-(const Point& p, const Hyperplane& h) {
    Vec diff;
    diff.deltas[h.fixed_coord] = p.coords[h.fixed_coord] - h.fixed_coord_value;
    return diff;
  }

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

};  // namespace aoc2021

#endif  // UTIL_DIMENSION_GRID_H_
