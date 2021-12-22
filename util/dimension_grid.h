#ifndef UTIL_DIMENSION_GRID_H_
#define UTIL_DIMENSION_GRID_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
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
  class MultiVolume;

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

    friend Vec operator-(const Point& a, const Point& b) {
      Vec diff;
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
      return *this + Vecs::Cardinal();
    }

    std::vector<Point> AdjacentWithDiagonal() const {
      return *this + Vecs::CardinalAndDiagonal();
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

    static const std::array<Vec, dim * 2>& Cardinal() {
      static const std::array<Vec, dim* 2>* const cardinal = ComputeCardinal();
      return *cardinal;
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

   private:
    static const std::array<Vec, 2 * dim>* ComputeCardinal() {
      auto cardinal = new std::array<Vec, 2 * dim>();
      for (size_t d = 0; d < dim; ++d) {
        (*cardinal)[d * 2].deltas[d] = -1;
        (*cardinal)[d * 2 + 1].deltas[d] = 1;
      }
      return cardinal;
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
      // Only 3D supported right now.
      static_assert(dim == 3);
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
      static_assert(dim == 3);

      static const std::array<Rotation, 24>* const orientations =
          ComputeAllOrientations();
      return *orientations;
    }

   private:
    static const std::array<Rotation, 24>* ComputeAllOrientations() {
      // Only 3D is currently supported.
      static_assert(dim == 3);
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

  struct Orthotope {
    Point min_point;
    Point max_point;

    class iterator {
     public:
      using iterator_category = std::input_iterator_tag;
      using difference_type = std::ptrdiff_t;
      using value_type = Point;
      using pointer = const Point*;
      using reference = const Point&;

      iterator() = default;

      reference operator*() const { return current_; }
      pointer operator->() const { return &current_; }

      iterator& operator++() {
        for (size_t d = 0; d < dim; ++d) {
          if (++current_.coords[d] <= rect_->max_point.coords[d]) {
            break;
          }
          if (d == dim - 1) break;
          current_.coords[d] = rect_->min_point.coords[d];
        }
        return *this;
      }

      iterator operator++(int) const {
        iterator tmp(*this);
        ++tmp;
        return tmp;
      }

      bool operator==(const iterator& other) const {
        return current_ == other.current_ && rect_ == other.rect_;
      }

     private:
      friend class Orthotope;

      explicit iterator(Point current, const Orthotope* rect)
          : current_(current), rect_(rect) {}

      Point current_;
      const Orthotope* rect_ = nullptr;
    };
    using const_iterator = iterator;

    bool operator==(const Orthotope& other) const {
      return min_point == other.min_point && max_point == other.max_point;
    }

    Vec Diagonal() const { return max_point - min_point; }

    int64_t HyperVolume() const {
      const Vec diag = Diagonal();
      return std::accumulate(diag.deltas.begin(), diag.deltas.end(), int64_t{1},
                             std::multiplies<int64_t>());
    }

    int64_t HyperVolumeInclusive() const {
      Vec diag = Diagonal();
      for (size_t d = 0; d < dim; ++d) {
        ++diag.deltas[d];
      }
      return std::accumulate(diag.deltas.begin(), diag.deltas.end(), int64_t{1},
                             std::multiplies<int64_t>());
    }

    bool Contains(const Point& p) const {
      for (size_t d = 0; d < dim; ++d) {
        if (p.coords[d] < min_point.coords[d]) return false;
        if (p.coords[d] > max_point.coords[d]) return false;
      }
      return true;
    }

    template <typename PointContainer>
    std::vector<Point> FilterContains(const PointContainer& points) const {
      std::vector<Point> filtered;
      for (const Point& point : points) {
        if (Contains(point)) filtered.emplace_back(point);
      }
      return filtered;
    }

    std::optional<Orthotope> Intersection(const Orthotope& other) const {
      Orthotope intersection;
      for (size_t d = 0; d < dim; ++d) {
        intersection.min_point.coords[d] =
            std::max(min_point.coords[d], other.min_point.coords[d]);
        intersection.max_point.coords[d] =
            std::min(max_point.coords[d], other.max_point.coords[d]);
        if (intersection.min_point.coords[d] > intersection.max_point.coords[d])
          return std::nullopt;
      }
      return intersection;
    }

    MultiVolume SpatialDifference(const Orthotope& sub) const {
      std::optional<Orthotope> intersection = Intersection(sub);
      if (!intersection.has_value()) return MultiVolume(*this);
      if (*intersection == *this) return MultiVolume();

      Orthotope remainder = *this;
      MultiVolume slices;
      for (size_t d = 0; d < dim; ++d) {
        if (remainder.min_point.coords[d] < intersection->min_point.coords[d]) {
          Orthotope slice(remainder);
          slice.max_point.coords[d] = intersection->min_point.coords[d] - 1;
          slices.regions_.emplace_back(std::move(slice));
          remainder.min_point.coords[d] = intersection->min_point.coords[d];
        }
        if (remainder.HyperVolumeInclusive() <= 0) break;
        if (remainder.max_point.coords[d] > intersection->max_point.coords[d]) {
          Orthotope slice(remainder);
          slice.min_point.coords[d] = intersection->max_point.coords[d] + 1;
          slices.regions_.emplace_back(std::move(slice));
          remainder.max_point.coords[d] = intersection->max_point.coords[d];
        }
        if (remainder.HyperVolumeInclusive() <= 0) break;
      }

      if (remainder.HyperVolumeInclusive() > 0) {
        CHECK(remainder == *intersection);
      }
      return slices;
    }

    MultiVolume Union(const Orthotope& other) const {
      MultiVolume result;

      std::optional<Orthotope> intersection = Intersection(other);
      if (!intersection.has_value()) {
        result.regions_ = {*this, other};
        return result;
      }

      const bool other_larger =
          other.HyperVolumeInclusive() > HyperVolumeInclusive();
      if (other_larger) {
        result.regions_.emplace_back(other);
      } else {
        result.regions_.emplace_back(*this);
      }
      const Orthotope& smaller = other_larger ? *this : other;
      MultiVolume diff = smaller.SpatialDifference(*intersection);
      result.regions_.insert(result.regions_.end(),
                             std::make_move_iterator(diff.regions_.begin()),
                             std::make_move_iterator(diff.regions_.end()));
      return result;
    }

    const_iterator cbegin() const { return const_iterator(min_point, this); }
    const_iterator begin() const { return cbegin(); }

    const_iterator cend() const {
      Point beyond;
      for (size_t d = 0; d < dim; ++d) {
        if (d == dim - 1) {
          beyond.coords[d] = max_point.coords[d] + 1;
        } else {
          beyond.coords[d] = min_point.coords[d];
        }
      }
      return const_iterator(beyond, this);
    }
    const_iterator end() const { return cend(); }
  };

  class MultiVolume {
   public:
    MultiVolume() = default;
    explicit MultiVolume(Orthotope ortho) : regions_{std::move(ortho)} {}

    bool empty() const { return regions_.empty(); }

    MultiVolume& operator+=(const Orthotope& ortho) {
      if (regions_.empty()) {
        regions_.emplace_back(ortho);
        return *this;
      }
      std::vector<Orthotope> remainder{ortho};
      for (const Orthotope& region : regions_) {
        std::vector<Orthotope> next_remainder;
        for (const Orthotope& remainder_region : remainder) {
          MultiVolume diff = remainder_region.SpatialDifference(region);
          next_remainder.insert(next_remainder.end(),
                                std::make_move_iterator(diff.regions_.begin()),
                                std::make_move_iterator(diff.regions_.end()));
        }
        remainder = std::move(next_remainder);
        if (remainder.empty()) break;
      }
      regions_.insert(regions_.end(),
                      std::make_move_iterator(remainder.begin()),
                      std::make_move_iterator(remainder.end()));
      return *this;
    }

    MultiVolume operator+(const Orthotope& ortho) const {
      MultiVolume result(*this);
      result += ortho;
      return result;
    }

    MultiVolume operator-(const Orthotope& ortho) const {
      MultiVolume result;
      for (const Orthotope& region : regions_) {
        MultiVolume diff = region.SpatialDifference(ortho);
        result.regions_.insert(result.regions_.end(),
                               std::make_move_iterator(diff.regions_.begin()),
                               std::make_move_iterator(diff.regions_.end()));
      }
      return result;
    }

    MultiVolume& operator-=(const Orthotope& ortho) {
      *this = *this - ortho;
      return *this;
    }

    int64_t HyperVolume() const {
      return std::accumulate(regions_.begin(), regions_.end(), int64_t{0},
                             [](int64_t sum, const Orthotope& region) {
                               return sum + region.HyperVolume();
                             });
    }

    int64_t HyperVolumeInclusive() const {
      return std::accumulate(regions_.begin(), regions_.end(), int64_t{0},
                             [](int64_t sum, const Orthotope& region) {
                               return sum + region.HyperVolumeInclusive();
                             });
    }

    bool Contains(const Point& p) const {
      for (const Orthotope& region : regions_) {
        if (region.Contains(p)) return true;
      }
      return false;
    }

    template <typename PointContainer>
    std::vector<Point> FilterContains(const PointContainer& points) const {
      std::vector<Point> filtered;
      for (const Point& point : points) {
        if (Contains(point)) filtered.emplace_back(point);
      }
      return filtered;
    }

   private:
    friend class Orthotope;

    std::vector<Orthotope> regions_;
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

  template <typename PointContainer>
  static Orthotope BoundingBox(const PointContainer& container) {
    Orthotope box;
    if (container.empty()) return box;
    box.min_point = MinDimensions(container);
    box.max_point = MaxDimensions(container);
    return box;
  }
};

};  // namespace aoc2021

#endif  // UTIL_DIMENSION_GRID_H_
