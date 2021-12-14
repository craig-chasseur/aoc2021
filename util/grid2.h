#ifndef UTIL_GRID2_H_
#define UTIL_GRID2_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "util/check.h"

namespace aoc2021::grid2 {

struct Vec;

// 2-dimensional point.
struct Point {
  int64_t x = 0;
  int64_t y = 0;

  bool operator==(const Point& other) const {
    return x == other.x && y == other.y;
  }

  bool operator!=(const Point& other) const {
    return x != other.x || y != other.y;
  }

  Point& operator+=(const Vec& vec);
  Point& operator-=(const Vec& vec);

  std::vector<Point> AdjacentCardinal() const;
  std::vector<Point> AdjacentWithDiagonal() const;

  template <typename H>
  friend H AbslHashValue(H h, const Point& point) {
    return H::combine(std::move(h), point.x, point.y);
  }
};

// Constants for 2-dimensional points.
class Points {
 public:
  Points() = delete;

  inline static constexpr Point kOrigin{.x = 0, .y = 0};
};

// Two-dimensional vector. Can represent the difference between Points, and can
// be added to a Point to yield another point.
struct Vec {
  int64_t dx = 0;
  int64_t dy = 0;

  bool operator==(const Vec& other) const {
    return dx == other.dx && dy == other.dy;
  }

  bool operator!=(const Vec& other) const {
    return dx != other.dx || dy != other.dy;
  }

  Vec operator-() const { return Vec{.dx = -dx, .dy = -dy}; }

  Vec& operator+=(const Vec& other) {
    dx += other.dx;
    dy += other.dy;
    return *this;
  }

  Vec operator+(const Vec& other) const {
    Vec result(*this);
    result += other;
    return result;
  }

  Vec& operator-=(const Vec& other) {
    dx -= other.dx;
    dy -= other.dy;
    return *this;
  }

  Vec operator-(const Vec& other) const {
    Vec result(*this);
    result -= other;
    return result;
  }

  // Scaling
  Vec& operator*=(const int64_t factor) {
    dx *= factor;
    dy *= factor;
    return *this;
  }

  Vec operator*(const int64_t factor) const {
    Vec result(*this);
    result *= factor;
    return result;
  }

  int64_t ManhattanDistance() const { return std::abs(dx) + std::abs(dy); }

  double Magnitude() const { return std::sqrt(dx * dx + dy * dy); }

  bool Horizontal() const { return dy == 0; }

  bool Vertical() const { return dx == 0; }

  bool Diag45() const { return std::abs(dx) == std::abs(dy); }

  // Iff Vec is horizontal, vertical, or a 45-degree diagonal, this gives a
  // unit-step in the direction of the vector.
  Vec UnitStep() const;

  template <typename H>
  friend H AbslHashValue(H h, const Vec& vec) {
    return H::combine(std::move(h), vec.dx, vec.dy);
  }
};

// Constants for 2-dimensional vectors.
class Vecs {
 public:
  Vecs() = delete;

  inline static constexpr Vec kLeft1{.dx = -1, .dy = 0};
  inline static constexpr Vec kRight1{.dx = 1, .dy = 0};
  inline static constexpr Vec kUp1{.dx = 0, .dy = -1};
  inline static constexpr Vec kDown1{.dx = 0, .dy = 1};

  inline static constexpr std::array<Vec, 4> kCardinal = {kLeft1, kRight1, kUp1,
                                                          kDown1};
  inline static constexpr std::array<Vec, 4> kDiagonal = {
      Vec{.dx = -1, .dy = -1}, Vec{.dx = -1, .dy = 1}, Vec{.dx = 1, .dy = -1},
      Vec{.dx = 1, .dy = 1}};
  inline static constexpr std::array<Vec, 8> kAdjacent8 = {
      kLeft1,
      kRight1,
      kUp1,
      kDown1,
      Vec{.dx = -1, .dy = -1},
      Vec{.dx = -1, .dy = 1},
      Vec{.dx = 1, .dy = -1},
      Vec{.dx = 1, .dy = 1}};
};

// Represents an infinite line. Currently only handles lines that are horizontal
// or vertical.
struct Line {
  enum class FixedCoord : bool {
    kX,
    kY
  };

  FixedCoord fixed_coord = FixedCoord::kX;
  int64_t fixed_coord_value = 0;

  static Line FixedX(const int64_t x) {
    return Line{.fixed_coord = FixedCoord::kX, .fixed_coord_value = x};
  }

  static Line FixedY(const int64_t y) {
    return Line{.fixed_coord = FixedCoord::kY, .fixed_coord_value = y};
  }

  Point Reflect(const Point& p) const;
};

// Arithmetic between Points and Vecs.

inline Point& Point::operator+=(const Vec& vec) {
  x += vec.dx;
  y += vec.dy;
  return *this;
}

inline Point& Point::operator-=(const Vec& vec) {
  *this += -vec;
  return *this;
}

inline Vec operator-(const Point& a, const Point& b) {
  return Vec{.dx = a.x - b.x, .dy = a.y - b.y};
}

inline Point operator+(const Point& p, const Vec& v) {
  return Point{.x = p.x + v.dx, .y = p.y + v.dy};
}

inline Point operator-(const Point& p, const Vec& v) { return p + (-v); }

template <typename VecContainer>
std::enable_if_t<std::is_same_v<Vec, typename VecContainer::value_type>,
                 std::vector<Point>>
operator+(const Point& p, const VecContainer& vecs) {
  std::vector<Point> result;
  for (const Vec& v : vecs) {
    result.emplace_back(p + v);
  }
  return result;
}

// Arithmetic between Points and Lines.

// Returns the shortest vector from `l` to `p`, which is also guaranteed to be
// normal to `l`.
inline Vec operator-(const Point& p, const Line& l) {
  switch (l.fixed_coord) {
    case Line::FixedCoord::kX:
      return Vec{.dx = p.x - l.fixed_coord_value, .dy = 0};
    case Line::FixedCoord::kY:
      return Vec{.dx = 0, .dy = p.y - l.fixed_coord_value};
  }
}

// Out-of-line override for scaling a Vec with the scaling factor first.
inline Vec operator*(const int64_t factor, const Vec& vec) {
  return vec * factor;
}

inline Point Line::Reflect(const Point &p) const {
  return p - 2 * (p - *this);
}

// A dense two-dimensional grid of 'T' values. Cells are addressable as points.
template <typename T = int>
class Grid {
 public:
  // Iterator implementation.
  template <typename ValueType, typename OuterIt, typename InnerIt>
  class IteratorImpl {
   public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = ValueType;
    using pointer = ValueType*;
    using reference = ValueType&;

    IteratorImpl() = default;

    reference operator*() const { return *inner_; }
    pointer operator->() const { return &*inner_; }

    IteratorImpl& operator++() {
      if (++inner_ == outer_->end()) {
        if (++outer_ == outer_end_) {
          inner_ = InnerIt();
        } else {
          inner_ = outer_->begin();
        }
      }
      return *this;
    }

    IteratorImpl operator++(int) const {
      IteratorImpl tmp = *this;
      ++tmp;
      return tmp;
    }

    bool operator==(const IteratorImpl& other) const {
      return outer_ == other.outer_ && inner_ == other.inner_;
    }

    bool operator!=(const IteratorImpl& other) const {
      return !(*this == other);
    }

   protected:
    OuterIt outer_;
    OuterIt outer_end_;
    InnerIt inner_;

   private:
    friend class Grid<T>;

    explicit IteratorImpl(OuterIt outer, OuterIt outer_end, InnerIt inner)
        : outer_(std::move(outer)),
          outer_end_(std::move(outer_end)),
          inner_(std::move(inner)) {}
  };

  // Iterator over the values held in cells.
  using iterator =
      IteratorImpl<T, typename std::vector<std::vector<T>>::iterator,
                   typename std::vector<T>::iterator>;

  // Const-iterator over the values held in cells.
  class const_iterator
      : public IteratorImpl<
            const T, typename std::vector<std::vector<T>>::const_iterator,
            typename std::vector<T>::const_iterator> {
   public:
    const_iterator(const iterator& it) {
      this->outer_ = it.outer_;
      this->outer_end_ = it.outer_end_;
      this->inner_ = it.inner_;
    }
  };

  // Iterator over Point coordinates in the grid.
  class PointIterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Point;
    using pointer = const Point*;
    using reference = const Point&;

    PointIterator() = default;

    reference operator*() const { return current_; }
    pointer operator->() const { return &current_; }

    PointIterator& operator++() {
      if (++current_.y == grid_->XSize()) {
        ++current_.x;
        current_.y = 0;
      }
      return *this;
    }

    PointIterator operator++(int) const {
      PointIterator tmp = *this;
      ++tmp;
      return tmp;
    }

    bool operator==(const PointIterator& other) const {
      return current_ == other.current_ && grid_ == other.grid_;
    }

    bool operator!=(const PointIterator& other) const {
      return !(*this == other);
    }

   private:
    friend class Grid;

    explicit PointIterator(Point current, const Grid* grid)
        : current_(current), grid_(grid) {}

    Point current_;
    const Grid* grid_ = nullptr;
  };

  // A view over the Points in range of the Grid.
  class PointView {
   public:
    using iterator = PointIterator;
    using const_iterator = PointIterator;

    const_iterator cbegin() const { return grid_->point_begin(); }

    const_iterator cend() const { return grid_->point_end(); }

    const_iterator begin() const { return cbegin(); }
    const_iterator end() const { return cend(); }

   private:
    friend class Grid;

    explicit PointView(const Grid* grid) : grid_(grid) {}

    const Grid* grid_ = nullptr;
  };

  explicit Grid(const size_t x_dim, const size_t y_dim)
      : values_(x_dim, std::vector<T>(y_dim)) {}

  // Parses a grid from a representation where each cell is an ASCII digit.
  // Technically X and Y are reversed from their intuitive direction in the
  // input, but thus far that hasn't been a problem.
  static Grid ReadFromDigits(const std::vector<std::string>& input) {
    Grid grid;
    for (const std::string& line : input) {
      grid.values_.emplace_back();
      for (const char c : line) {
        const T value = c - '0';
        CHECK(value >= 0);
        CHECK(value <= 9);
        grid.values_.back().emplace_back(value);
      }
    }
    return grid;
  }

  // Returns the grid's size in the X dimension.
  size_t XSize() const { return values_.size(); }

  // Returns the grid's size in the Y dimension.
  size_t YSize(size_t x = 0) const {
    CHECK(x < values_.size());
    return values_[x].size();
  }

  // Returns true iff `point` is within the range covered by the Grid.
  bool InRange(const Point& point) const {
    return point.x >= 0 && point.x < values_.size() && point.y >= 0 &&
           point.y < values_[point.x].size();
  }

  template <typename PointContainer>
  std::vector<Point> FilterInRange(const PointContainer& points) const {
    std::vector<Point> filtered;
    for (const Point& point : points) {
      if (InRange(point)) filtered.emplace_back(point);
    }
    return filtered;
  }

  // Accesses the value held in the cell at `point`.
  T& operator[](const Point& point) {
    CHECK(InRange(point));
    return values_[point.x][point.y];
  }

  // Accesses the value held in the cell at `point`.
  const T& operator[](const Point& point) const {
    CHECK(InRange(point));
    return values_[point.x][point.y];
  }

  // Returns a view of all the points in the Grid. Can be useful for iterating
  // over all grid points.
  PointView Points() const { return PointView(this); }

  // Iterators over values held in cells.
  iterator begin() {
    return iterator(values_.begin(), values_.end(), values_.front().begin());
  }

  iterator end() {
    return iterator(values_.end(), values_.end(),
                    typename std::vector<T>::iterator());
  }

  const_iterator cbegin() const {
    return const_iterator(const_cast<Grid*>(this)->begin());
  }

  const_iterator cend() const {
    return const_iterator(const_cast<Grid*>(this)->end());
  }

  const_iterator begin() const { return cbegin(); }

  const_iterator end() const { return cend(); }

  // Iterators over points on grid.
  PointIterator point_begin() const {
    return PointIterator(Points::kOrigin, this);
  }

  PointIterator point_end() const {
    return PointIterator({static_cast<int64_t>(XSize()), 0}, this);
  }

 private:
  Grid() = default;

  std::vector<std::vector<T>> values_;
};

}  // namespace aoc2021::grid2

#endif  // UTIL_GRID2_H_
