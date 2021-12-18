#include <cmath>
#include <memory>
#include <optional>
#include <variant>

#include "absl/strings/numbers.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/io.h"

namespace {

class Number {
 public:
  static std::unique_ptr<Number> Parse(absl::string_view rep,
                                       Number* parent = nullptr) {
    auto num = std::make_unique<Number>();
    num->parent_ = parent;

    CHECK(rep.front() == '[');
    int depth = 0;
    int comma = -1;
    for (int pos = 0; pos < rep.size(); ++pos) {
      switch (rep[pos]) {
        case '[':
          ++depth;
          break;
        case ']':
          --depth;
          break;
        case ',':
          if (depth == 1) {
            comma = pos;
          }
          break;
      }
    }
    CHECK(depth == 0);
    CHECK(comma != -1);

    absl::string_view left_rep = rep.substr(1, comma - 1);
    absl::string_view right_rep = rep.substr(comma + 1);
    right_rep.remove_suffix(1);

    int left_int = 0;
    if (absl::SimpleAtoi(left_rep, &left_int)) {
      num->left_ = left_int;
    } else {
      num->left_ = Parse(left_rep, num.get());
    }

    int right_int = 0;
    if (absl::SimpleAtoi(right_rep, &right_int)) {
      num->right_ = right_int;
    } else {
      num->right_ = Parse(right_rep, num.get());
    }

    return num;
  }

  static std::unique_ptr<Number> Add(std::unique_ptr<Number> left, std::unique_ptr<Number> right) {
    auto num = std::make_unique<Number>();
    left->parent_ = num.get();
    right->parent_ = num.get();
    num->left_ = std::move(left);
    num->right_ = std::move(right);
    return num;
  }

  void Reduce() {
    for (;;) {
      if (MaybeExplode(1)) continue;
      if (MaybeSplit()) continue;
      return;
    }
  }

  int Magnitude() const {
    int mag = 0;
    if (left_.index() == 0) {
      mag += 3 * std::get<0>(left_);
    } else {
      mag += 3 * std::get<1>(left_)->Magnitude();
    }
    if (right_.index() == 0) {
      mag += 2 * std::get<0>(right_);
    } else {
      mag += 2 * std::get<1>(right_)->Magnitude();
    }
    return mag;
  }

  void DebugPrint() const {
    std::cout << "[";
    if (left_.index() == 0) {
      std::cout << std::get<0>(left_);
    } else {
      std::get<1>(left_)->DebugPrint();
    }
    std::cout << ",";
    if (right_.index() == 0) {
      std::cout << std::get<0>(right_);
    } else {
      std::get<1>(right_)->DebugPrint();
    }
    std::cout << "]";
  }

 private:
  int Explode() const {
    CHECK(std::holds_alternative<int>(left_));
    CHECK(std::holds_alternative<int>(right_));

    int* left_of_pair = FindLeft();
    if (left_of_pair != nullptr) *left_of_pair += std::get<0>(left_);

    int* right_of_pair = FindRight();
    if (right_of_pair != nullptr) *right_of_pair += std::get<0>(right_);

    return 0;
  }

  std::unique_ptr<Number> Split(int value) {
    auto num = std::make_unique<Number>();
    num->parent_ = this;
    num->left_ = static_cast<int>(std::floor(static_cast<double>(value) / 2.0));
    num->right_ = static_cast<int>(std::ceil(static_cast<double>(value) / 2.0));
    return num;
  }

  bool MaybeExplode(int depth) {
    if (depth == 4) {
      if (left_.index() == 1) {
        left_ = std::get<1>(left_)->Explode();
        return true;
      }
      if (right_.index() == 1) {
        right_ = std::get<1>(right_)->Explode();
        return true;
      }
    }

    if (left_.index() == 1 && std::get<1>(left_)->MaybeExplode(depth + 1))
      return true;
    if (right_.index() == 1 && std::get<1>(right_)->MaybeExplode(depth + 1))
      return true;
    return false;
  }

  bool MaybeSplit() {
    if (left_.index() == 0 && std::get<0>(left_) >= 10) {
      left_ = Split(std::get<0>(left_));
      return true;
    }
    if (left_.index() == 1 && std::get<1>(left_)->MaybeSplit()) return true;

    if (right_.index() == 0 && std::get<0>(right_) >= 10) {
      right_ = Split(std::get<0>(right_));
      return true;
    }
    if (right_.index() == 1 && std::get<1>(right_)->MaybeSplit()) return true;

    return false;
  }

  int* LeftMost() {
    if (left_.index() == 0) return &std::get<0>(left_);
    return std::get<1>(left_)->LeftMost();
  }

  int* RightMost() {
    if (right_.index() == 0) return &std::get<0>(right_);
    return std::get<1>(right_)->RightMost();
  }

  int* FindLeft() const {
    if (parent_ == nullptr) return nullptr;
    if (parent_->left_.index() == 1 &&
        std::get<1>(parent_->left_).get() == this) {
      return parent_->FindLeft();
    }
    CHECK(parent_->right_.index() == 1 &&
          std::get<1>(parent_->right_).get() == this);
    if (parent_->left_.index() == 0) return &std::get<0>(parent_->left_);
    return std::get<1>(parent_->left_)->RightMost();
  }

  int* FindRight() const {
    if (parent_ == nullptr) return nullptr;
    if (parent_->right_.index() == 1 &&
        std::get<1>(parent_->right_).get() == this) {
      return parent_->FindRight();
    }
    CHECK(parent_->left_.index() == 1 &&
          std::get<1>(parent_->left_).get() == this);
    if (parent_->right_.index() == 0) return &std::get<0>(parent_->right_);
    return std::get<1>(parent_->right_)->LeftMost();
  }

  Number* parent_ = nullptr;
  std::variant<int, std::unique_ptr<Number>> left_;
  std::variant<int, std::unique_ptr<Number>> right_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);

  std::unique_ptr<Number> result = Number::Parse(lines.front());
  for (auto line_it = lines.begin() + 1; line_it != lines.end(); ++line_it) {
    std::unique_ptr<Number> next = Number::Parse(*line_it);
    result = Number::Add(std::move(result), std::move(next));
    result->Reduce();
  }

  std::cout << result->Magnitude() << "\n";

  return 0;
}
