#include <bits/stdint-intn.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/dimension_grid.h"
#include "util/io.h"

namespace {
/*
using MachineState = aoc2021::DimensionGrid<4>::Point;

size_t Coord(char var) {
  switch (var) {
    case 'w':
      return 0;
    case 'x':
      return 1;
    case 'y':
      return 2;
    case 'z':
      return 3;
    default:
      CHECK_FAIL();
  }
}

MachineState ExecInstruction(const MachineState& before, absl::string_view inst,
                             std::deque<int>& inputs) {
  const absl::string_view mnemonic = inst.substr(0, 3);
  const char reg_0 = inst[4];

  if (mnemonic == "inp") {
    CHECK(!inputs.empty());
    MachineState after(before);
    after.coords[Coord(reg_0)] = inputs.front();
    inputs.pop_front();
    return after;
  }

  int64_t arg_1_value = 0;
  if (!absl::SimpleAtoi(inst.substr(5), &arg_1_value)) {
    CHECK(inst.size() == 7);
    arg_1_value = before.coords[Coord(inst.back())];
  }

  if (mnemonic == "add") {
    MachineState after(before);
    after.coords[Coord(reg_0)] = before.coords[Coord(reg_0)] + arg_1_value;
    return after;
  }

  if (mnemonic == "mul") {
    MachineState after(before);
    after.coords[Coord(reg_0)] = before.coords[Coord(reg_0)] * arg_1_value;
    return after;
  }

  if (mnemonic == "div") {
    CHECK(arg_1_value != 0);
    MachineState after(before);
    after.coords[Coord(reg_0)] = before.coords[Coord(reg_0)] / arg_1_value;
    return after;
  }

  if (mnemonic == "mod") {
    CHECK(arg_1_value != 0);
    MachineState after(before);
    after.coords[Coord(reg_0)] = before.coords[Coord(reg_0)] % arg_1_value;
    return after;
  }

  if (mnemonic == "eql") {
    MachineState after(before);
    after.coords[Coord(reg_0)] =
        static_cast<int64_t>(before.coords[Coord(reg_0)] == arg_1_value);
    return after;
  }

  CHECK_FAIL();
}

// Returns true iff valid, i.e. z == 0 at end of execution.
bool RunProgram(const std::vector<std::string>& program,
                std::deque<int> inputs) {
  MachineState state;
  for (absl::string_view inst : program) {
    state = ExecInstruction(state, inst, inputs);
  }
  return state.coords[Coord('z')] == 0;
}

std::deque<int> GetNextLowerInputs(std::deque<int> inputs) {
  for (auto it = inputs.rbegin(); it != inputs.rend(); ++it) {
    if (--*it > 0) return inputs;
    *it = 9;
  }
  CHECK_FAIL();
}

std::deque<int> FindLargestModelNum(const std::vector<std::string>& program) {
  std::deque<int> inputs(14, 9);
  int64_t counter = 0;
  while (!RunProgram(program, inputs)) {
    if ((++counter % 1000) == 0) {
      std::cout << "Tried " << counter << " model numbers\n";
    }
    inputs = GetNextLowerInputs(inputs);
  }
  return inputs;
}
*/

constexpr std::array<int64_t, 14> kTens{
  INT64_C(10000000000000),
  INT64_C(1000000000000),
  INT64_C(100000000000),
  INT64_C(10000000000),
  INT64_C(1000000000),
  INT64_C(100000000),
  INT64_C(10000000),
  INT64_C(1000000),
  INT64_C(100000),
  INT64_C(10000),
  INT64_C(1000),
  INT64_C(100),
  INT64_C(10),
  INT64_C(1)};

struct Num {
  std::array<int8_t, 14> digits;

  Num& operator--() {
    for (int pos = 13; pos >= 0; --pos) {
      if (--digits[pos] >= 1) return *this;
      digits[pos] = 9;
    }
    CHECK_FAIL();
  }

  Num Mask(uint16_t mask) const {
    Num masked = *this;
    for (int pos = 0; pos < 14; ++pos) {
      if ((mask & (1 << pos)) == 0) {
        masked.digits[pos] = 0;
      }
    }
    return masked;
  }

  bool operator==(const Num& other) const { return digits == other.digits; }

  template <typename H>
  friend H AbslHashValue(H h, const Num& num) {
    return H::combine(std::move(h), num.digits);
  }
};

std::ostream& operator<<(std::ostream& os, const Num& num) {
  for (int8_t digit : num.digits) {
    os << static_cast<int>(digit);
  }
  return os;
}

constexpr Num kLargest{.digits = {9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9}};

class LiteralExpr;
class InputExpr;
class AddExpr;
class MulExpr;
class DivExpr;
class ModExpr;
class EqlExpr;

class Expr {
 public:
  virtual ~Expr() = default;

  //  virtual std::unique_ptr<Expr> Clone() const = 0;

  virtual std::string DebugStr(bool paren = false) const = 0;

  virtual const LiteralExpr* AsLiteral() const { return nullptr; }
  virtual const InputExpr* AsInput() const { return nullptr; }
  virtual const AddExpr* AsAdd() const { return nullptr; }
  virtual const MulExpr* AsMul() const { return nullptr; }
  virtual const DivExpr* AsDiv() const { return nullptr; }
  virtual const ModExpr* AsMod() const { return nullptr; }
  virtual const EqlExpr* AsEql() const { return nullptr; }

  virtual std::shared_ptr<Expr> Simplify() const { return nullptr; }

  virtual bool TreeEqual(const Expr& other) const = 0;

  int64_t Eval(const Num& input) const {
    Num masked = input.Mask(digit_mask_);
    auto [it, inserted] =
        masked_input_to_value_.try_emplace(std::move(masked), -1);
    if (inserted) {
      it->second = EvalImpl(input);
    }
    return it->second;
  }

  virtual int64_t EvalImpl(const Num& input) const = 0;

  virtual void CollectChildren(
      absl::flat_hash_set<std::shared_ptr<Expr>>& all_nodes) const {}

  std::uint16_t DigitMask() const { return digit_mask_; }

  virtual std::optional<std::pair<int64_t, int64_t>> ValueRange() const = 0;

  const absl::flat_hash_set<int64_t>& PossibleValues() const {
    if (possible_values_memo_.has_value()) return *possible_values_memo_;
    possible_values_memo_ = PossibleValuesImpl();
    return *possible_values_memo_;
  }

  virtual absl::flat_hash_set<int64_t> PossibleValuesImpl() const = 0;

  void PurgePossibleValuesMemo() { possible_values_memo_ = std::nullopt; }

 protected:
  mutable absl::flat_hash_map<Num, int64_t> masked_input_to_value_;

  mutable std::optional<absl::flat_hash_set<int64_t>> possible_values_memo_;

  std::uint16_t digit_mask_ = 0;
};

void TrySimplify(std::shared_ptr<Expr>& expr) {
//  return;
  std::shared_ptr<Expr> simplified = expr->Simplify();
  if (simplified != nullptr) {
    expr = std::move(simplified);
  }
}

class LiteralExpr : public Expr {
 public:
  explicit LiteralExpr(int64_t value) : value_(value) {}

  ~LiteralExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    return absl::StrCat(value_);
  }

  const LiteralExpr* AsLiteral() const override { return this; }

  bool TreeEqual(const Expr& other) const override {
    const LiteralExpr* other_lit = other.AsLiteral();
    return (other_lit != nullptr && value_ == other_lit->value_);
  }

  int64_t value() const { return value_; }

  int64_t EvalImpl(const Num& input) const override { return value_; }

  std::optional<std::pair<int64_t, int64_t>> ValueRange() const override {
    return std::make_pair(value_, value_);
  }

  absl::flat_hash_set<int64_t> PossibleValuesImpl() const override {
    return {value_};
  }

 private:
  int64_t value_ = 0;
};

class InputExpr : public Expr {
 public:
  explicit InputExpr(int input_idx) : input_idx_(input_idx) {
    digit_mask_ = 1 << input_idx_;
  }

  ~InputExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    return absl::StrCat("input_", input_idx_);
  }

  virtual const InputExpr* AsInput() const override { return this; }

  bool TreeEqual(const Expr& other) const override {
    const InputExpr* other_input = other.AsInput();
    return (other_input != nullptr && input_idx_ == other_input->input_idx_);
  }

  int64_t EvalImpl(const Num& input) const override {
    return input.digits[input_idx_];
  }

  std::optional<std::pair<int64_t, int64_t>> ValueRange() const override {
    return std::pair<int64_t, int64_t>(1, 9);
  }

  absl::flat_hash_set<int64_t> PossibleValuesImpl() const override {
    if (tentative_value_.has_value()) return {*tentative_value_};
    return {1, 2, 3, 4, 5, 6, 7, 8, 9};
  }

  void SetTentativeValue(std::optional<int64_t> tentative_value) {
    tentative_value_ = tentative_value;
  }

  int64_t FinalizeTentative() const { return tentative_value_.value_or(9); }

 private:
  int input_idx_ = 0;
  std::optional<int64_t> tentative_value_;
};

class AddExpr : public Expr {
 public:
  AddExpr(std::shared_ptr<Expr> left, std::shared_ptr<Expr> right) {
    digit_mask_ = left->DigitMask() | right->DigitMask();

    const LiteralExpr* left_lit = left->AsLiteral();
    const LiteralExpr* right_lit = right->AsLiteral();
    if (left_lit != nullptr) {
      if (right_lit != nullptr) {
        literal_ = std::make_shared<LiteralExpr>(left_lit->value() +
                                                 right_lit->value());
        return;
      }
      literal_ = std::dynamic_pointer_cast<LiteralExpr>(std::move(left));
      non_literal_.emplace_back(std::move(right));
      return;
    }
    if (right_lit != nullptr) {
      literal_ = std::dynamic_pointer_cast<LiteralExpr>(std::move(right));
      non_literal_.emplace_back(std::move(left));
      return;
    }
    non_literal_.emplace_back(std::move(left));
    non_literal_.emplace_back(std::move(right));
  }

  ~AddExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    std::string debugstr;
    if (paren) debugstr.push_back('(');
    if (literal_ != nullptr) {
      absl::StrAppend(&debugstr, literal_->DebugStr(true), " + ");
    }
    for (auto it = non_literal_.begin(); it != non_literal_.end(); ++it) {
      absl::StrAppend(&debugstr, (*it)->DebugStr(true));
      if (it != non_literal_.end() - 1) {
        absl::StrAppend(&debugstr, " + ");
      }
    }
    if (paren) debugstr.push_back(')');
    return debugstr;
  }

  virtual const AddExpr* AsAdd() const override { return this; }

  std::shared_ptr<Expr> Simplify() const override {
    if (literal_ != nullptr && non_literal_.empty()) {
      return literal_;
    }

    std::vector<std::shared_ptr<Expr>> non_add_args;
    std::vector<const AddExpr*> add_args;
    for (const std::shared_ptr<Expr>& non_literal : non_literal_) {
      const AddExpr* add_arg = non_literal->AsAdd();
      if (add_arg != nullptr) {
        add_args.emplace_back(add_arg);
      } else {
        non_add_args.emplace_back(non_literal);
      }
    }

    if (add_args.empty()) {
      if (literal_ != nullptr && literal_->value() == 0) {
        switch (non_add_args.size()) {
          case 0:
            CHECK_FAIL();
          case 1:
            return non_add_args.front();
          default:
            auto reduced_add = std::shared_ptr<AddExpr>(new AddExpr);
            reduced_add->non_literal_ = std::move(non_add_args);
            reduced_add->CalculateDigitMask();
            return reduced_add;
        }
      }
      return nullptr;
    }

    auto fused_add = std::shared_ptr<AddExpr>(new AddExpr);
    fused_add->digit_mask_ = digit_mask_;
    fused_add->literal_ = literal_;
    fused_add->non_literal_ = std::move(non_add_args);
    for (const AddExpr* add_arg : add_args) {
      if (add_arg->literal_ != nullptr) {
        if (fused_add->literal_ == nullptr) {
          fused_add->literal_ = add_arg->literal_;
        } else {
          fused_add->literal_ = std::make_shared<LiteralExpr>(
              fused_add->literal_->value() + add_arg->literal_->value());
        }
      }
      fused_add->non_literal_.insert(fused_add->non_literal_.end(),
                                     add_arg->non_literal_.begin(),
                                     add_arg->non_literal_.end());
    }
    if (fused_add->literal_ != nullptr && fused_add->literal_->value() == 0) {
      fused_add->literal_ = nullptr;
    }
    fused_add->CalculateDigitMask();
    return fused_add;
  }

  bool TreeEqual(const Expr& other) const override {
    if (this == &other) return true;

    const AddExpr* other_add = other.AsAdd();
    if (other_add == nullptr) return false;

    if (literal_ != nullptr) {
      if (other_add->literal_ == nullptr) return false;
      if (!literal_->TreeEqual(*other_add->literal_)) return false;
    } else {
      if (other_add->literal_ != nullptr) return false;
    }

    if (non_literal_.size() != other_add->non_literal_.size()) return false;

    // Could try all permutations, for now just compare in order.
    for (size_t i = 0; i < non_literal_.size(); ++i) {
      if (!non_literal_[i]->TreeEqual(*other_add->non_literal_[i]))
        return false;
    }
    return true;
  }

  int64_t EvalImpl(const Num& input) const override {
    int64_t sum = literal_ == nullptr ? 0 : literal_->Eval(input);
    for (const auto& non_literal : non_literal_) {
      sum += non_literal->Eval(input);
    }
    return sum;
  }

  void CollectChildren(
      absl::flat_hash_set<std::shared_ptr<Expr>>& all_nodes) const override {
    if (literal_ != nullptr) {
      if (all_nodes.emplace(literal_).second) {
        literal_->CollectChildren(all_nodes);
      }
    }
    for (const auto& non_literal : non_literal_) {
      if (all_nodes.emplace(non_literal).second) {
        non_literal->CollectChildren(all_nodes);
      }
    }
  }

  std::optional<std::pair<int64_t, int64_t>> ValueRange() const override {
    std::int64_t min_sum = 0;
    std::int64_t max_sum = 0;
    if (literal_) {
      min_sum = max_sum = literal_->value();
    }
    for (const auto& non_literal : non_literal_) {
      auto range = non_literal->ValueRange();
      if (!range.has_value()) return std::nullopt;
      min_sum += range->first;
      max_sum += range->second;
    }
    CHECK(min_sum <= max_sum);
    return std::make_pair(min_sum, max_sum);
  }

  absl::flat_hash_set<int64_t> PossibleValuesImpl() const override {
    absl::flat_hash_set<int64_t> sums = {
        literal_ == nullptr ? 0 : literal_->value()};
    for (const auto& non_literal : non_literal_) {
      absl::flat_hash_set<int64_t> sums_expanded;
      for (const int64_t val : non_literal->PossibleValues()) {
        for (const int64_t previous_sum : sums) {
          sums_expanded.insert(previous_sum + val);
        }
      }
      sums = std::move(sums_expanded);
    }
    return sums;
  }

  std::shared_ptr<Expr> DistributeMul(int64_t factor) const;
  std::shared_ptr<Expr> DistributeDiv(int64_t factor) const;
  std::shared_ptr<Expr> DistributeMod(int64_t mod) const;

 private:
  AddExpr() = default;

  void CalculateDigitMask() {
    digit_mask_ = 0;
    for (const auto& non_literal : non_literal_) {
      digit_mask_ |= non_literal->DigitMask();
    }
  }

  std::shared_ptr<LiteralExpr> literal_;  // May be null.
  std::vector<std::shared_ptr<Expr>> non_literal_;
};

class MulExpr : public Expr {
 public:
  MulExpr(std::shared_ptr<Expr> left, std::shared_ptr<Expr> right) {
    digit_mask_ = left->DigitMask() | right->DigitMask();

    const LiteralExpr* left_lit = left->AsLiteral();
    const LiteralExpr* right_lit = right->AsLiteral();
    if (left_lit != nullptr) {
      if (right_lit != nullptr) {
        literal_ = std::make_shared<LiteralExpr>(left_lit->value() *
                                                 right_lit->value());
        return;
      }
      literal_ = std::dynamic_pointer_cast<LiteralExpr>(std::move(left));
      non_literal_.emplace_back(std::move(right));
      return;
    }
    if (right_lit != nullptr) {
      literal_ = std::dynamic_pointer_cast<LiteralExpr>(std::move(right));
      non_literal_.emplace_back(std::move(left));
      return;
    }
    non_literal_.emplace_back(std::move(left));
    non_literal_.emplace_back(std::move(right));
  }

  ~MulExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    std::string debugstr;
    if (paren) debugstr.push_back('(');
    if (literal_ != nullptr) {
      absl::StrAppend(&debugstr, literal_->DebugStr(true), " * ");
    }
    for (auto it = non_literal_.begin(); it != non_literal_.end(); ++it) {
      absl::StrAppend(&debugstr, (*it)->DebugStr(true));
      if (it != non_literal_.end() - 1) {
        absl::StrAppend(&debugstr, " * ");
      }
    }
    if (paren) debugstr.push_back(')');
    return debugstr;
  }

  virtual const MulExpr* AsMul() const override { return this; }

  std::shared_ptr<Expr> Simplify() const override {
    if (literal_ != nullptr) {
      if (literal_->value() == 0 || non_literal_.empty()) return literal_;

      // Distribute multiplication by a constant over addition.
      if (non_literal_.size() == 1) {
        const AddExpr* right_add = non_literal_.front()->AsAdd();
        if (right_add != nullptr) {
          return right_add->DistributeMul(literal_->value());
        }
      }
    }

    std::vector<std::shared_ptr<Expr>> non_mul_args;
    std::vector<const MulExpr*> mul_args;
    for (const std::shared_ptr<Expr>& non_literal : non_literal_) {
      const MulExpr* mul_arg = non_literal->AsMul();
      if (mul_arg != nullptr) {
        mul_args.emplace_back(mul_arg);
      } else {
        non_mul_args.emplace_back(non_literal);
      }
    }

    if (mul_args.empty()) {
      if (literal_ != nullptr && literal_->value() == 1) {
        switch (non_mul_args.size()) {
          case 0:
            CHECK_FAIL();
          case 1:
            return non_mul_args.front();
          default:
            auto reduced_mul = std::shared_ptr<MulExpr>(new MulExpr);
            reduced_mul->non_literal_ = std::move(non_mul_args);
            reduced_mul->CalculateDigitMask();
            return reduced_mul;
        }
      }
      return nullptr;
    }

    auto fused_mul = std::shared_ptr<MulExpr>(new MulExpr);
    fused_mul->digit_mask_ = digit_mask_;
    fused_mul->literal_ = literal_;
    fused_mul->non_literal_ = std::move(non_mul_args);
    for (const MulExpr* mul_arg : mul_args) {
      if (mul_arg->literal_ != nullptr) {
        if (fused_mul->literal_ == nullptr) {
          fused_mul->literal_ = mul_arg->literal_;
        } else {
          fused_mul->literal_ = std::make_shared<LiteralExpr>(
              fused_mul->literal_->value() * mul_arg->literal_->value());
        }
      }
      fused_mul->non_literal_.insert(fused_mul->non_literal_.end(),
                                     mul_arg->non_literal_.begin(),
                                     mul_arg->non_literal_.end());
    }
    if (fused_mul->literal_ != nullptr) {
      if (fused_mul->literal_->value() == 0) {
        return fused_mul->literal_;
      }
      if (fused_mul->literal_->value() == 1) {
        fused_mul->literal_ = nullptr;
      }
    }
    fused_mul->CalculateDigitMask();
    return fused_mul;
  }

  bool TreeEqual(const Expr& other) const override {
    if (this == &other) return true;

    const MulExpr* other_mul = other.AsMul();
    if (other_mul == nullptr) return false;

    if (literal_ != nullptr) {
      if (other_mul->literal_ == nullptr) return false;
      if (!literal_->TreeEqual(*other_mul->literal_)) return false;
    } else {
      if (other_mul->literal_ != nullptr) return false;
    }

    if (non_literal_.size() != other_mul->non_literal_.size()) return false;

    // Could try all permutations, for now just compare in order.
    for (size_t i = 0; i < non_literal_.size(); ++i) {
      if (!non_literal_[i]->TreeEqual(*other_mul->non_literal_[i]))
        return false;
    }
    return true;
  }

  int64_t EvalImpl(const Num& input) const override {
    int64_t prod = literal_ == nullptr ? 1 : literal_->Eval(input);
    for (const auto& non_literal : non_literal_) {
      prod *= non_literal->Eval(input);
    }
    return prod;
  }

  void CollectChildren(
      absl::flat_hash_set<std::shared_ptr<Expr>>& all_nodes) const override {
    if (literal_ != nullptr) {
      if (all_nodes.emplace(literal_).second) {
        literal_->CollectChildren(all_nodes);
      }
    }
    for (const auto& non_literal : non_literal_) {
      if (all_nodes.emplace(non_literal).second) {
        non_literal->CollectChildren(all_nodes);
      }
    }
  }

  std::optional<std::pair<int64_t, int64_t>> ValueRange() const override {
    int64_t min_prod = 1;
    int64_t max_prod = 1;

    if (literal_ != nullptr) {
      min_prod = literal_->value();
      max_prod = literal_->value();
    }
    for (const auto& non_literal : non_literal_) {
      auto non_literal_range = non_literal->ValueRange();
      if (!non_literal_range.has_value()) return std::nullopt;
      std::array<int64_t, 4> prods = {min_prod * non_literal_range->first,
                                      min_prod * non_literal_range->second,
                                      max_prod * non_literal_range->first,
                                      max_prod * non_literal_range->second};
      auto minmax_elt = std::minmax_element(prods.begin(), prods.end());
      min_prod = *minmax_elt.first;
      max_prod = *minmax_elt.second;
    }
    return std::make_pair(min_prod, max_prod);
  }

  absl::flat_hash_set<int64_t> PossibleValuesImpl() const override {
    absl::flat_hash_set<int64_t> prods = {
        literal_ == nullptr ? 1 : literal_->value()};
    for (const auto& non_literal : non_literal_) {
      absl::flat_hash_set<int64_t> prods_expanded;
      for (const int64_t val : non_literal->PossibleValues()) {
        for (const int64_t previous_prod : prods) {
          prods_expanded.insert(previous_prod * val);
        }
      }
      prods = std::move(prods_expanded);
    }
    return prods;
  }

  std::shared_ptr<Expr> DivideLiteral(int64_t divisor) const {
    if (literal_ == nullptr) return nullptr;
    if (literal_->value() % divisor != 0) return nullptr;
    const int64_t new_lit = literal_->value() / divisor;

    std::shared_ptr<Expr> new_expr;
    switch (new_lit) {
      case 0:
        new_expr = std::make_shared<LiteralExpr>(0);
        break;
      case 1:
        if (non_literal_.size() == 1) {
          new_expr = non_literal_.front();
        } else {
          std::shared_ptr<MulExpr> new_mul(new MulExpr);
          new_mul->non_literal_ = non_literal_;
          new_mul->CalculateDigitMask();
          new_expr = std::move(new_mul);
        }
        break;
      default: {
        std::shared_ptr<MulExpr> new_mul(new MulExpr);
        new_mul->literal_ = std::make_shared<LiteralExpr>(new_lit);
        new_mul->non_literal_ = non_literal_;
        new_mul->CalculateDigitMask();
        new_expr = std::move(new_mul);
      }
    }
    CHECK(new_expr != nullptr);
    std::shared_ptr<Expr> new_expr_simplified = new_expr->Simplify();
    if (new_expr_simplified != nullptr) return new_expr_simplified;
    return new_expr;
  }

  const std::shared_ptr<LiteralExpr>& LiteralFactor() const { return literal_; }

 private:
  MulExpr() = default;

  void CalculateDigitMask() {
    digit_mask_ = 0;
    for (const auto& non_literal : non_literal_) {
      digit_mask_ |= non_literal->DigitMask();
    }
  }

  std::shared_ptr<LiteralExpr> literal_;  // May be null.
  std::vector<std::shared_ptr<Expr>> non_literal_;
};

std::shared_ptr<Expr> AddExpr::DistributeMul(int64_t factor) const {
  std::shared_ptr<AddExpr> distributed(new AddExpr);

  if (literal_ != nullptr) {
    distributed->literal_ =
        std::make_shared<LiteralExpr>(literal_->value() * factor);
  }

  for (const auto& non_literal : non_literal_) {
    auto mul = std::make_shared<MulExpr>(non_literal,
                                         std::make_shared<LiteralExpr>(factor));
    std::shared_ptr<Expr> mul_simplified = mul->Simplify();
    if (mul_simplified != nullptr) {
      const LiteralExpr* mul_simplified_lit = mul_simplified->AsLiteral();
      if (mul_simplified_lit != nullptr) {
        const int64_t original_lit = distributed->literal_ == nullptr
                                         ? 0
                                         : distributed->literal_->value();
        const int64_t sum = original_lit + mul_simplified_lit->value();
        if (sum == 0) {
          distributed->literal_ = nullptr;
        } else {
          distributed->literal_ = std::make_shared<LiteralExpr>(sum);
        }
      } else {
        distributed->non_literal_.emplace_back(std::move(mul_simplified));
      }
    } else {
      distributed->non_literal_.emplace_back(std::move(mul));
    }
  }
  distributed->CalculateDigitMask();

  std::shared_ptr<Expr> distributed_simplified = distributed->Simplify();
  if (distributed_simplified != nullptr) return distributed_simplified;
  return distributed;
}

class DivExpr : public Expr {
 public:
  DivExpr(std::shared_ptr<Expr> left, std::shared_ptr<Expr> right)
      : left_(std::move(left)), right_(std::move(right)) {
    digit_mask_ = left_->DigitMask() | right_->DigitMask();
  }

  ~DivExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    if (paren) {
      return absl::StrCat("(", left_->DebugStr(true), " / ",
                          right_->DebugStr(true), ")");
    }
    return absl::StrCat(left_->DebugStr(true), " / ", right_->DebugStr(true));
  }

  virtual const DivExpr* AsDiv() const override { return this; }

  std::shared_ptr<Expr> Simplify() const override {
    if (left_->TreeEqual(*right_)) {
      return std::make_shared<LiteralExpr>(1);
    }

    const LiteralExpr* left_lit = left_->AsLiteral();
    const LiteralExpr* right_lit = right_->AsLiteral();
    if (left_lit != nullptr) {
      if (right_lit != nullptr) {
        return std::make_shared<LiteralExpr>(left_lit->value() /
                                             right_lit->value());
      }
    }
    if (right_lit != nullptr) {
      if (right_lit->value() == 0) {
        CHECK_FAIL();
      } else if (right_lit->value() == 1) {
        return left_;
      }
      const MulExpr* left_mul = left_->AsMul();
      if (left_mul != nullptr) {
        std::shared_ptr<Expr> divided =
            left_mul->DivideLiteral(right_lit->value());
        if (divided != nullptr) return divided;
      }
      const AddExpr* left_add = left_->AsAdd();
      if (left_add != nullptr) {
        return left_add->DistributeDiv(right_lit->value());
      }
    }

    const MulExpr* left_mul = left_->AsMul();

    const DivExpr* left_div = left_->AsDiv();
    if (left_div != nullptr) {
      std::shared_ptr<Expr> new_divisor =
          std::make_shared<MulExpr>(left_div->right_, right_);
      TrySimplify(new_divisor);
      return std::make_shared<DivExpr>(left_div->left_, std::move(new_divisor));
    }
    return nullptr;
  }

  bool TreeEqual(const Expr& other) const override {
    const DivExpr* other_div = other.AsDiv();
    if (other_div == nullptr) return false;
    if (left_->TreeEqual(*other_div->left_) &&
        right_->TreeEqual(*other_div->right_)) {
      return true;
    }
    return false;
  }

  int64_t EvalImpl(const Num& input) const override {
    return left_->Eval(input) / right_->Eval(input);
  }

  void CollectChildren(
      absl::flat_hash_set<std::shared_ptr<Expr>>& all_nodes) const override {
    if (all_nodes.emplace(left_).second) {
      left_->CollectChildren(all_nodes);
    }
    if (all_nodes.emplace(right_).second) {
      right_->CollectChildren(all_nodes);
    }
  }

  std::optional<std::pair<int64_t, int64_t>> ValueRange() const override {
    /*
    auto left_range = left_->ValueRange();
    if (!left_range.has_value()) return std::nullopt;

    auto right_range = right_->ValueRange();
    if (!right_range.has_value()) return std::nullopt;
    if (right_range->first == 0) ++right_range->first;
    if (right_range->second == 0) --right_range->second;
    CHECK(right_range->first <= right_range->second);
    */

    // TODO... need to be careful of signs.
    return std::nullopt;
  }

  absl::flat_hash_set<int64_t> PossibleValuesImpl() const override {
    absl::flat_hash_set<int64_t> quots;
    for (int64_t left : left_->PossibleValues()) {
      for (int64_t right : right_->PossibleValues()) {
        if (right == 0) continue;
        quots.insert(left / right);
      }
    }
    return quots;
  }

 private:
  std::shared_ptr<Expr> left_;
  std::shared_ptr<Expr> right_;
};

std::shared_ptr<Expr> AddExpr::DistributeDiv(int64_t factor) const {
  std::vector<std::shared_ptr<LiteralExpr>> literal_terms;
  std::vector<std::shared_ptr<Expr>> divided_terms;
  std::vector<std::shared_ptr<Expr>> indivisible_terms;
  if (literal_ != nullptr) {
    literal_terms.emplace_back(literal_);
  }

  for (const auto& non_literal : non_literal_) {
    auto non_literal_mul = non_literal->AsMul();
    if (non_literal_mul == nullptr) {
      indivisible_terms.emplace_back(non_literal);
      continue;
    }
    std::shared_ptr<Expr> divided = non_literal_mul->DivideLiteral(factor);
    if (divided == nullptr) {
      indivisible_terms.emplace_back(non_literal);
      continue;
    }
    if (divided->AsLiteral() != nullptr) {
      literal_terms.emplace_back(
          std::dynamic_pointer_cast<LiteralExpr>(std::move(divided)));
    } else {
      divided_terms.emplace_back(divided);
    }
  }

  int64_t literal_sum = 0;
  for (const auto& literal_term : literal_terms) {
    literal_sum += literal_term->value();
  }

  std::shared_ptr<LiteralExpr> divided_literal;
  std::shared_ptr<LiteralExpr> pushdown_literal;
  if (literal_sum == 0) {
    // Do nothing.
  } else if (literal_sum % factor == 0 || indivisible_terms.empty()) {
    divided_literal = std::make_shared<LiteralExpr>(literal_sum / factor);
  } else {
    pushdown_literal = std::make_shared<LiteralExpr>(literal_sum);
  }

  std::shared_ptr<AddExpr> pushdown_add(new AddExpr);
  pushdown_add->literal_ = pushdown_literal;
  pushdown_add->non_literal_ = std::move(indivisible_terms);
  pushdown_add->CalculateDigitMask();
  std::shared_ptr<Expr> pushdown_add_simplified = pushdown_add->Simplify();
  if (pushdown_add_simplified == nullptr) {
    pushdown_add_simplified = pushdown_add;
  }

  auto pushdown_div =
      std::make_shared<DivExpr>(std::move(pushdown_add_simplified),
                                std::make_shared<LiteralExpr>(factor));
  divided_terms.emplace_back(std::move(pushdown_div));

  std::shared_ptr<AddExpr> distributed(new AddExpr);
  distributed->literal_ = std::move(divided_literal);
  distributed->non_literal_ = std::move(divided_terms);
  distributed->CalculateDigitMask();
  std::shared_ptr<Expr> distributed_simplified = distributed->Simplify();
  if (distributed_simplified != nullptr) return distributed_simplified;
  return distributed;
}

class ModExpr : public Expr {
 public:
  ModExpr(std::shared_ptr<Expr> left, std::shared_ptr<Expr> right)
      : left_(std::move(left)), right_(std::move(right)) {
    digit_mask_ = left_->DigitMask() | right_->DigitMask();
  }

  ~ModExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    if (paren) {
      return absl::StrCat("(", left_->DebugStr(true), " % ",
                          right_->DebugStr(true), ")");
    }
    return absl::StrCat(left_->DebugStr(true), " % ", right_->DebugStr(true));
  }

  virtual const ModExpr* AsMod() const override { return this; }

  std::shared_ptr<Expr> Simplify() const override {
    const LiteralExpr* left_lit = left_->AsLiteral();
    const LiteralExpr* right_lit = right_->AsLiteral();
    if (left_lit != nullptr) {
      if (right_lit != nullptr) {
        return std::make_shared<LiteralExpr>(left_lit->value() %
                                             right_lit->value());
      }
      if (left_lit->value() == 0) {
        return std::make_shared<LiteralExpr>(0);
      }
    }
    if (right_lit != nullptr) {
      if (right_lit->value() == 0) {
        CHECK_FAIL();
      } else if (right_lit->value() == 1) {
        return std::make_shared<LiteralExpr>(0);
      }

      const MulExpr* left_mul = left_->AsMul();
      if (left_mul != nullptr) {
        const auto& left_mul_literal_factor = left_mul->LiteralFactor();
        if (left_mul_literal_factor != nullptr) {
          if (left_mul_literal_factor->value() % right_lit->value() == 0) {
            return std::make_shared<LiteralExpr>(0);
          }
        }
      }

      const AddExpr* left_add = left_->AsAdd();
      if (left_add != nullptr) {
        return left_add->DistributeMod(right_lit->value());
      }
    }
    return nullptr;
  }

  bool TreeEqual(const Expr& other) const override {
    const ModExpr* other_mod = other.AsMod();
    if (other_mod == nullptr) return false;
    if (left_->TreeEqual(*other_mod->left_) &&
        right_->TreeEqual(*other_mod->right_)) {
      return true;
    }
    return false;
  }

  int64_t EvalImpl(const Num& input) const override {
    return left_->Eval(input) % right_->Eval(input);
  }

  void CollectChildren(
      absl::flat_hash_set<std::shared_ptr<Expr>>& all_nodes) const override {
    if (all_nodes.emplace(left_).second) {
      left_->CollectChildren(all_nodes);
    }
    if (all_nodes.emplace(right_).second) {
      right_->CollectChildren(all_nodes);
    }
  }

  std::optional<std::pair<int64_t, int64_t>> ValueRange() const override {
    // TODO: could narrow range if left range is fully within mod.
    auto right_range = right_->ValueRange();
    if (!right_range.has_value()) return std::nullopt;

    auto left_range = left_->ValueRange();
    if (!left_range.has_value() || left_range->second >= right_range->second) {
      return std::make_pair(int64_t{0}, right_range->second - 1);
    }
    return std::make_pair(int64_t{0}, left_range->second);
  }

  absl::flat_hash_set<int64_t> PossibleValuesImpl() const override {
    absl::flat_hash_set<int64_t> mods;
    for (int64_t left : left_->PossibleValues()) {
      if (left < 0) continue;
      for (int64_t right : right_->PossibleValues()) {
        if (right <= 0) continue;
        mods.insert(left % right);
      }
    }
    return mods;
  }

 private:
  std::shared_ptr<Expr> left_;
  std::shared_ptr<Expr> right_;
};

std::shared_ptr<Expr> AddExpr::DistributeMod(int64_t mod) const {
  int64_t lit_sum = 0;
  if (literal_ != nullptr) {
    lit_sum += literal_->value() % mod;
  }

  auto literal_mod = std::make_shared<LiteralExpr>(mod);
  std::vector<std::shared_ptr<Expr>> non_literal_modded;
  for (const auto& non_literal : non_literal_) {
    auto modded = std::make_shared<ModExpr>(non_literal, literal_mod);
    auto modded_simplified = modded->Simplify();
    if (modded_simplified != nullptr) {
      const LiteralExpr* modded_simplified_lit = modded_simplified->AsLiteral();
      if (modded_simplified_lit != nullptr) {
        lit_sum += modded_simplified_lit->value();
      } else {
        non_literal_modded.emplace_back(std::move(modded_simplified));
      }
    } else {
      non_literal_modded.emplace_back(std::move(modded));
    }
  }

  lit_sum %= mod;
  if (non_literal_modded.empty()) {
    return std::make_shared<LiteralExpr>(lit_sum);
  }

  if (lit_sum == 0 && non_literal_modded.size() == 1) {
    return non_literal_modded.front();
  }

  std::shared_ptr<LiteralExpr> lit_expr;
  if (lit_sum != 0) {
    lit_expr = std::make_shared<LiteralExpr>(lit_sum % mod);
  }

  std::shared_ptr<AddExpr> new_add(new AddExpr);
  new_add->literal_ = std::move(lit_expr);
  new_add->non_literal_ = std::move(non_literal_modded);
  new_add->CalculateDigitMask();
  std::shared_ptr<Expr> add_simplified = new_add->Simplify();
  if (add_simplified == nullptr) add_simplified = new_add;

  return std::make_shared<ModExpr>(std::move(add_simplified),
                                   std::move(literal_mod));
}

class EqlExpr : public Expr {
 public:
  EqlExpr(std::shared_ptr<Expr> left, std::shared_ptr<Expr> right,
          bool invert = false)
      : left_(std::move(left)), right_(std::move(right)), invert_(invert) {
    digit_mask_ = left_->DigitMask() | right_->DigitMask();
  }

  ~EqlExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    if (paren) {
      return absl::StrCat("(", left_->DebugStr(true),
                          invert_ ? " != " : " == ", right_->DebugStr(true),
                          ")");
    }
    return absl::StrCat(left_->DebugStr(true),
                        invert_ ? " != " : " == ", right_->DebugStr(true));
  }

  virtual const EqlExpr* AsEql() const override { return this; }

  std::shared_ptr<Expr> Simplify() const override {
    if (!invert_ && left_->TreeEqual(*right_)) {
      return std::make_shared<LiteralExpr>(1);
    }

    auto left_vals = left_->PossibleValues();
    auto right_vals = right_->PossibleValues();
    int64_t overlap_count = 0;
    for (const int64_t left_val : left_vals) {
      overlap_count += right_vals.contains(left_val);
    }
    if (!invert_) {
      if (overlap_count == 0) {
        return std::make_shared<LiteralExpr>(0);
      }
      if (overlap_count == left_vals.size() &&
          overlap_count == right_vals.size()) {
        return std::make_shared<LiteralExpr>(1);
      }
    } else {
      if (overlap_count == 0) {
        return std::make_shared<LiteralExpr>(1);
      } else if (overlap_count == left_vals.size() &&
                 overlap_count == right_vals.size()) {
        return std::make_shared<LiteralExpr>(0);
      }
    }

    const LiteralExpr* left_lit = left_->AsLiteral();
    const LiteralExpr* right_lit = right_->AsLiteral();
    const EqlExpr* left_eql = left_->AsEql();
    const EqlExpr* right_eql = right_->AsEql();
    if (left_lit != nullptr) {
      if (right_lit != nullptr) {
        if (invert_) {
          return std::make_shared<LiteralExpr>(left_lit->value() !=
                                               right_lit->value());
        }
        return std::make_shared<LiteralExpr>(left_lit->value() ==
                                             right_lit->value());
      }
      if (right_eql != nullptr) {
        switch (left_lit->value()) {
          case 0:
            return std::make_shared<EqlExpr>(right_eql->left_,
                                             right_eql->right_,
                                             invert_ == right_eql->invert_);
          case 1:
            return std::make_shared<EqlExpr>(right_eql->left_,
                                             right_eql->right_,
                                             invert_ != right_eql->invert_);
          default:
            return std::make_shared<LiteralExpr>(invert_ ? 1 : 0);
        }
      }
    }
    if (right_lit != nullptr && left_eql != nullptr) {
      switch (right_lit->value()) {
        case 0:
          return std::make_shared<EqlExpr>(left_eql->left_, left_eql->right_,
                                           invert_ == left_eql->invert_);
        case 1:
          return std::make_shared<EqlExpr>(left_eql->left_, left_eql->right_,
                                           invert_ != left_eql->invert_);
        default:
          return std::make_shared<LiteralExpr>(invert_ ? 1 : 0);
      }
    }
    return nullptr;
  }

  bool TreeEqual(const Expr& other) const override {
    const EqlExpr* other_eql = other.AsEql();
    if (other_eql == nullptr) return false;
    if (invert_ != other_eql->invert_) return false;
    if (left_->TreeEqual(*other_eql->left_) &&
        right_->TreeEqual(*other_eql->right_)) {
      return true;
    }
    if (left_->TreeEqual(*other_eql->right_) &&
        right_->TreeEqual(*other_eql->left_)) {
      return true;
    }
    return false;
  }

  int64_t EvalImpl(const Num& input) const override {
    return invert_ ? (left_->Eval(input) != right_->Eval(input))
                   : left_->Eval(input) == right_->Eval(input);
  }

  void CollectChildren(
      absl::flat_hash_set<std::shared_ptr<Expr>>& all_nodes) const override {
    if (all_nodes.emplace(left_).second) {
      left_->CollectChildren(all_nodes);
    }
    if (all_nodes.emplace(right_).second) {
      right_->CollectChildren(all_nodes);
    }
  }

  std::optional<std::pair<int64_t, int64_t>> ValueRange() const override {
    return std::pair<int64_t, int64_t>(0, 1);
  }

  absl::flat_hash_set<int64_t> PossibleValuesImpl() const override {
    absl::flat_hash_set<int64_t> results;
    for (int64_t left : left_->PossibleValues()) {
      for (int64_t right : right_->PossibleValues()) {
        if (left == right) {
          results.insert(invert_ ? 0 : 1);
        } else {
          results.insert(invert_ ? 1 : 0);
        }
        if (results.size() == 2) return results;
      }
    }
    return results;
  }

 private:
  std::shared_ptr<Expr> left_;
  std::shared_ptr<Expr> right_;
  bool invert_ = false;
};

class MachineState {
 public:
  MachineState() = default;

  const Expr& w() const { return *w_; }
  const Expr& x() const { return *x_; }
  const Expr& y() const { return *y_; }
  const Expr& z() const { return *z_; }

  std::shared_ptr<Expr> z_shared() const { return z_; }

  void RunInstruction(absl::string_view inst) {
    const absl::string_view mnemonic = inst.substr(0, 3);
    std::shared_ptr<Expr>& reg_0 = Register(inst[4]);

    if (mnemonic == "inp") {
      auto input_expr = std::make_shared<InputExpr>(inputs_.size());
      inputs_.emplace_back(input_expr);
      reg_0 = std::move(input_expr);
      return;
    }

    std::shared_ptr<Expr> arg_1;
    int64_t arg_1_int = 0;
    if (absl::SimpleAtoi(inst.substr(5), &arg_1_int)) {
      arg_1 = std::make_shared<LiteralExpr>(arg_1_int);
    } else {
      arg_1 = Register(inst.back());
    }

    if (mnemonic == "add") {
      reg_0 = std::make_shared<AddExpr>(std::move(reg_0), std::move(arg_1));
      TrySimplify(reg_0);
      return;
    }

    if (mnemonic == "mul") {
      reg_0 = std::make_shared<MulExpr>(std::move(reg_0), std::move(arg_1));
      TrySimplify(reg_0);
      return;
    }

    if (mnemonic == "div") {
      reg_0 = std::make_shared<DivExpr>(std::move(reg_0), std::move(arg_1));
      TrySimplify(reg_0);
      return;
    }

    if (mnemonic == "mod") {
      reg_0 = std::make_shared<ModExpr>(std::move(reg_0), std::move(arg_1));
      TrySimplify(reg_0);
      return;
    }

    if (mnemonic == "eql") {
      reg_0 = std::make_shared<EqlExpr>(std::move(reg_0), std::move(arg_1));
      TrySimplify(reg_0);
      return;
    }

    CHECK_FAIL();
  }

  void RunProgram(const std::vector<std::string>& program) {
    for (int pc = 0; pc < program.size(); ++pc) {
      RunInstruction(program[pc]);
    }
    std::cout << "Program evaluated\n";
  }

  int64_t SearchProgram() {
    absl::flat_hash_set<std::shared_ptr<Expr>> collected{z_};
    z_->CollectChildren(collected);

    std::map<int, int> deps_distro;
    for (const auto& expr : collected) {
      ++deps_distro[std::popcount(expr->DigitMask())];
    }
    for (const auto [deps, count] : deps_distro) {
      std::cout << deps << " dependencies -> " << count << " nodes\n";
    }

    std::vector<std::vector<std::shared_ptr<Expr>>> input_dependents(
        inputs_.size());
    for (int i = 0; i < inputs_.size(); ++i) {
      for (const auto& expr : collected) {
        if ((expr->DigitMask() & (1 << i)) != 0) {
          input_dependents[i].emplace_back(expr);
        }
      }
      std::cout << input_dependents[i].size() << " nodes depend on input " << i
                << "\n";
    }
    CHECK(SearchProgramImpl(input_dependents));

    int64_t result = 0;
    for (const auto& input : inputs_) {
      result = result * 10 + input->FinalizeTentative();
    }
    return result;
  }

 private:
  bool SearchProgramImpl(
      const std::vector<std::vector<std::shared_ptr<Expr>>>& input_dependents,
      int pos = 0, int64_t selected = 0) {
    CHECK(pos < inputs_.size());
    for (int tentative = 1; tentative <= 9; ++tentative) {
      for (const auto& dependent : input_dependents[pos]) {
        dependent->PurgePossibleValuesMemo();
      }
      inputs_[pos]->SetTentativeValue(tentative);
      const auto& possible = z_->PossibleValues();
      if (pos < 5) {
        for (int space = 0; space < pos; ++space) {
          std::cout << "  ";
        }
        std::cout << "input_" << pos << " = " << tentative << " -> "
                  << possible.size() << "\n";
      }
      if (!possible.contains(0)) continue;
      if (possible.size() == 1) return true;
      if (SearchProgramImpl(input_dependents, pos + 1)) return true;
    }
    for (const auto& dependent : input_dependents[pos]) {
      dependent->PurgePossibleValuesMemo();
    }
    inputs_[pos]->SetTentativeValue(std::nullopt);
    return false;
  }

  std::shared_ptr<Expr>& Register(char c) {
    switch (c) {
      case 'w':
        return w_;
      case 'x':
        return x_;
      case 'y':
        return y_;
      case 'z':
        return z_;
      default:
        CHECK_FAIL();
    }
  }

  std::shared_ptr<Expr> w_ = std::make_shared<LiteralExpr>(0);
  std::shared_ptr<Expr> x_ = std::make_shared<LiteralExpr>(0);
  std::shared_ptr<Expr> y_ = std::make_shared<LiteralExpr>(0);
  std::shared_ptr<Expr> z_ = std::make_shared<LiteralExpr>(0);
  std::vector<std::shared_ptr<InputExpr>> inputs_;
};

Num FindLargestModelNum(const std::vector<std::string>& program) {
  MachineState machine;
  machine.RunProgram(program);

  std::cout << "Program evaluated\n";

  Num input = kLargest;
  int64_t counter = 0;
  while (machine.z().Eval(input) != 0) {
    if ((++counter & 0xFFFF) == 0) {
      std::cout << "Trying " << input << "\n";
    }
    --input;
  }
  return input;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> program = aoc2021::ReadLinesFromFile(argv[1]);

  MachineState machine;
  machine.RunProgram(program);
  std::cout << machine.SearchProgram() << "\n";
  //  std::cout << "w = " << machine.w().DebugStr() << "\n\n";
  //  std::cout << "x = " << machine.x().DebugStr() << "\n\n";
  //  std::cout << "y = " << machine.y().DebugStr() << "\n\n";
  //  std::cout << "z = " << machine.z().DebugStr() << "\n\n";

  //  std::cout << FindLargestModelNum(program) << "\n";
  return 0;
}
