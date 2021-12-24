#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
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

 protected:
  mutable absl::flat_hash_map<Num, int64_t> masked_input_to_value_;

  std::uint16_t digit_mask_ = 0;
};

void TrySimplify(std::shared_ptr<Expr>& expr) {
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

 private:
  int input_idx_ = 0;
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

 private:
  AddExpr() = default;

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
    if (literal_ != nullptr && non_literal_.empty()) {
      return literal_;
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

 private:
  MulExpr() = default;

  std::shared_ptr<LiteralExpr> literal_;  // May be null.
  std::vector<std::shared_ptr<Expr>> non_literal_;
};

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
    }

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

 private:
  std::shared_ptr<Expr> left_;
  std::shared_ptr<Expr> right_;
};

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

 private:
  std::shared_ptr<Expr> left_;
  std::shared_ptr<Expr> right_;
};

class EqlExpr : public Expr {
 public:
  EqlExpr(std::shared_ptr<Expr> left, std::shared_ptr<Expr> right)
      : left_(std::move(left)), right_(std::move(right)) {
    digit_mask_ = left_->DigitMask() | right_->DigitMask();
  }

  ~EqlExpr() override = default;

  std::string DebugStr(bool paren = false) const override {
    if (paren) {
      return absl::StrCat("(", left_->DebugStr(true),
                          " == ", right_->DebugStr(true), ")");
    }
    return absl::StrCat(left_->DebugStr(true), " == ", right_->DebugStr(true));
  }

  virtual const EqlExpr* AsEql() const override { return this; }

  std::shared_ptr<Expr> Simplify() const override {
    const LiteralExpr* left_lit = left_->AsLiteral();
    const LiteralExpr* right_lit = right_->AsLiteral();
    if (left_lit != nullptr && right_lit != nullptr) {
      return std::make_shared<LiteralExpr>(left_lit->value() ==
                                           right_lit->value());
    }
    if (left_->TreeEqual(*right_)) {
      return std::make_shared<LiteralExpr>(1);
    }
    return nullptr;
  }

  bool TreeEqual(const Expr& other) const override {
    const EqlExpr* other_eql = other.AsEql();
    if (other_eql == nullptr) return false;
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
    return left_->Eval(input) == right_->Eval(input);
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

 private:
  std::shared_ptr<Expr> left_;
  std::shared_ptr<Expr> right_;
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
      reg_0 = std::make_shared<InputExpr>(input_counter_++);
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
      std::cout << "PC=" << pc << " : " << program[pc] << "\n";
      RunInstruction(program[pc]);
    }
  }

 private:
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
  int input_counter_ = 0;
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
  std::cout << "w = " << machine.w().DebugStr() << "\n\n";
  std::cout << "x = " << machine.x().DebugStr() << "\n\n";
  std::cout << "y = " << machine.y().DebugStr() << "\n\n";
  //  std::cout << "z = " << machine.z().DebugStr() << "\n\n";

//  std::cout << FindLargestModelNum(program) << "\n";
  return 0;
}
