#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <variant>
#include <vector>

#include "util/check.h"
#include "util/io.h"

namespace {

std::vector<bool> Binify(const std::string& hex_string) {
  std::vector<bool> binary;
  for (char c : hex_string) {
    uint8_t val = 0;
    if (c >= '0' && c <= '9') {
      val = c - '0';
    } else if (c >= 'A' && c <= 'F') {
      val = c - 'A' + 10;
    } else if (c == '\n') {
      continue;
    } else {
      CHECK_FAIL();
    }
    CHECK(val < 16);
    binary.push_back((val & 8) != 0);
    binary.push_back((val & 4) != 0);
    binary.push_back((val & 2) != 0);
    binary.push_back((val & 1) != 0);
  }
  return binary;
}

struct Packet {
  uint8_t version_ = 0;
  uint8_t type_ = 0;

  std::variant<uint64_t, std::vector<std::unique_ptr<Packet>>> contents_;

  uint64_t TreeVersion() const {
    uint64_t sum = version_;
    if (std::holds_alternative<std::vector<std::unique_ptr<Packet>>>(
            contents_)) {
      for (const std::unique_ptr<Packet>& subpacket :
           std::get<std::vector<std::unique_ptr<Packet>>>(contents_)) {
        sum += subpacket->TreeVersion();
      }
    }
    return sum;
  }

  uint64_t Value() const {
    switch (type_) {
      case 0: {
        uint64_t sum = 0;
        for (const auto& subpacket : std::get<1>(contents_)) {
          sum += subpacket->Value();
        }
        return sum;
      }
      case 1: {
        uint64_t product = 1;
        for (const auto& subpacket : std::get<1>(contents_)) {
          product *= subpacket->Value();
        }
        return product;
      }
      case 2: {
        uint64_t min = std::numeric_limits<uint64_t>::max();
        for (const auto& subpacket : std::get<1>(contents_)) {
          min = std::min(min, subpacket->Value());
        }
        return min;
      }
      case 3: {
        uint64_t max = std::numeric_limits<uint64_t>::min();
        for (const auto& subpacket : std::get<1>(contents_)) {
          max = std::max(max, subpacket->Value());
        }
        return max;
      }
      case 4: {
        return std::get<0>(contents_);
      }
      case 5: {
        CHECK(std::get<1>(contents_).size() == 2);
        return std::get<1>(contents_).front()->Value() >
                       std::get<1>(contents_).back()->Value()
                   ? 1
                   : 0;
      }
      case 6: {
        CHECK(std::get<1>(contents_).size() == 2);
        return std::get<1>(contents_).front()->Value() <
                       std::get<1>(contents_).back()->Value()
                   ? 1
                   : 0;
      }
      case 7: {
        CHECK(std::get<1>(contents_).size() == 2);
        return std::get<1>(contents_).front()->Value() ==
                       std::get<1>(contents_).back()->Value()
                   ? 1
                   : 0;
      }
      default:
        CHECK_FAIL();
    }
  }
};

uint64_t ParseBits(std::vector<bool>::const_iterator begin,
                   std::vector<bool>::const_iterator end) {
  uint64_t result = 0;
  for (; begin != end; ++begin) {
    result <<= 1;
    if (*begin) result |= 1;
  }
  return result;
}

std::pair<uint64_t, std::vector<bool>::const_iterator> ParseLit(
    std::vector<bool>::const_iterator begin,
    std::vector<bool>::const_iterator end) {
  uint64_t value = 0;
  for (;;) {
    CHECK(end - begin >= 5);
    const bool terminate = !*(begin++);
    value = (value << 4) | ParseBits(begin, begin + 4);
    begin += 4;
    if (terminate) return {value, begin};
  }
}

struct PacketAndIterator {
  std::unique_ptr<Packet> packet;
  std::vector<bool>::const_iterator after;
};

PacketAndIterator Parse(std::vector<bool>::const_iterator begin,
                        std::vector<bool>::const_iterator end) {
  PacketAndIterator result;
  result.packet = std::make_unique<Packet>();

  CHECK(end - begin >= 3);
  result.packet->version_ = ParseBits(begin, begin + 3);
  begin += 3;

  CHECK(end - begin >= 3);
  result.packet->type_ = ParseBits(begin, begin + 3);
  begin += 3;

  if (result.packet->type_ == 4) {
    auto [literal, after_literal] = ParseLit(begin, end);
    result.packet->contents_.emplace<uint64_t>(literal);
    begin = after_literal;
  } else {
    std::vector<std::unique_ptr<Packet>> subpackets;
    CHECK(begin != end);
    const bool length_as_subs = *(begin++);
    if (!length_as_subs) {
      CHECK(end - begin >= 15);
      const uint64_t subpacket_bits = ParseBits(begin, begin + 15);
      begin += 15;
      CHECK(end - begin >= subpacket_bits);
      const std::vector<bool>::const_iterator end_of_subpackets =
          begin + subpacket_bits;
      while (begin != end_of_subpackets) {
        PacketAndIterator subpacket_result = Parse(begin, end_of_subpackets);
        subpackets.emplace_back(std::move(subpacket_result.packet));
        begin = subpacket_result.after;
      }
    } else {
      CHECK(end - begin >= 11);
      const uint64_t num_subpackets = ParseBits(begin, begin + 11);
      begin += 11;
      for (uint64_t subpacket_num = 0; subpacket_num < num_subpackets;
           ++subpacket_num) {
        PacketAndIterator subpacket_result = Parse(begin, end);
        subpackets.emplace_back(std::move(subpacket_result.packet));
        begin = subpacket_result.after;
      }
    }
    result.packet->contents_.emplace<std::vector<std::unique_ptr<Packet>>>(
        std::move(subpackets));
  }

  result.after = begin;
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  std::string hex_string = aoc2021::ReadFile(argv[1]);
  std::vector<bool> binary = Binify(hex_string);
  PacketAndIterator result = Parse(binary.begin(), binary.end());
  std::cout << "Part 1: " << result.packet->TreeVersion() << "\n";
  std::cout << "Part 2: " << result.packet->Value() << "\n";
  return 0;
}
