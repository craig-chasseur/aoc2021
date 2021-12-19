#include <deque>
#include <iterator>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "util/check.h"
#include "util/dimension_grid.h"
#include "util/io.h"

namespace {

struct SensorData {
  SensorData() = default;

  SensorData(const SensorData& other)
      : position(other.position), beacon_signals(other.beacon_signals) {}
  SensorData(SensorData&& other)
      : position(std::move(other.position)),
        beacon_signals(std::move(other.beacon_signals)),
        all_orientations_(std::move(other.all_orientations_)) {}

  SensorData& operator=(const SensorData& other) {
    if (&other == this) return *this;
    position = other.position;
    beacon_signals = other.beacon_signals;
    all_orientations_.reset();
    return *this;
  }

  SensorData& operator=(SensorData&& other) {
    if (&other == this) return *this;
    position = std::move(other.position);
    beacon_signals = std::move(other.beacon_signals);
    all_orientations_ = std::move(other.all_orientations_);
    return *this;
  }

  std::optional<aoc2021::DimensionGrid<3>::Point> position;
  std::vector<aoc2021::DimensionGrid<3>::Vec> beacon_signals;

  const std::vector<SensorData>& AllOrientations() const {
    if (all_orientations_ != nullptr) return *all_orientations_;
    all_orientations_ = std::make_unique<std::vector<SensorData>>();
    for (const auto& rot :
         aoc2021::DimensionGrid<3>::Rotations::AllOrientations()) {
      SensorData oriented;
      oriented.position = position;
      oriented.beacon_signals = beacon_signals * rot;
      all_orientations_->emplace_back(std::move(oriented));
    }
    CHECK(all_orientations_->size() == 24);
    return *all_orientations_;
  }

 private:
  mutable std::unique_ptr<std::vector<SensorData>> all_orientations_;
};

std::optional<SensorData> MatchBeacons(const SensorData& known,
                                       const SensorData& candidate) {
  for (const auto& candidate_orientation : candidate.AllOrientations()) {
    absl::flat_hash_map<aoc2021::DimensionGrid<3>::Vec, int> offset_counts;
    for (const auto& outer_signal : known.beacon_signals) {
      for (const auto& inner_signal : candidate_orientation.beacon_signals) {
        ++offset_counts[outer_signal - inner_signal];
      }
    }
    for (auto& [offset, count] : offset_counts) {
      if (count >= 12) {
        SensorData solved_orientation(candidate_orientation);
        solved_orientation.position = known.position.value() + offset;
        return solved_orientation;
      }
    }
  }
  return std::nullopt;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<std::vector<std::string>> scanner_lines =
      aoc2021::SplitByEmptyStrings(std::move(lines));
  std::list<SensorData> sensor_data;
  for (const std::vector<std::string>& scanner_desc : scanner_lines) {
    SensorData sensor;
    for (auto line_it = scanner_desc.begin() + 1; line_it != scanner_desc.end();
         ++line_it) {
      sensor.beacon_signals.emplace_back(
          aoc2021::DimensionGrid<3>::Point::ParseCommaSeparated(*line_it) -
          aoc2021::DimensionGrid<3>::Points::kOrigin);
    }
    sensor_data.emplace_back(std::move(sensor));
  }

  std::vector<SensorData> exhausted;
  std::deque<SensorData> resolved;
  resolved.emplace_back(std::move(sensor_data.front()));
  sensor_data.erase(sensor_data.begin());
  resolved.front().position = aoc2021::DimensionGrid<3>::Points::kOrigin;
  while (!sensor_data.empty()) {
    std::cout << "Resolved " << (resolved.size() + exhausted.size())
              << " sensors (" << exhausted.size() << " exhausted)\n";
    SensorData current_resolved = std::move(resolved.front());
    resolved.pop_front();
    std::vector<std::list<SensorData>::iterator> newly_resolved;
    for (auto unresolved_it = sensor_data.begin();
         unresolved_it != sensor_data.end(); ++unresolved_it) {
      if (auto maybe_resolved = MatchBeacons(current_resolved, *unresolved_it);
          maybe_resolved.has_value()) {
        resolved.emplace_back(*maybe_resolved);
        newly_resolved.emplace_back(unresolved_it);
      }
    }
    exhausted.emplace_back(std::move(current_resolved));
    for (const auto to_remove : newly_resolved) {
      sensor_data.erase(to_remove);
    }
  }
  std::cout << "Resolved " << (resolved.size() + exhausted.size())
            << " sensors (" << exhausted.size() << " exhausted)\n";

  exhausted.insert(exhausted.end(), std::make_move_iterator(resolved.begin()),
                   std::make_move_iterator(resolved.end()));
  absl::flat_hash_set<aoc2021::DimensionGrid<3>::Point> beacons;
  for (const SensorData& sensor : exhausted) {
    for (const aoc2021::DimensionGrid<3>::Vec& signal : sensor.beacon_signals) {
      beacons.emplace(sensor.position.value() + signal);
    }
  }
  std::cout << "Part 1 (total beacons): " << beacons.size() << "\n";

  std::int64_t max_manhattan = -1;
  for (const SensorData& outer : exhausted) {
    for (const SensorData& inner : exhausted) {
      max_manhattan = std::max(max_manhattan,
                               (outer.position.value() - inner.position.value())
                                   .ManhattanDistance());
    }
  }
  std::cout << "Part 2 (max manhattan distance): " << max_manhattan << "\n";

  return 0;
}
