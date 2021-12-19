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

  explicit SensorData(
      std::vector<aoc2021::DimensionGrid<3>::Vec> beacon_signals)
      : beacon_signals_(std::move(beacon_signals)) {}

  SensorData(const SensorData& other)
      : position_(other.position_), beacon_signals_(other.beacon_signals_) {}
  SensorData(SensorData&& other)
      : position_(std::move(other.position_)),
        beacon_signals_(std::move(other.beacon_signals_)),
        all_orientations_(std::move(other.all_orientations_)) {}

  SensorData& operator=(const SensorData& other) {
    if (&other == this) return *this;
    position_ = other.position_;
    beacon_signals_ = other.beacon_signals_;
    all_orientations_.reset();
    return *this;
  }

  SensorData& operator=(SensorData&& other) {
    if (&other == this) return *this;
    position_ = std::move(other.position_);
    beacon_signals_ = std::move(other.beacon_signals_);
    all_orientations_ = std::move(other.all_orientations_);
    return *this;
  }

  void SetPosition(aoc2021::DimensionGrid<3>::Point position) {
    position_ = std::move(position);
  }

  bool MatchToKnown(const SensorData& known) {
    CHECK(known.position_.has_value());
    CHECK(!position_.has_value());
    for (const auto& candidate_orientation : AllOrientations()) {
      absl::flat_hash_map<aoc2021::DimensionGrid<3>::Vec, int> offset_counts;
      for (const auto& outer_signal : known.beacon_signals_) {
        for (const auto& inner_signal : candidate_orientation) {
          ++offset_counts[outer_signal - inner_signal];
        }
      }
      for (auto& [offset, count] : offset_counts) {
        if (count >= 12) {
          position_ = known.position_.value() + offset;
          beacon_signals_ = candidate_orientation;
          return true;
        }
      }
    }
    return false;
  }

  aoc2021::DimensionGrid<3>::Point position() const {
    CHECK(position_.has_value());
    return *position_;
  }

  std::vector<aoc2021::DimensionGrid<3>::Point> AbsoluteBeaconPositions()
      const {
    CHECK(position_.has_value());
    return *position_ + beacon_signals_;
  }

 private:
  const std::vector<std::vector<aoc2021::DimensionGrid<3>::Vec>>&
  AllOrientations() {
    if (all_orientations_ != nullptr) return *all_orientations_;
    all_orientations_ = std::make_unique<
        std::vector<std::vector<aoc2021::DimensionGrid<3>::Vec>>>();
    for (const auto& rot :
         aoc2021::DimensionGrid<3>::Rotations::AllOrientations()) {
      all_orientations_->emplace_back(beacon_signals_ * rot);
    }
    CHECK(all_orientations_->size() == 24);
    return *all_orientations_;
  }

  std::optional<aoc2021::DimensionGrid<3>::Point> position_;
  std::vector<aoc2021::DimensionGrid<3>::Vec> beacon_signals_;
  std::unique_ptr<std::vector<std::vector<aoc2021::DimensionGrid<3>::Vec>>>
      all_orientations_;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> lines = aoc2021::ReadLinesFromFile(argv[1]);
  std::vector<std::vector<std::string>> scanner_lines =
      aoc2021::SplitByEmptyStrings(std::move(lines));
  std::list<SensorData> sensor_data;
  for (const std::vector<std::string>& scanner_desc : scanner_lines) {
    std::vector<aoc2021::DimensionGrid<3>::Vec> beacon_signals;
    for (auto line_it = scanner_desc.begin() + 1; line_it != scanner_desc.end();
         ++line_it) {
      beacon_signals.emplace_back(
          aoc2021::DimensionGrid<3>::Point::ParseCommaSeparated(*line_it) -
          aoc2021::DimensionGrid<3>::Points::kOrigin);
    }
    sensor_data.emplace_back(std::move(beacon_signals));
  }

  std::vector<SensorData> exhausted;
  std::deque<SensorData> resolved;
  resolved.emplace_back(std::move(sensor_data.front()));
  sensor_data.erase(sensor_data.begin());
  resolved.front().SetPosition(aoc2021::DimensionGrid<3>::Points::kOrigin);
  while (!sensor_data.empty()) {
    std::cout << "Resolved " << (resolved.size() + exhausted.size())
              << " sensors (" << exhausted.size() << " exhausted)\n";
    SensorData current_resolved = std::move(resolved.front());
    resolved.pop_front();
    std::vector<std::list<SensorData>::iterator> newly_resolved;
    for (auto unresolved_it = sensor_data.begin();
         unresolved_it != sensor_data.end(); ++unresolved_it) {
      if (unresolved_it->MatchToKnown(current_resolved)) {
        resolved.emplace_back(std::move(*unresolved_it));
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
    auto beacons_from_sensor = sensor.AbsoluteBeaconPositions();
    beacons.insert(beacons_from_sensor.begin(), beacons_from_sensor.end());
  }
  std::cout << "Part 1 (total beacons): " << beacons.size() << "\n";

  std::int64_t max_manhattan = -1;
  for (const SensorData& outer : exhausted) {
    for (const SensorData& inner : exhausted) {
      max_manhattan =
          std::max(max_manhattan,
                   (outer.position() - inner.position()).ManhattanDistance());
    }
  }
  std::cout << "Part 2 (max manhattan distance): " << max_manhattan << "\n";

  return 0;
}
