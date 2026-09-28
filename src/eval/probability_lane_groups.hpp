#pragma once

#include <cstdint>
#include <cstring>

#include "eval_query.hpp"
#include "exact_plan_types.hpp"
#include "../leaf/dist_kind.hpp"

namespace accumulatr::eval::detail {

// A call has one variant and probability measure. Only parameters and interval
// endpoints vary between requests; component weights are applied afterwards.
// Grouping is rebuilt each call; numerical results never cross particles.
struct ProbabilityLaneGroups {
  ObservationLaneBatchView prepare(const ExactVariantPlan &plan,
                                   const ObservationLaneBatchView input,
                                   const double *lower = nullptr,
                                   const double *upper = nullptr) {
    grouped = false;
    if (input.size < 2U) return input;
    std::size_t width = lower == nullptr ? 0U : 2U;
    for (const auto &leaf : plan.leaf_descriptors)
      width += 3U + leaf::dist_param_count(static_cast<leaf::DistKind>(leaf.dist_kind));
    keys.resize(input.size * width);
    hashes.resize(input.size);
    destinations.resize(input.size);
    representatives.clear();
    std::size_t slots = 2U;
    while (slots < 2U * input.size) slots *= 2U;
    table.assign(slots, input.size);

    for (std::size_t lane = 0; lane < input.size; ++lane) {
      auto &hash = hashes[lane];
      hash = 14695981039346656037ULL;
      std::size_t field = lane * width;
      const auto append = [&](const double value) {
        std::uint64_t bits = 0;
        if (value != 0.0) std::memcpy(&bits, &value, sizeof(bits));
        keys[field++] = bits;
        hash = (hash ^ bits) * 1099511628211ULL;
      };
      for (std::size_t leaf = 0; leaf < plan.leaf_descriptors.size(); ++leaf) {
        const auto &descriptor = plan.leaf_descriptors[leaf];
        const auto row = input.physical_row(leaf, lane);
        const auto columns = 2U + leaf::dist_param_count(
            static_cast<leaf::DistKind>(descriptor.dist_kind));
        for (std::size_t column = 0; column < columns; ++column)
          append(input.matrix->base[column * input.matrix->nrow + row]);
        append(input.matrix->onset == nullptr
                   ? descriptor.onset_abs_value : input.matrix->onset[row]);
      }
      if (lower != nullptr) {
        append(lower[lane]);
        append(upper[lane]);
      }
      auto slot = static_cast<std::size_t>(hash) & (slots - 1U);
      while (table[slot] != input.size) {
        const auto previous = table[slot];
        if (hashes[previous] == hash &&
            std::memcmp(keys.data() + previous * width, keys.data() + lane * width,
                        width * sizeof(std::uint64_t)) == 0) break;
        slot = (slot + 1U) & (slots - 1U);
      }
      if (table[slot] == input.size) {
        table[slot] = lane;
        destinations[lane] = representatives.size();
        representatives.push_back(lane);
      } else {
        destinations[lane] = destinations[table[slot]];
      }
    }
    grouped = representatives.size() != input.size;
    if (!grouped) return input;
    lanes.clear();
    lower_bounds.clear();
    upper_bounds.clear();
    for (const auto lane : representatives) {
      lanes.emplace_back(input.row_maps[lane], input.row_offsets[lane],
                         input.observed_times[lane]);
      if (lower != nullptr) {
        lower_bounds.push_back(lower[lane]);
        upper_bounds.push_back(upper[lane]);
      }
    }
    return lanes.view(*input.matrix);
  }

  void expand(std::vector<double> *values) {
    if (!grouped) return;
    expanded.resize(destinations.size());
    for (std::size_t lane = 0; lane < destinations.size(); ++lane)
      expanded[lane] = (*values)[destinations[lane]];
    values->swap(expanded);
  }

  bool grouped{false};
  ObservationLaneBatch lanes;
  std::vector<std::uint64_t> keys, hashes;
  std::vector<std::size_t> table, destinations, representatives;
  std::vector<double> lower_bounds, upper_bounds, expanded;
};

} // namespace accumulatr::eval::detail
