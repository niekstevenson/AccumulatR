#pragma once

#include "exact_oracle.hpp"

namespace accumulatr::eval {
namespace detail {

inline leaf::EventChannels resolve_forced_source_channels(
    const ExactVariantPlan &plan,
    ExactSourceOracle *oracle,
    const std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> &forced,
    const double oracle_time,
    const ExactSourceKey key) {
  const auto it = forced.find(key);
  if (it != forced.end()) {
    return forced_channels(it->second);
  }
  (void) plan;
  return oracle->source_channels(key.kind, key.index, oracle_time);
}

} // namespace detail
} // namespace accumulatr::eval
