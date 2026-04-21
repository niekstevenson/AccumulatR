#pragma once

#include "exact_oracle.hpp"

namespace accumulatr::eval {
namespace detail {

inline leaf::EventChannels resolve_forced_source_channels(
    const ExactVariantPlan &plan,
    ExactSourceOracle *oracle,
    const std::vector<ExactRelation> *forced,
    const double oracle_time,
    const ExactSourceKey key) {
  if (forced != nullptr) {
    const auto source_id = source_ordinal(plan, key);
    if (source_id != semantic::kInvalidIndex) {
      const auto relation = (*forced)[static_cast<std::size_t>(source_id)];
      if (relation != ExactRelation::Unknown) {
        return forced_channels(relation);
      }
    }
  }
  return oracle->source_channels(key.kind, key.index, oracle_time);
}

} // namespace detail
} // namespace accumulatr::eval
