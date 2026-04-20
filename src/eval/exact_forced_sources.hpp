#pragma once

#include "exact_oracle.hpp"

namespace accumulatr::eval {
namespace detail {

inline const std::vector<semantic::Index> &source_support_for(
    const ExactVariantPlan &plan,
    const ExactSourceKey key) {
  if (key.kind == semantic::SourceKind::Leaf) {
    return plan.leaf_supports[static_cast<std::size_t>(key.index)];
  }
  return plan.pool_supports[static_cast<std::size_t>(key.index)];
}

inline void validate_forced_source_compatibility(
    const ExactVariantPlan &plan,
    const std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> &forced,
    const ExactSourceKey key) {
  const auto &support = source_support_for(plan, key);
  for (const auto &[forced_key, relation] : forced) {
    (void) relation;
    if (forced_key == key) {
      continue;
    }
    if (key.kind == semantic::SourceKind::Leaf &&
        forced_key.kind == semantic::SourceKind::Pool &&
        supports_overlap(support, source_support_for(plan, forced_key))) {
      throw std::runtime_error(
          "exact kernel cannot condition a leaf on an overlapping forced pool state yet");
    }
    if (key.kind == semantic::SourceKind::Pool &&
        forced_key.kind == semantic::SourceKind::Pool &&
        supports_overlap(support, source_support_for(plan, forced_key))) {
      throw std::runtime_error(
          "exact kernel cannot condition overlapping pool states yet");
    }
  }
}

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
  validate_forced_source_compatibility(plan, forced, key);
  return oracle->source_channels(key.kind, key.index, oracle_time);
}

} // namespace detail
} // namespace accumulatr::eval
