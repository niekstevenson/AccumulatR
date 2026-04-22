#pragma once

#include "exact_oracle.hpp"

namespace accumulatr::eval {
namespace detail {

inline leaf::EventChannels resolve_forced_source_channels(
    ExactSourceOracle *oracle,
    const RelationView &relations,
    const semantic::Index source_id) {
  const auto relation = relations.relation_for(source_id);
  if (relation != ExactRelation::Unknown) {
    return forced_channels(relation);
  }
  return oracle->conditional_source(source_id);
}

} // namespace detail
} // namespace accumulatr::eval
