#pragma once

#include <cstdint>
#include <vector>

#include "exact_common_types.hpp"

namespace accumulatr::eval {
namespace detail {

enum class ExactTransitionGuardKind : std::uint8_t {
  SourceBefore = 0,
  SourceAfter = 1,
  ExprBefore = 2,
  ExprAfter = 3
};

struct ExactTransitionGuard {
  ExactTransitionGuardKind kind{ExactTransitionGuardKind::SourceBefore};
  semantic::Index subject_id{semantic::kInvalidIndex};
};

struct ExactTransitionGuardSet {
  std::vector<ExactTransitionGuard> guards;
  [[nodiscard]] bool empty() const noexcept {
    return guards.empty();
  }
};

struct ExactSymbolicTransitionTime {
  semantic::Index release_source_id{semantic::kInvalidIndex};
  ExactTransitionGuardSet readiness;
  ExactTransitionGuardSet guards;
  std::vector<ExactSourceOrderFact> source_order_facts;
  semantic::Index source_view_id{0};
  ExactRelationTemplate relation_template;
};

struct ExactSymbolicTransitionRelation {
  bool competitor_can_strictly_precede{false};
  bool competitor_can_positively_coincide{false};
};

inline ExactSymbolicTransitionRelation exact_symbolic_transition_relation(
    const ExactSymbolicTransitionTime &target,
    const ExactSymbolicTransitionTime &competitor) {
  const auto target_source_id = target.release_source_id;
  const auto competitor_source_id = competitor.release_source_id;
  if (target_source_id == semantic::kInvalidIndex ||
      competitor_source_id == semantic::kInvalidIndex) {
    return {};
  }
  if (target_source_id == competitor_source_id) {
    return ExactSymbolicTransitionRelation{false, true};
  }
  return ExactSymbolicTransitionRelation{true, false};
}

struct ExactSymbolicTransitionScenario {
  ExactSymbolicTransitionTime transition;
  semantic::Index probability_root_id{semantic::kInvalidIndex};
};

} // namespace detail
} // namespace accumulatr::eval
