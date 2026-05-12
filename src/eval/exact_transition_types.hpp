#pragma once

#include <cstdint>
#include <vector>

#include "exact_common_types.hpp"

namespace accumulatr::eval {
namespace detail {

enum class ExactSymbolicTransitionTimeKind : std::uint8_t {
  SourceRelease = 0
};

struct ExactSymbolicTransitionTimeExpr {
  ExactSymbolicTransitionTimeKind kind{
      ExactSymbolicTransitionTimeKind::SourceRelease};
  semantic::Index source_id{semantic::kInvalidIndex};
};

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
  double empty_value{1.0};

  [[nodiscard]] bool empty() const noexcept {
    return guards.empty();
  }
};

struct ExactSymbolicReadinessTimeExpr {
  ExactTransitionGuardSet requirements;

  [[nodiscard]] bool present() const noexcept {
    return !requirements.empty();
  }
};

struct ExactSymbolicTransitionOrderRegion {
  std::vector<ExactSourceOrderFact> source_order_facts;
};

struct ExactSymbolicTransitionTime {
  ExactSymbolicTransitionTimeExpr transition_time_expr;
  ExactSymbolicReadinessTimeExpr readiness_time_expr;
  std::vector<semantic::Index> active_sources;
  ExactTransitionGuardSet guards;
  ExactSymbolicTransitionOrderRegion order_region;
  semantic::Index source_view_id{0};
  ExactRelationTemplate relation_template;
};

inline semantic::Index exact_symbolic_transition_release_source_id(
    const ExactSymbolicTransitionTime &transition) noexcept {
  return transition.transition_time_expr.kind ==
                 ExactSymbolicTransitionTimeKind::SourceRelease
             ? transition.transition_time_expr.source_id
             : semantic::kInvalidIndex;
}

inline bool exact_symbolic_transition_has_readiness(
    const ExactSymbolicTransitionTime &transition) noexcept {
  return transition.readiness_time_expr.present();
}

struct ExactSymbolicTransitionRelation {
  bool competitor_can_strictly_precede{false};
  bool competitor_can_positively_coincide{false};
};

inline ExactSymbolicTransitionRelation exact_symbolic_transition_relation(
    const ExactSymbolicTransitionTime &target,
    const ExactSymbolicTransitionTime &competitor) {
  const auto target_source_id =
      exact_symbolic_transition_release_source_id(target);
  const auto competitor_source_id =
      exact_symbolic_transition_release_source_id(competitor);
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
  semantic::Index visible_outcome{semantic::kInvalidIndex};
  ExactSymbolicTransitionTime transition;
  semantic::Index probability_root_id{semantic::kInvalidIndex};
};

} // namespace detail
} // namespace accumulatr::eval
