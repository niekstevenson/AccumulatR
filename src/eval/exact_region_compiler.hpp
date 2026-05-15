#pragma once

#include "exact_transition_lowering.hpp"
#include "exact_types.hpp"
#include "exact_region.hpp"
#include "exact_compiled_math_lowering.hpp"
#include "exact_projection_planner.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool exact_order_region_expr_before(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_not_before(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_after(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_before_or_at(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_at_time(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_satisfied_at_time(
    const ExactVariantBuildState &plan,
    semantic::Index expr_id,
    semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_not_satisfied_at_time(
    const ExactVariantBuildState &plan,
    semantic::Index expr_id,
    semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_relation_can_collapse(
    const ExactVariantBuildState &plan,
    semantic::Index expr_id);

inline bool exact_expr_completion_monotone(
    const ExactVariantBuildState &plan,
    semantic::Index expr_id);

inline bool exact_order_region_expand_relation_factor(
    const ExactVariantBuildState &plan,
    const ExactOrderRegionExprValueFactor &factor,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline ExactOrderRegionExpr exact_order_region_subtract_region(
    ExactOrderRegionExpr domain,
    const ExactOrderRegionExpr &covered);

inline bool exact_order_region_truth_at_time(
    const ExactVariantBuildState &plan,
    const ExactTransitionGuardSet &guards,
    semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_guard_before_or_at(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_transition_requirements_fail_at_times(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &transition,
    const ExactTransitionGuardSet &guards,
    semantic::Index readiness_time_id,
    semantic::Index guard_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_apply_source_order_facts(
    ExactOrderRegionExpr *region,
    const std::vector<ExactSourceOrderFact> &facts,
    ExactOrderRegionBuilder *builder);

inline bool exact_order_region_expr_relation_factor(
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool before,
    const bool inclusive,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex ||
      time_id == semantic::kInvalidIndex) {
    return false;
  }
  ExactOrderRegionExpr expr = exact_order_region_one();
  exact_order_region_append_expr_factor(
      &expr.terms.back(), expr_id, time_id, before, false, inclusive);
  *out = std::move(expr);
  return true;
}

inline bool exact_order_region_symbolic_transition_at_time(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &transition,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto release_source_id =
      exact_symbolic_transition_release_source_id(transition);
  if (release_source_id == semantic::kInvalidIndex) {
    return false;
  }
  ExactOrderRegionExpr region = exact_order_region_one();
  exact_order_region_append_exact(
      &region.terms.back(), release_source_id, time_id);

  ExactOrderRegionExpr readiness;
  if (!exact_order_region_truth_at_time(
          plan,
          transition.readiness_time_expr.requirements,
          time_id,
          builder,
          &readiness)) {
    return false;
  }
  region = exact_order_region_conjoin(std::move(region), std::move(readiness));

  ExactOrderRegionExpr guards;
  if (!exact_order_region_truth_at_time(
          plan,
          transition.guards,
          time_id,
          builder,
          &guards)) {
    return false;
  }
  region = exact_order_region_conjoin(std::move(region), std::move(guards));

  if (!exact_order_region_apply_source_order_facts(
          &region,
          transition.order_region.source_order_facts,
          builder)) {
    return false;
  }
  *out = std::move(region);
  return true;
}

inline bool exact_order_region_expr_transition_at_time(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex ||
      time_id == semantic::kInvalidIndex) {
    return false;
  }
  ExactOrderRegionExpr combined = exact_order_region_zero();
  auto scenarios = build_expr_transition_scenarios(plan, expr_id);
  for (const auto &scenario : scenarios) {
    ExactOrderRegionExpr branch;
    if (!exact_order_region_symbolic_transition_at_time(
            plan,
            scenario.transition,
            time_id,
            builder,
            &branch)) {
      return false;
    }
    combined = exact_order_region_union(
        std::move(combined), std::move(branch));
  }
  *out = exact_order_region_minimize_positive_union(
      exact_order_region_simplify(std::move(combined)));
  return true;
}

inline bool exact_order_region_expr_transition_before_or_at(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool inclusive,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto transition_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr at_time;
  if (!exact_order_region_expr_transition_at_time(
          plan,
          expr_id,
          transition_time_id,
          builder,
          &at_time)) {
    return false;
  }
  for (auto &term : at_time.terms) {
    exact_order_region_append_time_order(
        &term, transition_time_id, time_id, !inclusive);
  }
  *out = exact_order_region_simplify(std::move(at_time));
  return true;
}

inline bool exact_order_region_expr_transition_not_before_or_after(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool inclusive,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr combined = exact_order_region_one();
  auto scenarios = build_expr_transition_scenarios(plan, expr_id);
  for (const auto &scenario : scenarios) {
    const auto release_source_id =
        exact_symbolic_transition_release_source_id(scenario.transition);
    if (release_source_id == semantic::kInvalidIndex) {
      return false;
    }

    ExactOrderRegionExpr release_not_before = exact_order_region_one();
    exact_order_region_append_lower(
        &release_not_before.terms.back(),
        release_source_id,
        time_id,
        inclusive);

    const auto release_time_id = exact_order_region_new_time(builder);
    ExactOrderRegionExpr release_before = exact_order_region_one();
    exact_order_region_append_exact(
        &release_before.terms.back(),
        release_source_id,
        release_time_id);
    exact_order_region_append_time_order(
        &release_before.terms.back(),
        release_time_id,
        time_id,
        inclusive);

    ExactOrderRegionExpr requirements_fail;
    if (!exact_order_region_transition_requirements_fail_at_times(
            plan,
            scenario.transition,
            scenario.transition.guards,
            release_time_id,
            release_time_id,
            builder,
            &requirements_fail)) {
      return false;
    }
    ExactOrderRegionExpr scenario_not_before =
        exact_order_region_union(
            std::move(release_not_before),
            exact_order_region_conjoin(
                std::move(release_before),
                std::move(requirements_fail)));
    combined =
        exact_order_region_conjoin(
            std::move(combined),
            std::move(scenario_not_before));
  }
  *out = exact_order_region_minimize_positive_union(
      exact_order_region_simplify(std::move(combined)));
  return true;
}

inline bool exact_order_region_expr_satisfied_at_time(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto &kernel = plan.expr_kernels[pos];
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
    *out = exact_order_region_zero();
    return true;
  case semantic::ExprKind::TrueExpr:
    *out = exact_order_region_one();
    return true;
  case semantic::ExprKind::Event: {
    ExactOrderRegionExpr expr = exact_order_region_one();
    exact_order_region_append_upper(
        &expr.terms.back(), kernel.event_source_id, time_id, true);
    *out = std::move(expr);
    return true;
  }
  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    return exact_order_region_expr_not_before(
        plan, child, time_id, builder, out);
  }
  case semantic::ExprKind::And: {
    ExactOrderRegionExpr combined = exact_order_region_one();
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      ExactOrderRegionExpr child_satisfied;
      if (!exact_order_region_expr_satisfied_at_time(
              plan, child, time_id, builder, &child_satisfied)) {
        return false;
      }
      combined = exact_order_region_conjoin(
          std::move(combined), std::move(child_satisfied));
    }
    *out = std::move(combined);
    return true;
  }
  case semantic::ExprKind::Or: {
    ExactOrderRegionExpr combined = exact_order_region_zero();
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      ExactOrderRegionExpr child_satisfied;
      if (!exact_order_region_expr_satisfied_at_time(
              plan, child, time_id, builder, &child_satisfied)) {
        return false;
      }
      combined = exact_order_region_union(
          std::move(combined), std::move(child_satisfied));
    }
    *out = std::move(combined);
    return true;
  }
  case semantic::ExprKind::Guard:
    return exact_order_region_guard_before_or_at(
        plan, kernel, time_id, builder, out);
  }
  return false;
}

inline bool exact_order_region_expr_not_satisfied_at_time(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr satisfied;
  if (!exact_order_region_expr_satisfied_at_time(
          plan, expr_id, time_id, builder, &satisfied)) {
    return false;
  }
  *out = exact_order_region_subtract_region(
      exact_order_region_one(), satisfied);
  return true;
}

inline bool exact_order_region_guard_allowed_at_time(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto &program = plan.program;
  ExactOrderRegionExpr allowed;
  if (!exact_order_region_expr_after(
          plan,
          kernel.guard_blocker_expr_id,
          time_id,
          builder,
          &allowed)) {
    return false;
  }
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    const auto child =
        program.expr_args[
            static_cast<std::size_t>(kernel.children.offset + i)];
    ExactOrderRegionExpr unless_before;
    if (!exact_order_region_expr_before(
            plan, child, time_id, builder, &unless_before)) {
      return false;
    }
    allowed =
        allowed.terms.empty()
            ? std::move(unless_before)
            : exact_order_region_union(allowed, unless_before);
  }
  *out = std::move(allowed);
  return true;
}

inline bool exact_order_region_guard_before_or_at(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto ref_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr ref_at_time;
  if (!exact_order_region_expr_at_time(
          plan, kernel.guard_ref_expr_id, ref_time_id, builder, &ref_at_time)) {
    return false;
  }
  for (auto &term : ref_at_time.terms) {
    exact_order_region_append_time_order(&term, ref_time_id, time_id, false);
  }
  ExactOrderRegionExpr allowed;
  if (!exact_order_region_guard_allowed_at_time(
          plan, kernel, ref_time_id, builder, &allowed)) {
    return false;
  }
  *out =
      exact_order_region_conjoin(
          std::move(ref_at_time),
          std::move(allowed));
  return true;
}

inline bool exact_order_region_guard_blocked_at_time(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr blocked;
  if (!exact_order_region_expr_before_or_at(
          plan,
          kernel.guard_blocker_expr_id,
          time_id,
          builder,
          &blocked)) {
    return false;
  }
  const auto &program = plan.program;
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    const auto child =
        program.expr_args[
            static_cast<std::size_t>(kernel.children.offset + i)];
    ExactOrderRegionExpr unless_not_before;
    if (!exact_order_region_expr_not_before(
            plan, child, time_id, builder, &unless_not_before)) {
      return false;
    }
    blocked =
        exact_order_region_conjoin(
            std::move(blocked),
            std::move(unless_not_before));
  }
  *out = std::move(blocked);
  return true;
}

inline bool exact_order_region_guard_before(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  (void)expr_id;
  const auto ref_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr ref_at_time;
  if (!exact_order_region_expr_at_time(
          plan, kernel.guard_ref_expr_id, ref_time_id, builder, &ref_at_time)) {
    return false;
  }
  for (auto &term : ref_at_time.terms) {
    exact_order_region_append_time_order(&term, ref_time_id, time_id);
  }
  ExactOrderRegionExpr allowed;
  if (!exact_order_region_guard_allowed_at_time(
          plan, kernel, ref_time_id, builder, &allowed)) {
    return false;
  }
  *out =
      exact_order_region_conjoin(
          std::move(ref_at_time),
          std::move(allowed));
  return true;
}

inline bool exact_order_region_guard_not_before(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  (void)expr_id;
  ExactOrderRegionExpr ref_not_before;
  if (!exact_order_region_expr_not_before(
          plan, kernel.guard_ref_expr_id, time_id, builder, &ref_not_before)) {
    return false;
  }

  const auto ref_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr ref_at_time;
  if (!exact_order_region_expr_at_time(
          plan, kernel.guard_ref_expr_id, ref_time_id, builder, &ref_at_time)) {
    return false;
  }
  for (auto &term : ref_at_time.terms) {
    exact_order_region_append_time_order(&term, ref_time_id, time_id);
  }
  ExactOrderRegionExpr blocked;
  if (!exact_order_region_guard_blocked_at_time(
          plan, kernel, ref_time_id, builder, &blocked)) {
    return false;
  }
  *out =
      exact_order_region_union(
          std::move(ref_not_before),
          exact_order_region_conjoin(
              std::move(ref_at_time),
              std::move(blocked)));
  return true;
}

inline bool exact_order_region_guard_after(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr ref_after;
  if (!exact_order_region_expr_after(
          plan, kernel.guard_ref_expr_id, time_id, builder, &ref_after)) {
    return false;
  }

  const auto ref_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr ref_at_time;
  if (!exact_order_region_expr_at_time(
          plan, kernel.guard_ref_expr_id, ref_time_id, builder, &ref_at_time)) {
    return false;
  }
  for (auto &term : ref_at_time.terms) {
    exact_order_region_append_time_order(&term, ref_time_id, time_id, false);
  }
  ExactOrderRegionExpr blocked;
  if (!exact_order_region_guard_blocked_at_time(
          plan, kernel, ref_time_id, builder, &blocked)) {
    return false;
  }
  *out =
      exact_order_region_union(
          std::move(ref_after),
          exact_order_region_conjoin(
              std::move(ref_at_time),
          std::move(blocked)));
  return true;
}

inline bool exact_order_region_expr_relation_prefers_factor(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id) {
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_kernels.size() ||
      !exact_order_region_expr_relation_can_collapse(plan, expr_id)) {
    return false;
  }
  const auto kind =
      plan.expr_kernels[static_cast<std::size_t>(expr_id)].kind;
  return kind == semantic::ExprKind::And ||
         kind == semantic::ExprKind::Or ||
         kind == semantic::ExprKind::Guard;
}

inline bool exact_order_region_try_expr_relation_factor(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool before,
    const bool inclusive,
    ExactOrderRegionExpr *out) {
  if (!exact_order_region_expr_relation_prefers_factor(plan, expr_id)) {
    return false;
  }
  return exact_order_region_expr_relation_factor(
      expr_id, time_id, before, inclusive, out);
}

inline bool exact_order_region_expr_after(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (exact_order_region_try_expr_relation_factor(
          plan, expr_id, time_id, false, false, out)) {
    return true;
  }
  if (expr_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(expr_id) < plan.expr_kernels.size()) {
    const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
    if (kernel.kind == semantic::ExprKind::Guard) {
      return exact_order_region_guard_after(
          plan, kernel, time_id, builder, out);
    }
  }
  return exact_order_region_expr_transition_not_before_or_after(
      plan, expr_id, time_id, false, builder, out);
}

inline bool exact_order_region_expr_before_or_at(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (exact_order_region_try_expr_relation_factor(
          plan, expr_id, time_id, true, true, out)) {
    return true;
  }
  if (expr_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(expr_id) < plan.expr_kernels.size()) {
    const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
    if (kernel.kind == semantic::ExprKind::Guard) {
      return exact_order_region_guard_before_or_at(
          plan, kernel, time_id, builder, out);
    }
  }
  return exact_order_region_expr_transition_before_or_at(
      plan, expr_id, time_id, true, builder, out);
}

inline bool exact_order_region_guard_at_time(
    const ExactVariantBuildState &plan,
    const ExactExprKernel &kernel,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr ref_at_time;
  if (!exact_order_region_expr_at_time(
          plan, kernel.guard_ref_expr_id, time_id, builder, &ref_at_time)) {
    return false;
  }
  ExactOrderRegionExpr allowed;
  if (!exact_order_region_guard_allowed_at_time(
          plan, kernel, time_id, builder, &allowed)) {
    return false;
  }
  *out =
      exact_order_region_conjoin(
          std::move(ref_at_time),
          std::move(allowed));
  return true;
}

inline bool exact_order_region_expr_before(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (exact_order_region_try_expr_relation_factor(
          plan, expr_id, time_id, true, false, out)) {
    return true;
  }
  if (expr_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(expr_id) < plan.expr_kernels.size()) {
    const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
    if (kernel.kind == semantic::ExprKind::Guard) {
      return exact_order_region_guard_before(
          plan, kernel, expr_id, time_id, builder, out);
    }
  }
  return exact_order_region_expr_transition_before_or_at(
      plan, expr_id, time_id, false, builder, out);
}

inline bool exact_order_region_expr_not_before(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (exact_order_region_try_expr_relation_factor(
          plan, expr_id, time_id, false, true, out)) {
    return true;
  }
  if (expr_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(expr_id) < plan.expr_kernels.size()) {
    const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
    if (kernel.kind == semantic::ExprKind::Guard) {
      return exact_order_region_guard_not_before(
          plan, kernel, expr_id, time_id, builder, out);
    }
  }
  return exact_order_region_expr_transition_not_before_or_after(
      plan, expr_id, time_id, true, builder, out);
}

inline bool exact_order_region_expr_at_time(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(expr_id) < plan.expr_kernels.size()) {
    const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
    if (kernel.kind == semantic::ExprKind::Guard) {
      return exact_order_region_guard_at_time(
          plan, kernel, time_id, builder, out);
    }
  }
  return exact_order_region_expr_transition_at_time(
      plan, expr_id, time_id, builder, out);
}

inline bool exact_order_region_append_source_order_fact(
    ExactRegionCell *term,
    const ExactSourceOrderFact &fact,
    ExactOrderRegionBuilder *builder) {
  const auto before_time =
      exact_order_region_exact_time_for_source(*term, fact.before_source_id);
  const auto after_time =
      exact_order_region_exact_time_for_source(*term, fact.after_source_id);
  if (before_time != semantic::kInvalidIndex &&
      after_time != semantic::kInvalidIndex) {
    exact_order_region_append_time_order(term, before_time, after_time);
    return !term->impossible;
  }
  if (before_time != semantic::kInvalidIndex) {
    exact_order_region_append_lower(term, fact.after_source_id, before_time);
    return !term->impossible;
  }
  if (after_time != semantic::kInvalidIndex) {
    exact_order_region_append_upper(term, fact.before_source_id, after_time);
    return !term->impossible;
  }
  const auto upper_bounds = exact_region_upper_source_atoms(*term);
  const auto has_after_upper =
      std::any_of(
          upper_bounds.begin(),
          upper_bounds.end(),
          [&](const ExactOrderRegionSourceTime &bound) {
            return bound.source_id == fact.after_source_id;
          });
  if (has_after_upper && builder != nullptr) {
    const auto ordered_after_time = exact_order_region_new_time(builder);
    exact_order_region_append_exact(
        term, fact.after_source_id, ordered_after_time);
    exact_order_region_append_upper(
        term, fact.before_source_id, ordered_after_time);
    return !term->impossible;
  }
  return true;
}

inline bool exact_order_region_target_scenario(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &transition,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  return exact_order_region_symbolic_transition_at_time(
      plan, transition, observed_time_id, builder, out);
}

inline bool exact_order_region_transition_requirements_at_times(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &transition,
    const ExactTransitionGuardSet &guards,
    const semantic::Index release_time_id,
    const semantic::Index readiness_time_id,
    const semantic::Index guard_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto release_source_id =
      exact_symbolic_transition_release_source_id(transition);
  ExactOrderRegionExpr region = exact_order_region_one();
  exact_order_region_append_exact(
      &region.terms.back(), release_source_id, release_time_id);

  ExactOrderRegionExpr readiness;
  if (!exact_order_region_truth_at_time(
          plan,
          transition.readiness_time_expr.requirements,
          readiness_time_id,
          builder,
          &readiness)) {
    return false;
  }
  region = exact_order_region_conjoin(region, readiness);

  ExactOrderRegionExpr blockers;
  if (!exact_order_region_truth_at_time(
          plan,
          guards,
          guard_time_id,
          builder,
          &blockers)) {
    return false;
  }
  region = exact_order_region_conjoin(region, blockers);

  *out = std::move(region);
  return true;
}

inline bool exact_order_region_readiness_transition_at_time(
    const ExactVariantBuildState &plan,
    const ExactTransitionGuardSet &readiness,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr combined = exact_order_region_zero();

  for (const auto &active_guard : readiness.guards) {
    if (active_guard.kind != ExactTransitionGuardKind::SourceBefore) {
      continue;
    }
    ExactOrderRegionExpr branch = exact_order_region_one();
    exact_order_region_append_exact(
        &branch.terms.back(), active_guard.subject_id, time_id);
    for (const auto &guard : readiness.guards) {
      if (guard.kind == ExactTransitionGuardKind::SourceBefore &&
          guard.subject_id == active_guard.subject_id) {
        continue;
      }
      if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
        exact_order_region_append_upper(
            &branch.terms.back(), guard.subject_id, time_id, true);
      } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
        ExactOrderRegionExpr before;
        if (!exact_order_region_expr_satisfied_at_time(
                plan, guard.subject_id, time_id, builder, &before)) {
          return false;
        }
        branch = exact_order_region_conjoin(branch, before);
      } else {
        return false;
      }
    }
    combined =
        exact_order_region_union(
            std::move(combined),
            std::move(branch));
  }

  for (const auto &active_guard : readiness.guards) {
    if (active_guard.kind != ExactTransitionGuardKind::ExprBefore) {
      continue;
    }
    ExactOrderRegionExpr branch;
    if (!exact_order_region_expr_at_time(
            plan, active_guard.subject_id, time_id, builder, &branch)) {
      return false;
    }
    for (const auto &guard : readiness.guards) {
      if (guard.kind == ExactTransitionGuardKind::ExprBefore &&
          guard.subject_id == active_guard.subject_id) {
        continue;
      }
      if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
        for (auto &term : branch.terms) {
          exact_order_region_append_upper(
              &term, guard.subject_id, time_id, true);
        }
      } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
        ExactOrderRegionExpr before;
        if (!exact_order_region_expr_satisfied_at_time(
                plan, guard.subject_id, time_id, builder, &before)) {
          return false;
        }
        branch = exact_order_region_conjoin(branch, before);
      } else {
        return false;
      }
    }
    combined =
        exact_order_region_union(
            std::move(combined),
            std::move(branch));
  }

  *out = std::move(combined);
  return true;
}

inline bool exact_order_region_guard_at_time(
    const ExactVariantBuildState &plan,
    const ExactTransitionGuard &guard,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr region = exact_order_region_one();
  if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
    exact_order_region_append_upper(
        &region.terms.back(), guard.subject_id, time_id, true);
  } else if (guard.kind == ExactTransitionGuardKind::SourceAfter) {
    exact_order_region_append_lower(
        &region.terms.back(), guard.subject_id, time_id);
  } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
    if (exact_order_region_expr_relation_can_collapse(
            plan, guard.subject_id)) {
      return exact_order_region_expr_relation_factor(
          guard.subject_id, time_id, true, true, out);
    }
    return exact_order_region_expr_satisfied_at_time(
        plan, guard.subject_id, time_id, builder, out);
  } else if (guard.kind == ExactTransitionGuardKind::ExprAfter) {
    if (exact_order_region_expr_relation_can_collapse(
            plan, guard.subject_id)) {
      return exact_order_region_expr_relation_factor(
          guard.subject_id, time_id, false, true, out);
    }
    return exact_order_region_expr_after(
        plan, guard.subject_id, time_id, builder, out);
  }
  *out = std::move(region);
  return true;
}

inline bool exact_order_region_truth_at_time(
    const ExactVariantBuildState &plan,
    const ExactTransitionGuardSet &guards,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr region =
      guards.empty_value == 0.0
          ? exact_order_region_zero()
          : exact_order_region_one();
  for (const auto &guard : guards.guards) {
    ExactOrderRegionExpr guard_region;
    if (!exact_order_region_guard_at_time(
            plan, guard, time_id, builder, &guard_region)) {
      return false;
    }
    region = exact_order_region_conjoin(region, guard_region);
  }
  *out = std::move(region);
  return true;
}

inline bool exact_order_region_apply_source_order_facts(
    ExactOrderRegionExpr *region,
    const std::vector<ExactSourceOrderFact> &facts,
    ExactOrderRegionBuilder *builder) {
  for (const auto &fact : facts) {
    for (auto &term : region->terms) {
      exact_order_region_append_source_order_fact(&term, fact, builder);
    }
  }
  region->terms.erase(
      std::remove_if(
          region->terms.begin(),
          region->terms.end(),
          [](const ExactRegionCell &term) {
            return term.impossible || term.sign == 0.0;
          }),
      region->terms.end());
  return true;
}

struct ExactOrderRegionTargetBranch {
  ExactOrderRegionExpr expr;
  semantic::Index readiness_time_id{semantic::kInvalidIndex};
};

inline bool exact_order_region_target_branches(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &transition,
    const bool expose_readiness_time,
    ExactOrderRegionBuilder *builder,
    std::vector<ExactOrderRegionTargetBranch> *out) {
  if (!expose_readiness_time ||
      !exact_symbolic_transition_has_readiness(transition)) {
    ExactOrderRegionExpr target;
    if (!exact_order_region_target_scenario(
            plan, transition, builder, &target)) {
      return false;
    }
    out->push_back(
        ExactOrderRegionTargetBranch{
            std::move(target),
            semantic::kInvalidIndex});
    return true;
  }

  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto zero_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero);
  const auto release_source_id =
      exact_symbolic_transition_release_source_id(transition);

  ExactOrderRegionExpr base = exact_order_region_one();
  exact_order_region_append_exact(
      &base.terms.back(), release_source_id, observed_time_id);
  ExactOrderRegionExpr tail;
  if (!exact_order_region_truth_at_time(
          plan, transition.guards, observed_time_id, builder, &tail)) {
    return false;
  }
  base = exact_order_region_conjoin(base, tail);

  auto append_branch =
      [&](ExactOrderRegionExpr branch,
      const semantic::Index readiness_time_id) -> bool {
        if (!exact_order_region_apply_source_order_facts(
                &branch,
                transition.order_region.source_order_facts,
                builder)) {
          return false;
        }
        out->push_back(
            ExactOrderRegionTargetBranch{std::move(branch), readiness_time_id});
        return true;
      };

  ExactOrderRegionExpr initial_ready;
  if (!exact_order_region_truth_at_time(
          plan,
          transition.readiness_time_expr.requirements,
          zero_time_id,
          builder,
          &initial_ready)) {
    return false;
  }
  if (!initial_ready.terms.empty()) {
    if (!append_branch(
            exact_order_region_conjoin(base, initial_ready),
            zero_time_id)) {
      return false;
    }
  }

  const auto readiness_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr readiness;
  if (!exact_order_region_readiness_transition_at_time(
          plan,
          transition.readiness_time_expr.requirements,
          readiness_time_id,
          builder,
          &readiness)) {
    return false;
  }
  for (auto &readiness_term : readiness.terms) {
    exact_order_region_append_time_order(
        &readiness_term, readiness_time_id, observed_time_id);
  }
  if (!append_branch(
          exact_order_region_conjoin(base, readiness),
          readiness_time_id)) {
    return false;
  }
  return true;
}

inline bool exact_order_region_needs_target_readiness_time(
    const ExactOutcomeRegionCompileContext &outcome_context,
    const ExactSymbolicTransitionTime &transition) {
  if (!exact_symbolic_transition_has_readiness(transition)) {
    return false;
  }
  for (const auto &competitor : outcome_context.competitors) {
    for (const auto &scenario : competitor.scenarios) {
      const auto relation =
          exact_symbolic_transition_relation(transition, scenario.transition);
      if (relation.competitor_can_positively_coincide) {
        return true;
      }
    }
  }
  return false;
}

inline bool exact_transition_guard_equal(
    const ExactTransitionGuard &lhs,
    const ExactTransitionGuard &rhs) noexcept {
  return lhs.kind == rhs.kind && lhs.subject_id == rhs.subject_id;
}

inline void exact_order_region_erase_transition_guard_once(
    std::vector<ExactTransitionGuard> *guards,
    const ExactTransitionGuard &guard) {
  const auto it = std::find_if(
      guards->begin(),
      guards->end(),
      [&](const ExactTransitionGuard &existing) {
        return exact_transition_guard_equal(existing, guard);
      });
  if (it != guards->end()) {
    guards->erase(it);
  }
}

inline ExactTransitionGuardSet exact_order_region_guards_excluding(
    ExactTransitionGuardSet guards,
    const ExactTransitionGuardSet &covered_guards) {
  for (const auto &guard : covered_guards.guards) {
    exact_order_region_erase_transition_guard_once(
        &guards.guards, guard);
  }
  return guards;
}

inline void exact_order_region_append_union(
    ExactOrderRegionExpr *dst,
    ExactOrderRegionExpr src) {
  if (src.terms.empty()) {
    return;
  }
  if (dst->terms.empty()) {
    *dst = std::move(src);
    return;
  }
  *dst = exact_order_region_union(*dst, src);
}

inline bool exact_order_region_guard_not_at_time(
    const ExactVariantBuildState &plan,
    const ExactTransitionGuard &guard,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr failure = exact_order_region_one();
  if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
    exact_order_region_append_lower(
        &failure.terms.back(), guard.subject_id, time_id, true);
  } else if (guard.kind == ExactTransitionGuardKind::SourceAfter) {
    exact_order_region_append_upper(
        &failure.terms.back(), guard.subject_id, time_id, true);
  } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
    if (exact_order_region_expr_relation_can_collapse(
            plan, guard.subject_id)) {
      return exact_order_region_expr_relation_factor(
          guard.subject_id, time_id, false, true, out);
    }
    return exact_order_region_expr_not_satisfied_at_time(
        plan, guard.subject_id, time_id, builder, out);
  } else if (guard.kind == ExactTransitionGuardKind::ExprAfter) {
    if (exact_order_region_expr_relation_can_collapse(
            plan, guard.subject_id)) {
      return exact_order_region_expr_relation_factor(
          guard.subject_id, time_id, true, true, out);
    }
    return exact_order_region_expr_before_or_at(
        plan, guard.subject_id, time_id, builder, out);
  }
  *out = std::move(failure);
  return true;
}

inline bool exact_order_region_truth_not_at_time(
    const ExactVariantBuildState &plan,
    const ExactTransitionGuardSet &guards,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (guards.guards.empty()) {
    *out = guards.empty_value == 0.0
               ? exact_order_region_one()
               : exact_order_region_zero();
    return true;
  }
  ExactOrderRegionExpr combined = exact_order_region_zero();
  for (const auto &guard : guards.guards) {
    ExactOrderRegionExpr failure;
    if (!exact_order_region_guard_not_at_time(
            plan, guard, time_id, builder, &failure)) {
      return false;
    }
    exact_order_region_append_union(&combined, std::move(failure));
  }
  *out = std::move(combined);
  return true;
}

inline bool exact_order_region_transition_requirements_fail_at_times(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &transition,
    const ExactTransitionGuardSet &guards,
    const semantic::Index readiness_time_id,
    const semantic::Index guard_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr failure = exact_order_region_zero();
  if (exact_symbolic_transition_has_readiness(transition)) {
    ExactOrderRegionExpr readiness_fail;
    if (!exact_order_region_truth_not_at_time(
            plan,
            transition.readiness_time_expr.requirements,
            readiness_time_id,
            builder,
            &readiness_fail)) {
      return false;
    }
    exact_order_region_append_union(&failure, std::move(readiness_fail));
  }

  ExactOrderRegionExpr blocker_fail;
  if (!exact_order_region_truth_not_at_time(
          plan,
          guards,
          guard_time_id,
          builder,
          &blocker_fail)) {
    return false;
  }
  exact_order_region_append_union(&failure, std::move(blocker_fail));
  *out = std::move(failure);
  return true;
}

inline ExactOrderRegionExpr exact_order_region_subtract_region(
    ExactOrderRegionExpr domain,
    const ExactOrderRegionExpr &covered) {
  domain = exact_order_region_simplify(std::move(domain));
  for (const auto &term : exact_order_region_simplify(covered).terms) {
    domain = exact_order_region_subtract_expr(domain, term);
    if (domain.terms.empty()) {
      break;
    }
  }
  return exact_order_region_simplify(std::move(domain));
}

inline bool exact_order_region_expr_has_positive_measure(
    ExactOrderRegionExpr region) {
  region =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(std::move(region)));
  for (const auto &term : region.terms) {
    if (term.sign > 0.0 &&
        exact_order_region_cell_has_positive_measure(term)) {
      return true;
    }
  }
  return false;
}

inline bool exact_order_region_expr_is_tautology(
    ExactOrderRegionExpr region) {
  return !exact_order_region_expr_has_positive_measure(
      exact_order_region_subtract_region(
          exact_order_region_one(),
          exact_order_region_simplify(std::move(region))));
}

inline bool exact_order_region_materialize_relation_factors(
    const ExactVariantBuildState &plan,
    ExactOrderRegionExpr region,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const ExactProjectionRelationOps ops{
      nullptr,
      nullptr,
      exact_order_region_expand_relation_factor};
  region =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(std::move(region)));
  bool changed = true;
  while (changed) {
    changed = false;
    ExactOrderRegionExpr next = exact_order_region_zero();
    for (const auto &term : region.terms) {
      const auto factors = exact_projection_materializable_relation_factors(term);
      bool expanded = false;
      for (const auto &factor : factors) {
        ExactOrderRegionExpr materialized;
        if (!exact_projection_materialize_relation_factor(
                plan,
                term,
                factor,
                ops,
                builder,
                &materialized)) {
          continue;
        }
        exact_order_region_append_expr(&next, std::move(materialized));
        changed = true;
        expanded = true;
        break;
      }
      if (!expanded) {
        next.terms.push_back(term);
      }
    }
    region =
        exact_order_region_minimize_positive_union(
            exact_order_region_simplify(std::move(next)));
  }
  *out = std::move(region);
  return true;
}

inline bool exact_order_region_apply_transition_order_facts(
    ExactOrderRegionExpr *region,
    const ExactSymbolicTransitionTime &transition,
    ExactOrderRegionBuilder *builder) {
  return exact_order_region_apply_source_order_facts(
      region,
      transition.order_region.source_order_facts,
      builder);
}

inline bool exact_order_region_transition_strict_precedence(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &competitor,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto competitor_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr wins;
  if (!exact_order_region_symbolic_transition_at_time(
          plan,
          competitor,
          competitor_time_id,
          builder,
          &wins)) {
    return false;
  }
  for (auto &term : wins.terms) {
    exact_order_region_append_time_order(
        &term, competitor_time_id, observed_time_id);
  }
  *out =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(std::move(wins)));
  return true;
}

inline bool exact_order_region_transition_coincident_precedence(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &competitor,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto readiness_time_id =
      target_readiness_time_id == semantic::kInvalidIndex
          ? static_cast<semantic::Index>(CompiledMathTimeSlot::Zero)
          : target_readiness_time_id;
  ExactOrderRegionExpr wins;
  if (!exact_order_region_transition_requirements_at_times(
          plan,
          competitor,
          competitor.guards,
          observed_time_id,
          readiness_time_id,
          observed_time_id,
          builder,
          &wins)) {
    return false;
  }
  if (!exact_order_region_apply_transition_order_facts(
          &wins, competitor, builder)) {
    return false;
  }
  *out =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(std::move(wins)));
  return true;
}

inline bool exact_order_region_transition_precedence(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &target,
    const ExactSymbolicTransitionTime &competitor,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto relation =
      exact_symbolic_transition_relation(target, competitor);
  ExactOrderRegionExpr wins = exact_order_region_zero();
  if (relation.competitor_can_strictly_precede) {
    ExactOrderRegionExpr strict_wins;
    if (!exact_order_region_transition_strict_precedence(
            plan,
            competitor,
            builder,
            &strict_wins)) {
      return false;
    }
    wins = exact_order_region_union(std::move(wins), std::move(strict_wins));
  }
  if (relation.competitor_can_positively_coincide) {
    ExactOrderRegionExpr coincident_wins;
    if (!exact_order_region_transition_coincident_precedence(
            plan,
            competitor,
            target_readiness_time_id,
            builder,
            &coincident_wins)) {
      return false;
    }
    wins =
        exact_order_region_union(
            std::move(wins),
            std::move(coincident_wins));
  }
  *out =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(std::move(wins)));
  return true;
}

inline bool exact_order_region_competitor_scenario_can_win_on_branch(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionScenario &target,
    const ExactSymbolicTransitionScenario &competitor,
    const ExactOrderRegionExpr &target_branch,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder builder,
    bool *can_win) {
  ExactOrderRegionExpr wins;
  if (!exact_order_region_transition_precedence(
          plan,
          target.transition,
          competitor.transition,
          target_readiness_time_id,
          &builder,
          &wins)) {
    return false;
  }
  ExactOrderRegionExpr overlap =
      exact_order_region_conjoin(target_branch, wins);
  if (!exact_order_region_materialize_relation_factors(
          plan,
          std::move(overlap),
          &builder,
          &overlap)) {
    return false;
  }
  *can_win = exact_order_region_expr_has_positive_measure(std::move(overlap));
  return true;
}

inline bool exact_order_region_transition_no_strict_precedence(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionRelation &relation,
    const ExactSymbolicTransitionTime &competitor,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (!relation.competitor_can_strictly_precede) {
    *out = exact_order_region_one();
    return true;
  }
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto competitor_source_id =
      exact_symbolic_transition_release_source_id(competitor);
  ExactOrderRegionExpr release_not_before = exact_order_region_one();
  exact_order_region_append_lower(
      &release_not_before.terms.back(),
      competitor_source_id,
      observed_time_id,
      true);

  const auto competitor_time_id = exact_order_region_new_time(builder);
  ExactOrderRegionExpr release_before = exact_order_region_one();
  exact_order_region_append_exact(
      &release_before.terms.back(),
      competitor_source_id,
      competitor_time_id);
  exact_order_region_append_time_order(
      &release_before.terms.back(), competitor_time_id, observed_time_id);

  ExactOrderRegionExpr requirements_fail;
  if (!exact_order_region_transition_requirements_fail_at_times(
          plan,
          competitor,
          competitor.guards,
          competitor_time_id,
          competitor_time_id,
          builder,
          &requirements_fail)) {
    return false;
  }
  for (auto &term : requirements_fail.terms) {
    exact_order_region_append_time_order(
        &term, competitor_time_id, observed_time_id);
  }

  *out = exact_order_region_union(
      release_not_before,
      exact_order_region_conjoin(
          std::move(release_before),
          std::move(requirements_fail)));
  return true;
}

inline bool exact_order_region_transition_no_coincident_precedence(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &target,
    const ExactSymbolicTransitionRelation &relation,
    const ExactSymbolicTransitionTime &competitor,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (!relation.competitor_can_positively_coincide) {
    *out = exact_order_region_one();
    return true;
  }
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto readiness_time_id =
      target_readiness_time_id == semantic::kInvalidIndex
          ? static_cast<semantic::Index>(CompiledMathTimeSlot::Zero)
          : target_readiness_time_id;
  const auto guards =
      exact_order_region_guards_excluding(
          competitor.guards,
          target.guards);
  return exact_order_region_transition_requirements_fail_at_times(
      plan,
      competitor,
      guards,
      readiness_time_id,
      observed_time_id,
      builder,
      out);
}

inline bool exact_order_region_transition_not_before_target(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &target,
    const ExactSymbolicTransitionTime &competitor,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto relation =
      exact_symbolic_transition_relation(target, competitor);
  ExactOrderRegionExpr no_strict_precedence;
  if (!exact_order_region_transition_no_strict_precedence(
          plan,
          relation,
          competitor,
          builder,
          &no_strict_precedence)) {
    return false;
  }
  ExactOrderRegionExpr no_coincident_precedence;
  if (!exact_order_region_transition_no_coincident_precedence(
          plan,
          target,
          relation,
          competitor,
          target_readiness_time_id,
          builder,
          &no_coincident_precedence)) {
    return false;
  }
  *out = exact_order_region_conjoin(
      no_strict_precedence, no_coincident_precedence);
  return true;
}

inline bool exact_order_region_competitor_scenario_non_win(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionScenario &target,
    const ExactSymbolicTransitionScenario &competitor,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  return exact_order_region_transition_not_before_target(
      plan,
      target.transition,
      competitor.transition,
      target_readiness_time_id,
      builder,
      out);
}

inline bool exact_order_region_transition_support_overlaps_expr(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionTime &transition,
    const semantic::Index expr_id) {
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_supports.size()) {
    return true;
  }
  const auto &expr_support =
      plan.expr_supports[static_cast<std::size_t>(expr_id)];
  const auto source_overlaps = [&](const semantic::Index source_id) {
    return source_id != semantic::kInvalidIndex &&
           supports_overlap(
               planned_source_support_for_id(plan, source_id),
               expr_support);
  };
  if (source_overlaps(exact_symbolic_transition_release_source_id(transition))) {
    return true;
  }
  for (const auto source_id : transition.active_sources) {
    if (source_overlaps(source_id)) {
      return true;
    }
  }
  const auto guard_overlaps =
      [&](const ExactTransitionGuard &guard) {
        if (guard.kind == ExactTransitionGuardKind::SourceBefore ||
            guard.kind == ExactTransitionGuardKind::SourceAfter) {
          return source_overlaps(guard.subject_id);
        }
        return expr_supports_overlap(plan, guard.subject_id, expr_id);
      };
  for (const auto &guard :
       transition.readiness_time_expr.requirements.guards) {
    if (guard_overlaps(guard)) {
      return true;
    }
  }
  for (const auto &guard : transition.guards.guards) {
    if (guard_overlaps(guard)) {
      return true;
    }
  }
  return false;
}

inline bool exact_order_region_competitor_plan_non_win(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionScenario &target,
    const ExactCompetitorRegionPlan &competitor,
    const ExactOrderRegionExpr &target_branch,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (competitor.expr_root != semantic::kInvalidIndex &&
      exact_order_region_expr_relation_can_collapse(
          plan, competitor.expr_root) &&
      !exact_order_region_transition_support_overlaps_expr(
          plan, target.transition, competitor.expr_root)) {
    const auto observed_time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
    return exact_order_region_expr_relation_factor(
        competitor.expr_root, observed_time_id, false, true, out);
  }

  ExactOrderRegionExpr non_win = exact_order_region_one();
  for (const auto &scenario : competitor.scenarios) {
    bool can_win = true;
    if (!exact_order_region_competitor_scenario_can_win_on_branch(
            plan,
            target,
            scenario,
            target_branch,
            target_readiness_time_id,
            *builder,
            &can_win)) {
      return false;
    }
    if (!can_win) {
      continue;
    }
    ExactOrderRegionExpr scenario_non_win;
    if (!exact_order_region_competitor_scenario_non_win(
            plan,
            target,
            scenario,
            target_readiness_time_id,
            builder,
            &scenario_non_win)) {
      return false;
    }
    non_win = exact_order_region_conjoin(non_win, scenario_non_win);
  }
  *out = std::move(non_win);
  return true;
}

inline ExactOrderRegionExpr exact_order_region_gate_non_win(
    ExactOrderRegionExpr non_win,
    const std::vector<semantic::Index> &outcome_indices) {
  if (outcome_indices.empty() ||
      exact_order_region_expr_is_tautology(non_win)) {
    return non_win;
  }
  return exact_order_region_union(
      exact_order_region_with_outcome_used_gate(
          exact_order_region_one(),
          outcome_indices),
      exact_order_region_with_outcome_gate(
          std::move(non_win),
          outcome_indices));
}

inline bool exact_order_region_competitor_non_win(
    const ExactVariantBuildState &plan,
    const ExactOutcomeRegionCompileContext &outcome_context,
    const ExactSymbolicTransitionScenario &target,
    const ExactOrderRegionExpr &target_branch,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr non_win = exact_order_region_one();
  for (const auto &competitor : outcome_context.competitors) {
    ExactOrderRegionExpr competitor_non_win;
    if (!exact_order_region_competitor_plan_non_win(
            plan,
            target,
            competitor,
            target_branch,
            target_readiness_time_id,
            builder,
            &competitor_non_win)) {
      return false;
    }
    non_win =
        exact_order_region_conjoin(
            non_win,
            exact_order_region_gate_non_win(
                std::move(competitor_non_win),
                competitor.outcome_indices));
  }
  *out = std::move(non_win);
  return true;
}

inline bool exact_order_region_factor_context_overlaps_expr(
    const ExactVariantBuildState &plan,
    const ExactRegionCell &term,
    const semantic::Index expr_id) {
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_supports.size()) {
    return false;
  }
  const auto &support = plan.expr_supports[static_cast<std::size_t>(expr_id)];
  const auto source_overlaps =
      [&](const ExactOrderRegionSourceTime &source) {
        return support_contains_source(support, source.source_id);
      };
  const auto exact_sources = exact_region_exact_source_atoms(term);
  if (std::any_of(exact_sources.begin(), exact_sources.end(), source_overlaps)) {
    return true;
  }
  const auto lower_sources = exact_region_lower_source_atoms(term);
  if (std::any_of(lower_sources.begin(), lower_sources.end(), source_overlaps)) {
    return true;
  }
  const auto upper_sources = exact_region_upper_source_atoms(term);
  if (std::any_of(upper_sources.begin(), upper_sources.end(), source_overlaps)) {
    return true;
  }
  for (const auto &factor : exact_region_expr_atoms(term)) {
    if (factor.expr_id != expr_id &&
        expr_supports_overlap(plan, expr_id, factor.expr_id)) {
      return true;
    }
  }
  return false;
}

inline bool exact_order_region_expr_relation_can_collapse(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id) {
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_kernels.size()) {
    return false;
  }
  const auto &program = plan.program;
  const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::TrueExpr:
  case semantic::ExprKind::Event:
    return true;
  case semantic::ExprKind::Not:
    return false;
  case semantic::ExprKind::And:
  case semantic::ExprKind::Or:
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      if (!exact_order_region_expr_relation_can_collapse(plan, child)) {
        return false;
      }
    }
    return true;
  case semantic::ExprKind::Guard:
    if (!exact_order_region_expr_relation_can_collapse(
            plan, kernel.guard_ref_expr_id) ||
        !exact_order_region_expr_relation_can_collapse(
            plan, kernel.guard_blocker_expr_id)) {
      return false;
    }
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      if (!exact_order_region_expr_relation_can_collapse(plan, child)) {
        return false;
      }
    }
    return true;
  }
  return false;
}

inline bool exact_expr_completion_monotone(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id) {
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_kernels.size()) {
    return false;
  }
  const auto &program = plan.program;
  const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::TrueExpr:
  case semantic::ExprKind::Event:
    return true;
  case semantic::ExprKind::And:
  case semantic::ExprKind::Or:
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      if (!exact_expr_completion_monotone(plan, child)) {
        return false;
      }
    }
    return true;
  case semantic::ExprKind::Guard:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool exact_order_region_monotone_expr_before_or_at(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool inclusive,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  (void)builder;
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_kernels.size()) {
    return false;
  }
  const auto &program = plan.program;
  const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
    *out = exact_order_region_zero();
    return true;
  case semantic::ExprKind::TrueExpr:
    *out = exact_order_region_one();
    return true;
  case semantic::ExprKind::Event: {
    ExactOrderRegionExpr expr = exact_order_region_one();
    exact_order_region_append_upper(
        &expr.terms.back(), kernel.event_source_id, time_id, inclusive);
    *out = std::move(expr);
    return true;
  }
  case semantic::ExprKind::And: {
    ExactOrderRegionExpr combined = exact_order_region_one();
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      ExactOrderRegionExpr child_region;
      if (!exact_order_region_monotone_expr_before_or_at(
              plan, child, time_id, inclusive, builder, &child_region)) {
        return false;
      }
      combined =
          exact_order_region_conjoin(
              std::move(combined), std::move(child_region));
    }
    *out = exact_order_region_minimize_positive_union(
        exact_order_region_simplify(std::move(combined)));
    return true;
  }
  case semantic::ExprKind::Or: {
    ExactOrderRegionExpr combined = exact_order_region_zero();
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      ExactOrderRegionExpr child_region;
      if (!exact_order_region_monotone_expr_before_or_at(
              plan, child, time_id, inclusive, builder, &child_region)) {
        return false;
      }
      combined =
          exact_order_region_union(
              std::move(combined), std::move(child_region));
    }
    *out = exact_order_region_minimize_positive_union(
        exact_order_region_simplify(std::move(combined)));
    return true;
  }
  case semantic::ExprKind::Guard:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool exact_order_region_monotone_expr_not_before_or_after(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool inclusive,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr before;
  if (!exact_order_region_monotone_expr_before_or_at(
          plan, expr_id, time_id, !inclusive, builder, &before)) {
    return false;
  }
  *out = exact_order_region_minimize_positive_union(
      exact_order_region_subtract_region(exact_order_region_one(), before));
  return true;
}

inline bool exact_order_region_expand_relation_factor(
    const ExactVariantBuildState &plan,
    const ExactOrderRegionExprValueFactor &factor,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (factor.density) {
    return false;
  }
  if (exact_expr_completion_monotone(plan, factor.expr_id)) {
    if (factor.before) {
      return exact_order_region_monotone_expr_before_or_at(
          plan,
          factor.expr_id,
          factor.time_id,
          factor.inclusive,
          builder,
          out);
    }
    return exact_order_region_monotone_expr_not_before_or_after(
        plan,
        factor.expr_id,
        factor.time_id,
        factor.inclusive,
        builder,
        out);
  }
  if (factor.expr_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(factor.expr_id) < plan.expr_kernels.size()) {
    const auto &kernel =
        plan.expr_kernels[static_cast<std::size_t>(factor.expr_id)];
    if (kernel.kind == semantic::ExprKind::Guard) {
      if (factor.before && factor.inclusive) {
        return exact_order_region_guard_before_or_at(
            plan, kernel, factor.time_id, builder, out);
      }
      if (factor.before) {
        return exact_order_region_guard_before(
            plan, kernel, factor.expr_id, factor.time_id, builder, out);
      }
      if (factor.inclusive) {
        return exact_order_region_guard_not_before(
            plan, kernel, factor.expr_id, factor.time_id, builder, out);
      }
      return exact_order_region_guard_after(
          plan, kernel, factor.time_id, builder, out);
    }
  }
  if (factor.before && factor.inclusive) {
    return exact_order_region_expr_transition_before_or_at(
        plan, factor.expr_id, factor.time_id, true, builder, out);
  }
  if (factor.before) {
    return exact_order_region_expr_transition_before_or_at(
        plan, factor.expr_id, factor.time_id, false, builder, out);
  }
  if (factor.inclusive) {
    return exact_order_region_expr_transition_not_before_or_after(
        plan, factor.expr_id, factor.time_id, true, builder, out);
  }
  return exact_order_region_expr_transition_not_before_or_after(
      plan, factor.expr_id, factor.time_id, false, builder, out);
}

inline void exact_complexity_observe_region(
    ExactVariantBuildState *plan,
    const ExactOrderRegionExpr &region) {
  auto &metrics = plan->complexity;
  ++metrics.symbolic_region_count;
  const auto cell_count =
      static_cast<semantic::Index>(region.terms.size());
  metrics.symbolic_cell_count += cell_count;
  metrics.max_symbolic_cells_per_region =
      std::max(metrics.max_symbolic_cells_per_region, cell_count);
  for (std::size_t i = 0; i < region.terms.size(); ++i) {
    const auto &term = region.terms[i];
    if (term.sign < 0.0) {
      ++metrics.negative_symbolic_cell_count;
    }
    for (const auto &factor : exact_region_expr_atoms(term)) {
      if (!factor.density) {
        ++metrics.expr_relation_atom_count;
      }
    }
    for (std::size_t j = i + 1U; j < region.terms.size(); ++j) {
      auto intersection =
          exact_order_region_intersect_terms(term, region.terms[j]);
      exact_order_region_canonicalize_term(&intersection);
      if (exact_order_region_cell_has_positive_measure(intersection)) {
        ++metrics.overlapping_symbolic_cell_pair_count;
      }
    }
  }
}

inline bool exact_order_region_probability_root(
    ExactVariantBuildState *plan,
    const ExactOutcomeRegionCompileContext &outcome_context,
    const ExactSymbolicTransitionScenario &formula,
    semantic::Index *out_root_id) {
  ExactOrderRegionBuilder builder;
  std::vector<ExactOrderRegionTargetBranch> target_branches;
  if (!exact_order_region_target_branches(
          *plan,
          formula.transition,
          exact_order_region_needs_target_readiness_time(
              outcome_context, formula.transition),
          &builder,
          &target_branches)) {
    throw std::runtime_error("exact order-region target branch lowering failed");
  }
  std::vector<semantic::Index> terms;
  for (const auto &branch : target_branches) {
    ExactOrderRegionExpr non_win;
    if (!exact_order_region_competitor_non_win(
            *plan,
            outcome_context,
            formula,
            branch.expr,
            branch.readiness_time_id,
            &builder,
            &non_win)) {
      throw std::runtime_error("exact order-region competitor non-win lowering failed");
    }
    auto region =
        exact_order_region_simplify(
            exact_order_region_conjoin(branch.expr, non_win));
    region = exact_order_region_minimize_positive_union(std::move(region));
    const ExactProjectionRelationOps projection_ops{
        exact_order_region_factor_context_overlaps_expr,
        exact_order_region_expr_relation_can_collapse,
        exact_order_region_expand_relation_factor};
    ExactOrderRegionExpr planned_metric_region;
    for (const auto &term : region.terms) {
      ExactProjectionPlan projection_plan;
      if (!exact_projection_plan_cell(
              *plan,
              term,
              {},
              {},
              &projection_ops,
              builder,
              &projection_plan)) {
        throw std::runtime_error("exact order-region projection planning failed");
      }
      builder = projection_plan.builder_after;
      exact_projection_collect_metric_cells(
          projection_plan, &planned_metric_region);
    }
    region =
        exact_order_region_minimize_positive_union(
            exact_order_region_simplify(std::move(planned_metric_region)));
    exact_complexity_observe_region(plan, region);
    for (const auto &term : region.terms) {
      semantic::Index term_root{semantic::kInvalidIndex};
      if (!exact_order_region_lower_term_root(
              plan, term, 0, 0, &builder, &projection_ops, &term_root)) {
        throw std::runtime_error("exact order-region term lowering failed");
      }
      terms.push_back(
          compiled_math_root_node_id(plan->compiled_math, term_root));
    }
  }
  const auto node =
      terms.empty()
          ? compiled_math_constant(&plan->compiled_math, 0.0)
          : compiled_math_algebra_node(
                &plan->compiled_math,
                CompiledMathNodeKind::CleanSignedSum,
                std::move(terms),
                CompiledMathValueKind::Scalar);
  *out_root_id = compiled_math_make_root(&plan->compiled_math, node);
  return true;
}

} // namespace detail
} // namespace accumulatr::eval
