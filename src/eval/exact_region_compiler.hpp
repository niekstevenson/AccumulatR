#pragma once

#include "exact_types.hpp"
#include "exact_region.hpp"
#include "exact_compiled_math_lowering.hpp"
#include "exact_projection_planner.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool exact_order_region_expr_before(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_not_before(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_after(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_before_or_at(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_expr_at_time(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out);

inline bool exact_order_region_truth_at_time(
    const ExactVariantPlan &plan,
    const ExactTransitionGuardSet &guards,
    semantic::Index time_id,
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

inline bool exact_order_region_expr_span_conjunction(
    const ExactVariantPlan &plan,
    const ExactIndexSpan span,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out,
    const bool before) {
  ExactOrderRegionExpr combined = exact_order_region_one();
  const auto &program = plan.lowered.program;
  (void)builder;
  for (semantic::Index i = 0; i < span.size; ++i) {
    ExactOrderRegionExpr child;
    const auto child_id =
        program.expr_args[static_cast<std::size_t>(span.offset + i)];
    if (!exact_order_region_expr_relation_factor(
            child_id, time_id, before, !before, &child)) {
      return false;
    }
    combined = exact_order_region_conjoin(combined, child);
  }
  *out = std::move(combined);
  return true;
}

inline bool exact_order_region_expr_span_union(
    const ExactVariantPlan &plan,
    const ExactIndexSpan span,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out,
    const bool before) {
  ExactOrderRegionExpr combined = exact_order_region_zero();
  const auto &program = plan.lowered.program;
  (void)builder;
  bool first = true;
  for (semantic::Index i = 0; i < span.size; ++i) {
    ExactOrderRegionExpr child;
    const auto child_id =
        program.expr_args[static_cast<std::size_t>(span.offset + i)];
    if (!exact_order_region_expr_relation_factor(
            child_id, time_id, before, !before, &child)) {
      return false;
    }
    combined = first ? std::move(child)
                     : exact_order_region_union(combined, child);
    first = false;
  }
  *out = std::move(combined);
  return true;
}

inline bool exact_order_region_expr_span_before_or_at_conjunction(
    const ExactVariantPlan &plan,
    const ExactIndexSpan span,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr combined = exact_order_region_one();
  const auto &program = plan.lowered.program;
  (void)builder;
  for (semantic::Index i = 0; i < span.size; ++i) {
    ExactOrderRegionExpr child;
    const auto child_id =
        program.expr_args[static_cast<std::size_t>(span.offset + i)];
    if (!exact_order_region_expr_relation_factor(
            child_id, time_id, true, true, &child)) {
      return false;
    }
    combined = exact_order_region_conjoin(combined, child);
  }
  *out = std::move(combined);
  return true;
}

inline bool exact_order_region_expr_span_before_or_at_union(
    const ExactVariantPlan &plan,
    const ExactIndexSpan span,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr combined = exact_order_region_zero();
  const auto &program = plan.lowered.program;
  (void)builder;
  bool first = true;
  for (semantic::Index i = 0; i < span.size; ++i) {
    ExactOrderRegionExpr child;
    const auto child_id =
        program.expr_args[static_cast<std::size_t>(span.offset + i)];
    if (!exact_order_region_expr_relation_factor(
            child_id, time_id, true, true, &child)) {
      return false;
    }
    combined = first ? std::move(child)
                     : exact_order_region_union(combined, child);
    first = false;
  }
  *out = std::move(combined);
  return true;
}

inline bool exact_order_region_guard_allowed_at_time(
    const ExactVariantPlan &plan,
    const ExactExprKernel &kernel,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto &program = plan.lowered.program;
  (void)builder;
  ExactOrderRegionExpr allowed;
  if (!exact_order_region_expr_relation_factor(
          kernel.guard_blocker_expr_id,
          time_id,
          false,
          false,
          &allowed)) {
    return false;
  }
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    const auto child =
        program.expr_args[
            static_cast<std::size_t>(kernel.children.offset + i)];
    ExactOrderRegionExpr unless_before;
    if (!exact_order_region_expr_relation_factor(
            child, time_id, true, false, &unless_before)) {
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
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
    const ExactExprKernel &kernel,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  (void)builder;
  ExactOrderRegionExpr blocked;
  if (!exact_order_region_expr_relation_factor(
          kernel.guard_blocker_expr_id,
          time_id,
          true,
          true,
          &blocked)) {
    return false;
  }
  const auto &program = plan.lowered.program;
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    const auto child =
        program.expr_args[
            static_cast<std::size_t>(kernel.children.offset + i)];
    ExactOrderRegionExpr unless_not_before;
    if (!exact_order_region_expr_relation_factor(
            child, time_id, false, true, &unless_not_before)) {
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
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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

inline bool exact_order_region_expr_after(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto &kernel = plan.expr_kernels[pos];
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
    *out = exact_order_region_one();
    return true;
  case semantic::ExprKind::TrueExpr:
    *out = exact_order_region_zero();
    return true;
  case semantic::ExprKind::Event: {
    ExactOrderRegionExpr expr = exact_order_region_one();
    exact_order_region_append_lower(
        &expr.terms.back(), kernel.event_source_id, time_id);
    *out = std::move(expr);
    return true;
  }
  case semantic::ExprKind::And: {
    ExactOrderRegionExpr combined = exact_order_region_zero();
    bool first = true;
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      ExactOrderRegionExpr child_after;
      if (!exact_order_region_expr_relation_factor(
              child, time_id, false, false, &child_after)) {
        return false;
      }
      combined = first ? std::move(child_after)
                       : exact_order_region_union(combined, child_after);
      first = false;
    }
    *out = std::move(combined);
    return true;
  }
  case semantic::ExprKind::Or: {
    ExactOrderRegionExpr combined = exact_order_region_one();
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child =
          program.expr_args[
              static_cast<std::size_t>(kernel.children.offset + i)];
      ExactOrderRegionExpr child_after;
      if (!exact_order_region_expr_relation_factor(
              child, time_id, false, false, &child_after)) {
        return false;
      }
      combined =
          exact_order_region_conjoin(
              std::move(combined),
              std::move(child_after));
    }
    *out = std::move(combined);
    return true;
  }
  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    return exact_order_region_expr_before(
        plan, child, time_id, builder, out);
  }
  case semantic::ExprKind::Guard:
    return exact_order_region_guard_after(
        plan, kernel, time_id, builder, out);
  }
  return false;
}

inline bool exact_order_region_expr_before_or_at(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
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
  case semantic::ExprKind::And:
    return exact_order_region_expr_span_before_or_at_conjunction(
        plan, kernel.children, time_id, builder, out);
  case semantic::ExprKind::Or:
    return exact_order_region_expr_span_before_or_at_union(
        plan, kernel.children, time_id, builder, out);
  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    return exact_order_region_expr_not_before(
        plan, child, time_id, builder, out);
  }
  case semantic::ExprKind::Guard: {
    *out = exact_order_region_zero();
    if (!exact_order_region_guard_before_or_at(
            plan, kernel, time_id, builder, out)) {
      return false;
    }
    return true;
  }
  }
  return false;
}

inline bool exact_order_region_guard_at_time(
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
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
        &expr.terms.back(), kernel.event_source_id, time_id);
    *out = std::move(expr);
    return true;
  }
  case semantic::ExprKind::And:
    return exact_order_region_expr_span_conjunction(
        plan, kernel.children, time_id, builder, out, true);
  case semantic::ExprKind::Or:
    return exact_order_region_expr_span_union(
        plan, kernel.children, time_id, builder, out, true);
  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    return exact_order_region_expr_not_before(
        plan, child, time_id, builder, out);
  }
  case semantic::ExprKind::Guard:
    return exact_order_region_guard_before(
        plan, kernel, expr_id, time_id, builder, out);
  }
  return false;
}

inline bool exact_order_region_expr_not_before(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto &kernel = plan.expr_kernels[pos];
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
    *out = exact_order_region_one();
    return true;
  case semantic::ExprKind::TrueExpr:
    *out = exact_order_region_zero();
    return true;
  case semantic::ExprKind::Event: {
    ExactOrderRegionExpr expr = exact_order_region_one();
    exact_order_region_append_lower(
        &expr.terms.back(), kernel.event_source_id, time_id, true);
    *out = std::move(expr);
    return true;
  }
  case semantic::ExprKind::And:
    return exact_order_region_expr_span_union(
        plan, kernel.children, time_id, builder, out, false);
  case semantic::ExprKind::Or:
    return exact_order_region_expr_span_conjunction(
        plan, kernel.children, time_id, builder, out, false);
  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    return exact_order_region_expr_before(
        plan, child, time_id, builder, out);
  }
  case semantic::ExprKind::Guard:
    return exact_order_region_guard_not_before(
        plan, kernel, expr_id, time_id, builder, out);
  }
  return false;
}

inline bool exact_order_region_expr_at_time(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto &kernel = plan.expr_kernels[pos];
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::TrueExpr:
    *out = exact_order_region_zero();
    return true;
  case semantic::ExprKind::Event: {
    ExactOrderRegionExpr expr = exact_order_region_one();
    exact_order_region_append_exact(
        &expr.terms.back(), kernel.event_source_id, time_id);
    *out = std::move(expr);
    return true;
  }
  case semantic::ExprKind::And: {
    ExactOrderRegionExpr combined = exact_order_region_zero();
    const auto begin = kernel.children.offset;
    const auto end = begin + kernel.children.size;
    for (semantic::Index active = begin; active < end; ++active) {
      const auto active_child =
          program.expr_args[static_cast<std::size_t>(active)];
      ExactOrderRegionExpr branch;
      if (!exact_order_region_expr_at_time(
              plan, active_child, time_id, builder, &branch)) {
        return false;
      }
      for (semantic::Index other = begin; other < end; ++other) {
        if (other == active) {
          continue;
        }
        ExactOrderRegionExpr before_or_at;
        if (!exact_order_region_expr_relation_factor(
                program.expr_args[static_cast<std::size_t>(other)],
                time_id,
                true,
                true,
                &before_or_at)) {
          return false;
        }
        branch = exact_order_region_conjoin(branch, before_or_at);
      }
      combined =
          exact_order_region_union(
              std::move(combined),
              std::move(branch));
    }
    *out = std::move(combined);
    return true;
  }
  case semantic::ExprKind::Or: {
    ExactOrderRegionExpr combined = exact_order_region_zero();
    const auto begin = kernel.children.offset;
    const auto end = begin + kernel.children.size;
    for (semantic::Index active = begin; active < end; ++active) {
      const auto active_child =
          program.expr_args[static_cast<std::size_t>(active)];
      ExactOrderRegionExpr branch;
      if (!exact_order_region_expr_at_time(
              plan, active_child, time_id, builder, &branch)) {
        return false;
      }
      for (semantic::Index other = begin; other < end; ++other) {
        if (other == active) {
          continue;
        }
        ExactOrderRegionExpr not_before;
        if (!exact_order_region_expr_relation_factor(
                program.expr_args[static_cast<std::size_t>(other)],
                time_id,
                false,
                true,
                &not_before)) {
          return false;
        }
        branch = exact_order_region_conjoin(branch, not_before);
      }
      combined =
          exact_order_region_union(
              std::move(combined),
              std::move(branch));
    }
    *out = std::move(combined);
    return true;
  }
  case semantic::ExprKind::Not:
    *out = exact_order_region_zero();
    return true;
  case semantic::ExprKind::Guard:
    return exact_order_region_guard_at_time(
        plan, kernel, time_id, builder, out);
  }
  return false;
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
    const ExactVariantPlan &plan,
    const ExactSymbolicTransitionTime &transition,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto release_source_id =
      exact_symbolic_transition_release_source_id(transition);
  ExactOrderRegionExpr region = exact_order_region_one();
  exact_order_region_append_exact(
      &region.terms.back(), release_source_id, observed_time_id);
  for (const auto &guard :
       transition.readiness_time_expr.requirements.guards) {
    if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
      exact_order_region_append_upper(
          &region.terms.back(), guard.subject_id, observed_time_id);
      continue;
    }
    if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
      ExactOrderRegionExpr before;
      if (!exact_order_region_expr_relation_factor(
              guard.subject_id, observed_time_id, true, false, &before)) {
        return false;
      }
      region = exact_order_region_conjoin(region, before);
      continue;
    }
    return false;
  }
  for (const auto &guard : transition.guards.guards) {
    if (guard.kind == ExactTransitionGuardKind::SourceAfter) {
      exact_order_region_append_lower(
          &region.terms.back(), guard.subject_id, observed_time_id, true);
      continue;
    }
    if (guard.kind == ExactTransitionGuardKind::ExprAfter) {
      ExactOrderRegionExpr not_before;
      if (!exact_order_region_expr_relation_factor(
              guard.subject_id, observed_time_id, false, true, &not_before)) {
        return false;
      }
      region = exact_order_region_conjoin(region, not_before);
      continue;
    }
    return false;
  }
  for (const auto &fact : transition.order_region.source_order_facts) {
    for (auto &term : region.terms) {
      if (!exact_order_region_append_source_order_fact(
              &term, fact, builder)) {
        return false;
      }
    }
  }
  *out = std::move(region);
  return true;
}

inline bool exact_order_region_transition_requirements_at_times(
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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
        if (!exact_order_region_expr_relation_factor(
                guard.subject_id, time_id, true, false, &before)) {
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
        if (!exact_order_region_expr_relation_factor(
                guard.subject_id, time_id, true, false, &before)) {
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
    const ExactVariantPlan &plan,
    const ExactTransitionGuard &guard,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  (void)plan;
  (void)builder;
  ExactOrderRegionExpr region = exact_order_region_one();
  if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
    exact_order_region_append_upper(
        &region.terms.back(), guard.subject_id, time_id, true);
  } else if (guard.kind == ExactTransitionGuardKind::SourceAfter) {
    exact_order_region_append_lower(
        &region.terms.back(), guard.subject_id, time_id);
  } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
    ExactOrderRegionExpr before;
    if (!exact_order_region_expr_relation_factor(
            guard.subject_id, time_id, true, false, &before)) {
      return false;
    }
    *out = std::move(before);
    return true;
  } else if (guard.kind == ExactTransitionGuardKind::ExprAfter) {
    ExactOrderRegionExpr not_before;
    if (!exact_order_region_expr_relation_factor(
            guard.subject_id, time_id, false, true, &not_before)) {
      return false;
    }
    *out = std::move(not_before);
    return true;
  }
  *out = std::move(region);
  return true;
}

inline bool exact_order_region_truth_at_time(
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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
    const ExactRuntimeOutcomeCompileContext &runtime_outcome,
    const ExactSymbolicTransitionTime &transition) {
  if (!exact_symbolic_transition_has_readiness(transition)) {
    return false;
  }
  for (const auto &block : runtime_outcome.competitor_blocks) {
    for (const auto &subset : block.subsets) {
      for (const auto &scenario : subset.scenarios) {
        const auto relation =
            exact_symbolic_transition_relation(transition, scenario.transition);
        if (relation.competitor_can_positively_coincide) {
          return true;
        }
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
    const ExactVariantPlan &plan,
    const ExactTransitionGuard &guard,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  (void)plan;
  (void)builder;
  ExactOrderRegionExpr failure = exact_order_region_one();
  if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
    exact_order_region_append_lower(
        &failure.terms.back(), guard.subject_id, time_id, true);
  } else if (guard.kind == ExactTransitionGuardKind::SourceAfter) {
    exact_order_region_append_upper(
        &failure.terms.back(), guard.subject_id, time_id, true);
  } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
    if (!exact_order_region_expr_relation_factor(
            guard.subject_id, time_id, false, true, &failure)) {
      return false;
    }
  } else if (guard.kind == ExactTransitionGuardKind::ExprAfter) {
    if (!exact_order_region_expr_relation_factor(
            guard.subject_id, time_id, true, true, &failure)) {
      return false;
    }
  }
  *out = std::move(failure);
  return true;
}

inline bool exact_order_region_truth_not_at_time(
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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

inline bool exact_order_region_transition_no_strict_precedence(
    const ExactVariantPlan &plan,
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

  ExactOrderRegionExpr valid_before;
  if (!exact_order_region_transition_requirements_at_times(
          plan,
          competitor,
          competitor.guards,
          competitor_time_id,
          competitor_time_id,
          competitor_time_id,
          builder,
          &valid_before)) {
    return false;
  }
  for (auto &term : valid_before.terms) {
    exact_order_region_append_time_order(
        &term, competitor_time_id, observed_time_id);
  }

  *out = exact_order_region_union(
      release_not_before,
      exact_order_region_subtract_region(
          std::move(release_before),
          valid_before));
  return true;
}

inline bool exact_order_region_transition_no_coincident_precedence(
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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

inline bool exact_order_region_competitor_subset_non_win(
    const ExactVariantPlan &plan,
    const ExactSymbolicTransitionScenario &target,
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr non_win = exact_order_region_one();
  for (const auto &scenario : subset.scenarios) {
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
  if (outcome_indices.empty()) {
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

inline bool exact_order_region_competitor_block_non_win(
    const ExactVariantPlan &plan,
    const ExactSymbolicTransitionScenario &target,
    const ExactRuntimeCompetitorBlockPlan &block,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr block_non_win = exact_order_region_one();
  for (const auto &subset : block.subsets) {
    if (subset.inclusion_sign < 0 ||
        subset.expr_roots.size() != 1U) {
      continue;
    }
    ExactOrderRegionExpr subset_non_win;
    if (!exact_order_region_competitor_subset_non_win(
            plan,
            target,
            subset,
            target_readiness_time_id,
            builder,
            &subset_non_win)) {
      return false;
    }
    block_non_win =
        exact_order_region_conjoin(
            block_non_win,
            exact_order_region_gate_non_win(
                std::move(subset_non_win),
                subset.outcome_indices));
  }
  *out = std::move(block_non_win);
  return true;
}

inline bool exact_order_region_competitor_non_win(
    const ExactVariantPlan &plan,
    const ExactRuntimeOutcomeCompileContext &runtime_outcome,
    const ExactSymbolicTransitionScenario &target,
    const semantic::Index target_readiness_time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  ExactOrderRegionExpr non_win = exact_order_region_one();
  for (const auto &block : runtime_outcome.competitor_blocks) {
    ExactOrderRegionExpr block_non_win;
    if (!exact_order_region_competitor_block_non_win(
            plan,
            target,
            block,
            target_readiness_time_id,
            builder,
            &block_non_win)) {
      return false;
    }
    non_win = exact_order_region_conjoin(non_win, block_non_win);
  }
  *out = std::move(non_win);
  return true;
}

inline bool exact_order_region_factor_context_overlaps_expr(
    const ExactVariantPlan &plan,
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

inline bool exact_order_region_find_materialized_relation(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    ExactOrderRegionExprValueFactor *out) {
  for (const auto &factor : exact_region_expr_atoms(term)) {
    if (factor.density) {
      continue;
    }
    if (exact_order_region_factor_context_overlaps_expr(
            plan, term, factor.expr_id)) {
      *out = factor;
      return true;
    }
  }
  return false;
}

inline bool exact_order_region_expand_relation_factor(
    const ExactVariantPlan &plan,
    const ExactOrderRegionExprValueFactor &factor,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (factor.density) {
    return false;
  }
  if (factor.before && factor.inclusive) {
    return exact_order_region_expr_before_or_at(
        plan, factor.expr_id, factor.time_id, builder, out);
  }
  if (factor.before) {
    return exact_order_region_expr_before(
        plan, factor.expr_id, factor.time_id, builder, out);
  }
  if (factor.inclusive) {
    return exact_order_region_expr_not_before(
        plan, factor.expr_id, factor.time_id, builder, out);
  }
  return exact_order_region_expr_after(
      plan, factor.expr_id, factor.time_id, builder, out);
}

inline bool exact_order_region_materialize_coupled_relations(
    const ExactVariantPlan &plan,
    ExactOrderRegionExpr expr,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  expr = exact_order_region_simplify(std::move(expr));
  while (true) {
    bool changed = false;
    ExactOrderRegionExpr next = exact_order_region_zero();
    for (auto term : expr.terms) {
      exact_order_region_canonicalize_term(&term);
      if (term.impossible || term.sign == 0.0) {
        continue;
      }
      ExactOrderRegionExprValueFactor factor;
      if (!exact_order_region_find_materialized_relation(plan, term, &factor)) {
        next.terms.push_back(std::move(term));
        continue;
      }
      changed = true;
      ExactRegionCell residual = std::move(term);
      if (!exact_order_region_remove_expr_value_atom(&residual, factor)) {
        return false;
      }
      exact_order_region_canonicalize_term(&residual);
      if (residual.impossible || residual.sign == 0.0) {
        continue;
      }
      ExactOrderRegionExpr residual_expr;
      residual_expr.terms.push_back(std::move(residual));
      ExactOrderRegionExpr relation;
      if (!exact_order_region_expand_relation_factor(
              plan, factor, builder, &relation)) {
        return false;
      }
      auto materialized =
          exact_order_region_simplify(
              exact_order_region_conjoin(residual_expr, relation));
      exact_order_region_append_expr(&next, std::move(materialized));
    }
    next = exact_order_region_simplify(std::move(next));
    if (!changed) {
      *out = std::move(next);
      return true;
    }
    expr = std::move(next);
  }
}

inline void exact_complexity_observe_region(
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
    const ExactRuntimeOutcomeCompileContext &runtime_outcome,
    const ExactSymbolicTransitionScenario &formula,
    semantic::Index *out_root_id) {
  ExactOrderRegionBuilder builder;
  std::vector<ExactOrderRegionTargetBranch> target_branches;
  if (!exact_order_region_target_branches(
          *plan,
          formula.transition,
          exact_order_region_needs_target_readiness_time(
              runtime_outcome, formula.transition),
          &builder,
          &target_branches)) {
    throw std::runtime_error("exact order-region target branch lowering failed");
  }
  std::vector<semantic::Index> terms;
  for (const auto &branch : target_branches) {
    ExactOrderRegionExpr non_win;
    if (!exact_order_region_competitor_non_win(
            *plan,
            runtime_outcome,
            formula,
            branch.readiness_time_id,
            &builder,
            &non_win)) {
      throw std::runtime_error("exact order-region competitor non-win lowering failed");
    }
    auto region =
        exact_order_region_simplify(
            exact_order_region_conjoin(branch.expr, non_win));
    if (!exact_order_region_materialize_coupled_relations(
            *plan, std::move(region), &builder, &region)) {
      throw std::runtime_error(
          "exact order-region relation materialization failed");
    }
    region = exact_order_region_minimize_positive_union(std::move(region));
    exact_complexity_observe_region(plan, region);
    for (const auto &term : region.terms) {
      semantic::Index term_root{semantic::kInvalidIndex};
      if (!exact_order_region_lower_term_root(
              plan, term, 0, &builder, &term_root)) {
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
