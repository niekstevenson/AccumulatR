#pragma once

#include <algorithm>
#include <limits>

#include "exact_oracle.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool relation_view_with_overlay(const RelationView &base,
                                       const ExactRelationTemplate &overlay,
                                       RelationView *out) {
  for (std::size_t i = 0; i < overlay.source_ids.size(); ++i) {
    const auto source_id = overlay.source_ids[i];
    const auto relation = overlay.relations[i];
    const auto inherited = base.relation_for(source_id);
    if (inherited != ExactRelation::Unknown && inherited != relation) {
      return false;
    }
  }
  if (out != nullptr) {
    *out = overlay.empty() ? base : base.with_overlay(&overlay);
  }
  return true;
}

struct ExactScenarioRuntimeView {
  const ExactTransitionScenario &scenario;
  const semantic::Index *before_sources{nullptr};
  const semantic::Index *after_sources{nullptr};
  const semantic::Index *ready_exprs{nullptr};
  const semantic::Index *tail_exprs{nullptr};
  semantic::Index before_source_count{0};
  semantic::Index after_source_count{0};
  semantic::Index ready_expr_count{0};
  semantic::Index tail_expr_count{0};

  [[nodiscard]] bool has_readiness() const noexcept {
    return before_source_count > 0 || ready_expr_count > 0;
  }
};

inline const semantic::Index *scenario_span_data(
    const std::vector<semantic::Index> &arena,
    const ExactIndexSpan span) {
  if (span.empty()) {
    return nullptr;
  }
  return arena.data() + static_cast<std::size_t>(span.offset);
}

inline ExactScenarioRuntimeView make_exact_scenario_runtime_view(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario) {
  return ExactScenarioRuntimeView{
      scenario,
      scenario_span_data(plan.scenario_source_ids, scenario.before_source_span),
      scenario_span_data(plan.scenario_source_ids, scenario.after_source_span),
      scenario_span_data(plan.scenario_expr_ids, scenario.ready_expr_span),
      scenario_span_data(plan.scenario_expr_ids, scenario.tail_expr_span),
      scenario.before_source_span.size,
      scenario.after_source_span.size,
      scenario.ready_expr_span.size,
      scenario.tail_expr_span.size};
}

class ForcedExprEvaluator {
public:
  explicit ForcedExprEvaluator(const ExactVariantPlan &plan)
      : plan_(plan),
        program_(plan.lowered.program),
        cdf_cache_(program_.expr_kind.size(), 0.0),
        cdf_time_(program_.expr_kind.size(), 0.0),
        cdf_epoch_(program_.expr_kind.size(), 0U),
        density_cache_(program_.expr_kind.size(), 0.0),
        density_time_(program_.expr_kind.size(), 0.0),
        density_epoch_(program_.expr_kind.size(), 0U) {}

  void reset(ExactSourceOracle *oracle, const RelationView relation_view = {}) {
    oracle_ = oracle;
    relation_view_ = relation_view;
    invalidate_cache();
  }

  void invalidate_cache() {
    ++epoch_;
    if (epoch_ == 0U) {
      epoch_ = 1U;
      std::fill(cdf_epoch_.begin(), cdf_epoch_.end(), 0U);
      std::fill(density_epoch_.begin(), density_epoch_.end(), 0U);
    }
  }

  double expr_cdf(const semantic::Index expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    const double current_time = oracle_->conditional_time();
    const bool ignoring_guard_upper =
        ignored_guard_upper_ref_ != semantic::kInvalidIndex;
    const bool cacheable =
        ignored_expr_upper_ != expr_idx && !ignoring_guard_upper;
    if (cacheable && cdf_epoch_[pos] == epoch_ &&
        cdf_time_[pos] == current_time) {
      return cdf_cache_[pos];
    }
    if (cacheable) {
      if (const double *upper = oracle_->expr_upper_bound_for(expr_idx)) {
        const double normalizer = oracle_->expr_upper_bound_normalizer(expr_idx);
        if (!(normalizer > 0.0)) {
          return 0.0;
        }
        if (current_time >= *upper) {
          return 1.0;
        }
        const auto previous_ignore = ignored_expr_upper_;
        ignored_expr_upper_ = expr_idx;
        const double raw = expr_cdf(expr_idx);
        ignored_expr_upper_ = previous_ignore;
        return clamp_probability(raw / normalizer);
      }
    }
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    double value = 0.0;
    switch (kind) {
    case semantic::ExprKind::Impossible:
      value = 0.0;
      break;
    case semantic::ExprKind::TrueExpr:
      value = 1.0;
      break;
    case semantic::ExprKind::Event:
      value = source_channels(program_.expr_source_ids[pos]).cdf;
      break;
    case semantic::ExprKind::And: {
      double cdf = 1.0;
      for (semantic::Index i =
               program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        cdf *= expr_cdf(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      value = clamp_probability(cdf);
      break;
    }
    case semantic::ExprKind::Or: {
      std::vector<semantic::Index> children;
      for (semantic::Index i =
               program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        children.push_back(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      semantic::Index absorbed_source_id{semantic::kInvalidIndex};
      if (or_absorbed_event_source(children, &absorbed_source_id)) {
        value = source_cdf(absorbed_source_id);
      } else if ((oracle_->has_guard_upper_bound_overlay() ||
           oracle_->has_source_order_overlay()) &&
          expr_children_overlap(children)) {
        value = expr_union_cdf(children);
      } else {
        double surv = 1.0;
        for (const auto child : children) {
          surv *= expr_survival(child);
        }
        value = clamp_probability(1.0 - surv);
      }
      break;
    }
    case semantic::ExprKind::Guard: {
      if (!(current_time > 0.0)) {
        value = 0.0;
        break;
      }
      const auto ref =
          program_.expr_ref_child[static_cast<std::size_t>(expr_idx)];
      const auto blocker =
          program_.expr_blocker_child[static_cast<std::size_t>(expr_idx)];
      semantic::Index ref_source_id{semantic::kInvalidIndex};
      semantic::Index blocker_source_id{semantic::kInvalidIndex};
      const bool simple_guard =
          simple_guard_source_ids(ref, blocker, &ref_source_id, &blocker_source_id);
      if (simple_guard &&
          oracle_->source_order_known_before(blocker_source_id, ref_source_id)) {
        value = 0.0;
        break;
      }
      if (simple_guard &&
          !guard_upper_bound_ignored(ref_source_id, blocker_source_id)) {
        if (const auto *upper = oracle_->guard_upper_bound_for(
                ref_source_id, blocker_source_id)) {
          if (!(upper->normalizer > 0.0)) {
            value = 0.0;
            break;
          }
          if (current_time >= upper->time) {
            value = 1.0;
            break;
          }
          const auto previous_ref = ignored_guard_upper_ref_;
          const auto previous_blocker = ignored_guard_upper_blocker_;
          ignored_guard_upper_ref_ = ref_source_id;
          ignored_guard_upper_blocker_ = blocker_source_id;
          const double raw = expr_cdf(expr_idx);
          ignored_guard_upper_ref_ = previous_ref;
          ignored_guard_upper_blocker_ = previous_blocker;
          value = clamp_probability(raw / upper->normalizer);
          break;
        }
      }
      const bool has_unless =
          program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)] !=
          program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
      const auto ref_kind = static_cast<semantic::ExprKind>(
          program_.expr_kind[static_cast<std::size_t>(ref)]);
      const auto blocker_kind = static_cast<semantic::ExprKind>(
          program_.expr_kind[static_cast<std::size_t>(blocker)]);
      if (!has_unless && ref_kind == semantic::ExprKind::Event) {
        const auto ref_source_id =
            program_.expr_source_ids[static_cast<std::size_t>(ref)];
        if (const double *ref_time = oracle_->exact_time_for_source(ref_source_id)) {
          if (!(current_time >= *ref_time)) {
            value = 0.0;
          } else {
            const auto guard = oracle_->conditional_time_guard(*ref_time);
            value = expr_survival(blocker);
          }
          break;
        }
      }
      if (!has_unless && ref_kind == semantic::ExprKind::And) {
        std::vector<semantic::Index> remaining_ref_children;
        double exact_ref_time = -std::numeric_limits<double>::infinity();
        bool has_exact_ref_child = false;
        for (semantic::Index i =
                 program_.expr_arg_offsets[static_cast<std::size_t>(ref)];
             i < program_.expr_arg_offsets[static_cast<std::size_t>(ref + 1)];
             ++i) {
          const auto child = program_.expr_args[static_cast<std::size_t>(i)];
          const auto child_kind = static_cast<semantic::ExprKind>(
              program_.expr_kind[static_cast<std::size_t>(child)]);
          if (child_kind == semantic::ExprKind::Event) {
            const auto child_source_id =
                program_.expr_source_ids[static_cast<std::size_t>(child)];
            if (const double *child_time =
                    oracle_->exact_time_for_source(child_source_id)) {
              exact_ref_time = std::max(exact_ref_time, *child_time);
              has_exact_ref_child = true;
              continue;
            }
          }
          remaining_ref_children.push_back(child);
        }
        if (has_exact_ref_child) {
          if (!(current_time >= exact_ref_time)) {
            value = 0.0;
            break;
          }
          {
            const auto guard = oracle_->conditional_time_guard(exact_ref_time);
            const double ready_at_exact =
                remaining_ref_children.empty()
                    ? 1.0
                    : conjunction_cdf(remaining_ref_children);
            value = ready_at_exact * expr_survival(blocker);
          }
          if (!remaining_ref_children.empty() &&
              current_time > exact_ref_time) {
            value += quadrature::integrate_finite_default(
                [&](const double u) {
                  const auto guard = oracle_->conditional_time_guard(u);
                  double term = conjunction_density(remaining_ref_children);
                  if (!std::isfinite(term) || term == 0.0) {
                    return 0.0;
                  }
                  term *= expr_survival(blocker);
                  return std::isfinite(term) ? term : 0.0;
                },
                exact_ref_time,
                current_time);
          }
          value = clamp_probability(value);
          break;
        }

        bool ref_has_current_exact_child = false;
        for (semantic::Index i =
                 program_.expr_arg_offsets[static_cast<std::size_t>(ref)];
             i < program_.expr_arg_offsets[static_cast<std::size_t>(ref + 1)];
             ++i) {
          const auto child = program_.expr_args[static_cast<std::size_t>(i)];
          const auto child_kind = static_cast<semantic::ExprKind>(
              program_.expr_kind[static_cast<std::size_t>(child)]);
          if (child_kind != semantic::ExprKind::Event) {
            continue;
          }
          const auto child_source_id =
              program_.expr_source_ids[static_cast<std::size_t>(child)];
          const double *child_time =
              oracle_->exact_time_for_source(child_source_id);
          if (child_time != nullptr &&
              std::fabs(*child_time - current_time) <= 1e-12) {
            ref_has_current_exact_child = true;
            break;
          }
        }
        if (ref_has_current_exact_child) {
          value = clamp_probability(expr_cdf(ref) * expr_survival(blocker));
          break;
        }
      }
      if (!has_unless && blocker_kind == semantic::ExprKind::Event) {
        const auto blocker_source_id =
            program_.expr_source_ids[static_cast<std::size_t>(blocker)];
        if (const double *blocker_time =
                oracle_->exact_time_for_source(blocker_source_id)) {
          const double upper = std::min(current_time, *blocker_time);
          if (!(upper > 0.0)) {
            value = 0.0;
          } else {
            const auto guard = oracle_->conditional_time_guard(upper);
            value = expr_cdf(ref);
          }
          break;
        }
      }
      value = clamp_probability(quadrature::integrate_finite_default(
          [&](const double u) {
            const auto guard = oracle_->conditional_time_guard(u);
            return expr_density(expr_idx);
          },
          0.0,
          current_time));
      break;
    }
    case semantic::ExprKind::Not:
      value = clamp_probability(
          1.0 - expr_cdf(program_.expr_args[static_cast<std::size_t>(
                    program_.expr_arg_offsets[pos])]));
      break;
    }
    if (cacheable) {
      cdf_epoch_[pos] = epoch_;
      cdf_time_[pos] = current_time;
      cdf_cache_[pos] = value;
    }
    return value;
  }

  double expr_survival(const semantic::Index expr_idx) {
    return clamp_probability(1.0 - expr_cdf(expr_idx));
  }

  double source_cdf(const semantic::Index source_id) {
    return clamp_probability(source_channels(source_id).cdf);
  }

  double source_pdf(const semantic::Index source_id) {
    return safe_density(source_channels(source_id).pdf);
  }

  double source_survival(const semantic::Index source_id) {
    return clamp_probability(source_channels(source_id).survival);
  }

  ExactSourceOracle *oracle() const {
    return oracle_;
  }

  const ExactVariantPlan &plan() const {
    return plan_;
  }

  const RelationView &relation_view() const {
    return relation_view_;
  }

  double expr_density(const semantic::Index expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    const double current_time = oracle_->conditional_time();
    const bool ignoring_guard_upper =
        ignored_guard_upper_ref_ != semantic::kInvalidIndex;
    const bool cacheable =
        ignored_expr_upper_ != expr_idx && !ignoring_guard_upper;
    if (cacheable && density_epoch_[pos] == epoch_ &&
        density_time_[pos] == current_time) {
      return density_cache_[pos];
    }
    if (cacheable) {
      if (const double *upper = oracle_->expr_upper_bound_for(expr_idx)) {
        const double normalizer = oracle_->expr_upper_bound_normalizer(expr_idx);
        if (!(normalizer > 0.0) || current_time >= *upper) {
          return 0.0;
        }
        const auto previous_ignore = ignored_expr_upper_;
        ignored_expr_upper_ = expr_idx;
        const double raw = expr_density(expr_idx);
        ignored_expr_upper_ = previous_ignore;
        return safe_density(raw / normalizer);
      }
    }
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    double value = 0.0;
    switch (kind) {
    case semantic::ExprKind::Impossible:
    case semantic::ExprKind::TrueExpr:
      value = 0.0;
      break;
    case semantic::ExprKind::Event:
      value = source_channels(program_.expr_source_ids[pos]).pdf;
      break;
    case semantic::ExprKind::And: {
      double total = 0.0;
      for (semantic::Index i =
               program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        const auto child = program_.expr_args[static_cast<std::size_t>(i)];
        double term = expr_density(child);
        if (!std::isfinite(term) || term == 0.0) {
          continue;
        }
        for (semantic::Index j =
                 program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
             j <
             program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
             ++j) {
          if (i == j) {
            continue;
          }
          term *= expr_cdf(program_.expr_args[static_cast<std::size_t>(j)]);
          if (!std::isfinite(term) || term == 0.0) {
            break;
          }
        }
        total += term;
      }
      value = clean_signed_value(total);
      break;
    }
    case semantic::ExprKind::Or: {
      std::vector<semantic::Index> children;
      for (semantic::Index i =
               program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        children.push_back(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      semantic::Index absorbed_source_id{semantic::kInvalidIndex};
      if (or_absorbed_event_source(children, &absorbed_source_id)) {
        value = source_pdf(absorbed_source_id);
      } else if ((oracle_->has_guard_upper_bound_overlay() ||
           oracle_->has_source_order_overlay()) &&
          expr_children_overlap(children)) {
        value = expr_union_density(children);
      } else {
        double total = 0.0;
        for (std::size_t i = 0; i < children.size(); ++i) {
          double term = expr_density(children[i]);
          if (!std::isfinite(term) || term == 0.0) {
            continue;
          }
          for (std::size_t j = 0; j < children.size(); ++j) {
            if (i == j) {
              continue;
            }
            term *= expr_survival(children[j]);
            if (!std::isfinite(term) || term == 0.0) {
              break;
            }
          }
          total += term;
        }
        value = clean_signed_value(total);
      }
      break;
    }
    case semantic::ExprKind::Guard: {
      const auto ref = program_.expr_ref_child[static_cast<std::size_t>(expr_idx)];
      const auto blocker =
          program_.expr_blocker_child[static_cast<std::size_t>(expr_idx)];
      semantic::Index ref_source_id{semantic::kInvalidIndex};
      semantic::Index blocker_source_id{semantic::kInvalidIndex};
      const bool simple_guard =
          simple_guard_source_ids(ref, blocker, &ref_source_id, &blocker_source_id);
      if (simple_guard &&
          oracle_->source_order_known_before(blocker_source_id, ref_source_id)) {
        value = 0.0;
        break;
      }
      if (simple_guard &&
          !guard_upper_bound_ignored(ref_source_id, blocker_source_id)) {
        if (const auto *upper = oracle_->guard_upper_bound_for(
                ref_source_id, blocker_source_id)) {
          if (!(upper->normalizer > 0.0) || current_time >= upper->time) {
            value = 0.0;
            break;
          }
          const auto previous_ref = ignored_guard_upper_ref_;
          const auto previous_blocker = ignored_guard_upper_blocker_;
          ignored_guard_upper_ref_ = ref_source_id;
          ignored_guard_upper_blocker_ = blocker_source_id;
          const double raw = expr_density(expr_idx);
          ignored_guard_upper_ref_ = previous_ref;
          ignored_guard_upper_blocker_ = previous_blocker;
          value = safe_density(raw / upper->normalizer);
          break;
        }
      }
      double blocker_survival = 0.0;
      if (program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)] ==
          program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)]) {
        blocker_survival = expr_survival(blocker);
      } else {
        blocker_survival = clamp_probability(
            1.0 - quadrature::integrate_finite_default(
                      [&](const double u) {
                        const auto guard = oracle_->conditional_time_guard(u);
                        double term = expr_density(blocker);
                        if (!std::isfinite(term) || term == 0.0) {
                          return 0.0;
                        }
                        for (semantic::Index i = program_.expr_arg_offsets
                                                    [static_cast<std::size_t>(
                                                        expr_idx)];
                             i < program_.expr_arg_offsets[static_cast<std::size_t>(
                                      expr_idx + 1)];
                             ++i) {
                          term *= expr_survival(
                              program_.expr_args[static_cast<std::size_t>(i)]);
                          if (!std::isfinite(term) || term == 0.0) {
                            return 0.0;
                          }
                        }
                        return term;
                      },
                      0.0,
                      current_time));
      }
      value = clean_signed_value(expr_density(ref) * blocker_survival);
      break;
    }
    case semantic::ExprKind::Not:
      value = clean_signed_value(
          -expr_density(program_.expr_args[static_cast<std::size_t>(
               program_.expr_arg_offsets[pos])]));
      break;
    }
    if (cacheable) {
      density_epoch_[pos] = epoch_;
      density_time_[pos] = current_time;
      density_cache_[pos] = value;
    }
    return value;
  }

  double conjunction_density(const std::vector<semantic::Index> &exprs) {
    if (exprs.empty()) {
      return 0.0;
    }
    if (exprs.size() == 1U) {
      return expr_density(exprs.front());
    }
    double total = 0.0;
    for (std::size_t i = 0; i < exprs.size(); ++i) {
      double term = expr_density(exprs[i]);
      if (!std::isfinite(term) || term == 0.0) {
        continue;
      }
      for (std::size_t j = 0; j < exprs.size(); ++j) {
        if (i == j) {
          continue;
        }
        term *= expr_cdf(exprs[j]);
        if (!std::isfinite(term) || term == 0.0) {
          break;
        }
      }
      total += term;
    }
    return clean_signed_value(total);
  }

  double conjunction_cdf(const std::vector<semantic::Index> &exprs) {
    if (exprs.empty()) {
      return 0.0;
    }
    if (exprs.size() == 1U) {
      return expr_cdf(exprs.front());
    }
    const double current_time = oracle_->conditional_time();
    if (!(current_time > 0.0)) {
      return 0.0;
    }
    return clamp_probability(quadrature::integrate_finite_default(
        [&](const double u) {
          const auto guard = oracle_->conditional_time_guard(u);
          return conjunction_density(exprs);
        },
        0.0,
        current_time));
  }

private:
  bool simple_guard_source_ids(const semantic::Index ref,
                               const semantic::Index blocker,
                               semantic::Index *ref_source_id,
                               semantic::Index *blocker_source_id) const {
    const auto ref_kind = static_cast<semantic::ExprKind>(
        program_.expr_kind[static_cast<std::size_t>(ref)]);
    const auto blocker_kind = static_cast<semantic::ExprKind>(
        program_.expr_kind[static_cast<std::size_t>(blocker)]);
    if (ref_kind != semantic::ExprKind::Event ||
        blocker_kind != semantic::ExprKind::Event) {
      return false;
    }
    if (ref_source_id != nullptr) {
      *ref_source_id = program_.expr_source_ids[static_cast<std::size_t>(ref)];
    }
    if (blocker_source_id != nullptr) {
      *blocker_source_id =
          program_.expr_source_ids[static_cast<std::size_t>(blocker)];
    }
    return true;
  }

  bool guard_upper_bound_ignored(const semantic::Index ref_source_id,
                                 const semantic::Index blocker_source_id) const {
    return ignored_guard_upper_ref_ == ref_source_id &&
           ignored_guard_upper_blocker_ == blocker_source_id;
  }

  bool expr_equivalent_event_source(const semantic::Index expr_idx,
                                    semantic::Index *source_id) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    if (kind == semantic::ExprKind::Event) {
      if (source_id != nullptr) {
        *source_id = program_.expr_source_ids[pos];
      }
      return true;
    }
    if (kind != semantic::ExprKind::And) {
      return false;
    }

    semantic::Index candidate{semantic::kInvalidIndex};
    bool found_candidate = false;
    const auto begin = program_.expr_arg_offsets[pos];
    const auto end = program_.expr_arg_offsets[pos + 1U];
    for (semantic::Index i = begin; i < end; ++i) {
      const auto child = program_.expr_args[static_cast<std::size_t>(i)];
      semantic::Index child_source{semantic::kInvalidIndex};
      if (expr_equivalent_event_source(child, &child_source)) {
        const double child_cdf = expr_cdf(child);
        if (child_cdf < 1.0 - 1e-10) {
          if (found_candidate) {
            return false;
          }
          candidate = child_source;
          found_candidate = true;
          continue;
        }
      }
      if (expr_cdf(child) < 1.0 - 1e-10) {
        return false;
      }
    }
    if (!found_candidate) {
      return false;
    }
    if (source_id != nullptr) {
      *source_id = candidate;
    }
    return true;
  }

  bool expr_subset_of_event_source(const semantic::Index expr_idx,
                                   const semantic::Index source_id) const {
    if (source_id == semantic::kInvalidIndex) {
      return false;
    }
    const auto pos = static_cast<std::size_t>(expr_idx);
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    switch (kind) {
    case semantic::ExprKind::Event:
      return program_.expr_source_ids[pos] == source_id;
    case semantic::ExprKind::And: {
      const auto begin = program_.expr_arg_offsets[pos];
      const auto end = program_.expr_arg_offsets[pos + 1U];
      for (semantic::Index i = begin; i < end; ++i) {
        if (expr_subset_of_event_source(
                program_.expr_args[static_cast<std::size_t>(i)], source_id)) {
          return true;
        }
      }
      return false;
    }
    case semantic::ExprKind::Or: {
      const auto begin = program_.expr_arg_offsets[pos];
      const auto end = program_.expr_arg_offsets[pos + 1U];
      if (begin == end) {
        return false;
      }
      for (semantic::Index i = begin; i < end; ++i) {
        if (!expr_subset_of_event_source(
                program_.expr_args[static_cast<std::size_t>(i)], source_id)) {
          return false;
        }
      }
      return true;
    }
    case semantic::ExprKind::Guard: {
      const auto ref = program_.expr_ref_child[pos];
      return expr_subset_of_event_source(ref, source_id);
    }
    case semantic::ExprKind::Impossible:
      return true;
    case semantic::ExprKind::TrueExpr:
    case semantic::ExprKind::Not:
      return false;
    }
    return false;
  }

  bool or_absorbed_event_source(const std::vector<semantic::Index> &children,
                                semantic::Index *source_id) {
    for (const auto child : children) {
      semantic::Index candidate{semantic::kInvalidIndex};
      if (!expr_equivalent_event_source(child, &candidate)) {
        continue;
      }
      bool absorbs_all = true;
      for (const auto other : children) {
        if (!expr_subset_of_event_source(other, candidate)) {
          absorbs_all = false;
          break;
        }
      }
      if (!absorbs_all) {
        continue;
      }
      if (source_id != nullptr) {
        *source_id = candidate;
      }
      return true;
    }
    return false;
  }

  bool expr_children_overlap(
      const std::vector<semantic::Index> &children) const {
    for (std::size_t i = 0; i < children.size(); ++i) {
      const auto &lhs =
          plan_.expr_supports[static_cast<std::size_t>(children[i])];
      for (std::size_t j = i + 1U; j < children.size(); ++j) {
        const auto &rhs =
            plan_.expr_supports[static_cast<std::size_t>(children[j])];
        if (supports_overlap(lhs, rhs)) {
          return true;
        }
      }
    }
    return false;
  }

  double expr_union_cdf(const std::vector<semantic::Index> &children) {
    if (children.empty()) {
      return 0.0;
    }
    if (children.size() == 1U) {
      return expr_cdf(children.front());
    }
    double total = 0.0;
    const std::size_t max_mask = std::size_t{1} << children.size();
    for (std::size_t mask = 1U; mask < max_mask; ++mask) {
      std::vector<semantic::Index> subset;
      subset.reserve(children.size());
      for (std::size_t bit = 0; bit < children.size(); ++bit) {
        if ((mask & (std::size_t{1} << bit)) != 0U) {
          subset.push_back(children[bit]);
        }
      }
      const double value = subset.size() == 1U
                               ? expr_cdf(subset.front())
                               : conjunction_cdf(subset);
      total += (subset.size() % 2U == 1U) ? value : -value;
    }
    return clamp_probability(total);
  }

  double expr_union_density(const std::vector<semantic::Index> &children) {
    if (children.empty()) {
      return 0.0;
    }
    if (children.size() == 1U) {
      return expr_density(children.front());
    }
    double total = 0.0;
    const std::size_t max_mask = std::size_t{1} << children.size();
    for (std::size_t mask = 1U; mask < max_mask; ++mask) {
      std::vector<semantic::Index> subset;
      subset.reserve(children.size());
      for (std::size_t bit = 0; bit < children.size(); ++bit) {
        if ((mask & (std::size_t{1} << bit)) != 0U) {
          subset.push_back(children[bit]);
        }
      }
      const double value = subset.size() == 1U
                               ? expr_density(subset.front())
                               : conjunction_density(subset);
      total += (subset.size() % 2U == 1U) ? value : -value;
    }
    return clean_signed_value(total);
  }

  leaf::EventChannels source_channels(const semantic::Index source_id) const {
    if (oracle_->has_source_condition_overlay(source_id)) {
      return oracle_->conditional_source(source_id);
    }
    const auto relation = relation_view_.relation_for(source_id);
    if (relation != ExactRelation::Unknown) {
      return forced_channels(relation);
    }
    return oracle_->conditional_source(source_id);
  }

  const ExactVariantPlan &plan_;
  const runtime::ExactProgram &program_;
  ExactSourceOracle *oracle_{nullptr};
  RelationView relation_view_;
  std::vector<double> cdf_cache_;
  std::vector<double> cdf_time_;
  std::vector<std::uint32_t> cdf_epoch_;
  std::vector<double> density_cache_;
  std::vector<double> density_time_;
  std::vector<std::uint32_t> density_epoch_;
  std::uint32_t epoch_{1U};
  semantic::Index ignored_expr_upper_{semantic::kInvalidIndex};
  semantic::Index ignored_guard_upper_ref_{semantic::kInvalidIndex};
  semantic::Index ignored_guard_upper_blocker_{semantic::kInvalidIndex};
};

struct ForcedExprWorkspace {
  explicit ForcedExprWorkspace(const ExactVariantPlan &plan) : evaluator(plan) {}

  void reset(ExactSourceOracle *oracle, const RelationView relation_view = {}) {
    evaluator.reset(oracle, relation_view);
  }

  ForcedExprEvaluator evaluator;
};

struct ExactRuntimeTermCondition {
  semantic::Index exact_source_id{semantic::kInvalidIndex};
  std::vector<semantic::Index> upper_bound_source_ids;
  std::vector<semantic::Index> lower_bound_source_ids;
  std::vector<semantic::Index> upper_bound_expr_ids;
  std::vector<ExactSourceOrderFact> source_order_facts;
  std::vector<ExactGuardUpperBoundFact> guard_upper_bound_facts;
};

inline void append_runtime_condition_index(
    std::vector<semantic::Index> *sources,
    const semantic::Index index) {
  if (index == semantic::kInvalidIndex) {
    return;
  }
  if (std::find(sources->begin(), sources->end(), index) ==
      sources->end()) {
    sources->push_back(index);
  }
}

inline void append_runtime_source_order_fact(
    std::vector<ExactSourceOrderFact> *facts,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  if (before_source_id == semantic::kInvalidIndex ||
      after_source_id == semantic::kInvalidIndex ||
      before_source_id == after_source_id) {
    return;
  }
  for (const auto &fact : *facts) {
    if (fact.before_source_id == before_source_id &&
        fact.after_source_id == after_source_id) {
      return;
    }
  }
  facts->push_back(ExactSourceOrderFact{before_source_id, after_source_id});
}

inline void append_runtime_guard_upper_bound_fact(
    std::vector<ExactGuardUpperBoundFact> *facts,
    const semantic::Index expr_id,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  if (expr_id == semantic::kInvalidIndex ||
      ref_source_id == semantic::kInvalidIndex ||
      blocker_source_id == semantic::kInvalidIndex ||
      ref_source_id == blocker_source_id) {
    return;
  }
  for (const auto &fact : *facts) {
    if (fact.expr_id == expr_id &&
        fact.ref_source_id == ref_source_id &&
        fact.blocker_source_id == blocker_source_id) {
      return;
    }
  }
  facts->push_back(
      ExactGuardUpperBoundFact{expr_id, ref_source_id, blocker_source_id});
}

inline bool runtime_condition_empty(
    const ExactRuntimeTermCondition &condition) {
  return condition.exact_source_id == semantic::kInvalidIndex &&
         condition.upper_bound_source_ids.empty() &&
         condition.lower_bound_source_ids.empty() &&
         condition.upper_bound_expr_ids.empty() &&
         condition.source_order_facts.empty() &&
         condition.guard_upper_bound_facts.empty();
}

inline ForcedExprEvaluator *prepare_scenario_evaluator(
    const ExactScenarioRuntimeView &scenario_view,
    ForcedExprEvaluator *parent,
    ForcedExprWorkspace *workspace) {
  RelationView view;
  if (!relation_view_with_overlay(
          parent->relation_view(),
          scenario_view.scenario.relation_template,
          &view)) {
    return nullptr;
  }
  workspace->reset(parent->oracle(), view);
  return &workspace->evaluator;
}

inline ForcedExprEvaluator *prepare_runtime_scenario_evaluator(
    const ExactRuntimeScenarioFormula &scenario,
    ForcedExprEvaluator *parent,
    ForcedExprWorkspace *workspace) {
  RelationView view;
  if (!relation_view_with_overlay(
          parent->relation_view(),
          scenario.relation_template,
          &view)) {
    return nullptr;
  }
  workspace->reset(parent->oracle(), view);
  return &workspace->evaluator;
}

inline bool runtime_factors_empty(const ExactRuntimeFactors &factors) {
  return factors.source_pdf.empty() &&
         factors.source_cdf.empty() &&
         factors.source_survival.empty() &&
         factors.expr_density.empty() &&
         factors.expr_cdf.empty() &&
         factors.expr_survival.empty();
}

inline double evaluate_runtime_factors(const ExactRuntimeFactors &factors,
                                       ForcedExprEvaluator *parent,
                                       ForcedExprEvaluator *scenario) {
  double value = 1.0;
  auto multiply = [&](const double factor) {
    value *= factor;
    return std::isfinite(value) && value != 0.0;
  };

  for (const auto source_id : factors.source_pdf) {
    if (!multiply(parent->source_pdf(source_id))) {
      return 0.0;
    }
  }
  for (const auto source_id : factors.source_cdf) {
    if (!multiply(parent->source_cdf(source_id))) {
      return 0.0;
    }
  }
  for (const auto source_id : factors.source_survival) {
    if (!multiply(parent->source_survival(source_id))) {
      return 0.0;
    }
  }

  if (scenario == nullptr &&
      (!factors.expr_density.empty() ||
       !factors.expr_cdf.empty() ||
       !factors.expr_survival.empty())) {
    return 0.0;
  }
  for (const auto expr_id : factors.expr_density) {
    if (!multiply(scenario->expr_density(expr_id))) {
      return 0.0;
    }
  }
  for (const auto expr_id : factors.expr_cdf) {
    if (!multiply(scenario->expr_cdf(expr_id))) {
      return 0.0;
    }
  }
  for (const auto expr_id : factors.expr_survival) {
    if (!multiply(scenario->expr_survival(expr_id))) {
      return 0.0;
    }
  }

  return value;
}

inline semantic::Index runtime_product_term_source_density_id(
    const ExactRuntimeProductTerm &term) {
  if (term.factors.source_pdf.size() != 1U ||
      !term.factors.expr_density.empty()) {
    return semantic::kInvalidIndex;
  }
  return term.factors.source_pdf.front();
}

inline void append_runtime_factor_source_conditions(
    const ExactRuntimeProductTerm &term,
    ExactRuntimeTermCondition *condition) {
  for (const auto source_id : term.factors.source_cdf) {
    append_runtime_condition_index(
        &condition->upper_bound_source_ids,
        source_id);
  }
  for (const auto source_id : term.factors.source_survival) {
    append_runtime_condition_index(
        &condition->lower_bound_source_ids,
        source_id);
  }
}

inline bool simple_event_guard_sources(const ExactVariantPlan &plan,
                                       const semantic::Index expr_id,
                                       semantic::Index *ref_source_id,
                                       semantic::Index *blocker_source_id) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind != semantic::ExprKind::Guard ||
      program.expr_arg_offsets[pos] != program.expr_arg_offsets[pos + 1U]) {
    return false;
  }
  const auto ref = program.expr_ref_child[pos];
  const auto blocker = program.expr_blocker_child[pos];
  const auto ref_kind =
      static_cast<semantic::ExprKind>(program.expr_kind[static_cast<std::size_t>(ref)]);
  const auto blocker_kind = static_cast<semantic::ExprKind>(
      program.expr_kind[static_cast<std::size_t>(blocker)]);
  if (ref_kind != semantic::ExprKind::Event ||
      blocker_kind != semantic::ExprKind::Event) {
    return false;
  }
  *ref_source_id = program.expr_source_ids[static_cast<std::size_t>(ref)];
  *blocker_source_id =
      program.expr_source_ids[static_cast<std::size_t>(blocker)];
  return true;
}

inline void append_runtime_factor_expr_conditions(
    const ExactRuntimeProductTerm &term,
    const ExactVariantPlan &plan,
    ExactRuntimeTermCondition *condition) {
  for (const auto expr_id : term.factors.expr_cdf) {
    semantic::Index ref_source_id{semantic::kInvalidIndex};
    semantic::Index blocker_source_id{semantic::kInvalidIndex};
    if (simple_event_guard_sources(
            plan, expr_id, &ref_source_id, &blocker_source_id)) {
      append_runtime_guard_upper_bound_fact(
          &condition->guard_upper_bound_facts,
          expr_id,
          ref_source_id,
          blocker_source_id);
      append_runtime_source_order_fact(
          &condition->source_order_facts,
          ref_source_id,
          blocker_source_id);
      continue;
    }
    append_runtime_condition_index(&condition->upper_bound_expr_ids, expr_id);
  }
}

inline ExactRuntimeTermCondition runtime_product_term_condition(
    const ExactRuntimeProductTerm &term,
    const ExactVariantPlan &plan) {
  ExactRuntimeTermCondition condition;
  append_runtime_factor_source_conditions(term, &condition);
  append_runtime_factor_expr_conditions(term, plan, &condition);
  const auto source_id = runtime_product_term_source_density_id(term);
  if (source_id != semantic::kInvalidIndex) {
    condition.exact_source_id = source_id;
    return condition;
  }

  if (term.factors.expr_density.size() != 1U ||
      !term.factors.source_pdf.empty()) {
    return condition;
  }
  const auto expr_id = term.factors.expr_density.front();
  semantic::Index ref_source_id{semantic::kInvalidIndex};
  semantic::Index blocker_source_id{semantic::kInvalidIndex};
  if (!simple_event_guard_sources(
          plan, expr_id, &ref_source_id, &blocker_source_id)) {
    append_runtime_condition_index(&condition.upper_bound_expr_ids, expr_id);
    return condition;
  }

  append_runtime_guard_upper_bound_fact(
      &condition.guard_upper_bound_facts,
      expr_id,
      ref_source_id,
      blocker_source_id);
  condition.exact_source_id = ref_source_id;
  append_runtime_condition_index(
      &condition.lower_bound_source_ids,
      blocker_source_id);
  append_runtime_source_order_fact(
      &condition.source_order_facts,
      ref_source_id,
      blocker_source_id);
  return condition;
}

inline ExactRuntimeTermCondition runtime_scenario_tail_condition(
    const ExactRuntimeScenarioFormula &scenario_formula) {
  ExactRuntimeTermCondition condition;
  for (const auto source_id :
       scenario_formula.after_survival.product.source_survival) {
    append_runtime_condition_index(
        &condition.lower_bound_source_ids,
        source_id);
  }
  for (const auto &fact : scenario_formula.source_order_facts) {
    append_runtime_source_order_fact(
        &condition.source_order_facts,
        fact.before_source_id,
        fact.after_source_id);
  }
  return condition;
}

inline double evaluate_runtime_truth_formula(
    const ExactRuntimeTruthFormula &formula,
    const ExactRuntimeScenarioFormula &scenario_formula,
    ForcedExprEvaluator *parent,
    ForcedExprWorkspace *workspace) {
  ForcedExprEvaluator *scenario_evaluator = nullptr;
  if (formula.requires_scenario) {
    scenario_evaluator =
        prepare_runtime_scenario_evaluator(scenario_formula, parent, workspace);
    if (scenario_evaluator == nullptr) {
      return 0.0;
    }
  }

  if (!formula.sum_of_products) {
    const double value =
        runtime_factors_empty(formula.product)
            ? formula.empty_value
            : evaluate_runtime_factors(
                  formula.product, parent, scenario_evaluator);
    return clamp_probability(value);
  }

  double total = formula.sum_terms.empty() ? formula.empty_value : 0.0;
  for (const auto &term : formula.sum_terms) {
    const double value =
        evaluate_runtime_factors(term.factors, parent, scenario_evaluator);
    total += value;
  }
  return formula.clean_signed ? clean_signed_value(total) : total;
}

inline double runtime_readiness_cdf(
    const ExactRuntimeScenarioFormula &scenario_formula,
    ForcedExprEvaluator *evaluator,
    const double t,
    ForcedExprWorkspace *workspace) {
  if (t < 0.0) {
    return 0.0;
  }
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  return evaluate_runtime_truth_formula(
      scenario_formula.readiness_cdf,
      scenario_formula,
      evaluator,
      workspace);
}

inline double runtime_readiness_density(
    const ExactRuntimeScenarioFormula &scenario_formula,
    ForcedExprEvaluator *evaluator,
    const double t,
    ForcedExprWorkspace *workspace) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  return evaluate_runtime_truth_formula(
      scenario_formula.readiness_density,
      scenario_formula,
      evaluator,
      workspace);
}

inline double runtime_after_survival(
    const ExactRuntimeScenarioFormula &scenario_formula,
    ForcedExprEvaluator *evaluator,
    const double t,
    ForcedExprWorkspace *workspace) {
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  return evaluate_runtime_truth_formula(
      scenario_formula.after_survival,
      scenario_formula,
      evaluator,
      workspace);
}

inline double runtime_same_active_win_mass(
    const ExactRuntimeScenarioFormula &scenario_formula,
    ForcedExprEvaluator *evaluator,
    const double t,
    const double ready_upper,
    ForcedExprWorkspace *workspace) {
  if (!(ready_upper > 0.0)) {
    return 0.0;
  }
  const double tail =
      runtime_after_survival(scenario_formula, evaluator, t, workspace);
  if (!(tail > 0.0)) {
    return 0.0;
  }
  return tail *
         runtime_readiness_cdf(
             scenario_formula, evaluator, ready_upper, workspace);
}

inline double runtime_scenario_truth_cdf(
    const ExactRuntimeScenarioFormula &scenario_formula,
    ForcedExprEvaluator *evaluator,
    const double t,
    ForcedExprWorkspace *workspace) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  return clamp_probability(quadrature::integrate_finite_default(
      [&](const double u) {
        const auto guard = evaluator->oracle()->conditional_time_guard(u);
        double value = evaluator->source_pdf(scenario_formula.active_source_id);
        if (!(value > 0.0)) {
          return 0.0;
        }
        value *= runtime_after_survival(
            scenario_formula, evaluator, u, workspace);
        if (!(value > 0.0)) {
          return 0.0;
        }
        if (scenario_formula.has_readiness) {
          value *= runtime_readiness_cdf(
              scenario_formula, evaluator, u, workspace);
          if (!(value > 0.0)) {
            return 0.0;
          }
        }
        return value;
      },
      0.0,
      t));
}

inline double readiness_cdf(const ExactScenarioRuntimeView &scenario_view,
                            ForcedExprEvaluator *evaluator,
                            const double t,
                            ForcedExprWorkspace *workspace) {
  const bool has_before = scenario_view.before_source_count > 0;
  const bool has_ready = scenario_view.ready_expr_count > 0;
  if (!has_before && !has_ready) {
    return 1.0;
  }
  if (t < 0.0) {
    return 0.0;
  }
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  double value = 1.0;
  if (has_before) {
    for (semantic::Index i = 0; i < scenario_view.before_source_count; ++i) {
      const auto source_id =
          scenario_view.before_sources[static_cast<std::size_t>(i)];
      value *= evaluator->source_cdf(source_id);
      if (!(value > 0.0)) {
        return 0.0;
      }
    }
  }
  if (!has_ready) {
    return clamp_probability(value);
  }
  auto *scenario_evaluator =
      prepare_scenario_evaluator(scenario_view, evaluator, workspace);
  if (scenario_evaluator == nullptr) {
    return 0.0;
  }
  for (semantic::Index i = 0; i < scenario_view.ready_expr_count; ++i) {
    const auto expr_idx =
        scenario_view.ready_exprs[static_cast<std::size_t>(i)];
    value *= scenario_evaluator->expr_cdf(expr_idx);
    if (!(value > 0.0)) {
      return 0.0;
    }
  }
  return clamp_probability(value);
}

inline double readiness_density(const ExactScenarioRuntimeView &scenario_view,
                                ForcedExprEvaluator *evaluator,
                                const double t,
                                ForcedExprWorkspace *workspace) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  const bool has_before = scenario_view.before_source_count > 0;
  const bool has_ready = scenario_view.ready_expr_count > 0;
  if (!has_before && !has_ready) {
    return 0.0;
  }
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  ForcedExprEvaluator *scenario_evaluator = nullptr;
  if (has_ready) {
    scenario_evaluator =
        prepare_scenario_evaluator(scenario_view, evaluator, workspace);
    if (scenario_evaluator == nullptr) {
      return 0.0;
    }
  }
  double total = 0.0;
  if (has_before) {
    for (semantic::Index i = 0; i < scenario_view.before_source_count; ++i) {
      double term = evaluator->source_pdf(
          scenario_view.before_sources[static_cast<std::size_t>(i)]);
      if (!std::isfinite(term) || term == 0.0) {
        continue;
      }
      for (semantic::Index j = 0; j < scenario_view.before_source_count; ++j) {
        if (j == i) {
          continue;
        }
        term *= evaluator->source_cdf(
            scenario_view.before_sources[static_cast<std::size_t>(j)]);
        if (!std::isfinite(term) || term == 0.0) {
          break;
        }
      }
      if (!std::isfinite(term) || term == 0.0) {
        continue;
      }
      if (has_ready) {
        for (semantic::Index j = 0; j < scenario_view.ready_expr_count; ++j) {
          const auto expr_idx =
              scenario_view.ready_exprs[static_cast<std::size_t>(j)];
          term *= scenario_evaluator->expr_cdf(expr_idx);
          if (!std::isfinite(term) || term == 0.0) {
            break;
          }
        }
      }
      total += term;
    }
  }
  if (!has_ready) {
    return clean_signed_value(total);
  }
  for (semantic::Index i = 0; i < scenario_view.ready_expr_count; ++i) {
    double term = scenario_evaluator->expr_density(
        scenario_view.ready_exprs[static_cast<std::size_t>(i)]);
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    if (has_before) {
      for (semantic::Index j = 0; j < scenario_view.before_source_count; ++j) {
        const auto source_id =
            scenario_view.before_sources[static_cast<std::size_t>(j)];
        term *= evaluator->source_cdf(source_id);
        if (!std::isfinite(term) || term == 0.0) {
          break;
        }
      }
    }
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (semantic::Index j = 0; j < scenario_view.ready_expr_count; ++j) {
      if (i == j) {
        continue;
      }
      term *= scenario_evaluator->expr_cdf(
          scenario_view.ready_exprs[static_cast<std::size_t>(j)]);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    total += term;
  }
  return clean_signed_value(total);
}

inline double after_survival(const ExactScenarioRuntimeView &scenario_view,
                             ForcedExprEvaluator *evaluator,
                             const double t,
                             ForcedExprWorkspace *workspace) {
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  double value = 1.0;
  if (scenario_view.after_source_count > 0) {
    for (semantic::Index i = 0; i < scenario_view.after_source_count; ++i) {
      const auto source_id =
          scenario_view.after_sources[static_cast<std::size_t>(i)];
      value *= evaluator->source_survival(source_id);
      if (!(value > 0.0)) {
        return 0.0;
      }
    }
  }
  if (scenario_view.tail_expr_count == 0) {
    return clamp_probability(value);
  }
  auto *scenario_evaluator =
      prepare_scenario_evaluator(scenario_view, evaluator, workspace);
  if (scenario_evaluator == nullptr) {
    return 0.0;
  }
  for (semantic::Index i = 0; i < scenario_view.tail_expr_count; ++i) {
    const auto expr_idx =
        scenario_view.tail_exprs[static_cast<std::size_t>(i)];
    value *= scenario_evaluator->expr_survival(expr_idx);
    if (!(value > 0.0)) {
      return 0.0;
    }
  }
  return clamp_probability(value);
}

inline double same_active_win_mass(const ExactScenarioRuntimeView &scenario_view,
                                   ForcedExprEvaluator *evaluator,
                                   const double t,
                                   const double ready_upper,
                                   ForcedExprWorkspace *workspace) {
  if (!(ready_upper > 0.0)) {
    return 0.0;
  }
  const double tail = after_survival(scenario_view, evaluator, t, workspace);
  if (!(tail > 0.0)) {
    return 0.0;
  }
  return tail *
         readiness_cdf(scenario_view, evaluator, ready_upper, workspace);
}

inline double scenario_truth_cdf(const ExactScenarioRuntimeView &scenario_view,
                                 ForcedExprEvaluator *evaluator,
                                 const double t,
                                 ForcedExprWorkspace *workspace) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  return clamp_probability(quadrature::integrate_finite_default(
      [&](const double u) {
        const auto guard = evaluator->oracle()->conditional_time_guard(u);
        double value =
            evaluator->source_pdf(scenario_view.scenario.active_source_id);
        if (!(value > 0.0)) {
          return 0.0;
        }
        value *= after_survival(scenario_view, evaluator, u, workspace);
        if (!(value > 0.0)) {
          return 0.0;
        }
        if (scenario_view.has_readiness()) {
          value *= readiness_cdf(scenario_view, evaluator, u, workspace);
          if (!(value > 0.0)) {
            return 0.0;
          }
        }
        return value;
      },
      0.0,
      t));
}

} // namespace detail
} // namespace accumulatr::eval
