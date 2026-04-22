#pragma once

#include "exact_forced_sources.hpp"

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
    if (cdf_epoch_[pos] == epoch_ && cdf_time_[pos] == current_time) {
      return cdf_cache_[pos];
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
      double surv = 1.0;
      for (semantic::Index i =
               program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        surv *= expr_survival(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      value = clamp_probability(1.0 - surv);
      break;
    }
    case semantic::ExprKind::Guard: {
      if (!(current_time > 0.0)) {
        value = 0.0;
        break;
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
    cdf_epoch_[pos] = epoch_;
    cdf_time_[pos] = current_time;
    cdf_cache_[pos] = value;
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
    if (density_epoch_[pos] == epoch_ && density_time_[pos] == current_time) {
      return density_cache_[pos];
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
          term *= expr_survival(program_.expr_args[static_cast<std::size_t>(j)]);
          if (!std::isfinite(term) || term == 0.0) {
            break;
          }
        }
        total += term;
      }
      value = clean_signed_value(total);
      break;
    }
    case semantic::ExprKind::Guard: {
      const auto ref = program_.expr_ref_child[static_cast<std::size_t>(expr_idx)];
      const auto blocker =
          program_.expr_blocker_child[static_cast<std::size_t>(expr_idx)];
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
    density_epoch_[pos] = epoch_;
    density_time_[pos] = current_time;
    density_cache_[pos] = value;
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
  leaf::EventChannels source_channels(const semantic::Index source_id) const {
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
};

struct ForcedExprWorkspace {
  explicit ForcedExprWorkspace(const ExactVariantPlan &plan) : evaluator(plan) {}

  void reset(ExactSourceOracle *oracle, const RelationView relation_view = {}) {
    evaluator.reset(oracle, relation_view);
  }

  ForcedExprEvaluator evaluator;
};

inline ForcedExprEvaluator *prepare_scenario_evaluator(
    const ExactTransitionScenario &scenario,
    ForcedExprEvaluator *parent,
    ForcedExprWorkspace *workspace) {
  RelationView view;
  if (!relation_view_with_overlay(
          parent->relation_view(), scenario.relation_template, &view)) {
    return nullptr;
  }
  workspace->reset(parent->oracle(), view);
  return &workspace->evaluator;
}

inline double readiness_cdf(const ExactTransitionScenario &scenario,
                            ForcedExprEvaluator *evaluator,
                            const double t,
                            ForcedExprWorkspace *workspace) {
  const bool has_before = !scenario.before_source_ids.empty();
  const bool has_ready = !scenario.ready_exprs.empty();
  if (!has_before && !has_ready) {
    return 1.0;
  }
  if (t < 0.0) {
    return 0.0;
  }
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  double value = 1.0;
  if (has_before) {
    for (const auto source_id : scenario.before_source_ids) {
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
      prepare_scenario_evaluator(scenario, evaluator, workspace);
  if (scenario_evaluator == nullptr) {
    return 0.0;
  }
  for (const auto expr_idx : scenario.ready_exprs) {
    value *= scenario_evaluator->expr_cdf(expr_idx);
    if (!(value > 0.0)) {
      return 0.0;
    }
  }
  return clamp_probability(value);
}

inline double readiness_density(const ExactTransitionScenario &scenario,
                                ForcedExprEvaluator *evaluator,
                                const double t,
                                ForcedExprWorkspace *workspace) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  const bool has_before = !scenario.before_source_ids.empty();
  const bool has_ready = !scenario.ready_exprs.empty();
  if (!has_before && !has_ready) {
    return 0.0;
  }
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  ForcedExprEvaluator *scenario_evaluator = nullptr;
  if (has_ready) {
    scenario_evaluator = prepare_scenario_evaluator(scenario, evaluator, workspace);
    if (scenario_evaluator == nullptr) {
      return 0.0;
    }
  }
  double total = 0.0;
  if (has_before) {
    for (std::size_t i = 0; i < scenario.before_source_ids.size(); ++i) {
      double term = evaluator->source_pdf(scenario.before_source_ids[i]);
      if (!std::isfinite(term) || term == 0.0) {
        continue;
      }
      for (std::size_t j = 0; j < scenario.before_source_ids.size(); ++j) {
        if (j == i) {
          continue;
        }
        term *= evaluator->source_cdf(scenario.before_source_ids[j]);
        if (!std::isfinite(term) || term == 0.0) {
          break;
        }
      }
      if (!std::isfinite(term) || term == 0.0) {
        continue;
      }
      if (has_ready) {
        for (const auto expr_idx : scenario.ready_exprs) {
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
  for (std::size_t i = 0; i < scenario.ready_exprs.size(); ++i) {
    double term = scenario_evaluator->expr_density(scenario.ready_exprs[i]);
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    if (has_before) {
      for (const auto source_id : scenario.before_source_ids) {
        term *= evaluator->source_cdf(source_id);
        if (!std::isfinite(term) || term == 0.0) {
          break;
        }
      }
    }
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (std::size_t j = 0; j < scenario.ready_exprs.size(); ++j) {
      if (i == j) {
        continue;
      }
      term *= scenario_evaluator->expr_cdf(scenario.ready_exprs[j]);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    total += term;
  }
  return clean_signed_value(total);
}

inline double after_survival(const ExactTransitionScenario &scenario,
                             ForcedExprEvaluator *evaluator,
                             const double t,
                             ForcedExprWorkspace *workspace) {
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  double value = 1.0;
  if (!scenario.after_source_ids.empty()) {
    for (const auto source_id : scenario.after_source_ids) {
      value *= evaluator->source_survival(source_id);
      if (!(value > 0.0)) {
        return 0.0;
      }
    }
  }
  if (scenario.tail_exprs.empty()) {
    return clamp_probability(value);
  }
  auto *scenario_evaluator =
      prepare_scenario_evaluator(scenario, evaluator, workspace);
  if (scenario_evaluator == nullptr) {
    return 0.0;
  }
  for (const auto expr_idx : scenario.tail_exprs) {
    value *= scenario_evaluator->expr_survival(expr_idx);
    if (!(value > 0.0)) {
      return 0.0;
    }
  }
  return clamp_probability(value);
}

inline double same_active_win_mass(const ExactTransitionScenario &scenario,
                                   ForcedExprEvaluator *evaluator,
                                   const double t,
                                   const double ready_upper,
                                   ForcedExprWorkspace *workspace) {
  if (!(ready_upper > 0.0)) {
    return 0.0;
  }
  const double tail = after_survival(scenario, evaluator, t, workspace);
  if (!(tail > 0.0)) {
    return 0.0;
  }
  return tail * readiness_cdf(scenario, evaluator, ready_upper, workspace);
}

inline double scenario_truth_cdf(const ExactTransitionScenario &scenario,
                                 ForcedExprEvaluator *evaluator,
                                 const double t,
                                 ForcedExprWorkspace *workspace) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  return clamp_probability(quadrature::integrate_finite_default(
      [&](const double u) {
        const auto guard = evaluator->oracle()->conditional_time_guard(u);
        double value = evaluator->source_pdf(scenario.active_source_id);
        if (!(value > 0.0)) {
          return 0.0;
        }
        value *= after_survival(scenario, evaluator, u, workspace);
        if (!(value > 0.0)) {
          return 0.0;
        }
        if (!scenario.before_source_ids.empty() ||
            !scenario.ready_exprs.empty()) {
          value *= readiness_cdf(scenario, evaluator, u, workspace);
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
