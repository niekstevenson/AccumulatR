#pragma once

#include "exact_forced_sources.hpp"

namespace accumulatr::eval {
namespace detail {

class ForcedExprEvaluator {
public:
  ForcedExprEvaluator(const ExactVariantPlan &plan,
                      ExactSourceOracle *oracle,
                      const std::vector<ExactRelation> *forced_relations)
      : plan_(plan),
        program_(plan.lowered.program),
        oracle_(oracle),
        forced_relations_(forced_relations),
        cdf_cache_(program_.expr_kind.size(), 0.0),
        cdf_time_(program_.expr_kind.size(), 0.0),
        cdf_epoch_(program_.expr_kind.size(), 0U),
        density_cache_(program_.expr_kind.size(), 0.0),
        density_time_(program_.expr_kind.size(), 0.0),
        density_epoch_(program_.expr_kind.size(), 0U) {}

  void reset_forced_relations(const std::vector<ExactRelation> *forced_relations) {
    forced_relations_ = forced_relations;
    ++epoch_;
  }

  double expr_cdf(const semantic::Index expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    if (cdf_epoch_[pos] == epoch_ && cdf_time_[pos] == oracle_time_) {
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
      value = resolve_forced_source_channels(
                 plan_,
                 oracle_,
                 forced_relations_,
                 oracle_time_,
                 ExactSourceKey{child_event_source_kind(program_, expr_idx),
                                child_event_source_index(program_, expr_idx)})
                  .cdf;
      break;
    case semantic::ExprKind::And: {
      double cdf = 1.0;
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        cdf *= expr_cdf(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      value = clamp_probability(cdf);
      break;
    }
    case semantic::ExprKind::Or: {
      double surv = 1.0;
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        surv *= expr_survival(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      value = clamp_probability(1.0 - surv);
      break;
    }
    case semantic::ExprKind::Guard: {
      const double current_time = oracle_time_;
      if (!(current_time > 0.0)) {
        value = 0.0;
        break;
      }
      const double cdf = quadrature::integrate_finite_robust(
          [&](const double u) {
            oracle_time_ = u;
            return expr_density(expr_idx);
          },
          0.0,
          current_time);
      oracle_time_ = current_time;
      value = clamp_probability(cdf);
      break;
    }
    case semantic::ExprKind::Not:
      value = clamp_probability(
          1.0 - expr_cdf(program_.expr_args[static_cast<std::size_t>(
                    program_.expr_arg_offsets[pos])]));
      break;
    }
    cdf_epoch_[pos] = epoch_;
    cdf_time_[pos] = oracle_time_;
    cdf_cache_[pos] = value;
    return value;
  }

  double expr_survival(const semantic::Index expr_idx) {
    return clamp_probability(1.0 - expr_cdf(expr_idx));
  }

  double source_cdf(const ExactSourceKey key) {
    return clamp_probability(
        resolve_forced_source_channels(
            plan_, oracle_, forced_relations_, oracle_time_, key)
            .cdf);
  }

  double source_pdf(const ExactSourceKey key) {
    return safe_density(
        resolve_forced_source_channels(
            plan_, oracle_, forced_relations_, oracle_time_, key)
            .pdf);
  }

  double source_survival(const ExactSourceKey key) {
    return clamp_probability(
        resolve_forced_source_channels(plan_, oracle_, forced_relations_, oracle_time_, key)
            .survival);
  }

  ExactSourceOracle *oracle() const {
    return oracle_;
  }

  const ExactVariantPlan &plan() const {
    return plan_;
  }

  const std::vector<ExactRelation> *forced_relations() const {
    return forced_relations_;
  }

  double expr_density(const semantic::Index expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    if (density_epoch_[pos] == epoch_ && density_time_[pos] == oracle_time_) {
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
      value = resolve_forced_source_channels(
                 plan_,
                 oracle_,
                 forced_relations_,
                 oracle_time_,
                 ExactSourceKey{child_event_source_kind(program_, expr_idx),
                                child_event_source_index(program_, expr_idx)})
                  .pdf;
      break;
    case semantic::ExprKind::And: {
      double total = 0.0;
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        const auto child = program_.expr_args[static_cast<std::size_t>(i)];
        double term = expr_density(child);
        if (!std::isfinite(term) || term == 0.0) {
          continue;
        }
        for (semantic::Index j = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
             j < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
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
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        const auto child = program_.expr_args[static_cast<std::size_t>(i)];
        double term = expr_density(child);
        if (!std::isfinite(term) || term == 0.0) {
          continue;
        }
        for (semantic::Index j = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
             j < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
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
        const double current_time = oracle_time_;
        blocker_survival = clamp_probability(
            1.0 - quadrature::integrate_finite_robust(
                      [&](const double u) {
                        oracle_time_ = u;
                        double term = expr_density(blocker);
                        if (!std::isfinite(term) || term == 0.0) {
                          return 0.0;
                        }
                        for (semantic::Index i =
                                 program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
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
        oracle_time_ = current_time;
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
    density_time_[pos] = oracle_time_;
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
    const double current_time = oracle_time_;
    if (!(current_time > 0.0)) {
      return 0.0;
    }
    const double cdf = quadrature::integrate_finite_robust(
        [&](const double u) {
          oracle_time_ = u;
          return conjunction_density(exprs);
        },
        0.0,
        current_time);
    oracle_time_ = current_time;
    return clamp_probability(cdf);
  }

private:
  const ExactVariantPlan &plan_;
  const runtime::ExactProgram &program_;
  ExactSourceOracle *oracle_;
  const std::vector<ExactRelation> *forced_relations_;
  std::vector<double> cdf_cache_;
  std::vector<double> cdf_time_;
  std::vector<std::uint32_t> cdf_epoch_;
  std::vector<double> density_cache_;
  std::vector<double> density_time_;
  std::vector<std::uint32_t> density_epoch_;
  std::uint32_t epoch_{1U};

public:
  double oracle_time_{0.0};
};

struct ForcedExprWorkspace {
  explicit ForcedExprWorkspace(const ExactVariantPlan &plan,
                               ExactSourceOracle *oracle)
      : relations(static_cast<std::size_t>(plan.source_count),
                  ExactRelation::Unknown),
        evaluator(plan, oracle, nullptr) {}

  std::vector<ExactRelation> relations;
  ForcedExprEvaluator evaluator;
};

inline bool merge_forced_relations(const std::vector<ExactRelation> *base,
                                   const ExactTransitionScenario &scenario,
                                   std::vector<ExactRelation> *out) {
  if (base != nullptr) {
    *out = *base;
  } else {
    std::fill(out->begin(), out->end(), ExactRelation::Unknown);
  }
  for (std::size_t i = 0; i < scenario.forced_source_ids.size(); ++i) {
    const auto source_id = scenario.forced_source_ids[i];
    const auto relation = scenario.forced_source_relations[i];
    auto &slot = (*out)[static_cast<std::size_t>(source_id)];
    if (slot != ExactRelation::Unknown && slot != relation) {
      return false;
    }
    slot = relation;
  }
  return true;
}

inline double readiness_cdf(const ExactTransitionScenario &scenario,
                            ForcedExprEvaluator *evaluator,
                            const double t,
                            ForcedExprWorkspace *workspace) {
  if (scenario.before_keys.empty() && scenario.ready_exprs.empty()) {
    return 1.0;
  }
  if (t < 0.0) {
    return 0.0;
  }
  double value = 1.0;
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  for (const auto &key : scenario.before_keys) {
    value *= evaluator->source_cdf(key);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  if (!merge_forced_relations(
          evaluator->forced_relations(), scenario, &workspace->relations)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
  }
  workspace->evaluator.reset_forced_relations(&workspace->relations);
  workspace->evaluator.oracle_time_ = t;
  for (const auto expr_idx : scenario.ready_exprs) {
    value *= workspace->evaluator.expr_cdf(expr_idx);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  evaluator->oracle_time_ = current_time;
  return clamp_probability(value);
}

inline double readiness_density(const ExactTransitionScenario &scenario,
                                ForcedExprEvaluator *evaluator,
                                const double t,
                                ForcedExprWorkspace *workspace) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  if (scenario.before_keys.empty() && scenario.ready_exprs.empty()) {
    return 0.0;
  }
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  if (!merge_forced_relations(
          evaluator->forced_relations(), scenario, &workspace->relations)) {
    evaluator->oracle_time_ = current_time;
    return 0.0;
  }
  workspace->evaluator.reset_forced_relations(&workspace->relations);
  workspace->evaluator.oracle_time_ = t;
  double total = 0.0;
  for (std::size_t i = 0; i < scenario.before_keys.size(); ++i) {
    double term = evaluator->source_pdf(scenario.before_keys[i]);
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (std::size_t j = 0; j < scenario.before_keys.size(); ++j) {
      if (j == i) {
        continue;
      }
      term *= evaluator->source_cdf(scenario.before_keys[j]);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (const auto expr_idx : scenario.ready_exprs) {
      term *= workspace->evaluator.expr_cdf(expr_idx);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    total += term;
  }
  for (std::size_t i = 0; i < scenario.ready_exprs.size(); ++i) {
    double term = workspace->evaluator.expr_density(scenario.ready_exprs[i]);
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (const auto &key : scenario.before_keys) {
      term *= evaluator->source_cdf(key);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (std::size_t j = 0; j < scenario.ready_exprs.size(); ++j) {
      if (i == j) {
        continue;
      }
      term *= workspace->evaluator.expr_cdf(scenario.ready_exprs[j]);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    total += term;
  }
  evaluator->oracle_time_ = current_time;
  return clean_signed_value(total);
}

inline double after_survival(const ExactTransitionScenario &scenario,
                             ForcedExprEvaluator *evaluator,
                             const double t,
                             ForcedExprWorkspace *workspace) {
  double value = 1.0;
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  for (const auto &key : scenario.after_keys) {
    value *= evaluator->source_survival(key);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  if (!merge_forced_relations(
          evaluator->forced_relations(), scenario, &workspace->relations)) {
    evaluator->oracle_time_ = current_time;
    return 0.0;
  }
  workspace->evaluator.reset_forced_relations(&workspace->relations);
  workspace->evaluator.oracle_time_ = t;
  for (const auto expr_idx : scenario.tail_exprs) {
    value *= workspace->evaluator.expr_survival(expr_idx);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  evaluator->oracle_time_ = current_time;
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
  const double current_time = evaluator->oracle_time_;
  const double integrated = quadrature::integrate_finite_robust(
      [&](const double u) {
        evaluator->oracle_time_ = u;
        double value = evaluator->source_pdf(scenario.active_key);
        if (!(value > 0.0)) {
          return 0.0;
        }
        value *= after_survival(scenario, evaluator, u, workspace);
        if (!(value > 0.0)) {
          return 0.0;
        }
        if (!scenario.before_keys.empty() || !scenario.ready_exprs.empty()) {
          value *= readiness_cdf(scenario, evaluator, u, workspace);
          if (!(value > 0.0)) {
            return 0.0;
          }
        }
        return value;
      },
      0.0,
      t);
  evaluator->oracle_time_ = current_time;
  return clamp_probability(integrated);
}

} // namespace detail
} // namespace accumulatr::eval
