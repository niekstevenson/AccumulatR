#include "exact_outcome_density.h"

#include <Rcpp.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "competitor_cache.h"
#include "integrate.h"
#include "ranked_transitions.h"
#include "runtime_stats.h"

namespace {

struct ExactScenarioPointRef;

enum class RelativeSelfConstraintKind : std::uint8_t {
  kNone = 0u,
  kForcedComplete,
  kForcedSurvive,
  kExactAtEval,
  kBounded,
  kImpossible
};

inline bool exact_density_supported_for_outcome(
    const uuber::NativeContext &ctx, int outcome_idx) {
  if (outcome_idx < 0 ||
      outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return true;
  }
  const uuber::OutcomeContextInfo &outcome =
      ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
  return outcome.alias_sources.empty() && outcome.guess_donors.empty();
}

inline bool exact_forced_contains_id(const uuber::NativeContext &ctx,
                                     const uuber::BitsetState &bits,
                                     bool bits_valid, int source_id) {
  if (!bits_valid || source_id < 0) {
    return false;
  }
  auto bit_it = ctx.ir.label_id_to_bit_idx.find(source_id);
  if (bit_it == ctx.ir.label_id_to_bit_idx.end()) {
    return false;
  }
  return bits.test(bit_it->second);
}

inline ForcedStateView exact_make_point_forced_state(
    const uuber::NativeContext &ctx, const ExactScenarioPointRef &point);

inline bool exact_point_active(const std::vector<std::uint8_t> *active_mask,
                               std::size_t idx);

inline std::size_t exact_count_active_points(
    const std::vector<std::uint8_t> *active_mask, std::size_t point_count);

inline const uuber::TrialParamsSoA *exact_trial_params_soa_for_point(
    const uuber::NativeContext &ctx, const ExactScenarioPointRef &point,
    const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *fallback_trial_params_soa = nullptr);

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask = nullptr,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

template <typename PointVec>
bool exact_eval_simple_acc_event_batch(
    const uuber::NativeContext &ctx, const uuber::IrEvent &event,
    const PointVec &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask = nullptr,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_collect_deterministic_scenarios_batch_from_points_internal(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &aligned_points,
    std::vector<std::uint8_t> &active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

inline void exact_scenario_batch_append_point(ExactScenarioBatch &batch,
                                              ExactScenarioPoint &&point) {
  batch.t.push_back(point.t);
  batch.density_index.push_back(point.density_index);
  batch.weight.push_back(point.weight);
  batch.trial_params_soa.push_back(point.trial_params_soa);
  batch.forced_complete_bits.push_back(std::move(point.forced_complete_bits));
  batch.forced_survive_bits.push_back(std::move(point.forced_survive_bits));
  batch.forced_complete_bits_valid.push_back(point.forced_complete_bits_valid
                                                 ? 1u
                                                 : 0u);
  batch.forced_survive_bits_valid.push_back(point.forced_survive_bits_valid
                                                ? 1u
                                                : 0u);
  batch.time_constraints.push_back(std::move(point.time_constraints));
}

inline void exact_scenario_batch_append_point(ExactScenarioBatch &batch,
                                              const ExactScenarioPoint &point,
                                              double weight_override) {
  batch.t.push_back(point.t);
  batch.density_index.push_back(point.density_index);
  batch.weight.push_back(weight_override);
  batch.trial_params_soa.push_back(point.trial_params_soa);
  batch.forced_complete_bits.push_back(point.forced_complete_bits);
  batch.forced_survive_bits.push_back(point.forced_survive_bits);
  batch.forced_complete_bits_valid.push_back(point.forced_complete_bits_valid
                                                 ? 1u
                                                 : 0u);
  batch.forced_survive_bits_valid.push_back(point.forced_survive_bits_valid
                                                ? 1u
                                                : 0u);
  batch.time_constraints.push_back(point.time_constraints);
}

inline void exact_scenario_batch_append_row(ExactScenarioBatch &dst,
                                            ExactScenarioBatch &src,
                                            std::size_t idx,
                                            double weight_override) {
  dst.t.push_back(src.t[idx]);
  dst.density_index.push_back(src.density_index[idx]);
  dst.weight.push_back(weight_override);
  dst.trial_params_soa.push_back(src.trial_params_soa[idx]);
  dst.forced_complete_bits.push_back(
      std::move(src.forced_complete_bits[idx]));
  dst.forced_survive_bits.push_back(std::move(src.forced_survive_bits[idx]));
  dst.forced_complete_bits_valid.push_back(
      idx < src.forced_complete_bits_valid.size()
          ? src.forced_complete_bits_valid[idx]
          : 0u);
  dst.forced_survive_bits_valid.push_back(
      idx < src.forced_survive_bits_valid.size()
          ? src.forced_survive_bits_valid[idx]
          : 0u);
  dst.time_constraints.push_back(std::move(src.time_constraints[idx]));
}

inline ExactScenarioPoint exact_scenario_batch_point_at(
    const ExactScenarioBatch &batch, std::size_t idx) {
  ExactScenarioPoint point;
  point.t = batch.t[idx];
  point.density_index = batch.density_index[idx];
  point.weight = batch.weight[idx];
  point.trial_params_soa =
      idx < batch.trial_params_soa.size() ? batch.trial_params_soa[idx]
                                          : nullptr;
  point.forced_complete_bits = batch.forced_complete_bits[idx];
  point.forced_survive_bits = batch.forced_survive_bits[idx];
  point.forced_complete_bits_valid =
      idx < batch.forced_complete_bits_valid.size() &&
      batch.forced_complete_bits_valid[idx] != 0u;
  point.forced_survive_bits_valid =
      idx < batch.forced_survive_bits_valid.size() &&
      batch.forced_survive_bits_valid[idx] != 0u;
  point.time_constraints = batch.time_constraints[idx];
  return point;
}

inline void exact_scenario_batch_materialize_points(
    const ExactScenarioBatch &batch, std::vector<ExactScenarioPoint> &points) {
  points.resize(batch.size());
  for (std::size_t i = 0; i < batch.size(); ++i) {
    points[i] = exact_scenario_batch_point_at(batch, i);
  }
}

struct ExactScenarioPointRef {
  double t{0.0};
  std::size_t density_index{0u};
  double weight{0.0};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  const uuber::BitsetState *forced_complete_bits{nullptr};
  const uuber::BitsetState *forced_survive_bits{nullptr};
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  const TimeConstraintMap *time_constraints{nullptr};
};

inline ExactScenarioPointRef exact_point_ref(
    const std::vector<ExactScenarioPoint> &points, std::size_t idx) {
  const ExactScenarioPoint &point = points[idx];
  return ExactScenarioPointRef{point.t,
                               point.density_index,
                               point.weight,
                               point.trial_params_soa,
                               &point.forced_complete_bits,
                               &point.forced_survive_bits,
                               point.forced_complete_bits_valid,
                               point.forced_survive_bits_valid,
                               &point.time_constraints};
}

inline ExactScenarioPointRef exact_point_ref(const ExactScenarioBatch &points,
                                             std::size_t idx) {
  return ExactScenarioPointRef{
      points.t[idx],
      points.density_index[idx],
      points.weight[idx],
      idx < points.trial_params_soa.size() ? points.trial_params_soa[idx]
                                           : nullptr,
      &points.forced_complete_bits[idx],
      &points.forced_survive_bits[idx],
      idx < points.forced_complete_bits_valid.size() &&
          points.forced_complete_bits_valid[idx] != 0u,
      idx < points.forced_survive_bits_valid.size() &&
          points.forced_survive_bits_valid[idx] != 0u,
      &points.time_constraints[idx]};
}

inline void exact_apply_transition_mask_to_point(
    const uuber::NativeContext &ctx, std::vector<ExactScenarioPoint> &points,
    std::size_t idx, const uuber::CompetitorCompiledOp &op) {
  if (idx >= points.size()) {
    return;
  }
  ExactScenarioPoint &point = points[idx];
  apply_transition_mask_words(ctx, op.transition_mask_begin,
                              op.transition_mask_count,
                              op.transition_invalidate_slot,
                              point.forced_survive_bits,
                              point.forced_survive_bits_valid, nullptr, nullptr);
}

inline void exact_apply_transition_mask_to_point(
    const uuber::NativeContext &ctx, ExactScenarioBatch &points,
    std::size_t idx, const uuber::CompetitorCompiledOp &op) {
  if (idx >= points.size()) {
    return;
  }
  bool forced_survive_bits_valid =
      idx < points.forced_survive_bits_valid.size() &&
      points.forced_survive_bits_valid[idx] != 0u;
  apply_transition_mask_words(ctx, op.transition_mask_begin,
                              op.transition_mask_count,
                              op.transition_invalidate_slot,
                              points.forced_survive_bits[idx],
                              forced_survive_bits_valid, nullptr, nullptr);
  if (idx < points.forced_survive_bits_valid.size()) {
    points.forced_survive_bits_valid[idx] =
        forced_survive_bits_valid ? 1u : 0u;
  }
}

inline int exact_find_label_ref_node_idx(const uuber::NativeContext &ctx,
                                         const uuber::LabelRef &ref) {
  if (ref.acc_idx >= 0) {
    for (const uuber::IrEvent &event : ctx.ir.events) {
      if (event.acc_idx == ref.acc_idx && event.node_idx >= 0) {
        return event.node_idx;
      }
    }
  }
  if (ref.pool_idx >= 0) {
    for (const uuber::IrEvent &event : ctx.ir.events) {
      if (event.pool_idx == ref.pool_idx && event.node_idx >= 0) {
        return event.node_idx;
      }
    }
  }
  if (ref.label_id >= 0 && ref.label_id != NA_INTEGER) {
    for (const uuber::IrEvent &event : ctx.ir.events) {
      if (event.label_id == ref.label_id && event.node_idx >= 0) {
        return event.node_idx;
      }
    }
  }
  return -1;
}

template <typename PointVec>
inline void exact_scenario_batch_append_retimed_point(
    ExactScenarioBatch &batch, const PointVec &points, std::size_t idx,
    double t_new) {
  const ExactScenarioPointRef point = exact_point_ref(points, idx);
  batch.t.push_back(t_new);
  batch.density_index.push_back(point.density_index);
  batch.weight.push_back(point.weight);
  batch.trial_params_soa.push_back(point.trial_params_soa);
  batch.forced_complete_bits.push_back(*point.forced_complete_bits);
  batch.forced_survive_bits.push_back(*point.forced_survive_bits);
  batch.forced_complete_bits_valid.push_back(point.forced_complete_bits_valid
                                                 ? 1u
                                                 : 0u);
  batch.forced_survive_bits_valid.push_back(point.forced_survive_bits_valid
                                                ? 1u
                                                : 0u);
  batch.time_constraints.push_back(*point.time_constraints);
}

struct ExactRelativeOnsetBatchPlan {
  bool has_conditioning_bounds{false};
  bool need_density{false};
  bool need_cdf{false};
  bool has_exact{false};
  double exact_time{std::numeric_limits<double>::quiet_NaN()};
  double bound_lower{0.0};
  double bound_upper{std::numeric_limits<double>::infinity()};
  double upper_int{0.0};
  double source_mass{1.0};
  std::size_t lower_cdf_query_idx{std::numeric_limits<std::size_t>::max()};
  std::size_t upper_cdf_query_idx{std::numeric_limits<std::size_t>::max()};
  std::size_t quad_begin{0u};
  std::size_t quad_count{0u};
};

[[noreturn]] inline void exact_relative_onset_batch_invariant(
    const char *reason, const uuber::LabelRef &source_ref,
    std::size_t point_count) {
  Rcpp::stop(
      "Exact relative-onset batch invariant failed: %s "
      "(source_label_id=%d source_acc_idx=%d source_pool_idx=%d points=%d)",
      reason, source_ref.label_id, source_ref.acc_idx, source_ref.pool_idx,
      static_cast<int>(point_count));
}

template <typename PointVec>
bool exact_eval_label_ref_batch(
    const uuber::NativeContext &ctx, const uuber::LabelRef &source_ref,
    const PointVec &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa,
    uuber::KernelNodeBatchValues &out_values) {
  const int source_node_idx = exact_find_label_ref_node_idx(ctx, source_ref);
  if (source_node_idx >= 0) {
    return exact_eval_node_batch_with_points(
        ctx, source_node_idx, points, component_idx, trial_params,
        trial_type_key, need, out_values, active_mask,
        uniform_trial_params_soa);
  }
  if (source_ref.acc_idx >= 0) {
    uuber::IrEvent source_event;
    source_event.node_idx = -1;
    source_event.acc_idx = source_ref.acc_idx;
    source_event.pool_idx = -1;
    source_event.label_id = source_ref.label_id;
    source_event.outcome_idx = source_ref.outcome_idx;
    return exact_eval_simple_acc_event_batch(
        ctx, source_event, points, component_idx, trial_params, trial_type_key,
        need, out_values, active_mask, uniform_trial_params_soa);
  }
  return false;
}

template <typename PointVec>
bool exact_eval_relative_onset_base_batch(
    const uuber::NativeContext &ctx, const PointVec &points,
    const uuber::LabelRef &source_ref, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const uuber::AccDistParams &cfg,
    const LowerBoundTransform &lower_bound,
    const std::vector<double> &success_probs,
    const std::vector<double> &eval_times,
    const std::vector<std::uint8_t> &need_density_mask,
    const std::vector<std::uint8_t> &need_cdf_mask, double x_shift,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa,
    std::vector<double> &base_density_out, std::vector<double> &base_cdf_out) {
  const std::size_t point_count = points.size();
  base_density_out.assign(point_count, 0.0);
  base_cdf_out.assign(point_count, 0.0);
  if (success_probs.size() != point_count || eval_times.size() != point_count ||
      need_density_mask.size() != point_count ||
      need_cdf_mask.size() != point_count) {
    return false;
  }
  if (source_ref.acc_idx < 0 && source_ref.pool_idx < 0) {
    return true;
  }

  const int source_label_id = source_ref.label_id;
  std::vector<ExactRelativeOnsetBatchPlan> plans(point_count);
  ExactScenarioBatch cdf_query_points;
  ExactScenarioBatch density_query_points;
  std::vector<double> density_query_weights;
  std::vector<double> density_query_shifted_times;
  bool any_density_query = false;
  bool any_cdf_query = false;
  constexpr int kFusedSegments = 4;

  for (std::size_t i = 0; i < point_count; ++i) {
    if (!exact_point_active(active_mask, i) ||
        (need_density_mask[i] == 0u && need_cdf_mask[i] == 0u)) {
      continue;
    }
    ExactRelativeOnsetBatchPlan &plan = plans[i];
    plan.need_density = need_density_mask[i] != 0u;
    plan.need_cdf = need_cdf_mask[i] != 0u;

    const ExactScenarioPointRef point = exact_point_ref(points, i);
    const double eval_t = eval_times[i];
    const ForcedStateView forced_state = exact_make_point_forced_state(ctx, point);
    if (source_label_id >= 0 && source_label_id != NA_INTEGER &&
        evaluator_state_contains_survive_at(
            forced_state, point.time_constraints, source_label_id, eval_t)) {
      continue;
    }

    if (source_label_id >= 0 && source_label_id != NA_INTEGER) {
      if (forced_state_contains_complete(forced_state, source_label_id)) {
        plan.bound_upper = std::min(plan.bound_upper, eval_t);
        plan.has_conditioning_bounds = true;
      }
      bool source_has_bounds = false;
      double source_bound_lower = 0.0;
      double source_bound_upper = std::numeric_limits<double>::infinity();
      (void)evaluator_resolve_label_time_constraint(
          point.time_constraints, source_label_id, plan.has_exact,
          plan.exact_time, source_has_bounds, source_bound_lower,
          source_bound_upper);
      if (source_has_bounds) {
        plan.bound_lower = std::max(plan.bound_lower, source_bound_lower);
        plan.bound_upper = std::min(plan.bound_upper, source_bound_upper);
        plan.has_conditioning_bounds = true;
      }
    }

    if (plan.has_exact) {
      const double shifted_time = eval_t - plan.exact_time - x_shift;
      if (plan.need_density) {
        base_density_out[i] = safe_density(
            success_probs[i] *
            eval_pdf_single_with_lower_bound(cfg, shifted_time, lower_bound));
      }
      if (plan.need_cdf) {
        base_cdf_out[i] = clamp_probability(
            success_probs[i] *
            eval_cdf_single_with_lower_bound(cfg, shifted_time, lower_bound));
      }
      continue;
    }

    if (!(plan.bound_upper > plan.bound_lower)) {
      continue;
    }
    if (std::isfinite(eval_t) && (eval_t - x_shift <= plan.bound_lower)) {
      continue;
    }

    if (plan.has_conditioning_bounds) {
      if (std::isfinite(plan.bound_lower) && plan.bound_lower > 0.0) {
        plan.lower_cdf_query_idx = cdf_query_points.size();
        exact_scenario_batch_append_retimed_point(
            cdf_query_points, points, i, plan.bound_lower);
        any_cdf_query = true;
      }
      if (std::isfinite(plan.bound_upper)) {
        if (plan.bound_upper > 0.0) {
          plan.upper_cdf_query_idx = cdf_query_points.size();
          exact_scenario_batch_append_retimed_point(
              cdf_query_points, points, i, plan.bound_upper);
          any_cdf_query = true;
        }
      } else {
        plan.upper_cdf_query_idx = cdf_query_points.size();
        exact_scenario_batch_append_retimed_point(
            cdf_query_points, points, i, std::numeric_limits<double>::infinity());
        any_cdf_query = true;
      }
    } else if (!std::isfinite(eval_t) && plan.need_cdf) {
      plan.upper_cdf_query_idx = cdf_query_points.size();
      exact_scenario_batch_append_retimed_point(
          cdf_query_points, points, i, std::numeric_limits<double>::infinity());
      any_cdf_query = true;
    }

    if (!std::isfinite(eval_t)) {
      continue;
    }

    plan.upper_int =
        std::min(plan.bound_upper, std::max(plan.bound_lower, eval_t - x_shift));
    if (!(plan.upper_int > plan.bound_lower)) {
      continue;
    }

    plan.quad_begin = density_query_points.size();
    for (int seg = 0; seg < kFusedSegments; ++seg) {
      const double seg_lo =
          plan.bound_lower +
          ((plan.upper_int - plan.bound_lower) * static_cast<double>(seg) /
           static_cast<double>(kFusedSegments));
      const double seg_hi =
          plan.bound_lower +
          ((plan.upper_int - plan.bound_lower) * static_cast<double>(seg + 1) /
           static_cast<double>(kFusedSegments));
      const uuber::TimeBatch batch = uuber::build_time_batch(seg_lo, seg_hi);
      if (batch.nodes.empty() || batch.nodes.size() != batch.weights.size()) {
        exact_relative_onset_batch_invariant(
            "quadrature node construction failed", source_ref, point_count);
      }
      for (std::size_t node_idx = 0; node_idx < batch.nodes.size(); ++node_idx) {
        const double u = batch.nodes[node_idx];
        exact_scenario_batch_append_retimed_point(density_query_points, points, i,
                                                  u);
        density_query_weights.push_back(batch.weights[node_idx]);
        density_query_shifted_times.push_back(eval_t - u - x_shift);
      }
    }
    plan.quad_count = density_query_points.size() - plan.quad_begin;
    any_density_query = true;
  }

  std::vector<double> cdf_query_values;
  if (any_cdf_query) {
    uuber::KernelNodeBatchValues cdf_values;
    if (!exact_eval_label_ref_batch(
            ctx, source_ref, cdf_query_points, component_idx, trial_params,
            trial_type_key, EvalNeed::kCDF, nullptr, uniform_trial_params_soa,
            cdf_values) ||
        cdf_values.cdf.size() != cdf_query_points.size()) {
      exact_relative_onset_batch_invariant(
          "source CDF batch evaluation failed", source_ref, point_count);
    }
    cdf_query_values = std::move(cdf_values.cdf);
  }

  std::vector<double> density_query_values;
  if (any_density_query) {
    uuber::KernelNodeBatchValues source_density_values;
    if (!exact_eval_label_ref_batch(
            ctx, source_ref, density_query_points, component_idx, trial_params,
            trial_type_key, EvalNeed::kDensity, nullptr,
            uniform_trial_params_soa, source_density_values) ||
        source_density_values.density.size() != density_query_points.size()) {
      exact_relative_onset_batch_invariant(
          "source density batch evaluation failed", source_ref, point_count);
    }
    density_query_values = std::move(source_density_values.density);
  }

  std::vector<double> pdf_values;
  std::vector<double> cdf_values;
  if (any_density_query) {
    const bool need_any_density = std::any_of(
        plans.begin(), plans.end(),
        [](const ExactRelativeOnsetBatchPlan &plan) { return plan.need_density; });
    const bool need_any_cdf = std::any_of(
        plans.begin(), plans.end(),
        [](const ExactRelativeOnsetBatchPlan &plan) { return plan.need_cdf; });
    if (need_any_density) {
      pdf_values.assign(density_query_shifted_times.size(), 0.0);
      eval_pdf_vec_with_lower_bound(cfg, lower_bound,
                                    density_query_shifted_times.data(),
                                    density_query_shifted_times.size(),
                                    pdf_values.data());
    }
    if (need_any_cdf) {
      cdf_values.assign(density_query_shifted_times.size(), 0.0);
      eval_cdf_vec_with_lower_bound(cfg, lower_bound,
                                    density_query_shifted_times.data(),
                                    density_query_shifted_times.size(),
                                    cdf_values.data());
    }
  }

  std::vector<double> density_sums(point_count, 0.0);
  std::vector<double> cdf_sums(point_count, 0.0);
  for (std::size_t i = 0; i < point_count; ++i) {
    ExactRelativeOnsetBatchPlan &plan = plans[i];
    if (!exact_point_active(active_mask, i) ||
        (!plan.need_density && !plan.need_cdf) || plan.has_exact) {
      continue;
    }
    if (plan.has_conditioning_bounds) {
      const double lower_cdf =
          plan.lower_cdf_query_idx == std::numeric_limits<std::size_t>::max()
              ? 0.0
              : clamp_probability(cdf_query_values[plan.lower_cdf_query_idx]);
      const double upper_cdf =
          plan.upper_cdf_query_idx == std::numeric_limits<std::size_t>::max()
              ? 0.0
              : clamp_probability(cdf_query_values[plan.upper_cdf_query_idx]);
      plan.source_mass = upper_cdf - lower_cdf;
      if (!std::isfinite(plan.source_mass) || plan.source_mass <= 0.0) {
        continue;
      }
    }
    if (!std::isfinite(eval_times[i])) {
      if (plan.need_cdf) {
        const double cdf_total =
            plan.has_conditioning_bounds
                ? clamp_probability(success_probs[i])
                : clamp_probability(success_probs[i] *
                                    clamp_probability(cdf_query_values[plan.upper_cdf_query_idx]));
        base_cdf_out[i] = cdf_total;
      }
      continue;
    }
    if (plan.quad_count == 0u) {
      continue;
    }
    for (std::size_t query_idx = plan.quad_begin;
         query_idx < plan.quad_begin + plan.quad_count; ++query_idx) {
      const double fs = safe_density(density_query_values[query_idx]);
      if (!std::isfinite(fs) || fs <= 0.0) {
        continue;
      }
      const double normalized_fs =
          plan.has_conditioning_bounds ? (fs / plan.source_mass) : fs;
      if (!std::isfinite(normalized_fs) || normalized_fs <= 0.0) {
        continue;
      }
      if (plan.need_density && !pdf_values.empty()) {
        const double fx = safe_density(pdf_values[query_idx]);
        if (fx > 0.0) {
          density_sums[i] += density_query_weights[query_idx] *
                             normalized_fs * fx;
        }
      }
      if (plan.need_cdf && !cdf_values.empty()) {
        const double Fx = clamp_probability(cdf_values[query_idx]);
        if (Fx > 0.0) {
          cdf_sums[i] += density_query_weights[query_idx] *
                         normalized_fs * Fx;
        }
      }
    }
    if (plan.need_density) {
      base_density_out[i] =
          safe_density(success_probs[i] * density_sums[i]);
    }
    if (plan.need_cdf) {
      base_cdf_out[i] =
          clamp_probability(success_probs[i] * clamp_probability(cdf_sums[i]));
    }
  }

  return true;
}

bool classify_relative_self_constraint(const SourceTimeConstraint *constraint,
                                       double t, bool forced_complete,
                                       bool forced_survive,
                                       RelativeSelfConstraintKind &kind,
                                       bool &has_bounds, double &bound_lower,
                                       double &bound_upper) {
  kind = RelativeSelfConstraintKind::kNone;
  has_bounds = false;
  bound_lower = 0.0;
  bound_upper = std::numeric_limits<double>::infinity();
  const bool self_is_time_conditioned =
      constraint != nullptr && !time_constraint_empty(*constraint);
  if (!self_is_time_conditioned && forced_complete) {
    kind = RelativeSelfConstraintKind::kForcedComplete;
    return true;
  }
  if (forced_survive) {
    kind = self_is_time_conditioned ? RelativeSelfConstraintKind::kImpossible
                                    : RelativeSelfConstraintKind::kForcedSurvive;
    return true;
  }
  if (!self_is_time_conditioned) {
    return true;
  }
  if (constraint->has_exact) {
    if (!time_constraint_same_time(constraint->exact_time, t)) {
      return false;
    }
    if (constraint->has_lower) {
      bound_lower = constraint->lower;
      has_bounds = true;
    }
    if (constraint->has_upper) {
      bound_upper = constraint->upper;
      has_bounds = true;
    }
    kind = RelativeSelfConstraintKind::kExactAtEval;
    return true;
  }
  if (constraint->has_lower || constraint->has_upper) {
    if (constraint->has_lower) {
      bound_lower = constraint->lower;
      has_bounds = true;
    }
    if (constraint->has_upper) {
      bound_upper = constraint->upper;
      has_bounds = true;
    }
    kind = RelativeSelfConstraintKind::kBounded;
    return true;
  }
  return true;
}

template <typename PointVec>
bool exact_eval_simple_acc_event_batch(
    const uuber::NativeContext &ctx, const uuber::IrEvent &event,
    const PointVec &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  out_values = uuber::KernelNodeBatchValues{};
  const std::size_t point_count = points.size();
  out_values.density.assign(point_count, 0.0);
  out_values.survival.assign(point_count, 1.0);
  out_values.cdf.assign(point_count, 0.0);
  if (point_count == 0u || event.acc_idx < 0 ||
      event.acc_idx >= static_cast<int>(ctx.accumulators.size())) {
    return false;
  }
  if (needs_density(need) &&
      !exact_density_supported_for_outcome(ctx, event.outcome_idx)) {
    return false;
  }

  const uuber::NativeAccumulator &acc =
      ctx.accumulators[static_cast<std::size_t>(event.acc_idx)];
  const TrialAccumulatorParams *override =
      evaluator_get_trial_param_entry(trial_params, event.acc_idx);
  if (!evaluator_component_active_idx(acc, component_idx, override)) {
    return true;
  }
  std::size_t anchor_idx = point_count;
  for (std::size_t i = 0; i < point_count; ++i) {
    if (exact_point_active(active_mask, i)) {
      anchor_idx = i;
      break;
    }
  }
  if (anchor_idx >= point_count) {
    return true;
  }
  const uuber::TrialParamsSoA *trial_params_soa =
      uniform_trial_params_soa != nullptr
          ? uniform_trial_params_soa
          : exact_trial_params_soa_for_point(ctx, exact_point_ref(points, anchor_idx),
                                             trial_params);
  if (!trial_params_soa) {
    return false;
  }

  const int onset_kind = override ? override->onset_kind : acc.onset_kind;
  const int onset_source_acc_idx =
      override ? override->onset_source_acc_idx : acc.onset_source_acc_idx;
  const int onset_source_pool_idx =
      override ? override->onset_source_pool_idx : acc.onset_source_pool_idx;
  const double onset_lag = override ? override->onset_lag : acc.onset_lag;

  double onset = 0.0;
  double q = 0.0;
  uuber::AccDistParams cfg;
  evaluator_resolve_event_numeric_params(
      acc, event.acc_idx, override, trial_params_soa, onset, q, cfg);
  const LowerBoundTransform lower_bound = default_lower_bound_transform(cfg);
  const double onset_eff = total_onset_with_t0(onset, cfg);
  const double success_prob = clamp_probability(1.0 - q);
  const int label_id = event.label_id;
  if (label_id == NA_INTEGER || label_id < 0) {
    return false;
  }
  const bool needs_relative_onset =
      ctx.has_chained_onsets && onset_kind != uuber::ONSET_ABSOLUTE;
  const uuber::LabelRef onset_source_ref = evaluator_make_onset_source_ref(
      ctx, onset_kind, onset_source_acc_idx, onset_source_pool_idx);
  const int onset_source_label_id = onset_source_ref.label_id;

  std::vector<RelativeSelfConstraintKind> kinds(
      point_count, RelativeSelfConstraintKind::kNone);
  std::vector<std::uint8_t> self_has_bounds(point_count, 0u);
  std::vector<double> self_bound_lower(point_count, 0.0);
  std::vector<double> self_bound_upper(
      point_count, std::numeric_limits<double>::infinity());
  std::vector<double> shifted_times(point_count, 0.0);
  std::vector<double> eval_times(point_count, 0.0);
  std::vector<double> success_probs(point_count, success_prob);
  std::vector<std::uint8_t> direct_shifted_eval(point_count, 0u);
  std::vector<std::uint8_t> need_base_density_mask(point_count, 0u);
  std::vector<std::uint8_t> need_base_cdf_mask(point_count, 0u);
  std::vector<std::uint8_t> bounded_lower_need_cdf(point_count, 0u);
  std::vector<std::uint8_t> bounded_upper_need_cdf(point_count, 0u);
  bool need_pdf_vec = false;
  bool need_cdf_vec = false;

  for (std::size_t i = 0; i < point_count; ++i) {
    if (active_mask &&
        (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    const ExactScenarioPointRef point = exact_point_ref(points, i);
    if (uniform_trial_params_soa == nullptr) {
      const uuber::TrialParamsSoA *point_trial_params_soa =
          exact_trial_params_soa_for_point(ctx, point, trial_params,
                                           trial_params_soa);
      if (!point_trial_params_soa ||
          !same_acc_batch_params_except_q(*trial_params_soa,
                                          *point_trial_params_soa,
                                          event.acc_idx)) {
        return false;
      }
      if (point_trial_params_soa->valid && event.acc_idx >= 0 &&
          event.acc_idx < point_trial_params_soa->n_acc &&
          static_cast<std::size_t>(event.acc_idx) <
              point_trial_params_soa->q.size()) {
        success_probs[i] = clamp_probability(
            1.0 -
            point_trial_params_soa->q[static_cast<std::size_t>(event.acc_idx)]);
      }
    }
    const bool forced_complete = exact_forced_contains_id(
        ctx, *point.forced_complete_bits, point.forced_complete_bits_valid,
        label_id);
    const bool forced_survive = exact_forced_contains_id(
        ctx, *point.forced_survive_bits, point.forced_survive_bits_valid,
        label_id);
    const SourceTimeConstraint *constraint =
        time_constraints_find(point.time_constraints, label_id);
    bool has_bounds = false;
    double bound_lower = 0.0;
    double bound_upper = std::numeric_limits<double>::infinity();
    if (!classify_relative_self_constraint(constraint, point.t, forced_complete,
                                          forced_survive, kinds[i], has_bounds,
                                          bound_lower, bound_upper)) {
      return false;
    }
    self_has_bounds[i] = has_bounds ? 1u : 0u;
    self_bound_lower[i] = bound_lower;
    self_bound_upper[i] = bound_upper;
    eval_times[i] = point.t;
    const bool need_base_density =
        needs_density(need) &&
        (kinds[i] == RelativeSelfConstraintKind::kNone ||
         kinds[i] == RelativeSelfConstraintKind::kBounded);
    const bool need_base_cdf =
        (needs_cdf(need) || needs_survival(need) ||
         kinds[i] == RelativeSelfConstraintKind::kBounded) &&
        (kinds[i] == RelativeSelfConstraintKind::kNone ||
         kinds[i] == RelativeSelfConstraintKind::kBounded);
    need_base_density_mask[i] = need_base_density ? 1u : 0u;
    need_base_cdf_mask[i] = need_base_cdf ? 1u : 0u;
    if (!needs_relative_onset) {
      direct_shifted_eval[i] = 1u;
      shifted_times[i] = point.t - onset_eff;
    } else {
      const SourceTimeConstraint *source_constraint =
          time_constraints_find(point.time_constraints, onset_source_label_id);
      if (source_constraint != nullptr && source_constraint->has_exact &&
          std::isfinite(source_constraint->exact_time) &&
          !exact_forced_contains_id(ctx, *point.forced_survive_bits,
                                    point.forced_survive_bits_valid,
                                    onset_source_label_id)) {
        direct_shifted_eval[i] = 1u;
        shifted_times[i] =
            point.t - source_constraint->exact_time - onset_lag - onset_eff;
      }
    }
    if (kinds[i] == RelativeSelfConstraintKind::kBounded &&
        direct_shifted_eval[i] == 0u) {
      bounded_lower_need_cdf[i] = std::isfinite(bound_lower) ? 1u : 0u;
      bounded_upper_need_cdf[i] = std::isfinite(bound_upper) ? 1u : 0u;
    }
    if (direct_shifted_eval[i] != 0u && need_base_density) {
      need_pdf_vec = true;
    }
    if (direct_shifted_eval[i] != 0u && need_base_cdf) {
      need_cdf_vec = true;
    }
  }

  std::vector<double> pdf_values;
  std::vector<double> cdf_values;
  if (need_pdf_vec) {
    pdf_values.assign(point_count, 0.0);
    eval_pdf_vec_with_lower_bound(cfg, lower_bound, shifted_times.data(),
                                  point_count, pdf_values.data());
  }
  if (need_cdf_vec) {
    cdf_values.assign(point_count, 0.0);
    eval_cdf_vec_with_lower_bound(cfg, lower_bound, shifted_times.data(),
                                  point_count, cdf_values.data());
  }

  auto any_mask_set = [](const std::vector<std::uint8_t> &mask) -> bool {
    return std::any_of(mask.begin(), mask.end(),
                       [](std::uint8_t value) { return value != 0u; });
  };
  std::vector<double> relative_base_density;
  std::vector<double> relative_base_cdf;
  std::vector<double> relative_lower_cdf;
  std::vector<double> relative_upper_cdf;
  std::vector<double> relative_edge_density_unused;
  if (needs_relative_onset) {
    if (!exact_eval_relative_onset_base_batch(
            ctx, points, onset_source_ref, component_idx, trial_params,
            trial_type_key, cfg, lower_bound, success_probs,
            eval_times, need_base_density_mask, need_base_cdf_mask,
            onset_lag + onset_eff, active_mask, uniform_trial_params_soa,
            relative_base_density, relative_base_cdf)) {
      return false;
    }
    const std::vector<std::uint8_t> zero_need(point_count, 0u);
    if (any_mask_set(bounded_lower_need_cdf) &&
        !exact_eval_relative_onset_base_batch(
            ctx, points, onset_source_ref, component_idx, trial_params,
            trial_type_key, cfg, lower_bound, success_probs,
            self_bound_lower, zero_need, bounded_lower_need_cdf,
            onset_lag + onset_eff, active_mask, uniform_trial_params_soa,
            relative_edge_density_unused, relative_lower_cdf)) {
      return false;
    }
    if (any_mask_set(bounded_upper_need_cdf) &&
        !exact_eval_relative_onset_base_batch(
            ctx, points, onset_source_ref, component_idx, trial_params,
            trial_type_key, cfg, lower_bound, success_probs,
            self_bound_upper, zero_need, bounded_upper_need_cdf,
            onset_lag + onset_eff, active_mask, uniform_trial_params_soa,
            relative_edge_density_unused, relative_upper_cdf)) {
      return false;
    }
  }

  for (std::size_t i = 0; i < point_count; ++i) {
    if (active_mask &&
        (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    const bool need_base_density = need_base_density_mask[i] != 0u;
    const bool need_base_cdf = need_base_cdf_mask[i] != 0u;
    double base_density = 0.0;
    double base_cdf = 0.0;
    if (direct_shifted_eval[i] != 0u) {
      if (need_base_density) {
        base_density = safe_density(success_probs[i] * pdf_values[i]);
      }
      if (need_base_cdf) {
        base_cdf = clamp_probability(success_probs[i] * cdf_values[i]);
      }
    } else if (needs_relative_onset) {
      if (need_base_density && i < relative_base_density.size()) {
        base_density = relative_base_density[i];
      }
      if (need_base_cdf && i < relative_base_cdf.size()) {
        base_cdf = relative_base_cdf[i];
      }
    }
    switch (kinds[i]) {
    case RelativeSelfConstraintKind::kNone: {
      if (needs_density(need)) {
        out_values.density[i] = safe_density(base_density);
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] = clamp_probability(base_cdf);
      }
      if (needs_survival(need)) {
        out_values.survival[i] = clamp_probability(1.0 - base_cdf);
      }
      break;
    }
    case RelativeSelfConstraintKind::kForcedComplete:
    case RelativeSelfConstraintKind::kExactAtEval:
      if (self_has_bounds[i] != 0u &&
          (!(exact_point_ref(points, i).t > self_bound_lower[i]) ||
           exact_point_ref(points, i).t > self_bound_upper[i])) {
        break;
      }
      if (needs_survival(need)) {
        out_values.survival[i] = 0.0;
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] = 1.0;
      }
      break;
    case RelativeSelfConstraintKind::kForcedSurvive:
      if (needs_survival(need)) {
        out_values.survival[i] = 1.0;
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] = 0.0;
      }
      break;
    case RelativeSelfConstraintKind::kImpossible:
      if (needs_density(need)) {
        out_values.density[i] = 0.0;
      }
      if (needs_survival(need)) {
        out_values.survival[i] = 0.0;
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] = 0.0;
      }
      break;
    case RelativeSelfConstraintKind::kBounded: {
      const double lower = self_bound_lower[i];
      const double upper = self_bound_upper[i];
      const double lower_cdf =
          std::isfinite(lower)
              ? (direct_shifted_eval[i] != 0u
                     ? clamp_probability(success_probs[i] * eval_cdf_single_with_lower_bound(
                                               cfg, lower - (exact_point_ref(points, i).t - shifted_times[i]),
                                               lower_bound))
                     : (i < relative_lower_cdf.size()
                            ? clamp_probability(relative_lower_cdf[i])
                            : 0.0))
              : 0.0;
      const double upper_cdf =
          std::isfinite(upper)
              ? (direct_shifted_eval[i] != 0u
                     ? clamp_probability(success_probs[i] * eval_cdf_single_with_lower_bound(
                                               cfg, upper - (exact_point_ref(points, i).t - shifted_times[i]),
                                               lower_bound))
                     : (i < relative_upper_cdf.size()
                            ? clamp_probability(relative_upper_cdf[i])
                            : clamp_probability(success_probs[i])))
              : clamp_probability(success_probs[i]);
      const double condition_mass = clamp_probability(upper_cdf - lower_cdf);
      if (!(upper > lower) || !(condition_mass > 0.0)) {
        break;
      }
      if (needs_density(need)) {
        if (std::isfinite(exact_point_ref(points, i).t) &&
            exact_point_ref(points, i).t > lower &&
            exact_point_ref(points, i).t <= upper) {
          out_values.density[i] =
              safe_density(base_density / condition_mass);
        } else {
          out_values.density[i] = 0.0;
        }
      }
      if (needs_cdf(need) || needs_survival(need)) {
        double cdf_cond = 0.0;
        if (!std::isfinite(exact_point_ref(points, i).t) ||
            exact_point_ref(points, i).t >= upper) {
          cdf_cond = 1.0;
        } else if (!(exact_point_ref(points, i).t > lower)) {
          cdf_cond = 0.0;
        } else {
          cdf_cond =
              clamp_probability((base_cdf - lower_cdf) / condition_mass);
        }
        if (needs_cdf(need)) {
          out_values.cdf[i] = cdf_cond;
        }
        if (needs_survival(need)) {
          out_values.survival[i] = 1.0 - cdf_cond;
        }
      }
      break;
    }
    }
  }

  return true;
}

inline bool time_constraint_equal(const SourceTimeConstraint &lhs,
                                  const SourceTimeConstraint &rhs) {
  return lhs.has_exact == rhs.has_exact &&
         lhs.has_lower == rhs.has_lower &&
         lhs.has_upper == rhs.has_upper &&
         (!lhs.has_exact ||
          time_constraint_same_time(lhs.exact_time, rhs.exact_time)) &&
         (!lhs.has_lower ||
          time_constraint_same_time(lhs.lower, rhs.lower)) &&
         (!lhs.has_upper ||
          time_constraint_same_time(lhs.upper, rhs.upper));
}

inline bool filtered_forced_contains_id(const uuber::BitsetState &bits,
                                        bool bits_valid,
                                        int source_id,
                                        const uuber::NativeContext &ctx) {
  if (!bits_valid || source_id < 0) {
    return false;
  }
  auto bit_it = ctx.ir.label_id_to_bit_idx.find(source_id);
  if (bit_it == ctx.ir.label_id_to_bit_idx.end()) {
    return false;
  }
  return bits.test(bit_it->second);
}

inline void build_filtered_time_constraints(
    const TimeConstraintMap &source,
    const std::vector<int> &relevant_source_ids,
    TimeConstraintMap &filtered) {
  filtered.clear();
  for (int source_id : relevant_source_ids) {
    auto it = source.find(source_id);
    if (it != source.end()) {
      filtered.emplace(it->first, it->second);
    }
  }
}

inline void build_filtered_forced_bits(
    const uuber::NativeContext &ctx,
    const uuber::BitsetState &source_bits, bool source_bits_valid,
    const std::vector<int> &relevant_source_ids,
    uuber::BitsetState &filtered_bits, bool &filtered_bits_valid) {
  filtered_bits_valid = false;
  if (!source_bits_valid) {
    return;
  }
  (void)ensure_forced_bitset_capacity(ctx, filtered_bits, filtered_bits_valid);
  filtered_bits.clear_all();
  for (int source_id : relevant_source_ids) {
    if (filtered_forced_contains_id(source_bits, source_bits_valid, source_id,
                                    ctx)) {
      set_forced_id_bit_strict(ctx, source_id, filtered_bits,
                               filtered_bits_valid);
    }
  }
}

inline std::vector<int> exact_relevant_source_ids_for_node(
    const uuber::NativeContext &ctx, int node_idx) {
  static const std::vector<int> kEmpty;
  if (node_idx < 0 || static_cast<std::size_t>(node_idx) >= ctx.ir.nodes.size()) {
    return kEmpty;
  }
  const uuber::IrNode &node = ir_node_required(ctx, node_idx);

  std::vector<int> ids = ensure_source_ids(ctx, node);

  if ((node.op == uuber::IrNodeOp::EventAcc ||
       node.op == uuber::IrNodeOp::EventPool) &&
      node.event_idx >= 0 &&
      node.event_idx < static_cast<int>(ctx.ir.events.size())) {
    const uuber::IrEvent &event =
        ctx.ir.events[static_cast<std::size_t>(node.event_idx)];
    if (event.acc_idx >= 0 &&
        event.acc_idx < static_cast<int>(ctx.accumulators.size())) {
      const uuber::NativeAccumulator &acc =
          ctx.accumulators[static_cast<std::size_t>(event.acc_idx)];
      if (acc.onset_kind == uuber::ONSET_AFTER_ACCUMULATOR &&
          acc.onset_source_acc_idx >= 0 &&
          acc.onset_source_acc_idx < static_cast<int>(ctx.accumulators.size())) {
        ids.push_back(accumulator_label_id_of(ctx, acc.onset_source_acc_idx));
      } else if (acc.onset_kind == uuber::ONSET_AFTER_POOL &&
                 acc.onset_source_pool_idx >= 0 &&
                 acc.onset_source_pool_idx < static_cast<int>(ctx.pools.size())) {
        ids.push_back(pool_label_id_of(ctx, acc.onset_source_pool_idx));
      }
    }
  }

  if (node.op == uuber::IrNodeOp::Guard && node.blocker_idx >= 0) {
    std::vector<int> blocker_sources =
        evaluator_gather_blocker_sources(ctx, node_idx);
    ids.insert(ids.end(), blocker_sources.begin(), blocker_sources.end());
  }

  sort_unique(ids);
  return ids;
}

struct ExactSharedStateStorage {
  TimeConstraintMap time_constraints;
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
};

template <typename PointVec>
inline bool exact_points_match_relevant_state(
    const PointVec &points, std::size_t lhs_idx, std::size_t rhs_idx,
    const std::vector<int> &relevant_source_ids,
    const uuber::NativeContext &ctx) {
  const ExactScenarioPointRef lhs = exact_point_ref(points, lhs_idx);
  const ExactScenarioPointRef rhs = exact_point_ref(points, rhs_idx);
  for (int source_id : relevant_source_ids) {
    if (filtered_forced_contains_id(*lhs.forced_complete_bits,
                                    lhs.forced_complete_bits_valid, source_id,
                                    ctx) !=
            filtered_forced_contains_id(*rhs.forced_complete_bits,
                                        rhs.forced_complete_bits_valid,
                                        source_id, ctx) ||
        filtered_forced_contains_id(*lhs.forced_survive_bits,
                                    lhs.forced_survive_bits_valid, source_id,
                                    ctx) !=
            filtered_forced_contains_id(*rhs.forced_survive_bits,
                                        rhs.forced_survive_bits_valid,
                                        source_id, ctx)) {
      return false;
    }
    const SourceTimeConstraint *lhs_constraint =
        time_constraints_find(lhs.time_constraints, source_id);
    const SourceTimeConstraint *rhs_constraint =
        time_constraints_find(rhs.time_constraints, source_id);
    const SourceTimeConstraint lhs_value =
        lhs_constraint ? *lhs_constraint : SourceTimeConstraint{};
    const SourceTimeConstraint rhs_value =
        rhs_constraint ? *rhs_constraint : SourceTimeConstraint{};
    if (!time_constraint_equal(lhs_value, rhs_value)) {
      return false;
    }
  }
  return true;
}

template <typename PointVec>
bool exact_points_share_relevant_state(
    const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, const uuber::NativeContext &ctx,
    std::size_t &anchor_idx);

template <typename PointVec>
inline bool exact_prepare_shared_state_storage(
    const uuber::NativeContext &ctx,
    const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, std::size_t &anchor_idx,
    ExactSharedStateStorage &storage) {
  anchor_idx = points.size();
  if (!exact_points_share_relevant_state(points, active_mask,
                                         relevant_source_ids, ctx,
                                         anchor_idx)) {
    return false;
  }
  const ExactScenarioPointRef anchor = exact_point_ref(points, anchor_idx);
  build_filtered_time_constraints(*anchor.time_constraints, relevant_source_ids,
                                  storage.time_constraints);
  build_filtered_forced_bits(
      ctx, *anchor.forced_complete_bits, anchor.forced_complete_bits_valid,
      relevant_source_ids, storage.forced_complete_bits,
      storage.forced_complete_bits_valid);
  build_filtered_forced_bits(
      ctx, *anchor.forced_survive_bits, anchor.forced_survive_bits_valid,
      relevant_source_ids, storage.forced_survive_bits,
      storage.forced_survive_bits_valid);
  return true;
}

template <typename PointVec, typename SharedFn>
bool exact_eval_shared_relevant_state_group(
    const uuber::NativeContext &ctx,
    const PointVec &points,
    const std::vector<std::uint8_t> *active_mask, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::size_t anchor_idx, const ExactSharedStateStorage &storage,
    SharedFn &&shared_eval,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  if (anchor_idx >= points.size()) {
    return false;
  }
  const ExactScenarioPointRef anchor = exact_point_ref(points, anchor_idx);
  NodeEvalState shared_state(
      ctx, anchor.t, component_idx, trial_params, trial_type_key, false, -1,
      &storage.time_constraints, nullptr,
      storage.forced_complete_bits_valid ? &storage.forced_complete_bits
                                         : nullptr,
      storage.forced_complete_bits_valid,
      storage.forced_survive_bits_valid ? &storage.forced_survive_bits
                                        : nullptr,
      storage.forced_survive_bits_valid, nullptr);
  shared_state.trial_params_soa =
      uniform_trial_params_soa != nullptr
          ? uniform_trial_params_soa
          : exact_trial_params_soa_for_point(ctx, anchor, trial_params,
                                             shared_state.trial_params_soa);
  if (!shared_state.trial_params_soa) {
    return false;
  }
  std::vector<const uuber::TrialParamsSoA *> point_trial_params_soa;
  if (uniform_trial_params_soa == nullptr) {
    std::size_t first_trial_param_mismatch = points.size();
    const uuber::TrialParamsSoA *mismatch_trial_params_soa = nullptr;
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (!exact_point_active(active_mask, i)) {
        continue;
      }
      const uuber::TrialParamsSoA *point_soa =
          exact_trial_params_soa_for_point(ctx, exact_point_ref(points, i),
                                           trial_params,
                                           shared_state.trial_params_soa);
      if (point_soa != shared_state.trial_params_soa) {
        first_trial_param_mismatch = i;
        mismatch_trial_params_soa = point_soa;
        break;
      }
    }
    if (first_trial_param_mismatch < points.size()) {
      point_trial_params_soa.assign(points.size(),
                                    shared_state.trial_params_soa);
      point_trial_params_soa[first_trial_param_mismatch] =
          mismatch_trial_params_soa;
      for (std::size_t i = first_trial_param_mismatch + 1; i < points.size();
           ++i) {
        if (!exact_point_active(active_mask, i)) {
          continue;
        }
        point_trial_params_soa[i] =
            exact_trial_params_soa_for_point(ctx, exact_point_ref(points, i),
                                             trial_params,
                                             shared_state.trial_params_soa);
      }
      shared_state.trial_params_soa_batch = &point_trial_params_soa;
    }
  }
  return shared_eval(shared_state, active_mask);
}

template <typename PointVec>
inline bool exact_partition_points_by_relevant_state(
    const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, const uuber::NativeContext &ctx,
    std::vector<std::vector<std::uint8_t>> &group_masks) {
  group_masks.clear();
  std::vector<std::size_t> anchor_indices;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (!exact_point_active(active_mask, i)) {
      continue;
    }
    bool assigned = false;
    for (std::size_t group_idx = 0; group_idx < anchor_indices.size();
         ++group_idx) {
      if (exact_points_match_relevant_state(points, anchor_indices[group_idx], i,
                                            relevant_source_ids, ctx)) {
        group_masks[group_idx][i] = 1u;
        assigned = true;
        break;
      }
    }
    if (assigned) {
      continue;
    }
    anchor_indices.push_back(i);
    group_masks.emplace_back(points.size(), 0u);
    group_masks.back()[i] = 1u;
  }
  return !group_masks.empty();
}

template <typename PointVec, typename SharedFn, typename FallbackFn>
bool exact_eval_with_shared_relevant_state(
    const uuber::NativeContext &ctx,
    const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    SharedFn &&shared_eval, FallbackFn &&fallback_eval,
    bool fallback_on_singleton_partition = false,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  if (exact_count_active_points(active_mask, points.size()) == 0u) {
    return true;
  }
  std::size_t anchor_idx = points.size();
  ExactSharedStateStorage storage;
  if (!exact_prepare_shared_state_storage(ctx, points, active_mask,
                                          relevant_source_ids, anchor_idx,
                                          storage)) {
    const std::size_t active_point_count =
        exact_count_active_points(active_mask, points.size());
    std::vector<std::vector<std::uint8_t>> group_masks;
    if (!exact_partition_points_by_relevant_state(points, active_mask,
                                                  relevant_source_ids, ctx,
                                                  group_masks) ||
        group_masks.size() <= 1u) {
      return fallback_eval();
    }
    if (fallback_on_singleton_partition &&
        group_masks.size() >= active_point_count) {
      return fallback_eval();
    }
    record_unified_outcome_exact_shared_state_partition(
        static_cast<std::uint64_t>(active_point_count),
        static_cast<std::uint64_t>(group_masks.size()));
    for (const std::vector<std::uint8_t> &group_mask : group_masks) {
      std::size_t group_anchor_idx = points.size();
      ExactSharedStateStorage group_storage;
      if (!exact_prepare_shared_state_storage(
              ctx, points, &group_mask, relevant_source_ids, group_anchor_idx,
              group_storage) ||
          !exact_eval_shared_relevant_state_group(
              ctx, points, &group_mask, component_idx, trial_params,
              trial_type_key, group_anchor_idx, group_storage, shared_eval,
              uniform_trial_params_soa)) {
        return fallback_eval();
      }
    }
    return true;
  }
  return exact_eval_shared_relevant_state_group(
      ctx, points, active_mask, component_idx, trial_params, trial_type_key,
      anchor_idx, storage, shared_eval, uniform_trial_params_soa);
}

inline bool exact_point_active(const std::vector<std::uint8_t> *active_mask,
                               std::size_t idx) {
  return !active_mask ||
         (idx < active_mask->size() && (*active_mask)[idx] != 0u);
}

inline std::size_t exact_count_active_points(
    const std::vector<std::uint8_t> *active_mask, std::size_t point_count) {
  if (active_mask == nullptr) {
    return point_count;
  }
  std::size_t active_count = 0u;
  const std::size_t limit = std::min(point_count, active_mask->size());
  for (std::size_t i = 0; i < limit; ++i) {
    active_count += ((*active_mask)[i] != 0u) ? 1u : 0u;
  }
  return active_count;
}

thread_local int g_exact_node_batch_depth = 0;

struct ExactNodeBatchDepthGuard {
  ExactNodeBatchDepthGuard() { ++g_exact_node_batch_depth; }
  ~ExactNodeBatchDepthGuard() { --g_exact_node_batch_depth; }

  bool recursive_call() const noexcept { return g_exact_node_batch_depth > 1; }
};

inline void exact_init_batch_values(std::size_t point_count,
                                    uuber::KernelNodeBatchValues &out_values) {
  out_values.density.assign(point_count, 0.0);
  out_values.survival.assign(point_count, 1.0);
  out_values.cdf.assign(point_count, 0.0);
}

inline bool exact_batch_values_size_matches(
    const uuber::KernelNodeBatchValues &values, std::size_t point_count) {
  return values.density.size() == point_count &&
         values.survival.size() == point_count &&
         values.cdf.size() == point_count;
}

inline void exact_merge_batch_values_subset(
    const uuber::KernelNodeBatchValues &src_values,
    const std::vector<std::uint8_t> *active_mask,
    uuber::KernelNodeBatchValues &dst_values) {
  const std::size_t point_count = dst_values.density.size();
  if (!exact_batch_values_size_matches(src_values, point_count) ||
      !exact_batch_values_size_matches(dst_values, point_count)) {
    return;
  }
  if (active_mask == nullptr) {
    dst_values = src_values;
    return;
  }
  for (std::size_t i = 0; i < point_count; ++i) {
    if (!exact_point_active(active_mask, i)) {
      continue;
    }
    dst_values.density[i] = src_values.density[i];
    dst_values.survival[i] = src_values.survival[i];
    dst_values.cdf[i] = src_values.cdf[i];
  }
}

inline EvalNeed exact_eval_need_from_kernel_need(
    const uuber::KernelEvalNeed &kneed) {
  EvalNeed need = static_cast<EvalNeed>(0u);
  if (kneed.density) {
    need = need | EvalNeed::kDensity;
  }
  if (kneed.survival) {
    need = need | EvalNeed::kSurvival;
  }
  if (kneed.cdf) {
    need = need | EvalNeed::kCDF;
  }
  return need;
}

inline ForcedStateView exact_make_point_forced_state(
    const uuber::NativeContext &ctx, const ExactScenarioPointRef &point) {
  return make_forced_state_view(
      nullptr, point.forced_complete_bits_valid ? point.forced_complete_bits
                                                : nullptr,
      &point.forced_complete_bits_valid,
      point.forced_survive_bits_valid ? point.forced_survive_bits : nullptr,
      &point.forced_survive_bits_valid,
      ctx.ir.label_id_to_bit_idx.empty() ? nullptr : &ctx.ir.label_id_to_bit_idx);
}

inline const uuber::TrialParamsSoA *exact_trial_params_soa_for_point(
    const uuber::NativeContext &ctx, const ExactScenarioPointRef &point,
    const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *fallback_trial_params_soa) {
  if (point.trial_params_soa != nullptr) {
    return point.trial_params_soa;
  }
  if (fallback_trial_params_soa != nullptr) {
    return fallback_trial_params_soa;
  }
  return resolve_trial_params_soa(ctx, trial_params);
}

inline bool exact_guard_batch_inputs_valid(
    const std::vector<double> &batch_times, std::size_t point_count,
    const uuber::KernelNodeBatchValues &reference_values,
    const uuber::KernelNodeBatchValues &blocker_values) {
  return batch_times.size() == point_count &&
         reference_values.density.size() == point_count &&
         reference_values.survival.size() == point_count &&
         reference_values.cdf.size() == point_count &&
         blocker_values.density.size() == point_count &&
         blocker_values.survival.size() == point_count &&
         blocker_values.cdf.size() == point_count;
}

inline bool exact_point_active(const std::vector<std::uint8_t> *active_mask,
                               std::size_t idx);
inline bool exact_batch_values_size_matches(
    const uuber::KernelNodeBatchValues &values, std::size_t point_count);

inline void exact_seed_guard_density(
    const uuber::KernelNodeBatchValues &reference_values,
    const uuber::KernelNodeBatchValues &blocker_values,
    uuber::KernelNodeBatchValues &guard_out,
    const std::vector<std::uint8_t> *active_mask) {
  const std::size_t point_count = guard_out.density.size();
  for (std::size_t i = 0; i < point_count; ++i) {
    if (!exact_point_active(active_mask, i)) {
      continue;
    }
    const double ref_density = safe_density(reference_values.density[i]);
    const double blocker_survival =
        clamp_probability(blocker_values.survival[i]);
    guard_out.density[i] = safe_density(ref_density * blocker_survival);
  }
}

inline bool exact_fill_guard_prepared_batch(
    const uuber::NativeContext &ctx, const std::vector<double> &batch_times,
    int component_idx, const TrialParamSet *trial_params,
    NodeEvalState &shared_state, const uuber::KernelOp &op,
    const std::vector<std::uint8_t> *group_mask,
    uuber::KernelNodeBatchValues &guard_out) {
  const std::size_t point_count = batch_times.size();
  if (!exact_batch_values_size_matches(guard_out, point_count)) {
    return false;
  }
  std::vector<const uuber::TrialParamsSoA *> subgroup_params;
  std::vector<std::vector<std::size_t>> subgroup_indices;
  subgroup_params.reserve(point_count);
  subgroup_indices.reserve(point_count);
  for (std::size_t i = 0; i < point_count; ++i) {
    if (!exact_point_active(group_mask, i)) {
      continue;
    }
    const uuber::TrialParamsSoA *point_soa =
        (shared_state.trial_params_soa_batch != nullptr &&
         i < shared_state.trial_params_soa_batch->size())
            ? (*shared_state.trial_params_soa_batch)[i]
            : shared_state.trial_params_soa;
    std::size_t subgroup_idx = subgroup_params.size();
    for (std::size_t j = 0; j < subgroup_params.size(); ++j) {
      if (subgroup_params[j] == point_soa) {
        subgroup_idx = j;
        break;
      }
    }
    if (subgroup_idx == subgroup_params.size()) {
      subgroup_params.push_back(point_soa);
      subgroup_indices.emplace_back();
    }
    subgroup_indices[subgroup_idx].push_back(i);
  }
  for (std::size_t subgroup_idx = 0; subgroup_idx < subgroup_params.size();
       ++subgroup_idx) {
    const std::vector<std::size_t> &indices = subgroup_indices[subgroup_idx];
    if (indices.empty()) {
      continue;
    }
    std::vector<double> subgroup_times(indices.size(), 0.0);
    for (std::size_t j = 0; j < indices.size(); ++j) {
      subgroup_times[j] = batch_times[indices[j]];
    }
    GuardEvalInput guard_input = evaluator_make_guard_input(
        ctx, op.node_idx, component_idx, &shared_state.trial_type_key,
        trial_params, subgroup_params[subgroup_idx], &shared_state.time_constraints,
        nullptr, nullptr, nullptr, nullptr, &shared_state.forced_state);
    std::vector<double> subgroup_cdf;
    if (!evaluator_guard_cdf_batch_prepared(guard_input, subgroup_times,
                                            subgroup_cdf) ||
        subgroup_cdf.size() != indices.size()) {
      return false;
    }
    record_unified_outcome_exact_guard_prepared_batch_call(
        static_cast<std::uint64_t>(indices.size()));
    for (std::size_t j = 0; j < indices.size(); ++j) {
      const std::size_t point_idx = indices[j];
      guard_out.cdf[point_idx] = clamp_probability(subgroup_cdf[j]);
      guard_out.survival[point_idx] =
          clamp_probability(1.0 - guard_out.cdf[point_idx]);
    }
  }
  return true;
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa);
bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa);

template <typename PointVec>
bool exact_eval_node_batch_common(
    const uuber::NativeContext &ctx, int node_idx, const PointVec &points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, EvalNeed need,
    uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa);

template <typename PointVec, typename DirectEvalFn>
bool exact_eval_event_batch_common(
    const uuber::NativeContext &ctx, const PointVec &points,
    const std::vector<double> &batch_times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    int event_idx,
    const uuber::KernelEvalNeed &kneed,
    uuber::KernelNodeBatchValues &event_out,
    const std::vector<std::uint8_t> *active_mask, DirectEvalFn &&direct_eval,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  const std::size_t point_count = points.size();
  if (batch_times.size() != point_count ||
      event_idx < 0 ||
      event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return false;
  }
  const uuber::IrEvent &event =
      ctx.ir.events[static_cast<std::size_t>(event_idx)];
  exact_init_batch_values(point_count, event_out);
  if (exact_count_active_points(active_mask, point_count) == 0u) {
    return true;
  }
  const EvalNeed event_need = exact_eval_need_from_kernel_need(kneed);

  if (direct_eval(event, event_need, event_out)) {
    return true;
  }
  if (event.node_idx < 0) {
    return false;
  }
  const uuber::IrNode &event_node = ir_node_required(ctx, event.node_idx);
  if (event_node.op == uuber::IrNodeOp::EventAcc) {
    return exact_eval_node_batch_with_points(
        ctx, event.node_idx, points, component_idx, trial_params,
        trial_type_key, event_need, event_out, active_mask,
        uniform_trial_params_soa);
  }
  if (event_node.op == uuber::IrNodeOp::EventPool) {
    const std::vector<int> event_source_ids =
        exact_relevant_source_ids_for_node(ctx, event.node_idx);
    return exact_eval_with_shared_relevant_state(
        ctx, points, active_mask, event_source_ids, component_idx,
        trial_params, trial_type_key,
        [&](NodeEvalState &shared_state,
            const std::vector<std::uint8_t> *group_mask) -> bool {
          uuber::KernelEventBatchEvalFn shared_event_eval =
              evaluator_make_kernel_event_eval_batch(shared_state);
          uuber::KernelNodeBatchValues shared_event_out;
          if (!shared_event_eval(event_idx, batch_times, kneed,
                                 shared_event_out)) {
            return false;
          }
          exact_merge_batch_values_subset(shared_event_out, group_mask,
                                          event_out);
          return true;
        },
      [&]() -> bool {
        Rcpp::stop(
            "Exact batch invariant failed for EventPool node: "
            "event_idx=%d node_idx=%d active_points=%d",
            event_idx, event.node_idx,
              static_cast<int>(exact_count_active_points(active_mask,
                                                         point_count)));
          return false;
        },
        false, uniform_trial_params_soa);
  }
  return false;
}

template <typename PointVec>
bool exact_eval_guard_batch_common(
    const uuber::NativeContext &ctx, const PointVec &points,
    const std::vector<double> &batch_times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const uuber::KernelOp &op,
    const uuber::KernelNodeBatchValues &reference_values,
    const uuber::KernelNodeBatchValues &blocker_values,
    const uuber::KernelEvalNeed &kneed,
    uuber::KernelNodeBatchValues &guard_out,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  const std::size_t point_count = points.size();
  if (!exact_guard_batch_inputs_valid(batch_times, point_count, reference_values,
                                      blocker_values)) {
    return false;
  }
  exact_init_batch_values(point_count, guard_out);
  if (kneed.density) {
    exact_seed_guard_density(reference_values, blocker_values, guard_out,
                             active_mask);
  }
  if (exact_count_active_points(active_mask, point_count) == 0u) {
    return true;
  }
  if (op.node_idx < 0) {
    return false;
  }
  const std::vector<int> guard_source_ids =
      exact_relevant_source_ids_for_node(ctx, op.node_idx);
  return exact_eval_with_shared_relevant_state(
      ctx, points, active_mask, guard_source_ids, component_idx,
      trial_params, trial_type_key,
      [&](NodeEvalState &shared_state,
          const std::vector<std::uint8_t> *group_mask) -> bool {
        return exact_fill_guard_prepared_batch(
            ctx, batch_times, component_idx, trial_params, shared_state, op,
            group_mask, guard_out);
      },
      [&]() -> bool {
        Rcpp::stop(
            "Exact batch invariant failed for Guard node: "
            "node_idx=%d active_points=%d",
            op.node_idx,
            static_cast<int>(exact_count_active_points(active_mask,
                                                       point_count)));
        return false;
      },
      false, uniform_trial_params_soa);
}

template <typename PointVec>
bool exact_eval_node_batch_incremental_common(
    const uuber::NativeContext &ctx, int node_idx, const PointVec &points,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  uuber::KernelBatchRuntimeState batch_runtime;
  uuber::KernelEventBatchEvalFn event_eval_batch =
      [&](int event_idx, const std::vector<double> &batch_times,
          const uuber::KernelEvalNeed &kneed,
          uuber::KernelNodeBatchValues &event_out) -> bool {
    return exact_eval_event_batch_common(
        ctx, points, batch_times, component_idx, trial_params, trial_type_key,
        event_idx, kneed, event_out, active_mask,
        [&](const uuber::IrEvent &event, EvalNeed event_need,
            uuber::KernelNodeBatchValues &out) -> bool {
          return exact_eval_simple_acc_event_batch(
              ctx, event, points, component_idx, trial_params, trial_type_key,
              event_need, out, active_mask, uniform_trial_params_soa);
        },
        uniform_trial_params_soa);
  };

  uuber::KernelGuardBatchEvalFn guard_eval_batch =
      [&](const uuber::KernelOp &op, const std::vector<double> &batch_times,
          const uuber::KernelNodeBatchValues &reference_values,
          const uuber::KernelNodeBatchValues &blocker_values,
          const uuber::KernelEvalNeed &kneed,
          uuber::KernelNodeBatchValues &guard_out) -> bool {
    return exact_eval_guard_batch_common(
        ctx, points, batch_times, component_idx, trial_params, trial_type_key,
        op, reference_values, blocker_values, kneed, guard_out, active_mask,
        uniform_trial_params_soa);
  };

  uuber::KernelEvalNeed kernel_need;
  kernel_need.density = needs_density(need);
  kernel_need.survival = needs_survival(need);
  kernel_need.cdf = needs_cdf(need);
  if (!kernel_need.density && !kernel_need.survival && !kernel_need.cdf) {
    kernel_need.density = true;
    kernel_need.survival = true;
    kernel_need.cdf = true;
  }
  return uuber::eval_kernel_node_batch_incremental(
      ctx.kernel_program, batch_runtime, node_idx, times, kernel_need,
      event_eval_batch, guard_eval_batch, out_values);
}

template <typename PointVec>
bool exact_points_share_relevant_state(
    const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, const uuber::NativeContext &ctx,
    std::size_t &anchor_idx) {
  anchor_idx = points.size();
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (active_mask && (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    anchor_idx = i;
    break;
  }
  if (anchor_idx >= points.size()) {
    return false;
  }
  for (std::size_t i = anchor_idx + 1u; i < points.size(); ++i) {
    if (active_mask && (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    if (!exact_points_match_relevant_state(points, anchor_idx, i,
                                           relevant_source_ids, ctx)) {
      return false;
    }
  }
  return true;
}

template <typename PointVec>
bool exact_eval_node_batch_common(
    const uuber::NativeContext &ctx, int node_idx,
    const PointVec &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  if (points.empty()) {
    out_values = uuber::KernelNodeBatchValues{};
    return true;
  }
  ExactNodeBatchDepthGuard depth_guard;
  const std::size_t active_point_count =
      exact_count_active_points(active_mask, points.size());
  record_unified_outcome_exact_node_batch_call(
      static_cast<std::uint64_t>(active_point_count),
      depth_guard.recursive_call());

  std::vector<double> times(points.size(), 0.0);
  for (std::size_t i = 0; i < points.size(); ++i) {
    times[i] = exact_point_ref(points, i).t;
  }

  if (node_idx >= 0 && node_idx < static_cast<int>(ctx.ir.nodes.size())) {
    const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
    if ((node.op == uuber::IrNodeOp::EventAcc ||
         node.op == uuber::IrNodeOp::EventPool) &&
        node.event_idx >= 0 &&
        node.event_idx < static_cast<int>(ctx.ir.events.size())) {
      const uuber::IrEvent &event =
          ctx.ir.events[static_cast<std::size_t>(node.event_idx)];
      if (exact_eval_simple_acc_event_batch(
              ctx, event, points, component_idx, trial_params, trial_type_key,
              need, out_values, active_mask, uniform_trial_params_soa)) {
        return true;
      }
      const std::vector<int> event_source_ids =
          exact_relevant_source_ids_for_node(ctx, node_idx);
      exact_init_batch_values(points.size(), out_values);
      return exact_eval_with_shared_relevant_state(
          ctx, points, active_mask, event_source_ids, component_idx,
          trial_params, trial_type_key,
          [&](NodeEvalState &shared_state,
              const std::vector<std::uint8_t> *group_mask) -> bool {
            uuber::KernelBatchRuntimeState batch_runtime;
            uuber::KernelNodeBatchValues shared_out;
            if (!evaluator_eval_node_batch_with_state_dense(
                    node_idx, times, shared_state, need, batch_runtime,
                    shared_out)) {
              return false;
            }
            exact_merge_batch_values_subset(shared_out, group_mask, out_values);
            return true;
          },
          [&]() -> bool {
            Rcpp::stop(
                "Exact batch invariant failed for event node: "
                "node_idx=%d event_idx=%d active_points=%d",
                node_idx, node.event_idx,
                static_cast<int>(active_point_count));
            return false;
          },
          false, uniform_trial_params_soa);
    }
    if (node.op == uuber::IrNodeOp::Guard &&
        node.reference_idx >= 0 &&
        node.blocker_idx >= 0) {
      uuber::KernelNodeBatchValues reference_values;
      if (!exact_eval_node_batch_with_points(
              ctx, node.reference_idx, points, component_idx, trial_params,
              trial_type_key, EvalNeed::kDensity, reference_values,
              active_mask, uniform_trial_params_soa)) {
        return false;
      }
      uuber::KernelNodeBatchValues blocker_values;
      if (!exact_eval_node_batch_with_points(
              ctx, node.blocker_idx, points, component_idx, trial_params,
              trial_type_key, EvalNeed::kSurvival, blocker_values,
              active_mask, uniform_trial_params_soa)) {
        return false;
      }
      uuber::KernelOp guard_op;
      guard_op.code = uuber::KernelOpCode::Guard;
      guard_op.node_idx = node_idx;
      guard_op.flags = node.flags;
      uuber::KernelEvalNeed kernel_need;
      kernel_need.density = needs_density(need);
      kernel_need.survival = needs_survival(need);
      kernel_need.cdf = needs_cdf(need);
      if (!kernel_need.density && !kernel_need.survival && !kernel_need.cdf) {
        kernel_need.density = true;
        kernel_need.survival = true;
        kernel_need.cdf = true;
      }
      return exact_eval_guard_batch_common(
          ctx, points, times, component_idx, trial_params, trial_type_key,
          guard_op, reference_values, blocker_values, kernel_need,
          out_values, active_mask, uniform_trial_params_soa);
    }
  }

  const std::vector<int> relevant_source_ids =
      exact_relevant_source_ids_for_node(ctx, node_idx);
  if (node_idx >= 0 && node_idx < static_cast<int>(ctx.ir.nodes.size())) {
    const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
    if (node.op == uuber::IrNodeOp::And ||
        node.op == uuber::IrNodeOp::Or ||
        node.op == uuber::IrNodeOp::Not) {
      std::size_t anchor_idx = points.size();
      if (!exact_points_share_relevant_state(points, active_mask,
                                             relevant_source_ids, ctx,
                                             anchor_idx)) {
        exact_init_batch_values(points.size(), out_values);
        return exact_eval_node_batch_incremental_common(
            ctx, node_idx, points, times, component_idx, trial_params,
            trial_type_key, need, out_values, active_mask,
            uniform_trial_params_soa);
      }
    }
  }
  exact_init_batch_values(points.size(), out_values);
  return exact_eval_with_shared_relevant_state(
      ctx, points, active_mask, relevant_source_ids, component_idx,
      trial_params, trial_type_key,
      [&](NodeEvalState &shared_state,
          const std::vector<std::uint8_t> *group_mask) -> bool {
        uuber::KernelBatchRuntimeState batch_runtime;
        uuber::KernelNodeBatchValues shared_out;
        if (!evaluator_eval_node_batch_with_state_dense(
                node_idx, times, shared_state, need, batch_runtime,
                shared_out)) {
          return false;
        }
        exact_merge_batch_values_subset(shared_out, group_mask, out_values);
        return true;
      },
      [&]() -> bool {
        const uuber::IrNode &node = ir_node_required(ctx, node_idx);
        Rcpp::stop(
            "Exact batch invariant failed for generic node: "
            "node_idx=%d op=%d event_idx=%d active_points=%d",
            node_idx, static_cast<int>(node.op), node.event_idx,
            static_cast<int>(active_point_count));
        return false;
      },
      true, uniform_trial_params_soa);
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_eval_node_batch_common(ctx, node_idx, points, component_idx,
                                      trial_params, trial_type_key, need,
                                      out_values, active_mask,
                                      uniform_trial_params_soa);
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_eval_node_batch_common(ctx, node_idx, points, component_idx,
                                      trial_params, trial_type_key, need,
                                      out_values, active_mask,
                                      uniform_trial_params_soa);
}

bool exact_collect_scenarios_batch_aligned_from_points_impl(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, ExactScenarioBatch &scenario_batch,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  scenario_batch.clear();
  const RankedNodeTransitionPlan &plan = compiler.plan_for_node(node_idx);
  if (!plan.valid || plan.transitions.empty()) {
    return false;
  }

  std::vector<ExactScenarioPoint> aligned_points;
  std::vector<std::uint8_t> aligned_active_mask;
  if (exact_collect_deterministic_scenarios_batch_from_points_internal(
          compiler, ctx, node_idx, seed_points, component_idx, trial_params,
          trial_type_key, aligned_points, aligned_active_mask,
          uniform_trial_params_soa)) {
    scenario_batch.reserve(seed_points.size());
    for (std::size_t i = 0; i < aligned_points.size(); ++i) {
      if (i < aligned_active_mask.size() && aligned_active_mask[i] != 0u &&
          std::isfinite(aligned_points[i].weight) &&
          aligned_points[i].weight > 0.0) {
        exact_scenario_batch_append_point(scenario_batch,
                                          std::move(aligned_points[i]));
      }
    }
    return !scenario_batch.empty();
  }

  scenario_batch.reserve(seed_points.size() * plan.transitions.size());
  ExactScenarioBatch transition_batch;
  transition_batch.reserve(seed_points.size());
  std::vector<std::uint8_t> active;
  active.reserve(seed_points.size());
  std::vector<double> transition_weights;
  transition_weights.reserve(seed_points.size());
  uuber::KernelNodeBatchValues step_values;
  uuber::KernelNodeBatchValues blocker_values;
  for (const RankedTransitionTemplate &transition : plan.transitions) {
    transition_batch.clear();
    active.assign(seed_points.size(), 0u);
    transition_weights.assign(seed_points.size(), 0.0);
    std::size_t active_count = 0u;
    for (std::size_t i = 0; i < seed_points.size(); ++i) {
      const ExactScenarioPoint &point = seed_points[i];
      const bool is_active = std::isfinite(point.t) && point.t >= 0.0 &&
                             std::isfinite(point.weight) && point.weight > 0.0;
      active[i] = is_active ? 1u : 0u;
      transition_weights[i] = is_active ? point.weight : 0.0;
      active_count += is_active ? 1u : 0u;
    }
    if (active_count == 0u) {
      continue;
    }

    bool points_materialized = false;
    auto materialize_transition_batch = [&]() {
      if (points_materialized) {
        return;
      }
      transition_batch.reserve(seed_points.size());
      for (std::size_t i = 0; i < seed_points.size(); ++i) {
        exact_scenario_batch_append_point(transition_batch, seed_points[i],
                                          transition_weights[i]);
      }
      points_materialized = true;
    };

    bool transition_valid = true;
    auto deactivate_point = [&](std::size_t idx) {
      if (idx < active.size() && active[idx] != 0u) {
        active[idx] = 0u;
        if (active_count > 0u) {
          --active_count;
        }
      }
      if (idx < transition_weights.size()) {
        transition_weights[idx] = 0.0;
      }
      if (points_materialized && idx < transition_batch.weight.size()) {
        transition_batch.weight[idx] = 0.0;
      }
    };
    for (const RankedTransitionStep &step : transition.steps) {
      if (step.kind == RankedTransitionStepKind::EvalDensityNode ||
          step.kind == RankedTransitionStepKind::EvalCDFNode) {
        const EvalNeed step_need =
            (step.kind == RankedTransitionStepKind::EvalDensityNode)
                ? EvalNeed::kDensity
                : EvalNeed::kCDF;
        const bool step_ok =
            points_materialized
                ? exact_eval_node_batch_with_points(
                      ctx, step.node_idx, transition_batch, component_idx,
                      trial_params, trial_type_key, step_need, step_values,
                      &active, uniform_trial_params_soa)
                : exact_eval_node_batch_with_points(
                      ctx, step.node_idx, seed_points, component_idx,
                      trial_params, trial_type_key, step_need, step_values,
                      &active, uniform_trial_params_soa);
        if (!step_ok ||
            !exact_batch_values_size_matches(step_values, seed_points.size())) {
          transition_valid = false;
          break;
        }
        for (std::size_t i = 0; i < seed_points.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          const double factor =
              (step.kind == RankedTransitionStepKind::EvalDensityNode)
                  ? safe_density(step_values.density[i])
                  : clamp_probability(step_values.cdf[i]);
          if (!std::isfinite(factor) || factor <= 0.0) {
            deactivate_point(i);
            continue;
          }
          transition_weights[i] *= factor;
          if (points_materialized && i < transition_batch.weight.size()) {
            transition_batch.weight[i] = transition_weights[i];
          }
          if (!std::isfinite(transition_weights[i]) ||
              transition_weights[i] <= 0.0) {
            deactivate_point(i);
          }
        }
      } else if (step.kind == RankedTransitionStepKind::EvalGuardEffective) {
        const uuber::IrNode &guard_node = ir_node_required(ctx, step.node_idx);
        if (guard_node.blocker_idx < 0) {
          transition_valid = false;
          break;
        }
        const bool blocker_ok =
            points_materialized
                ? exact_eval_node_batch_with_points(
                      ctx, guard_node.blocker_idx, transition_batch,
                      component_idx, trial_params, trial_type_key,
                      EvalNeed::kSurvival, blocker_values, &active,
                      uniform_trial_params_soa)
                : exact_eval_node_batch_with_points(
                      ctx, guard_node.blocker_idx, seed_points, component_idx,
                      trial_params, trial_type_key, EvalNeed::kSurvival,
                      blocker_values, &active, uniform_trial_params_soa);
        if (!blocker_ok ||
            !exact_batch_values_size_matches(blocker_values,
                                             seed_points.size())) {
          transition_valid = false;
          break;
        }
        for (std::size_t i = 0; i < seed_points.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          const double factor = clamp_probability(blocker_values.survival[i]);
          if (!std::isfinite(factor) || factor <= 0.0) {
            deactivate_point(i);
            continue;
          }
          transition_weights[i] *= factor;
          if (points_materialized && i < transition_batch.weight.size()) {
            transition_batch.weight[i] = transition_weights[i];
          }
          if (!std::isfinite(transition_weights[i]) ||
              transition_weights[i] <= 0.0) {
            deactivate_point(i);
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddCompleteSources) {
        materialize_transition_batch();
        for (std::size_t i = 0; i < transition_batch.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          bool ok = true;
          for (int id : step.source_ids) {
            if (!time_constraints_mark_complete(
                    id, transition_batch.t[i], step.bind_exact_current_time,
                    transition_batch.time_constraints[i])) {
              ok = false;
              break;
            }
          }
          if (!ok) {
            deactivate_point(i);
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddSurviveSources) {
        materialize_transition_batch();
        for (std::size_t i = 0; i < transition_batch.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          bool ok = true;
          for (int id : step.source_ids) {
            if (!time_constraints_mark_survive(
                    id, transition_batch.t[i],
                    transition_batch.time_constraints[i])) {
              ok = false;
              break;
            }
          }
          if (!ok) {
            deactivate_point(i);
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddOrWitnessFromSources) {
        materialize_transition_batch();
        for (std::size_t i = 0; i < transition_batch.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          int witness = NA_INTEGER;
          bool all_forced = true;
          for (std::size_t src_i = 0; src_i < step.source_ids.size(); ++src_i) {
            const int id = step.source_ids[src_i];
            if (id == NA_INTEGER || id < 0) {
              continue;
            }
            const int bit_idx =
                (src_i < step.source_bits.size()) ? step.source_bits[src_i] : -1;
            const bool is_forced =
                ((i < transition_batch.forced_complete_bits_valid.size()) &&
                 transition_batch.forced_complete_bits_valid[i] != 0u &&
                 bit_idx >= 0 &&
                 transition_batch.forced_complete_bits[i].test(bit_idx)) ||
                time_constraints_contains_complete_at(
                    &transition_batch.time_constraints[i], id,
                    transition_batch.t[i]);
            if (!is_forced) {
              witness = id;
              all_forced = false;
              break;
            }
          }
          if (all_forced || witness == NA_INTEGER ||
              !time_constraints_mark_survive(
                  witness, transition_batch.t[i],
                  transition_batch.time_constraints[i])) {
            deactivate_point(i);
          }
        }
      }

      if (active_count == 0u) {
        transition_valid = false;
        break;
      }
    }

    if (!transition_valid) {
      continue;
    }
    if (points_materialized) {
      for (std::size_t i = 0; i < transition_batch.size(); ++i) {
        if (active[i] == 0u || !std::isfinite(transition_weights[i]) ||
            transition_weights[i] <= 0.0) {
          continue;
        }
        exact_scenario_batch_append_row(scenario_batch, transition_batch, i,
                                        transition_weights[i]);
      }
    } else {
      for (std::size_t i = 0; i < seed_points.size(); ++i) {
        if (active[i] == 0u || !std::isfinite(transition_weights[i]) ||
            transition_weights[i] <= 0.0) {
          continue;
        }
        exact_scenario_batch_append_point(scenario_batch, seed_points[i],
                                          transition_weights[i]);
      }
    }
  }

  return !scenario_batch.empty();
}

bool exact_collect_deterministic_scenarios_batch_from_points_internal(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &aligned_points,
    std::vector<std::uint8_t> &active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  aligned_points.clear();
  active_mask.clear();
  const RankedNodeTransitionPlan &plan = compiler.plan_for_node(node_idx);
  if (!plan.valid || plan.transitions.size() != 1u) {
    return false;
  }

  const RankedTransitionTemplate &transition = plan.transitions.front();
  if (transition.steps.empty() ||
      transition.steps.front().kind != RankedTransitionStepKind::EvalDensityNode) {
    return false;
  }
  for (std::size_t step_idx = 1; step_idx < transition.steps.size(); ++step_idx) {
    if (transition.steps[step_idx].kind !=
            RankedTransitionStepKind::AddCompleteSources &&
        transition.steps[step_idx].kind !=
            RankedTransitionStepKind::AddSurviveSources) {
      return false;
    }
  }

  uuber::KernelNodeBatchValues step_values;
  if (!exact_eval_node_batch_with_points(
          ctx, transition.steps.front().node_idx, seed_points, component_idx,
          trial_params, trial_type_key, EvalNeed::kDensity, step_values, nullptr,
          uniform_trial_params_soa) ||
      step_values.density.size() != seed_points.size()) {
    return false;
  }

  aligned_points = seed_points;
  active_mask.assign(seed_points.size(), 0u);
  for (std::size_t i = 0; i < seed_points.size(); ++i) {
    const ExactScenarioPoint &seed_point = seed_points[i];
    if (!(std::isfinite(seed_point.t) && seed_point.t >= 0.0 &&
          std::isfinite(seed_point.weight) && seed_point.weight > 0.0)) {
      continue;
    }
    ExactScenarioPoint &point = aligned_points[i];
    const double factor = safe_density(step_values.density[i]);
    if (!std::isfinite(factor) || factor <= 0.0) {
      point.weight = 0.0;
      continue;
    }
    point.weight *= factor;
    if (!std::isfinite(point.weight) || point.weight <= 0.0) {
      point.weight = 0.0;
      continue;
    }
    bool valid = true;
    for (std::size_t step_idx = 1; step_idx < transition.steps.size(); ++step_idx) {
      const RankedTransitionStep &step = transition.steps[step_idx];
      if (step.kind == RankedTransitionStepKind::AddCompleteSources) {
        for (int id : step.source_ids) {
          if (!time_constraints_mark_complete(
                  id, point.t, step.bind_exact_current_time,
                  point.time_constraints)) {
            valid = false;
            break;
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddSurviveSources) {
        for (int id : step.source_ids) {
          if (!time_constraints_mark_survive(id, point.t,
                                             point.time_constraints)) {
            valid = false;
            break;
          }
        }
      }
      if (!valid) {
        break;
      }
    }
    if (valid) {
      active_mask[i] = 1u;
    } else {
      point.weight = 0.0;
    }
  }
  return true;
}

bool exact_collect_scenarios_batch_impl(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    ExactScenarioBatch &scenario_points) {
  std::vector<ExactScenarioPoint> seed_points(times.size());
  for (std::size_t i = 0; i < times.size(); ++i) {
    ExactScenarioPoint &point = seed_points[i];
    point.t = times[i];
    point.density_index = i;
    point.weight =
        (std::isfinite(times[i]) && times[i] >= 0.0) ? 1.0 : 0.0;
    if (state.trial_params_soa_batch != nullptr &&
        i < state.trial_params_soa_batch->size()) {
      point.trial_params_soa = (*state.trial_params_soa_batch)[i];
    } else {
      point.trial_params_soa = state.trial_params_soa;
    }
    point.forced_complete_bits = state.forced_complete_bits;
    point.forced_survive_bits = state.forced_survive_bits;
    point.forced_complete_bits_valid = state.forced_complete_bits_valid;
    point.forced_survive_bits_valid = state.forced_survive_bits_valid;
    point.time_constraints = state.time_constraints;
  }
  const uuber::TrialParamsSoA *uniform_trial_params_soa =
      state.trial_params_soa_batch == nullptr ? state.trial_params_soa
                                              : nullptr;
  return exact_collect_scenarios_batch_aligned_from_points(
      compiler, ctx, node_idx, seed_points, state.component_idx,
      state.trial_params, state.trial_type_key, scenario_points,
      uniform_trial_params_soa);
}

inline std::unordered_map<std::uint64_t,
                          std::unique_ptr<RankedTransitionCompiler>> &
exact_scenario_compiler_cache_registry() {
  thread_local std::unordered_map<std::uint64_t,
                                  std::unique_ptr<RankedTransitionCompiler>>
      cache;
  return cache;
}

RankedTransitionCompiler &
exact_scenario_compiler_for_ctx(const uuber::NativeContext &ctx) {
  std::unique_ptr<RankedTransitionCompiler> &entry =
      exact_scenario_compiler_cache_registry()[ctx.runtime_cache_instance_id];
  if (!entry) {
    entry = std::make_unique<RankedTransitionCompiler>(ctx);
  }
  return *entry;
}

inline bool exact_competitor_guard_fastpath_eligible(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache) {
  if (competitor_cache.compiled_ops.empty()) {
    return false;
  }
  for (const uuber::CompetitorCompiledOp &op : competitor_cache.compiled_ops) {
    if (op.target_node_indices.size() != 1u) {
      return false;
    }
    const int node_idx = op.target_node_indices.front();
    if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      return false;
    }
    const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
    if (node.op != uuber::IrNodeOp::Guard) {
      return false;
    }
  }
  return true;
}

inline bool exact_competitor_single_node_fastpath_eligible(
    const uuber::CompetitorClusterCacheEntry &competitor_cache) {
  if (competitor_cache.mutates_forced_survive ||
      competitor_cache.compiled_ops.size() != 1u) {
    return false;
  }
  const uuber::CompetitorCompiledOp &op = competitor_cache.compiled_ops.front();
  return op.target_node_indices.size() == 1u &&
         op.transition_mask_begin < 0 &&
         op.transition_mask_count <= 0;
}

[[noreturn]] inline void exact_competitor_batch_invariant_failure(
    const char *reason, std::size_t point_count, std::size_t compiled_op_count) {
  Rcpp::stop("Exact competitor batch invariant failed: %s "
             "(points=%d compiled_ops=%d)",
             reason, static_cast<int>(point_count),
             static_cast<int>(compiled_op_count));
}

inline void exact_competitor_batch_require(bool condition, const char *reason,
                                           std::size_t point_count,
                                           std::size_t compiled_op_count) {
  if (!condition) {
    exact_competitor_batch_invariant_failure(reason, point_count,
                                             compiled_op_count);
  }
}

template <typename PointVec>
void exact_competitor_survival_batch_impl(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const PointVec &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  survival_out.assign(points.size(), 1.0);
  if (points.empty() || competitor_cache.compiled_ops.empty()) {
    return;
  }
  record_unified_outcome_exact_competitor_batch_call(
      static_cast<std::uint64_t>(points.size()),
      static_cast<std::uint64_t>(competitor_cache.compiled_ops.size()));

  PointVec mutable_points;
  PointVec *eval_points_ptr = nullptr;
  if (competitor_cache.mutates_forced_survive) {
    mutable_points = points;
    eval_points_ptr = &mutable_points;
  }
  const PointVec &eval_points =
      eval_points_ptr ? *eval_points_ptr : points;
  std::vector<double> times(points.size(), 0.0);
  for (std::size_t i = 0; i < points.size(); ++i) {
    times[i] = exact_point_ref(points, i).t;
  }

  if (exact_competitor_single_node_fastpath_eligible(competitor_cache)) {
    uuber::KernelNodeBatchValues node_values;
    const int node_idx =
        competitor_cache.compiled_ops.front().target_node_indices.front();
    exact_competitor_batch_require(
        exact_eval_node_batch_with_points(
            ctx, node_idx, points, component_idx, trial_params, trial_type_key,
            EvalNeed::kSurvival, node_values, nullptr,
            uniform_trial_params_soa),
        "single-node exact batch evaluation failed", points.size(),
        competitor_cache.compiled_ops.size());
    exact_competitor_batch_require(
        node_values.survival.size() == survival_out.size(),
        "single-node exact batch size mismatch", points.size(),
        competitor_cache.compiled_ops.size());
    for (std::size_t i = 0; i < survival_out.size(); ++i) {
      survival_out[i] = clamp_probability(node_values.survival[i]);
    }
    return;
  }

  if (exact_competitor_guard_fastpath_eligible(ctx, competitor_cache)) {
    record_unified_outcome_exact_competitor_guard_fastpath_call();
    std::vector<std::uint8_t> active(eval_points.size(), 0u);
    for (std::size_t i = 0; i < eval_points.size(); ++i) {
      active[i] = (std::isfinite(exact_point_ref(eval_points, i).t) &&
                   exact_point_ref(eval_points, i).t >= 0.0 &&
                   survival_out[i] > 0.0)
                      ? 1u
                      : 0u;
    }
    uuber::KernelNodeBatchValues guard_values;
    for (const uuber::CompetitorCompiledOp &op : competitor_cache.compiled_ops) {
      exact_competitor_batch_require(
          exact_eval_node_batch_with_points(
              ctx, op.target_node_indices.front(), eval_points, component_idx,
              trial_params, trial_type_key, EvalNeed::kSurvival, guard_values,
              &active, uniform_trial_params_soa),
          "guard competitor exact batch evaluation failed", points.size(),
          competitor_cache.compiled_ops.size());
      for (std::size_t i = 0; i < survival_out.size(); ++i) {
        if (active[i] == 0u) {
          continue;
        }
        const double surv = clamp_probability(guard_values.survival[i]);
        if (!std::isfinite(surv) || surv <= 0.0) {
          survival_out[i] = 0.0;
          active[i] = 0u;
          continue;
        }
        survival_out[i] *= surv;
        if (!std::isfinite(survival_out[i]) || survival_out[i] <= 0.0) {
          survival_out[i] = 0.0;
          active[i] = 0u;
        }
      }
      if (op.transition_mask_begin >= 0 && op.transition_mask_count > 0 &&
          eval_points_ptr != nullptr) {
        for (std::size_t i = 0; i < eval_points_ptr->size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          exact_apply_transition_mask_to_point(ctx, *eval_points_ptr, i, op);
        }
      }
    }
    return;
  }

  uuber::KernelBatchRuntimeState batch_runtime;

  uuber::KernelEventBatchEvalFn event_eval_batch =
      [&](int event_idx, const std::vector<double> &batch_times,
          const uuber::KernelEvalNeed &kneed,
          uuber::KernelNodeBatchValues &out_values) -> bool {
    return exact_eval_event_batch_common(
        ctx, eval_points, batch_times, component_idx, trial_params,
        trial_type_key, event_idx, kneed, out_values, nullptr,
        [&](const uuber::IrEvent &event, EvalNeed need,
            uuber::KernelNodeBatchValues &out) -> bool {
          return exact_eval_simple_acc_event_batch(
              ctx, event, eval_points, component_idx, trial_params,
              trial_type_key, need, out, nullptr, uniform_trial_params_soa);
        },
        uniform_trial_params_soa);
  };

  uuber::KernelGuardBatchEvalFn guard_eval_batch =
      [&](const uuber::KernelOp &op, const std::vector<double> &batch_times,
          const uuber::KernelNodeBatchValues &reference_values,
          const uuber::KernelNodeBatchValues &blocker_values,
          const uuber::KernelEvalNeed &kneed,
          uuber::KernelNodeBatchValues &out_values) -> bool {
    return exact_eval_guard_batch_common(
        ctx, eval_points, batch_times, component_idx, trial_params,
        trial_type_key, op, reference_values, blocker_values, kneed,
        out_values, nullptr, uniform_trial_params_soa);
  };

  uuber::KernelBatchTransitionApplyFn apply_transition_batch =
      [&](const uuber::CompetitorCompiledOp &op,
          uuber::KernelBatchRuntimeState &runtime) {
    if (op.transition_mask_begin < 0 || op.transition_mask_count <= 0) {
      return;
    }
    if (eval_points_ptr == nullptr) {
      return;
    }
    for (std::size_t i = 0; i < eval_points_ptr->size(); ++i) {
      exact_apply_transition_mask_to_point(ctx, *eval_points_ptr, i, op);
    }
    if (op.transition_invalidate_slot >= 0) {
      uuber::invalidate_kernel_batch_runtime_from_slot(
          runtime, op.transition_invalidate_slot);
    }
  };

  exact_competitor_batch_require(
      uuber::eval_kernel_competitor_product_batch_incremental(
          ctx.kernel_program, batch_runtime, competitor_cache.compiled_ops,
          times, event_eval_batch, guard_eval_batch, apply_transition_batch,
          survival_out),
      "generic kernel competitor batch evaluation failed", points.size(),
      competitor_cache.compiled_ops.size());
}

} // namespace

namespace uuber {

void clear_exact_outcome_runtime_caches(
    std::uint64_t runtime_cache_instance_id) noexcept {
  if (runtime_cache_instance_id == 0) {
    return;
  }
  exact_scenario_compiler_cache_registry().erase(runtime_cache_instance_id);
}

} // namespace uuber

bool exact_eval_node_batch_from_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_eval_node_batch_with_points(ctx, node_idx, points, component_idx,
                                           trial_params, trial_type_key, need,
                                           out_values, nullptr,
                                           uniform_trial_params_soa);
}

bool exact_collect_scenarios_batch_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &scenario_points,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  ExactScenarioBatch scenario_batch;
  if (!exact_collect_scenarios_batch_aligned_from_points_impl(
          compiler, ctx, node_idx, seed_points, component_idx, trial_params,
          trial_type_key, scenario_batch, uniform_trial_params_soa)) {
    scenario_points.clear();
    return false;
  }
  exact_scenario_batch_materialize_points(scenario_batch, scenario_points);
  return !scenario_points.empty();
}

bool exact_collect_scenarios_batch_aligned_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, ExactScenarioBatch &scenario_batch,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_collect_scenarios_batch_aligned_from_points_impl(
      compiler, ctx, node_idx, seed_points, component_idx, trial_params,
      trial_type_key, scenario_batch, uniform_trial_params_soa);
}

bool exact_collect_deterministic_scenarios_batch_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &aligned_points,
    std::vector<std::uint8_t> &active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_collect_deterministic_scenarios_batch_from_points_internal(
      compiler, ctx, node_idx, seed_points, component_idx, trial_params,
      trial_type_key, aligned_points, active_mask,
      uniform_trial_params_soa);
}

bool exact_collect_deterministic_scenarios_batch_aligned_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, ExactScenarioBatch &aligned_batch,
    std::vector<std::uint8_t> &active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  const RankedNodeTransitionPlan &plan = compiler.plan_for_node(node_idx);
  if (!plan.valid || plan.transitions.size() != 1u) {
    aligned_batch.clear();
    return false;
  }

  const RankedTransitionTemplate &transition = plan.transitions.front();
  if (transition.steps.empty() ||
      transition.steps.front().kind != RankedTransitionStepKind::EvalDensityNode) {
    aligned_batch.clear();
    return false;
  }
  for (std::size_t step_idx = 1; step_idx < transition.steps.size(); ++step_idx) {
    if (transition.steps[step_idx].kind !=
            RankedTransitionStepKind::AddCompleteSources &&
        transition.steps[step_idx].kind !=
            RankedTransitionStepKind::AddSurviveSources) {
      aligned_batch.clear();
      return false;
    }
  }

  uuber::KernelNodeBatchValues step_values;
  if (!exact_eval_node_batch_with_points(
          ctx, transition.steps.front().node_idx, seed_points, component_idx,
          trial_params, trial_type_key, EvalNeed::kDensity, step_values, nullptr,
          uniform_trial_params_soa) ||
      step_values.density.size() != seed_points.size()) {
    aligned_batch.clear();
    return false;
  }

  aligned_batch.clear();
  aligned_batch.reserve(seed_points.size());
  active_mask.assign(seed_points.size(), 0u);
  for (const ExactScenarioPoint &point : seed_points) {
    exact_scenario_batch_append_point(aligned_batch, point, point.weight);
  }
  for (std::size_t i = 0; i < seed_points.size(); ++i) {
    const ExactScenarioPoint &seed_point = seed_points[i];
    if (!(std::isfinite(seed_point.t) && seed_point.t >= 0.0 &&
          std::isfinite(seed_point.weight) && seed_point.weight > 0.0)) {
      continue;
    }
    const double factor = safe_density(step_values.density[i]);
    if (!std::isfinite(factor) || factor <= 0.0) {
      aligned_batch.weight[i] = 0.0;
      continue;
    }
    aligned_batch.weight[i] *= factor;
    if (!std::isfinite(aligned_batch.weight[i]) ||
        aligned_batch.weight[i] <= 0.0) {
      aligned_batch.weight[i] = 0.0;
      continue;
    }
    bool valid = true;
    for (std::size_t step_idx = 1; step_idx < transition.steps.size();
         ++step_idx) {
      const RankedTransitionStep &step = transition.steps[step_idx];
      if (step.kind == RankedTransitionStepKind::AddCompleteSources) {
        for (int id : step.source_ids) {
          if (!time_constraints_mark_complete(
                  id, aligned_batch.t[i], step.bind_exact_current_time,
                  aligned_batch.time_constraints[i])) {
            valid = false;
            break;
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddSurviveSources) {
        for (int id : step.source_ids) {
          if (!time_constraints_mark_survive(id, aligned_batch.t[i],
                                             aligned_batch.time_constraints[i])) {
            valid = false;
            break;
          }
        }
      }
      if (!valid) {
        break;
      }
    }
    if (valid) {
      active_mask[i] = 1u;
    } else {
      aligned_batch.weight[i] = 0.0;
    }
  }
  return true;
}

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const std::vector<ExactScenarioPoint> &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  exact_competitor_survival_batch_impl(
      ctx, competitor_cache, component_idx, trial_params, trial_type_key,
      points, survival_out, uniform_trial_params_soa);
}

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const ExactScenarioBatch &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  exact_competitor_survival_batch_impl(
      ctx, competitor_cache, component_idx, trial_params, trial_type_key,
      points, survival_out, uniform_trial_params_soa);
}

bool exact_outcome_density_batch_from_state(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    const std::vector<double> &times, std::vector<double> &density_out,
    const uuber::CompetitorClusterCacheEntry *competitor_cache) {
  record_unified_outcome_exact_batch_density_call();
  density_out.assign(times.size(), 0.0);
  if (times.empty()) {
    return true;
  }
  RankedTransitionCompiler &compiler = exact_scenario_compiler_for_ctx(ctx);
  const double saved_t = state.t;
  ExactScenarioBatch scenario_points;
  if (!exact_collect_scenarios_batch_impl(compiler, ctx, node_idx, times, state,
                                          scenario_points)) {
    state.t = saved_t;
    invalidate_kernel_runtime_root(state);
    return true;
  }

  if (competitor_ids.empty()) {
    for (std::size_t i = 0; i < scenario_points.size(); ++i) {
      const double weight = scenario_points.weight[i];
      const std::size_t density_index = scenario_points.density_index[i];
      if (!std::isfinite(weight) || weight <= 0.0 ||
          density_index >= density_out.size()) {
        continue;
      }
      density_out[density_index] += weight;
    }
    for (double &value : density_out) {
      value = safe_density(value);
    }
  } else if (!scenario_points.empty()) {
    const uuber::CompetitorClusterCacheEntry &cache =
        competitor_cache ? *competitor_cache
                         : fetch_competitor_cluster_cache(ctx, competitor_ids);
    const uuber::TrialParamsSoA *uniform_trial_params_soa =
        state.trial_params_soa_batch == nullptr ? state.trial_params_soa
                                                : nullptr;
    std::vector<double> survival_out;
    exact_competitor_survival_batch(
        ctx, cache, state.component_idx, state.trial_params,
        state.trial_type_key, scenario_points, survival_out,
        uniform_trial_params_soa);
    for (std::size_t i = 0; i < scenario_points.size(); ++i) {
      const double surv = clamp_probability(survival_out[i]);
      if (!std::isfinite(surv) || surv <= 0.0) {
        continue;
      }
      const double weight = scenario_points.weight[i];
      const std::size_t density_index = scenario_points.density_index[i];
      if (!std::isfinite(weight) || weight <= 0.0 ||
          density_index >= density_out.size()) {
        continue;
      }
      density_out[density_index] += weight * surv;
    }
    for (double &value : density_out) {
      value = safe_density(value);
    }
  }
  state.t = saved_t;
  invalidate_kernel_runtime_root(state);
  return true;
}
