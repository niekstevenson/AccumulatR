#include "exact_outcome_density.h"

#include <Rcpp.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "competitor_cache.h"
#include "integrate.h"
#include "ranked_transitions.h"
#include "vector_lane_utils.h"

namespace {

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

inline bool exact_node_allowed_in_component(
    const uuber::NativeContext &ctx, int node_idx, int component_idx) {
  if (component_idx < 0 || node_idx < 0 ||
      node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return true;
  }
  const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  if (node.component_mask_offset >= 0) {
    return ir_mask_has_component(ctx, node.component_mask_offset,
                                 component_idx);
  }
  auto out_it = ctx.ir.node_idx_to_outcomes.find(node_idx);
  if (out_it == ctx.ir.node_idx_to_outcomes.end() || out_it->second.empty()) {
    return true;
  }
  for (int outcome_idx : out_it->second) {
    if (ir_outcome_allows_component(ctx, outcome_idx, component_idx)) {
      return true;
    }
  }
  return false;
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

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask = nullptr,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const uuber::ExactScenarioLaneViewBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask = nullptr,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

template <typename PointVec>
bool exact_eval_simple_acc_event_batch(
    const uuber::NativeContext &ctx, const uuber::IrEvent &event,
    const PointVec &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask = nullptr,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

inline void exact_build_seed_batch_from_state(
    const std::vector<double> &times, const NodeEvalState &state,
    ExactScenarioBatch &seed_batch,
    const uuber::TrialParamsSoA *&uniform_trial_params_soa);

template <typename ProcessBatchFn>
ExactTransitionExecutionResult exact_execute_compiled_transition_plan_process_impl(
    const RankedNodeBatchPlan &plan, const uuber::NativeContext &ctx,
    const ExactScenarioBatch &seed_batch, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    ProcessBatchFn &&process_batch,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

template <typename PointVec>
void exact_competitor_survival_batch_impl(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const PointVec &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

int exact_resolve_transition_node_idx_impl(
    const uuber::NativeContext &ctx, int node_idx, int component_idx) {
  int current_idx = node_idx;
  while (current_idx >= 0 &&
         current_idx < static_cast<int>(ctx.ir.nodes.size())) {
    const uuber::IrNode &node =
        ctx.ir.nodes[static_cast<std::size_t>(current_idx)];
    if (node.op != uuber::IrNodeOp::Guard) {
      break;
    }
    if (!exact_node_allowed_in_component(ctx, node.reference_idx,
                                         component_idx)) {
      return NA_INTEGER;
    }
    if (!exact_node_allowed_in_component(ctx, node.blocker_idx,
                                         component_idx)) {
      current_idx = node.reference_idx;
      continue;
    }
    break;
  }
  return current_idx;
}

inline void exact_scenario_batch_append_lane_copy(
    ExactScenarioBatch &batch, const uuber::VectorLaneRef &point, double t,
    double weight) {
  batch.t.push_back(t);
  batch.density_index.push_back(point.density_index);
  batch.weight.push_back(weight);
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

inline void exact_scenario_view_batch_append_lane(
    uuber::ExactScenarioLaneViewBatch &batch,
    const uuber::VectorLaneRef &point, double weight,
    const uuber::TrialParamsSoA *resolved_trial_params_soa = nullptr) {
  batch.t.push_back(point.t);
  batch.density_index.push_back(point.density_index);
  batch.weight.push_back(weight);
  batch.trial_params_soa.push_back(resolved_trial_params_soa != nullptr
                                       ? resolved_trial_params_soa
                                       : point.trial_params_soa);
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

inline void exact_build_scenario_view_batch_from_seed(
    const ExactScenarioBatch &seed_batch, const std::vector<double> &weights,
    uuber::ExactScenarioLaneViewBatch &batch) {
  batch.clear();
  batch.reserve(seed_batch.size());
  for (std::size_t i = 0; i < seed_batch.size(); ++i) {
    exact_scenario_view_batch_append_lane(
        batch, uuber::vector_lane_ref(seed_batch, i),
        i < weights.size() ? weights[i] : 0.0);
  }
}

inline void exact_apply_transition_mask_to_point(
    const uuber::NativeContext &ctx, ExactScenarioBatch &points,
    std::size_t idx, const uuber::TreeCompetitorOp &op) {
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
                              forced_survive_bits_valid);
  if (idx < points.forced_survive_bits_valid.size()) {
    points.forced_survive_bits_valid[idx] =
        forced_survive_bits_valid ? 1u : 0u;
  }
}

inline void exact_apply_transition_mask_to_point(
    const uuber::NativeContext &ctx, uuber::ExactScenarioLaneViewBatch &points,
    std::size_t idx, const uuber::TreeCompetitorOp &op) {
  if (idx >= points.size() || idx >= points.mutable_forced_survive_bits.size()) {
    Rcpp::stop(
        "Exact competitor lane-view invariant failed: missing mutable "
        "forced-survive storage");
  }
  bool forced_survive_bits_valid =
      idx < points.forced_survive_bits_valid.size() &&
      points.forced_survive_bits_valid[idx] != 0u;
  apply_transition_mask_words(
      ctx, op.transition_mask_begin, op.transition_mask_count,
      op.transition_invalidate_slot, points.mutable_forced_survive_bits[idx],
      forced_survive_bits_valid);
  if (idx < points.forced_survive_bits.size()) {
    points.forced_survive_bits[idx] = &points.mutable_forced_survive_bits[idx];
  }
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
  const uuber::VectorLaneRef point = uuber::vector_lane_ref(points, idx);
  exact_scenario_batch_append_lane_copy(batch, point, t_new, point.weight);
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
    uuber::TreeNodeBatchValues &out_values) {
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
    if (!uuber::vector_lane_active(active_mask, i) ||
        (need_density_mask[i] == 0u && need_cdf_mask[i] == 0u)) {
      continue;
    }
    ExactRelativeOnsetBatchPlan &plan = plans[i];
    plan.need_density = need_density_mask[i] != 0u;
    plan.need_cdf = need_cdf_mask[i] != 0u;

    const uuber::VectorLaneRef point = uuber::vector_lane_ref(points, i);
    const double eval_t = eval_times[i];
    const ForcedStateView forced_state =
        uuber::make_vector_lane_forced_state_view(ctx, point);
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
    uuber::TreeNodeBatchValues cdf_values;
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
    uuber::TreeNodeBatchValues source_density_values;
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
    if (!uuber::vector_lane_active(active_mask, i) ||
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
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  out_values = uuber::TreeNodeBatchValues{};
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
    if (uuber::vector_lane_active(active_mask, i)) {
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
          : uuber::resolve_vector_lane_trial_params_soa(
                ctx, uuber::vector_lane_ref(points, anchor_idx), trial_params);
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
    const uuber::VectorLaneRef point = uuber::vector_lane_ref(points, i);
    if (uniform_trial_params_soa == nullptr) {
      const uuber::TrialParamsSoA *point_trial_params_soa =
          uuber::resolve_vector_lane_trial_params_soa(ctx, point, trial_params,
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
          (!(points.t[i] > self_bound_lower[i]) ||
           points.t[i] > self_bound_upper[i])) {
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
                                               cfg, lower - (points.t[i] - shifted_times[i]),
                                               lower_bound))
                     : (i < relative_lower_cdf.size()
                            ? clamp_probability(relative_lower_cdf[i])
                            : 0.0))
              : 0.0;
      const double upper_cdf =
          std::isfinite(upper)
              ? (direct_shifted_eval[i] != 0u
                     ? clamp_probability(success_probs[i] * eval_cdf_single_with_lower_bound(
                                               cfg, upper - (points.t[i] - shifted_times[i]),
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
        if (std::isfinite(points.t[i]) &&
            points.t[i] > lower &&
            points.t[i] <= upper) {
          out_values.density[i] =
              safe_density(base_density / condition_mass);
        } else {
          out_values.density[i] = 0.0;
        }
      }
      if (needs_cdf(need) || needs_survival(need)) {
        double cdf_cond = 0.0;
        if (!std::isfinite(points.t[i]) ||
            points.t[i] >= upper) {
          cdf_cond = 1.0;
        } else if (!(points.t[i] > lower)) {
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

thread_local int g_exact_node_batch_depth = 0;

struct ExactNodeBatchDepthGuard {
  ExactNodeBatchDepthGuard() { ++g_exact_node_batch_depth; }
  ~ExactNodeBatchDepthGuard() { --g_exact_node_batch_depth; }

  bool recursive_call() const noexcept { return g_exact_node_batch_depth > 1; }
};

inline void exact_init_batch_values(std::size_t point_count,
                                    uuber::TreeNodeBatchValues &out_values) {
  out_values.density.assign(point_count, 0.0);
  out_values.survival.assign(point_count, 1.0);
  out_values.cdf.assign(point_count, 0.0);
}

inline bool exact_batch_values_size_matches(
    const uuber::TreeNodeBatchValues &values, std::size_t point_count) {
  return values.density.size() == point_count &&
         values.survival.size() == point_count &&
         values.cdf.size() == point_count;
}

template <typename PointVec>
inline void exact_collect_active_subset(
    const PointVec &points, const std::vector<std::uint8_t> *active_mask,
    std::vector<std::size_t> &subset_indices,
    std::vector<double> &subset_times) {
  subset_indices.clear();
  subset_times.clear();
  subset_indices.reserve(points.size());
  subset_times.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (!uuber::vector_lane_active(active_mask, i)) {
      continue;
    }
    subset_indices.push_back(i);
    subset_times.push_back(points.t[i]);
  }
}

inline bool exact_merge_batch_values_indices(
    const uuber::TreeNodeBatchValues &src_values,
    const std::vector<std::size_t> &dst_indices,
    uuber::TreeNodeBatchValues &dst_values) {
  const std::size_t point_count = dst_indices.size();
  if (!exact_batch_values_size_matches(src_values, point_count) ||
      !exact_batch_values_size_matches(dst_values, dst_values.density.size())) {
    return false;
  }
  for (std::size_t i = 0; i < point_count; ++i) {
    const std::size_t dst_idx = dst_indices[i];
    if (dst_idx >= dst_values.density.size()) {
      return false;
    }
    dst_values.density[dst_idx] = src_values.density[i];
    dst_values.survival[dst_idx] = src_values.survival[i];
    dst_values.cdf[dst_idx] = src_values.cdf[i];
  }
  return true;
}

template <typename PointVec, typename EvalFn>
inline bool exact_with_shared_state_compact_subset(
    NodeEvalState &shared_state, const PointVec &points,
    const std::vector<std::uint8_t> *shared_group_mask,
    std::vector<std::size_t> &subset_indices,
    std::vector<double> &subset_times,
    std::vector<const uuber::TrialParamsSoA *> &subset_trial_params_soa,
    EvalFn &&eval_fn) {
  exact_collect_active_subset(points, shared_group_mask, subset_indices,
                              subset_times);
  if (subset_indices.empty()) {
    return true;
  }

  const auto *saved_trial_params_soa_batch = shared_state.trial_params_soa_batch;
  const std::vector<double> *group_times = &subset_times;
  if (subset_indices.size() == points.size()) {
    group_times = &points.t;
  } else if (saved_trial_params_soa_batch != nullptr) {
    subset_trial_params_soa.resize(subset_indices.size());
    for (std::size_t i = 0; i < subset_indices.size(); ++i) {
      const std::size_t src_idx = subset_indices[i];
      subset_trial_params_soa[i] =
          src_idx < saved_trial_params_soa_batch->size()
              ? (*saved_trial_params_soa_batch)[src_idx]
              : shared_state.trial_params_soa;
    }
    shared_state.trial_params_soa_batch = &subset_trial_params_soa;
  }

  const bool ok = eval_fn(*group_times);
  shared_state.trial_params_soa_batch = saved_trial_params_soa_batch;
  return ok;
}

template <typename PointVec, typename EvalFn>
inline bool exact_eval_shared_state_compact_subset(
    NodeEvalState &shared_state, const PointVec &points,
    const std::vector<std::uint8_t> *shared_group_mask, EvalFn &&eval_fn,
    uuber::TreeNodeBatchValues &shared_out,
    std::vector<std::size_t> &subset_indices,
    std::vector<double> &subset_times,
    std::vector<const uuber::TrialParamsSoA *> &subset_trial_params_soa) {
  const bool ok = exact_with_shared_state_compact_subset(
      shared_state, points, shared_group_mask, subset_indices, subset_times,
      subset_trial_params_soa,
      [&](const std::vector<double> &group_times) -> bool {
        if (group_times.empty()) {
          shared_out = uuber::TreeNodeBatchValues{};
          return true;
        }
        return eval_fn(group_times, shared_state, shared_out);
      });
  return ok && (subset_indices.empty() ||
                exact_batch_values_size_matches(shared_out,
                                                subset_indices.size()));
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa);

template <typename PointVec>
bool exact_eval_node_batch_common(
    const uuber::NativeContext &ctx, int node_idx, const PointVec &points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, EvalNeed need,
    uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa);

template <typename PointVec>
bool exact_eval_node_batch_common(
    const uuber::NativeContext &ctx, int node_idx,
    const PointVec &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  if (points.empty()) {
    out_values = uuber::TreeNodeBatchValues{};
    return true;
  }
  ExactNodeBatchDepthGuard depth_guard;
  const std::size_t active_point_count =
      uuber::count_active_vector_lanes(active_mask, points.size());
  (void)depth_guard;

  std::vector<std::size_t> subset_indices;
  std::vector<double> subset_times;
  std::vector<const uuber::TrialParamsSoA *> subset_trial_params_soa;
  auto eval_partitioned_node_batch =
      [&](const std::vector<int> &relevant_source_ids) -> bool {
    exact_init_batch_values(points.size(), out_values);
    return eval_vector_lanes_with_shared_node_state(
        ctx, points, active_mask, relevant_source_ids, component_idx,
        trial_params, trial_type_key,
        [&](NodeEvalState &shared_state,
            const std::vector<std::uint8_t> *shared_group_mask) -> bool {
          uuber::TreeNodeBatchValues shared_out;
          if (!exact_eval_shared_state_compact_subset(
                  shared_state, points, shared_group_mask,
                  [&](const std::vector<double> &group_times,
                      NodeEvalState &group_state,
                      uuber::TreeNodeBatchValues &group_out) -> bool {
                    return evaluator_eval_node_batch_with_state_dense(
                        node_idx, group_times, group_state, need, group_out);
                  },
                  shared_out, subset_indices, subset_times,
                  subset_trial_params_soa)) {
            return false;
          }
          return exact_merge_batch_values_indices(shared_out, subset_indices,
                                                 out_values);
        },
        [&]() -> bool {
          const uuber::IrNode &node = ir_node_required(ctx, node_idx);
          Rcpp::stop(
              "Exact batch invariant failed for node: "
              "node_idx=%d op=%d event_idx=%d active_points=%d",
              node_idx, static_cast<int>(node.op), node.event_idx,
              static_cast<int>(active_point_count));
          return false;
        },
        [&](std::uint64_t active_count, std::uint64_t group_count) {
          (void)active_count;
          (void)group_count;
        },
        false, uniform_trial_params_soa);
  };

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
      return eval_partitioned_node_batch(
          relevant_source_ids_for_batch_node(ctx, node_idx));
    }
  }

  return eval_partitioned_node_batch(
      relevant_source_ids_for_batch_node(ctx, node_idx));
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_eval_node_batch_common(ctx, node_idx, points, component_idx,
                                      trial_params, trial_type_key, need,
                                      out_values, active_mask,
                                      uniform_trial_params_soa);
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const uuber::ExactScenarioLaneViewBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_eval_node_batch_common(ctx, node_idx, points, component_idx,
                                      trial_params, trial_type_key, need,
                                      out_values, active_mask,
                                      uniform_trial_params_soa);
}

inline EvalNeed exact_ranked_eval_need(RankedProgramEvalKind kind) {
  switch (kind) {
  case RankedProgramEvalKind::Density:
    return EvalNeed::kDensity;
  case RankedProgramEvalKind::CDF:
    return EvalNeed::kCDF;
  case RankedProgramEvalKind::Survival:
    return EvalNeed::kSurvival;
  }
  return EvalNeed::kDensity;
}

inline double
exact_ranked_eval_factor(RankedProgramEvalKind kind,
                         const uuber::TreeNodeBatchValues &values,
                         std::size_t idx) {
  switch (kind) {
  case RankedProgramEvalKind::Density:
    return safe_density(values.density[idx]);
  case RankedProgramEvalKind::CDF:
    return clamp_probability(values.cdf[idx]);
  case RankedProgramEvalKind::Survival:
    return clamp_probability(values.survival[idx]);
  }
  return 0.0;
}

inline uuber::BitsetState *exact_mutable_forced_survive_bits_for_lane(
    uuber::ExactScenarioLaneViewBatch &batch, std::size_t idx) {
  if (idx >= batch.size()) {
    return nullptr;
  }
  if (batch.mutable_forced_survive_bits.size() < batch.size()) {
    batch.mutable_forced_survive_bits.resize(batch.size());
  }
  uuber::BitsetState *mutable_bits = &batch.mutable_forced_survive_bits[idx];
  if (idx < batch.forced_survive_bits.size() &&
      batch.forced_survive_bits[idx] != mutable_bits) {
    if (batch.forced_survive_bits[idx] != nullptr) {
      *mutable_bits = *batch.forced_survive_bits[idx];
    } else {
      mutable_bits->clear_all();
    }
    batch.forced_survive_bits[idx] = mutable_bits;
  }
  return mutable_bits;
}

inline TimeConstraintMap *
exact_mutable_time_constraints_for_lane(uuber::ExactScenarioLaneViewBatch &batch,
                                        std::size_t idx) {
  if (idx >= batch.size()) {
    return nullptr;
  }
  if (batch.mutable_time_constraints.size() < batch.size()) {
    batch.mutable_time_constraints.resize(batch.size());
  }
  TimeConstraintMap *mutable_constraints =
      &batch.mutable_time_constraints[idx];
  if (idx < batch.time_constraints.size() &&
      batch.time_constraints[idx] != mutable_constraints) {
    if (batch.time_constraints[idx] != nullptr) {
      *mutable_constraints = *batch.time_constraints[idx];
    } else {
      mutable_constraints->clear();
    }
    batch.time_constraints[idx] = mutable_constraints;
  }
  return mutable_constraints;
}

inline bool exact_apply_ranked_state_delta_row(
    uuber::ExactScenarioLaneViewBatch &batch, std::size_t idx,
    const RankedStateDelta &delta) {
  uuber::BitsetState *forced_survive_bits =
      exact_mutable_forced_survive_bits_for_lane(batch, idx);
  TimeConstraintMap *time_constraints =
      exact_mutable_time_constraints_for_lane(batch, idx);
  if (forced_survive_bits == nullptr || time_constraints == nullptr) {
    return false;
  }
  bool forced_survive_bits_valid =
      idx < batch.forced_survive_bits_valid.size() &&
      batch.forced_survive_bits_valid[idx] != 0u;
  const bool applied = apply_ranked_state_delta_raw(
      batch.t[idx], delta, *batch.forced_complete_bits[idx],
      idx < batch.forced_complete_bits_valid.size() &&
          batch.forced_complete_bits_valid[idx] != 0u,
      *forced_survive_bits, forced_survive_bits_valid,
      *time_constraints);
  if (idx < batch.forced_survive_bits_valid.size()) {
    batch.forced_survive_bits_valid[idx] =
        forced_survive_bits_valid ? 1u : 0u;
  }
  return applied;
}

struct ExactScenarioCollapseKey {
  std::size_t density_index{0u};
  std::uint64_t time_bits{0u};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  SequenceStateKey state_key{};

  bool operator==(const ExactScenarioCollapseKey &other) const noexcept {
    return density_index == other.density_index &&
           time_bits == other.time_bits &&
           trial_params_soa == other.trial_params_soa &&
           state_key == other.state_key;
  }
};

struct ExactScenarioCollapseKeyHash {
  std::size_t operator()(const ExactScenarioCollapseKey &key) const noexcept {
    std::size_t seed = std::hash<std::size_t>{}(key.density_index);
    auto combine = [](std::size_t &seed_ref, std::size_t value) {
      seed_ref ^= value + 0x9e3779b97f4a7c15ULL + (seed_ref << 6) +
                  (seed_ref >> 2);
    };
    combine(seed, static_cast<std::size_t>(key.time_bits));
    combine(seed, std::hash<const uuber::TrialParamsSoA *>{}(key.trial_params_soa));
    combine(seed, SequenceStateKeyHash{}(key.state_key));
    return seed;
  }
};

inline ExactScenarioCollapseKey
exact_make_scenario_collapse_key(const uuber::VectorLaneRef &point) {
  ExactScenarioCollapseKey key;
  key.density_index = point.density_index;
  key.time_bits = canonical_double_bits(point.t);
  key.trial_params_soa = point.trial_params_soa;
  key.state_key = sequence_state_key(point.forced_complete_bits,
                                     point.forced_complete_bits_valid,
                                     point.forced_survive_bits,
                                     point.forced_survive_bits_valid,
                                     point.time_constraints);
  return key;
}

inline void exact_collapse_active_scenario_lane_view(
    uuber::ExactScenarioLaneViewBatch &points, std::vector<std::uint8_t> &active,
    std::vector<double> &weights) {
  if (points.size() < 2u || active.size() != points.size() ||
      weights.size() != points.size()) {
    return;
  }
  std::unordered_map<ExactScenarioCollapseKey, std::size_t,
                     ExactScenarioCollapseKeyHash>
      first_index_by_key;
  first_index_by_key.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (active[i] == 0u || !std::isfinite(weights[i]) || weights[i] <= 0.0) {
      continue;
    }
    const ExactScenarioCollapseKey key =
        exact_make_scenario_collapse_key(uuber::vector_lane_ref(points, i));
    const auto it = first_index_by_key.find(key);
    if (it == first_index_by_key.end()) {
      first_index_by_key.emplace(key, i);
      continue;
    }
    const std::size_t dst_idx = it->second;
    weights[dst_idx] += weights[i];
    if (dst_idx < points.weight.size()) {
      points.weight[dst_idx] = weights[dst_idx];
    }
    weights[i] = 0.0;
    if (i < points.weight.size()) {
      points.weight[i] = 0.0;
    }
    active[i] = 0u;
  }
}

inline void exact_build_seed_batch_from_state(
    const std::vector<double> &times, const NodeEvalState &state,
    ExactScenarioBatch &seed_batch,
    const uuber::TrialParamsSoA *&uniform_trial_params_soa) {
  seed_batch.clear();
  seed_batch.reserve(times.size());
  for (std::size_t i = 0; i < times.size(); ++i) {
    seed_batch.t.push_back(times[i]);
    seed_batch.density_index.push_back(i);
    seed_batch.weight.push_back((std::isfinite(times[i]) && times[i] >= 0.0)
                                    ? 1.0
                                    : 0.0);
    if (state.trial_params_soa_batch != nullptr &&
        i < state.trial_params_soa_batch->size()) {
      seed_batch.trial_params_soa.push_back((*state.trial_params_soa_batch)[i]);
    } else {
      seed_batch.trial_params_soa.push_back(state.trial_params_soa);
    }
    seed_batch.forced_complete_bits.push_back(state.forced_complete_bits);
    seed_batch.forced_survive_bits.push_back(state.forced_survive_bits);
    seed_batch.forced_complete_bits_valid.push_back(
        state.forced_complete_bits_valid ? 1u : 0u);
    seed_batch.forced_survive_bits_valid.push_back(
        state.forced_survive_bits_valid ? 1u : 0u);
    seed_batch.time_constraints.push_back(state.time_constraints);
  }
  uniform_trial_params_soa =
      state.trial_params_soa_batch == nullptr ? state.trial_params_soa : nullptr;
}

inline std::size_t exact_initialize_seed_transition_state(
    const ExactScenarioBatch &seed_batch, std::vector<std::uint8_t> &active,
    std::vector<double> &transition_weights) {
  active.assign(seed_batch.size(), 0u);
  transition_weights.assign(seed_batch.size(), 0.0);
  std::size_t active_count = 0u;
  for (std::size_t i = 0; i < seed_batch.size(); ++i) {
    const uuber::VectorLaneRef point = uuber::vector_lane_ref(seed_batch, i);
    const bool is_active = std::isfinite(point.t) && point.t >= 0.0 &&
                           std::isfinite(point.weight) && point.weight > 0.0;
    active[i] = is_active ? 1u : 0u;
    transition_weights[i] = is_active ? point.weight : 0.0;
    active_count += is_active ? 1u : 0u;
  }
  return active_count;
}

template <typename PointVec>
bool exact_eval_ranked_batch_plan_slice_fused(
    const uuber::NativeContext &ctx, const RankedProgramSliceRef &slice,
    const PointVec &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<std::uint8_t> &active,
    const uuber::TrialParamsSoA *uniform_trial_params_soa,
    std::vector<uuber::TreeNodeBatchValues> &step_values_out) {
  step_values_out.clear();
  if (slice.evals.empty()) {
    return true;
  }

  if (slice.evals.size() == 1u) {
    step_values_out.resize(1u);
    return exact_eval_node_batch_with_points(
               ctx, slice.evals.front().node_idx, points, component_idx,
               trial_params, trial_type_key,
               exact_ranked_eval_need(slice.evals.front().kind),
               step_values_out.front(), &active, uniform_trial_params_soa) &&
           exact_batch_values_size_matches(step_values_out.front(), points.size());
  }

  step_values_out.resize(slice.evals.size());
  for (uuber::TreeNodeBatchValues &values : step_values_out) {
    exact_init_batch_values(points.size(), values);
  }

  std::vector<std::size_t> subset_indices;
  std::vector<double> subset_times;
  std::vector<const uuber::TrialParamsSoA *> subset_trial_params_soa;
  std::vector<uuber::TreeNodeBatchValues> shared_values;
  shared_values.reserve(slice.evals.size());

  return eval_vector_lanes_with_shared_node_state(
      ctx, points, &active, slice.relevant_source_ids, component_idx,
      trial_params, trial_type_key,
      [&](NodeEvalState &shared_state,
          const std::vector<std::uint8_t> *shared_group_mask) -> bool {
        return exact_with_shared_state_compact_subset(
            shared_state, points, shared_group_mask, subset_indices,
            subset_times, subset_trial_params_soa,
            [&](const std::vector<double> &group_times) -> bool {
              if (group_times.empty()) {
                return true;
              }
              if (!evaluator_eval_nodes_batch_with_state_dense(
                      slice.node_indices, group_times, shared_state,
                      shared_values) ||
                  shared_values.size() != slice.evals.size()) {
                return false;
              }
              for (std::size_t eval_idx = 0; eval_idx < shared_values.size();
                   ++eval_idx) {
                if (!exact_merge_batch_values_indices(
                        shared_values[eval_idx], subset_indices,
                        step_values_out[eval_idx])) {
                  return false;
                }
              }
              return true;
            });
      },
      [&]() -> bool {
        step_values_out.clear();
        step_values_out.reserve(slice.evals.size());
        for (const RankedProgramEvalRef &eval_ref : slice.evals) {
          uuber::TreeNodeBatchValues step_values;
          if (!exact_eval_node_batch_with_points(
                  ctx, eval_ref.node_idx, points, component_idx, trial_params,
                  trial_type_key, exact_ranked_eval_need(eval_ref.kind),
                  step_values, &active, uniform_trial_params_soa) ||
              !exact_batch_values_size_matches(step_values, points.size())) {
            return false;
          }
          step_values_out.push_back(std::move(step_values));
        }
        return true;
      },
      [&](std::uint64_t active_count, std::uint64_t group_count) {
        (void)active_count;
        (void)group_count;
      },
      false, uniform_trial_params_soa);
}

template <typename ProcessBatchFn>
ExactTransitionExecutionResult exact_execute_compiled_transition_plan_process_impl(
    const RankedNodeBatchPlan &plan, const uuber::NativeContext &ctx,
    const ExactScenarioBatch &seed_batch, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    ProcessBatchFn &&process_batch,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  if (!plan.valid || plan.batch_plans.empty()) {
    return ExactTransitionExecutionResult::NoEmission;
  }

  uuber::ExactScenarioLaneViewBatch transition_batch;
  transition_batch.reserve(seed_batch.size());
  std::vector<std::uint8_t> seed_active;
  seed_active.reserve(seed_batch.size());
  std::vector<double> seed_transition_weights;
  seed_transition_weights.reserve(seed_batch.size());
  const std::size_t seed_active_count = exact_initialize_seed_transition_state(
      seed_batch, seed_active, seed_transition_weights);
  if (seed_active_count == 0u) {
    return ExactTransitionExecutionResult::NoEmission;
  }
  uuber::ExactScenarioLaneViewBatch base_transition_batch;
  std::vector<std::uint8_t> active;
  std::vector<std::uint8_t> base_active;
  std::vector<double> transition_weights;
  std::vector<double> base_transition_weights;
  std::vector<uuber::TreeNodeBatchValues> step_values_by_eval;
  bool emitted_any = false;
  const bool use_grouped_plans =
      plan.has_shared_slice_groups && !plan.batch_plan_groups.empty();
  const std::size_t group_count =
      use_grouped_plans ? plan.batch_plan_groups.size() : plan.batch_plans.size();
  for (std::size_t group_idx = 0; group_idx < group_count; ++group_idx) {
    const RankedBatchPlanGroup *batch_plan_group =
        use_grouped_plans ? &plan.batch_plan_groups[group_idx] : nullptr;
    const std::size_t leader_index =
        batch_plan_group != nullptr ? batch_plan_group->leader_index : group_idx;
    const RankedBatchPlan &leader_plan = plan.batch_plans[leader_index];
    active = seed_active;
    transition_weights = seed_transition_weights;
    std::size_t active_count = seed_active_count;
    const bool eval_ok = exact_eval_ranked_batch_plan_slice_fused(
        ctx, leader_plan.slice, seed_batch, component_idx, trial_params,
        trial_type_key, active, uniform_trial_params_soa, step_values_by_eval);
    if (!eval_ok || step_values_by_eval.size() != leader_plan.slice.evals.size()) {
      continue;
    }
    bool base_valid = true;
    auto deactivate_base_point = [&](std::size_t idx) {
      if (idx < active.size() && active[idx] != 0u) {
        active[idx] = 0u;
        if (active_count > 0u) {
          --active_count;
        }
      }
      if (idx < transition_weights.size()) {
        transition_weights[idx] = 0.0;
      }
    };
    for (std::size_t eval_idx = 0; eval_idx < leader_plan.slice.evals.size();
         ++eval_idx) {
      const RankedProgramEvalRef &eval_ref = leader_plan.slice.evals[eval_idx];
      const uuber::TreeNodeBatchValues &step_values =
          step_values_by_eval[eval_idx];
      if (!exact_batch_values_size_matches(step_values, seed_batch.size())) {
        base_valid = false;
        break;
      }
      for (std::size_t i = 0; i < seed_batch.size(); ++i) {
        if (active[i] == 0u) {
          continue;
        }
        const double factor =
            exact_ranked_eval_factor(eval_ref.kind, step_values, i);
        if (!std::isfinite(factor) || factor <= 0.0) {
          deactivate_base_point(i);
          continue;
        }
        transition_weights[i] *= factor;
        if (!std::isfinite(transition_weights[i]) ||
            transition_weights[i] <= 0.0) {
          deactivate_base_point(i);
        }
      }
    }
    if (!base_valid || active_count == 0u) {
      continue;
    }

    base_active = active;
    base_transition_weights = transition_weights;
    const std::size_t base_active_count = active_count;
    const bool group_has_deltas =
        batch_plan_group != nullptr ? batch_plan_group->has_deltas
                                    : !leader_plan.deltas.empty();
    if (group_has_deltas) {
      exact_build_scenario_view_batch_from_seed(
          seed_batch, base_transition_weights, base_transition_batch);
    }

    const std::size_t group_plan_count =
        batch_plan_group != nullptr ? batch_plan_group->plan_indices.size() : 1u;
    for (std::size_t plan_pos = 0; plan_pos < group_plan_count; ++plan_pos) {
      const std::size_t batch_plan_index =
          batch_plan_group != nullptr ? batch_plan_group->plan_indices[plan_pos]
                                      : group_idx;
      const RankedBatchPlan &batch_plan = plan.batch_plans[batch_plan_index];
      if (batch_plan.deltas.empty()) {
        if (!std::forward<ProcessBatchFn>(process_batch)(
                seed_batch, &base_active, base_transition_weights)) {
          return ExactTransitionExecutionResult::Error;
        }
        emitted_any = true;
        continue;
      }

      transition_batch = base_transition_batch;
      active = base_active;
      transition_weights = base_transition_weights;
      active_count = base_active_count;
      bool transition_valid = true;
      auto deactivate_transition_point = [&](std::size_t idx) {
        if (idx < active.size() && active[idx] != 0u) {
          active[idx] = 0u;
          if (active_count > 0u) {
            --active_count;
          }
        }
        if (idx < transition_weights.size()) {
          transition_weights[idx] = 0.0;
        }
        if (idx < transition_batch.weight.size()) {
          transition_batch.weight[idx] = 0.0;
        }
      };
      for (const RankedStateDelta &delta : batch_plan.deltas) {
        for (std::size_t i = 0; i < transition_batch.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          if (!exact_apply_ranked_state_delta_row(transition_batch, i, delta)) {
            deactivate_transition_point(i);
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
      exact_collapse_active_scenario_lane_view(transition_batch, active,
                                               transition_weights);
      if (!std::forward<ProcessBatchFn>(process_batch)(transition_batch, &active,
                                                       transition_weights)) {
        return ExactTransitionExecutionResult::Error;
      }
      emitted_any = true;
    }
  }

  return emitted_any ? ExactTransitionExecutionResult::Emitted
                     : ExactTransitionExecutionResult::NoEmission;
}

ExactTransitionExecutionResult exact_execute_compiled_transition_plan_impl(
    const RankedNodeBatchPlan &plan, const uuber::NativeContext &ctx,
    const ExactScenarioBatch &seed_batch, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    ExactTransitionReducer &reducer,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_execute_compiled_transition_plan_process_impl(
      plan, ctx, seed_batch, component_idx, trial_params,
      trial_type_key,
      [&](const auto &points, const std::vector<std::uint8_t> *active_mask,
          const std::vector<double> &weights) -> bool {
        return reducer.consume(points, active_mask, weights);
      },
      uniform_trial_params_soa);
}

inline bool exact_competitor_guard_fastpath_eligible(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache) {
  if (competitor_cache.compiled_ops.empty()) {
    return false;
  }
  for (const uuber::TreeCompetitorOp &op : competitor_cache.compiled_ops) {
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
  const uuber::TreeCompetitorOp &op = competitor_cache.compiled_ops.front();
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
inline void exact_competitor_init_active_mask(const PointVec &points,
                                              std::vector<std::uint8_t> &active) {
  active.assign(points.size(), 0u);
  for (std::size_t i = 0; i < points.size(); ++i) {
    active[i] = (std::isfinite(points.t[i]) && points.t[i] >= 0.0) ? 1u : 0u;
  }
}

inline void exact_prepare_mutable_competitor_points(
    const ExactScenarioBatch &points, ExactScenarioBatch &mutable_points) {
  mutable_points = points;
}

inline void exact_prepare_mutable_competitor_points(
    const uuber::ExactScenarioLaneViewBatch &points,
    uuber::ExactScenarioLaneViewBatch &mutable_points) {
  mutable_points = points;
  mutable_points.mutable_forced_survive_bits.clear();
  mutable_points.mutable_forced_survive_bits.reserve(points.size());
  mutable_points.forced_survive_bits.clear();
  mutable_points.forced_survive_bits.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i < points.forced_survive_bits.size() &&
        points.forced_survive_bits[i] != nullptr) {
      mutable_points.mutable_forced_survive_bits.push_back(
          *points.forced_survive_bits[i]);
    } else {
      mutable_points.mutable_forced_survive_bits.emplace_back();
    }
    mutable_points.forced_survive_bits.push_back(
        &mutable_points.mutable_forced_survive_bits.back());
  }
}

template <typename PointVec>
inline const uuber::TrialParamsSoA *exact_resolve_competitor_point_soa(
    const uuber::NativeContext &ctx, const PointVec &points, std::size_t idx,
    const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return uuber::resolve_vector_lane_trial_params_soa(
      ctx, uuber::vector_lane_ref(points, idx), trial_params,
      uniform_trial_params_soa);
}

struct ExactTimeTrialParamsSoAKey {
  std::uint64_t time_bits{0u};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};

  bool operator==(const ExactTimeTrialParamsSoAKey &other) const noexcept {
    return time_bits == other.time_bits &&
           trial_params_soa == other.trial_params_soa;
  }
};

struct ExactTimeTrialParamsSoAKeyHash {
  std::size_t operator()(const ExactTimeTrialParamsSoAKey &key) const noexcept {
    std::uint64_t hash = kFNV64Offset;
    hash_append_u64(hash, key.time_bits);
    hash_append_u64(
        hash, static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(
                  key.trial_params_soa)));
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

template <typename PointVec>
bool exact_eval_node_survival_batch_compressed(
    const uuber::NativeContext &ctx, int node_idx, const PointVec &points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const std::vector<std::uint8_t> &active,
    const uuber::TrialParamsSoA *uniform_trial_params_soa,
    uuber::TreeNodeBatchValues &node_values) {
  const std::size_t point_count = points.size();
  const std::size_t active_count =
      uuber::count_active_vector_lanes(&active, point_count);
  if (active_count == 0u || active_count < 8u) {
    return exact_eval_node_batch_with_points(
        ctx, node_idx, points, component_idx, trial_params, trial_type_key,
        EvalNeed::kSurvival, node_values, &active, uniform_trial_params_soa);
  }

  const std::vector<int> relevant_source_ids =
      relevant_source_ids_for_batch_node(ctx, node_idx);
  std::size_t anchor_idx = point_count;
  if (!uuber::vector_lanes_share_relevant_state(points, &active,
                                                relevant_source_ids, ctx,
                                                anchor_idx)) {
    return exact_eval_node_batch_with_points(
        ctx, node_idx, points, component_idx, trial_params, trial_type_key,
        EvalNeed::kSurvival, node_values, &active, uniform_trial_params_soa);
  }

  std::unordered_map<ExactTimeTrialParamsSoAKey, std::size_t,
                     ExactTimeTrialParamsSoAKeyHash>
      index_by_key;
  index_by_key.reserve(active_count);
  std::vector<std::size_t> inverse_indices(point_count,
                                           std::numeric_limits<std::size_t>::max());
  uuber::ExactScenarioLaneViewBatch unique_points;
  unique_points.reserve(active_count);
  bool any_duplicate = false;

  for (std::size_t i = 0; i < point_count; ++i) {
    if (!uuber::vector_lane_active(&active, i)) {
      continue;
    }
    const uuber::TrialParamsSoA *point_soa =
        exact_resolve_competitor_point_soa(ctx, points, i, trial_params,
                                           uniform_trial_params_soa);
    if (point_soa == nullptr) {
      return false;
    }
    const ExactTimeTrialParamsSoAKey key{canonical_double_bits(points.t[i]),
                                         point_soa};
    auto it = index_by_key.find(key);
    if (it != index_by_key.end()) {
      inverse_indices[i] = it->second;
      any_duplicate = true;
      continue;
    }
    const std::size_t unique_idx = unique_points.size();
    inverse_indices[i] = unique_idx;
    index_by_key.emplace(key, unique_idx);
    exact_scenario_view_batch_append_lane(unique_points,
                                          uuber::vector_lane_ref(points, i),
                                          uuber::vector_lane_ref(points, i).weight,
                                          point_soa);
  }

  if (!any_duplicate || unique_points.size() >= active_count) {
    return exact_eval_node_batch_with_points(
        ctx, node_idx, points, component_idx, trial_params, trial_type_key,
        EvalNeed::kSurvival, node_values, &active, uniform_trial_params_soa);
  }

  uuber::TreeNodeBatchValues unique_values;
  if (!exact_eval_node_batch_with_points(ctx, node_idx, unique_points,
                                         component_idx, trial_params,
                                         trial_type_key, EvalNeed::kSurvival,
                                         unique_values, nullptr, nullptr) ||
      unique_values.survival.size() != unique_points.size()) {
    return false;
  }

  node_values = uuber::TreeNodeBatchValues{};
  node_values.density.assign(point_count, 0.0);
  node_values.survival.assign(point_count, 1.0);
  node_values.cdf.assign(point_count, 0.0);
  for (std::size_t i = 0; i < point_count; ++i) {
    if (!uuber::vector_lane_active(&active, i)) {
      continue;
    }
    const std::size_t unique_idx = inverse_indices[i];
    if (unique_idx >= unique_values.survival.size()) {
      return false;
    }
    node_values.survival[i] = unique_values.survival[unique_idx];
  }
  return true;
}

template <typename PointVec>
bool exact_eval_guard_survival_batch_direct(
    const uuber::NativeContext &ctx, int node_idx, const PointVec &points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const std::vector<std::uint8_t> &active,
    const uuber::TrialParamsSoA *uniform_trial_params_soa,
    uuber::TreeNodeBatchValues &node_values) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  if (node.op != uuber::IrNodeOp::Guard) {
    return false;
  }

  ExactNodeBatchDepthGuard depth_guard;
  const std::size_t active_point_count =
      uuber::count_active_vector_lanes(&active, points.size());
  (void)depth_guard;

  exact_init_batch_values(points.size(), node_values);
  if (active_point_count == 0u) {
    return true;
  }

  const std::vector<int> relevant_source_ids =
      relevant_source_ids_for_batch_node(ctx, node_idx);
  return uuber::eval_vector_lanes_with_shared_state(
      ctx, points, &active, relevant_source_ids,
      [&](std::size_t /*group_anchor_idx*/,
          const uuber::VectorRelevantStateStorage &group_storage,
          const std::vector<std::uint8_t> *group_mask) -> bool {
        std::vector<const uuber::TrialParamsSoA *> subgroup_params;
        std::vector<std::vector<std::size_t>> subgroup_indices;
        subgroup_params.reserve(active_point_count);
        subgroup_indices.reserve(active_point_count);

        for (std::size_t i = 0; i < points.size(); ++i) {
          if (!uuber::vector_lane_active(group_mask, i)) {
            continue;
          }
          const uuber::TrialParamsSoA *point_soa =
              uniform_trial_params_soa != nullptr
                  ? uniform_trial_params_soa
                  : exact_resolve_competitor_point_soa(
                        ctx, points, i, trial_params, uniform_trial_params_soa);
          if (point_soa == nullptr) {
            return false;
          }
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

        const bool forced_complete_bits_valid =
            group_storage.forced_complete_bits_valid;
        const bool forced_survive_bits_valid =
            group_storage.forced_survive_bits_valid;
        const ForcedStateView forced_state = make_forced_state_view(
            nullptr,
            forced_complete_bits_valid ? &group_storage.forced_complete_bits
                                       : nullptr,
            &forced_complete_bits_valid,
            forced_survive_bits_valid ? &group_storage.forced_survive_bits
                                      : nullptr,
            &forced_survive_bits_valid,
            ctx.ir.label_id_to_bit_idx.empty() ? nullptr
                                               : &ctx.ir.label_id_to_bit_idx);

        for (std::size_t subgroup_idx = 0; subgroup_idx < subgroup_params.size();
             ++subgroup_idx) {
          const std::vector<std::size_t> &indices =
              subgroup_indices[subgroup_idx];
          if (indices.empty()) {
            continue;
          }
          std::vector<double> subgroup_times(indices.size(), 0.0);
          for (std::size_t j = 0; j < indices.size(); ++j) {
            subgroup_times[j] = points.t[indices[j]];
          }
          const GuardEvalInput guard_input = make_guard_input_forced_state(
              ctx, node_idx, component_idx, &trial_type_key, trial_params,
              subgroup_params[subgroup_idx], &group_storage.time_constraints,
              forced_state);
          std::vector<double> subgroup_cdf;
          if (!evaluator_guard_cdf_batch_prepared(guard_input, subgroup_times,
                                                  subgroup_cdf) ||
              subgroup_cdf.size() != indices.size()) {
            return false;
          }
          for (std::size_t j = 0; j < indices.size(); ++j) {
            const std::size_t point_idx = indices[j];
            node_values.cdf[point_idx] = clamp_probability(subgroup_cdf[j]);
            node_values.survival[point_idx] =
                clamp_probability(1.0 - node_values.cdf[point_idx]);
          }
        }
        return true;
      },
      [&]() -> bool {
        return exact_eval_node_batch_with_points(
            ctx, node_idx, points, component_idx, trial_params, trial_type_key,
            EvalNeed::kSurvival, node_values, &active,
            uniform_trial_params_soa);
      },
      [&](std::uint64_t active_count, std::uint64_t group_count) {
        (void)active_count;
        (void)group_count;
      },
      false);
}

inline void exact_competitor_copy_survival_values(
    const std::vector<double> &node_survival, std::vector<double> &survival_out,
    std::size_t point_count, std::size_t compiled_op_count,
    const char *reason) {
  exact_competitor_batch_require(node_survival.size() == survival_out.size(),
                                 reason, point_count, compiled_op_count);
  for (std::size_t i = 0; i < survival_out.size(); ++i) {
    survival_out[i] = clamp_probability(node_survival[i]);
  }
}

inline void exact_competitor_apply_survival_values(
    const std::vector<double> &node_survival, std::vector<std::uint8_t> &active,
    std::vector<double> &survival_out, std::size_t point_count,
    std::size_t compiled_op_count, const char *reason) {
  exact_competitor_batch_require(node_survival.size() == survival_out.size(),
                                 reason, point_count, compiled_op_count);
  for (std::size_t i = 0; i < survival_out.size(); ++i) {
    if (active[i] == 0u) {
      continue;
    }
    const double surv = clamp_probability(node_survival[i]);
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
}

template <typename PointVec>
inline void exact_competitor_apply_transition_masks(
    const uuber::NativeContext &ctx, PointVec *eval_points_ptr,
    const std::vector<std::uint8_t> &active,
    const uuber::TreeCompetitorOp &op) {
  if (eval_points_ptr == nullptr || op.transition_mask_begin < 0 ||
      op.transition_mask_count <= 0) {
    return;
  }
  for (std::size_t i = 0; i < eval_points_ptr->size(); ++i) {
    if (active[i] == 0u) {
      continue;
    }
    exact_apply_transition_mask_to_point(ctx, *eval_points_ptr, i, op);
  }
}

inline void exact_finalize_density_output(std::vector<double> &density_out) {
  for (double &value : density_out) {
    value = safe_density(value);
  }
}

template <typename PointVec>
inline void exact_accumulate_scenario_density_raw(
    const PointVec &scenario_points, std::vector<double> &density_out,
    const std::vector<double> *survival_out = nullptr) {
  for (std::size_t i = 0; i < scenario_points.size(); ++i) {
    const uuber::VectorLaneRef point = uuber::vector_lane_ref(scenario_points, i);
    double weight = point.weight;
    if (survival_out != nullptr) {
      const double surv = clamp_probability((*survival_out)[i]);
      if (!std::isfinite(surv) || surv <= 0.0) {
        continue;
      }
      weight *= surv;
    }
    if (!std::isfinite(weight) || weight <= 0.0 ||
        point.density_index >= density_out.size()) {
      continue;
    }
    density_out[point.density_index] += weight;
  }
}

template <typename PointVec>
inline void exact_accumulate_scenario_density(
    const PointVec &scenario_points, std::vector<double> &density_out,
    const std::vector<double> *survival_out = nullptr) {
  exact_accumulate_scenario_density_raw(scenario_points, density_out,
                                        survival_out);
  exact_finalize_density_output(density_out);
}

template <typename PointVec>
inline void exact_flush_competitor_scenario_chunk(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, PointVec &scenario_chunk,
    std::vector<double> &density_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  if (scenario_chunk.empty()) {
    return;
  }
  std::vector<double> survival_out;
  exact_competitor_survival_batch_impl(ctx, competitor_cache, component_idx,
                                       trial_params, trial_type_key,
                                       scenario_chunk, survival_out,
                                       uniform_trial_params_soa);
  exact_accumulate_scenario_density_raw(scenario_chunk, density_out,
                                        &survival_out);
  scenario_chunk.clear();
}

class ExactDensityAccumulatorReducer final : public ExactTransitionReducer {
public:
  ExactDensityAccumulatorReducer(
      const uuber::NativeContext &ctx_,
      const uuber::CompetitorClusterCacheEntry *competitor_cache_,
      int component_idx_, const TrialParamSet *trial_params_,
      const std::string &trial_type_key_, std::vector<double> &density_out_,
      const uuber::TrialParamsSoA *uniform_trial_params_soa_)
      : ctx(ctx_), competitor_cache(competitor_cache_),
        component_idx(component_idx_), trial_params(trial_params_),
        trial_type_key(trial_type_key_), density_out(density_out_),
        uniform_trial_params_soa(uniform_trial_params_soa_) {
    if (competitor_cache != nullptr) {
      constexpr std::size_t kExactCompetitorChunkSize = 512u;
      scenario_chunk.reserve(std::min<std::size_t>(
          std::max<std::size_t>(density_out.size(), 1u),
          kExactCompetitorChunkSize));
    }
  }

  bool consume(const ExactScenarioBatch &points,
               const std::vector<std::uint8_t> *active_mask,
               const std::vector<double> &weights) override {
    return consume_impl(points, active_mask, weights);
  }

  bool consume(const uuber::ExactScenarioLaneViewBatch &points,
               const std::vector<std::uint8_t> *active_mask,
               const std::vector<double> &weights) override {
    return consume_impl(points, active_mask, weights);
  }

  void finish() {
    if (emitted_any) {
      exact_finalize_density_output(density_out);
    }
  }

private:
  template <typename PointVec>
  bool consume_impl(const PointVec &points,
                    const std::vector<std::uint8_t> *active_mask,
                    const std::vector<double> &weights) {
    if (competitor_cache == nullptr) {
      for (std::size_t i = 0; i < points.size(); ++i) {
        if (!uuber::vector_lane_active(active_mask, i) || i >= weights.size() ||
            !std::isfinite(weights[i]) || weights[i] <= 0.0) {
          continue;
        }
        const uuber::VectorLaneRef point = uuber::vector_lane_ref(points, i);
        if (point.density_index >= density_out.size()) {
          continue;
        }
        density_out[point.density_index] += weights[i];
        emitted_any = true;
      }
      return true;
    }

    constexpr std::size_t kExactCompetitorChunkSize = 512u;
    scenario_chunk.clear();
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (!uuber::vector_lane_active(active_mask, i) || i >= weights.size() ||
          !std::isfinite(weights[i]) || weights[i] <= 0.0) {
        continue;
      }
      exact_scenario_view_batch_append_lane(
          scenario_chunk, uuber::vector_lane_ref(points, i), weights[i]);
      emitted_any = true;
      if (scenario_chunk.size() >= kExactCompetitorChunkSize) {
        exact_flush_competitor_scenario_chunk(
            ctx, *competitor_cache, component_idx, trial_params, trial_type_key,
            scenario_chunk, density_out, uniform_trial_params_soa);
      }
    }
    exact_flush_competitor_scenario_chunk(
        ctx, *competitor_cache, component_idx, trial_params, trial_type_key,
        scenario_chunk, density_out, uniform_trial_params_soa);
    return true;
  }

  const uuber::NativeContext &ctx;
  const uuber::CompetitorClusterCacheEntry *competitor_cache;
  int component_idx;
  const TrialParamSet *trial_params;
  const std::string &trial_type_key;
  std::vector<double> &density_out;
  const uuber::TrialParamsSoA *uniform_trial_params_soa;
  uuber::ExactScenarioLaneViewBatch scenario_chunk;
  bool emitted_any{false};
};

template <typename PointVec>
void exact_competitor_survival_batch_impl(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const PointVec &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  survival_out.assign(points.size(), 1.0);
  if (points.empty() || competitor_cache.compiled_ops.empty()) {
    return;
  }
  PointVec mutable_points;
  PointVec *eval_points_ptr = nullptr;
  if (competitor_cache.mutates_forced_survive) {
    exact_prepare_mutable_competitor_points(points, mutable_points);
    eval_points_ptr = &mutable_points;
  }
  const PointVec &eval_points = eval_points_ptr ? *eval_points_ptr : points;
  std::vector<std::uint8_t> active;
  exact_competitor_init_active_mask(eval_points, active);

  if (exact_competitor_single_node_fastpath_eligible(competitor_cache)) {
    uuber::TreeNodeBatchValues node_values;
    const int node_idx =
        competitor_cache.compiled_ops.front().target_node_indices.front();
    exact_competitor_batch_require(
        ((node_idx >= 0 &&
          node_idx < static_cast<int>(ctx.ir.nodes.size()) &&
          ctx.ir.nodes[static_cast<std::size_t>(node_idx)].op ==
              uuber::IrNodeOp::Guard)
             ? exact_eval_guard_survival_batch_direct(
                   ctx, node_idx, eval_points, component_idx, trial_params,
                   trial_type_key, active, uniform_trial_params_soa,
                   node_values)
             : exact_eval_node_survival_batch_compressed(
                   ctx, node_idx, eval_points, component_idx, trial_params,
                   trial_type_key, active, uniform_trial_params_soa,
                   node_values)),
        "single-node exact batch evaluation failed", points.size(),
        competitor_cache.compiled_ops.size());
    exact_competitor_copy_survival_values(
        node_values.survival, survival_out, points.size(),
        competitor_cache.compiled_ops.size(),
        "single-node exact batch size mismatch");
    return;
  }

  if (exact_competitor_guard_fastpath_eligible(ctx, competitor_cache)) {
    uuber::TreeNodeBatchValues guard_values;
    for (const uuber::TreeCompetitorOp &op : competitor_cache.compiled_ops) {
      exact_competitor_batch_require(
          exact_eval_guard_survival_batch_direct(
              ctx, op.target_node_indices.front(), eval_points, component_idx,
              trial_params, trial_type_key, active, uniform_trial_params_soa,
              guard_values),
          "guard competitor exact batch evaluation failed", points.size(),
          competitor_cache.compiled_ops.size());
      exact_competitor_apply_survival_values(
          guard_values.survival, active, survival_out, points.size(),
          competitor_cache.compiled_ops.size(),
          "guard competitor exact batch size mismatch");
      exact_competitor_apply_transition_masks(ctx, eval_points_ptr, active, op);
    }
    return;
  }

  uuber::TreeNodeBatchValues node_values;
  for (const uuber::TreeCompetitorOp &op : competitor_cache.compiled_ops) {
    for (int target_node_idx : op.target_node_indices) {
      exact_competitor_batch_require(
          exact_eval_node_survival_batch_compressed(
              ctx, target_node_idx, eval_points, component_idx, trial_params,
              trial_type_key, active, uniform_trial_params_soa, node_values),
          "generic tree vector competitor batch evaluation failed", points.size(),
          competitor_cache.compiled_ops.size());
      exact_competitor_apply_survival_values(
          node_values.survival, active, survival_out, points.size(),
          competitor_cache.compiled_ops.size(),
          "generic tree vector competitor batch size mismatch");
    }
    exact_competitor_apply_transition_masks(ctx, eval_points_ptr, active, op);
  }
}

} // namespace

int exact_resolve_transition_node_idx(const uuber::NativeContext &ctx,
                                      int node_idx, int component_idx) {
  return exact_resolve_transition_node_idx_impl(ctx, node_idx, component_idx);
}

namespace uuber {

} // namespace uuber

bool exact_eval_node_batch_from_batch(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_eval_node_batch_with_points(ctx, node_idx, points, component_idx,
                                           trial_params, trial_type_key, need,
                                           out_values, nullptr,
                                           uniform_trial_params_soa);
}

ExactTransitionExecutionResult exact_execute_compiled_transition_plan_from_batch(
    const RankedNodeBatchPlan &plan, const uuber::NativeContext &ctx,
    const ExactScenarioBatch &seed_batch, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    ExactTransitionReducer &reducer,
    const uuber::TrialParamsSoA *uniform_trial_params_soa) {
  return exact_execute_compiled_transition_plan_impl(
      plan, ctx, seed_batch, component_idx, trial_params,
      trial_type_key, reducer, uniform_trial_params_soa);
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

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const uuber::ExactScenarioLaneViewBatch &points,
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
  density_out.assign(times.size(), 0.0);
  if (times.empty()) {
    return true;
  }
  RankedBatchPlanner &planner = ranked_transition_planner_for_ctx(ctx);
  const int exec_node_idx =
      exact_resolve_transition_node_idx_impl(ctx, node_idx, state.component_idx);
  if (exec_node_idx == NA_INTEGER) {
    return true;
  }
  const RankedNodeBatchPlan &transition_plan = planner.plan_for_node(exec_node_idx);
  const double saved_t = state.t;
  ExactScenarioBatch seed_batch;
  const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr;
  exact_build_seed_batch_from_state(times, state, seed_batch,
                                    uniform_trial_params_soa);

  const uuber::CompetitorClusterCacheEntry *cache_ptr = nullptr;
  if (!competitor_ids.empty()) {
    cache_ptr = competitor_cache != nullptr
                    ? competitor_cache
                    : &fetch_competitor_cluster_cache(ctx, competitor_ids);
  }
  ExactDensityAccumulatorReducer reducer(
      ctx, cache_ptr, state.component_idx, state.trial_params,
      state.trial_type_key, density_out, uniform_trial_params_soa);
  const ExactTransitionExecutionResult result =
      exact_execute_compiled_transition_plan_impl(
          transition_plan, ctx, seed_batch, state.component_idx,
      state.trial_params, state.trial_type_key, reducer,
      uniform_trial_params_soa);
  if (result == ExactTransitionExecutionResult::Error) {
    state.t = saved_t;
    return false;
  }
  if (result == ExactTransitionExecutionResult::Emitted) {
    reducer.finish();
  }
  state.t = saved_t;
  return true;
}
