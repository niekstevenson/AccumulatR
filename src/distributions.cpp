// [[Rcpp::depends(Rcpp, BH)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "accumulator.h"
#include "bitset_state.h"
#include "competitor_cache.h"
#include "context.h"
#include "distribution_core.h"
#include "dist_vector.h"
#include "exact_outcome_density.h"
#include "evaluator_internal.h"
#include "forced_state.h"
#include "integrate.h"
#include "native_utils.h"
#include "pool_math.h"
#include "proto.h"
#include "quadrature_batch.h"
#include "ranked_transitions.h"
#include "trial_params.h"

using uuber::AccDistParams;
using uuber::ComponentMap;
using uuber::LabelRef;

namespace {

constexpr int kOutcomeIdxNA = -2;

inline int component_index_of(const uuber::NativeContext &ctx,
                              const std::string &component) {
  if (component.empty() || component == "__default__")
    return -1;
  auto it = ctx.component_index.find(component);
  if (it == ctx.component_index.end())
    return -1;
  return it->second;
}

inline const std::string &
component_label_by_index_or_empty(const uuber::NativeContext &ctx,
                                  int component_idx) {
  static const std::string kEmptyComponentLabel;
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.components.ids.size())) {
    return kEmptyComponentLabel;
  }
  return ctx.components.ids[static_cast<std::size_t>(component_idx)];
}

inline bool ir_mask_has_component(const uuber::NativeContext &ctx,
                                  int mask_offset, int component_idx) {
  if (component_idx < 0)
    return true;
  if (!ctx.ir.valid || ctx.ir.component_mask_words <= 0 || mask_offset < 0) {
    return true;
  }
  int word_idx = component_idx / 64;
  int bit = component_idx % 64;
  if (word_idx < 0 || word_idx >= ctx.ir.component_mask_words) {
    return false;
  }
  std::size_t idx = static_cast<std::size_t>(mask_offset + word_idx);
  if (idx >= ctx.ir.component_masks.size())
    return false;
  std::uint64_t word = ctx.ir.component_masks[idx];
  return (word & (1ULL << bit)) != 0ULL;
}

inline bool ir_outcome_allows_component(const uuber::NativeContext &ctx,
                                        int outcome_idx, int component_idx) {
  if (!ctx.ir.valid || outcome_idx < 0 ||
      outcome_idx >= static_cast<int>(ctx.ir.outcomes.size())) {
    return true;
  }
  const uuber::IrOutcome &out =
      ctx.ir.outcomes[static_cast<std::size_t>(outcome_idx)];
  return ir_mask_has_component(ctx, out.allowed_component_mask_offset,
                               component_idx);
}

inline bool
outcome_allows_component_idx(const uuber::NativeContext &ctx,
                             const uuber::OutcomeContextInfo & /*info*/,
                             int outcome_idx, int component_idx) {
  return ir_outcome_allows_component(ctx, outcome_idx, component_idx);
}

inline bool competitor_node_allowed_idx(const uuber::NativeContext &ctx,
                                        int node_id, int component_idx) {
  if (component_idx < 0)
    return true;
  int dense_node = -1;
  auto node_it = ctx.ir.id_to_node_idx.find(node_id);
  if (node_it != ctx.ir.id_to_node_idx.end()) {
    dense_node = node_it->second;
  } else if (node_id >= 0 &&
             node_id < static_cast<int>(ctx.ir.nodes.size())) {
    dense_node = node_id;
  }
  if (dense_node < 0)
    return true;
  if (ctx.ir.valid &&
      static_cast<std::size_t>(dense_node) < ctx.ir.nodes.size()) {
    const uuber::IrNode &node =
        ctx.ir.nodes[static_cast<std::size_t>(dense_node)];
    if (node.component_mask_offset >= 0) {
      return ir_mask_has_component(ctx, node.component_mask_offset,
                                   component_idx);
    }
  }
  auto out_it = ctx.ir.node_idx_to_outcomes.find(dense_node);
  if (out_it == ctx.ir.node_idx_to_outcomes.end() || out_it->second.empty())
    return true;
  for (int out_idx : out_it->second) {
    if (ir_outcome_allows_component(ctx, out_idx, component_idx))
      return true;
  }
  return false;
}

inline const std::vector<int> &
filter_competitor_ids(const uuber::NativeContext &ctx,
                      const std::vector<int> &competitor_ids,
                      int component_idx, std::vector<int> &scratch) {
  if (competitor_ids.empty() || component_idx < 0)
    return competitor_ids;
  bool all_allowed = true;
  for (int node_id : competitor_ids) {
    if (!competitor_node_allowed_idx(ctx, node_id, component_idx)) {
      all_allowed = false;
      break;
    }
  }
  if (all_allowed)
    return competitor_ids;
  scratch.clear();
  scratch.reserve(competitor_ids.size());
  for (int node_id : competitor_ids) {
    if (competitor_node_allowed_idx(ctx, node_id, component_idx)) {
      scratch.push_back(node_id);
    }
  }
  return scratch;
}

inline bool state_allows_exact_batch_dispatch(const uuber::NativeContext &ctx,
                                              const NodeEvalState &state) {
  if (state.include_na_donors || state.outcome_idx < 0 ||
      state.outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return !state.include_na_donors;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(state.outcome_idx)];
  return info.alias_sources.empty() && info.guess_donors.empty();
}

inline bool shared_trigger_mask_batch_supported(
    const uuber::NativeContext &ctx, bool include_na_donors,
    int outcome_idx_context) {
  if (include_na_donors) {
    return false;
  }
  if (outcome_idx_context < 0 ||
      outcome_idx_context >= static_cast<int>(ctx.outcome_info.size())) {
    return true;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(outcome_idx_context)];
  return info.alias_sources.empty() && info.guess_donors.empty();
}

inline bool shared_trigger_density_batch_supported(
    const uuber::NativeContext & /*ctx*/, bool /*include_na_donors*/,
    int /*outcome_idx_context*/) {
  return true;
}

inline bool simple_direct_outcome_fastpath_supported(
    const uuber::NativeContext &ctx, bool include_na_donors,
    int outcome_idx_context) {
  if (include_na_donors) {
    return false;
  }
  if (outcome_idx_context < 0 ||
      outcome_idx_context >= static_cast<int>(ctx.outcome_info.size())) {
    return true;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(outcome_idx_context)];
  if (!info.alias_sources.empty()) {
    return false;
  }
  for (const auto &donor : info.guess_donors) {
    if (donor.rt_policy != "keep" || donor.outcome_idx < 0 ||
        donor.outcome_idx >= static_cast<int>(ctx.outcome_info.size()) ||
        donor.outcome_idx == outcome_idx_context ||
        !shared_trigger_mask_batch_supported(ctx, false, donor.outcome_idx)) {
      return false;
    }
  }
  return true;
}

inline std::uint64_t outcome_coupling_mix_lookup_hash(std::uint64_t hash,
                                                      int value) {
  const std::uint32_t v = static_cast<std::uint32_t>(value);
  for (int i = 0; i < 4; ++i) {
    const std::uint8_t b =
        static_cast<std::uint8_t>((v >> (8 * i)) & 0xFFU);
    hash ^= static_cast<std::uint64_t>(b);
    hash *= kFNV64Prime;
  }
  return hash;
}

inline std::uint64_t outcome_coupling_lookup_key(
    int node_idx, const std::vector<int> &competitors) {
  std::uint64_t hash = kFNV64Offset;
  hash = outcome_coupling_mix_lookup_hash(hash, node_idx);
  hash = outcome_coupling_mix_lookup_hash(hash,
                                          static_cast<int>(competitors.size()));
  for (int comp : competitors) {
    hash = outcome_coupling_mix_lookup_hash(hash, comp);
  }
  return hash;
}

inline bool ir_node_source_masks_overlap(const uuber::NativeContext &ctx,
                                         const uuber::IrNode &a,
                                         const uuber::IrNode &b) {
  if (!ctx.ir.valid || a.source_mask_begin < 0 || a.source_mask_count <= 0 ||
      b.source_mask_begin < 0 || b.source_mask_count <= 0 ||
      a.source_mask_count != b.source_mask_count) {
    return false;
  }
  for (int i = 0; i < a.source_mask_count; ++i) {
    const std::uint64_t wa = ctx.ir.node_source_masks[static_cast<std::size_t>(
        a.source_mask_begin + i)];
    const std::uint64_t wb = ctx.ir.node_source_masks[static_cast<std::size_t>(
        b.source_mask_begin + i)];
    if ((wa & wb) != 0u) {
      return true;
    }
  }
  return false;
}

inline bool generic_overlap_requires_exact_scenario_eval(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids) {
  if (!ctx.ir.valid || node_idx < 0 ||
      node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  const uuber::IrNode &target = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  for (int comp_node_id : competitor_ids) {
    const int comp_idx = resolve_dense_node_idx(ctx, comp_node_id);
    if (comp_idx < 0 || comp_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      continue;
    }
    const uuber::IrNode &competitor =
        ctx.ir.nodes[static_cast<std::size_t>(comp_idx)];
    if (ir_node_source_masks_overlap(ctx, target, competitor)) {
      return true;
    }
  }
  return false;
}

inline const uuber::IrOutcomeCouplingOp *find_outcome_coupling_spec(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids) {
  if (!ctx.ir.valid || node_idx < 0 ||
      node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return nullptr;
  }
  std::vector<int> competitor_dense;
  competitor_dense.reserve(competitor_ids.size());
  for (int comp_node_id : competitor_ids) {
    const int comp_idx = resolve_dense_node_idx(ctx, comp_node_id);
    if (comp_idx < 0) {
      return nullptr;
    }
    competitor_dense.push_back(comp_idx);
  }
  const std::uint64_t key =
      outcome_coupling_lookup_key(node_idx, competitor_dense);
  auto it = ctx.ir.outcome_coupling_lookup.find(key);
  if (it == ctx.ir.outcome_coupling_lookup.end()) {
    return nullptr;
  }
  const int spec_idx = it->second;
  if (spec_idx < 0 ||
      spec_idx >= static_cast<int>(ctx.ir.outcome_coupling_ops.size())) {
    return nullptr;
  }
  return &ctx.ir.outcome_coupling_ops[static_cast<std::size_t>(spec_idx)];
}

inline bool should_use_exact_density_program(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, const NodeEvalState &state,
    const uuber::IrOutcomeCouplingOp *coupling_spec) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size()) ||
      !state_allows_exact_batch_dispatch(ctx, state)) {
    return false;
  }
  if (coupling_spec != nullptr &&
      coupling_spec->kind == uuber::IrOutcomeCouplingKind::GenericNodeIntegral) {
    return coupling_spec->requires_exact_scenario_eval;
  }
  return generic_overlap_requires_exact_scenario_eval(ctx, node_idx,
                                                      competitor_ids);
}

inline const uuber::TrialParamsSoA *
trial_params_soa_for_batch_point(const NodeEvalState &state,
                                 std::size_t point_idx) {
  if (state.trial_params_soa_batch != nullptr &&
      point_idx < state.trial_params_soa_batch->size()) {
    return (*state.trial_params_soa_batch)[point_idx];
  }
  return state.trial_params_soa;
}

inline bool try_eval_conditioned_simple_acc_event_batch(
    const NodeEvalState &state, const uuber::IrEvent &event,
    const std::vector<double> &times, EvalNeed need,
    uuber::TreeNodeBatchValues &out_values);

inline bool eval_event_ref_batch(const NodeEvalState &state,
                                 const LabelRef &label_ref,
                                 std::uint32_t node_flags,
                                 const std::vector<double> &times,
                                 EvalNeed need,
                                 uuber::TreeNodeBatchValues &out_values);

inline bool eval_accumulator_base_batch(const NodeEvalState &state,
                                        const LabelRef &label_ref,
                                        std::uint32_t node_flags,
                                        const std::vector<double> &times,
                                        EvalNeed need,
                                        uuber::TreeNodeBatchValues &out_values);

inline bool component_guess_density_applicable(
    const uuber::NativeContext &ctx, int target_label_id,
    int target_outcome_idx, bool target_is_guess, int component_idx);
inline bool outcome_alias_density_applicable(
    const uuber::NativeContext &ctx, int target_outcome_idx,
    bool include_na_donors);
inline bool accumulate_component_guess_density_batch_idx(
    const uuber::NativeContext &ctx, int target_label_id,
    int target_outcome_idx, bool target_is_guess,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::vector<double> &density_out,
    const SharedTriggerPlan *trigger_plan,
    bool use_shared_trigger_eval,
    std::vector<double> *donor_density_scratch);
inline bool accumulate_outcome_alias_density_batch_idx(
    const uuber::NativeContext &ctx, int target_outcome_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, std::vector<double> &density_out,
    const SharedTriggerPlan *trigger_plan,
    bool use_shared_trigger_eval,
    std::vector<double> *donor_density_scratch);
inline bool eval_node_batch_from_guard_input(
    const GuardEvalInput &input, int node_idx, const std::vector<double> &times,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values);
struct DensityBatchWorkspace;
inline bool integrate_outcome_probability_from_density_batch_idx(
    const uuber::NativeContext &ctx, int node_id, double upper,
    int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, double rel_tol, double abs_tol,
    const std::vector<const TrialParamSet *> &trial_params_batch,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context, std::vector<double> &out_probabilities,
    const std::vector<const uuber::TrialParamsSoA *>
        *trial_params_soa_batch_override,
    const std::vector<const PreparedTrialParamsRuntime *>
        *prepared_runtime_batch,
    DensityBatchWorkspace *workspace);

inline bool try_evaluate_simple_overlap_event_batch_values(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    uuber::TreeNodeBatchValues &out_values);

// Shared-blocker specialization only needs simple accumulator metadata.
struct SimpleAccEventBatchInfo {
  const uuber::TrialParamsSoA *base_soa{nullptr};
  int acc_idx{-1};
  double onset_eff{0.0};
  double q0{0.0};
  AccDistParams cfg{};
  LowerBoundTransform lower_bound{};
  bool component_ok{false};
};

inline bool resolve_simple_acc_event_batch_info(
    const uuber::NativeContext &ctx, int node_idx,
    const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    int component_idx, SimpleAccEventBatchInfo &out,
    bool require_simple_outcome_support = true);

inline bool evaluate_overlap_event_batch_values(
    const uuber::NativeContext &ctx, int event_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    uuber::TreeNodeBatchValues &out_values,
    const TimeConstraintMap *time_constraints = nullptr,
    const ForcedStateView *forced_state_view = nullptr) {
  if (event_idx < 0 || event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return false;
  }
  const uuber::IrEvent &event =
      ctx.ir.events[static_cast<std::size_t>(event_idx)];
  if (event.node_idx < 0 ||
      event.node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }

  const bool has_conditioning =
      time_constraints_any(time_constraints) ||
      (forced_state_view != nullptr &&
       (forced_state_complete_valid(*forced_state_view) ||
        forced_state_survive_valid(*forced_state_view)));

  if (has_conditioning) {
    NodeEvalState state(ctx, 0.0, component_idx, trial_params, trial_type_key,
                        false, -1, time_constraints, nullptr, nullptr, false,
                        nullptr, false, forced_state_view);
    state.trial_params_soa_batch = trial_params_soa_batch;
    if (try_eval_conditioned_simple_acc_event_batch(state, event, times,
                                                    EvalNeed::kDensity |
                                                        EvalNeed::kSurvival |
                                                        EvalNeed::kCDF,
                                                    out_values)) {
      return true;
    }
    if (!evaluator_eval_node_batch_with_state_dense(
            event.node_idx, times, state,
            EvalNeed::kDensity | EvalNeed::kSurvival | EvalNeed::kCDF,
            out_values) ||
        out_values.density.size() != times.size() ||
        out_values.survival.size() != times.size() ||
        out_values.cdf.size() != times.size()) {
      return false;
    }
  } else if (try_evaluate_simple_overlap_event_batch_values(
                 ctx, event.node_idx, times, component_idx, trial_params,
                 trial_params_soa_batch, out_values)) {
    return true;
  } else {
    NodeEvalState state(ctx, 0.0, component_idx, trial_params, trial_type_key,
                        false, -1);
    state.trial_params_soa_batch = trial_params_soa_batch;
    if (!evaluator_eval_node_batch_with_state_dense(
            event.node_idx, times, state,
            EvalNeed::kDensity | EvalNeed::kSurvival | EvalNeed::kCDF,
            out_values) ||
        out_values.density.size() != times.size() ||
        out_values.survival.size() != times.size() ||
        out_values.cdf.size() != times.size()) {
      return false;
    }
  }
  for (std::size_t i = 0; i < times.size(); ++i) {
    out_values.density[i] = safe_density(out_values.density[i]);
    out_values.survival[i] = clamp_probability(out_values.survival[i]);
    out_values.cdf[i] = clamp_probability(out_values.cdf[i]);
  }
  return true;
}

struct TimeIntegralPointPlan {
  std::size_t begin{0u};
  std::size_t count{0u};
};

struct DensityBatchWorkspace;

inline double max_positive_finite_upper(const std::vector<double> &upper_bounds) {
  double max_upper = 0.0;
  for (double upper : upper_bounds) {
    if (std::isfinite(upper) && upper > max_upper) {
      max_upper = upper;
    }
  }
  return max_upper;
}

inline bool build_time_integral_query_plan(
    const std::vector<double> &upper_bounds, double lower, double upper,
    int segments, std::vector<TimeIntegralPointPlan> &plans,
    std::vector<double> &query_times, std::vector<double> &query_weights,
    const std::vector<const uuber::TrialParamsSoA *> *point_trial_params_soa,
    std::vector<const uuber::TrialParamsSoA *> *query_trial_params_soa,
    const std::vector<std::uint8_t> *active_points = nullptr) {
  if (!(std::isfinite(lower) && std::isfinite(upper)) || !(upper > lower) ||
      segments < 1) {
    return false;
  }

  plans.assign(upper_bounds.size(), TimeIntegralPointPlan{});
  query_times.clear();
  query_weights.clear();
  if (query_trial_params_soa != nullptr) {
    query_trial_params_soa->clear();
  }
  if (point_trial_params_soa != nullptr &&
      point_trial_params_soa->size() != upper_bounds.size()) {
    return false;
  }

  const double window_width = upper - lower;
  std::size_t estimated_node_count = 0u;
  for (std::size_t i = 0; i < upper_bounds.size(); ++i) {
    if (active_points != nullptr &&
        (i >= active_points->size() || (*active_points)[i] == 0u)) {
      continue;
    }
    const double point_upper = upper_bounds[i];
    const double interval_hi =
        std::isfinite(point_upper) ? std::min(upper, point_upper) : upper;
    if (interval_hi > lower) {
      const double interval_width = interval_hi - lower;
      const int local_segments =
          std::max(1, static_cast<int>(std::ceil(
                          static_cast<double>(segments) *
                          interval_width / std::max(window_width, interval_width))));
      estimated_node_count +=
          static_cast<std::size_t>(local_segments * 15);
    }
  }
  query_times.reserve(estimated_node_count);
  query_weights.reserve(estimated_node_count);
  if (query_trial_params_soa != nullptr) {
    query_trial_params_soa->reserve(estimated_node_count);
  }

  for (std::size_t i = 0; i < upper_bounds.size(); ++i) {
    TimeIntegralPointPlan &plan = plans[i];
    plan.begin = query_times.size();
    if (active_points != nullptr &&
        (i >= active_points->size() || (*active_points)[i] == 0u)) {
      continue;
    }
    const double point_upper = upper_bounds[i];
    const double interval_hi =
        std::isfinite(point_upper) ? std::min(upper, point_upper) : upper;
    if (!(interval_hi > lower)) {
      continue;
    }
    const double interval_width = interval_hi - lower;
    const int local_segments =
        std::max(1, static_cast<int>(std::ceil(
                        static_cast<double>(segments) *
                        interval_width / std::max(window_width, interval_width))));
    const uuber::TrialParamsSoA *point_trial_params =
        point_trial_params_soa != nullptr ? (*point_trial_params_soa)[i] : nullptr;
    if (query_trial_params_soa != nullptr && point_trial_params == nullptr) {
      return false;
    }
    for (int seg = 0; seg < local_segments; ++seg) {
      const double seg_lo =
          lower + interval_width * static_cast<double>(seg) /
                      static_cast<double>(local_segments);
      const double seg_hi =
          lower + interval_width * static_cast<double>(seg + 1) /
                      static_cast<double>(local_segments);
      const uuber::TimeBatch batch = uuber::build_time_batch(seg_lo, seg_hi);
      if (batch.nodes.empty() || batch.nodes.size() != batch.weights.size()) {
        return false;
      }
      query_times.insert(query_times.end(), batch.nodes.begin(), batch.nodes.end());
      query_weights.insert(query_weights.end(), batch.weights.begin(),
                           batch.weights.end());
      if (query_trial_params_soa != nullptr) {
        query_trial_params_soa->insert(query_trial_params_soa->end(),
                                       batch.nodes.size(), point_trial_params);
      }
    }
    plan.count = query_times.size() - plan.begin;
  }
  return true;
}

inline bool integrate_query_plan_density_per_point(
    const std::vector<TimeIntegralPointPlan> &plans,
    const std::vector<double> &query_weights,
    const std::vector<double> &density, std::vector<double> &out) {
  if (query_weights.size() != density.size()) {
    return false;
  }
  out.assign(plans.size(), 0.0);
  for (std::size_t point_idx = 0; point_idx < plans.size(); ++point_idx) {
    const TimeIntegralPointPlan &plan = plans[point_idx];
    if (plan.begin + plan.count > density.size()) {
      return false;
    }
    double total = 0.0;
    for (std::size_t query_idx = plan.begin;
         query_idx < plan.begin + plan.count; ++query_idx) {
      const double weight = query_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      const double value = safe_density(density[query_idx]);
      if (!std::isfinite(value) || value <= 0.0) {
        continue;
      }
      total += weight * value;
    }
    out[point_idx] = safe_density(total);
  }
  return true;
}

inline bool integrate_query_plan_density_per_point_slice(
    const std::vector<TimeIntegralPointPlan> &plans,
    const std::vector<double> &query_weights, const double *density_begin,
    std::size_t density_count, std::vector<double> &out) {
  if (density_begin == nullptr || query_weights.size() != density_count) {
    return false;
  }
  out.assign(plans.size(), 0.0);
  for (std::size_t point_idx = 0; point_idx < plans.size(); ++point_idx) {
    const TimeIntegralPointPlan &plan = plans[point_idx];
    if (plan.begin + plan.count > density_count) {
      return false;
    }
    double total = 0.0;
    for (std::size_t query_idx = plan.begin;
         query_idx < plan.begin + plan.count; ++query_idx) {
      const double weight = query_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      const double value = safe_density(density_begin[query_idx]);
      if (!std::isfinite(value) || value <= 0.0) {
        continue;
      }
      total += weight * value;
    }
    out[point_idx] = safe_density(total);
  }
  return true;
}

inline bool resolve_point_trial_params_soa_batch(
    const uuber::NativeContext &ctx, std::size_t point_count,
    const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<const uuber::TrialParamsSoA *> &storage,
    const std::vector<const uuber::TrialParamsSoA *> *&out_ptr) {
  if (trial_params_soa_batch != nullptr) {
    if (trial_params_soa_batch->size() != point_count) {
      return false;
    }
    for (const uuber::TrialParamsSoA *point_soa : *trial_params_soa_batch) {
      if (point_soa == nullptr || !point_soa->valid) {
        return false;
      }
    }
    out_ptr = trial_params_soa_batch;
    return true;
  }

  const uuber::TrialParamsSoA *default_trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);
  if (default_trial_params_soa == nullptr || !default_trial_params_soa->valid) {
    return false;
  }
  storage.assign(point_count, default_trial_params_soa);
  out_ptr = &storage;
  return true;
}

inline int find_event_idx_for_label_ref(const uuber::NativeContext &ctx,
                                        const LabelRef &label_ref) {
  for (std::size_t event_idx = 0; event_idx < ctx.ir.events.size();
       ++event_idx) {
    const uuber::IrEvent &event = ctx.ir.events[event_idx];
    if (label_ref.label_id >= 0 && event.label_id == label_ref.label_id) {
      return static_cast<int>(event_idx);
    }
    if (label_ref.acc_idx >= 0 && event.acc_idx == label_ref.acc_idx) {
      return static_cast<int>(event_idx);
    }
    if (label_ref.pool_idx >= 0 && event.pool_idx == label_ref.pool_idx) {
      return static_cast<int>(event_idx);
    }
  }
  return -1;
}

inline bool eval_event_ref_batch_with_params(
    const NodeEvalState &state, const LabelRef &label_ref,
    std::uint32_t node_flags, const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values) {
  NodeEvalState local_state = state;
  local_state.trial_params_soa_batch = &trial_params_soa_batch;
  if (!trial_params_soa_batch.empty()) {
    local_state.trial_params_soa = trial_params_soa_batch.front();
  }
  if (!times.empty()) {
    local_state.t = times.front();
  }
  return eval_event_ref_batch(local_state, label_ref, node_flags, times, need,
                              out_values);
}

inline bool eval_accumulator_base_batch_with_params(
    const NodeEvalState &state, const LabelRef &label_ref,
    std::uint32_t node_flags, const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values) {
  NodeEvalState local_state = state;
  local_state.trial_params_soa_batch = &trial_params_soa_batch;
  if (!trial_params_soa_batch.empty()) {
    local_state.trial_params_soa = trial_params_soa_batch.front();
  }
  if (!times.empty()) {
    local_state.t = times.front();
  }
  return eval_accumulator_base_batch(local_state, label_ref, node_flags, times,
                                     need, out_values);
}

inline bool eval_accumulator_base_batch(
    const NodeEvalState &state, const LabelRef &label_ref,
    std::uint32_t node_flags, const std::vector<double> &times,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values) {
  out_values = uuber::TreeNodeBatchValues{};
  const std::size_t point_count = times.size();
  out_values.density.assign(point_count, 0.0);
  out_values.survival.assign(point_count, 1.0);
  out_values.cdf.assign(point_count, 0.0);
  if (point_count == 0u) {
    return true;
  }

  const int acc_idx = label_ref.acc_idx;
  if (acc_idx < 0 || acc_idx >= static_cast<int>(state.ctx.accumulators.size())) {
    return true;
  }
  const uuber::NativeAccumulator &acc =
      state.ctx.accumulators[static_cast<std::size_t>(acc_idx)];
  const TrialAccumulatorParams *override =
      evaluator_get_trial_param_entry(state.trial_params, acc_idx);
  if (!evaluator_component_active_idx(acc, state.component_idx, override)) {
    return true;
  }

  const int onset_kind = override ? override->onset_kind : acc.onset_kind;
  const int onset_source_acc_idx =
      override ? override->onset_source_acc_idx : acc.onset_source_acc_idx;
  const int onset_source_pool_idx =
      override ? override->onset_source_pool_idx : acc.onset_source_pool_idx;

  if (!state.ctx.has_chained_onsets ||
      onset_kind == uuber::ONSET_ABSOLUTE) {
    for (std::size_t i = 0; i < point_count; ++i) {
      double onset = 0.0;
      double q = 0.0;
      AccDistParams cfg;
      evaluator_resolve_event_numeric_params(
          acc, acc_idx, override, trial_params_soa_for_batch_point(state, i),
          onset, q, cfg);
      if (needs_density(need)) {
        out_values.density[i] = acc_density_from_cfg(times[i], onset, q, cfg);
      }
      if (needs_survival(need)) {
        out_values.survival[i] =
            acc_survival_from_cfg(times[i], onset, q, cfg);
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] = acc_cdf_success_from_cfg(times[i], onset, q, cfg);
      }
    }
    return true;
  }

  const LabelRef source_ref = evaluator_make_onset_source_ref(
      state.ctx, onset_kind, onset_source_acc_idx, onset_source_pool_idx);
  if (source_ref.acc_idx < 0 && source_ref.pool_idx < 0) {
    return true;
  }

  const int source_label_id = source_ref.label_id;
  bool source_has_exact = false;
  double source_exact_time = std::numeric_limits<double>::quiet_NaN();
  bool source_has_bounds = false;
  double source_bound_lower = 0.0;
  double source_bound_upper = std::numeric_limits<double>::infinity();
  if (source_label_id >= 0 && source_label_id != NA_INTEGER) {
    (void)evaluator_resolve_label_time_constraint(
        &state.time_constraints, source_label_id, source_has_exact,
        source_exact_time, source_has_bounds, source_bound_lower,
        source_bound_upper);
  }

  struct ChainedPointMeta {
    const uuber::TrialParamsSoA *trial_params_soa{nullptr};
    double t{0.0};
    double q{0.0};
    double x_shift{0.0};
    double bound_lower{0.0};
    double bound_upper{std::numeric_limits<double>::infinity()};
    double upper_int{0.0};
    double source_mass{1.0};
    AccDistParams cfg{};
    LowerBoundTransform lower_bound{};
    bool has_conditioning_bounds{false};
    bool needs_source_finish{false};
    bool needs_integral{false};
  };

  std::vector<ChainedPointMeta> meta(point_count);
  std::vector<std::uint8_t> active_integral(point_count, 0u);
  std::vector<double> source_lower_cdf(point_count, 0.0);
  std::vector<double> source_upper_cdf(point_count, 0.0);
  std::vector<double> source_finish_cdf(point_count, 0.0);

  for (std::size_t i = 0; i < point_count; ++i) {
    ChainedPointMeta &point = meta[i];
    point.trial_params_soa = trial_params_soa_for_batch_point(state, i);
    point.t = times[i];

    double onset = 0.0;
    evaluator_resolve_event_numeric_params(
        acc, acc_idx, override, point.trial_params_soa, onset, point.q,
        point.cfg);
    point.lower_bound = default_lower_bound_transform(point.cfg);
    point.x_shift = onset + (override ? override->onset_lag : acc.onset_lag) +
                    point.cfg.t0;

    if (source_label_id >= 0 && source_label_id != NA_INTEGER &&
        evaluator_state_contains_survive_at(state.forced_state,
                                            &state.time_constraints,
                                            source_label_id, point.t)) {
      continue;
    }

    point.bound_lower = 0.0;
    point.bound_upper = std::numeric_limits<double>::infinity();
    if (source_label_id >= 0 && source_label_id != NA_INTEGER) {
      if (forced_state_contains_complete(state.forced_state, source_label_id)) {
        point.bound_upper = std::min(point.bound_upper, point.t);
        point.has_conditioning_bounds = true;
      }
      if (source_has_bounds) {
        point.bound_lower = std::max(point.bound_lower, source_bound_lower);
        point.bound_upper = std::min(point.bound_upper, source_bound_upper);
        point.has_conditioning_bounds = true;
      }
    }

    if (source_has_exact) {
      const double x = point.t - source_exact_time - point.x_shift;
      const double cdf_success = clamp_probability(
          eval_cdf_single_with_lower_bound(point.cfg, x, point.lower_bound));
      if (needs_density(need)) {
        out_values.density[i] =
            safe_density((1.0 - point.q) *
                         eval_pdf_single_with_lower_bound(point.cfg, x,
                                                          point.lower_bound));
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] =
            clamp_probability((1.0 - point.q) * cdf_success);
      }
      if (needs_survival(need)) {
        out_values.survival[i] =
            clamp_probability(point.q + (1.0 - point.q) * (1.0 - cdf_success));
      }
      continue;
    }

    if (!(point.bound_upper > point.bound_lower)) {
      continue;
    }
    if (std::isfinite(point.t) &&
        (point.t - point.x_shift <= point.bound_lower)) {
      continue;
    }

    if (!std::isfinite(point.t)) {
      if (point.has_conditioning_bounds) {
        const double cdf_total = clamp_probability(1.0 - point.q);
        if (needs_cdf(need)) {
          out_values.cdf[i] = cdf_total;
        }
        if (needs_survival(need)) {
          out_values.survival[i] = clamp_probability(1.0 - cdf_total);
        }
      } else {
        point.needs_source_finish = needs_cdf(need) || needs_survival(need);
      }
      continue;
    }

    point.upper_int = std::min(
        point.bound_upper,
        std::max(point.bound_lower, point.t - point.x_shift));
    if (!(point.upper_int > point.bound_lower)) {
      continue;
    }
    point.needs_integral = true;
    active_integral[i] = 1u;
  }

  std::vector<double> cdf_query_times;
  std::vector<const uuber::TrialParamsSoA *> cdf_query_params;
  std::vector<std::size_t> cdf_query_point_idx;
  std::vector<std::uint8_t> cdf_query_slot;
  cdf_query_times.reserve(point_count * 2u);
  cdf_query_params.reserve(point_count * 2u);
  cdf_query_point_idx.reserve(point_count * 2u);
  cdf_query_slot.reserve(point_count * 2u);

  for (std::size_t i = 0; i < point_count; ++i) {
    const ChainedPointMeta &point = meta[i];
    if (point.has_conditioning_bounds) {
      if (std::isfinite(point.bound_lower) && point.bound_lower > 0.0) {
        cdf_query_times.push_back(point.bound_lower);
        cdf_query_params.push_back(point.trial_params_soa);
        cdf_query_point_idx.push_back(i);
        cdf_query_slot.push_back(0u);
      }
      cdf_query_times.push_back(point.bound_upper);
      cdf_query_params.push_back(point.trial_params_soa);
      cdf_query_point_idx.push_back(i);
      cdf_query_slot.push_back(1u);
    } else if (point.needs_source_finish) {
      cdf_query_times.push_back(std::numeric_limits<double>::infinity());
      cdf_query_params.push_back(point.trial_params_soa);
      cdf_query_point_idx.push_back(i);
      cdf_query_slot.push_back(2u);
    }
  }

  if (!cdf_query_times.empty()) {
    uuber::TreeNodeBatchValues source_cdf_values;
    if (!eval_event_ref_batch_with_params(
            state, source_ref, 0u, cdf_query_times, cdf_query_params,
            EvalNeed::kCDF, source_cdf_values) ||
        source_cdf_values.cdf.size() != cdf_query_times.size()) {
      return false;
    }
    for (std::size_t query_idx = 0; query_idx < cdf_query_times.size();
         ++query_idx) {
      const std::size_t point_idx = cdf_query_point_idx[query_idx];
      const double value = clamp_probability(source_cdf_values.cdf[query_idx]);
      switch (cdf_query_slot[query_idx]) {
      case 0u:
        source_lower_cdf[point_idx] = value;
        break;
      case 1u:
        source_upper_cdf[point_idx] = value;
        break;
      default:
        source_finish_cdf[point_idx] = value;
        break;
      }
    }
  }

  for (std::size_t i = 0; i < point_count; ++i) {
    ChainedPointMeta &point = meta[i];
    if (point.has_conditioning_bounds) {
      point.source_mass = source_upper_cdf[i] - source_lower_cdf[i];
      if (!std::isfinite(point.source_mass) || point.source_mass <= 0.0) {
        point.needs_integral = false;
        active_integral[i] = 0u;
        continue;
      }
    }
    if (point.needs_source_finish) {
      const double cdf_total =
          clamp_probability((1.0 - point.q) * source_finish_cdf[i]);
      if (needs_cdf(need)) {
        out_values.cdf[i] = cdf_total;
      }
      if (needs_survival(need)) {
        out_values.survival[i] = clamp_probability(1.0 - cdf_total);
      }
    }
  }

  std::vector<double> upper_bounds(point_count, 0.0);
  bool any_integral = false;
  for (std::size_t i = 0; i < point_count; ++i) {
    upper_bounds[i] = meta[i].upper_int;
    any_integral = any_integral || active_integral[i] != 0u;
  }
  if (!any_integral) {
    return true;
  }

  std::vector<const uuber::TrialParamsSoA *> point_trial_params_soa_storage;
  const std::vector<const uuber::TrialParamsSoA *> *point_trial_params_soa =
      nullptr;
  if (!resolve_point_trial_params_soa_batch(
          state.ctx, point_count, state.trial_params, state.trial_params_soa_batch,
          point_trial_params_soa_storage, point_trial_params_soa)) {
    return false;
  }

  std::vector<TimeIntegralPointPlan> point_plans;
  std::vector<double> query_times;
  std::vector<double> query_weights;
  std::vector<const uuber::TrialParamsSoA *> query_trial_params_soa;
  const double integration_lower =
      source_has_bounds ? source_bound_lower : 0.0;
  if (!build_time_integral_query_plan(
          upper_bounds, integration_lower,
          max_positive_finite_upper(upper_bounds), 4, point_plans, query_times,
          query_weights, point_trial_params_soa, &query_trial_params_soa,
          &active_integral)) {
    return false;
  }

  if (query_times.empty()) {
    return true;
  }

  uuber::TreeNodeBatchValues source_density_values;
  if (!eval_event_ref_batch_with_params(
          state, source_ref, 0u, query_times, query_trial_params_soa,
          EvalNeed::kDensity, source_density_values) ||
      source_density_values.density.size() != query_times.size()) {
    return false;
  }

  std::vector<double> shifted;
  std::vector<double> kernel_density;
  std::vector<double> kernel_cdf;
  for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
    if (active_integral[point_idx] == 0u) {
      continue;
    }
    const ChainedPointMeta &point = meta[point_idx];
    const TimeIntegralPointPlan &plan = point_plans[point_idx];
    if (plan.count == 0u) {
      continue;
    }

    shifted.resize(plan.count);
    for (std::size_t offset = 0; offset < plan.count; ++offset) {
      const std::size_t query_idx = plan.begin + offset;
      shifted[offset] = point.t - query_times[query_idx] - point.x_shift;
    }

    if (needs_density(need)) {
      kernel_density.assign(plan.count, 0.0);
      eval_pdf_vec_with_lower_bound(point.cfg, point.lower_bound,
                                    shifted.data(), plan.count,
                                    kernel_density.data());
    }
    if (needs_cdf(need) || needs_survival(need)) {
      kernel_cdf.assign(plan.count, 0.0);
      eval_cdf_vec_with_lower_bound(point.cfg, point.lower_bound,
                                    shifted.data(), plan.count,
                                    kernel_cdf.data());
    }

    double density_total = 0.0;
    double cdf_total_cond = 0.0;
    for (std::size_t offset = 0; offset < plan.count; ++offset) {
      const std::size_t query_idx = plan.begin + offset;
      const double weight = query_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      double source_density =
          safe_density(source_density_values.density[query_idx]);
      if (!(source_density > 0.0)) {
        continue;
      }
      if (point.has_conditioning_bounds) {
        source_density /= point.source_mass;
      }
      if (needs_density(need)) {
        const double fx = safe_density(kernel_density[offset]);
        if (fx > 0.0) {
          density_total += weight * source_density * fx;
        }
      }
      if (needs_cdf(need) || needs_survival(need)) {
        const double Fx = clamp_probability(kernel_cdf[offset]);
        if (Fx > 0.0) {
          cdf_total_cond += weight * source_density * Fx;
        }
      }
    }

    if (needs_density(need)) {
      out_values.density[point_idx] =
          safe_density((1.0 - point.q) * density_total);
    }
    if (needs_cdf(need) || needs_survival(need)) {
      const double cdf_total =
          clamp_probability((1.0 - point.q) *
                            clamp_probability(cdf_total_cond));
      if (needs_cdf(need)) {
        out_values.cdf[point_idx] = cdf_total;
      }
      if (needs_survival(need)) {
        out_values.survival[point_idx] = clamp_probability(1.0 - cdf_total);
      }
    }
  }

  return true;
}

inline bool eval_event_ref_batch(const NodeEvalState &state,
                                 const LabelRef &label_ref,
                                 std::uint32_t node_flags,
                                 const std::vector<double> &times,
                                 EvalNeed need,
                                 uuber::TreeNodeBatchValues &out_values) {
  out_values = uuber::TreeNodeBatchValues{};
  const std::size_t point_count = times.size();
  out_values.density.assign(point_count, 0.0);
  out_values.survival.assign(point_count, 1.0);
  out_values.cdf.assign(point_count, 0.0);
  if (point_count == 0u) {
    return true;
  }

  const int event_idx = find_event_idx_for_label_ref(state.ctx, label_ref);
  if (event_idx >= 0 &&
      event_idx < static_cast<int>(state.ctx.ir.events.size())) {
    const uuber::IrEvent &event =
        state.ctx.ir.events[static_cast<std::size_t>(event_idx)];
    if (try_eval_conditioned_simple_acc_event_batch(
            state, event, times, need, out_values)) {
      return true;
    }
    if (node_flags == 0u && event.node_idx >= 0 &&
        event.node_idx < static_cast<int>(state.ctx.ir.nodes.size())) {
      node_flags =
          state.ctx.ir.nodes[static_cast<std::size_t>(event.node_idx)].flags;
    }
  }

  const bool is_special_deadline =
      (node_flags & uuber::IR_NODE_FLAG_SPECIAL_DEADLINE) != 0u;
  const bool is_special_guess =
      (node_flags & uuber::IR_NODE_FLAG_SPECIAL_GUESS) != 0u;
  if (is_special_deadline) {
    for (std::size_t i = 0; i < point_count; ++i) {
      if (needs_survival(need) || needs_cdf(need)) {
        if (std::isfinite(times[i]) && times[i] < 0.0) {
          out_values.survival[i] = 1.0;
          out_values.cdf[i] = 0.0;
        } else {
          out_values.survival[i] = 0.0;
          out_values.cdf[i] = 1.0;
        }
      }
    }
    return true;
  }

  std::vector<double> donor_density(point_count, 0.0);
  if (needs_density(need)) {
    int target_outcome_idx = (state.outcome_idx >= 0) ? state.outcome_idx
                                                      : label_ref.outcome_idx;
    int target_label_id = label_ref.label_id;
    if (state.outcome_idx >= 0 &&
        state.outcome_idx <
            static_cast<int>(state.ctx.outcome_label_ids.size())) {
      target_label_id = state.ctx.outcome_label_ids
          [static_cast<std::size_t>(state.outcome_idx)];
    }
    const bool needs_component_guess_density =
        component_guess_density_applicable(
            state.ctx, target_label_id, target_outcome_idx, is_special_guess,
            state.component_idx);
    const bool needs_outcome_alias_density = outcome_alias_density_applicable(
        state.ctx, target_outcome_idx, state.include_na_donors);
    std::vector<double> density_scratch;
    if (needs_component_guess_density &&
        !accumulate_component_guess_density_batch_idx(
            state.ctx, target_label_id, target_outcome_idx, is_special_guess,
            times, state.component_idx, state.trial_params, state.trial_type_key,
            donor_density, nullptr, false, &density_scratch)) {
      return false;
    }
    if (needs_outcome_alias_density) {
      std::vector<double> alias_density;
      if (!accumulate_outcome_alias_density_batch_idx(
              state.ctx, target_outcome_idx, times, state.component_idx,
              state.trial_params, state.trial_type_key,
              state.include_na_donors, alias_density, nullptr, false,
              &density_scratch) ||
          alias_density.size() != point_count) {
        return false;
      }
      for (std::size_t i = 0; i < point_count; ++i) {
        donor_density[i] =
            safe_density(donor_density[i] + alias_density[i]);
      }
    }
  }

  if (is_special_guess) {
    if (needs_density(need)) {
      out_values.density = std::move(donor_density);
    }
    if (needs_survival(need)) {
      std::fill(out_values.survival.begin(), out_values.survival.end(), 0.0);
    }
    if (needs_cdf(need)) {
      std::fill(out_values.cdf.begin(), out_values.cdf.end(), 1.0);
    }
    return true;
  }

  int label_idx = label_ref.label_id;
  bool self_has_exact_source_time = false;
  double self_exact_source_time = std::numeric_limits<double>::quiet_NaN();
  bool self_has_source_bounds = false;
  double self_bound_lower = 0.0;
  double self_bound_upper = std::numeric_limits<double>::infinity();
  if (label_idx >= 0 && label_idx != NA_INTEGER) {
    (void)evaluator_resolve_label_time_constraint(
        &state.time_constraints, label_idx, self_has_exact_source_time,
        self_exact_source_time, self_has_source_bounds, self_bound_lower,
        self_bound_upper);
  }
  const bool self_is_time_conditioned =
      self_has_exact_source_time || self_has_source_bounds;
  if (label_idx >= 0 && label_idx != NA_INTEGER) {
    if (!self_is_time_conditioned &&
        forced_state_contains_complete(state.forced_state, label_idx)) {
      if (needs_cdf(need)) {
        std::fill(out_values.cdf.begin(), out_values.cdf.end(), 1.0);
      }
      if (needs_survival(need)) {
        std::fill(out_values.survival.begin(), out_values.survival.end(), 0.0);
      }
      if (needs_density(need)) {
        out_values.density = donor_density;
      }
      return true;
    }
    if (forced_state_contains_survive(state.forced_state, label_idx)) {
      if (self_is_time_conditioned) {
        if (needs_survival(need)) {
          std::fill(out_values.survival.begin(), out_values.survival.end(), 0.0);
        }
      } else {
        if (needs_survival(need)) {
          std::fill(out_values.survival.begin(), out_values.survival.end(), 1.0);
        }
      }
      if (needs_cdf(need)) {
        std::fill(out_values.cdf.begin(), out_values.cdf.end(), 0.0);
      }
      if (needs_density(need)) {
        out_values.density = donor_density;
      }
      return true;
    }
  }

  if (label_ref.acc_idx >= 0 &&
      label_ref.acc_idx < static_cast<int>(state.ctx.accumulators.size())) {
    const int acc_idx = label_ref.acc_idx;
    const uuber::NativeAccumulator &acc =
        state.ctx.accumulators[static_cast<std::size_t>(acc_idx)];
    const TrialAccumulatorParams *override =
        evaluator_get_trial_param_entry(state.trial_params, acc_idx);
    if (!evaluator_component_active_idx(acc, state.component_idx, override)) {
      if (needs_density(need)) {
        out_values.density = donor_density;
      }
      return true;
    }

    uuber::TreeNodeBatchValues base_values;
    EvalNeed base_need = static_cast<EvalNeed>(0u);
    if (needs_density(need)) {
      base_need = base_need | EvalNeed::kDensity;
    }
    if (needs_survival(need)) {
      base_need = base_need | EvalNeed::kSurvival;
    }
    if (needs_cdf(need)) {
      base_need = base_need | EvalNeed::kCDF;
    }
    if (self_is_time_conditioned && !self_has_exact_source_time) {
      base_need = base_need | EvalNeed::kCDF;
      if (needs_density(need)) {
        base_need = base_need | EvalNeed::kDensity;
      }
    }
    if (!eval_accumulator_base_batch(state, label_ref, node_flags, times,
                                     base_need, base_values)) {
      return false;
    }

    if (self_is_time_conditioned) {
      if (self_has_exact_source_time) {
        if (self_has_source_bounds &&
            (!(self_exact_source_time > self_bound_lower) ||
             self_exact_source_time > self_bound_upper)) {
          if (needs_density(need)) {
            out_values.density = donor_density;
          }
          if (needs_survival(need)) {
            std::fill(out_values.survival.begin(), out_values.survival.end(), 0.0);
          }
          if (needs_cdf(need)) {
            std::fill(out_values.cdf.begin(), out_values.cdf.end(), 0.0);
          }
          return true;
        }
        for (std::size_t i = 0; i < point_count; ++i) {
          if (needs_density(need)) {
            out_values.density[i] = donor_density[i];
          }
          const bool before_exact =
              std::isfinite(times[i]) && times[i] < self_exact_source_time;
          const double cdf_exact = before_exact ? 0.0 : 1.0;
          if (needs_cdf(need)) {
            out_values.cdf[i] = cdf_exact;
          }
          if (needs_survival(need)) {
            out_values.survival[i] = 1.0 - cdf_exact;
          }
        }
        return true;
      }

      std::vector<double> bound_query_times;
      std::vector<const uuber::TrialParamsSoA *> bound_query_params;
      std::vector<std::size_t> bound_query_point_idx;
      std::vector<std::uint8_t> bound_query_slot;
      bound_query_times.reserve(point_count * 2u);
      bound_query_params.reserve(point_count * 2u);
      bound_query_point_idx.reserve(point_count * 2u);
      bound_query_slot.reserve(point_count * 2u);
      for (std::size_t i = 0; i < point_count; ++i) {
        const uuber::TrialParamsSoA *point_soa =
            trial_params_soa_for_batch_point(state, i);
        if (std::isfinite(self_bound_lower)) {
          bound_query_times.push_back(self_bound_lower);
          bound_query_params.push_back(point_soa);
          bound_query_point_idx.push_back(i);
          bound_query_slot.push_back(0u);
        }
        bound_query_times.push_back(self_bound_upper);
        bound_query_params.push_back(point_soa);
        bound_query_point_idx.push_back(i);
        bound_query_slot.push_back(1u);
      }

      std::vector<double> lower_cdf(point_count, 0.0);
      std::vector<double> upper_cdf(point_count, 1.0);
      if (!bound_query_times.empty()) {
        uuber::TreeNodeBatchValues bound_values;
        if (!eval_accumulator_base_batch_with_params(
                state, label_ref, node_flags, bound_query_times,
                bound_query_params, EvalNeed::kCDF, bound_values) ||
            bound_values.cdf.size() != bound_query_times.size()) {
          return false;
        }
        for (std::size_t query_idx = 0; query_idx < bound_query_times.size();
             ++query_idx) {
          const std::size_t point_idx = bound_query_point_idx[query_idx];
          const double value = clamp_probability(bound_values.cdf[query_idx]);
          if (bound_query_slot[query_idx] == 0u) {
            lower_cdf[point_idx] = value;
          } else {
            upper_cdf[point_idx] = value;
          }
        }
      }

      for (std::size_t i = 0; i < point_count; ++i) {
        const double condition_mass = upper_cdf[i] - lower_cdf[i];
        if (!(self_bound_upper > self_bound_lower) ||
            !std::isfinite(condition_mass) || condition_mass <= 0.0) {
          if (needs_density(need)) {
            out_values.density[i] = donor_density[i];
          }
          if (needs_survival(need)) {
            out_values.survival[i] = 0.0;
          }
          if (needs_cdf(need)) {
            out_values.cdf[i] = 0.0;
          }
          continue;
        }

        if (needs_density(need)) {
          if (std::isfinite(times[i]) && times[i] > self_bound_lower &&
              times[i] <= self_bound_upper) {
            out_values.density[i] = safe_density(
                safe_density(base_values.density[i]) / condition_mass);
          } else {
            out_values.density[i] = 0.0;
          }
          out_values.density[i] =
              safe_density(out_values.density[i] + donor_density[i]);
        }

        double cdf_cond = 0.0;
        if (!std::isfinite(times[i]) || times[i] >= self_bound_upper) {
          cdf_cond = 1.0;
        } else if (!(times[i] > self_bound_lower)) {
          cdf_cond = 0.0;
        } else {
          cdf_cond = clamp_probability(
              (clamp_probability(base_values.cdf[i]) - lower_cdf[i]) /
              condition_mass);
        }
        if (needs_cdf(need)) {
          out_values.cdf[i] = cdf_cond;
        }
        if (needs_survival(need)) {
          out_values.survival[i] = clamp_probability(1.0 - cdf_cond);
        }
      }
      return true;
    }

    for (std::size_t i = 0; i < point_count; ++i) {
      if (needs_density(need)) {
        out_values.density[i] =
            safe_density(base_values.density[i] + donor_density[i]);
      }
      if (needs_survival(need)) {
        out_values.survival[i] = clamp_probability(base_values.survival[i]);
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] = clamp_probability(base_values.cdf[i]);
      }
    }
    return true;
  }

  if (label_ref.pool_idx >= 0 &&
      label_ref.pool_idx < static_cast<int>(state.ctx.pools.size())) {
    const uuber::NativePool &pool =
        state.ctx.pools[static_cast<std::size_t>(label_ref.pool_idx)];
    std::vector<std::size_t> active_members;
    active_members.reserve(pool.member_refs.size());
    for (std::size_t member_idx = 0; member_idx < pool.member_refs.size();
         ++member_idx) {
      const LabelRef &member_ref = pool.member_refs[member_idx];
      if (member_ref.acc_idx >= 0 &&
          member_ref.acc_idx <
              static_cast<int>(state.ctx.accumulators.size())) {
        const uuber::NativeAccumulator &member_acc =
            state.ctx.accumulators[static_cast<std::size_t>(member_ref.acc_idx)];
        const TrialAccumulatorParams *override =
            evaluator_get_trial_param_entry(state.trial_params,
                                            member_ref.acc_idx);
        if (!evaluator_component_active_idx(member_acc, state.component_idx,
                                            override)) {
          continue;
        }
      }
      active_members.push_back(member_idx);
    }
    if (active_members.empty()) {
      if (needs_density(need)) {
        out_values.density = donor_density;
      }
      return true;
    }

    EvalNeed child_need = needs_density(need)
                              ? (EvalNeed::kDensity | EvalNeed::kSurvival)
                              : EvalNeed::kSurvival;
    std::vector<std::vector<double>> member_density;
    std::vector<std::vector<double>> member_survival;
    if (needs_density(need)) {
      member_density.resize(active_members.size());
    }
    member_survival.resize(active_members.size());

    for (std::size_t pos = 0; pos < active_members.size(); ++pos) {
      const LabelRef &member_ref = pool.member_refs[active_members[pos]];
      uuber::TreeNodeBatchValues child_values;
      if (!eval_event_ref_batch(state, member_ref, 0u, times, child_need,
                                child_values) ||
          child_values.survival.size() != point_count ||
          (needs_density(need) && child_values.density.size() != point_count)) {
        return false;
      }
      member_survival[pos] = std::move(child_values.survival);
      if (needs_density(need)) {
        member_density[pos] = std::move(child_values.density);
      }
    }

    std::vector<double> density_scratch(active_members.size(), 0.0);
    std::vector<double> survival_scratch(active_members.size(), 1.0);
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      for (std::size_t member_pos = 0; member_pos < active_members.size();
           ++member_pos) {
        survival_scratch[member_pos] =
            clamp_probability(member_survival[member_pos][point_idx]);
        if (needs_density(need)) {
          density_scratch[member_pos] =
              safe_density(member_density[member_pos][point_idx]);
        }
      }
      if (needs_density(need)) {
        out_values.density[point_idx] =
            pool_density_fast(density_scratch, survival_scratch, pool.k);
        if (!std::isfinite(out_values.density[point_idx]) ||
            out_values.density[point_idx] < 0.0) {
          out_values.density[point_idx] = 0.0;
        }
        out_values.density[point_idx] =
            safe_density(out_values.density[point_idx] + donor_density[point_idx]);
      }
      if (needs_survival(need) || needs_cdf(need)) {
        out_values.survival[point_idx] =
            pool_survival_fast(survival_scratch, pool.k);
        if (!std::isfinite(out_values.survival[point_idx])) {
          out_values.survival[point_idx] = 0.0;
        }
      }
      if (needs_cdf(need)) {
        out_values.cdf[point_idx] =
            clamp_probability(1.0 - out_values.survival[point_idx]);
      }
    }
    return true;
  }

  if (needs_density(need)) {
    out_values.density = donor_density;
  }
  return true;
}

inline bool evaluate_specialized_pair_overlap_density_batch(
    const uuber::NativeContext &ctx, const uuber::IrOutcomeCouplingOp &spec,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints = nullptr,
    const ForcedStateView *forced_state_view = nullptr) {
  density_out.assign(times.size(), 0.0);
  if (spec.target_event_idx < 0 || spec.gate_event_idx < 0 ||
      spec.aux_event_count != 1 || spec.aux_event_begin < 0 ||
      spec.aux_event_begin >=
          static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
    return false;
  }
  const int y_event_idx = ctx.ir.outcome_coupling_event_indices
      [static_cast<std::size_t>(spec.aux_event_begin)];
  if (y_event_idx < 0 || y_event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return false;
  }

  std::vector<double> eval_times(times);
  std::vector<std::uint8_t> valid_points(times.size(), 0u);
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (std::isfinite(times[i]) && times[i] >= 0.0) {
      valid_points[i] = 1u;
    } else {
      eval_times[i] = 0.0;
    }
  }

  uuber::TreeNodeBatchValues gate_values;
  uuber::TreeNodeBatchValues target_values;
  uuber::TreeNodeBatchValues other_values;
  if (!evaluate_overlap_event_batch_values(
          ctx, spec.gate_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          gate_values, time_constraints, forced_state_view) ||
      !evaluate_overlap_event_batch_values(
          ctx, spec.target_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          target_values, time_constraints, forced_state_view) ||
      !evaluate_overlap_event_batch_values(
          ctx, y_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          other_values, time_constraints, forced_state_view)) {
    return false;
  }

  std::vector<const uuber::TrialParamsSoA *> point_trial_params_soa_storage;
  const std::vector<const uuber::TrialParamsSoA *> *point_trial_params_soa =
      nullptr;
  if (!resolve_point_trial_params_soa_batch(
          ctx, times.size(), trial_params, trial_params_soa_batch,
          point_trial_params_soa_storage, point_trial_params_soa)) {
    return false;
  }

  std::vector<TimeIntegralPointPlan> quadrature_plans;
  std::vector<double> quadrature_times;
  std::vector<double> quadrature_weights;
  std::vector<const uuber::TrialParamsSoA *> quadrature_trial_params_soa;
  quadrature_plans.assign(times.size(), TimeIntegralPointPlan{});
  const double quadrature_upper = max_positive_finite_upper(eval_times);
  if (quadrature_upper > 0.0 &&
      !build_time_integral_query_plan(
          eval_times, 0.0, quadrature_upper, 4, quadrature_plans,
          quadrature_times, quadrature_weights, point_trial_params_soa,
          &quadrature_trial_params_soa, &valid_points)) {
    return false;
  }

  uuber::TreeNodeBatchValues quadrature_target_values;
  uuber::TreeNodeBatchValues quadrature_other_values;
  if (!quadrature_times.empty() &&
      (!evaluate_overlap_event_batch_values(
           ctx, spec.target_event_idx, quadrature_times, component_idx,
           trial_params, trial_type_key, &quadrature_trial_params_soa,
           quadrature_target_values, time_constraints, forced_state_view) ||
       !evaluate_overlap_event_batch_values(
           ctx, y_event_idx, quadrature_times, component_idx, trial_params,
           trial_type_key, &quadrature_trial_params_soa,
           quadrature_other_values, time_constraints, forced_state_view))) {
    return false;
  }

  for (std::size_t i = 0; i < times.size(); ++i) {
    if (valid_points[i] == 0u) {
      continue;
    }
    double int_fx_fy = 0.0;
    const TimeIntegralPointPlan &plan = quadrature_plans[i];
    for (std::size_t query_idx = plan.begin;
         query_idx < plan.begin + plan.count; ++query_idx) {
      const double weight = quadrature_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      const double term = safe_density(
          quadrature_target_values.density[query_idx] *
          quadrature_other_values.cdf[query_idx]);
      if (term > 0.0) {
        int_fx_fy += weight * term;
      }
    }
    const double density = safe_density(
        gate_values.density[i] * target_values.cdf[i] +
        target_values.density[i] * gate_values.cdf[i] -
        gate_values.density[i] * int_fx_fy -
        gate_values.cdf[i] * safe_density(target_values.density[i] *
                                          other_values.cdf[i]));
    density_out[i] = density;
  }
  return true;
}

inline bool evaluate_specialized_guarded_pair_overlap_density_batch(
    const uuber::NativeContext &ctx, const uuber::IrOutcomeCouplingOp &spec,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints = nullptr,
    const ForcedStateView *forced_state_view = nullptr) {
  density_out.assign(times.size(), 0.0);
  if (spec.target_event_idx < 0 || spec.gate_event_idx < 0 ||
      spec.blocker_event_idx < 0 || spec.aux_event_count != 1 ||
      spec.aux_event_begin < 0 ||
      spec.aux_event_begin >=
          static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
    return false;
  }
  const int other_event_idx = ctx.ir.outcome_coupling_event_indices
      [static_cast<std::size_t>(spec.aux_event_begin)];
  if (other_event_idx < 0 ||
      other_event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return false;
  }

  std::vector<double> eval_times(times);
  std::vector<std::uint8_t> valid_points(times.size(), 0u);
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (std::isfinite(times[i]) && times[i] >= 0.0) {
      valid_points[i] = 1u;
    } else {
      eval_times[i] = 0.0;
    }
  }

  uuber::TreeNodeBatchValues gate_values;
  uuber::TreeNodeBatchValues target_values;
  uuber::TreeNodeBatchValues other_values;
  uuber::TreeNodeBatchValues blocker_values;
  if (!evaluate_overlap_event_batch_values(
          ctx, spec.gate_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          gate_values, time_constraints, forced_state_view) ||
      !evaluate_overlap_event_batch_values(
          ctx, spec.target_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          target_values, time_constraints, forced_state_view) ||
      !evaluate_overlap_event_batch_values(
          ctx, other_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          other_values, time_constraints, forced_state_view) ||
      !evaluate_overlap_event_batch_values(
          ctx, spec.blocker_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          blocker_values, time_constraints, forced_state_view)) {
    return false;
  }

  std::vector<const uuber::TrialParamsSoA *> point_trial_params_soa_storage;
  const std::vector<const uuber::TrialParamsSoA *> *point_trial_params_soa =
      nullptr;
  if (!resolve_point_trial_params_soa_batch(
          ctx, times.size(), trial_params, trial_params_soa_batch,
          point_trial_params_soa_storage, point_trial_params_soa)) {
    return false;
  }

  std::vector<TimeIntegralPointPlan> quadrature_plans;
  std::vector<double> quadrature_times;
  std::vector<double> quadrature_weights;
  std::vector<const uuber::TrialParamsSoA *> quadrature_trial_params_soa;
  quadrature_plans.assign(times.size(), TimeIntegralPointPlan{});
  const double quadrature_upper = max_positive_finite_upper(eval_times);
  if (quadrature_upper > 0.0 &&
      !build_time_integral_query_plan(
          eval_times, 0.0, quadrature_upper, 2, quadrature_plans,
          quadrature_times, quadrature_weights, point_trial_params_soa,
          &quadrature_trial_params_soa, &valid_points)) {
    return false;
  }

  uuber::TreeNodeBatchValues quadrature_gate_values;
  uuber::TreeNodeBatchValues quadrature_target_values;
  uuber::TreeNodeBatchValues quadrature_other_values;
  uuber::TreeNodeBatchValues quadrature_blocker_values;
  if (!quadrature_times.empty() &&
      (!evaluate_overlap_event_batch_values(
           ctx, spec.gate_event_idx, quadrature_times, component_idx,
           trial_params, trial_type_key, &quadrature_trial_params_soa,
           quadrature_gate_values, time_constraints, forced_state_view) ||
       !evaluate_overlap_event_batch_values(
           ctx, spec.target_event_idx, quadrature_times, component_idx,
           trial_params, trial_type_key, &quadrature_trial_params_soa,
           quadrature_target_values, time_constraints, forced_state_view) ||
       !evaluate_overlap_event_batch_values(
           ctx, other_event_idx, quadrature_times, component_idx, trial_params,
           trial_type_key, &quadrature_trial_params_soa,
           quadrature_other_values, time_constraints, forced_state_view) ||
       !evaluate_overlap_event_batch_values(
           ctx, spec.blocker_event_idx, quadrature_times, component_idx,
           trial_params, trial_type_key, &quadrature_trial_params_soa,
           quadrature_blocker_values, time_constraints, forced_state_view))) {
    return false;
  }

  for (std::size_t i = 0; i < times.size(); ++i) {
    if (valid_points[i] == 0u) {
      continue;
    }
    double target_before_guarded_mass = 0.0;
    double blocked_gate_order_mass = 0.0;
    double blocked_gate_mass = 0.0;
    double blocked_branch_mass = 0.0;
    const TimeIntegralPointPlan &plan = quadrature_plans[i];
    for (std::size_t query_idx = plan.begin;
         query_idx < plan.begin + plan.count; ++query_idx) {
      const double weight = quadrature_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      const double target_density = quadrature_target_values.density[query_idx];
      const double other_survival =
          quadrature_other_values.survival[query_idx];
      const double other_cdf = quadrature_other_values.cdf[query_idx];
      const double gate_density = quadrature_gate_values.density[query_idx];
      const double gate_cdf = quadrature_gate_values.cdf[query_idx];
      const double blocker_cdf = quadrature_blocker_values.cdf[query_idx];

      const double before_guarded =
          safe_density(target_density * other_survival);
      if (before_guarded > 0.0) {
        target_before_guarded_mass += weight * before_guarded;
      }
      if (!spec.target_branch_guarded) {
        const double gate_order = safe_density(target_density * other_cdf);
        if (gate_order > 0.0) {
          blocked_gate_order_mass += weight * gate_order;
        }
        const double gate_blocked =
            safe_density(gate_density * blocker_cdf * other_cdf);
        if (gate_blocked > 0.0) {
          blocked_gate_mass += weight * gate_blocked;
        }
        const double branch_blocked =
            safe_density(quadrature_other_values.density[query_idx] *
                         blocker_cdf * gate_cdf);
        if (branch_blocked > 0.0) {
          blocked_branch_mass += weight * branch_blocked;
        }
      }
    }

    double density = 0.0;
    if (spec.target_branch_guarded) {
      density = safe_density(
          gate_values.density[i] * blocker_values.survival[i] *
              target_before_guarded_mass +
          gate_values.cdf[i] *
              safe_density(target_values.density[i] * other_values.survival[i] *
                           blocker_values.survival[i]));
    } else {
      density = safe_density(
          gate_values.density[i] * target_before_guarded_mass +
          gate_values.density[i] * blocker_values.cdf[i] *
              blocked_gate_order_mass +
          target_values.density[i] * gate_values.cdf[i] *
              other_values.survival[i] +
          target_values.density[i] * blocked_gate_mass +
          target_values.density[i] * blocked_branch_mass);
    }
    density_out[i] = density;
  }
  return true;
}

inline bool evaluate_specialized_shared_blocker_pair_density_batch(
    const uuber::NativeContext &ctx, const uuber::IrOutcomeCouplingOp &spec,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints = nullptr,
    const ForcedStateView *forced_state_view = nullptr) {
  density_out.assign(times.size(), 0.0);
  if (spec.target_event_idx < 0 || spec.blocker_event_idx < 0 ||
      spec.aux_event_count != 1 || spec.aux_event_begin < 0 ||
      spec.aux_event_begin >=
          static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
    return false;
  }
  const int other_event_idx = ctx.ir.outcome_coupling_event_indices
      [static_cast<std::size_t>(spec.aux_event_begin)];
  if (other_event_idx < 0 ||
      other_event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return false;
  }

  std::vector<double> eval_times(times);
  std::vector<std::uint8_t> valid_points(times.size(), 0u);
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (std::isfinite(times[i]) && times[i] >= 0.0) {
      valid_points[i] = 1u;
    } else {
      eval_times[i] = 0.0;
    }
  }

  uuber::TreeNodeBatchValues target_values;
  uuber::TreeNodeBatchValues blocker_values;
  if (!evaluate_overlap_event_batch_values(
          ctx, spec.target_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch, target_values,
          time_constraints, forced_state_view) ||
      !evaluate_overlap_event_batch_values(
          ctx, spec.blocker_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch, blocker_values,
          time_constraints, forced_state_view)) {
    return false;
  }

  if (spec.target_branch_guarded) {
    for (std::size_t i = 0; i < times.size(); ++i) {
      if (valid_points[i] == 0u) {
        continue;
      }
      density_out[i] = safe_density(target_values.density[i] *
                                    blocker_values.survival[i]);
    }
    return true;
  }

  uuber::TreeNodeBatchValues guard_values;
  if (!evaluate_overlap_event_batch_values(
          ctx, other_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch, guard_values,
          time_constraints, forced_state_view)) {
    return false;
  }

  std::vector<const uuber::TrialParamsSoA *> point_trial_params_soa_storage;
  const std::vector<const uuber::TrialParamsSoA *> *point_trial_params_soa =
      nullptr;
  if (!resolve_point_trial_params_soa_batch(
          ctx, times.size(), trial_params, trial_params_soa_batch,
          point_trial_params_soa_storage, point_trial_params_soa)) {
    return false;
  }

  SimpleAccEventBatchInfo blocker_info;
  const uuber::IrEvent &blocker_event =
      ctx.ir.events[static_cast<std::size_t>(spec.blocker_event_idx)];
  if (!resolve_simple_acc_event_batch_info(
          ctx, blocker_event.node_idx, trial_params, point_trial_params_soa,
          component_idx, blocker_info, false)) {
    return false;
  }
  const double quadrature_lower =
      blocker_info.component_ok ? std::max(0.0, blocker_info.onset_eff) : 0.0;

  std::vector<TimeIntegralPointPlan> quadrature_plans;
  std::vector<double> quadrature_times;
  std::vector<double> quadrature_weights;
  std::vector<const uuber::TrialParamsSoA *> quadrature_trial_params_soa;
  quadrature_plans.assign(times.size(), TimeIntegralPointPlan{});
  const double quadrature_upper = max_positive_finite_upper(eval_times);
  if (quadrature_upper > quadrature_lower &&
      !build_time_integral_query_plan(
          eval_times, quadrature_lower, quadrature_upper, 2, quadrature_plans,
          quadrature_times, quadrature_weights, point_trial_params_soa,
          &quadrature_trial_params_soa, &valid_points)) {
    return false;
  }

  uuber::TreeNodeBatchValues quadrature_blocker_values;
  uuber::TreeNodeBatchValues quadrature_guard_values;
  if (!quadrature_times.empty() &&
      (!evaluate_overlap_event_batch_values(
           ctx, spec.blocker_event_idx, quadrature_times, component_idx,
           trial_params, trial_type_key, &quadrature_trial_params_soa,
           quadrature_blocker_values, time_constraints, forced_state_view) ||
       !evaluate_overlap_event_batch_values(
           ctx, other_event_idx, quadrature_times, component_idx, trial_params,
           trial_type_key, &quadrature_trial_params_soa,
           quadrature_guard_values, time_constraints, forced_state_view))) {
    return false;
  }

  for (std::size_t i = 0; i < times.size(); ++i) {
    if (valid_points[i] == 0u) {
      continue;
    }
    double blocked_mass = 0.0;
    const TimeIntegralPointPlan &plan = quadrature_plans[i];
    for (std::size_t query_idx = plan.begin;
         query_idx < plan.begin + plan.count; ++query_idx) {
      const double weight = quadrature_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      const double term = safe_density(
          quadrature_blocker_values.density[query_idx] *
          quadrature_guard_values.survival[query_idx]);
      if (term > 0.0) {
        blocked_mass += weight * term;
      }
    }
    density_out[i] = safe_density(
        blocker_values.density[i] * target_values.cdf[i] *
            guard_values.survival[i] +
        target_values.density[i] * blocked_mass);
  }
  return true;
}

inline bool evaluate_specialized_nway_overlap_density_batch(
    const uuber::NativeContext &ctx, const uuber::IrOutcomeCouplingOp &spec,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints = nullptr,
    const ForcedStateView *forced_state_view = nullptr) {
  density_out.assign(times.size(), 0.0);
  if (spec.target_event_idx < 0 || spec.gate_event_idx < 0 ||
      spec.aux_event_count <= 0 || spec.aux_event_begin < 0 ||
      spec.aux_event_begin + spec.aux_event_count >
          static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
    return false;
  }

  std::vector<int> competitor_event_indices;
  competitor_event_indices.reserve(static_cast<std::size_t>(spec.aux_event_count));
  for (int i = 0; i < spec.aux_event_count; ++i) {
    const int event_idx = ctx.ir.outcome_coupling_event_indices
        [static_cast<std::size_t>(spec.aux_event_begin + i)];
    if (event_idx < 0 || event_idx >= static_cast<int>(ctx.ir.events.size())) {
      return false;
    }
    competitor_event_indices.push_back(event_idx);
  }

  std::vector<double> eval_times(times);
  std::vector<std::uint8_t> valid_points(times.size(), 0u);
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (std::isfinite(times[i]) && times[i] >= 0.0) {
      valid_points[i] = 1u;
    } else {
      eval_times[i] = 0.0;
    }
  }

  uuber::TreeNodeBatchValues gate_values;
  uuber::TreeNodeBatchValues target_values;
  if (!evaluate_overlap_event_batch_values(
          ctx, spec.gate_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          gate_values, time_constraints, forced_state_view) ||
      !evaluate_overlap_event_batch_values(
          ctx, spec.target_event_idx, eval_times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch,
          target_values, time_constraints, forced_state_view)) {
    return false;
  }

  std::vector<uuber::TreeNodeBatchValues> competitor_values(
      competitor_event_indices.size());
  for (std::size_t comp_idx = 0; comp_idx < competitor_event_indices.size();
       ++comp_idx) {
    if (!evaluate_overlap_event_batch_values(
            ctx, competitor_event_indices[comp_idx], eval_times, component_idx,
            trial_params, trial_type_key, trial_params_soa_batch,
            competitor_values[comp_idx], time_constraints,
            forced_state_view)) {
      return false;
    }
  }

  std::vector<const uuber::TrialParamsSoA *> point_trial_params_soa_storage;
  const std::vector<const uuber::TrialParamsSoA *> *point_trial_params_soa =
      nullptr;
  if (!resolve_point_trial_params_soa_batch(
          ctx, times.size(), trial_params, trial_params_soa_batch,
          point_trial_params_soa_storage, point_trial_params_soa)) {
    return false;
  }

  std::vector<TimeIntegralPointPlan> quadrature_plans;
  std::vector<double> quadrature_times;
  std::vector<double> quadrature_weights;
  std::vector<const uuber::TrialParamsSoA *> quadrature_trial_params_soa;
  const double quadrature_upper = max_positive_finite_upper(eval_times);
  if (quadrature_upper > 0.0 &&
      !build_time_integral_query_plan(
          eval_times, 0.0, quadrature_upper, 2, quadrature_plans,
          quadrature_times, quadrature_weights, point_trial_params_soa,
          &quadrature_trial_params_soa, &valid_points)) {
    return false;
  }

  uuber::TreeNodeBatchValues quadrature_target_values;
  std::vector<uuber::TreeNodeBatchValues> quadrature_competitor_values(
      competitor_event_indices.size());
  if (!quadrature_times.empty() &&
      !evaluate_overlap_event_batch_values(
          ctx, spec.target_event_idx, quadrature_times, component_idx,
          trial_params, trial_type_key, &quadrature_trial_params_soa,
          quadrature_target_values, time_constraints, forced_state_view)) {
    return false;
  }
  for (std::size_t comp_idx = 0; comp_idx < competitor_event_indices.size();
       ++comp_idx) {
    if (!quadrature_times.empty() &&
        !evaluate_overlap_event_batch_values(
            ctx, competitor_event_indices[comp_idx], quadrature_times,
            component_idx, trial_params, trial_type_key,
            &quadrature_trial_params_soa,
            quadrature_competitor_values[comp_idx], time_constraints,
            forced_state_view)) {
      return false;
    }
  }

  for (std::size_t i = 0; i < times.size(); ++i) {
    if (valid_points[i] == 0u) {
      continue;
    }
    double current_order_density = target_values.density[i];
    for (const uuber::TreeNodeBatchValues &comp_values : competitor_values) {
      current_order_density = safe_density(current_order_density *
                                           comp_values.survival[i]);
    }
    double order_mass = 0.0;
    const TimeIntegralPointPlan &plan = quadrature_plans[i];
    for (std::size_t query_idx = plan.begin;
         query_idx < plan.begin + plan.count; ++query_idx) {
      const double weight = quadrature_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      double term = quadrature_target_values.density[query_idx];
      for (const uuber::TreeNodeBatchValues &comp_values :
           quadrature_competitor_values) {
        term = safe_density(term * comp_values.survival[query_idx]);
        if (term <= 0.0) {
          break;
        }
      }
      if (term > 0.0) {
        order_mass += weight * term;
      }
    }
    density_out[i] = safe_density(gate_values.density[i] * order_mass +
                                  gate_values.cdf[i] * current_order_density);
  }
  return true;
}

inline bool evaluate_specialized_overlap_density_batch(
    const uuber::NativeContext &ctx,
    const uuber::IrOutcomeCouplingOp *coupling_spec,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints = nullptr,
    const ForcedStateView *forced_state_view = nullptr) {
  if (coupling_spec == nullptr) {
    return false;
  }
  switch (coupling_spec->kind) {
  case uuber::IrOutcomeCouplingKind::Pair:
    return evaluate_specialized_pair_overlap_density_batch(
        ctx, *coupling_spec, times, component_idx, trial_params, trial_type_key,
        trial_params_soa_batch, density_out, time_constraints,
        forced_state_view);
  case uuber::IrOutcomeCouplingKind::GuardedPair:
    return evaluate_specialized_guarded_pair_overlap_density_batch(
        ctx, *coupling_spec, times, component_idx, trial_params, trial_type_key,
        trial_params_soa_batch, density_out, time_constraints,
        forced_state_view);
  case uuber::IrOutcomeCouplingKind::SharedBlockerPair:
    return evaluate_specialized_shared_blocker_pair_density_batch(
        ctx, *coupling_spec, times, component_idx, trial_params, trial_type_key,
        trial_params_soa_batch, density_out, time_constraints,
        forced_state_view);
  case uuber::IrOutcomeCouplingKind::NWay:
    return evaluate_specialized_nway_overlap_density_batch(
        ctx, *coupling_spec, times, component_idx, trial_params, trial_type_key,
        trial_params_soa_batch, density_out, time_constraints,
        forced_state_view);
  default:
    return false;
  }
}

struct TimeTrialParamsSoAKey {
  std::uint64_t time_bits{0ULL};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};

  bool operator==(const TimeTrialParamsSoAKey &other) const noexcept {
    return time_bits == other.time_bits &&
           trial_params_soa == other.trial_params_soa;
  }
};

struct TimeTrialParamsSoAKeyHash {
  std::size_t operator()(const TimeTrialParamsSoAKey &key) const noexcept {
    const std::size_t h1 = std::hash<std::uint64_t>{}(key.time_bits);
    const std::size_t h2 = std::hash<std::uintptr_t>{}(
        reinterpret_cast<std::uintptr_t>(key.trial_params_soa));
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

inline void compress_time_trial_params_soa_batch(
    const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    std::vector<double> &compressed_times,
    std::vector<const uuber::TrialParamsSoA *> &compressed_trial_params_soa,
    std::vector<std::size_t> &inverse_indices) {
  compressed_times.clear();
  compressed_trial_params_soa.clear();
  inverse_indices.assign(times.size(), 0u);
  if (times.empty()) {
    return;
  }
  std::unordered_map<TimeTrialParamsSoAKey, std::size_t,
                     TimeTrialParamsSoAKeyHash>
      index_by_key;
  index_by_key.reserve(times.size());
  for (std::size_t i = 0; i < times.size(); ++i) {
    const uuber::TrialParamsSoA *point_soa =
        i < trial_params_soa_batch.size() ? trial_params_soa_batch[i] : nullptr;
    const TimeTrialParamsSoAKey key{canonical_double_bits(times[i]), point_soa};
    auto it = index_by_key.find(key);
    if (it == index_by_key.end()) {
      const std::size_t next_idx = compressed_times.size();
      compressed_times.push_back(times[i]);
      compressed_trial_params_soa.push_back(point_soa);
      inverse_indices[i] = next_idx;
      index_by_key.emplace(key, next_idx);
    } else {
      inverse_indices[i] = it->second;
    }
  }
}

inline bool times_have_duplicates(const std::vector<double> &times) {
  if (times.size() < 2u) {
    return false;
  }
  std::unordered_set<std::uint64_t> seen;
  seen.reserve(times.size());
  for (double time_value : times) {
    const std::uint64_t bits = canonical_double_bits(time_value);
    if (!seen.insert(bits).second) {
      return true;
    }
  }
  return false;
}

struct RankedSequenceQueryKey {
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  const double *time_data{nullptr};
  std::size_t time_count{0u};

  bool operator==(const RankedSequenceQueryKey &other) const noexcept {
    if (trial_params_soa != other.trial_params_soa ||
        time_count != other.time_count) {
      return false;
    }
    if (time_data == other.time_data) {
      return true;
    }
    for (std::size_t i = 0; i < time_count; ++i) {
      if (canonical_double_bits(time_data[i]) !=
          canonical_double_bits(other.time_data[i])) {
        return false;
      }
    }
    return true;
  }
};

struct RankedSequenceQueryKeyHash {
  std::size_t operator()(const RankedSequenceQueryKey &key) const noexcept {
    std::uint64_t hash = kFNV64Offset;
    hash_append_u64(hash, static_cast<std::uint64_t>(key.time_count));
    hash_append_u64(
        hash, static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(
                  key.trial_params_soa)));
    for (std::size_t i = 0; i < key.time_count; ++i) {
      hash_append_u64(hash, canonical_double_bits(key.time_data[i]));
    }
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

inline void compress_ranked_sequence_query_batch(
    const std::vector<const double *> &times_by_query,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_by_query,
    std::size_t time_count, std::vector<const double *> &compressed_times_by_query,
    std::vector<const uuber::TrialParamsSoA *> &compressed_trial_params_soa,
    std::vector<std::size_t> &inverse_indices) {
  compressed_times_by_query.clear();
  compressed_trial_params_soa.clear();
  inverse_indices.assign(times_by_query.size(), 0u);
  if (times_by_query.empty()) {
    return;
  }
  std::unordered_map<RankedSequenceQueryKey, std::size_t,
                     RankedSequenceQueryKeyHash>
      index_by_key;
  index_by_key.reserve(times_by_query.size());
  for (std::size_t i = 0; i < times_by_query.size(); ++i) {
    const RankedSequenceQueryKey key{
        (i < trial_params_soa_by_query.size()) ? trial_params_soa_by_query[i]
                                               : nullptr,
        times_by_query[i], time_count};
    auto it = index_by_key.find(key);
    if (it == index_by_key.end()) {
      const std::size_t next_idx = compressed_times_by_query.size();
      compressed_times_by_query.push_back(times_by_query[i]);
      compressed_trial_params_soa.push_back(key.trial_params_soa);
      inverse_indices[i] = next_idx;
      index_by_key.emplace(key, next_idx);
    } else {
      inverse_indices[i] = it->second;
    }
  }
}

inline int resolve_outcome_index_ir(const uuber::NativeContext &ctx,
                                    int label_id, int component_idx);

inline int resolve_outcome_index_ir(const uuber::NativeContext &ctx,
                                    int label_id, int component_idx) {
  if (!ctx.ir.valid || label_id < 0)
    return -1;
  auto it = ctx.ir.label_id_to_outcomes.find(label_id);
  if (it == ctx.ir.label_id_to_outcomes.end() || it->second.empty()) {
    return -1;
  }
  const std::vector<int> &indices = it->second;
  if (indices.size() == 1 || component_idx < 0) {
    return indices[0];
  }
  int fallback_idx = -1;
  for (int idx : indices) {
    if (ir_outcome_allows_component(ctx, idx, component_idx)) {
      const uuber::OutcomeContextInfo &info =
          ctx.outcome_info[static_cast<std::size_t>(idx)];
      if (!info.allowed_components.empty())
        return idx;
      if (fallback_idx < 0)
        fallback_idx = idx;
    }
  }
  return fallback_idx;
}

struct ComponentCacheEntry;
double component_keep_weight(const uuber::NativeContext &ctx, int component_idx,
                             int outcome_idx);

double native_outcome_probability_impl_idx(
    const uuber::NativeContext &ctx, int node_id, double upper, int component_idx,
    SEXP forced_complete, SEXP forced_survive,
    const Rcpp::IntegerVector &competitor_ids, double rel_tol, double abs_tol,
    const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context);

struct DensityBatchWorkspace {
  std::vector<TimeIntegralPointPlan> point_plans;
  std::vector<double> query_times;
  std::vector<double> query_weights;
  std::vector<const uuber::TrialParamsSoA *> query_params_soa;
  std::vector<std::uint8_t> active_points;
  std::vector<std::uint8_t> tail_saw_positive;
  std::vector<int> tail_small_segments;
  std::vector<double> expanded_times;
  std::vector<const uuber::TrialParamsSoA *> expanded_params_soa;
  std::vector<double> compressed_times;
  std::vector<const uuber::TrialParamsSoA *> compressed_params_soa;
  std::vector<std::size_t> expanded_to_compressed;
  std::vector<double> batch_density;
  std::vector<double> segment_probabilities;
  std::vector<double> point_probabilities;
  std::vector<const uuber::TrialParamsSoA *> unique_params_soa;
  std::vector<const TrialParamSet *> unique_trial_params;
  std::vector<std::size_t> point_to_eval_index;
};

bool node_density_with_competitors_batch_internal(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<double> &times, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context, std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints = nullptr,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch =
        nullptr,
    DensityBatchWorkspace *workspace = nullptr);

bool node_density_entry_batch_idx(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<double> &times, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan, bool use_shared_trigger_eval,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints = nullptr,
    const std::vector<const uuber::TrialParamsSoA *>
        *trial_params_soa_batch_override = nullptr,
    const SharedTriggerMaskSoABatch *shared_trigger_mask_batch_override =
        nullptr,
    DensityBatchWorkspace *workspace = nullptr);

struct WeightedNodeDensityTerm {
  int node_id{-1};
  int component_idx{-1};
  int outcome_idx_context{-1};
  std::string trial_type_key;
  std::vector<int> competitor_ids;
  bool include_na_donors{false};
  double scaled_weight{0.0};

  bool same_key(const WeightedNodeDensityTerm &other) const {
    return node_id == other.node_id &&
           component_idx == other.component_idx &&
           outcome_idx_context == other.outcome_idx_context &&
           trial_type_key == other.trial_type_key &&
           competitor_ids == other.competitor_ids &&
           include_na_donors == other.include_na_donors;
  }

  bool same_eval_signature(const WeightedNodeDensityTerm &other) const {
    return component_idx == other.component_idx &&
           outcome_idx_context == other.outcome_idx_context &&
           trial_type_key == other.trial_type_key &&
           competitor_ids == other.competitor_ids &&
           include_na_donors == other.include_na_donors;
  }
};

inline void append_weighted_node_density_term(
    std::vector<WeightedNodeDensityTerm> &terms,
    int node_id, int component_idx, int outcome_idx_context,
    const std::string &trial_type_key, const std::vector<int> &competitor_ids,
    bool include_na_donors, double scaled_weight) {
  if (node_id < 0 || !std::isfinite(scaled_weight) || scaled_weight <= 0.0) {
    return;
  }
  WeightedNodeDensityTerm candidate;
  candidate.node_id = node_id;
  candidate.component_idx = component_idx;
  candidate.outcome_idx_context = outcome_idx_context;
  candidate.trial_type_key = trial_type_key;
  candidate.competitor_ids = competitor_ids;
  candidate.include_na_donors = include_na_donors;
  for (WeightedNodeDensityTerm &term : terms) {
    if (term.same_key(candidate)) {
      term.scaled_weight += scaled_weight;
      return;
    }
  }
  candidate.scaled_weight = scaled_weight;
  terms.push_back(std::move(candidate));
}

inline bool evaluate_node_density_multi_generic_batch_idx(
    const uuber::NativeContext &ctx, const std::vector<int> &node_ids,
    const std::vector<double> &times, int component_idx,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context, std::vector<double> &density_matrix_out,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch =
        nullptr) {
  density_matrix_out.assign(node_ids.size() * times.size(), 0.0);
  if (node_ids.empty() || times.empty()) {
    return true;
  }
  if (trial_params_soa_batch != nullptr &&
      trial_params_soa_batch->size() != times.size()) {
    return false;
  }

  std::vector<double> eval_times(times);
  std::vector<std::uint8_t> valid_points(times.size(), 0u);
  for (std::size_t i = 0; i < times.size(); ++i) {
    const double t = times[i];
    if (std::isfinite(t) && t >= 0.0) {
      valid_points[i] = 1u;
    } else {
      eval_times[i] = 0.0;
    }
  }

  NodeEvalState state(ctx, 0.0, component_idx, trial_params, trial_type_key,
                      include_na_donors, outcome_idx_context, nullptr, nullptr,
                      nullptr, false, nullptr, false);
  state.trial_params_soa_batch = trial_params_soa_batch;

  std::vector<int> node_indices;
  node_indices.reserve(node_ids.size());
  for (int node_id : node_ids) {
    const int node_idx = resolve_dense_node_idx_required(ctx, node_id);
    if (node_idx < 0) {
      return false;
    }
    const uuber::IrOutcomeCouplingOp *coupling_spec =
        find_outcome_coupling_spec(ctx, node_idx, competitor_ids);
    if (should_use_exact_density_program(ctx, node_idx, competitor_ids, state,
                                         coupling_spec)) {
      return false;
    }
    node_indices.push_back(node_idx);
  }

  std::vector<uuber::TreeNodeBatchValues> node_values;
  if (!evaluator_eval_nodes_batch_with_state_dense(node_indices, eval_times, state,
                                                   node_values) ||
      node_values.size() != node_ids.size()) {
    return false;
  }

  std::vector<double> survival_product(times.size(), 1.0);
  const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr;
  if (!competitor_ids.empty()) {
    competitor_cache = &fetch_competitor_cluster_cache(ctx, competitor_ids);
  }
  const bool has_competitors =
      competitor_cache != nullptr && !competitor_cache->compiled_ops.empty();
  if (has_competitors) {
    state.include_na_donors = false;
    state.outcome_idx = -1;
    competitor_survival_batch_from_state_compiled_ops(
        ctx, competitor_ids, state, eval_times, survival_product,
        competitor_cache);
    if (survival_product.size() != times.size()) {
      return false;
    }
  }

  for (std::size_t node_pos = 0; node_pos < node_values.size(); ++node_pos) {
    const uuber::TreeNodeBatchValues &values = node_values[node_pos];
    if (values.density.size() != times.size()) {
      return false;
    }
    double *density_row = density_matrix_out.data() +
                          node_pos * static_cast<std::ptrdiff_t>(times.size());
    for (std::size_t point_idx = 0; point_idx < times.size(); ++point_idx) {
      if (!valid_points[point_idx]) {
        density_row[point_idx] = 0.0;
        continue;
      }
      double density = safe_density(values.density[point_idx]);
      if (!(std::isfinite(density) && density > 0.0)) {
        density_row[point_idx] = 0.0;
        continue;
      }
      if (has_competitors) {
        const double surv = clamp_probability(survival_product[point_idx]);
        if (!(std::isfinite(surv) && surv > 0.0)) {
          density_row[point_idx] = 0.0;
          continue;
        }
        density *= surv;
      }
      density_row[point_idx] = safe_density(density);
    }
  }
  return true;
}

inline bool accumulate_weighted_node_density_terms_batch_idx(
    const uuber::NativeContext &ctx,
    const std::vector<WeightedNodeDensityTerm> &terms,
    const std::vector<double> &times, const TrialParamSet *trial_params,
    const SharedTriggerPlan *trigger_plan, bool use_shared_trigger_eval,
    std::vector<double> &density_out,
    std::vector<double> *density_scratch = nullptr) {
  density_out.assign(times.size(), 0.0);
  if (terms.empty()) {
    return true;
  }

  const bool allow_multi_node =
      !(use_shared_trigger_eval && trigger_plan != nullptr &&
        shared_trigger_count(*trigger_plan) > 0u);
  std::vector<double> density_local;
  std::vector<double> &density =
      density_scratch ? *density_scratch : density_local;
  std::vector<double> density_matrix;
  std::vector<std::uint8_t> consumed(terms.size(), 0u);

  for (std::size_t term_idx = 0; term_idx < terms.size(); ++term_idx) {
    if (consumed[term_idx] != 0u) {
      continue;
    }
    const WeightedNodeDensityTerm &term = terms[term_idx];
    if (!(std::isfinite(term.scaled_weight) && term.scaled_weight > 0.0)) {
      consumed[term_idx] = 1u;
      continue;
    }

    std::vector<std::size_t> group_indices;
    group_indices.reserve(terms.size() - term_idx);
    std::vector<int> group_node_ids;
    group_node_ids.reserve(terms.size() - term_idx);
    for (std::size_t other_idx = term_idx; other_idx < terms.size(); ++other_idx) {
      if (consumed[other_idx] != 0u ||
          !terms[other_idx].same_eval_signature(term)) {
        continue;
      }
      group_indices.push_back(other_idx);
      group_node_ids.push_back(terms[other_idx].node_id);
      consumed[other_idx] = 1u;
    }

    bool grouped = false;
    if (allow_multi_node && group_node_ids.size() > 1u) {
      grouped = evaluate_node_density_multi_generic_batch_idx(
          ctx, group_node_ids, times, term.component_idx, term.competitor_ids,
          trial_params, term.trial_type_key, term.include_na_donors,
          term.outcome_idx_context, density_matrix, nullptr);
    }

    if (grouped) {
      for (std::size_t group_pos = 0; group_pos < group_indices.size(); ++group_pos) {
        const WeightedNodeDensityTerm &group_term = terms[group_indices[group_pos]];
        const double *density_row = density_matrix.data() +
                                    group_pos *
                                        static_cast<std::ptrdiff_t>(times.size());
        for (std::size_t i = 0; i < density_out.size(); ++i) {
          const double value = density_row[i];
          if (std::isfinite(value) && value > 0.0) {
            density_out[i] = safe_density(
                density_out[i] + group_term.scaled_weight * value);
          }
        }
      }
      continue;
    }

    for (std::size_t group_pos = 0; group_pos < group_indices.size(); ++group_pos) {
      const WeightedNodeDensityTerm &group_term = terms[group_indices[group_pos]];
      density.clear();
      if (!node_density_entry_batch_idx(
              ctx, group_term.node_id, times, group_term.component_idx, nullptr,
              false, nullptr, false, group_term.competitor_ids, trial_params,
              group_term.trial_type_key, group_term.include_na_donors,
              group_term.outcome_idx_context, trigger_plan,
              use_shared_trigger_eval, density, nullptr) ||
          density.size() != density_out.size()) {
        return false;
      }
      for (std::size_t i = 0; i < density_out.size(); ++i) {
        const double value = density[i];
        if (std::isfinite(value) && value > 0.0) {
          density_out[i] = safe_density(
              density_out[i] + group_term.scaled_weight * value);
        }
      }
    }
  }
  return true;
}

inline bool component_guess_density_applicable(
    const uuber::NativeContext &ctx, int target_label_id,
    int target_outcome_idx, bool target_is_guess, int component_idx) {
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.component_info.size())) {
    return false;
  }
  const uuber::ComponentGuessPolicy &guess =
      ctx.component_info[static_cast<std::size_t>(component_idx)].guess;
  if (!guess.valid) {
    return false;
  }
  if (guess.target_is_guess) {
    if (!target_is_guess) {
      return false;
    }
  } else if (guess.target_outcome_idx >= 0) {
    if (target_outcome_idx != guess.target_outcome_idx) {
      return false;
    }
  } else if (guess.target_label_id != NA_INTEGER) {
    if (target_label_id == NA_INTEGER || target_label_id != guess.target_label_id) {
      return false;
    }
  } else if (!guess.target.empty()) {
    return false;
  }
  for (const auto &entry : guess.keep_weights) {
    if (entry.first >= 0 && (1.0 - entry.second) > 0.0) {
      return true;
    }
  }
  return false;
}

inline bool outcome_alias_density_applicable(
    const uuber::NativeContext &ctx, int target_outcome_idx,
    bool include_na_donors) {
  if (target_outcome_idx < 0 ||
      target_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return false;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(target_outcome_idx)];
  if (!info.alias_sources.empty()) {
    return true;
  }
  for (const auto &donor : info.guess_donors) {
    if (donor.rt_policy == "na" && !include_na_donors) {
      continue;
    }
    if (donor.outcome_idx >= 0 && donor.weight > 0.0) {
      return true;
    }
  }
  return false;
}

inline bool accumulate_component_guess_density_batch_idx(
    const uuber::NativeContext &ctx, int target_label_id,
    int target_outcome_idx, bool target_is_guess,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::vector<double> &density_out,
    const SharedTriggerPlan *trigger_plan = nullptr,
    bool use_shared_trigger_eval = false,
    std::vector<double> *donor_density_scratch = nullptr) {
  density_out.assign(times.size(), 0.0);
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.component_info.size())) {
    return true;
  }
  const uuber::ComponentGuessPolicy &guess =
      ctx.component_info[static_cast<std::size_t>(component_idx)].guess;
  if (!guess.valid) {
    return true;
  }
  if (guess.target_is_guess) {
    if (!target_is_guess) {
      return true;
    }
  } else if (guess.target_outcome_idx >= 0) {
    if (target_outcome_idx != guess.target_outcome_idx) {
      return true;
    }
  } else if (guess.target_label_id != NA_INTEGER) {
    if (target_label_id == NA_INTEGER || target_label_id != guess.target_label_id) {
      return true;
    }
  } else if (!guess.target.empty()) {
    return true;
  }
  std::vector<WeightedNodeDensityTerm> terms;
  terms.reserve(guess.keep_weights.size());
  for (const auto &entry : guess.keep_weights) {
    const int donor_outcome_idx = entry.first;
    const double release = 1.0 - entry.second;
    if (donor_outcome_idx < 0 || !(release > 0.0)) {
      continue;
    }
    if (donor_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
      continue;
    }
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(donor_outcome_idx)];
    if (info.node_id < 0) {
      continue;
    }
    std::vector<int> filtered_competitors;
    const std::vector<int> &competitors = filter_competitor_ids(
        ctx, info.competitor_ids, component_idx, filtered_competitors);
    append_weighted_node_density_term(terms, info.node_id, component_idx,
                                      donor_outcome_idx, trial_type_key,
                                      competitors, false, release);
  }
  return accumulate_weighted_node_density_terms_batch_idx(
      ctx, terms, times, trial_params, trigger_plan, use_shared_trigger_eval,
      density_out, donor_density_scratch);
}

inline bool accumulate_outcome_alias_density_batch_idx(
    const uuber::NativeContext &ctx, int target_outcome_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, std::vector<double> &density_out,
    const SharedTriggerPlan *trigger_plan = nullptr,
    bool use_shared_trigger_eval = false,
    std::vector<double> *donor_density_scratch = nullptr) {
  density_out.assign(times.size(), 0.0);
  if (target_outcome_idx < 0 ||
      target_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return true;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(target_outcome_idx)];
  std::vector<WeightedNodeDensityTerm> terms;
  terms.reserve(info.alias_sources.size() + info.guess_donors.size());
  for (int source_idx : info.alias_sources) {
    if (source_idx < 0 ||
        source_idx >= static_cast<int>(ctx.outcome_info.size())) {
      continue;
    }
    const uuber::OutcomeContextInfo &source_info =
        ctx.outcome_info[static_cast<std::size_t>(source_idx)];
    if (source_info.node_id < 0) {
      continue;
    }
    std::vector<int> filtered_competitors;
    const std::vector<int> &competitors = filter_competitor_ids(
        ctx, source_info.competitor_ids, component_idx, filtered_competitors);
    append_weighted_node_density_term(terms, source_info.node_id, component_idx,
                                      source_idx, trial_type_key, competitors,
                                      include_na_donors, 1.0);
  }
  for (const auto &donor : info.guess_donors) {
    if (donor.rt_policy == "na" && !include_na_donors) {
      continue;
    }
    if (donor.outcome_idx < 0 || !(donor.weight > 0.0)) {
      continue;
    }
    if (donor.outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
      continue;
    }
    const uuber::OutcomeContextInfo &donor_info =
        ctx.outcome_info[static_cast<std::size_t>(donor.outcome_idx)];
    if (donor_info.node_id < 0) {
      continue;
    }
    std::vector<int> filtered_competitors;
    const std::vector<int> &competitors = filter_competitor_ids(
        ctx, donor_info.competitor_ids, component_idx, filtered_competitors);
    append_weighted_node_density_term(terms, donor_info.node_id, component_idx,
                                      donor.outcome_idx, trial_type_key,
                                      competitors, include_na_donors,
                                      donor.weight);
  }
  return accumulate_weighted_node_density_terms_batch_idx(
      ctx, terms, times, trial_params, trigger_plan, use_shared_trigger_eval,
      density_out, donor_density_scratch);
}

inline const TrialAccumulatorParams *
get_trial_param_entry(const TrialParamSet *trial_params, int acc_index) {
  if (!trial_params)
    return nullptr;
  if (acc_index < 0 ||
      acc_index >= static_cast<int>(trial_params->acc_params.size())) {
    return nullptr;
  }
  const TrialAccumulatorParams &entry = trial_params->acc_params[acc_index];
  if (!entry.has_override)
    return nullptr;
  return &entry;
}

inline void resolve_event_numeric_params(
    const uuber::NativeAccumulator &acc, int acc_index,
    const TrialAccumulatorParams *override,
    const uuber::TrialParamsSoA *trial_params_soa, double &onset_out,
    double &q_out, AccDistParams &cfg_out) {
  if (trial_params_soa && trial_params_soa->valid && acc_index >= 0 &&
      acc_index < trial_params_soa->n_acc &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->dist_code.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->onset.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->q.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->t0.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p1.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p2.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p3.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p4.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p5.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p6.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p7.size() &&
      static_cast<std::size_t>(acc_index) < trial_params_soa->p8.size()) {
    const std::size_t idx = static_cast<std::size_t>(acc_index);
    onset_out = trial_params_soa->onset[idx];
    q_out = trial_params_soa->q[idx];
    cfg_out.code = trial_params_soa->dist_code[idx];
    cfg_out.t0 = trial_params_soa->t0[idx];
    cfg_out.p1 = trial_params_soa->p1[idx];
    cfg_out.p2 = trial_params_soa->p2[idx];
    cfg_out.p3 = trial_params_soa->p3[idx];
    cfg_out.p4 = trial_params_soa->p4[idx];
    cfg_out.p5 = trial_params_soa->p5[idx];
    cfg_out.p6 = trial_params_soa->p6[idx];
    cfg_out.p7 = trial_params_soa->p7[idx];
    cfg_out.p8 = trial_params_soa->p8[idx];
    return;
  }
  onset_out = override ? override->onset : acc.onset;
  q_out = override ? override->q : acc.q;
  cfg_out = override ? override->dist_cfg : acc.dist_cfg;
}

bool component_active_idx(const uuber::NativeAccumulator &acc, int component_idx,
                          const TrialAccumulatorParams *override = nullptr) {
  if (component_idx < 0) {
    return true;
  }
  const std::vector<int> *comps_idx = nullptr;
  if (override && override->has_components &&
      !override->component_indices.empty()) {
    comps_idx = &override->component_indices;
  } else if (!acc.component_indices.empty()) {
    comps_idx = &acc.component_indices;
  }
  if (comps_idx) {
    return std::find(comps_idx->begin(), comps_idx->end(), component_idx) !=
           comps_idx->end();
  }
  return true;
}

struct PoolEvalScratch {
  std::vector<std::size_t> active_indices;
  std::vector<double> density;
  std::vector<double> survival;
};

struct PoolEvalScratchStack {
  std::vector<PoolEvalScratch> stack;
  std::size_t depth{0};
};

inline PoolEvalScratchStack &pool_eval_scratch_stack() {
  thread_local PoolEvalScratchStack scratch;
  return scratch;
}

class PoolEvalScratchGuard {
public:
  PoolEvalScratchGuard()
      : stack(pool_eval_scratch_stack()), index(stack.depth++) {
    if (index >= stack.stack.size()) {
      stack.stack.emplace_back();
    }
  }

  ~PoolEvalScratchGuard() {
    if (stack.depth > 0) {
      --stack.depth;
    }
  }

  PoolEvalScratch &scratch() { return stack.stack[index]; }

private:
  PoolEvalScratchStack &stack;
  std::size_t index;
};

struct ForcedKey {
  int complete_id{0};
  int survive_id{0};
};

struct TimeKey {
  double value{0.0};
};

std::string component_cache_key(const uuber::NativeContext &ctx,
                                int component_idx,
                                const std::string &trial_key);

struct ComponentCacheEntry {
  std::string trial_type_key;
};

ComponentCacheEntry default_component_cache_entry_idx(
    const uuber::NativeContext &ctx, int component_idx) {
  ComponentCacheEntry entry;
  entry.trial_type_key = component_cache_key(ctx, component_idx, std::string());
  return entry;
}

std::vector<ComponentCacheEntry> build_component_cache_entries_from_indices(
    const uuber::NativeContext &ctx, const std::vector<int> &component_indices) {
  std::vector<ComponentCacheEntry> entries;
  entries.reserve(component_indices.size());
  for (int component_idx : component_indices) {
    entries.push_back(default_component_cache_entry_idx(ctx, component_idx));
  }
  return entries;
}

std::string component_cache_key(const uuber::NativeContext &ctx,
                                int component_idx,
                                const std::string &trial_key = std::string()) {
  if (!trial_key.empty()) {
    return trial_key;
  }
  const std::string &component = component_label_by_index_or_empty(ctx, component_idx);
  return component.empty() ? std::string("__default__") : component;
}

double component_keep_weight(const uuber::NativeContext &ctx, int component_idx,
                             int outcome_idx) {
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.component_info.size())) {
    return 1.0;
  }
  const uuber::ComponentGuessPolicy &guess =
      ctx.component_info[static_cast<std::size_t>(component_idx)].guess;
  if (!guess.valid)
    return 1.0;
  if (outcome_idx == kOutcomeIdxNA) {
    return guess.has_keep_weight_na ? guess.keep_weight_na : 1.0;
  }
  for (const auto &entry : guess.keep_weights) {
    if (entry.first == outcome_idx)
      return entry.second;
  }
  return 1.0;
}

inline const std::string &
guard_trial_type_key(const GuardEvalInput &input) {
  static const std::string kEmptyTrialTypeKey;
  return input.trial_type_key ? *input.trial_type_key : kEmptyTrialTypeKey;
}

bool guard_cdf_batch_prepared_internal(const GuardEvalInput &input,
                                       const std::vector<double> &times,
                                       std::vector<double> &cdf_out);

bool eval_node_batch_with_state_dense(
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values);

GuardEvalInput make_guard_input(const uuber::NativeContext &ctx,
                                int node_idx,
                                int component_idx,
                                const std::string *trial_type_key,
                                const TrialParamSet *trial_params,
                                const uuber::TrialParamsSoA *trial_params_soa,
                                const TimeConstraintMap *time_constraints,
                                const ForcedScopeFilter *forced_scope_filter,
                                const uuber::BitsetState *forced_complete_bits,
                                const uuber::BitsetState *forced_survive_bits,
                                const std::unordered_map<int, int>
                                    *forced_label_id_to_bit_idx,
                                const ForcedStateView *forced_state_view) {
  const uuber::IrNode &node = ir_node_required(ctx, node_idx);
  static const std::string kEmptyTrialTypeKey;
  const std::string *resolved_trial_type_key =
      trial_type_key ? trial_type_key : &kEmptyTrialTypeKey;
  bool forced_complete_bits_valid_fallback = (forced_complete_bits != nullptr);
  bool forced_survive_bits_valid_fallback = (forced_survive_bits != nullptr);
  const ForcedScopeFilter *resolved_scope_filter = forced_scope_filter;
  const uuber::BitsetState *resolved_forced_complete_bits =
      forced_complete_bits;
  const uuber::BitsetState *resolved_forced_survive_bits =
      forced_survive_bits;
  const std::unordered_map<int, int> *resolved_label_id_to_bit_idx =
      forced_label_id_to_bit_idx;
  const bool *resolved_forced_complete_bits_valid =
      &forced_complete_bits_valid_fallback;
  const bool *resolved_forced_survive_bits_valid =
      &forced_survive_bits_valid_fallback;
  if (forced_state_view) {
    if (forced_state_view->scope_filter) {
      resolved_scope_filter = forced_state_view->scope_filter;
    }
    if (forced_state_view->forced_complete_bits) {
      resolved_forced_complete_bits = forced_state_view->forced_complete_bits;
    }
    if (forced_state_view->forced_survive_bits) {
      resolved_forced_survive_bits = forced_state_view->forced_survive_bits;
    }
    if (forced_state_view->label_id_to_bit_idx) {
      resolved_label_id_to_bit_idx = forced_state_view->label_id_to_bit_idx;
    }
    if (forced_state_view->forced_complete_bits_valid) {
      resolved_forced_complete_bits_valid =
          forced_state_view->forced_complete_bits_valid;
    }
    if (forced_state_view->forced_survive_bits_valid) {
      resolved_forced_survive_bits_valid =
          forced_state_view->forced_survive_bits_valid;
    }
  }
  GuardEvalInput input{ctx,
                       node_idx,
                       -1,
                       nullptr,
                       component_idx,
                       resolved_trial_type_key,
                       {},
                       resolved_forced_complete_bits,
                       false,
                       resolved_forced_survive_bits,
                       false,
                       resolved_label_id_to_bit_idx,
                       nullptr,
                       {},
                       false,
                       trial_params,
                       trial_params_soa,
                       time_constraints};
  if (node.source_id_begin >= 0 && node.source_id_count > 0 &&
      node.source_id_begin + node.source_id_count <=
          static_cast<int>(ctx.ir.node_source_label_ids.size())) {
    input.local_scope_filter.source_ids_data =
        &ctx.ir.node_source_label_ids[static_cast<std::size_t>(
            node.source_id_begin)];
    input.local_scope_filter.source_ids_count = node.source_id_count;
  }
  if (node.source_mask_begin >= 0 && node.source_mask_count > 0 &&
      node.source_mask_begin + node.source_mask_count <=
          static_cast<int>(ctx.ir.node_source_masks.size())) {
    input.local_scope_filter.source_mask_words =
        &ctx.ir.node_source_masks[static_cast<std::size_t>(
            node.source_mask_begin)];
    input.local_scope_filter.source_mask_count = node.source_mask_count;
    input.local_scope_filter.label_id_to_bit_idx =
        resolved_label_id_to_bit_idx;
  }
  input.local_scope_filter.parent = resolved_scope_filter;
  input.forced_scope_filter = &input.local_scope_filter;
  input.forced_complete_bits_valid =
      resolved_forced_complete_bits_valid &&
      *resolved_forced_complete_bits_valid;
  input.forced_survive_bits_valid =
      resolved_forced_survive_bits_valid &&
      *resolved_forced_survive_bits_valid;
  input.forced_state = make_forced_state_view(
      input.forced_scope_filter, input.forced_complete_bits,
      &input.forced_complete_bits_valid, input.forced_survive_bits,
      &input.forced_survive_bits_valid, input.forced_label_id_to_bit_idx);
  input.has_scoped_forced =
      forced_state_intersects_scope_complete(input.forced_state) ||
      forced_state_intersects_scope_survive(input.forced_state);
  if (node.op == uuber::IrNodeOp::Guard) {
    std::vector<int> relevant_source_ids = ensure_source_ids(ctx, node);
    std::vector<int> blocker_sources =
        evaluator_gather_blocker_sources(ctx, node_idx);
    relevant_source_ids.insert(relevant_source_ids.end(),
                               blocker_sources.begin(),
                               blocker_sources.end());
    build_compact_time_constraint_lookup(time_constraints, relevant_source_ids,
                                         input.local_time_constraint_lookup);
    input.time_constraint_lookup = &input.local_time_constraint_lookup;
    if (!ctx.tree_runtime.valid ||
        node_idx < 0 ||
        node_idx >= static_cast<int>(
                        ctx.tree_runtime.node_guard_transition_idx.size())) {
      Rcpp::stop("IR guard transition metadata missing for guard node_idx=%d",
                 node_idx);
    }
    const int tr_idx = ctx.tree_runtime.node_guard_transition_idx
        [static_cast<std::size_t>(node_idx)];
    if (tr_idx < 0 ||
        tr_idx >=
            static_cast<int>(ctx.tree_runtime.guard_transitions.size())) {
      Rcpp::stop("IR guard transition mapping missing for guard node_idx=%d",
                 node_idx);
    }
    const uuber::TreeGuardTransition &tr =
        ctx.tree_runtime.guard_transitions[static_cast<std::size_t>(
            tr_idx)];
    if (tr.node_idx != node_idx) {
      Rcpp::stop("IR guard transition node mismatch for guard node_idx=%d",
                 node_idx);
    }
    input.guard_transition_idx = tr_idx;
    input.guard_transition = &tr;
  } else {
    build_compact_time_constraint_lookup(time_constraints, ensure_source_ids(ctx, node),
                                         input.local_time_constraint_lookup);
    input.time_constraint_lookup = &input.local_time_constraint_lookup;
  }
  return input;
}

uuber::TreeEventBatchEvalFn make_tree_event_eval_batch(NodeEvalState &state) {
  return [&](int event_idx, const std::vector<double> &times,
             const uuber::TreeEvalNeed &kneed,
             uuber::TreeNodeBatchValues &out) -> bool {
    out.density.assign(times.size(), 0.0);
    out.survival.assign(times.size(), 1.0);
    out.cdf.assign(times.size(), 0.0);
    if (event_idx < 0 ||
        event_idx >= static_cast<int>(state.ctx.ir.events.size())) {
      return true;
    }
    const uuber::IrEvent &event =
        state.ctx.ir.events[static_cast<std::size_t>(event_idx)];
    LabelRef ref;
    ref.label_id = event.label_id;
    ref.acc_idx = event.acc_idx;
    ref.pool_idx = event.pool_idx;
    ref.outcome_idx = event.outcome_idx;
    std::uint32_t node_flags = 0u;
    if (event.node_idx >= 0 &&
        event.node_idx < static_cast<int>(state.ctx.ir.nodes.size())) {
      node_flags =
          state.ctx.ir.nodes[static_cast<std::size_t>(event.node_idx)].flags;
    }
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
    if (!eval_event_ref_batch(state, ref, node_flags, times, need, out)) {
      return false;
    }
    for (std::size_t i = 0; i < times.size(); ++i) {
      out.density[i] = safe_density(out.density[i]);
      out.survival[i] = clamp_probability(out.survival[i]);
      out.cdf[i] = clamp_probability(out.cdf[i]);
    }
    return true;
  };
}

uuber::TreeGuardBatchEvalFn make_tree_guard_eval_batch(NodeEvalState &state) {
  return [&](const uuber::TreeEvalOp &op, const std::vector<double> &times,
             const uuber::TreeNodeBatchValues &reference_values,
             const uuber::TreeNodeBatchValues &blocker_values,
             const uuber::TreeEvalNeed &kneed,
             uuber::TreeNodeBatchValues &out) -> bool {
    const std::size_t point_count = times.size();
    if (reference_values.density.size() != point_count ||
        reference_values.survival.size() != point_count ||
        reference_values.cdf.size() != point_count ||
        blocker_values.density.size() != point_count ||
        blocker_values.survival.size() != point_count ||
        blocker_values.cdf.size() != point_count) {
      return false;
    }
    out.density.assign(point_count, 0.0);
    out.survival.assign(point_count, 1.0);
    out.cdf.assign(point_count, 0.0);
    if (state.trial_params_soa_batch == nullptr) {
      GuardEvalInput guard_input = make_guard_input_forced_state(
          state.ctx, op.node_idx, state.component_idx, &state.trial_type_key,
          state.trial_params, state.trial_params_soa, &state.time_constraints,
          state.forced_state);
      for (std::size_t i = 0; i < point_count; ++i) {
        if (kneed.density) {
          const double ref_density = safe_density(reference_values.density[i]);
          const double blocker_survival =
              clamp_probability(blocker_values.survival[i]);
          out.density[i] = safe_density(ref_density * blocker_survival);
        }
      }
      if (kneed.cdf || kneed.survival) {
        std::vector<double> cdf_values;
        if (!guard_cdf_batch_prepared_internal(guard_input, times, cdf_values) ||
            cdf_values.size() != point_count) {
          return false;
        }
        for (std::size_t i = 0; i < point_count; ++i) {
          out.cdf[i] = clamp_probability(cdf_values[i]);
          out.survival[i] = clamp_probability(1.0 - out.cdf[i]);
        }
      }
    } else {
      std::vector<const uuber::TrialParamsSoA *> subgroup_params;
      std::vector<std::vector<std::size_t>> subgroup_indices;
      subgroup_params.reserve(point_count);
      subgroup_indices.reserve(point_count);
      for (std::size_t i = 0; i < point_count; ++i) {
        const uuber::TrialParamsSoA *point_trial_params_soa =
            trial_params_soa_for_batch_point(state, i);
        std::size_t subgroup_idx = subgroup_params.size();
        for (std::size_t j = 0; j < subgroup_params.size(); ++j) {
          if (subgroup_params[j] == point_trial_params_soa) {
            subgroup_idx = j;
            break;
          }
        }
        if (subgroup_idx == subgroup_params.size()) {
          subgroup_params.push_back(point_trial_params_soa);
          subgroup_indices.emplace_back();
        }
        subgroup_indices[subgroup_idx].push_back(i);
      }
      for (std::size_t i = 0; i < point_count; ++i) {
        if (kneed.density) {
          const double ref_density = safe_density(reference_values.density[i]);
          const double blocker_survival =
              clamp_probability(blocker_values.survival[i]);
          out.density[i] = safe_density(ref_density * blocker_survival);
        }
      }
      if (kneed.cdf || kneed.survival) {
        for (std::size_t subgroup_idx = 0; subgroup_idx < subgroup_params.size();
             ++subgroup_idx) {
          const std::vector<std::size_t> &indices = subgroup_indices[subgroup_idx];
          if (indices.empty()) {
            continue;
          }
          std::vector<double> subgroup_times(indices.size(), 0.0);
          for (std::size_t j = 0; j < indices.size(); ++j) {
            subgroup_times[j] = times[indices[j]];
          }
          GuardEvalInput guard_input = make_guard_input_forced_state(
              state.ctx, op.node_idx, state.component_idx, &state.trial_type_key,
              state.trial_params, subgroup_params[subgroup_idx],
              &state.time_constraints, state.forced_state);
          std::vector<double> subgroup_cdf;
          if (!guard_cdf_batch_prepared_internal(guard_input, subgroup_times,
                                                  subgroup_cdf) ||
              subgroup_cdf.size() != indices.size()) {
            return false;
          }
          for (std::size_t j = 0; j < indices.size(); ++j) {
            const std::size_t point_idx = indices[j];
            out.cdf[point_idx] = clamp_probability(subgroup_cdf[j]);
            out.survival[point_idx] =
                clamp_probability(1.0 - out.cdf[point_idx]);
          }
        }
      }
    }
    return true;
  };
}

inline const uuber::VectorProgram &
tree_vector_program_required(const uuber::NativeContext &ctx) {
  if (!ctx.tree_program || !ctx.tree_program->valid ||
      ctx.tree_program->domain != uuber::VectorProgramDomain::Tree) {
    Rcpp::stop("IR tree vector program unavailable for node evaluation");
  }
  return *ctx.tree_program;
}

inline void init_tree_batch_values(std::size_t point_count,
                                     uuber::TreeNodeBatchValues &values) {
  values.density.assign(point_count, 0.0);
  values.survival.assign(point_count, 1.0);
  values.cdf.assign(point_count, 0.0);
}

inline bool tree_batch_values_size_matches(
    const uuber::TreeNodeBatchValues &values, std::size_t point_count) {
  return values.density.size() == point_count &&
         values.survival.size() == point_count &&
         values.cdf.size() == point_count;
}

inline void clamp_tree_batch_values(uuber::TreeNodeBatchValues &values) {
  for (std::size_t i = 0; i < values.density.size(); ++i) {
    values.density[i] = safe_density(values.density[i]);
    values.survival[i] = clamp_probability(values.survival[i]);
    values.cdf[i] = clamp_probability(values.cdf[i]);
  }
}

inline uuber::TreeEvalOp
vector_tree_guard_op(const uuber::VectorOp &op) {
  uuber::TreeEvalOp guard_op;
  guard_op.code = uuber::TreeEvalOpCode::Guard;
  guard_op.node_idx = op.node_idx;
  guard_op.out_slot = op.dst_slot;
  guard_op.event_idx = op.event_idx;
  guard_op.child_begin = op.child_begin;
  guard_op.child_count = op.child_count;
  guard_op.ref_slot = op.ref_slot;
  guard_op.blocker_slot = op.blocker_slot;
  guard_op.flags = op.flags;
  return guard_op;
}

inline bool eval_tree_vector_slots_batch(
    const uuber::VectorProgram &program, int target_slot,
    const std::vector<double> &times, const uuber::TreeEvalNeed &need,
    const uuber::TreeEventBatchEvalFn &event_eval_batch,
    const uuber::TreeGuardBatchEvalFn &guard_eval_batch,
    std::vector<uuber::TreeNodeBatchValues> &slots_out) {
  if (!program.valid || program.domain != uuber::VectorProgramDomain::Tree ||
      target_slot < 0 || !event_eval_batch) {
    return false;
  }
  if (target_slot >= static_cast<int>(program.ops.size())) {
    return false;
  }

  const std::size_t point_count = times.size();
  const uuber::TreeEvalNeed full_need{true, true, true};
  slots_out.assign(program.ops.size(), uuber::TreeNodeBatchValues{});
  std::vector<double> child_primary(
      static_cast<std::size_t>(std::max(0, program.max_child_count)), 0.0);
  std::vector<double> child_density(
      static_cast<std::size_t>(std::max(0, program.max_child_count)), 0.0);
  std::vector<double> prefix(
      static_cast<std::size_t>(std::max(0, program.max_child_count) + 1), 1.0);
  std::vector<double> suffix(
      static_cast<std::size_t>(std::max(0, program.max_child_count) + 1), 1.0);

  for (int op_idx = 0; op_idx <= target_slot; ++op_idx) {
    const uuber::VectorOp &op = program.ops[static_cast<std::size_t>(op_idx)];
    uuber::TreeNodeBatchValues out;
    init_tree_batch_values(point_count, out);
    switch (op.code) {
    case uuber::VectorOpCode::TreeEvent: {
      if (!event_eval_batch(op.event_idx, times, full_need, out) ||
          !tree_batch_values_size_matches(out, point_count)) {
        return false;
      }
      clamp_tree_batch_values(out);
      break;
    }
    case uuber::VectorOpCode::TreeAnd: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        return false;
      }
      const int n = op.child_count;
      for (int i = 0; i < n; ++i) {
        const int slot =
            program.children[static_cast<std::size_t>(op.child_begin + i)];
        if (slot < 0 || slot > target_slot) {
          child_primary[static_cast<std::size_t>(i)] = 0.0;
          child_density[static_cast<std::size_t>(i)] = 0.0;
          continue;
        }
        child_primary[static_cast<std::size_t>(i)] =
            clamp_probability(slots_out[static_cast<std::size_t>(slot)].cdf[0]);
      }
      for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
        for (int i = 0; i < n; ++i) {
          const int slot =
              program.children[static_cast<std::size_t>(op.child_begin + i)];
          if (slot < 0 || slot > target_slot) {
            child_primary[static_cast<std::size_t>(i)] = 0.0;
            child_density[static_cast<std::size_t>(i)] = 0.0;
            continue;
          }
          child_primary[static_cast<std::size_t>(i)] = clamp_probability(
              slots_out[static_cast<std::size_t>(slot)].cdf[point_idx]);
          child_density[static_cast<std::size_t>(i)] = safe_density(
              slots_out[static_cast<std::size_t>(slot)].density[point_idx]);
        }
        prefix[0] = 1.0;
        suffix[static_cast<std::size_t>(n)] = 1.0;
        for (int i = 0; i < n; ++i) {
          prefix[static_cast<std::size_t>(i + 1)] =
              prefix[static_cast<std::size_t>(i)] *
              child_primary[static_cast<std::size_t>(i)];
        }
        for (int i = n - 1; i >= 0; --i) {
          suffix[static_cast<std::size_t>(i)] =
              suffix[static_cast<std::size_t>(i + 1)] *
              child_primary[static_cast<std::size_t>(i)];
        }
        out.cdf[point_idx] =
            clamp_probability(prefix[static_cast<std::size_t>(n)]);
        out.survival[point_idx] = clamp_probability(1.0 - out.cdf[point_idx]);
        double density_sum = 0.0;
        for (int i = 0; i < n; ++i) {
          const double others =
              prefix[static_cast<std::size_t>(i)] *
              suffix[static_cast<std::size_t>(i + 1)];
          density_sum +=
              child_density[static_cast<std::size_t>(i)] *
              clamp_probability(others);
        }
        out.density[point_idx] = safe_density(density_sum);
      }
      break;
    }
    case uuber::VectorOpCode::TreeOr: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        return false;
      }
      const int n = op.child_count;
      for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
        for (int i = 0; i < n; ++i) {
          const int slot =
              program.children[static_cast<std::size_t>(op.child_begin + i)];
          if (slot < 0 || slot > target_slot) {
            child_primary[static_cast<std::size_t>(i)] = 1.0;
            child_density[static_cast<std::size_t>(i)] = 0.0;
            continue;
          }
          child_primary[static_cast<std::size_t>(i)] = clamp_probability(
              slots_out[static_cast<std::size_t>(slot)].survival[point_idx]);
          child_density[static_cast<std::size_t>(i)] = safe_density(
              slots_out[static_cast<std::size_t>(slot)].density[point_idx]);
        }
        prefix[0] = 1.0;
        suffix[static_cast<std::size_t>(n)] = 1.0;
        for (int i = 0; i < n; ++i) {
          prefix[static_cast<std::size_t>(i + 1)] =
              prefix[static_cast<std::size_t>(i)] *
              child_primary[static_cast<std::size_t>(i)];
        }
        for (int i = n - 1; i >= 0; --i) {
          suffix[static_cast<std::size_t>(i)] =
              suffix[static_cast<std::size_t>(i + 1)] *
              child_primary[static_cast<std::size_t>(i)];
        }
        out.survival[point_idx] =
            clamp_probability(prefix[static_cast<std::size_t>(n)]);
        out.cdf[point_idx] = clamp_probability(1.0 - out.survival[point_idx]);
        double density_sum = 0.0;
        for (int i = 0; i < n; ++i) {
          const double others =
              prefix[static_cast<std::size_t>(i)] *
              suffix[static_cast<std::size_t>(i + 1)];
          density_sum +=
              child_density[static_cast<std::size_t>(i)] *
              clamp_probability(others);
        }
        out.density[point_idx] = safe_density(density_sum);
      }
      break;
    }
    case uuber::VectorOpCode::TreeNot: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin >= static_cast<int>(program.children.size())) {
        return false;
      }
      const int child_slot =
          program.children[static_cast<std::size_t>(op.child_begin)];
      if (child_slot < 0 || child_slot > target_slot) {
        return false;
      }
      const uuber::TreeNodeBatchValues &child =
          slots_out[static_cast<std::size_t>(child_slot)];
      for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
        out.cdf[point_idx] = clamp_probability(1.0 - child.cdf[point_idx]);
        out.survival[point_idx] = clamp_probability(child.cdf[point_idx]);
        out.density[point_idx] = 0.0;
      }
      break;
    }
    case uuber::VectorOpCode::TreeGuard: {
      uuber::TreeNodeBatchValues ref_values;
      uuber::TreeNodeBatchValues blocker_values;
      init_tree_batch_values(point_count, ref_values);
      init_tree_batch_values(point_count, blocker_values);
      if (op.ref_slot >= 0 && op.ref_slot <= target_slot) {
        ref_values = slots_out[static_cast<std::size_t>(op.ref_slot)];
      }
      if (op.blocker_slot >= 0 && op.blocker_slot <= target_slot) {
        blocker_values = slots_out[static_cast<std::size_t>(op.blocker_slot)];
      }
      const uuber::TreeEvalOp guard_op = vector_tree_guard_op(op);
      if (guard_eval_batch) {
        if (!guard_eval_batch(guard_op, times, ref_values, blocker_values,
                              full_need, out) ||
            !tree_batch_values_size_matches(out, point_count)) {
          return false;
        }
        clamp_tree_batch_values(out);
      } else {
        for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
          const double ref_density = safe_density(ref_values.density[point_idx]);
          const double blocker_survival =
              clamp_probability(blocker_values.survival[point_idx]);
          out.density[point_idx] = safe_density(ref_density * blocker_survival);
          out.cdf[point_idx] = clamp_probability(
              ref_values.cdf[point_idx] * blocker_survival);
          out.survival[point_idx] =
              clamp_probability(1.0 - out.cdf[point_idx]);
        }
      }
      break;
    }
    default:
      return false;
    }
    slots_out[static_cast<std::size_t>(op_idx)] = std::move(out);
  }

  return true;
}

inline bool eval_tree_vector_target_slot(const uuber::VectorProgram &program,
                                         int target_node_idx,
                                         int &target_slot) {
  if (target_node_idx < 0 ||
      target_node_idx >=
          static_cast<int>(program.outputs.node_idx_to_slot.size())) {
    return false;
  }
  target_slot =
      program.outputs.node_idx_to_slot[static_cast<std::size_t>(target_node_idx)];
  return target_slot >= 0 && target_slot < static_cast<int>(program.ops.size());
}

inline bool eval_tree_vector_node_batch(
    const uuber::VectorProgram &program, int target_node_idx,
    const std::vector<double> &times, const uuber::TreeEvalNeed &need,
    const uuber::TreeEventBatchEvalFn &event_eval_batch,
    const uuber::TreeGuardBatchEvalFn &guard_eval_batch,
    uuber::TreeNodeBatchValues &out_values) {
  out_values = uuber::TreeNodeBatchValues{};
  int target_slot = -1;
  if (!eval_tree_vector_target_slot(program, target_node_idx, target_slot)) {
    return false;
  }
  std::vector<uuber::TreeNodeBatchValues> slots;
  if (!eval_tree_vector_slots_batch(program, target_slot, times, need,
                                    event_eval_batch, guard_eval_batch, slots)) {
    return false;
  }
  out_values = slots[static_cast<std::size_t>(target_slot)];
  return tree_batch_values_size_matches(out_values, times.size());
}

inline bool eval_tree_vector_nodes_batch(
    const uuber::VectorProgram &program, const std::vector<int> &target_node_indices,
    const std::vector<double> &times, const uuber::TreeEvalNeed &need,
    const uuber::TreeEventBatchEvalFn &event_eval_batch,
    const uuber::TreeGuardBatchEvalFn &guard_eval_batch,
    std::vector<uuber::TreeNodeBatchValues> &out_values) {
  out_values.clear();
  if (target_node_indices.empty()) {
    return true;
  }
  int max_target_slot = -1;
  std::vector<int> target_slots;
  target_slots.reserve(target_node_indices.size());
  for (int target_node_idx : target_node_indices) {
    int target_slot = -1;
    if (!eval_tree_vector_target_slot(program, target_node_idx, target_slot)) {
      return false;
    }
    target_slots.push_back(target_slot);
    max_target_slot = std::max(max_target_slot, target_slot);
  }
  std::vector<uuber::TreeNodeBatchValues> slots;
  if (!eval_tree_vector_slots_batch(program, max_target_slot, times, need,
                                    event_eval_batch, guard_eval_batch, slots)) {
    return false;
  }
  out_values.reserve(target_slots.size());
  for (int target_slot : target_slots) {
    out_values.push_back(slots[static_cast<std::size_t>(target_slot)]);
    if (!tree_batch_values_size_matches(out_values.back(), times.size())) {
      return false;
    }
  }
  return true;
}

bool eval_node_batch_with_state_dense(
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values) {
  out_values = uuber::TreeNodeBatchValues{};
  if (times.empty()) {
    return true;
  }
  uuber::TreeEvalNeed tree_need;
  tree_need.density = needs_density(need);
  tree_need.survival = needs_survival(need);
  tree_need.cdf = needs_cdf(need);
  if (!tree_need.density && !tree_need.survival && !tree_need.cdf) {
    tree_need.density = true;
    tree_need.survival = true;
    tree_need.cdf = true;
  }
  uuber::TreeEventBatchEvalFn event_eval_batch_cb =
      make_tree_event_eval_batch(state);
  uuber::TreeGuardBatchEvalFn guard_eval_batch_cb =
      make_tree_guard_eval_batch(state);
  return eval_tree_vector_node_batch(tree_vector_program_required(state.ctx),
                                     node_idx, times, tree_need,
                                     event_eval_batch_cb,
                                     guard_eval_batch_cb, out_values);
}

inline const uuber::TreeGuardTransition &
guard_transition_required(const GuardEvalInput &input,
                          bool require_blocker_idx = false) {
  if (input.guard_transition != nullptr) {
    if (input.guard_transition->node_idx != input.node_idx) {
      Rcpp::stop("IR guard transition node mismatch for guard node_idx=%d",
                 input.node_idx);
    }
    if (require_blocker_idx && input.guard_transition->blocker_node_idx < 0) {
      Rcpp::stop(
          "IR guard transition blocker metadata invalid for guard node_idx=%d",
          input.node_idx);
    }
    return *input.guard_transition;
  }
  if (!input.ctx.tree_runtime.valid ||
      input.node_idx < 0 ||
      input.node_idx >= static_cast<int>(
                            input.ctx.tree_runtime.node_guard_transition_idx
                                .size())) {
    Rcpp::stop("IR guard transition metadata missing for guard node_idx=%d",
               input.node_idx);
  }
  const int tr_idx = input.ctx.tree_runtime.node_guard_transition_idx
      [static_cast<std::size_t>(input.node_idx)];
  if (tr_idx < 0 ||
      tr_idx >=
          static_cast<int>(input.ctx.tree_runtime.guard_transitions.size())) {
    Rcpp::stop("IR guard transition mapping missing for guard node_idx=%d",
               input.node_idx);
  }
  const uuber::TreeGuardTransition &tr =
      input.ctx.tree_runtime.guard_transitions[static_cast<std::size_t>(
          tr_idx)];
  if (tr.node_idx != input.node_idx) {
    Rcpp::stop("IR guard transition node mismatch for guard node_idx=%d",
               input.node_idx);
  }
  if (require_blocker_idx && tr.blocker_node_idx < 0) {
    Rcpp::stop("IR guard transition blocker metadata invalid for guard node_idx=%d",
               input.node_idx);
  }
  return tr;
}

// --- Nested Guard Optimization ---

struct LinearGuardChain {
  const int *reference_indices{nullptr}; // Dense IR node indices [Ref0, Ref1, ...]
  std::size_t reference_count{0u};
  int leaf_blocker_idx{-1};          // Dense IR node index
  bool valid{false};
};

struct FastEventInfo {
  const uuber::NativeAccumulator *acc{nullptr};
  const TrialAccumulatorParams *override{nullptr};
  int acc_idx{-1};
  double onset{0.0};
  double q{0.0};
  AccDistParams cfg{};
  LowerBoundTransform lower_bound{};
  int outcome_idx{-1};
  bool component_ok{false};
};

struct PreparedLinearGuardChainInfo {
  LinearGuardChain chain{};
  std::vector<FastEventInfo> ref_infos;
  std::vector<std::uint8_t> ref_info_valid;
  std::vector<bool> ref_density_ok;
  FastEventInfo leaf_info{};
  bool leaf_info_valid{false};
  bool valid{false};
};

struct GuardPreparedConditioningTemplate {
  LinearGuardChain chain{};
  std::vector<int> conditioning_label_ids;
  bool valid{false};
};

struct GuardPreparedCacheKey {
  int node_idx{-1};
  int component_idx{-1};
  std::uintptr_t trial_params_bits{0u};
  std::uintptr_t trial_params_soa_bits{0u};
  std::vector<std::uint8_t> conditioning_pattern;

  bool operator==(const GuardPreparedCacheKey &other) const noexcept {
    return node_idx == other.node_idx &&
           component_idx == other.component_idx &&
           trial_params_bits == other.trial_params_bits &&
           trial_params_soa_bits == other.trial_params_soa_bits &&
           conditioning_pattern == other.conditioning_pattern;
  }
};

struct GuardPreparedCacheKeyHash {
  std::size_t operator()(const GuardPreparedCacheKey &key) const noexcept {
    std::uint64_t hash = kFNV64Offset;
    hash_append_u64(hash, static_cast<std::uint64_t>(key.node_idx));
    hash_append_u64(hash, static_cast<std::uint64_t>(key.component_idx));
    hash_append_u64(hash, static_cast<std::uint64_t>(key.trial_params_bits));
    hash_append_u64(hash, static_cast<std::uint64_t>(key.trial_params_soa_bits));
    for (std::uint8_t state : key.conditioning_pattern) {
      hash_append_u64(hash, static_cast<std::uint64_t>(state));
    }
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

using GuardPreparedConditioningTemplateMap =
    std::unordered_map<int, GuardPreparedConditioningTemplate>;
using GuardPreparedCacheMap =
    std::unordered_map<GuardPreparedCacheKey, PreparedLinearGuardChainInfo,
                       GuardPreparedCacheKeyHash>;

inline const SourceTimeConstraint *guard_time_constraint_find(
    const GuardEvalInput &input, int label_id) {
  if (label_id < 0 || label_id == NA_INTEGER) {
    return nullptr;
  }
  const SourceTimeConstraint *constraint =
      compact_time_constraint_find(input.time_constraint_lookup, label_id);
  if (constraint != nullptr) {
    if (!time_constraint_empty(*constraint)) {
      return constraint;
    }
    if (input.time_constraints == nullptr) {
      return constraint;
    }
  }
  return time_constraints_find(input.time_constraints, label_id);
}

inline bool guard_label_has_time_constraint(const GuardEvalInput &input,
                                            int label_id) {
  const SourceTimeConstraint *constraint =
      guard_time_constraint_find(input, label_id);
  return constraint != nullptr &&
         (constraint->has_exact || constraint->has_lower || constraint->has_upper);
}

inline bool guard_label_is_forced(const GuardEvalInput &input, int label_id) {
  if (label_id < 0 || label_id == NA_INTEGER) {
    return false;
  }
  return forced_state_contains_complete(input.forced_state, label_id) ||
         forced_state_contains_survive(input.forced_state, label_id);
}

inline bool fast_event_density_supported(const GuardEvalInput &input,
                                         const FastEventInfo &info);
LinearGuardChain linear_guard_chain_from_transition(
    const uuber::NativeContext &ctx, const uuber::TreeGuardTransition &tr);
bool resolve_prepared_linear_guard_chain_info(
    const GuardEvalInput &input, PreparedLinearGuardChainInfo &scratch,
    const PreparedLinearGuardChainInfo *&prepared);

bool fast_event_info_dense_idx(const GuardEvalInput &input, int node_idx,
                               FastEventInfo &info) {
  const uuber::IrNode &node = ir_node_required(input.ctx, node_idx);
  if (!(node.op == uuber::IrNodeOp::EventAcc || node.op == uuber::IrNodeOp::EventPool)) {
    return false;
  }
  if ((node.flags & uuber::IR_NODE_FLAG_SPECIAL_DEADLINE) != 0u ||
      (node.flags & uuber::IR_NODE_FLAG_SPECIAL_GUESS) != 0u) {
    return false;
  }
  LabelRef ref = node_label_ref(input.ctx, node);
  if (guard_label_is_forced(input, ref.label_id) ||
      guard_label_has_time_constraint(input, ref.label_id)) {
    return false;
  }

  int acc_idx = ref.acc_idx;
  if (acc_idx < 0 ||
      acc_idx >= static_cast<int>(input.ctx.accumulators.size())) {
    return false;
  }
  info.acc_idx = acc_idx;
  info.acc = &input.ctx.accumulators[acc_idx];
  info.override = get_trial_param_entry(input.trial_params, acc_idx);
  resolve_event_numeric_params(*info.acc, acc_idx, info.override,
                               input.trial_params_soa, info.onset, info.q,
                               info.cfg);
  const uuber::LabelRef onset_source_ref = evaluator_make_onset_source_ref(
      input.ctx, info.acc->onset_kind, info.acc->onset_source_acc_idx,
      info.acc->onset_source_pool_idx);
  if (guard_label_is_forced(input, onset_source_ref.label_id) ||
      guard_label_has_time_constraint(input, onset_source_ref.label_id)) {
    return false;
  }
  info.lower_bound = default_lower_bound_transform(info.cfg);
  info.outcome_idx = ref.outcome_idx;
  info.component_ok =
      component_active_idx(*info.acc, input.component_idx, info.override);
  return true;
}

inline bool fast_event_density_supported(const GuardEvalInput &input,
                                         const FastEventInfo &info) {
  if (info.outcome_idx < 0 ||
      info.outcome_idx >= static_cast<int>(input.ctx.outcome_info.size())) {
    return true;
  }
  const uuber::OutcomeContextInfo &outcome =
      input.ctx.outcome_info[static_cast<std::size_t>(info.outcome_idx)];
  return outcome.alias_sources.empty() && outcome.guess_donors.empty();
}

inline bool eval_node_batch_from_guard_input(
    const GuardEvalInput &input, int node_idx, const std::vector<double> &times,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values) {
  if (times.empty()) {
    out_values = uuber::TreeNodeBatchValues{};
    return true;
  }
  NodeEvalState local(input.ctx, times.front(), input.component_idx,
                      input.trial_params, guard_trial_type_key(input), false,
                      -1, input.time_constraints,
                      nullptr, nullptr, false, nullptr, false,
                      &input.forced_state);
  local.trial_params_soa = input.trial_params_soa;
  return eval_node_batch_with_state_dense(node_idx, times, local, need,
                                          out_values);
}

LinearGuardChain linear_guard_chain_from_transition(
    const uuber::NativeContext &ctx, const uuber::TreeGuardTransition &tr) {
  LinearGuardChain chain;
  if (tr.eval_mode != uuber::TreeGuardEvalMode::LinearChainODE) {
    return chain;
  }
  if (tr.linear_chain_begin < 0 || tr.linear_chain_count <= 0 ||
      tr.linear_chain_begin + tr.linear_chain_count >
          static_cast<int>(ctx.tree_runtime.guard_linear_chain_nodes.size()) ||
      tr.linear_chain_leaf_idx < 0 ||
      tr.linear_chain_leaf_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return {};
  }
  const int begin = tr.linear_chain_begin;
  chain.reference_indices =
      &ctx.tree_runtime
           .guard_linear_chain_nodes[static_cast<std::size_t>(begin)];
  chain.reference_count = static_cast<std::size_t>(tr.linear_chain_count);
  chain.leaf_blocker_idx = tr.linear_chain_leaf_idx;
  if (chain.reference_count == 0u) {
    return {};
  }
  chain.valid = true;
  return chain;
}

inline void append_guard_conditioning_labels_for_node(
    const uuber::NativeContext &ctx, int node_idx, std::vector<int> &out) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return;
  }
  const uuber::IrNode &node = ir_node_required(ctx, node_idx);
  if (!(node.op == uuber::IrNodeOp::EventAcc ||
        node.op == uuber::IrNodeOp::EventPool)) {
    return;
  }
  const LabelRef ref = node_label_ref(ctx, node);
  if (ref.label_id >= 0 && ref.label_id != NA_INTEGER) {
    out.push_back(ref.label_id);
  }
  if (ref.acc_idx < 0 ||
      ref.acc_idx >= static_cast<int>(ctx.accumulators.size())) {
    return;
  }
  const uuber::NativeAccumulator &acc =
      ctx.accumulators[static_cast<std::size_t>(ref.acc_idx)];
  const LabelRef onset_source_ref = evaluator_make_onset_source_ref(
      ctx, acc.onset_kind, acc.onset_source_acc_idx, acc.onset_source_pool_idx);
  if (onset_source_ref.label_id >= 0 &&
      onset_source_ref.label_id != NA_INTEGER) {
    out.push_back(onset_source_ref.label_id);
  }
}

inline bool build_guard_prepared_conditioning_template(
    const GuardEvalInput &input, GuardPreparedConditioningTemplate &out) {
  out = GuardPreparedConditioningTemplate{};
  const uuber::TreeGuardTransition &tr = guard_transition_required(input);
  out.chain = linear_guard_chain_from_transition(input.ctx, tr);
  if (!out.chain.valid) {
    return false;
  }
  for (std::size_t i = 0; i < out.chain.reference_count; ++i) {
    append_guard_conditioning_labels_for_node(input.ctx,
                                              out.chain.reference_indices[i],
                                              out.conditioning_label_ids);
  }
  append_guard_conditioning_labels_for_node(input.ctx, out.chain.leaf_blocker_idx,
                                            out.conditioning_label_ids);
  sort_unique(out.conditioning_label_ids);
  out.valid = true;
  return true;
}

inline std::unordered_map<std::uint64_t, GuardPreparedConditioningTemplateMap> &
guard_prepared_conditioning_template_registry() {
  thread_local std::unordered_map<std::uint64_t,
                                  GuardPreparedConditioningTemplateMap>
      cache_by_runtime_id;
  return cache_by_runtime_id;
}

inline std::unordered_map<std::uint64_t, GuardPreparedCacheMap> &
guard_prepared_cache_registry() {
  thread_local std::unordered_map<std::uint64_t, GuardPreparedCacheMap>
      cache_by_runtime_id;
  return cache_by_runtime_id;
}

inline GuardPreparedConditioningTemplateMap &
guard_prepared_conditioning_template_cache_for_ctx(
    const uuber::NativeContext &ctx) {
  return guard_prepared_conditioning_template_registry()
      [ctx.runtime_cache_instance_id];
}

inline GuardPreparedCacheMap &
guard_prepared_cache_for_ctx(const uuber::NativeContext &ctx) {
  return guard_prepared_cache_registry()[ctx.runtime_cache_instance_id];
}

inline bool guard_prepared_conditioning_template_for_input(
    const GuardEvalInput &input, GuardPreparedConditioningTemplate &out) {
  auto &template_cache = guard_prepared_conditioning_template_cache_for_ctx(
      input.ctx);
  auto it = template_cache.find(input.node_idx);
  if (it != template_cache.end()) {
    out = it->second;
    return out.valid;
  }
  GuardPreparedConditioningTemplate built;
  if (!build_guard_prepared_conditioning_template(input, built)) {
    return false;
  }
  auto inserted = template_cache.emplace(input.node_idx, std::move(built));
  out = inserted.first->second;
  return out.valid;
}

inline bool build_guard_prepared_cache_key(
    const GuardEvalInput &input,
    const GuardPreparedConditioningTemplate &conditioning_template,
    GuardPreparedCacheKey &out) {
  if (!conditioning_template.valid) {
    return false;
  }
  out = GuardPreparedCacheKey{};
  out.node_idx = input.node_idx;
  out.component_idx = input.component_idx;
  out.trial_params_bits =
      reinterpret_cast<std::uintptr_t>(input.trial_params);
  out.trial_params_soa_bits =
      reinterpret_cast<std::uintptr_t>(input.trial_params_soa);
  out.conditioning_pattern.reserve(
      conditioning_template.conditioning_label_ids.size());
  for (int label_id : conditioning_template.conditioning_label_ids) {
    std::uint8_t state = 0u;
    if (guard_label_is_forced(input, label_id)) {
      state |= 0x1u;
    }
    if (guard_label_has_time_constraint(input, label_id)) {
      state |= 0x2u;
    }
    out.conditioning_pattern.push_back(state);
  }
  return true;
}

inline bool prepare_linear_guard_chain_info_uncached(
    const GuardEvalInput &input,
    const GuardPreparedConditioningTemplate &conditioning_template,
    PreparedLinearGuardChainInfo &out) {
  out = PreparedLinearGuardChainInfo{};
  out.chain = conditioning_template.chain;
  if (!out.chain.valid) {
    return false;
  }
  const std::size_t depth = out.chain.reference_count;
  out.ref_infos.resize(depth);
  out.ref_info_valid.assign(depth, 0u);
  out.ref_density_ok.assign(depth, false);
  for (std::size_t i = 0; i < depth; ++i) {
    if (fast_event_info_dense_idx(input, out.chain.reference_indices[i],
                                  out.ref_infos[i])) {
      out.ref_info_valid[i] = 1u;
      out.ref_density_ok[i] =
          fast_event_density_supported(input, out.ref_infos[i]);
    }
  }
  out.leaf_info_valid =
      fast_event_info_dense_idx(input, out.chain.leaf_blocker_idx, out.leaf_info);
  out.valid = true;
  return true;
}

bool resolve_prepared_linear_guard_chain_info(
    const GuardEvalInput &input, PreparedLinearGuardChainInfo &scratch,
    const PreparedLinearGuardChainInfo *&prepared) {
  prepared = nullptr;
  GuardPreparedConditioningTemplate conditioning_template;
  if (!guard_prepared_conditioning_template_for_input(input,
                                                      conditioning_template)) {
    return false;
  }
  GuardPreparedCacheKey cache_key;
  if (!build_guard_prepared_cache_key(input, conditioning_template, cache_key)) {
    return false;
  }
  auto &cache = guard_prepared_cache_for_ctx(input.ctx);
  auto it = cache.find(cache_key);
  if (it != cache.end()) {
    prepared = &it->second;
    return true;
  }
  if (!prepare_linear_guard_chain_info_uncached(input, conditioning_template,
                                                scratch)) {
    return false;
  }
  auto inserted =
      cache.emplace(std::move(cache_key), std::move(scratch));
  prepared = &inserted.first->second;
  return true;
}

constexpr int kLinearChainBaseSteps = 88;
constexpr int kLinearChainDepthStep = 24;
constexpr int kLinearChainTolBonus = 16;
constexpr int kLinearChainMaxSteps = 240;

inline int linear_chain_step_budget(std::size_t depth,
                                    const CanonicalIntegrationSettings &canonical) {
  int steps = kLinearChainBaseSteps;
  if (depth > 3u) {
    steps += static_cast<int>(depth - 3u) * kLinearChainDepthStep;
  }
  if (canonical.effective_tol <= 1e-7) {
    steps += kLinearChainTolBonus;
  }
  if (canonical.effective_tol <= 1e-8) {
    steps += kLinearChainTolBonus;
  }
  return std::min(steps, kLinearChainMaxSteps);
}

inline bool eval_fast_event_density_batch_from_info(
    const FastEventInfo &info, const std::vector<double> &times,
    std::vector<double> &shifted, std::vector<double> &values,
    std::vector<double> &density_out) {
  density_out.assign(times.size(), 0.0);
  if (!info.component_ok) {
    return true;
  }
  const double onset_eff = total_onset_with_t0(info.onset, info.cfg);
  const double success_prob = clamp(1.0 - info.q, 0.0, 1.0);
  shifted.assign(times.size(), 0.0);
  values.assign(times.size(), 0.0);
  for (std::size_t i = 0; i < times.size(); ++i) {
    shifted[i] = times[i] - onset_eff;
  }
  eval_pdf_vec_with_lower_bound(info.cfg, info.lower_bound, shifted.data(),
                                times.size(), values.data());
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (!std::isfinite(times[i]) || times[i] < onset_eff ||
        !(success_prob > 0.0)) {
      continue;
    }
    density_out[i] = safe_density(success_prob * safe_density(values[i]));
  }
  return true;
}

inline bool eval_fast_event_survival_batch_from_info(
    const FastEventInfo &info, const std::vector<double> &times,
    std::vector<double> &shifted, std::vector<double> &values,
    std::vector<double> &survival_out) {
  survival_out.assign(times.size(), 1.0);
  if (!info.component_ok) {
    return true;
  }
  const double onset_eff = total_onset_with_t0(info.onset, info.cfg);
  const double success_prob = clamp(1.0 - info.q, 0.0, 1.0);
  shifted.assign(times.size(), 0.0);
  values.assign(times.size(), 0.0);
  for (std::size_t i = 0; i < times.size(); ++i) {
    shifted[i] = times[i] - onset_eff;
  }
  eval_cdf_vec_with_lower_bound(info.cfg, info.lower_bound, shifted.data(),
                                times.size(), values.data());
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (!std::isfinite(times[i]) || times[i] < onset_eff) {
      continue;
    }
    survival_out[i] =
        clamp_probability(info.q + success_prob * (1.0 - clamp_probability(values[i])));
  }
  return true;
}

inline bool guard_cdf_simple_pair_batch_prepared(
    const GuardEvalInput &input, const PreparedLinearGuardChainInfo &prepared,
    const std::vector<double> &times, std::vector<double> &cdf_out) {
  cdf_out.assign(times.size(), 0.0);
  if (times.empty() || !prepared.valid || !prepared.chain.valid ||
      prepared.chain.reference_count != 1u) {
    return false;
  }

  std::vector<std::uint8_t> valid_points(times.size(), 0u);
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (std::isfinite(times[i]) && times[i] > 0.0) {
      valid_points[i] = 1u;
    }
  }
  std::vector<TimeIntegralPointPlan> quadrature_plans;
  std::vector<double> quadrature_times;
  std::vector<double> quadrature_weights;
  const double quadrature_upper = max_positive_finite_upper(times);
  if (quadrature_upper > 0.0 &&
      !build_time_integral_query_plan(
          times, 0.0, quadrature_upper, 4, quadrature_plans,
          quadrature_times, quadrature_weights, nullptr, nullptr,
          &valid_points)) {
    return false;
  }

  std::vector<double> shifted;
  std::vector<double> values;
  std::vector<double> reference_density;
  std::vector<double> blocker_survival;
  const bool fast_ref_ok =
      !prepared.ref_infos.empty() &&
      !prepared.ref_info_valid.empty() && prepared.ref_info_valid[0] != 0u &&
      !prepared.ref_density_ok.empty() && prepared.ref_density_ok[0];
  const bool fast_leaf_ok = prepared.leaf_info_valid;
  auto try_eval_conditioned_guard_event_batch =
      [&](int node_idx, EvalNeed need,
          uuber::TreeNodeBatchValues &values_out) -> bool {
    if (node_idx < 0 ||
        node_idx >= static_cast<int>(input.ctx.ir.nodes.size())) {
      return false;
    }
    const uuber::IrNode &node =
        input.ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
    if ((node.op != uuber::IrNodeOp::EventAcc &&
         node.op != uuber::IrNodeOp::EventPool) ||
        node.event_idx < 0 ||
        node.event_idx >= static_cast<int>(input.ctx.ir.events.size())) {
      return false;
    }
    NodeEvalState local(input.ctx, quadrature_times.empty() ? 0.0
                                                            : quadrature_times.front(),
                        input.component_idx, input.trial_params,
                        guard_trial_type_key(input), false, -1,
                        input.time_constraints, nullptr, nullptr, false,
                        nullptr, false, &input.forced_state);
    local.trial_params_soa = input.trial_params_soa;
    return try_eval_conditioned_simple_acc_event_batch(
        local, input.ctx.ir.events[static_cast<std::size_t>(node.event_idx)],
        quadrature_times, need, values_out);
  };

  if (fast_ref_ok) {
    if (!eval_fast_event_density_batch_from_info(prepared.ref_infos[0],
                                                 quadrature_times, shifted,
                                                 values, reference_density)) {
      return false;
    }
  } else {
    uuber::TreeNodeBatchValues reference_values;
    if ((!try_eval_conditioned_guard_event_batch(
             prepared.chain.reference_indices[0], EvalNeed::kDensity,
             reference_values) &&
         !eval_node_batch_from_guard_input(
             input, prepared.chain.reference_indices[0], quadrature_times,
             EvalNeed::kDensity, reference_values)) ||
        reference_values.density.size() != quadrature_times.size()) {
      return false;
    }
    reference_density.resize(reference_values.density.size());
    for (std::size_t i = 0; i < reference_values.density.size(); ++i) {
      reference_density[i] = safe_density(reference_values.density[i]);
    }
  }

  if (fast_leaf_ok) {
    if (!eval_fast_event_survival_batch_from_info(prepared.leaf_info,
                                                  quadrature_times, shifted,
                                                  values, blocker_survival)) {
      return false;
    }
  } else {
    uuber::TreeNodeBatchValues blocker_values;
    if ((!try_eval_conditioned_guard_event_batch(
             prepared.chain.leaf_blocker_idx, EvalNeed::kSurvival,
             blocker_values) &&
         !eval_node_batch_from_guard_input(
             input, prepared.chain.leaf_blocker_idx, quadrature_times,
             EvalNeed::kSurvival, blocker_values)) ||
        blocker_values.survival.size() != quadrature_times.size()) {
      return false;
    }
    blocker_survival.resize(blocker_values.survival.size());
    for (std::size_t i = 0; i < blocker_values.survival.size(); ++i) {
      blocker_survival[i] = clamp_probability(blocker_values.survival[i]);
    }
  }

  for (std::size_t i = 0; i < times.size(); ++i) {
    if (valid_points[i] == 0u) {
      continue;
    }
    const TimeIntegralPointPlan &plan = quadrature_plans[i];
    double total = 0.0;
    for (std::size_t query_idx = plan.begin;
         query_idx < plan.begin + plan.count; ++query_idx) {
      const double weight = quadrature_weights[query_idx];
      if (!std::isfinite(weight) || weight <= 0.0) {
        continue;
      }
      const double term = safe_density(reference_density[query_idx] *
                                       blocker_survival[query_idx]);
      if (term > 0.0) {
        total += weight * term;
      }
    }
    cdf_out[i] = clamp_probability(total);
  }
  return true;
}

inline std::uint64_t guard_time_bits(double value) {
  std::uint64_t bits = 0u;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline bool guard_cdf_general_batch_prepared(
    const GuardEvalInput &input, const PreparedLinearGuardChainInfo &prepared,
    const std::vector<double> &times, std::vector<double> &cdf_out) {
  cdf_out.assign(times.size(), 0.0);
  if (times.empty()) {
    return true;
  }
  if (!prepared.valid || !prepared.chain.valid ||
      prepared.chain.reference_count < 1u) {
    return false;
  }

  const CanonicalIntegrationSettings canonical =
      canonicalize_integration_settings(IntegrationSettings{});
  const std::size_t depth = prepared.chain.reference_count;
  const int n_steps = linear_chain_step_budget(depth, canonical);
  if (n_steps < 1) {
    return false;
  }
  const std::size_t grid_points = static_cast<std::size_t>(2 * n_steps + 1);
  const std::size_t value_width = depth + 1u;

  std::vector<int> time_to_unique(times.size(), -1);
  std::vector<double> unique_times;
  unique_times.reserve(times.size());
  std::unordered_map<std::uint64_t, int> unique_time_lookup;
  unique_time_lookup.reserve(times.size());
  for (std::size_t i = 0; i < times.size(); ++i) {
    const double t = times[i];
    if (!std::isfinite(t)) {
      cdf_out[i] = 1.0;
      continue;
    }
    if (!(t > 0.0)) {
      continue;
    }
    const std::uint64_t bits = guard_time_bits(t);
    auto it = unique_time_lookup.find(bits);
    if (it == unique_time_lookup.end()) {
      const int unique_idx = static_cast<int>(unique_times.size());
      unique_times.push_back(t);
      unique_time_lookup.emplace(bits, unique_idx);
      time_to_unique[i] = unique_idx;
    } else {
      time_to_unique[i] = it->second;
    }
  }
  if (unique_times.empty()) {
    return true;
  }

  std::vector<double> query_times(unique_times.size() * grid_points, 0.0);
  for (std::size_t unique_idx = 0; unique_idx < unique_times.size();
       ++unique_idx) {
    const double h = unique_times[unique_idx] / static_cast<double>(n_steps);
    const std::size_t base = unique_idx * grid_points;
    for (std::size_t grid_idx = 0; grid_idx < grid_points; ++grid_idx) {
      query_times[base + grid_idx] =
          0.5 * h * static_cast<double>(grid_idx);
    }
  }

  std::vector<double> grid_values(query_times.size() * value_width, 0.0);
  std::vector<double> shifted;
  std::vector<double> values;
  auto grid_slot = [&](std::size_t query_idx) {
    return grid_values.data() + query_idx * value_width;
  };

  for (std::size_t ref_idx = 0; ref_idx < depth; ++ref_idx) {
    const bool fast_ref_ok =
        ref_idx < prepared.ref_infos.size() &&
        ref_idx < prepared.ref_info_valid.size() &&
        prepared.ref_info_valid[ref_idx] != 0u &&
        ref_idx < prepared.ref_density_ok.size() &&
        prepared.ref_density_ok[ref_idx];
    if (fast_ref_ok) {
      std::vector<double> density_values;
      if (!eval_fast_event_density_batch_from_info(
              prepared.ref_infos[ref_idx], query_times, shifted, values,
              density_values) ||
          density_values.size() != query_times.size()) {
        return false;
      }
      for (std::size_t query_idx = 0; query_idx < query_times.size();
           ++query_idx) {
        grid_slot(query_idx)[ref_idx] = safe_density(density_values[query_idx]);
      }
    } else {
      uuber::TreeNodeBatchValues ref_batch;
      if (!eval_node_batch_from_guard_input(
              input, prepared.chain.reference_indices[ref_idx], query_times,
              EvalNeed::kDensity, ref_batch) ||
          ref_batch.density.size() != query_times.size()) {
        return false;
      }
      for (std::size_t query_idx = 0; query_idx < query_times.size();
           ++query_idx) {
        grid_slot(query_idx)[ref_idx] =
            safe_density(ref_batch.density[query_idx]);
      }
    }
  }

  if (prepared.leaf_info_valid) {
    std::vector<double> survival_values;
    if (!eval_fast_event_survival_batch_from_info(
            prepared.leaf_info, query_times, shifted, values, survival_values) ||
        survival_values.size() != query_times.size()) {
      return false;
    }
    for (std::size_t query_idx = 0; query_idx < query_times.size();
         ++query_idx) {
      grid_slot(query_idx)[depth] =
          clamp_probability(survival_values[query_idx]);
    }
  } else {
    uuber::TreeNodeBatchValues leaf_batch;
    if (!eval_node_batch_from_guard_input(
            input, prepared.chain.leaf_blocker_idx, query_times,
            EvalNeed::kSurvival, leaf_batch) ||
        leaf_batch.survival.size() != query_times.size()) {
      return false;
    }
    for (std::size_t query_idx = 0; query_idx < query_times.size();
         ++query_idx) {
      grid_slot(query_idx)[depth] =
          clamp_probability(leaf_batch.survival[query_idx]);
    }
  }

  std::vector<double> state(depth, 0.0);
  std::vector<double> next(depth, 0.0);
  std::vector<double> k1(depth, 0.0);
  std::vector<double> k2(depth, 0.0);
  std::vector<double> k3(depth, 0.0);
  std::vector<double> k4(depth, 0.0);
  std::vector<double> tmp(depth, 0.0);
  auto deriv_from_vals = [&](const double *vals, const std::vector<double> &src,
                             std::vector<double> &dst) {
    std::fill(dst.begin(), dst.end(), 0.0);
    for (std::size_t i = depth; i-- > 0;) {
      const double f_ref = vals[i];
      if (!std::isfinite(f_ref) || f_ref <= 0.0) {
        continue;
      }
      const double s_down = (i + 1 < depth)
                                ? clamp_probability(1.0 - src[i + 1])
                                : vals[depth];
      if (!std::isfinite(s_down) || s_down <= 0.0) {
        continue;
      }
      const double val = f_ref * s_down;
      if (std::isfinite(val) && val > 0.0) {
        dst[i] = val;
      }
    }
  };

  std::vector<double> unique_cdf(unique_times.size(), 0.0);
  for (std::size_t unique_idx = 0; unique_idx < unique_times.size();
       ++unique_idx) {
    const double h = unique_times[unique_idx] / static_cast<double>(n_steps);
    const std::size_t base = unique_idx * grid_points;
    std::fill(state.begin(), state.end(), 0.0);
    std::fill(next.begin(), next.end(), 0.0);
    for (int step = 0; step < n_steps; ++step) {
      const std::size_t g0 = base + static_cast<std::size_t>(2 * step);
      const double *vals_x = grid_slot(g0);
      const double *vals_half = grid_slot(g0 + 1u);
      const double *vals_next = grid_slot(g0 + 2u);
      deriv_from_vals(vals_x, state, k1);
      for (std::size_t j = 0; j < depth; ++j) {
        tmp[j] = state[j] + 0.5 * h * k1[j];
      }
      deriv_from_vals(vals_half, tmp, k2);
      for (std::size_t j = 0; j < depth; ++j) {
        tmp[j] = state[j] + 0.5 * h * k2[j];
      }
      deriv_from_vals(vals_half, tmp, k3);
      for (std::size_t j = 0; j < depth; ++j) {
        tmp[j] = state[j] + h * k3[j];
      }
      deriv_from_vals(vals_next, tmp, k4);
      for (std::size_t j = 0; j < depth; ++j) {
        next[j] = state[j] + (h / 6.0) *
                                 (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
      }
      state.swap(next);
      for (double value : state) {
        if (!std::isfinite(value)) {
          return false;
        }
      }
    }
    unique_cdf[unique_idx] = clamp_probability(state.front());
  }

  for (std::size_t i = 0; i < times.size(); ++i) {
    if (time_to_unique[i] >= 0) {
      cdf_out[i] = unique_cdf[static_cast<std::size_t>(time_to_unique[i])];
    }
  }
  return true;
}

bool guard_cdf_batch_prepared_internal(const GuardEvalInput &input,
                                       const std::vector<double> &times,
                                       std::vector<double> &cdf_out) {
  cdf_out.assign(times.size(), 0.0);
  if (times.empty()) {
    return true;
  }
  PreparedLinearGuardChainInfo scratch;
  const PreparedLinearGuardChainInfo *prepared = nullptr;
  if (!resolve_prepared_linear_guard_chain_info(input, scratch, prepared) ||
      prepared == nullptr) {
    return false;
  }
  if (guard_cdf_simple_pair_batch_prepared(input, *prepared, times, cdf_out)) {
    return true;
  }
  return guard_cdf_general_batch_prepared(input, *prepared, times, cdf_out);
}

std::vector<int> gather_blocker_sources(const uuber::NativeContext &ctx,
                                        int guard_node_idx) {
  std::vector<int> sources;
  const uuber::IrNode &guard = ir_node_required(ctx, guard_node_idx);
  if (guard.blocker_idx < 0) {
    return sources;
  }
  const uuber::IrNode &blocker = ir_node_required(ctx, guard.blocker_idx);
  sources = ensure_source_ids(ctx, blocker);
  if (blocker.op == uuber::IrNodeOp::EventAcc ||
      blocker.op == uuber::IrNodeOp::EventPool) {
    LabelRef blocker_ref = node_label_ref(ctx, blocker);
    int label_idx = blocker_ref.label_id;
    if (label_idx < 0) {
      label_idx = NA_INTEGER;
    }
    sources.push_back(label_idx);
  }
  sort_unique(sources);
  return sources;
}

inline bool resolve_simple_acc_event_batch_info(
    const uuber::NativeContext &ctx, int node_idx,
    const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    int component_idx, SimpleAccEventBatchInfo &out,
    bool require_simple_outcome_support) {
  out = SimpleAccEventBatchInfo{};
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  if (node.op != uuber::IrNodeOp::EventAcc || node.event_idx < 0 ||
      node.event_idx >= static_cast<int>(ctx.ir.events.size()) ||
      (node.flags & (uuber::IR_NODE_FLAG_SPECIAL_DEADLINE |
                     uuber::IR_NODE_FLAG_SPECIAL_GUESS)) != 0u) {
    return false;
  }
  const uuber::IrEvent &event =
      ctx.ir.events[static_cast<std::size_t>(node.event_idx)];
  if (event.acc_idx < 0 || event.pool_idx >= 0 ||
      (require_simple_outcome_support &&
       !simple_direct_outcome_fastpath_supported(ctx, false,
                                                 event.outcome_idx))) {
    return false;
  }

  const int acc_idx = event.acc_idx;
  const uuber::NativeAccumulator &acc =
      ctx.accumulators[static_cast<std::size_t>(acc_idx)];
  const TrialAccumulatorParams *override =
      get_trial_param_entry(trial_params, acc_idx);
  if (!component_active_idx(acc, component_idx, override)) {
    out.acc_idx = acc_idx;
    out.component_ok = false;
    return true;
  }
  const int onset_kind = override ? override->onset_kind : acc.onset_kind;
  if (ctx.has_chained_onsets && onset_kind != uuber::ONSET_ABSOLUTE) {
    return false;
  }

  const uuber::TrialParamsSoA *base_soa =
      trial_params_soa_batch != nullptr && !trial_params_soa_batch->empty()
          ? (*trial_params_soa_batch)[0]
          : resolve_trial_params_soa(ctx, trial_params);
  if (!base_soa) {
    return false;
  }

  double onset = 0.0;
  double q0 = 0.0;
  AccDistParams cfg{};
  resolve_event_numeric_params(acc, acc_idx, override, base_soa, onset, q0, cfg);
  if (trial_params_soa_batch != nullptr) {
    for (const uuber::TrialParamsSoA *point_soa : *trial_params_soa_batch) {
      if (!point_soa || !same_acc_batch_params_except_q(*base_soa, *point_soa,
                                                        acc_idx)) {
        return false;
      }
    }
  }

  out.base_soa = base_soa;
  out.acc_idx = acc_idx;
  out.q0 = q0;
  out.cfg = cfg;
  out.lower_bound = default_lower_bound_transform(cfg);
  out.onset_eff = total_onset_with_t0(onset, cfg);
  out.component_ok = true;
  return true;
}

inline bool eval_simple_acc_event_density_batch_from_info(
    const SimpleAccEventBatchInfo &info, const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &shifted, std::vector<double> &values,
    std::vector<double> &density_out);

inline bool eval_simple_acc_event_survival_batch_from_info(
    const SimpleAccEventBatchInfo &info, const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &shifted, std::vector<double> &values,
    std::vector<double> &survival_out);

struct SimpleAccEventBatchNodeCacheEntry {
  int node_id{-1};
  SimpleAccEventBatchInfo info{};
  std::vector<double> density;
  std::vector<double> survival;
  bool density_ready{false};
  bool survival_ready{false};
};

inline bool eval_simple_acc_event_competing_density_batch_from_info(
    const uuber::NativeContext &ctx,
    const SimpleAccEventBatchInfo &target_info,
    const std::vector<int> &competitor_ids, const std::vector<double> &times,
    int component_idx, const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    int outcome_idx_context, std::vector<double> &density_out);

inline bool eval_simple_acc_event_competing_density_batch_from_info(
    const uuber::NativeContext &ctx,
    const SimpleAccEventBatchInfo &target_info,
    const std::vector<int> &competitor_ids, const std::vector<double> &times,
    int component_idx, const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    int outcome_idx_context, std::vector<double> &density_out) {
  density_out.assign(times.size(), 0.0);
  if (!target_info.component_ok) {
    return true;
  }

  std::vector<double> shifted;
  std::vector<double> values;
  if (!eval_simple_acc_event_density_batch_from_info(
          target_info, times, trial_params_soa_batch, shifted, values,
          density_out)) {
    return false;
  }

  std::vector<SimpleAccEventBatchNodeCacheEntry> node_cache;
  struct SimpleDonorBatchTerms {
    double weight{0.0};
    std::size_t density_index{0u};
    std::vector<std::size_t> survival_indices;
  };
  std::vector<SimpleDonorBatchTerms> donor_terms;
  std::vector<std::size_t> target_survival_indices;
  std::vector<int> filtered_competitors;

  auto find_or_create_entry_index = [&](int node_id) -> std::size_t {
    for (std::size_t idx = 0; idx < node_cache.size(); ++idx) {
      if (node_cache[idx].node_id == node_id) {
        return idx;
      }
    }
    const int node_idx = resolve_dense_node_idx_required(ctx, node_id);
    SimpleAccEventBatchNodeCacheEntry entry;
    entry.node_id = node_id;
    if (!resolve_simple_acc_event_batch_info(
            ctx, node_idx, trial_params, trial_params_soa_batch,
            component_idx, entry.info, false)) {
      return std::numeric_limits<std::size_t>::max();
    }
    node_cache.push_back(std::move(entry));
    return node_cache.size() - 1u;
  };

  auto ensure_density = [&](std::size_t entry_idx) -> bool {
    SimpleAccEventBatchNodeCacheEntry &entry = node_cache[entry_idx];
    if (entry.density_ready) {
      return true;
    }
    entry.density_ready = eval_simple_acc_event_density_batch_from_info(
        entry.info, times, trial_params_soa_batch, shifted, values,
        entry.density);
    return entry.density_ready;
  };

  auto ensure_survival = [&](std::size_t entry_idx) -> bool {
    SimpleAccEventBatchNodeCacheEntry &entry = node_cache[entry_idx];
    if (entry.survival_ready) {
      return true;
    }
    entry.survival_ready = eval_simple_acc_event_survival_batch_from_info(
        entry.info, times, trial_params_soa_batch, shifted, values,
        entry.survival);
    return entry.survival_ready;
  };

  if (outcome_idx_context >= 0 &&
      outcome_idx_context < static_cast<int>(ctx.outcome_info.size())) {
    const uuber::OutcomeContextInfo &target_outcome =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx_context)];
    if (!target_outcome.alias_sources.empty()) {
      return false;
    }
    for (const auto &donor : target_outcome.guess_donors) {
      if (donor.rt_policy != "keep" || donor.weight == 0.0 ||
          donor.outcome_idx < 0 ||
          donor.outcome_idx >= static_cast<int>(ctx.outcome_info.size()) ||
          donor.outcome_idx == outcome_idx_context ||
          !shared_trigger_mask_batch_supported(ctx, false, donor.outcome_idx)) {
        return false;
      }
      const uuber::OutcomeContextInfo &donor_outcome =
          ctx.outcome_info[static_cast<std::size_t>(donor.outcome_idx)];
      const std::size_t donor_density_index =
          find_or_create_entry_index(donor_outcome.node_id);
      if (donor_density_index == std::numeric_limits<std::size_t>::max() ||
          !ensure_density(donor_density_index)) {
        return false;
      }
      const std::vector<int> &donor_competitors = filter_competitor_ids(
          ctx, donor_outcome.competitor_ids, component_idx,
          filtered_competitors);
      SimpleDonorBatchTerms donor_term;
      donor_term.weight = donor.weight;
      donor_term.density_index = donor_density_index;
      donor_term.survival_indices.reserve(donor_competitors.size());
      for (int donor_competitor_id : donor_competitors) {
        const std::size_t comp_index =
            find_or_create_entry_index(donor_competitor_id);
        if (comp_index == std::numeric_limits<std::size_t>::max() ||
            !ensure_survival(comp_index)) {
          return false;
        }
        donor_term.survival_indices.push_back(comp_index);
      }
      donor_terms.push_back(std::move(donor_term));
    }
  }

  for (int competitor_id : competitor_ids) {
    const std::size_t comp_index = find_or_create_entry_index(competitor_id);
    if (comp_index == std::numeric_limits<std::size_t>::max() ||
        !ensure_survival(comp_index)) {
      return false;
    }
    target_survival_indices.push_back(comp_index);
  }

  if (!donor_terms.empty() || !target_survival_indices.empty()) {
    for (std::size_t i = 0; i < density_out.size(); ++i) {
      double total_density = density_out[i];
      for (const SimpleDonorBatchTerms &donor_term : donor_terms) {
        double donor_value =
            node_cache[donor_term.density_index].density[i];
        if (!(donor_value > 0.0)) {
          continue;
        }
        for (std::size_t survival_index : donor_term.survival_indices) {
          const double surv =
              clamp_probability(node_cache[survival_index].survival[i]);
          if (!std::isfinite(surv) || surv <= 0.0) {
            donor_value = 0.0;
            break;
          }
          donor_value = safe_density(donor_value * surv);
        }
        if (donor_value > 0.0) {
          total_density =
              safe_density(total_density + donor_term.weight * donor_value);
        }
      }
      if (!(total_density > 0.0)) {
        density_out[i] = 0.0;
        continue;
      }
      for (std::size_t survival_index : target_survival_indices) {
        const double surv =
            clamp_probability(node_cache[survival_index].survival[i]);
        if (!std::isfinite(surv) || surv <= 0.0) {
          total_density = 0.0;
          break;
        }
        total_density = safe_density(total_density * surv);
      }
      density_out[i] = total_density > 0.0 ? total_density : 0.0;
    }
  }
  return true;
}

inline double simple_acc_event_q_for_batch_point(
    const SimpleAccEventBatchInfo &info,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::size_t point_idx) {
  const uuber::TrialParamsSoA *point_soa =
      trial_params_soa_batch != nullptr ? (*trial_params_soa_batch)[point_idx]
                                        : info.base_soa;
  if (point_soa && point_soa->valid && info.acc_idx >= 0 &&
      info.acc_idx < point_soa->n_acc &&
      static_cast<std::size_t>(info.acc_idx) < point_soa->q.size()) {
    return point_soa->q[static_cast<std::size_t>(info.acc_idx)];
  }
  return info.q0;
}

inline bool try_eval_conditioned_simple_acc_event_batch(
    const NodeEvalState &state, const uuber::IrEvent &event,
    const std::vector<double> &times, EvalNeed need,
    uuber::TreeNodeBatchValues &out_values) {
  out_values = uuber::TreeNodeBatchValues{};
  const std::size_t point_count = times.size();
  out_values.density.assign(point_count, 0.0);
  out_values.survival.assign(point_count, 1.0);
  out_values.cdf.assign(point_count, 0.0);
  if (point_count == 0u || event.node_idx < 0) {
    return true;
  }

  const int outcome_idx_context =
      state.outcome_idx >= 0 ? state.outcome_idx : event.outcome_idx;
  if (needs_density(need) &&
      !simple_direct_outcome_fastpath_supported(
          state.ctx, state.include_na_donors, outcome_idx_context)) {
    return false;
  }

  SimpleAccEventBatchInfo info;
  if (!resolve_simple_acc_event_batch_info(
          state.ctx, event.node_idx, state.trial_params,
          state.trial_params_soa_batch, state.component_idx, info, false)) {
    return false;
  }

  const int label_id = event.label_id;
  const SourceTimeConstraint *constraint =
      time_constraints_find(&state.time_constraints, label_id);
  const bool self_is_time_conditioned =
      constraint != nullptr && !time_constraint_empty(*constraint);
  const bool forced_complete =
      forced_state_contains_complete(state.forced_state, label_id);
  const bool forced_survive =
      forced_state_contains_survive(state.forced_state, label_id);

  if (!self_is_time_conditioned && forced_complete) {
    if (needs_cdf(need)) {
      std::fill(out_values.cdf.begin(), out_values.cdf.end(), 1.0);
    }
    if (needs_survival(need)) {
      std::fill(out_values.survival.begin(), out_values.survival.end(), 0.0);
    }
    return true;
  }

  if (forced_survive) {
    if (self_is_time_conditioned) {
      if (needs_survival(need)) {
        std::fill(out_values.survival.begin(), out_values.survival.end(), 0.0);
      }
      if (needs_cdf(need)) {
        std::fill(out_values.cdf.begin(), out_values.cdf.end(), 0.0);
      }
    }
    return true;
  }

  if (!info.component_ok) {
    return true;
  }

  std::vector<double> shifted;
  std::vector<double> values;
  std::vector<double> base_density;
  std::vector<double> base_survival;
  if (needs_density(need) &&
      !eval_simple_acc_event_density_batch_from_info(
          info, times, state.trial_params_soa_batch, shifted, values,
          base_density)) {
    return false;
  }
  if ((needs_cdf(need) || needs_survival(need) || self_is_time_conditioned) &&
      !eval_simple_acc_event_survival_batch_from_info(
          info, times, state.trial_params_soa_batch, shifted, values,
          base_survival)) {
    return false;
  }

  if (!self_is_time_conditioned) {
    if (needs_density(need)) {
      out_values.density = std::move(base_density);
      for (double &value : out_values.density) {
        value = safe_density(value);
      }
    }
    if (needs_survival(need) || needs_cdf(need)) {
      for (std::size_t i = 0; i < point_count; ++i) {
        const double survival = clamp_probability(base_survival[i]);
        if (needs_survival(need)) {
          out_values.survival[i] = survival;
        }
        if (needs_cdf(need)) {
          out_values.cdf[i] = clamp_probability(1.0 - survival);
        }
      }
    }
    return true;
  }

  double bound_lower = 0.0;
  double bound_upper = std::numeric_limits<double>::infinity();
  bool has_exact = false;
  double exact_time = std::numeric_limits<double>::quiet_NaN();
  if (constraint->has_exact) {
    has_exact = true;
    exact_time = constraint->exact_time;
  }
  if (constraint->has_lower) {
    bound_lower = constraint->lower;
  }
  if (constraint->has_upper) {
    bound_upper = constraint->upper;
  }

  if (has_exact) {
    if ((constraint->has_lower &&
         (!(exact_time > bound_lower) &&
          !time_constraint_same_time(exact_time, bound_lower))) ||
        (constraint->has_upper && exact_time > bound_upper &&
         !time_constraint_same_time(exact_time, bound_upper))) {
      if (needs_survival(need)) {
        std::fill(out_values.survival.begin(), out_values.survival.end(), 0.0);
      }
      if (needs_cdf(need)) {
        std::fill(out_values.cdf.begin(), out_values.cdf.end(), 0.0);
      }
      return true;
    }
    for (std::size_t i = 0; i < point_count; ++i) {
      const bool before_exact =
          std::isfinite(times[i]) && times[i] < exact_time &&
          !time_constraint_same_time(times[i], exact_time);
      const double cdf_exact = before_exact ? 0.0 : 1.0;
      if (needs_cdf(need)) {
        out_values.cdf[i] = cdf_exact;
      }
      if (needs_survival(need)) {
        out_values.survival[i] = 1.0 - cdf_exact;
      }
    }
    return true;
  }

  const double lower_base_cdf =
      std::isfinite(bound_lower)
          ? clamp_probability(eval_cdf_single_with_lower_bound(
                info.cfg, bound_lower - info.onset_eff, info.lower_bound))
          : 0.0;
  const double upper_base_cdf =
      std::isfinite(bound_upper)
          ? clamp_probability(eval_cdf_single_with_lower_bound(
                info.cfg, bound_upper - info.onset_eff, info.lower_bound))
          : 1.0;

  for (std::size_t i = 0; i < point_count; ++i) {
    const double q =
        simple_acc_event_q_for_batch_point(info, state.trial_params_soa_batch, i);
    const double success_prob = clamp_probability(1.0 - q);
    const double lower_cdf = clamp_probability(success_prob * lower_base_cdf);
    const double upper_cdf = clamp_probability(success_prob * upper_base_cdf);
    const double condition_mass =
        clamp_probability(upper_cdf - lower_cdf);
    if (!(bound_upper > bound_lower) || !(condition_mass > 0.0)) {
      if (needs_survival(need)) {
        out_values.survival[i] = 0.0;
      }
      if (needs_cdf(need)) {
        out_values.cdf[i] = 0.0;
      }
      continue;
    }

    const double base_cdf =
        clamp_probability(1.0 - clamp_probability(base_survival[i]));
    if (needs_density(need)) {
      if (std::isfinite(times[i]) && times[i] > bound_lower &&
          times[i] <= bound_upper) {
        out_values.density[i] =
            safe_density(safe_density(base_density[i]) / condition_mass);
      } else {
        out_values.density[i] = 0.0;
      }
    }
    double cdf_cond = 0.0;
    if (!std::isfinite(times[i]) || times[i] >= bound_upper) {
      cdf_cond = 1.0;
    } else if (!(times[i] > bound_lower)) {
      cdf_cond = 0.0;
    } else {
      cdf_cond = clamp_probability((base_cdf - lower_cdf) / condition_mass);
    }
    if (needs_cdf(need)) {
      out_values.cdf[i] = cdf_cond;
    }
    if (needs_survival(need)) {
      out_values.survival[i] = clamp_probability(1.0 - cdf_cond);
    }
  }

  return true;
}

inline bool eval_simple_acc_event_density_batch_from_info(
    const SimpleAccEventBatchInfo &info, const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &shifted, std::vector<double> &values,
    std::vector<double> &density_out) {
  density_out.assign(times.size(), 0.0);
  if (!info.component_ok) {
    return true;
  }

  shifted.assign(times.size(), 0.0);
  values.assign(times.size(), 0.0);
  for (std::size_t i = 0; i < times.size(); ++i) {
    shifted[i] = times[i] - info.onset_eff;
  }
  eval_pdf_vec_with_lower_bound(info.cfg, info.lower_bound, shifted.data(),
                                times.size(), values.data());

  for (std::size_t i = 0; i < times.size(); ++i) {
    const double q =
        simple_acc_event_q_for_batch_point(info, trial_params_soa_batch, i);
    const double success_prob = 1.0 - q;
    if (!std::isfinite(times[i]) || times[i] < 0.0 ||
        times[i] < info.onset_eff ||
        !(success_prob > 0.0)) {
      continue;
    }
    density_out[i] = safe_density(success_prob * safe_density(values[i]));
  }
  return true;
}

inline bool eval_simple_acc_event_survival_batch_from_info(
    const SimpleAccEventBatchInfo &info, const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &shifted, std::vector<double> &values,
    std::vector<double> &survival_out) {
  survival_out.assign(times.size(), 1.0);
  if (!info.component_ok) {
    return true;
  }

  shifted.assign(times.size(), 0.0);
  values.assign(times.size(), 0.0);
  for (std::size_t i = 0; i < times.size(); ++i) {
    shifted[i] = times[i] - info.onset_eff;
  }
  eval_cdf_vec_with_lower_bound(info.cfg, info.lower_bound, shifted.data(),
                                times.size(), values.data());

  for (std::size_t i = 0; i < times.size(); ++i) {
    if (!std::isfinite(times[i]) || times[i] < 0.0 ||
        times[i] < info.onset_eff) {
      continue;
    }
    const double q =
        simple_acc_event_q_for_batch_point(info, trial_params_soa_batch, i);
    const double success_prob = clamp(1.0 - q, 0.0, 1.0);
    const double cdf = clamp_probability(values[i]);
    survival_out[i] = clamp_probability(q + success_prob * (1.0 - cdf));
  }
  return true;
}

inline bool try_evaluate_simple_overlap_event_batch_values(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    uuber::TreeNodeBatchValues &out_values) {
  SimpleAccEventBatchInfo info;
  if (!resolve_simple_acc_event_batch_info(
          ctx, node_idx, trial_params, trial_params_soa_batch, component_idx,
          info, false)) {
    return false;
  }
  std::vector<double> shifted;
  std::vector<double> values;
  if (!eval_simple_acc_event_density_batch_from_info(
          info, times, trial_params_soa_batch, shifted, values,
          out_values.density) ||
      !eval_simple_acc_event_survival_batch_from_info(
          info, times, trial_params_soa_batch, shifted, values,
          out_values.survival)) {
    return false;
  }
  out_values.cdf.assign(times.size(), 0.0);
  for (std::size_t i = 0; i < times.size(); ++i) {
    out_values.density[i] = safe_density(out_values.density[i]);
    out_values.survival[i] = clamp_probability(out_values.survival[i]);
    out_values.cdf[i] = clamp_probability(1.0 - out_values.survival[i]);
  }
  return true;
}

inline bool eval_simple_acc_event_density_batch(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    int outcome_idx_context, std::vector<double> &density_out) {
  SimpleAccEventBatchInfo info;
  if (!resolve_simple_acc_event_batch_info(
          ctx, node_idx, trial_params, trial_params_soa_batch, component_idx,
          info)) {
    return false;
  }
  const std::vector<int> no_competitors;
  if (!eval_simple_acc_event_competing_density_batch_from_info(
          ctx, info, no_competitors, times, component_idx, trial_params,
          trial_params_soa_batch, outcome_idx_context, density_out)) {
    return false;
  }
  return true;
}

inline bool eval_simple_acc_event_competing_density_batch(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, const std::vector<double> &times,
    int component_idx, const TrialParamSet *trial_params,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    int outcome_idx_context, std::vector<double> &density_out) {
  SimpleAccEventBatchInfo target_info;
  if (!resolve_simple_acc_event_batch_info(
          ctx, node_idx, trial_params, trial_params_soa_batch, component_idx,
          target_info)) {
    return false;
  }
  return eval_simple_acc_event_competing_density_batch_from_info(
      ctx, target_info, competitor_ids, times, component_idx, trial_params,
      trial_params_soa_batch, outcome_idx_context, density_out);
}

bool node_density_with_competitors_batch_internal(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<double> &times, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context, std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    DensityBatchWorkspace *workspace) {
  density_out.assign(times.size(), 0.0);
  if (times.empty()) {
    return true;
  }
  if (trial_params_soa_batch != nullptr &&
      trial_params_soa_batch->size() != times.size()) {
    return false;
  }

  NodeEvalState state(
      ctx, 0.0, component_idx, trial_params, trial_type_key, include_na_donors,
      outcome_idx_context, time_constraints, nullptr,
      forced_complete_bits_valid ? forced_complete_bits : nullptr,
      forced_complete_bits_valid,
      forced_survive_bits_valid ? forced_survive_bits : nullptr,
      forced_survive_bits_valid);
  state.trial_params_soa_batch = trial_params_soa_batch;

  const uuber::BitsetState forced_complete_seed = state.forced_complete_bits;
  const uuber::BitsetState forced_survive_seed = state.forced_survive_bits;
  const bool forced_complete_seed_valid = state.forced_complete_bits_valid;
  const bool forced_survive_seed_valid = state.forced_survive_bits_valid;

  const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr;
  if (!competitor_ids.empty()) {
    competitor_cache = &fetch_competitor_cluster_cache(ctx, competitor_ids);
  }
  const int node_idx = resolve_dense_node_idx_required(ctx, node_id);
  const uuber::IrOutcomeCouplingOp *coupling_spec =
      find_outcome_coupling_spec(ctx, node_idx, competitor_ids);
  const bool coupling_forced_empty =
      !forced_bits_any(forced_complete_bits, forced_complete_bits_valid) &&
      !forced_bits_any(forced_survive_bits, forced_survive_bits_valid) &&
      !time_constraints_any(time_constraints);
  const bool allow_conditioned_specialized_overlap =
      coupling_spec != nullptr &&
      coupling_spec->kind == uuber::IrOutcomeCouplingKind::SharedBlockerPair;
  if (!competitor_ids.empty() &&
      (coupling_forced_empty || allow_conditioned_specialized_overlap) &&
      evaluate_specialized_overlap_density_batch(
          ctx, coupling_spec, times, component_idx, trial_params,
          trial_type_key, trial_params_soa_batch, density_out,
          allow_conditioned_specialized_overlap ? time_constraints : nullptr,
          allow_conditioned_specialized_overlap ? &state.forced_state
                                               : nullptr)) {
    return true;
  }
  if (should_use_exact_density_program(ctx, node_idx, competitor_ids, state,
                                       coupling_spec)) {
    return exact_outcome_density_batch_from_state(
        ctx, node_idx, competitor_ids, state, times, density_out,
        competitor_cache);
  }
  const bool has_competitors =
      competitor_cache && !competitor_cache->compiled_ops.empty();
  const bool simple_direct_fastpath_base_eligible = coupling_forced_empty;
  if (simple_direct_fastpath_base_eligible) {
    const bool plain_simple_outcome_supported =
        simple_direct_outcome_fastpath_supported(ctx, include_na_donors,
                                                 outcome_idx_context);
    if (has_competitors &&
        plain_simple_outcome_supported &&
        eval_simple_acc_event_competing_density_batch(
            ctx, node_idx, competitor_ids, times, component_idx, trial_params,
            trial_params_soa_batch, outcome_idx_context, density_out)) {
      return true;
    }
    if (!has_competitors && plain_simple_outcome_supported &&
        eval_simple_acc_event_density_batch(
            ctx, node_idx, times, component_idx, trial_params,
            trial_params_soa_batch, outcome_idx_context, density_out)) {
      return true;
    }
  }
  std::vector<double> eval_times(times);
  std::vector<std::uint8_t> valid_points(times.size(), 0);
  for (std::size_t i = 0; i < times.size(); ++i) {
    const double t = times[i];
    if (std::isfinite(t) && t >= 0.0) {
      valid_points[i] = 1;
    } else {
      eval_times[i] = 0.0;
    }
  }
  state.include_na_donors = include_na_donors;
  state.outcome_idx = outcome_idx_context;
  state.forced_complete_bits = forced_complete_seed;
  state.forced_survive_bits = forced_survive_seed;
  state.forced_complete_bits_valid = forced_complete_seed_valid;
  state.forced_survive_bits_valid = forced_survive_seed_valid;
  uuber::TreeNodeBatchValues base_values;
  if (!eval_node_batch_with_state_dense(node_idx, eval_times, state,
                                        EvalNeed::kDensity, base_values) ||
      base_values.density.size() != times.size()) {
    Rcpp::stop("IR tree vector batch execution failed for node density");
  }
  for (std::size_t i = 0; i < density_out.size(); ++i) {
    if (!valid_points[i]) {
      continue;
    }
    density_out[i] = safe_density(base_values.density[i]);
  }

  if (has_competitors) {
    state.include_na_donors = false;
    state.outcome_idx = -1;
    state.forced_complete_bits = forced_complete_seed;
    state.forced_survive_bits = forced_survive_seed;
    state.forced_complete_bits_valid = forced_complete_seed_valid;
    state.forced_survive_bits_valid = forced_survive_seed_valid;
    std::vector<double> survival_product;
    competitor_survival_batch_from_state_compiled_ops(
        ctx, competitor_ids, state, eval_times, survival_product,
        competitor_cache);
    if (survival_product.size() != times.size()) {
      Rcpp::stop("IR competitor batch output size mismatch");
    }
    for (std::size_t i = 0; i < density_out.size(); ++i) {
      if (!valid_points[i] || density_out[i] <= 0.0) {
        density_out[i] = 0.0;
        continue;
      }
      const double surv = clamp_probability(survival_product[i]);
      if (!std::isfinite(surv) || surv <= 0.0) {
        density_out[i] = 0.0;
        continue;
      }
      density_out[i] = safe_density(density_out[i] * surv);
    }
  }

  return true;
}

inline bool node_density_entry_batch_idx(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<double> &times, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan, bool use_shared_trigger_eval,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints,
    const std::vector<const uuber::TrialParamsSoA *>
        *trial_params_soa_batch_override,
    const SharedTriggerMaskSoABatch *shared_trigger_mask_batch_override,
    DensityBatchWorkspace *workspace) {
  if (trial_params_soa_batch_override != nullptr) {
    return node_density_with_competitors_batch_internal(
        ctx, node_id, times, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, density_out, time_constraints,
        trial_params_soa_batch_override, workspace);
  }
  if (!use_shared_trigger_eval || !trial_params) {
    return node_density_with_competitors_batch_internal(
        ctx, node_id, times, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, density_out, time_constraints,
        nullptr, workspace);
  }

  SharedTriggerPlan local_plan;
  const SharedTriggerPlan *plan_ptr = trigger_plan;
  if (!plan_ptr || shared_trigger_count(*plan_ptr) == 0u) {
    local_plan = build_shared_trigger_plan(ctx, trial_params);
    plan_ptr = &local_plan;
  }
  if (shared_trigger_count(*plan_ptr) == 0u) {
    return node_density_with_competitors_batch_internal(
        ctx, node_id, times, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, density_out, time_constraints,
        nullptr, workspace);
  }
  if (!shared_trigger_density_batch_supported(ctx, include_na_donors,
                                              outcome_idx_context)) {
    return false;
  }
  SharedTriggerMaskSoABatch local_mask_batch;
  const SharedTriggerMaskSoABatch *mask_batch_ptr =
      shared_trigger_mask_batch_override;
  if (mask_batch_ptr == nullptr || mask_batch_ptr->mask_param_ptrs.empty()) {
    if (!build_shared_trigger_mask_soa_batch(ctx, trial_params, *plan_ptr,
                                             local_mask_batch) ||
        local_mask_batch.mask_param_ptrs.empty()) {
      return false;
    }
    mask_batch_ptr = &local_mask_batch;
  }
  std::vector<double> local_expanded_times;
  std::vector<const uuber::TrialParamsSoA *> local_expanded_params_soa;
  std::vector<double> local_compressed_times;
  std::vector<const uuber::TrialParamsSoA *> local_compressed_params_soa;
  std::vector<std::size_t> local_expanded_to_compressed;
  std::vector<double> local_batch_density;
  std::vector<double> &expanded_times =
      workspace ? workspace->expanded_times : local_expanded_times;
  std::vector<const uuber::TrialParamsSoA *> &expanded_params_soa =
      workspace ? workspace->expanded_params_soa : local_expanded_params_soa;
  std::vector<double> &compressed_times =
      workspace ? workspace->compressed_times : local_compressed_times;
  std::vector<const uuber::TrialParamsSoA *> &compressed_params_soa =
      workspace ? workspace->compressed_params_soa : local_compressed_params_soa;
  std::vector<std::size_t> &expanded_to_compressed =
      workspace ? workspace->expanded_to_compressed : local_expanded_to_compressed;
  std::vector<double> &batch_density =
      workspace ? workspace->batch_density : local_batch_density;
  expanded_times.clear();
  expanded_params_soa.clear();
  expanded_times.reserve(times.size() * mask_batch_ptr->mask_param_ptrs.size());
  expanded_params_soa.reserve(times.size() *
                              mask_batch_ptr->mask_param_ptrs.size());
  for (std::size_t mask_idx = 0; mask_idx < mask_batch_ptr->mask_param_ptrs.size();
       ++mask_idx) {
    const uuber::TrialParamsSoA *mask_params_soa =
        mask_batch_ptr->mask_param_ptrs[mask_idx];
    for (double t : times) {
      expanded_times.push_back(t);
      expanded_params_soa.push_back(mask_params_soa);
    }
  }

  const bool can_compress = times_have_duplicates(times);
  const std::vector<double> *batch_times_ptr = &expanded_times;
  const std::vector<const uuber::TrialParamsSoA *> *batch_params_ptr =
      &expanded_params_soa;
  if (can_compress) {
    compressed_times.clear();
    compressed_params_soa.clear();
    expanded_to_compressed.clear();
    compress_time_trial_params_soa_batch(
        expanded_times, expanded_params_soa, compressed_times,
        compressed_params_soa, expanded_to_compressed);
    batch_times_ptr = &compressed_times;
    batch_params_ptr = &compressed_params_soa;
  }

  batch_density.clear();
  const bool batched = node_density_with_competitors_batch_internal(
      ctx, node_id, *batch_times_ptr, component_idx, forced_complete_bits,
      forced_complete_bits_valid, forced_survive_bits,
      forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
      include_na_donors, outcome_idx_context, batch_density,
      time_constraints, batch_params_ptr, workspace);
  if (!batched || batch_density.size() != batch_times_ptr->size()) {
    return false;
  }

  density_out.assign(times.size(), 0.0);
  std::size_t offset = 0u;
  for (std::size_t mask_idx = 0; mask_idx < mask_batch_ptr->mask_weights.size();
       ++mask_idx) {
    const double weight = mask_batch_ptr->mask_weights[mask_idx];
    if (!(std::isfinite(weight) && weight > 0.0)) {
      offset += times.size();
      continue;
    }
    for (std::size_t time_idx = 0; time_idx < times.size(); ++time_idx) {
      const double density = can_compress
                                 ? batch_density[expanded_to_compressed[offset +
                                                                       time_idx]]
                                 : batch_density[offset + time_idx];
      if (std::isfinite(density) && density > 0.0) {
        density_out[time_idx] += weight * density;
      }
    }
    offset += times.size();
  }
  for (double &density : density_out) {
    density = safe_density(density);
  }
  return true;
}

inline double normalize_outcome_probability(double value);

inline double normalize_outcome_probability(double value) {
  const double clamped = clamp_probability(value);
  return clamped <= 1e-11 ? 0.0 : clamped;
}

inline bool integrate_outcome_probability_from_density_batch_idx(
    const uuber::NativeContext &ctx, int node_id, double upper,
    int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, double rel_tol, double abs_tol,
    const std::vector<const TrialParamSet *> &trial_params_batch,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context, std::vector<double> &out_probabilities,
    const std::vector<const uuber::TrialParamsSoA *>
        *trial_params_soa_batch_override = nullptr,
    const std::vector<const PreparedTrialParamsRuntime *>
        *prepared_runtime_batch = nullptr,
    DensityBatchWorkspace *workspace = nullptr) {
  const std::size_t point_count =
      trial_params_soa_batch_override ? trial_params_soa_batch_override->size()
      : (prepared_runtime_batch ? prepared_runtime_batch->size()
                                : trial_params_batch.size());
  out_probabilities.assign(point_count, 0.0);
  if (point_count == 0u || !(upper > 0.0)) {
    return point_count == 0u || upper <= 0.0;
  }

  std::vector<const uuber::TrialParamsSoA *> trial_params_soa_storage;
  std::vector<const TrialParamSet *> local_unique_trial_params;
  std::vector<const uuber::TrialParamsSoA *> local_unique_params_soa;
  std::vector<std::size_t> local_point_to_eval_index;
  std::vector<const TrialParamSet *> &unique_trial_params =
      workspace ? workspace->unique_trial_params : local_unique_trial_params;
  std::vector<const uuber::TrialParamsSoA *> &unique_params_soa =
      workspace ? workspace->unique_params_soa : local_unique_params_soa;
  std::vector<std::size_t> &point_to_eval_index =
      workspace ? workspace->point_to_eval_index : local_point_to_eval_index;
  point_to_eval_index.assign(point_count, 0u);

  const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch_ptr =
      trial_params_soa_batch_override;
  const std::vector<const TrialParamSet *> *trial_params_eval_ptr =
      &trial_params_batch;
  if (trial_params_soa_batch_ptr == nullptr && prepared_runtime_batch != nullptr) {
    if (prepared_runtime_batch->size() != point_count) {
      return false;
    }
    unique_trial_params.clear();
    unique_params_soa.clear();
    unique_trial_params.reserve(point_count);
    unique_params_soa.reserve(point_count);
    std::unordered_map<uuber::NAMapCacheKey, std::size_t, uuber::NAMapCacheKeyHash>
        runtime_index_by_key;
    runtime_index_by_key.reserve(point_count);
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      PreparedTrialParamsRuntime *runtime =
          const_cast<PreparedTrialParamsRuntime *>(
              (*prepared_runtime_batch)[point_idx]);
      if (runtime == nullptr || !ensure_prepared_trial_params_soa(ctx, *runtime)) {
        return false;
      }
      const auto insert_result = runtime_index_by_key.emplace(
          runtime->value_key, unique_params_soa.size());
      if (insert_result.second) {
        unique_trial_params.push_back(runtime->params);
        unique_params_soa.push_back(runtime->soa);
      }
      point_to_eval_index[point_idx] = insert_result.first->second;
    }
    trial_params_eval_ptr = &unique_trial_params;
    trial_params_soa_batch_ptr = &unique_params_soa;
  } else if (trial_params_soa_batch_ptr == nullptr) {
    trial_params_soa_storage.reserve(point_count);
    for (const TrialParamSet *params_ptr : trial_params_batch) {
      const uuber::TrialParamsSoA *trial_params_soa =
          resolve_trial_params_soa(ctx, params_ptr);
      if (!trial_params_soa || !trial_params_soa->valid) {
        return false;
      }
      trial_params_soa_storage.push_back(trial_params_soa);
    }
    trial_params_soa_batch_ptr = &trial_params_soa_storage;
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      point_to_eval_index[point_idx] = point_idx;
    }
  } else if (trial_params_soa_batch_ptr->size() != point_count) {
    return false;
  } else {
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      point_to_eval_index[point_idx] = point_idx;
    }
  }

  for (const uuber::TrialParamsSoA *trial_params_soa :
       *trial_params_soa_batch_ptr) {
    if (!trial_params_soa || !trial_params_soa->valid) {
      return false;
    }
  }

  TrialParamSet base_params_holder;
  const TrialParamSet *rep_params =
      trial_params_eval_ptr->empty() ? nullptr : trial_params_eval_ptr->front();
  if (!rep_params) {
    base_params_holder = build_base_paramset(ctx);
    rep_params = &base_params_holder;
  }

  const std::size_t eval_point_count = trial_params_soa_batch_ptr->size();
  SharedTriggerPlan local_trigger_plan;
  const SharedTriggerPlan *plan_ptr = nullptr;
  bool use_shared_trigger_eval_local = false;
  if (eval_point_count == 1u) {
    local_trigger_plan = build_shared_trigger_plan(ctx, rep_params);
    plan_ptr = &local_trigger_plan;
    use_shared_trigger_eval_local = shared_trigger_count(*plan_ptr) > 0u;
  }
  constexpr int kIntegralSegments = 16;
  constexpr int kInfiniteMaxSegments = 32;
  const bool finite_upper = std::isfinite(upper);
  std::vector<double> upper_bounds(eval_point_count, upper);
  std::vector<TimeIntegralPointPlan> local_point_plans;
  std::vector<double> local_query_times;
  std::vector<double> local_query_weights;
  std::vector<const uuber::TrialParamsSoA *> local_query_params_soa;
  std::vector<std::uint8_t> local_active_points;
  std::vector<std::uint8_t> local_tail_saw_positive;
  std::vector<int> local_tail_small_segments;
  std::vector<double> local_density_out;
  std::vector<double> local_segment_probabilities;
  std::vector<double> local_point_probabilities;
  std::vector<TimeIntegralPointPlan> &point_plans =
      workspace ? workspace->point_plans : local_point_plans;
  std::vector<double> &query_times =
      workspace ? workspace->query_times : local_query_times;
  std::vector<double> &query_weights =
      workspace ? workspace->query_weights : local_query_weights;
  std::vector<const uuber::TrialParamsSoA *> &query_params_soa =
      workspace ? workspace->query_params_soa : local_query_params_soa;
  std::vector<std::uint8_t> &active_points =
      workspace ? workspace->active_points : local_active_points;
  std::vector<std::uint8_t> &tail_saw_positive =
      workspace ? workspace->tail_saw_positive : local_tail_saw_positive;
  std::vector<int> &tail_small_segments =
      workspace ? workspace->tail_small_segments : local_tail_small_segments;
  std::vector<double> &density_out =
      workspace ? workspace->batch_density : local_density_out;
  std::vector<double> &segment_probabilities =
      workspace ? workspace->segment_probabilities : local_segment_probabilities;
  std::vector<double> &point_probabilities =
      workspace ? workspace->point_probabilities : local_point_probabilities;
  point_probabilities.assign(eval_point_count, 0.0);
  active_points.assign(eval_point_count, 1u);
  tail_saw_positive.assign(eval_point_count, 0u);
  tail_small_segments.assign(eval_point_count, 0);

  const auto eval_segment = [&](double lower, double upper_window) -> bool {
    if (!(std::isfinite(lower) && std::isfinite(upper_window)) ||
        !(upper_window > lower)) {
      segment_probabilities.assign(eval_point_count, 0.0);
      return true;
    }
    if (!build_time_integral_query_plan(
            upper_bounds, lower, upper_window, kIntegralSegments, point_plans,
            query_times, query_weights, trial_params_soa_batch_ptr,
            &query_params_soa, &active_points)) {
      return false;
    }
    if (query_times.empty()) {
      segment_probabilities.assign(eval_point_count, 0.0);
      return true;
    }
    density_out.clear();
    const std::vector<const uuber::TrialParamsSoA *> *query_params_soa_arg =
        use_shared_trigger_eval_local ? nullptr : &query_params_soa;
    const bool ok = node_density_entry_batch_idx(
        ctx, node_id, query_times, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, rep_params, trial_type_key,
        include_na_donors, outcome_idx_context, plan_ptr,
        use_shared_trigger_eval_local, density_out, nullptr,
        query_params_soa_arg, nullptr, workspace);
    if (!ok) {
      return false;
    }
    return integrate_query_plan_density_per_point(point_plans, query_weights,
                                                  density_out,
                                                  segment_probabilities);
  };

  if (finite_upper) {
    if (!eval_segment(0.0, upper)) {
      return false;
    }
    for (std::size_t i = 0; i < eval_point_count; ++i) {
      point_probabilities[i] =
          normalize_outcome_probability(segment_probabilities[i]);
    }
  } else {
    double lo = 0.0;
    double hi = 0.5;
    for (int seg = 0; seg < kInfiniteMaxSegments; ++seg) {
      if (!eval_segment(lo, hi)) {
        return false;
      }
      bool any_active = false;
      for (std::size_t point_idx = 0; point_idx < eval_point_count; ++point_idx) {
        if (point_plans[point_idx].count == 0u) {
          active_points[point_idx] = 0u;
          continue;
        }
        const double piece = segment_probabilities[point_idx];
        if (piece > 0.0) {
          tail_saw_positive[point_idx] = 1u;
        }
        point_probabilities[point_idx] += piece;
        const double tail_tol =
            std::max(std::max(0.0, abs_tol),
                     std::max(0.0, rel_tol) *
                         std::max(1.0, std::fabs(point_probabilities[point_idx])));
        if (piece <= tail_tol) {
          ++tail_small_segments[point_idx];
        } else {
          tail_small_segments[point_idx] = 0;
        }
        const bool done =
            tail_saw_positive[point_idx] != 0u && hi >= 4.0 &&
            tail_small_segments[point_idx] >= 3;
        active_points[point_idx] = done ? 0u : 1u;
        any_active = any_active || active_points[point_idx] != 0u;
      }
      if (!any_active) {
        break;
      }
      lo = hi;
      hi *= 2.0;
    }

    for (double &probability : point_probabilities) {
      probability = normalize_outcome_probability(probability);
    }
  }

  for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
    const std::size_t eval_idx = point_to_eval_index[point_idx];
    if (eval_idx >= point_probabilities.size()) {
      return false;
    }
    out_probabilities[point_idx] = point_probabilities[eval_idx];
  }
  return true;
}

inline bool integrate_outcome_probabilities_from_density_multi_idx(
    const uuber::NativeContext &ctx, const std::vector<int> &node_ids,
    double upper, int component_idx, const std::vector<int> &competitor_ids,
    double rel_tol, double abs_tol,
    const std::vector<const TrialParamSet *> &trial_params_batch,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context, std::vector<double> &out_probabilities_matrix,
    const std::vector<const uuber::TrialParamsSoA *>
        *trial_params_soa_batch_override = nullptr,
    const std::vector<const PreparedTrialParamsRuntime *>
        *prepared_runtime_batch = nullptr,
    DensityBatchWorkspace *workspace = nullptr) {
  const std::size_t node_count = node_ids.size();
  const std::size_t point_count =
      trial_params_soa_batch_override ? trial_params_soa_batch_override->size()
      : (prepared_runtime_batch ? prepared_runtime_batch->size()
                                : trial_params_batch.size());
  out_probabilities_matrix.assign(node_count * point_count, 0.0);
  if (node_count == 0u) {
    return true;
  }
  if (point_count == 0u || !(upper > 0.0)) {
    return point_count == 0u || upper <= 0.0;
  }

  std::vector<const uuber::TrialParamsSoA *> trial_params_soa_storage;
  std::vector<const TrialParamSet *> local_unique_trial_params;
  std::vector<const uuber::TrialParamsSoA *> local_unique_params_soa;
  std::vector<std::size_t> local_point_to_eval_index;
  std::vector<const TrialParamSet *> &unique_trial_params =
      workspace ? workspace->unique_trial_params : local_unique_trial_params;
  std::vector<const uuber::TrialParamsSoA *> &unique_params_soa =
      workspace ? workspace->unique_params_soa : local_unique_params_soa;
  std::vector<std::size_t> &point_to_eval_index =
      workspace ? workspace->point_to_eval_index : local_point_to_eval_index;
  point_to_eval_index.assign(point_count, 0u);

  const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch_ptr =
      trial_params_soa_batch_override;
  const std::vector<const TrialParamSet *> *trial_params_eval_ptr =
      &trial_params_batch;
  if (trial_params_soa_batch_ptr == nullptr && prepared_runtime_batch != nullptr) {
    if (prepared_runtime_batch->size() != point_count) {
      return false;
    }
    unique_trial_params.clear();
    unique_params_soa.clear();
    unique_trial_params.reserve(point_count);
    unique_params_soa.reserve(point_count);
    std::unordered_map<uuber::NAMapCacheKey, std::size_t, uuber::NAMapCacheKeyHash>
        runtime_index_by_key;
    runtime_index_by_key.reserve(point_count);
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      PreparedTrialParamsRuntime *runtime =
          const_cast<PreparedTrialParamsRuntime *>(
              (*prepared_runtime_batch)[point_idx]);
      if (runtime == nullptr || !ensure_prepared_trial_params_soa(ctx, *runtime)) {
        return false;
      }
      const auto insert_result = runtime_index_by_key.emplace(
          runtime->value_key, unique_params_soa.size());
      if (insert_result.second) {
        unique_trial_params.push_back(runtime->params);
        unique_params_soa.push_back(runtime->soa);
      }
      point_to_eval_index[point_idx] = insert_result.first->second;
    }
    trial_params_eval_ptr = &unique_trial_params;
    trial_params_soa_batch_ptr = &unique_params_soa;
  } else if (trial_params_soa_batch_ptr == nullptr) {
    trial_params_soa_storage.reserve(point_count);
    for (const TrialParamSet *params_ptr : trial_params_batch) {
      const uuber::TrialParamsSoA *trial_params_soa =
          resolve_trial_params_soa(ctx, params_ptr);
      if (!trial_params_soa || !trial_params_soa->valid) {
        return false;
      }
      trial_params_soa_storage.push_back(trial_params_soa);
    }
    trial_params_soa_batch_ptr = &trial_params_soa_storage;
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      point_to_eval_index[point_idx] = point_idx;
    }
  } else if (trial_params_soa_batch_ptr->size() != point_count) {
    return false;
  } else {
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      point_to_eval_index[point_idx] = point_idx;
    }
  }

  for (const uuber::TrialParamsSoA *trial_params_soa :
       *trial_params_soa_batch_ptr) {
    if (!trial_params_soa || !trial_params_soa->valid) {
      return false;
    }
  }

  TrialParamSet base_params_holder;
  const TrialParamSet *rep_params =
      trial_params_eval_ptr->empty() ? nullptr : trial_params_eval_ptr->front();
  if (!rep_params) {
    base_params_holder = build_base_paramset(ctx);
    rep_params = &base_params_holder;
  }

  const std::size_t eval_point_count = trial_params_soa_batch_ptr->size();
  SharedTriggerPlan local_trigger_plan;
  const SharedTriggerPlan *plan_ptr = nullptr;
  bool use_shared_trigger_eval_local = false;
  if (eval_point_count == 1u) {
    local_trigger_plan = build_shared_trigger_plan(ctx, rep_params);
    plan_ptr = &local_trigger_plan;
    use_shared_trigger_eval_local = shared_trigger_count(*plan_ptr) > 0u;
  }
  constexpr int kIntegralSegments = 16;
  constexpr int kInfiniteMaxSegments = 32;
  const bool finite_upper = std::isfinite(upper);
  std::vector<double> upper_bounds(eval_point_count, upper);
  std::vector<TimeIntegralPointPlan> local_point_plans;
  std::vector<double> local_query_times;
  std::vector<double> local_query_weights;
  std::vector<const uuber::TrialParamsSoA *> local_query_params_soa;
  std::vector<std::uint8_t> local_active_points;
  std::vector<std::uint8_t> local_tail_saw_positive;
  std::vector<int> local_tail_small_segments;
  std::vector<double> local_density_out;
  std::vector<double> local_point_probabilities_matrix;
  std::vector<double> local_segment_probabilities_matrix;
  std::vector<TimeIntegralPointPlan> &point_plans =
      workspace ? workspace->point_plans : local_point_plans;
  std::vector<double> &query_times =
      workspace ? workspace->query_times : local_query_times;
  std::vector<double> &query_weights =
      workspace ? workspace->query_weights : local_query_weights;
  std::vector<const uuber::TrialParamsSoA *> &query_params_soa =
      workspace ? workspace->query_params_soa : local_query_params_soa;
  std::vector<std::uint8_t> &active_points =
      workspace ? workspace->active_points : local_active_points;
  std::vector<std::uint8_t> &tail_saw_positive =
      workspace ? workspace->tail_saw_positive : local_tail_saw_positive;
  std::vector<int> &tail_small_segments =
      workspace ? workspace->tail_small_segments : local_tail_small_segments;
  std::vector<double> &density_out =
      workspace ? workspace->batch_density : local_density_out;
  std::vector<double> &segment_probabilities_matrix =
      workspace ? workspace->segment_probabilities
                : local_segment_probabilities_matrix;
  std::vector<double> &point_probabilities_matrix =
      workspace ? workspace->point_probabilities : local_point_probabilities_matrix;
  point_probabilities_matrix.assign(node_count * eval_point_count, 0.0);
  segment_probabilities_matrix.assign(node_count * eval_point_count, 0.0);
  active_points.assign(eval_point_count, 1u);
  tail_saw_positive.assign(eval_point_count, 0u);
  tail_small_segments.assign(eval_point_count, 0);

  const auto eval_segment = [&](double lower, double upper_window) -> bool {
    if (!(std::isfinite(lower) && std::isfinite(upper_window)) ||
        !(upper_window > lower)) {
      segment_probabilities_matrix.assign(node_count * eval_point_count, 0.0);
      return true;
    }
    if (!build_time_integral_query_plan(
            upper_bounds, lower, upper_window, kIntegralSegments, point_plans,
            query_times, query_weights, trial_params_soa_batch_ptr,
            &query_params_soa, &active_points)) {
      return false;
    }
    if (query_times.empty()) {
      segment_probabilities_matrix.assign(node_count * eval_point_count, 0.0);
      return true;
    }
    const std::vector<const uuber::TrialParamsSoA *> *query_params_soa_arg =
        use_shared_trigger_eval_local ? nullptr : &query_params_soa;

    const bool grouped =
        !use_shared_trigger_eval_local &&
        evaluate_node_density_multi_generic_batch_idx(
            ctx, node_ids, query_times, component_idx, competitor_ids,
            rep_params, trial_type_key, include_na_donors, outcome_idx_context,
            density_out, query_params_soa_arg);
    if (!grouped) {
      density_out.assign(node_count * query_times.size(), 0.0);
      for (std::size_t node_pos = 0; node_pos < node_count; ++node_pos) {
        std::vector<double> node_density;
        if (!node_density_entry_batch_idx(
                ctx, node_ids[node_pos], query_times, component_idx, nullptr,
                false, nullptr, false, competitor_ids, rep_params,
                trial_type_key, include_na_donors, outcome_idx_context, plan_ptr,
                use_shared_trigger_eval_local, node_density, nullptr,
                query_params_soa_arg, nullptr, workspace) ||
            node_density.size() != query_times.size()) {
          return false;
        }
        std::copy(node_density.begin(), node_density.end(),
                  density_out.begin() +
                      static_cast<std::ptrdiff_t>(node_pos * query_times.size()));
      }
    } else if (density_out.size() != node_count * query_times.size()) {
      return false;
    }

    for (std::size_t node_pos = 0; node_pos < node_count; ++node_pos) {
      std::vector<double> node_segment_probabilities;
      if (!integrate_query_plan_density_per_point_slice(
              point_plans, query_weights,
              density_out.data() +
                  static_cast<std::ptrdiff_t>(node_pos * query_times.size()),
              query_times.size(), node_segment_probabilities) ||
          node_segment_probabilities.size() != eval_point_count) {
        return false;
      }
      std::copy(node_segment_probabilities.begin(),
                node_segment_probabilities.end(),
                segment_probabilities_matrix.begin() +
                    static_cast<std::ptrdiff_t>(node_pos * eval_point_count));
    }
    return true;
  };

  if (finite_upper) {
    if (!eval_segment(0.0, upper)) {
      return false;
    }
    for (std::size_t node_pos = 0; node_pos < node_count; ++node_pos) {
      double *prob_row = point_probabilities_matrix.data() +
                         node_pos * static_cast<std::ptrdiff_t>(eval_point_count);
      const double *seg_row = segment_probabilities_matrix.data() +
                              node_pos * static_cast<std::ptrdiff_t>(eval_point_count);
      for (std::size_t i = 0; i < eval_point_count; ++i) {
        prob_row[i] = normalize_outcome_probability(seg_row[i]);
      }
    }
  } else {
    double lo = 0.0;
    double hi = 0.5;
    for (int seg = 0; seg < kInfiniteMaxSegments; ++seg) {
      if (!eval_segment(lo, hi)) {
        return false;
      }
      bool any_active = false;
      for (std::size_t point_idx = 0; point_idx < eval_point_count; ++point_idx) {
        if (point_plans[point_idx].count == 0u) {
          active_points[point_idx] = 0u;
          continue;
        }
        double max_piece = 0.0;
        double max_total = 0.0;
        for (std::size_t node_pos = 0; node_pos < node_count; ++node_pos) {
          double *prob_row = point_probabilities_matrix.data() +
                             node_pos * static_cast<std::ptrdiff_t>(eval_point_count);
          const double *seg_row =
              segment_probabilities_matrix.data() +
              node_pos * static_cast<std::ptrdiff_t>(eval_point_count);
          const double piece = seg_row[point_idx];
          prob_row[point_idx] += piece;
          max_piece = std::max(max_piece, piece);
          max_total = std::max(max_total, std::fabs(prob_row[point_idx]));
        }
        if (max_piece > 0.0) {
          tail_saw_positive[point_idx] = 1u;
        }
        const double tail_tol =
            std::max(std::max(0.0, abs_tol),
                     std::max(0.0, rel_tol) * std::max(1.0, max_total));
        if (max_piece <= tail_tol) {
          ++tail_small_segments[point_idx];
        } else {
          tail_small_segments[point_idx] = 0;
        }
        const bool done =
            tail_saw_positive[point_idx] != 0u && hi >= 4.0 &&
            tail_small_segments[point_idx] >= 3;
        active_points[point_idx] = done ? 0u : 1u;
        any_active = any_active || active_points[point_idx] != 0u;
      }
      if (!any_active) {
        break;
      }
      lo = hi;
      hi *= 2.0;
    }

    for (std::size_t node_pos = 0; node_pos < node_count; ++node_pos) {
      double *prob_row = point_probabilities_matrix.data() +
                         node_pos * static_cast<std::ptrdiff_t>(eval_point_count);
      for (std::size_t i = 0; i < eval_point_count; ++i) {
        prob_row[i] = normalize_outcome_probability(prob_row[i]);
      }
    }
  }

  for (std::size_t node_pos = 0; node_pos < node_count; ++node_pos) {
    double *out_row = out_probabilities_matrix.data() +
                      node_pos * static_cast<std::ptrdiff_t>(point_count);
    const double *prob_row =
        point_probabilities_matrix.data() +
        node_pos * static_cast<std::ptrdiff_t>(eval_point_count);
    for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
      const std::size_t eval_idx = point_to_eval_index[point_idx];
      if (eval_idx >= eval_point_count) {
        return false;
      }
      out_row[point_idx] = prob_row[eval_idx];
    }
  }
  return true;
}

double native_outcome_probability_bits_impl_idx(
    const uuber::NativeContext &ctx, int node_id, double upper, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids_raw, double rel_tol, double abs_tol,
    const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    bool include_na_donors,
    int outcome_idx_context) {
  if (upper <= 0.0) {
    return 0.0;
  }
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.components.ids.size())) {
    component_idx = -1;
  }
  std::vector<int> comp_vec_filtered;
  const std::vector<int> &comp_vec = filter_competitor_ids(
      ctx, competitor_ids_raw, component_idx, comp_vec_filtered);
  std::vector<const TrialParamSet *> trial_params_batch(1u, trial_params);
  std::vector<double> out_probabilities;
  if (!integrate_outcome_probability_from_density_batch_idx(
      ctx, node_id, upper, component_idx, forced_complete_bits,
      forced_complete_bits_valid, forced_survive_bits,
      forced_survive_bits_valid, comp_vec, rel_tol, abs_tol,
      trial_params_batch, trial_type_key, include_na_donors,
      outcome_idx_context, out_probabilities, nullptr, nullptr, nullptr) ||
      out_probabilities.size() != 1u) {
    return 0.0;
  }
  return out_probabilities[0];
}

double native_outcome_probability_bits_impl_idx(
    const uuber::NativeContext &ctx, int node_id, double upper, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const Rcpp::IntegerVector &competitor_ids, double rel_tol, double abs_tol,
    const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    bool include_na_donors,
    int outcome_idx_context) {
  std::vector<int> competitor_ids_raw = integer_vector_to_std(competitor_ids, false);
  return native_outcome_probability_bits_impl_idx(
      ctx, node_id, upper, component_idx,
      forced_complete_bits, forced_complete_bits_valid,
      forced_survive_bits, forced_survive_bits_valid,
      competitor_ids_raw, rel_tol, abs_tol, trial_params,
      trial_type_key, include_na_donors, outcome_idx_context);
}

double native_outcome_probability_impl_idx(
    const uuber::NativeContext &ctx, int node_id, double upper, int component_idx,
    SEXP forced_complete, SEXP forced_survive,
    const Rcpp::IntegerVector &competitor_ids, double rel_tol, double abs_tol,
    const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    bool include_na_donors,
    int outcome_idx_context) {
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid = false;
  bool forced_survive_bits_valid = false;
  build_forced_bitset_strict(ctx, fc_vec, forced_complete_bits,
                             forced_complete_bits_valid);
  build_forced_bitset_strict(ctx, fs_vec, forced_survive_bits,
                             forced_survive_bits_valid);
  return native_outcome_probability_bits_impl_idx(
      ctx, node_id, upper, component_idx,
      forced_complete_bits_valid ? &forced_complete_bits : nullptr,
      forced_complete_bits_valid,
      forced_survive_bits_valid ? &forced_survive_bits : nullptr,
      forced_survive_bits_valid, competitor_ids, rel_tol, abs_tol,
      trial_params, trial_type_key, include_na_donors, outcome_idx_context);
}

double native_outcome_probability_impl_idx(
    SEXP ctxSEXP, int node_id, double upper, int component_idx,
    SEXP forced_complete, SEXP forced_survive,
    const Rcpp::IntegerVector &competitor_ids, double rel_tol, double abs_tol,
    const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  return native_outcome_probability_impl_idx(
      *ctx, node_id, upper, component_idx, forced_complete, forced_survive,
      competitor_ids, rel_tol, abs_tol, trial_params, trial_type_key,
      include_na_donors, outcome_idx_context);
}

} // namespace

std::string evaluator_component_cache_key(const uuber::NativeContext &ctx,
                                          int component_idx,
                                          const std::string &trial_key) {
  return component_cache_key(ctx, component_idx, trial_key);
}

const TrialAccumulatorParams *
evaluator_get_trial_param_entry(const TrialParamSet *trial_params,
                                int acc_index) {
  return get_trial_param_entry(trial_params, acc_index);
}

void evaluator_resolve_event_numeric_params(
    const uuber::NativeAccumulator &acc, int acc_index,
    const TrialAccumulatorParams *override,
    const uuber::TrialParamsSoA *trial_params_soa, double &onset_out,
    double &q_out, uuber::AccDistParams &cfg_out) {
  resolve_event_numeric_params(acc, acc_index, override, trial_params_soa,
                               onset_out, q_out, cfg_out);
}

bool evaluator_component_active_idx(const uuber::NativeAccumulator &acc,
                                    int component_idx,
                                    const TrialAccumulatorParams *override) {
  return component_active_idx(acc, component_idx, override);
}

GuardEvalInput evaluator_make_guard_input(
    const uuber::NativeContext &ctx, int node_idx, int component_idx,
    const std::string *trial_type_key, const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *trial_params_soa,
    const TimeConstraintMap *time_constraints,
    const ForcedScopeFilter *forced_scope_filter,
    const uuber::BitsetState *forced_complete_bits,
    const uuber::BitsetState *forced_survive_bits,
    const std::unordered_map<int, int> *forced_label_id_to_bit_idx,
    const ForcedStateView *forced_state_view) {
  return make_guard_input(
      ctx, node_idx, component_idx, trial_type_key, trial_params,
      trial_params_soa, time_constraints,
      forced_scope_filter, forced_complete_bits, forced_survive_bits,
      forced_label_id_to_bit_idx, forced_state_view);
}

std::vector<int> evaluator_gather_blocker_sources(
    const uuber::NativeContext &ctx, int guard_node_idx) {
  return gather_blocker_sources(ctx, guard_node_idx);
}

bool evaluator_guard_cdf_batch_prepared(const GuardEvalInput &input,
                                        const std::vector<double> &times,
                                        std::vector<double> &cdf_out) {
  return guard_cdf_batch_prepared_internal(input, times, cdf_out);
}

uuber::TreeEventBatchEvalFn
evaluator_make_tree_event_eval_batch(NodeEvalState &state) {
  return make_tree_event_eval_batch(state);
}

uuber::TreeGuardBatchEvalFn
evaluator_make_tree_guard_eval_batch(NodeEvalState &state) {
  return make_tree_guard_eval_batch(state);
}

bool evaluator_eval_node_batch_with_state_dense(
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values) {
  return eval_node_batch_with_state_dense(node_idx, times, state, need,
                                          out_values);
}

bool evaluator_eval_nodes_batch_with_state_dense(
    const std::vector<int> &node_indices, const std::vector<double> &times,
    NodeEvalState &state, std::vector<uuber::TreeNodeBatchValues> &out_values) {
  out_values.clear();
  if (times.empty()) {
    return true;
  }
  uuber::TreeEventBatchEvalFn event_eval_batch_cb =
      make_tree_event_eval_batch(state);
  uuber::TreeGuardBatchEvalFn guard_eval_batch_cb =
      make_tree_guard_eval_batch(state);
  const uuber::TreeEvalNeed full_need{true, true, true};
  return eval_tree_vector_nodes_batch(tree_vector_program_required(state.ctx),
                                      node_indices, times, full_need,
                                      event_eval_batch_cb,
                                      guard_eval_batch_cb, out_values);
}

bool evaluator_node_density_with_competitors_batch_internal(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<double> &times, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid, const std::vector<int> &competitor_ids,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    std::vector<double> &density_out,
    const TimeConstraintMap *time_constraints,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch) {
  return node_density_with_competitors_batch_internal(
      ctx, node_id, times, component_idx, forced_complete_bits,
      forced_complete_bits_valid, forced_survive_bits,
      forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
      include_na_donors, outcome_idx_context, density_out, time_constraints,
      trial_params_soa_batch);
}

bool mix_outcome_mass_batch_idx(
    uuber::NativeContext &ctx,
    const std::vector<int> &nonresponse_outcome_indices,
    const std::vector<int> &component_indices,
    const std::vector<const std::vector<double> *> &component_weights_batch,
    const std::vector<const TrialParamSet *> &trial_params_batch, double rel_tol,
    double abs_tol, int max_depth, bool complement_probability,
    std::vector<double> &out_probabilities,
    const std::vector<std::size_t> *trial_to_point_index = nullptr,
    const std::vector<const uuber::TrialParamsSoA *>
        *trial_params_soa_batch_override = nullptr,
    const std::vector<const PreparedTrialParamsRuntime *>
        *prepared_runtime_batch = nullptr) {
  out_probabilities.clear();
  const std::size_t point_count =
      trial_params_soa_batch_override ? trial_params_soa_batch_override->size()
      : (prepared_runtime_batch ? prepared_runtime_batch->size()
                                : trial_params_batch.size());
  const std::size_t trial_count = component_weights_batch.size();
  if (trial_count == 0u) {
    return true;
  }
  if (point_count == 0u) {
    return false;
  }
  if (trial_to_point_index != nullptr &&
      trial_to_point_index->size() != trial_count) {
    return false;
  }
  if (trial_params_batch.size() < point_count) {
    return false;
  }

  std::vector<const uuber::TrialParamsSoA *> trial_params_soa_batch_storage;
  const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch_ptr =
      trial_params_soa_batch_override;
  if (trial_params_soa_batch_ptr == nullptr && prepared_runtime_batch != nullptr) {
    if (prepared_runtime_batch->size() != point_count) {
      return false;
    }
    for (const PreparedTrialParamsRuntime *runtime_const : *prepared_runtime_batch) {
      PreparedTrialParamsRuntime *runtime =
          const_cast<PreparedTrialParamsRuntime *>(runtime_const);
      if (runtime == nullptr || !ensure_prepared_trial_params_soa(ctx, *runtime)) {
        return false;
      }
    }
  } else if (trial_params_soa_batch_ptr == nullptr) {
    trial_params_soa_batch_storage.reserve(point_count);
    for (const TrialParamSet *params_ptr : trial_params_batch) {
      const uuber::TrialParamsSoA *trial_params_soa =
          resolve_trial_params_soa(ctx, params_ptr);
      if (!trial_params_soa || !trial_params_soa->valid) {
        return false;
      }
      trial_params_soa_batch_storage.push_back(trial_params_soa);
    }
    trial_params_soa_batch_ptr = &trial_params_soa_batch_storage;
  } else {
    for (const uuber::TrialParamsSoA *trial_params_soa :
         *trial_params_soa_batch_ptr) {
      if (!trial_params_soa || !trial_params_soa->valid) {
        return false;
      }
    }
  }

  out_probabilities.assign(trial_count, 0.0);
  DensityBatchWorkspace integration_workspace;
  const bool use_implicit_component = component_indices.empty();
  const std::size_t component_count =
      use_implicit_component ? 1u : component_indices.size();
  struct OutcomeMassComponentUse {
    std::size_t trial_pos{0u};
    std::size_t point_idx{0u};
    double mix_weight{0.0};
  };
  std::vector<std::vector<OutcomeMassComponentUse>> component_uses(
      component_count);
  for (std::size_t component_pos = 0; component_pos < component_count;
       ++component_pos) {
    std::vector<OutcomeMassComponentUse> &uses = component_uses[component_pos];
    uses.reserve(trial_count);
    for (std::size_t trial_pos = 0; trial_pos < trial_count; ++trial_pos) {
      const std::vector<double> *trial_weights = component_weights_batch[trial_pos];
      const double mix_weight = use_implicit_component
                                    ? 1.0
                                    : ((trial_weights &&
                                        component_pos < trial_weights->size())
                                           ? (*trial_weights)[component_pos]
                                           : 0.0);
      if (!std::isfinite(mix_weight) || mix_weight <= 0.0) {
        continue;
      }
      const std::size_t point_idx =
          (trial_to_point_index != nullptr) ? (*trial_to_point_index)[trial_pos]
                                            : trial_pos;
      if (point_idx >= point_count) {
        return false;
      }
      uses.push_back(
          OutcomeMassComponentUse{trial_pos, point_idx, mix_weight});
    }
  }

  std::vector<double> batch_probabilities;
  batch_probabilities.reserve(point_count);
  for (std::size_t outcome_pos = 0;
       outcome_pos < nonresponse_outcome_indices.size(); ++outcome_pos) {
    const int outcome_idx = nonresponse_outcome_indices[outcome_pos];
    if (outcome_idx < 0 ||
        outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
      continue;
    }

    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
    if (info.node_id < 0) {
      continue;
    }

    for (std::size_t component_pos = 0; component_pos < component_count;
         ++component_pos) {
      const std::vector<OutcomeMassComponentUse> &uses =
          component_uses[component_pos];
      if (uses.empty()) {
        continue;
      }
      const int component_idx =
          use_implicit_component ? -1 : component_indices[component_pos];
      if (!use_implicit_component &&
          !outcome_allows_component_idx(ctx, info, outcome_idx, component_idx)) {
        continue;
      }

      const double keep_weight = use_implicit_component
                                     ? 1.0
                                     : component_keep_weight(ctx, component_idx,
                                                             outcome_idx);
      if (!std::isfinite(keep_weight) || keep_weight <= 0.0) {
        continue;
      }

      std::vector<int> filtered_competitors;
      const std::vector<int> &competitors = filter_competitor_ids(
          ctx, info.competitor_ids, component_idx, filtered_competitors);

      batch_probabilities.clear();
      if (!integrate_outcome_probability_from_density_batch_idx(
              ctx, info.node_id, std::numeric_limits<double>::infinity(),
              component_idx, nullptr, false, nullptr, false, competitors,
              rel_tol, abs_tol, trial_params_batch, std::string(),
              false, outcome_idx, batch_probabilities,
              trial_params_soa_batch_ptr, prepared_runtime_batch,
              &integration_workspace) ||
          batch_probabilities.size() != point_count) {
        return false;
      }

      for (const OutcomeMassComponentUse &use : uses) {
        const double probability = batch_probabilities[use.point_idx];
        if (std::isfinite(probability) && probability > 0.0) {
          out_probabilities[use.trial_pos] +=
              use.mix_weight * keep_weight * probability;
        }
      }
    }
  }

  for (double &probability : out_probabilities) {
    if (!std::isfinite(probability)) {
      probability = 0.0;
      continue;
    }
    if (complement_probability) {
      probability = clamp_probability(1.0 - probability);
    } else if (probability <= 0.0) {
      probability = 0.0;
    }
  }

  return true;
}

enum class TrialKernelIntent : std::uint8_t {
  SequenceDensity = 0,
  OutcomeMass = 1
};

enum class TrialSequenceExecutionKind : std::uint8_t {
  LowerLayerDirect = 0,
  SequenceState = 1
};

enum class TrialProbabilityTransform : std::uint8_t {
  Identity = 0,
  Complement = 1
};

struct TrialContributionSpec {
  TrialKernelIntent intent{TrialKernelIntent::SequenceDensity};
  TrialSequenceExecutionKind sequence_execution{
      TrialSequenceExecutionKind::LowerLayerDirect};
  double scaled_weight{1.0};
  int component_idx{-1};
  const std::string *trial_type_key_ptr{nullptr};
  std::string trial_type_key_storage;
  std::vector<int> sequence_outcome_indices;
  std::vector<int> sequence_node_indices;
  const double *sequence_time_data{nullptr};
  std::vector<const std::vector<int> *> step_competitor_ids_ptrs;
  std::vector<std::vector<int>> step_competitor_ids_storage;
  std::vector<std::vector<int>> step_persistent_sources;
};

inline const std::string &empty_trial_type_key_ref() {
  static const std::string kEmptyTrialTypeKey;
  return kEmptyTrialTypeKey;
}

struct TrialEvalInput {
  bool valid{false};
  TrialProbabilityTransform probability_transform{
      TrialProbabilityTransform::Identity};
  bool enforce_component_outcome_gate{false};
  int single_label_storage{NA_INTEGER};
  double single_time_storage{std::numeric_limits<double>::quiet_NaN()};
  std::vector<int> sequence_label_storage;
  std::vector<double> sequence_time_storage;
  const int *sequence_label_data{nullptr};
  const double *sequence_time_data{nullptr};
  std::size_t sequence_length{0u};
  const std::vector<int> *nonresponse_outcome_indices{nullptr};
};

struct LogLikWorkload {
  R_xlen_t n_rows{0};
  int n_trials{0};
  bool ranked_mode{false};
  bool has_onset{false};
  bool enable_ranked_sequence_batch{false};
  std::vector<int> row_to_trial_index;
  std::vector<int> row_acc_indices;
  std::vector<double> onset_by_row;
  std::vector<int> forced_component_idx_by_trial;
  std::vector<TrialEvalInput> trial_eval_inputs;
  std::vector<int> nonresponse_outcome_indices;
  std::vector<int> default_component_indices;
  std::vector<ComponentCacheEntry> default_cache_entries;
  std::vector<std::vector<int>> leader_components_by_acc_idx;
  std::vector<std::uint8_t> ok_by_trial;
  std::vector<int> expand_trial_indices;
};

struct LogLikPreparedParams {
  std::vector<TrialParamSet> param_sets;
  std::vector<PreparedTrialParamsRuntime> prepared_trial_runtimes;
  std::vector<const PreparedTrialParamsRuntime *> prepared_trial_runtime_ptrs;
  std::vector<SharedTriggerPlan> trigger_plans;
  std::vector<int> trigger_plan_index_by_trial;
  std::vector<int> param_value_group_by_trial;
  std::vector<int> param_value_group_representative_trial;
  std::vector<std::vector<double>> weights_by_trial;
  std::vector<int> component_weight_group_by_trial;
  std::vector<int> component_weight_group_representative_trial;
};

struct OutcomeMassBatchKey {
  TrialProbabilityTransform probability_transform{
      TrialProbabilityTransform::Identity};
  bool valid{false};
  bool enforce_component_outcome_gate{false};
  const std::vector<int> *nonresponse_outcome_indices{nullptr};
  const std::vector<int> *component_indices{nullptr};

  bool operator==(const OutcomeMassBatchKey &other) const {
    return probability_transform == other.probability_transform &&
           valid == other.valid &&
           enforce_component_outcome_gate == other.enforce_component_outcome_gate &&
           nonresponse_outcome_indices == other.nonresponse_outcome_indices &&
           component_indices == other.component_indices;
  }
};

struct OutcomeMassBatchKeyHash {
  std::size_t operator()(const OutcomeMassBatchKey &key) const {
    std::uint64_t hash = kFNV64Offset;
    const int transform_code = static_cast<int>(key.probability_transform);
    hash_append_bytes(hash, &transform_code, sizeof(transform_code));
    hash_append_bool(hash, key.valid);
    hash_append_bool(hash, key.enforce_component_outcome_gate);
    hash_append_u64(
        hash, static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(
                  key.nonresponse_outcome_indices)));
    hash_append_u64(hash, static_cast<std::uint64_t>(
                              reinterpret_cast<std::uintptr_t>(
                                  key.component_indices)));
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

struct OutcomeMassSharedTriggerBatchKey {
  TrialProbabilityTransform probability_transform{
      TrialProbabilityTransform::Identity};
  bool valid{false};
  bool enforce_component_outcome_gate{false};
  int param_group_idx{-1};
  int trigger_plan_idx{-1};
  const std::vector<int> *nonresponse_outcome_indices{nullptr};
  const std::vector<int> *component_indices{nullptr};

  bool operator==(const OutcomeMassSharedTriggerBatchKey &other) const {
    return probability_transform == other.probability_transform &&
           valid == other.valid &&
           enforce_component_outcome_gate == other.enforce_component_outcome_gate &&
           param_group_idx == other.param_group_idx &&
           trigger_plan_idx == other.trigger_plan_idx &&
           nonresponse_outcome_indices == other.nonresponse_outcome_indices &&
           component_indices == other.component_indices;
  }
};

struct OutcomeMassSharedTriggerBatchKeyHash {
  std::size_t operator()(const OutcomeMassSharedTriggerBatchKey &key) const {
    std::uint64_t hash = kFNV64Offset;
    const int transform_code = static_cast<int>(key.probability_transform);
    hash_append_bytes(hash, &transform_code, sizeof(transform_code));
    hash_append_bool(hash, key.valid);
    hash_append_bool(hash, key.enforce_component_outcome_gate);
    hash_append_bytes(hash, &key.param_group_idx, sizeof(key.param_group_idx));
    hash_append_bytes(hash, &key.trigger_plan_idx, sizeof(key.trigger_plan_idx));
    hash_append_u64(
        hash, static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(
                  key.nonresponse_outcome_indices)));
    hash_append_u64(hash, static_cast<std::uint64_t>(
                              reinterpret_cast<std::uintptr_t>(
                                  key.component_indices)));
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

inline bool build_outcome_mass_batch_key(const TrialEvalInput &eval_input,
                                         const std::vector<int> &component_indices,
                                         OutcomeMassBatchKey &out) {
  out = OutcomeMassBatchKey{};
  if (!eval_input.valid || eval_input.sequence_length != 0u ||
      eval_input.nonresponse_outcome_indices == nullptr ||
      eval_input.nonresponse_outcome_indices->empty()) {
    return false;
  }
  out.probability_transform = eval_input.probability_transform;
  out.valid = eval_input.valid;
  out.enforce_component_outcome_gate = eval_input.enforce_component_outcome_gate;
  out.nonresponse_outcome_indices = eval_input.nonresponse_outcome_indices;
  out.component_indices = &component_indices;
  return true;
}

inline bool build_outcome_mass_shared_trigger_batch_key(
    const TrialEvalInput &eval_input, const std::vector<int> &component_indices,
    int param_group_idx, int trigger_plan_idx,
    OutcomeMassSharedTriggerBatchKey &out) {
  if (param_group_idx < 0 || trigger_plan_idx < 0) {
    return false;
  }
  out = OutcomeMassSharedTriggerBatchKey{};
  if (!eval_input.valid || eval_input.sequence_length != 0u ||
      eval_input.nonresponse_outcome_indices == nullptr ||
      eval_input.nonresponse_outcome_indices->empty()) {
    return false;
  }
  out.probability_transform = eval_input.probability_transform;
  out.valid = eval_input.valid;
  out.enforce_component_outcome_gate = eval_input.enforce_component_outcome_gate;
  out.param_group_idx = param_group_idx;
  out.trigger_plan_idx = trigger_plan_idx;
  out.nonresponse_outcome_indices = eval_input.nonresponse_outcome_indices;
  out.component_indices = &component_indices;
  return true;
}

struct OutcomeMassTrialSignature {
  std::size_t point_idx{0u};
  int component_weight_group_idx{-1};

  bool operator==(const OutcomeMassTrialSignature &other) const {
    return point_idx == other.point_idx &&
           component_weight_group_idx == other.component_weight_group_idx;
  }
};

struct OutcomeMassTrialSignatureHash {
  std::size_t operator()(const OutcomeMassTrialSignature &key) const {
    std::uint64_t hash = kFNV64Offset;
    hash_append_u64(hash, static_cast<std::uint64_t>(key.point_idx));
    hash_append_bytes(hash, &key.component_weight_group_idx,
                      sizeof(key.component_weight_group_idx));
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

bool build_trial_contributions_unified(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    std::vector<TrialContributionSpec> &contributions);

struct SingleStepDirectBatchSpec {
  int node_id{-1};
  int component_idx{-1};
  int outcome_idx_context{-1};
  int param_group_idx{-1};
  int trigger_plan_idx{-1};
  double time{std::numeric_limits<double>::quiet_NaN()};
  double scaled_weight{0.0};
  std::string trial_type_key;
  std::vector<int> competitor_ids;
};

struct SingleStepDirectBatchKey {
  int node_id{-1};
  int component_idx{-1};
  int outcome_idx_context{-1};
  int param_group_idx{-1};
  int trigger_plan_idx{-1};
  std::string trial_type_key;
  std::vector<int> competitor_ids;

  bool operator==(const SingleStepDirectBatchKey &other) const {
    return node_id == other.node_id && component_idx == other.component_idx &&
           outcome_idx_context == other.outcome_idx_context &&
           param_group_idx == other.param_group_idx &&
           trigger_plan_idx == other.trigger_plan_idx &&
           trial_type_key == other.trial_type_key &&
           competitor_ids == other.competitor_ids;
  }
};

struct SingleStepDirectBatchKeyHash {
  std::size_t operator()(const SingleStepDirectBatchKey &key) const {
    std::uint64_t hash = kFNV64Offset;
    hash_append_bytes(hash, &key.node_id, sizeof(key.node_id));
    hash_append_bytes(hash, &key.component_idx, sizeof(key.component_idx));
    hash_append_bytes(hash, &key.outcome_idx_context,
                      sizeof(key.outcome_idx_context));
    hash_append_bytes(hash, &key.param_group_idx,
                      sizeof(key.param_group_idx));
    hash_append_bytes(hash, &key.trigger_plan_idx,
                      sizeof(key.trigger_plan_idx));
    hash_append_u64(hash,
                    static_cast<std::uint64_t>(key.trial_type_key.size()));
    if (!key.trial_type_key.empty()) {
      hash_append_bytes(hash, key.trial_type_key.data(),
                        key.trial_type_key.size());
    }
    hash_append_u64(hash,
                    static_cast<std::uint64_t>(key.competitor_ids.size()));
    for (int competitor_id : key.competitor_ids) {
      hash_append_bytes(hash, &competitor_id, sizeof(competitor_id));
    }
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

struct DirectContributionTemplate {
  int node_id{-1};
  int component_idx{-1};
  int outcome_idx_context{-1};
  double scaled_weight{0.0};
  std::string trial_type_key;
  std::vector<int> competitor_ids;

  bool operator==(const DirectContributionTemplate &other) const {
    return node_id == other.node_id &&
           component_idx == other.component_idx &&
           outcome_idx_context == other.outcome_idx_context &&
           scaled_weight == other.scaled_weight &&
           trial_type_key == other.trial_type_key &&
           competitor_ids == other.competitor_ids;
  }
};

struct DirectTrialBatchSpec {
  std::vector<DirectContributionTemplate> contributions;
  int param_group_idx{-1};
  int trigger_plan_idx{-1};
  double time{std::numeric_limits<double>::quiet_NaN()};
};

struct DirectTrialBatchKey {
  std::vector<DirectContributionTemplate> contributions;
  int param_group_idx{-1};
  int trigger_plan_idx{-1};

  bool operator==(const DirectTrialBatchKey &other) const {
    return contributions == other.contributions &&
           param_group_idx == other.param_group_idx &&
           trigger_plan_idx == other.trigger_plan_idx;
  }
};

struct DirectTrialBatchKeyHash {
  std::size_t operator()(const DirectTrialBatchKey &key) const {
    std::uint64_t hash = kFNV64Offset;
    hash_append_u64(hash,
                    static_cast<std::uint64_t>(key.contributions.size()));
    for (const DirectContributionTemplate &contribution : key.contributions) {
      hash_append_bytes(hash, &contribution.node_id,
                        sizeof(contribution.node_id));
      hash_append_bytes(hash, &contribution.component_idx,
                        sizeof(contribution.component_idx));
      hash_append_bytes(hash, &contribution.outcome_idx_context,
                        sizeof(contribution.outcome_idx_context));
      hash_append_u64(hash, canonical_double_bits(contribution.scaled_weight));
      hash_append_u64(
          hash,
          static_cast<std::uint64_t>(contribution.trial_type_key.size()));
      if (!contribution.trial_type_key.empty()) {
        hash_append_bytes(hash, contribution.trial_type_key.data(),
                          contribution.trial_type_key.size());
      }
      hash_append_u64(
          hash,
          static_cast<std::uint64_t>(contribution.competitor_ids.size()));
      for (int competitor_id : contribution.competitor_ids) {
        hash_append_bytes(hash, &competitor_id, sizeof(competitor_id));
      }
    }
    hash_append_bytes(hash, &key.param_group_idx,
                      sizeof(key.param_group_idx));
    hash_append_bytes(hash, &key.trigger_plan_idx,
                      sizeof(key.trigger_plan_idx));
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

struct RankedSequenceContributionTemplate {
  int component_idx{-1};
  double scaled_weight{0.0};
  std::string trial_type_key;
  std::vector<int> sequence_outcome_indices;
  std::vector<int> sequence_node_indices;
  std::vector<std::vector<int>> step_competitor_ids;
  std::vector<const std::vector<int> *> step_competitor_ids_ptrs;
  std::vector<std::vector<int>> step_persistent_sources;

  void bind_step_competitor_ids_ptrs() {
    step_competitor_ids_ptrs.resize(step_competitor_ids.size(), nullptr);
    for (std::size_t i = 0; i < step_competitor_ids.size(); ++i) {
      step_competitor_ids_ptrs[i] = &step_competitor_ids[i];
    }
  }
};

struct RankedSequenceBatchTemplate {
  std::vector<RankedSequenceContributionTemplate> contributions;

  void bind_step_competitor_ids_ptrs() {
    for (RankedSequenceContributionTemplate &contribution : contributions) {
      contribution.bind_step_competitor_ids_ptrs();
    }
  }
};

enum class RankedSequenceTemplateRejectReason : std::uint8_t {
  None = 0,
  Contribution = 1,
  Shape = 2
};

struct RankedSequenceTemplateCacheEntry {
  bool buildable{false};
  RankedSequenceTemplateRejectReason reject_reason{
      RankedSequenceTemplateRejectReason::None};
  RankedSequenceBatchTemplate batch_template;
};

struct RankedSequenceBatchSpec {
  const RankedSequenceBatchTemplate *template_ptr{nullptr};
  int param_group_idx{-1};
  int trigger_plan_idx{-1};
  const double *time_data{nullptr};
  std::size_t time_count{0u};
};

struct RankedSequenceTemplateCacheKey {
  bool enforce_component_outcome_gate{false};
  std::vector<int> sequence_label_ids;
  std::vector<int> component_indices;
  std::vector<std::uint64_t> component_weight_bits;

  bool operator==(const RankedSequenceTemplateCacheKey &other) const {
    return enforce_component_outcome_gate ==
               other.enforce_component_outcome_gate &&
           sequence_label_ids == other.sequence_label_ids &&
           component_indices == other.component_indices &&
           component_weight_bits == other.component_weight_bits;
  }
};

struct RankedSequenceTemplateCacheKeyHash {
  std::size_t operator()(const RankedSequenceTemplateCacheKey &key) const {
    std::uint64_t hash = kFNV64Offset;
    hash_append_bool(hash, key.enforce_component_outcome_gate);
    auto append_int_vec = [&](const std::vector<int> &values) {
      hash_append_u64(hash, static_cast<std::uint64_t>(values.size()));
      for (int value : values) {
        hash_append_bytes(hash, &value, sizeof(value));
      }
    };
    auto append_u64_vec = [&](const std::vector<std::uint64_t> &values) {
      hash_append_u64(hash, static_cast<std::uint64_t>(values.size()));
      for (std::uint64_t value : values) {
        hash_append_u64(hash, value);
      }
    };
    append_int_vec(key.sequence_label_ids);
    append_int_vec(key.component_indices);
    append_u64_vec(key.component_weight_bits);
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

struct RankedSequenceBatchKey {
  const RankedSequenceBatchTemplate *template_ptr{nullptr};
  int param_group_idx{-1};
  int trigger_plan_idx{-1};

  bool operator==(const RankedSequenceBatchKey &other) const {
    return template_ptr == other.template_ptr &&
           param_group_idx == other.param_group_idx &&
           trigger_plan_idx == other.trigger_plan_idx;
  }
};

struct RankedSequenceBatchKeyHash {
  std::size_t operator()(const RankedSequenceBatchKey &key) const {
    std::uint64_t hash = kFNV64Offset;
    const std::uintptr_t template_bits =
        reinterpret_cast<std::uintptr_t>(key.template_ptr);
    hash_append_bytes(hash, &template_bits, sizeof(template_bits));
    hash_append_bytes(hash, &key.param_group_idx, sizeof(key.param_group_idx));
    hash_append_bytes(hash, &key.trigger_plan_idx,
                      sizeof(key.trigger_plan_idx));
    return static_cast<std::size_t>(mix_hash64(hash));
  }
};

using RankedSequenceTemplateCacheMap =
    std::unordered_map<RankedSequenceTemplateCacheKey,
                       RankedSequenceTemplateCacheEntry,
                       RankedSequenceTemplateCacheKeyHash>;

inline bool build_ranked_sequence_template_cache_key(
    const TrialEvalInput &eval_input, const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    RankedSequenceTemplateCacheKey &out) {
  if (eval_input.sequence_length <= 1u ||
      eval_input.sequence_label_data == nullptr) {
    return false;
  }
  out = RankedSequenceTemplateCacheKey{};
  out.enforce_component_outcome_gate = eval_input.enforce_component_outcome_gate;
  out.sequence_label_ids.assign(eval_input.sequence_label_data,
                                eval_input.sequence_label_data +
                                    eval_input.sequence_length);
  out.component_indices = component_indices;
  out.component_weight_bits.reserve(component_indices.size());
  for (std::size_t i = 0; i < component_indices.size(); ++i) {
    const double weight =
        (i < component_weights.size()) ? component_weights[i] : 0.0;
    out.component_weight_bits.push_back(canonical_double_bits(weight));
  }
  return true;
}

inline std::unordered_map<std::uint64_t, RankedSequenceTemplateCacheMap> &
ranked_sequence_template_cache_registry() {
  thread_local std::unordered_map<std::uint64_t, RankedSequenceTemplateCacheMap>
      cache_by_runtime_id;
  return cache_by_runtime_id;
}

inline RankedSequenceTemplateCacheMap &
ranked_sequence_template_cache_for_ctx(const uuber::NativeContext &ctx) {
  return ranked_sequence_template_cache_registry()[ctx.runtime_cache_instance_id];
}

namespace uuber {

void clear_ranked_distribution_runtime_caches(
    std::uint64_t runtime_cache_instance_id) noexcept {
  if (runtime_cache_instance_id == 0) {
    return;
  }
  ranked_sequence_template_cache_registry().erase(runtime_cache_instance_id);
  guard_prepared_conditioning_template_registry().erase(runtime_cache_instance_id);
  guard_prepared_cache_registry().erase(runtime_cache_instance_id);
}

} // namespace uuber

inline bool build_single_step_direct_batch_spec(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    SingleStepDirectBatchSpec &out) {
  out = SingleStepDirectBatchSpec{};
  if (!eval_input.valid ||
      eval_input.probability_transform != TrialProbabilityTransform::Identity ||
      eval_input.sequence_length != 1u || eval_input.sequence_time_data == nullptr) {
    return false;
  }

  std::vector<TrialContributionSpec> contributions;
  if (!build_trial_contributions_unified(
          ctx, eval_input, component_indices, component_weights, cache_entries,
          contributions) ||
      contributions.size() != 1u) {
    return false;
  }

  const TrialContributionSpec &spec = contributions.front();
  if (spec.intent != TrialKernelIntent::SequenceDensity ||
      spec.sequence_execution != TrialSequenceExecutionKind::LowerLayerDirect ||
      spec.sequence_node_indices.size() != 1u ||
      spec.sequence_outcome_indices.size() != 1u ||
      !std::isfinite(spec.scaled_weight) || spec.scaled_weight <= 0.0 ||
      spec.sequence_time_data == nullptr ||
      !std::isfinite(spec.sequence_time_data[0]) ||
      spec.sequence_time_data[0] < 0.0) {
    return false;
  }

  const int node_idx = spec.sequence_node_indices[0];
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  const int stable_node_id =
      ctx.ir.nodes[static_cast<std::size_t>(node_idx)].node_id;
  out.node_id = (stable_node_id >= 0) ? stable_node_id : node_idx;
  out.component_idx = spec.component_idx;
  out.outcome_idx_context = spec.sequence_outcome_indices[0];
  out.time = spec.sequence_time_data[0];
  out.scaled_weight = spec.scaled_weight;
  out.trial_type_key =
      spec.trial_type_key_ptr ? *spec.trial_type_key_ptr : std::string();
  if (!spec.step_competitor_ids_ptrs.empty() &&
      spec.step_competitor_ids_ptrs[0] != nullptr) {
    out.competitor_ids = *spec.step_competitor_ids_ptrs[0];
  }
  return true;
}

inline bool build_ranked_sequence_batch_spec(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    RankedSequenceBatchSpec &out) {
  out = RankedSequenceBatchSpec{};
  if (!eval_input.valid ||
      eval_input.probability_transform != TrialProbabilityTransform::Identity ||
      eval_input.sequence_length <= 1u || eval_input.sequence_time_data == nullptr) {
    return false;
  }

  RankedSequenceTemplateCacheKey cache_key;
  if (!build_ranked_sequence_template_cache_key(
          eval_input, component_indices, component_weights, cache_key)) {
    return false;
  }
  auto &template_cache = ranked_sequence_template_cache_for_ctx(ctx);
  auto template_it = template_cache.find(cache_key);
  if (template_it != template_cache.end()) {
    const RankedSequenceTemplateCacheEntry &entry = template_it->second;
    if (!entry.buildable) {
      return false;
    }
    out.template_ptr = &entry.batch_template;
    out.time_data = eval_input.sequence_time_data;
    out.time_count = eval_input.sequence_length;
    return true;
  }

  std::vector<TrialContributionSpec> contributions;
  if (!build_trial_contributions_unified(
          ctx, eval_input, component_indices, component_weights, cache_entries,
          contributions) ||
      contributions.empty()) {
    RankedSequenceTemplateCacheEntry entry;
    entry.buildable = false;
    entry.reject_reason = RankedSequenceTemplateRejectReason::Contribution;
    template_cache.emplace(std::move(cache_key), std::move(entry));
    return false;
  }

  for (const TrialContributionSpec &spec : contributions) {
    if (spec.intent != TrialKernelIntent::SequenceDensity ||
        spec.sequence_execution != TrialSequenceExecutionKind::SequenceState ||
        !std::isfinite(spec.scaled_weight) || spec.scaled_weight <= 0.0 ||
        spec.sequence_outcome_indices.empty() ||
        spec.sequence_outcome_indices.size() != spec.sequence_node_indices.size() ||
        spec.sequence_outcome_indices.size() != spec.step_persistent_sources.size() ||
        spec.sequence_outcome_indices.size() != spec.step_competitor_ids_ptrs.size() ||
        spec.sequence_time_data == nullptr) {
      RankedSequenceTemplateCacheEntry entry;
      entry.buildable = false;
      entry.reject_reason = RankedSequenceTemplateRejectReason::Shape;
      template_cache.emplace(std::move(cache_key), std::move(entry));
      return false;
    }
  }

  RankedSequenceTemplateCacheEntry entry;
  entry.buildable = true;
  RankedSequenceBatchTemplate &batch_template = entry.batch_template;
  out.time_data = eval_input.sequence_time_data;
  out.time_count = eval_input.sequence_length;
  batch_template.contributions.reserve(contributions.size());
  for (const TrialContributionSpec &spec : contributions) {
    RankedSequenceContributionTemplate contribution;
    contribution.component_idx = spec.component_idx;
    contribution.scaled_weight = spec.scaled_weight;
    contribution.trial_type_key =
        spec.trial_type_key_ptr ? *spec.trial_type_key_ptr : std::string();
    contribution.sequence_outcome_indices = spec.sequence_outcome_indices;
    contribution.sequence_node_indices = spec.sequence_node_indices;
    contribution.step_persistent_sources = spec.step_persistent_sources;
    contribution.step_competitor_ids.resize(spec.step_competitor_ids_ptrs.size());
    for (std::size_t i = 0; i < spec.step_competitor_ids_ptrs.size(); ++i) {
      if (spec.step_competitor_ids_ptrs[i] != nullptr) {
        contribution.step_competitor_ids[i] = *spec.step_competitor_ids_ptrs[i];
      }
    }
    batch_template.contributions.push_back(std::move(contribution));
  }
  auto inserted =
      template_cache.emplace(std::move(cache_key), std::move(entry));
  RankedSequenceTemplateCacheEntry &stored_entry = inserted.first->second;
  stored_entry.batch_template.bind_step_competitor_ids_ptrs();
  out.template_ptr = &stored_entry.batch_template;
  return true;
}

inline void reset_trial_eval_input(
    TrialEvalInput &eval_input,
    const std::vector<int> *nonresponse_outcome_indices) {
  eval_input = TrialEvalInput{};
  eval_input.nonresponse_outcome_indices = nonresponse_outcome_indices;
}

inline void initialize_ranked_trial_eval_input_stage(TrialEvalInput &eval_input,
                                                     std::size_t rank_width) {
  eval_input = TrialEvalInput{};
  eval_input.sequence_label_storage.assign(rank_width, -1);
  eval_input.sequence_time_storage.assign(
      rank_width, std::numeric_limits<double>::quiet_NaN());
}

inline void finalize_ranked_trial_eval_input_inline(
    TrialEvalInput &eval_input,
    const std::vector<int> *nonresponse_outcome_indices) {
  eval_input.valid = false;
  eval_input.probability_transform = TrialProbabilityTransform::Identity;
  eval_input.enforce_component_outcome_gate = false;
  eval_input.sequence_label_data = nullptr;
  eval_input.sequence_time_data = nullptr;
  eval_input.sequence_length = 0u;
  eval_input.nonresponse_outcome_indices = nonresponse_outcome_indices;

  if (eval_input.sequence_label_storage.size() !=
      eval_input.sequence_time_storage.size()) {
    return;
  }

  bool ranked_valid = true;
  bool seen_truncation = false;
  std::vector<int> seen_label_ids;
  seen_label_ids.reserve(eval_input.sequence_label_storage.size());
  std::size_t write_idx = 0u;
  for (std::size_t rank_i = 0; rank_i < eval_input.sequence_label_storage.size();
       ++rank_i) {
    const int observed_label_id = eval_input.sequence_label_storage[rank_i];
    const double rank_rt = eval_input.sequence_time_storage[rank_i];
    const bool has_label = (observed_label_id >= 0);
    const bool has_rt = std::isfinite(rank_rt);
    if (!has_label && !has_rt) {
      seen_truncation = true;
      continue;
    }
    if (seen_truncation || (has_label != has_rt) || rank_rt < 0.0 ||
        observed_label_id < 0) {
      ranked_valid = false;
      break;
    }
    if (write_idx > 0u &&
        !(rank_rt > eval_input.sequence_time_storage[write_idx - 1u])) {
      ranked_valid = false;
      break;
    }
    auto seen_it = std::lower_bound(seen_label_ids.begin(), seen_label_ids.end(),
                                    observed_label_id);
    if (seen_it != seen_label_ids.end() && *seen_it == observed_label_id) {
      ranked_valid = false;
      break;
    }
    seen_label_ids.insert(seen_it, observed_label_id);
    eval_input.sequence_label_storage[write_idx] = observed_label_id;
    eval_input.sequence_time_storage[write_idx] = rank_rt;
    ++write_idx;
  }

  if (!ranked_valid) {
    return;
  }

  eval_input.sequence_label_storage.resize(write_idx);
  eval_input.sequence_time_storage.resize(write_idx);
  if (write_idx == 0u) {
    eval_input.valid = true;
    eval_input.probability_transform = TrialProbabilityTransform::Complement;
    return;
  }

  eval_input.valid = true;
  eval_input.sequence_label_data = eval_input.sequence_label_storage.data();
  eval_input.sequence_time_data = eval_input.sequence_time_storage.data();
  eval_input.sequence_length = write_idx;
}

inline void prepare_observed_trial_eval_input(
    TrialEvalInput &eval_input, int outcome_label_id, double rt,
    const std::vector<int> *nonresponse_outcome_indices) {
  reset_trial_eval_input(eval_input, nonresponse_outcome_indices);
  eval_input.enforce_component_outcome_gate = true;
  if (outcome_label_id < 0) {
    eval_input.valid = true;
    eval_input.probability_transform = TrialProbabilityTransform::Complement;
    return;
  }
  if (!std::isfinite(rt) || rt < 0.0) {
    return;
  }
  eval_input.valid = true;
  eval_input.single_label_storage = outcome_label_id;
  eval_input.single_time_storage = rt;
  eval_input.sequence_label_data = &eval_input.single_label_storage;
  eval_input.sequence_time_data = &eval_input.single_time_storage;
  eval_input.sequence_length = 1u;
}

bool build_trial_contributions_unified(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    std::vector<TrialContributionSpec> &contributions) {
  contributions.clear();

  const bool nonresponse_mode = (eval_input.sequence_length == 0u ||
                                 eval_input.sequence_label_data == nullptr ||
                                 eval_input.sequence_time_data == nullptr);

  if (nonresponse_mode) {
    const bool use_precomputed = eval_input.nonresponse_outcome_indices &&
                                 !eval_input.nonresponse_outcome_indices->empty();
    const std::size_t n = use_precomputed
                              ? eval_input.nonresponse_outcome_indices->size()
                              : ctx.outcome_info.size();
    contributions.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      const int oi = use_precomputed ? (*eval_input.nonresponse_outcome_indices)[i]
                                     : static_cast<int>(i);
      if (oi < 0 || oi >= static_cast<int>(ctx.outcome_info.size())) {
        continue;
      }
      const uuber::OutcomeContextInfo &info =
          ctx.outcome_info[static_cast<std::size_t>(oi)];
      TrialContributionSpec spec;
      spec.intent = TrialKernelIntent::OutcomeMass;
      if (info.node_id >= 0) {
        contributions.push_back(std::move(spec));
      }
    }
    return true;
  }

  if (component_indices.empty()) {
    return true;
  }

  const std::size_t sequence_length = eval_input.sequence_length;
  if (sequence_length == 0u || !eval_input.sequence_label_data ||
      !eval_input.sequence_time_data) {
    return false;
  }

  contributions.reserve(component_indices.size());
  const TrialSequenceExecutionKind sequence_execution =
      (sequence_length > 1u) ? TrialSequenceExecutionKind::SequenceState
                             : TrialSequenceExecutionKind::LowerLayerDirect;
  const bool uses_sequence_state =
      (sequence_execution == TrialSequenceExecutionKind::SequenceState);
  for (std::size_t c = 0; c < component_indices.size(); ++c) {
    const double mix_w =
        (c < component_weights.size()) ? component_weights[c] : 0.0;
    if (!std::isfinite(mix_w) || mix_w <= 0.0) {
      continue;
    }
    const int comp_idx = component_indices[c];
    std::vector<int> outcome_indices;
    std::vector<int> node_indices;
    outcome_indices.reserve(sequence_length);
    node_indices.reserve(sequence_length);
    const bool enforce_unique_sequence_nodes = uses_sequence_state;
    std::vector<std::uint8_t> seen_node_bits(
        enforce_unique_sequence_nodes ? ctx.ir.nodes.size() : 0u, 0u);

    bool component_ok = true;
    double keep_mult = 1.0;

    for (std::size_t rank_i = 0; rank_i < sequence_length; ++rank_i) {
      const int label_id = eval_input.sequence_label_data[rank_i];
      const int comp_outcome_idx =
          resolve_outcome_index_ir(ctx, label_id, comp_idx);
      if (comp_outcome_idx < 0 ||
          comp_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
        component_ok = false;
        break;
      }
      const uuber::OutcomeContextInfo &info =
          ctx.outcome_info[static_cast<std::size_t>(comp_outcome_idx)];
      const int node_idx = resolve_dense_node_idx_required(ctx, info.node_id);
      if (node_idx < 0) {
        component_ok = false;
        break;
      }
      if (enforce_unique_sequence_nodes) {
        if (static_cast<std::size_t>(node_idx) >= seen_node_bits.size() ||
            seen_node_bits[static_cast<std::size_t>(node_idx)] != 0u) {
          component_ok = false;
          break;
        }
        seen_node_bits[static_cast<std::size_t>(node_idx)] = 1u;
      } else if (eval_input.enforce_component_outcome_gate &&
                 !outcome_allows_component_idx(ctx, info, comp_outcome_idx,
                                              comp_idx)) {
        component_ok = false;
        break;
      }

      const double keep_w = component_keep_weight(ctx, comp_idx, comp_outcome_idx);
      if (!std::isfinite(keep_w) || keep_w <= 0.0) {
        component_ok = false;
        break;
      }
      keep_mult *= keep_w;
      if (!std::isfinite(keep_mult) || keep_mult <= 0.0) {
        component_ok = false;
        break;
      }

      outcome_indices.push_back(comp_outcome_idx);
      node_indices.push_back(node_idx);
    }

    if (!component_ok) {
      continue;
    }

    const double scaled_weight = mix_w * keep_mult;
    if (!std::isfinite(scaled_weight) || scaled_weight <= 0.0) {
      continue;
    }

    TrialContributionSpec spec;
    spec.sequence_execution = sequence_execution;
    spec.scaled_weight = scaled_weight;
    spec.component_idx = comp_idx;
    if (c < cache_entries.size()) {
      spec.trial_type_key_ptr = &cache_entries[c].trial_type_key;
    } else if (eval_input.enforce_component_outcome_gate) {
      spec.trial_type_key_storage =
          component_cache_key(ctx, comp_idx, std::string());
      spec.trial_type_key_ptr = &spec.trial_type_key_storage;
    } else {
      spec.trial_type_key_ptr = &empty_trial_type_key_ref();
    }
    spec.intent = TrialKernelIntent::SequenceDensity;
    spec.sequence_outcome_indices = std::move(outcome_indices);
    spec.sequence_node_indices = std::move(node_indices);
    spec.sequence_time_data = eval_input.sequence_time_data;
    spec.step_competitor_ids_ptrs.assign(sequence_length, nullptr);
    spec.step_competitor_ids_storage.resize(sequence_length);
    spec.step_persistent_sources.resize(sequence_length);
    std::vector<std::uint8_t> observed_node_bits;
    std::vector<std::uint8_t> future_node_bits;
    if (uses_sequence_state) {
      observed_node_bits.assign(ctx.ir.nodes.size(), 0u);
      future_node_bits.assign(ctx.ir.nodes.size(), 0u);
      for (int node_idx : spec.sequence_node_indices) {
        if (node_idx >= 0 &&
            static_cast<std::size_t>(node_idx) < future_node_bits.size()) {
          future_node_bits[static_cast<std::size_t>(node_idx)] = 1u;
        }
      }
    }
    for (std::size_t rank_i = 0; rank_i < sequence_length; ++rank_i) {
      const int comp_outcome_idx = spec.sequence_outcome_indices[rank_i];
      const int node_idx = spec.sequence_node_indices[rank_i];
      if (comp_outcome_idx < 0 ||
          comp_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
        continue;
      }
      const uuber::OutcomeContextInfo &info =
          ctx.outcome_info[static_cast<std::size_t>(comp_outcome_idx)];
      std::vector<int> filtered_competitors;
      const std::vector<int> &base_competitors = filter_competitor_ids(
          ctx, info.competitor_ids, comp_idx, filtered_competitors);
      if (!uses_sequence_state) {
        if (&base_competitors == &info.competitor_ids) {
          spec.step_competitor_ids_ptrs[rank_i] = &info.competitor_ids;
        } else {
          spec.step_competitor_ids_storage[rank_i] =
              std::move(filtered_competitors);
          spec.step_competitor_ids_ptrs[rank_i] =
              &spec.step_competitor_ids_storage[rank_i];
        }
      } else {
        std::vector<int> &step_competitors =
            spec.step_competitor_ids_storage[rank_i];
        step_competitors.reserve(base_competitors.size());
        for (int node_id : base_competitors) {
          if (node_id == NA_INTEGER || node_id == info.node_id) {
            continue;
          }
          const int competitor_node_idx =
              resolve_dense_node_idx_required(ctx, node_id);
          if (competitor_node_idx >= 0 &&
              static_cast<std::size_t>(competitor_node_idx) < observed_node_bits.size() &&
              observed_node_bits[static_cast<std::size_t>(competitor_node_idx)] != 0u) {
            continue;
          }
          step_competitors.push_back(node_id);
        }
        sort_unique(step_competitors);
        spec.step_competitor_ids_ptrs[rank_i] = &step_competitors;
        std::vector<int> persistent_competitors;
        persistent_competitors.reserve(step_competitors.size());
        for (int node_id : step_competitors) {
          const int competitor_node_idx =
              resolve_dense_node_idx_required(ctx, node_id);
          if (competitor_node_idx < 0 ||
              static_cast<std::size_t>(competitor_node_idx) >= future_node_bits.size() ||
              future_node_bits[static_cast<std::size_t>(competitor_node_idx)] == 0u) {
            persistent_competitors.push_back(node_id);
          }
        }
        spec.step_persistent_sources[rank_i] =
            collect_competitor_sources(ctx, persistent_competitors);
        if (node_idx >= 0 &&
            static_cast<std::size_t>(node_idx) < observed_node_bits.size()) {
          observed_node_bits[static_cast<std::size_t>(node_idx)] = 1u;
          future_node_bits[static_cast<std::size_t>(node_idx)] = 0u;
        }
      }
    }
    contributions.push_back(std::move(spec));
  }

  return true;
}

inline bool build_unified_direct_trial_batch_spec(
    const uuber::NativeContext &ctx, const TrialContributionSpec &spec,
    SingleStepDirectBatchSpec &out) {
  out = SingleStepDirectBatchSpec{};
  if (spec.intent != TrialKernelIntent::SequenceDensity ||
      spec.sequence_execution != TrialSequenceExecutionKind::LowerLayerDirect ||
      spec.sequence_node_indices.size() != 1u ||
      spec.sequence_outcome_indices.size() != 1u ||
      spec.sequence_time_data == nullptr ||
      !std::isfinite(spec.sequence_time_data[0]) ||
      spec.sequence_time_data[0] < 0.0 ||
      !std::isfinite(spec.scaled_weight) || spec.scaled_weight <= 0.0) {
    return false;
  }

  const int node_idx = spec.sequence_node_indices[0];
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  const int stable_node_id =
      ctx.ir.nodes[static_cast<std::size_t>(node_idx)].node_id;
  out.node_id = stable_node_id >= 0 ? stable_node_id : node_idx;
  out.component_idx = spec.component_idx;
  out.outcome_idx_context = spec.sequence_outcome_indices[0];
  out.time = spec.sequence_time_data[0];
  out.scaled_weight = spec.scaled_weight;
  out.trial_type_key =
      spec.trial_type_key_ptr ? *spec.trial_type_key_ptr
                              : empty_trial_type_key_ref();
  if (!spec.step_competitor_ids_ptrs.empty() &&
      spec.step_competitor_ids_ptrs[0] != nullptr) {
    out.competitor_ids = *spec.step_competitor_ids_ptrs[0];
  }
  return true;
}

inline bool build_direct_trial_batch_spec(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    DirectTrialBatchSpec &out) {
  out = DirectTrialBatchSpec{};
  if (!eval_input.valid ||
      eval_input.probability_transform != TrialProbabilityTransform::Identity ||
      eval_input.sequence_length != 1u || eval_input.sequence_time_data == nullptr ||
      !std::isfinite(eval_input.sequence_time_data[0]) ||
      eval_input.sequence_time_data[0] < 0.0) {
    return false;
  }

  std::vector<TrialContributionSpec> contributions;
  if (!build_trial_contributions_unified(
          ctx, eval_input, component_indices, component_weights, cache_entries,
          contributions) ||
      contributions.size() <= 1u) {
    return false;
  }

  out.time = eval_input.sequence_time_data[0];
  out.contributions.reserve(contributions.size());
  for (const TrialContributionSpec &contribution : contributions) {
    SingleStepDirectBatchSpec direct_spec;
    if (!build_unified_direct_trial_batch_spec(ctx, contribution,
                                               direct_spec)) {
      return false;
    }
    DirectContributionTemplate stored;
    stored.node_id = direct_spec.node_id;
    stored.component_idx = direct_spec.component_idx;
    stored.outcome_idx_context = direct_spec.outcome_idx_context;
    stored.scaled_weight = direct_spec.scaled_weight;
    stored.trial_type_key = std::move(direct_spec.trial_type_key);
    stored.competitor_ids = std::move(direct_spec.competitor_ids);
    out.contributions.push_back(std::move(stored));
  }
  return true;
}

inline void validate_loglik_context_or_stop(
    const uuber::NativeContext &ctx) {
  if (!ctx.ir.valid) {
    Rcpp::stop("loglik requires a compiled IR context");
  }
  if (!ctx.tree_program || !ctx.tree_program->valid) {
    Rcpp::stop("loglik requires a compiled tree vector program");
  }
  if (!ctx.tree_runtime.valid) {
    Rcpp::stop("loglik requires compiled tree guard metadata");
  }
  if (ctx.outcome_label_ids.empty()) {
    Rcpp::stop("IR loglik requires native outcome labels in context");
  }
}

inline void build_loglik_workload(const uuber::NativeContext &ctx,
                                  Rcpp::DataFrame data_df,
                                  Rcpp::LogicalVector ok,
                                  Rcpp::IntegerVector expand,
                                  LogLikWorkload &out) {
  out = LogLikWorkload{};

  const ComponentMap &comp_map = ctx.components;
  const std::vector<std::string> &comp_ids = comp_map.ids;
  const std::vector<int> &outcome_label_ids = ctx.outcome_label_ids;

  auto decode_outcome_factor = [&](SEXP col_sexp) -> Rcpp::IntegerVector {
    const R_xlen_t n = Rf_xlength(col_sexp);
    Rcpp::IntegerVector decoded(n, NA_INTEGER);
    Rcpp::IntegerVector codes(col_sexp);
    for (R_xlen_t i = 0; i < n; ++i) {
      const int code = codes[i];
      if (code == NA_INTEGER) {
        continue;
      }
      decoded[i] = outcome_label_ids[static_cast<std::size_t>(code - 1)];
    }
    return decoded;
  };

  auto decode_component_factor = [&](SEXP col_sexp) -> Rcpp::IntegerVector {
    const R_xlen_t n = Rf_xlength(col_sexp);
    Rcpp::IntegerVector decoded(n, NA_INTEGER);
    Rcpp::IntegerVector codes(col_sexp);
    for (R_xlen_t i = 0; i < n; ++i) {
      const int code = codes[i];
      if (code == NA_INTEGER) {
        continue;
      }
      decoded[i] = code - 1;
    }
    return decoded;
  };

  Rcpp::NumericVector trial_col = data_df["trial"];
  out.n_rows = trial_col.size();
  out.row_to_trial_index.assign(static_cast<std::size_t>(out.n_rows), -1);
  out.row_acc_indices.assign(static_cast<std::size_t>(out.n_rows), -1);
  out.trial_eval_inputs.reserve(static_cast<std::size_t>(out.n_rows));
  out.forced_component_idx_by_trial.reserve(static_cast<std::size_t>(out.n_rows));

  const bool has_r = data_df.containsElementNamed("R");
  Rcpp::IntegerVector outcome_id_col =
      has_r ? decode_outcome_factor(data_df["R"]) : Rcpp::IntegerVector();
  Rcpp::NumericVector rt_col = data_df["rt"];
  std::vector<Rcpp::IntegerVector> rank_outcome_id_cols;
  std::vector<Rcpp::NumericVector> rank_rt_cols;
  rank_outcome_id_cols.push_back(outcome_id_col);
  rank_rt_cols.push_back(rt_col);
  for (int rank = 2;; ++rank) {
    const std::string r_name = "R" + std::to_string(rank);
    const std::string rt_name = "rt" + std::to_string(rank);
    const bool has_rank_label = data_df.containsElementNamed(r_name.c_str());
    const bool has_rank_rt = data_df.containsElementNamed(rt_name.c_str());
    if (!(has_rank_label && has_rank_rt)) {
      break;
    }
    rank_outcome_id_cols.push_back(decode_outcome_factor(data_df[r_name]));
    rank_rt_cols.push_back(Rcpp::NumericVector(data_df[rt_name]));
  }
  const int rank_width = static_cast<int>(rank_outcome_id_cols.size());
  out.ranked_mode = rank_width > 1;
  if (outcome_id_col.size() == 0) {
    Rcpp::stop("IR loglik requires factor column 'R'");
  }

  const bool has_component_label = data_df.containsElementNamed("component");
  Rcpp::IntegerVector component_idx_col =
      has_component_label
          ? decode_component_factor(data_df["component"])
          : Rcpp::IntegerVector();

  out.has_onset = data_df.containsElementNamed("onset");
  if (out.has_onset) {
    Rcpp::NumericVector onset_col(data_df["onset"]);
    out.onset_by_row.assign(onset_col.begin(), onset_col.end());
  }

  std::vector<int> all_acc_indices(ctx.accumulators.size());
  for (std::size_t i = 0; i < all_acc_indices.size(); ++i) {
    all_acc_indices[i] = static_cast<int>(i);
  }
  std::vector<std::vector<int>> comp_acc_indices(comp_ids.size());
  for (std::size_t c = 0; c < comp_ids.size(); ++c) {
    const std::string &cid = comp_ids[c];
    for (std::size_t a = 0; a < ctx.accumulators.size(); ++a) {
      const auto &comps = ctx.accumulators[a].components;
      if (std::find(comps.begin(), comps.end(), cid) != comps.end()) {
        comp_acc_indices[c].push_back(static_cast<int>(a));
      }
    }
    if (comp_acc_indices[c].empty()) {
      comp_acc_indices[c] = all_acc_indices;
    }
  }

  std::vector<int> outcome_label_id_by_trial;
  std::vector<double> rt_by_trial;
  outcome_label_id_by_trial.reserve(static_cast<std::size_t>(out.n_rows));
  rt_by_trial.reserve(static_cast<std::size_t>(out.n_rows));

  int current_trial_label = std::numeric_limits<int>::min();
  int trial_idx = -1;
  std::size_t acc_cursor = 0u;
  const std::vector<int> *current_acc_order = &all_acc_indices;

  for (R_xlen_t r = 0; r < out.n_rows; ++r) {
    const int trial_label = static_cast<int>(trial_col[r]);
    if (r == 0 || trial_label != current_trial_label) {
      current_trial_label = trial_label;
      ++trial_idx;
      out.trial_eval_inputs.push_back(TrialEvalInput{});
      out.forced_component_idx_by_trial.push_back(-1);
      outcome_label_id_by_trial.push_back(-1);
      rt_by_trial.push_back(std::numeric_limits<double>::quiet_NaN());
      if (out.ranked_mode) {
        initialize_ranked_trial_eval_input_stage(
            out.trial_eval_inputs.back(), static_cast<std::size_t>(rank_width));
      }
      acc_cursor = 0u;
      current_acc_order = &all_acc_indices;
    }

    out.row_to_trial_index[static_cast<std::size_t>(r)] = trial_idx;

    if (has_component_label &&
        out.forced_component_idx_by_trial[static_cast<std::size_t>(trial_idx)] < 0 &&
        r < component_idx_col.size()) {
      const int comp_idx_val = component_idx_col[r];
      if (comp_idx_val != NA_INTEGER) {
        if (comp_idx_val < 0 ||
            comp_idx_val >= static_cast<int>(comp_acc_indices.size())) {
          Rcpp::stop("IR loglik encountered an out-of-range component index");
        }
        out.forced_component_idx_by_trial[static_cast<std::size_t>(trial_idx)] =
            comp_idx_val;
        current_acc_order =
            &comp_acc_indices[static_cast<std::size_t>(comp_idx_val)];
      }
    }

    const int acc_idx =
        (acc_cursor < current_acc_order->size())
            ? (*current_acc_order)[acc_cursor]
            : all_acc_indices[acc_cursor % all_acc_indices.size()];
    out.row_acc_indices[static_cast<std::size_t>(r)] = acc_idx;
    ++acc_cursor;

    if (outcome_label_id_by_trial[static_cast<std::size_t>(trial_idx)] < 0 &&
        r < outcome_id_col.size()) {
      const int id_val = outcome_id_col[r];
      if (id_val != NA_INTEGER) {
        outcome_label_id_by_trial[static_cast<std::size_t>(trial_idx)] = id_val;
      }
    }
    if (!std::isfinite(rt_by_trial[static_cast<std::size_t>(trial_idx)]) &&
        r < rt_col.size()) {
      rt_by_trial[static_cast<std::size_t>(trial_idx)] =
          static_cast<double>(rt_col[r]);
    }

    if (out.ranked_mode) {
      TrialEvalInput &trial_eval_input =
          out.trial_eval_inputs[static_cast<std::size_t>(trial_idx)];
      std::vector<int> &trial_ranked_ids =
          trial_eval_input.sequence_label_storage;
      std::vector<double> &trial_ranked_rt =
          trial_eval_input.sequence_time_storage;
      for (int rank_idx = 0; rank_idx < rank_width; ++rank_idx) {
        if (trial_ranked_ids[static_cast<std::size_t>(rank_idx)] < 0 &&
            r < rank_outcome_id_cols[static_cast<std::size_t>(rank_idx)].size()) {
          const int id_val =
              rank_outcome_id_cols[static_cast<std::size_t>(rank_idx)][r];
          if (id_val != NA_INTEGER) {
            trial_ranked_ids[static_cast<std::size_t>(rank_idx)] = id_val;
          }
        }
        if (!std::isfinite(
                trial_ranked_rt[static_cast<std::size_t>(rank_idx)]) &&
            r < rank_rt_cols[static_cast<std::size_t>(rank_idx)].size()) {
          const double cand_rt =
              rank_rt_cols[static_cast<std::size_t>(rank_idx)][r];
          if (std::isfinite(cand_rt)) {
            trial_ranked_rt[static_cast<std::size_t>(rank_idx)] = cand_rt;
          }
        }
      }
    }
  }

  out.n_trials = static_cast<int>(out.trial_eval_inputs.size());

  out.nonresponse_outcome_indices.reserve(ctx.outcome_info.size());
  for (std::size_t oi = 0; oi < ctx.outcome_info.size(); ++oi) {
    bool is_deadline = false;
    const bool maps_to_na =
        ctx.outcome_info[static_cast<std::size_t>(oi)].maps_to_na;
    if (oi < ctx.ir.outcomes.size()) {
      const int node_idx = ctx.ir.outcomes[oi].node_idx;
      if (node_idx >= 0 && node_idx < static_cast<int>(ctx.ir.nodes.size())) {
        const auto flags =
            ctx.ir.nodes[static_cast<std::size_t>(node_idx)].flags;
        is_deadline = (flags & uuber::IR_NODE_FLAG_SPECIAL_DEADLINE) != 0u;
      }
    }
    if (!is_deadline && !maps_to_na) {
      out.nonresponse_outcome_indices.push_back(static_cast<int>(oi));
    }
  }

  for (int t = 0; t < out.n_trials; ++t) {
    TrialEvalInput &eval_input =
        out.trial_eval_inputs[static_cast<std::size_t>(t)];
    if (out.ranked_mode) {
      finalize_ranked_trial_eval_input_inline(eval_input,
                                              &out.nonresponse_outcome_indices);
      continue;
    }
    const int outcome_label_id =
        (t < static_cast<int>(outcome_label_id_by_trial.size()))
            ? outcome_label_id_by_trial[static_cast<std::size_t>(t)]
            : -1;
    const double rt =
        (t < static_cast<int>(rt_by_trial.size()))
            ? rt_by_trial[static_cast<std::size_t>(t)]
            : std::numeric_limits<double>::quiet_NaN();
    prepare_observed_trial_eval_input(eval_input, outcome_label_id, rt,
                                      &out.nonresponse_outcome_indices);
  }

  out.default_component_indices.reserve(comp_map.ids.size());
  for (const auto &comp_id : comp_map.ids) {
    out.default_component_indices.push_back(component_index_of(ctx, comp_id));
  }
  out.default_cache_entries = build_component_cache_entries_from_indices(
      ctx, out.default_component_indices);

  out.leader_components_by_acc_idx.assign(ctx.accumulators.size(), {});
  for (std::size_t comp_idx = 0; comp_idx < comp_map.leader_idx.size();
       ++comp_idx) {
    const int acc_idx = comp_map.leader_idx[comp_idx];
    if (acc_idx >= 0 &&
        acc_idx < static_cast<int>(out.leader_components_by_acc_idx.size())) {
      out.leader_components_by_acc_idx[static_cast<std::size_t>(acc_idx)]
          .push_back(static_cast<int>(comp_idx));
    }
  }

  out.ok_by_trial.assign(static_cast<std::size_t>(out.n_trials), 1u);
  if (ok.size() != 0 && ok.size() != out.n_trials) {
    Rcpp::stop("Length of ok must equal number of trials");
  }
  if (ok.size() != 0) {
    for (int t = 0; t < out.n_trials; ++t) {
      out.ok_by_trial[static_cast<std::size_t>(t)] =
          (ok[static_cast<std::size_t>(t)] == TRUE) ? 1u : 0u;
    }
  }

  out.expand_trial_indices.clear();
  if (expand.size() == 0) {
    out.expand_trial_indices.reserve(static_cast<std::size_t>(out.n_trials));
    for (int t = 0; t < out.n_trials; ++t) {
      out.expand_trial_indices.push_back(t);
    }
  } else {
    out.expand_trial_indices.reserve(static_cast<std::size_t>(expand.size()));
    for (R_xlen_t i = 0; i < expand.size(); ++i) {
      const int trial_idx_1based = expand[i];
      if (trial_idx_1based == NA_INTEGER || trial_idx_1based < 1 ||
          trial_idx_1based > out.n_trials) {
        Rcpp::stop("expand indices must reference trials");
      }
      out.expand_trial_indices.push_back(trial_idx_1based - 1);
    }
  }

  out.enable_ranked_sequence_batch = false;
  for (int t = 0; t < out.n_trials; ++t) {
    if (out.ok_by_trial[static_cast<std::size_t>(t)] == 0u) {
      continue;
    }
    const TrialEvalInput &eval_input =
        out.trial_eval_inputs[static_cast<std::size_t>(t)];
    if (eval_input.valid &&
        eval_input.probability_transform ==
            TrialProbabilityTransform::Identity &&
        eval_input.sequence_length > 1u &&
        eval_input.sequence_time_data != nullptr) {
      out.enable_ranked_sequence_batch = true;
      break;
    }
  }
}

inline void build_loglik_prepared_params_batch(
    const uuber::NativeContext &ctx, const LogLikWorkload &workload,
    const std::vector<Rcpp::NumericMatrix> &params_mats,
    std::vector<LogLikPreparedParams> &out) {
  out.clear();
  if (params_mats.empty()) {
    return;
  }

  struct ParamMatrixLayout {
    const Rcpp::NumericMatrix *params_mat{nullptr};
    std::array<int, 8> p_cols{};
  };

  const int q_col = 0;
  const int w_col = 1;
  const int t0_col = 2;
  const std::size_t n_trials = static_cast<std::size_t>(workload.n_trials);
  const std::size_t comp_count = workload.default_component_indices.size();
  const std::vector<TrialAccumulatorParams> base_acc_params =
      build_base_paramset(ctx).acc_params;

  std::vector<ParamMatrixLayout> layouts;
  layouts.reserve(params_mats.size());
  for (const Rcpp::NumericMatrix &params_mat : params_mats) {
    if (params_mat.nrow() != workload.n_rows) {
      Rcpp::stop(
          "Parameter matrix rows (%lld) must match nested likelihood rows (%lld)",
          static_cast<long long>(params_mat.nrow()),
          static_cast<long long>(workload.n_rows));
    }
    ParamMatrixLayout layout;
    layout.params_mat = &params_mat;
    layout.p_cols.fill(-1);
    const int n_param_cols = std::max(0, params_mat.ncol() - 3);
    const int n_slot_cols = std::min(8, n_param_cols);
    for (int i = 0; i < n_slot_cols; ++i) {
      layout.p_cols[static_cast<std::size_t>(i)] = 3 + i;
    }
    layouts.push_back(layout);
  }

  out.resize(params_mats.size());
  for (LogLikPreparedParams &prepared : out) {
    prepared.param_sets.resize(n_trials);
    for (TrialParamSet &params : prepared.param_sets) {
      params.acc_params = base_acc_params;
    }
    prepared.prepared_trial_runtimes.resize(n_trials);
    prepared.weights_by_trial.assign(
        n_trials,
        std::vector<double>(comp_count, std::numeric_limits<double>::quiet_NaN()));
    prepared.prepared_trial_runtime_ptrs.clear();
    prepared.trigger_plans.clear();
    prepared.trigger_plan_index_by_trial.clear();
    prepared.param_value_group_by_trial.clear();
    prepared.param_value_group_representative_trial.clear();
    prepared.component_weight_group_by_trial.clear();
    prepared.component_weight_group_representative_trial.clear();
  }

  for (R_xlen_t r = 0; r < workload.n_rows; ++r) {
    const int trial_idx = workload.row_to_trial_index[static_cast<std::size_t>(r)];
    if (trial_idx < 0 || trial_idx >= workload.n_trials) {
      Rcpp::stop("loglik workload trial mapping is invalid");
    }
    const int acc_idx = workload.row_acc_indices[static_cast<std::size_t>(r)];
    if (acc_idx < 0 || acc_idx >= static_cast<int>(base_acc_params.size())) {
      Rcpp::stop("loglik workload accumulator mapping is invalid");
    }

    const bool has_onset =
        workload.has_onset && static_cast<std::size_t>(r) < workload.onset_by_row.size();
    const double onset_val =
        has_onset ? workload.onset_by_row[static_cast<std::size_t>(r)]
                  : std::numeric_limits<double>::quiet_NaN();
    const std::vector<int> *leader_components_ptr =
        (acc_idx >= 0 &&
         acc_idx < static_cast<int>(workload.leader_components_by_acc_idx.size()))
            ? &workload.leader_components_by_acc_idx[static_cast<std::size_t>(
                  acc_idx)]
            : nullptr;

    for (std::size_t batch_idx = 0; batch_idx < out.size(); ++batch_idx) {
      LogLikPreparedParams &prepared = out[batch_idx];
      TrialParamSet &trial_params =
          prepared.param_sets[static_cast<std::size_t>(trial_idx)];
      TrialAccumulatorParams &tap =
          trial_params.acc_params[static_cast<std::size_t>(acc_idx)];
      const Rcpp::NumericMatrix &params_mat = *layouts[batch_idx].params_mat;

      const double q_val = clamp_probability(params_mat(r, q_col));
      if (!tap.shared_trigger_id.empty()) {
        tap.shared_q = q_val;
        tap.q = 0.0;
      } else {
        tap.q = q_val;
        tap.shared_q = q_val;
      }
      tap.dist_cfg.t0 = params_mat(r, t0_col);
      if (has_onset && std::isfinite(onset_val)) {
        tap.onset = onset_val;
      }
      tap.has_override = true;
      trial_params.has_any_override = true;
      trial_params.soa_cache_valid = false;

      const std::array<int, 8> &p_cols = layouts[batch_idx].p_cols;
      if (p_cols[0] >= 0) {
        tap.dist_cfg.p1 = params_mat(r, p_cols[0]);
      }
      if (p_cols[1] >= 0) {
        tap.dist_cfg.p2 = params_mat(r, p_cols[1]);
      }
      if (p_cols[2] >= 0) {
        tap.dist_cfg.p3 = params_mat(r, p_cols[2]);
      }
      if (p_cols[3] >= 0) {
        tap.dist_cfg.p4 = params_mat(r, p_cols[3]);
      }
      if (p_cols[4] >= 0) {
        tap.dist_cfg.p5 = params_mat(r, p_cols[4]);
      }
      if (p_cols[5] >= 0) {
        tap.dist_cfg.p6 = params_mat(r, p_cols[5]);
      }
      if (p_cols[6] >= 0) {
        tap.dist_cfg.p7 = params_mat(r, p_cols[6]);
      }
      if (p_cols[7] >= 0) {
        tap.dist_cfg.p8 = params_mat(r, p_cols[7]);
      }

      if (leader_components_ptr != nullptr) {
        std::vector<double> &weight_overrides =
            prepared.weights_by_trial[static_cast<std::size_t>(trial_idx)];
        for (int comp_idx : *leader_components_ptr) {
          if (comp_idx >= 0 &&
              static_cast<std::size_t>(comp_idx) < weight_overrides.size()) {
            weight_overrides[static_cast<std::size_t>(comp_idx)] =
                params_mat(r, w_col);
          }
        }
      }
    }
  }

  std::unordered_map<uuber::NAMapCacheKey, PreparedTrialParamsRuntime *,
                     uuber::NAMapCacheKeyHash>
      canonical_runtime_by_value_key;
  canonical_runtime_by_value_key.reserve(n_trials * out.size());
  std::unordered_map<uuber::NAMapCacheKey, SharedTriggerPlan,
                     uuber::NAMapCacheKeyHash>
      trigger_plan_by_source_key;
  trigger_plan_by_source_key.reserve(n_trials * out.size());

  for (LogLikPreparedParams &prepared : out) {
    for (std::size_t t = 0; t < n_trials; ++t) {
      TrialParamSet &trial_params = prepared.param_sets[t];
      PreparedTrialParamsRuntime &runtime = prepared.prepared_trial_runtimes[t];
      if (!prepare_trial_params_runtime(ctx, &trial_params, runtime)) {
        Rcpp::stop("Failed to prepare trial runtime metadata");
      }

      std::vector<double> &weights = prepared.weights_by_trial[t];
      if (weights.empty()) {
        continue;
      }
      double total = 0.0;
      for (std::size_t c = 0; c < weights.size(); ++c) {
        if (!std::isfinite(weights[c])) {
          weights[c] = ctx.components.base_weights[c];
        }
        total += weights[c];
      }
      if (!std::isfinite(total) || total <= 0.0) {
        const double uniform = 1.0 / static_cast<double>(weights.size());
        std::fill(weights.begin(), weights.end(), uniform);
      } else {
        for (double &value : weights) {
          value /= total;
        }
      }

      auto soa_insert = canonical_runtime_by_value_key.emplace(runtime.value_key,
                                                               &runtime);
      if (soa_insert.second) {
        if (!ensure_prepared_trial_params_soa(ctx, runtime)) {
          Rcpp::stop("Failed to prepare trial runtime SoA");
        }
      } else {
        PreparedTrialParamsRuntime *canonical_runtime = soa_insert.first->second;
        if (canonical_runtime == nullptr ||
            !ensure_prepared_trial_params_soa(ctx, *canonical_runtime)) {
          Rcpp::stop("Failed to resolve canonical trial runtime SoA");
        }
        runtime.soa = canonical_runtime->soa;
      }

      auto trigger_insert =
          trigger_plan_by_source_key.emplace(runtime.shared_trigger_source_key,
                                            SharedTriggerPlan{});
      if (trigger_insert.second) {
        if (!ensure_prepared_trial_params_trigger_plan(ctx, runtime)) {
          Rcpp::stop("Failed to prepare shared-trigger plan");
        }
        trigger_insert.first->second = runtime.trigger_plan;
      } else {
        runtime.trigger_plan = trigger_insert.first->second;
        runtime.trigger_plan_ready = true;
      }
    }
  }
}

struct LogLikPreparedBatchPoint {
  int param_batch_idx{-1};
  int trial_idx{-1};
};

inline Rcpp::NumericVector execute_loglik_workload_multi(
    uuber::NativeContext &ctx, const LogLikWorkload &workload,
    std::vector<LogLikPreparedParams> &prepared_batch, double min_ll,
    double rel_tol, double abs_tol, int max_depth) {
  ctx.na_map_cache.clear();
  ctx.na_cache_order.clear();

  const int n_trials = workload.n_trials;
  const int n_param_batches = static_cast<int>(prepared_batch.size());
  Rcpp::NumericVector total_loglik(n_param_batches, 0.0);
  if (n_param_batches == 0) {
    return total_loglik;
  }

  const std::vector<TrialEvalInput> &trial_eval_inputs =
      workload.trial_eval_inputs;
  const std::vector<int> &default_component_indices =
      workload.default_component_indices;
  const std::vector<ComponentCacheEntry> &default_cache_entries =
      workload.default_cache_entries;
  const std::vector<int> &forced_component_idx_by_trial =
      workload.forced_component_idx_by_trial;
  const std::vector<std::uint8_t> &ok_by_trial = workload.ok_by_trial;
  const std::vector<int> &expand_trial_indices = workload.expand_trial_indices;

  for (int param_batch_idx = 0; param_batch_idx < n_param_batches;
       ++param_batch_idx) {
    const LogLikPreparedParams &prepared =
        prepared_batch[static_cast<std::size_t>(param_batch_idx)];
    if (prepared.param_sets.size() != static_cast<std::size_t>(n_trials) ||
        prepared.prepared_trial_runtimes.size() !=
            static_cast<std::size_t>(n_trials) ||
        prepared.weights_by_trial.size() !=
            static_cast<std::size_t>(n_trials)) {
      Rcpp::stop("multi-loglik prepared params do not match workload trials");
    }
  }

  struct ForcedComponentBundle {
    std::vector<int> component_indices;
    std::vector<double> component_weights;
    std::vector<ComponentCacheEntry> cache_entries;
  };
  std::unordered_map<int, ForcedComponentBundle> forced_component_bundle_cache;
  forced_component_bundle_cache.reserve(static_cast<std::size_t>(n_trials));
  auto forced_bundle_for =
      [&](int forced_component_idx) -> const ForcedComponentBundle & {
    auto it = forced_component_bundle_cache.find(forced_component_idx);
    if (it != forced_component_bundle_cache.end()) {
      return it->second;
    }
    ForcedComponentBundle bundle;
    bundle.component_indices.push_back(forced_component_idx);
    bundle.component_weights.push_back(1.0);
    bundle.cache_entries.push_back(
        default_component_cache_entry_idx(ctx, forced_component_idx));
    auto inserted = forced_component_bundle_cache.emplace(forced_component_idx,
                                                          std::move(bundle));
    return inserted.first->second;
  };

  struct SharedTriggerMaskBatchCacheEntry {
    bool ready{false};
    bool valid{false};
    SharedTriggerMaskSoABatch batch;
  };

  auto flatten_point_index =
      [n_trials](int param_batch_idx, int trial_idx) -> std::size_t {
    return static_cast<std::size_t>(param_batch_idx) *
               static_cast<std::size_t>(n_trials) +
           static_cast<std::size_t>(trial_idx);
  };
  auto flat_idx_of_point =
      [&](const LogLikPreparedBatchPoint &point) -> std::size_t {
    return flatten_point_index(point.param_batch_idx, point.trial_idx);
  };
  auto first_flat_idx_of_group =
      [&](const std::vector<LogLikPreparedBatchPoint> &points) -> std::size_t {
    return flat_idx_of_point(points.front());
  };

  auto assign_group_loglik =
      [&](const std::vector<LogLikPreparedBatchPoint> &points,
          const std::vector<double> &probabilities,
          std::vector<double> &trial_loglik,
          std::vector<std::uint8_t> &trial_loglik_batched) -> void {
    for (std::size_t i = 0; i < points.size(); ++i) {
      const LogLikPreparedBatchPoint &point = points[i];
      const std::size_t flat_idx = flat_idx_of_point(point);
      const double prob = (i < probabilities.size()) ? probabilities[i] : 0.0;
      double ll_val = min_ll;
      if (std::isfinite(prob) && prob > 0.0) {
        const double lp = std::log(prob);
        if (std::isfinite(lp) && lp > min_ll) {
          ll_val = lp;
        }
      }
      trial_loglik[flat_idx] = ll_val;
      trial_loglik_batched[flat_idx] = 1u;
    }
  };

  const std::size_t total_point_count =
      static_cast<std::size_t>(n_trials) *
      static_cast<std::size_t>(n_param_batches);
  std::unordered_map<uuber::NAMapCacheKey, int, uuber::NAMapCacheKeyHash>
      global_param_group_by_value_key;
  global_param_group_by_value_key.reserve(total_point_count);
  auto global_param_group_for =
      [&](const PreparedTrialParamsRuntime &runtime) -> int {
    const auto insert_result = global_param_group_by_value_key.emplace(
        runtime.value_key,
        static_cast<int>(global_param_group_by_value_key.size()));
    return insert_result.first->second;
  };
  std::unordered_map<uuber::NAMapCacheKey, SharedTriggerMaskBatchCacheEntry,
                     uuber::NAMapCacheKeyHash>
      shared_trigger_mask_batch_by_value_key;
  shared_trigger_mask_batch_by_value_key.reserve(total_point_count);
  auto resolve_shared_trigger_mask_batch =
      [&](PreparedTrialParamsRuntime &runtime)
      -> const SharedTriggerMaskSoABatch * {
    if (!ensure_prepared_trial_params_trigger_plan(ctx, runtime) ||
        shared_trigger_count(runtime.trigger_plan) == 0u) {
      return nullptr;
    }
    SharedTriggerMaskBatchCacheEntry &entry =
        shared_trigger_mask_batch_by_value_key[runtime.value_key];
    if (!entry.ready) {
      entry.ready = true;
      entry.valid = build_shared_trigger_mask_soa_batch(
                        ctx, runtime.params, runtime.trigger_plan, entry.batch) &&
                    !entry.batch.mask_param_ptrs.empty();
      if (!entry.valid) {
        entry.batch = SharedTriggerMaskSoABatch{};
      }
    }
    return entry.valid ? &entry.batch : nullptr;
  };
  auto prepared_params_for_point =
      [&](const LogLikPreparedBatchPoint &point) -> LogLikPreparedParams & {
    return prepared_batch[static_cast<std::size_t>(point.param_batch_idx)];
  };
  auto prepared_runtime_for_point =
      [&](const LogLikPreparedBatchPoint &point)
      -> PreparedTrialParamsRuntime & {
    return prepared_params_for_point(point).prepared_trial_runtimes
        [static_cast<std::size_t>(point.trial_idx)];
  };
  auto for_each_prepared_query_lane =
      [&](const std::vector<LogLikPreparedBatchPoint> &points,
          bool shared_trigger_supported, const TrialParamSet **rep_params_out,
          auto &&visitor) -> bool {
    std::unordered_map<uuber::NAMapCacheKey, const uuber::TrialParamsSoA *,
                       uuber::NAMapCacheKeyHash>
        canonical_soa_by_value_key;
    canonical_soa_by_value_key.reserve(points.size());
    for (std::size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
      const LogLikPreparedBatchPoint &point = points[point_idx];
      PreparedTrialParamsRuntime &runtime = prepared_runtime_for_point(point);
      if (!ensure_prepared_trial_params_soa(ctx, runtime) || runtime.soa == nullptr) {
        return false;
      }
      if (rep_params_out != nullptr && *rep_params_out == nullptr) {
        *rep_params_out = runtime.params;
      }
      const SharedTriggerMaskSoABatch *mask_batch_ptr =
          shared_trigger_supported ? resolve_shared_trigger_mask_batch(runtime)
                                   : nullptr;
      if (mask_batch_ptr != nullptr) {
        for (std::size_t mask_idx = 0;
             mask_idx < mask_batch_ptr->mask_param_ptrs.size(); ++mask_idx) {
          const uuber::TrialParamsSoA *mask_soa =
              mask_batch_ptr->mask_param_ptrs[mask_idx];
          const double mask_weight =
              (mask_idx < mask_batch_ptr->mask_weights.size())
                  ? mask_batch_ptr->mask_weights[mask_idx]
                  : 0.0;
          if (mask_soa == nullptr || !mask_soa->valid ||
              !std::isfinite(mask_weight) || mask_weight <= 0.0) {
            continue;
          }
          if (!visitor(point_idx, point, runtime, mask_soa, mask_weight, true)) {
            return false;
          }
        }
        continue;
      }
      const auto insert_result =
          canonical_soa_by_value_key.emplace(runtime.value_key, runtime.soa);
      if (!visitor(point_idx, point, runtime, insert_result.first->second, 1.0,
                   false)) {
        return false;
      }
    }
    return true;
  };
  std::vector<double> trial_loglik(total_point_count, min_ll);
  std::vector<std::uint8_t> trial_loglik_batched(total_point_count, 0u);

  std::vector<SingleStepDirectBatchSpec> single_step_batch_specs(total_point_count);
  std::unordered_map<SingleStepDirectBatchKey,
                     std::vector<LogLikPreparedBatchPoint>,
                     SingleStepDirectBatchKeyHash>
      single_step_batch_groups;
  single_step_batch_groups.reserve(total_point_count);

  std::vector<DirectTrialBatchSpec> direct_trial_batch_specs(total_point_count);
  std::unordered_map<DirectTrialBatchKey, std::vector<LogLikPreparedBatchPoint>,
                     DirectTrialBatchKeyHash>
      direct_trial_batch_groups;
  direct_trial_batch_groups.reserve(total_point_count);

  std::vector<RankedSequenceBatchSpec> ranked_sequence_batch_specs(
      workload.enable_ranked_sequence_batch ? total_point_count : 0u);
  std::unordered_map<RankedSequenceBatchKey,
                     std::vector<LogLikPreparedBatchPoint>,
                     RankedSequenceBatchKeyHash>
      ranked_sequence_batch_groups;
  if (workload.enable_ranked_sequence_batch) {
    ranked_sequence_batch_groups.reserve(total_point_count);
  }

  std::unordered_map<const RankedSequenceBatchTemplate *, bool>
      ranked_shared_trigger_support_by_template;
  auto ranked_template_supports_shared_trigger =
      [&](const RankedSequenceBatchTemplate *batch_template) -> bool {
    if (batch_template == nullptr) {
      return false;
    }
    auto it = ranked_shared_trigger_support_by_template.find(batch_template);
    if (it != ranked_shared_trigger_support_by_template.end()) {
      return it->second;
    }
    bool supported = true;
    for (const RankedSequenceContributionTemplate &contribution :
         batch_template->contributions) {
      for (int outcome_idx : contribution.sequence_outcome_indices) {
        if (!shared_trigger_mask_batch_supported(ctx, false, outcome_idx)) {
          supported = false;
          break;
        }
      }
      if (!supported) {
        break;
      }
    }
    ranked_shared_trigger_support_by_template.emplace(batch_template, supported);
    return supported;
  };

  std::unordered_map<OutcomeMassBatchKey, std::vector<LogLikPreparedBatchPoint>,
                     OutcomeMassBatchKeyHash>
      outcome_mass_batch_groups;
  outcome_mass_batch_groups.reserve(total_point_count);

  auto resolve_point_components =
      [&](const LogLikPreparedBatchPoint &point,
          const std::vector<int> *&component_indices_ptr,
          const std::vector<ComponentCacheEntry> *&cache_entries_ptr,
          const std::vector<double> *&component_weights_ptr) -> void {
    const LogLikPreparedParams &prepared =
        prepared_batch[static_cast<std::size_t>(point.param_batch_idx)];
    component_indices_ptr = &default_component_indices;
    cache_entries_ptr = &default_cache_entries;
    component_weights_ptr =
        &prepared.weights_by_trial[static_cast<std::size_t>(point.trial_idx)];
    const int forced_component_idx =
        (point.trial_idx < static_cast<int>(forced_component_idx_by_trial.size()))
            ? forced_component_idx_by_trial[static_cast<std::size_t>(
                  point.trial_idx)]
            : -1;
    if (forced_component_idx >= 0) {
      const ForcedComponentBundle &forced_bundle =
          forced_bundle_for(forced_component_idx);
      component_indices_ptr = &forced_bundle.component_indices;
      cache_entries_ptr = &forced_bundle.cache_entries;
      component_weights_ptr = &forced_bundle.component_weights;
    }
  };

  for (int param_batch_idx = 0; param_batch_idx < n_param_batches;
       ++param_batch_idx) {
    for (int trial_idx = 0; trial_idx < n_trials; ++trial_idx) {
      if (trial_idx >= static_cast<int>(ok_by_trial.size()) ||
          ok_by_trial[static_cast<std::size_t>(trial_idx)] == 0u) {
        continue;
      }

      const LogLikPreparedBatchPoint point{param_batch_idx, trial_idx};
      const std::size_t flat_idx = flat_idx_of_point(point);
      const TrialEvalInput &eval_input =
          trial_eval_inputs[static_cast<std::size_t>(trial_idx)];
      const LogLikPreparedParams &prepared =
          prepared_batch[static_cast<std::size_t>(param_batch_idx)];
      const PreparedTrialParamsRuntime &runtime =
          prepared.prepared_trial_runtimes[static_cast<std::size_t>(trial_idx)];
      const int global_param_group_idx = global_param_group_for(runtime);

      const std::vector<int> *component_indices_ptr = nullptr;
      const std::vector<ComponentCacheEntry> *cache_entries_ptr = nullptr;
      const std::vector<double> *component_weights_ptr = nullptr;
      resolve_point_components(point, component_indices_ptr, cache_entries_ptr,
                               component_weights_ptr);

      SingleStepDirectBatchSpec single_step_spec;
      if (build_single_step_direct_batch_spec(
              ctx, eval_input, *component_indices_ptr, *component_weights_ptr,
              *cache_entries_ptr, single_step_spec)) {
        single_step_spec.param_group_idx = global_param_group_idx;
        single_step_spec.trigger_plan_idx = -1;
        single_step_batch_specs[flat_idx] = single_step_spec;

        SingleStepDirectBatchKey key;
        key.node_id = single_step_spec.node_id;
        key.component_idx = single_step_spec.component_idx;
        key.outcome_idx_context = single_step_spec.outcome_idx_context;
        key.param_group_idx = global_param_group_idx;
        key.trigger_plan_idx = -1;
        key.trial_type_key = single_step_spec.trial_type_key;
        key.competitor_ids = single_step_spec.competitor_ids;
        single_step_batch_groups[std::move(key)].push_back(point);
        continue;
      }

      DirectTrialBatchSpec direct_trial_spec;
      if (build_direct_trial_batch_spec(
              ctx, eval_input, *component_indices_ptr, *component_weights_ptr,
              *cache_entries_ptr, direct_trial_spec)) {
        direct_trial_spec.param_group_idx = global_param_group_idx;
        direct_trial_spec.trigger_plan_idx = -1;
        direct_trial_batch_specs[flat_idx] = std::move(direct_trial_spec);

        DirectTrialBatchKey key;
        key.contributions = direct_trial_batch_specs[flat_idx].contributions;
        key.param_group_idx = global_param_group_idx;
        key.trigger_plan_idx = -1;
        direct_trial_batch_groups[std::move(key)].push_back(point);
        continue;
      }

      if (!workload.enable_ranked_sequence_batch) {
        continue;
      }

      RankedSequenceBatchSpec ranked_spec;
      if (!build_ranked_sequence_batch_spec(
              ctx, eval_input, *component_indices_ptr, *component_weights_ptr,
              *cache_entries_ptr, ranked_spec)) {
        continue;
      }
      ranked_spec.param_group_idx = global_param_group_idx;
      ranked_spec.trigger_plan_idx = -1;
      ranked_sequence_batch_specs[flat_idx] = std::move(ranked_spec);

      RankedSequenceBatchKey key;
      key.template_ptr = ranked_sequence_batch_specs[flat_idx].template_ptr;
      key.param_group_idx = global_param_group_idx;
      key.trigger_plan_idx = -1;
      ranked_sequence_batch_groups[std::move(key)].push_back(point);
    }
  }

  for (int param_batch_idx = 0; param_batch_idx < n_param_batches;
       ++param_batch_idx) {
    for (int trial_idx = 0; trial_idx < n_trials; ++trial_idx) {
      if (trial_idx >= static_cast<int>(ok_by_trial.size()) ||
          ok_by_trial[static_cast<std::size_t>(trial_idx)] == 0u) {
        continue;
      }

      const LogLikPreparedBatchPoint point{param_batch_idx, trial_idx};
      const TrialEvalInput &eval_input =
          trial_eval_inputs[static_cast<std::size_t>(trial_idx)];
      const std::vector<int> *component_indices_ptr = nullptr;
      const std::vector<ComponentCacheEntry> *cache_entries_ptr = nullptr;
      const std::vector<double> *component_weights_ptr = nullptr;
      resolve_point_components(point, component_indices_ptr, cache_entries_ptr,
                               component_weights_ptr);
      (void)cache_entries_ptr;
      (void)component_weights_ptr;

      OutcomeMassBatchKey key;
      if (!build_outcome_mass_batch_key(eval_input, *component_indices_ptr, key)) {
        continue;
      }
      outcome_mass_batch_groups[std::move(key)].push_back(point);
    }
  }

  for (auto &kv : single_step_batch_groups) {
    const std::vector<LogLikPreparedBatchPoint> &points = kv.second;
    if (points.empty()) {
      continue;
    }
    const std::size_t first_flat_idx = first_flat_idx_of_group(points);
    const SingleStepDirectBatchSpec &group_spec =
        single_step_batch_specs[first_flat_idx];
    const bool shared_trigger_supported =
        shared_trigger_density_batch_supported(ctx, false,
                                              group_spec.outcome_idx_context);

    std::vector<double> query_times;
    std::vector<const uuber::TrialParamsSoA *> query_params_soa;
    std::vector<double> query_weights;
    std::vector<std::size_t> query_to_point_index;
    query_times.reserve(points.size());
    query_params_soa.reserve(points.size());
    query_weights.reserve(points.size());
    query_to_point_index.reserve(points.size());

    const TrialParamSet *rep_params = nullptr;
    if (!for_each_prepared_query_lane(
            points, shared_trigger_supported, &rep_params,
            [&](std::size_t point_idx, const LogLikPreparedBatchPoint &point,
                PreparedTrialParamsRuntime &, const uuber::TrialParamsSoA *soa,
                double expansion_weight, bool) -> bool {
              const SingleStepDirectBatchSpec &point_spec =
                  single_step_batch_specs[flatten_point_index(
                      point.param_batch_idx, point.trial_idx)];
              query_times.push_back(point_spec.time);
              query_params_soa.push_back(soa);
              query_weights.push_back(point_spec.scaled_weight *
                                      expansion_weight);
              query_to_point_index.push_back(point_idx);
              return true;
            })) {
      Rcpp::stop("Direct multi-parameter group preparation failed");
    }

    std::vector<double> point_probabilities(points.size(), 0.0);
    if (!query_times.empty()) {
      std::vector<double> density_out;
      DensityBatchWorkspace density_workspace;
      const bool batched = node_density_entry_batch_idx(
          ctx, group_spec.node_id, query_times, group_spec.component_idx, nullptr,
          false, nullptr, false, group_spec.competitor_ids, rep_params,
          group_spec.trial_type_key, false, group_spec.outcome_idx_context, nullptr,
          false, density_out, nullptr, &query_params_soa, nullptr,
          &density_workspace);
      if (!batched || density_out.size() != query_times.size()) {
        Rcpp::stop("Direct multi-parameter group batch evaluation failed");
      }
      for (std::size_t query_idx = 0; query_idx < density_out.size(); ++query_idx) {
        const double density = density_out[query_idx];
        const double weight = query_weights[query_idx];
        if (!(std::isfinite(density) && density > 0.0 && std::isfinite(weight) &&
              weight > 0.0)) {
          continue;
        }
        point_probabilities[query_to_point_index[query_idx]] += weight * density;
      }
    }
    for (double &prob : point_probabilities) {
      prob = safe_density(prob);
    }
    assign_group_loglik(points, point_probabilities, trial_loglik,
                        trial_loglik_batched);
  }

  for (auto &kv : direct_trial_batch_groups) {
    const std::vector<LogLikPreparedBatchPoint> &points = kv.second;
    if (points.empty()) {
      continue;
    }
    const std::size_t first_flat_idx = first_flat_idx_of_group(points);
    const DirectTrialBatchSpec &group_spec = direct_trial_batch_specs[first_flat_idx];

    std::vector<double> point_probabilities(points.size(), 0.0);
    DensityBatchWorkspace density_workspace;
    for (const DirectContributionTemplate &contribution : group_spec.contributions) {
      const bool shared_trigger_supported =
          shared_trigger_density_batch_supported(ctx, false,
                                                contribution.outcome_idx_context);
      std::vector<double> query_times;
      std::vector<const uuber::TrialParamsSoA *> query_params_soa;
      std::vector<double> query_weights;
      std::vector<std::size_t> query_to_point_index;
      query_times.reserve(points.size());
      query_params_soa.reserve(points.size());
      query_weights.reserve(points.size());
      query_to_point_index.reserve(points.size());

      const TrialParamSet *rep_params = nullptr;
      if (!for_each_prepared_query_lane(
              points, shared_trigger_supported, &rep_params,
              [&](std::size_t point_idx, const LogLikPreparedBatchPoint &point,
                  PreparedTrialParamsRuntime &, const uuber::TrialParamsSoA *soa,
                  double expansion_weight, bool) -> bool {
                const DirectTrialBatchSpec &point_spec =
                    direct_trial_batch_specs[flatten_point_index(
                        point.param_batch_idx, point.trial_idx)];
                query_times.push_back(point_spec.time);
                query_params_soa.push_back(soa);
                query_weights.push_back(contribution.scaled_weight *
                                        expansion_weight);
                query_to_point_index.push_back(point_idx);
                return true;
              })) {
        Rcpp::stop("Direct multi-parameter contribution preparation failed");
      }
      if (query_times.empty()) {
        continue;
      }

      std::vector<double> density_out;
      const bool batched = node_density_entry_batch_idx(
          ctx, contribution.node_id, query_times, contribution.component_idx,
          nullptr, false, nullptr, false, contribution.competitor_ids, rep_params,
          contribution.trial_type_key, false, contribution.outcome_idx_context,
          nullptr, false, density_out, nullptr, &query_params_soa, nullptr,
          &density_workspace);
      if (!batched || density_out.size() != query_times.size()) {
        Rcpp::stop("Direct multi-parameter contribution evaluation failed");
      }
      for (std::size_t query_idx = 0; query_idx < density_out.size(); ++query_idx) {
        const double density = density_out[query_idx];
        const double weight = query_weights[query_idx];
        if (!(std::isfinite(density) && density > 0.0 && std::isfinite(weight) &&
              weight > 0.0)) {
          continue;
        }
        const std::size_t point_idx = query_to_point_index[query_idx];
        point_probabilities[point_idx] = safe_density(
            point_probabilities[point_idx] + weight * density);
      }
    }
    assign_group_loglik(points, point_probabilities, trial_loglik,
                        trial_loglik_batched);
  }

  if (workload.enable_ranked_sequence_batch) {
    RankedBatchPlanner ranked_batch_transition_planner(ctx);
    for (auto &kv : ranked_sequence_batch_groups) {
      const std::vector<LogLikPreparedBatchPoint> &points = kv.second;
      if (points.empty()) {
        continue;
      }
      const std::size_t first_flat_idx = first_flat_idx_of_group(points);
      const RankedSequenceBatchSpec &group_spec =
          ranked_sequence_batch_specs[first_flat_idx];
      const RankedSequenceBatchTemplate *group_template = group_spec.template_ptr;
      if (group_template == nullptr || group_template->contributions.empty()) {
        continue;
      }

      const bool shared_trigger_ranked_supported =
          ranked_template_supports_shared_trigger(group_template);

      std::vector<const double *> query_times_by_point;
      std::vector<const uuber::TrialParamsSoA *> query_params_soa;
      std::vector<double> query_weights;
      std::vector<std::size_t> query_to_point_index;
      query_times_by_point.reserve(points.size());
      query_params_soa.reserve(points.size());
      query_weights.reserve(points.size());
      query_to_point_index.reserve(points.size());

      const TrialParamSet *rep_params = nullptr;
      if (!for_each_prepared_query_lane(
              points, shared_trigger_ranked_supported, &rep_params,
              [&](std::size_t point_idx, const LogLikPreparedBatchPoint &point,
                  PreparedTrialParamsRuntime &, const uuber::TrialParamsSoA *soa,
                  double expansion_weight, bool) -> bool {
                const RankedSequenceBatchSpec &point_spec =
                    ranked_sequence_batch_specs[flatten_point_index(
                        point.param_batch_idx, point.trial_idx)];
                query_times_by_point.push_back(point_spec.time_data);
                query_params_soa.push_back(soa);
                query_weights.push_back(expansion_weight);
                query_to_point_index.push_back(point_idx);
                return true;
              })) {
        Rcpp::stop("Ranked multi-parameter group preparation failed");
      }
      std::vector<const double *> compressed_query_times_by_point;
      std::vector<const uuber::TrialParamsSoA *> compressed_query_params_soa;
      std::vector<std::size_t> query_inverse_indices;
      compress_ranked_sequence_query_batch(
          query_times_by_point, query_params_soa, group_spec.time_count,
          compressed_query_times_by_point, compressed_query_params_soa,
          query_inverse_indices);

      std::vector<double> point_probabilities(points.size(), 0.0);
      for (const RankedSequenceContributionTemplate &contribution :
           group_template->contributions) {
        if (compressed_query_times_by_point.empty()) {
          continue;
        }
        std::vector<double> density_out;
        const bool batched = sequence_prefix_density_batch_resolved(
            ctx, contribution.sequence_outcome_indices,
            contribution.sequence_node_indices, compressed_query_times_by_point,
            contribution.component_idx, rep_params, contribution.trial_type_key,
            density_out, &compressed_query_params_soa,
            &ranked_batch_transition_planner,
            &contribution.step_competitor_ids_ptrs,
            &contribution.step_persistent_sources);
        if (!batched ||
            density_out.size() != compressed_query_times_by_point.size()) {
          Rcpp::stop("Ranked multi-parameter group evaluation failed");
        }
        std::vector<double> contribution_probabilities(points.size(), 0.0);
        for (std::size_t query_idx = 0; query_idx < query_inverse_indices.size();
             ++query_idx) {
          const std::size_t unique_idx = query_inverse_indices[query_idx];
          if (unique_idx >= density_out.size()) {
            continue;
          }
          const double density = density_out[unique_idx];
          const double mask_weight = query_weights[query_idx];
          if (!(std::isfinite(density) && density > 0.0 &&
                std::isfinite(mask_weight) && mask_weight > 0.0)) {
            continue;
          }
          contribution_probabilities[query_to_point_index[query_idx]] +=
              mask_weight * density;
        }
        for (std::size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
          const double contribution_density =
              safe_density(contribution_probabilities[point_idx]);
          if (std::isfinite(contribution_density) && contribution_density > 0.0) {
            point_probabilities[point_idx] = safe_density(
                point_probabilities[point_idx] +
                contribution.scaled_weight * contribution_density);
          }
        }
      }

      assign_group_loglik(points, point_probabilities, trial_loglik,
                          trial_loglik_batched);
    }
  }

  for (auto &kv : outcome_mass_batch_groups) {
    const std::vector<LogLikPreparedBatchPoint> &points = kv.second;
    if (points.empty()) {
      continue;
    }
    const OutcomeMassBatchKey &group_key = kv.first;
    const std::vector<int> *group_component_indices = group_key.component_indices;
    if (group_component_indices == nullptr ||
        group_key.nonresponse_outcome_indices == nullptr) {
      continue;
    }

    std::unordered_map<uuber::NAMapCacheKey, std::size_t,
                       uuber::NAMapCacheKeyHash>
        nonshared_point_index_by_value_key;
    std::unordered_map<const uuber::TrialParamsSoA *, std::size_t>
        shared_point_index_by_soa;
    nonshared_point_index_by_value_key.reserve(points.size());
    shared_point_index_by_soa.reserve(points.size());

    std::vector<const TrialParamSet *> unique_trial_params;
    std::vector<const uuber::TrialParamsSoA *> unique_params_soa;
    std::vector<const std::vector<double> *> component_weights_batch;
    std::vector<std::size_t> trial_to_point_index;
    std::vector<double> expansion_weights;
    std::vector<std::size_t> expansion_to_original_point;
    unique_trial_params.reserve(points.size());
    unique_params_soa.reserve(points.size());
    component_weights_batch.reserve(points.size());
    trial_to_point_index.reserve(points.size());
    expansion_weights.reserve(points.size());
    expansion_to_original_point.reserve(points.size());

    if (!for_each_prepared_query_lane(
            points, true, nullptr,
            [&](std::size_t point_idx, const LogLikPreparedBatchPoint &point,
                PreparedTrialParamsRuntime &runtime,
                const uuber::TrialParamsSoA *soa, double expansion_weight,
                bool shared_trigger_expanded) -> bool {
              const int forced_component_idx =
                  (point.trial_idx <
                   static_cast<int>(forced_component_idx_by_trial.size()))
                      ? forced_component_idx_by_trial[static_cast<std::size_t>(
                            point.trial_idx)]
                      : -1;
              const std::vector<double> *component_weights_ptr =
                  &prepared_params_for_point(point)
                       .weights_by_trial[static_cast<std::size_t>(
                           point.trial_idx)];
              if (forced_component_idx >= 0) {
                const ForcedComponentBundle &forced_bundle =
                    forced_bundle_for(forced_component_idx);
                component_weights_ptr = &forced_bundle.component_weights;
                if (&forced_bundle.component_indices !=
                    group_component_indices) {
                  return false;
                }
              } else if (&default_component_indices != group_component_indices) {
                return false;
              }

              std::size_t point_index = 0u;
              if (shared_trigger_expanded) {
                auto insert_result =
                    shared_point_index_by_soa.emplace(soa, unique_params_soa.size());
                point_index = insert_result.first->second;
                if (insert_result.second) {
                  unique_trial_params.push_back(runtime.params);
                  unique_params_soa.push_back(soa);
                }
              } else {
                auto insert_result = nonshared_point_index_by_value_key.emplace(
                    runtime.value_key, unique_params_soa.size());
                point_index = insert_result.first->second;
                if (insert_result.second) {
                  unique_trial_params.push_back(runtime.params);
                  unique_params_soa.push_back(soa);
                }
              }
              component_weights_batch.push_back(component_weights_ptr);
              trial_to_point_index.push_back(point_index);
              expansion_weights.push_back(expansion_weight);
              expansion_to_original_point.push_back(point_idx);
              return true;
            })) {
      Rcpp::stop("Outcome-mass multi-parameter group preparation failed");
    }

    std::vector<double> point_probabilities(points.size(), 0.0);
    if (!component_weights_batch.empty()) {
      std::vector<double> expanded_probabilities;
      if (!mix_outcome_mass_batch_idx(
              ctx, *group_key.nonresponse_outcome_indices,
              *group_component_indices, component_weights_batch,
              unique_trial_params, rel_tol, abs_tol, max_depth, false,
              expanded_probabilities, &trial_to_point_index, &unique_params_soa,
              nullptr) ||
          expanded_probabilities.size() != component_weights_batch.size()) {
        Rcpp::stop("Outcome-mass multi-parameter group evaluation failed");
      }
      for (std::size_t expansion_idx = 0; expansion_idx < expanded_probabilities.size();
           ++expansion_idx) {
        const double probability = expanded_probabilities[expansion_idx];
        const double weight = expansion_weights[expansion_idx];
        if (!(std::isfinite(probability) && probability > 0.0 &&
              std::isfinite(weight) && weight > 0.0)) {
          continue;
        }
        point_probabilities[expansion_to_original_point[expansion_idx]] +=
            weight * probability;
      }
    }
    for (double &probability : point_probabilities) {
      probability = safe_density(probability);
      if (group_key.probability_transform ==
          TrialProbabilityTransform::Complement) {
        probability = clamp_probability(1.0 - probability);
      }
    }
    assign_group_loglik(points, point_probabilities, trial_loglik,
                        trial_loglik_batched);
  }

  for (int param_batch_idx = 0; param_batch_idx < n_param_batches;
       ++param_batch_idx) {
    for (int trial_idx = 0; trial_idx < n_trials; ++trial_idx) {
      const std::size_t flat_idx = flatten_point_index(param_batch_idx, trial_idx);
      if (trial_idx >= static_cast<int>(ok_by_trial.size()) ||
          ok_by_trial[static_cast<std::size_t>(trial_idx)] == 0u) {
        trial_loglik[flat_idx] = min_ll;
        continue;
      }
      if (!trial_eval_inputs[static_cast<std::size_t>(trial_idx)].valid) {
        trial_loglik[flat_idx] = min_ll;
        continue;
      }
      if (trial_loglik_batched[flat_idx] == 0u) {
        Rcpp::stop("cpp_loglik_multiple left an eligible trial unassigned after "
                   "grouped batch evaluation");
      }
    }
  }

  for (int param_batch_idx = 0; param_batch_idx < n_param_batches;
       ++param_batch_idx) {
    double sum = 0.0;
    for (int trial_idx : expand_trial_indices) {
      sum += trial_loglik[flatten_point_index(param_batch_idx, trial_idx)];
    }
    total_loglik[param_batch_idx] = sum;
  }

  return total_loglik;
}

// [[Rcpp::export]]
double cpp_loglik(SEXP ctxSEXP, Rcpp::NumericMatrix params_mat,
                  Rcpp::DataFrame data_df, Rcpp::LogicalVector ok,
                  Rcpp::IntegerVector expand, double min_ll, double rel_tol,
                  double abs_tol, int max_depth) {
  Rcpp::XPtr<uuber::NativeContext> ctx_xptr(ctxSEXP);
  validate_loglik_context_or_stop(*ctx_xptr);
  LogLikWorkload workload;
  build_loglik_workload(*ctx_xptr, data_df, ok, expand, workload);
  std::vector<Rcpp::NumericMatrix> params_batch;
  params_batch.reserve(1);
  params_batch.push_back(params_mat);
  std::vector<LogLikPreparedParams> prepared_batch;
  prepared_batch.reserve(1);
  build_loglik_prepared_params_batch(*ctx_xptr, workload, params_batch,
                                     prepared_batch);
  const Rcpp::NumericVector out = execute_loglik_workload_multi(
      *ctx_xptr, workload, prepared_batch, min_ll, rel_tol, abs_tol,
      max_depth);
  return out.size() == 0 ? NA_REAL : out[0];
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_loglik_multiple(SEXP ctxSEXP, Rcpp::List params_list,
                                        Rcpp::DataFrame data_df,
                                        Rcpp::LogicalVector ok,
                                        Rcpp::IntegerVector expand,
                                        double min_ll, double rel_tol,
                                        double abs_tol, int max_depth) {
  Rcpp::XPtr<uuber::NativeContext> ctx_xptr(ctxSEXP);
  validate_loglik_context_or_stop(*ctx_xptr);
  LogLikWorkload workload;
  build_loglik_workload(*ctx_xptr, data_df, ok, expand, workload);
  Rcpp::NumericVector out(params_list.size());
  std::vector<int> active_param_indices;
  active_param_indices.reserve(params_list.size());
  std::vector<Rcpp::NumericMatrix> active_params_mats;
  active_params_mats.reserve(params_list.size());
  std::vector<LogLikPreparedParams> prepared_batch;
  prepared_batch.reserve(params_list.size());
  for (R_xlen_t i = 0; i < params_list.size(); ++i) {
    if (params_list[i] == R_NilValue) {
      out[i] = NA_REAL;
      continue;
    }
    active_param_indices.push_back(static_cast<int>(i));
    active_params_mats.emplace_back(params_list[i]);
  }
  if (!active_params_mats.empty()) {
    build_loglik_prepared_params_batch(*ctx_xptr, workload, active_params_mats,
                                       prepared_batch);
    const Rcpp::NumericVector active_out = execute_loglik_workload_multi(
        *ctx_xptr, workload, prepared_batch, min_ll, rel_tol, abs_tol,
        max_depth);
    for (std::size_t i = 0; i < active_param_indices.size(); ++i) {
      out[static_cast<R_xlen_t>(active_param_indices[i])] = active_out[i];
    }
  }
  return out;
}

// [[Rcpp::export]]
double native_outcome_probability_params_cpp_idx(
    SEXP ctxSEXP, int node_id, double upper, int component_idx,
    SEXP forced_complete, SEXP forced_survive,
    Rcpp::IntegerVector competitor_ids, double rel_tol, double abs_tol,
    Rcpp::Nullable<Rcpp::DataFrame> trial_rows) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::unique_ptr<TrialParamSet> params_holder =
      build_trial_params_from_df(*ctx, trial_rows);
  return native_outcome_probability_impl_idx(
      ctxSEXP, node_id, upper, component_idx, forced_complete, forced_survive,
      competitor_ids, rel_tol, abs_tol,
      params_holder ? params_holder.get() : nullptr, std::string(), false, -1);
}

// [[Rcpp::export]]
Rcpp::NumericVector native_outcome_probability_params_batch_cpp_idx(
    SEXP ctxSEXP, Rcpp::IntegerVector node_ids, double upper, int component_idx,
    Rcpp::List competitor_ids_list, double rel_tol, double abs_tol,
    Rcpp::Nullable<Rcpp::DataFrame> trial_rows) {
  if (competitor_ids_list.size() != node_ids.size()) {
    Rcpp::stop("competitor_ids_list must match node_ids length");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::unique_ptr<TrialParamSet> params_holder =
      build_trial_params_from_df(*ctx, trial_rows);
  Rcpp::NumericVector out(node_ids.size());
  std::vector<const TrialParamSet *> trial_params_batch(
      1u, params_holder ? params_holder.get() : nullptr);
  std::vector<std::uint8_t> consumed(static_cast<std::size_t>(node_ids.size()), 0u);
  DensityBatchWorkspace workspace;
  for (R_xlen_t i = 0; i < node_ids.size(); ++i) {
    if (consumed[static_cast<std::size_t>(i)] != 0u) {
      continue;
    }
    Rcpp::IntegerVector competitor_ids_vec =
        competitor_ids_list[i] == R_NilValue
            ? Rcpp::IntegerVector()
            : Rcpp::as<Rcpp::IntegerVector>(competitor_ids_list[i]);
    const std::vector<int> competitor_ids =
        integer_vector_to_std(competitor_ids_vec, false);
    std::vector<R_xlen_t> group_indices;
    std::vector<int> group_node_ids;
    for (R_xlen_t j = i; j < node_ids.size(); ++j) {
      if (consumed[static_cast<std::size_t>(j)] != 0u) {
        continue;
      }
      Rcpp::IntegerVector other_competitor_ids_vec =
          competitor_ids_list[j] == R_NilValue
              ? Rcpp::IntegerVector()
              : Rcpp::as<Rcpp::IntegerVector>(competitor_ids_list[j]);
      const std::vector<int> other_competitor_ids =
          integer_vector_to_std(other_competitor_ids_vec, false);
      if (other_competitor_ids != competitor_ids) {
        continue;
      }
      consumed[static_cast<std::size_t>(j)] = 1u;
      group_indices.push_back(j);
      group_node_ids.push_back(node_ids[j]);
    }

    std::vector<double> group_probabilities;
    if (group_node_ids.size() > 1u &&
        integrate_outcome_probabilities_from_density_multi_idx(
            *ctx, group_node_ids, upper, component_idx, competitor_ids, rel_tol,
            abs_tol, trial_params_batch, std::string(), false, -1,
            group_probabilities, nullptr, nullptr, &workspace) &&
        group_probabilities.size() == group_node_ids.size()) {
      for (std::size_t group_pos = 0; group_pos < group_indices.size(); ++group_pos) {
        out[group_indices[group_pos]] = group_probabilities[group_pos];
      }
      continue;
    }

    for (std::size_t group_pos = 0; group_pos < group_indices.size(); ++group_pos) {
      out[group_indices[group_pos]] = native_outcome_probability_bits_impl_idx(
          *ctx, group_node_ids[group_pos], upper, component_idx, nullptr, false,
          nullptr, false, competitor_ids, rel_tol, abs_tol,
          params_holder ? params_holder.get() : nullptr, std::string(), false,
          -1);
    }
  }
  return out;
}
