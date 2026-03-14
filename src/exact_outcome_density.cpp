#include "exact_outcome_density.h"

#include <Rcpp.h>

#include <algorithm>
#include <cstdlib>
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

constexpr double kDefaultRelTol = 1e-5;
constexpr double kDefaultAbsTol = 1e-6;
constexpr int kDefaultMaxDepth = 12;

inline bool exact_batch_debug_enabled() {
  static const bool enabled = []() {
    const char *raw = std::getenv("ACCUMULATR_DEBUG_EXACT_BATCH");
    return raw != nullptr && raw[0] != '\0' && raw[0] != '0';
  }();
  return enabled;
}

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
    const uuber::NativeContext &ctx, const ExactScenarioPoint &point);

bool exact_eval_relative_onset_base_point(
    const uuber::NativeContext &ctx, const ExactScenarioPoint &point,
    const uuber::LabelRef &source_ref, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const uuber::TrialParamsSoA *trial_params_soa,
    const uuber::AccDistParams &cfg,
    const LowerBoundTransform &lower_bound, double success_prob,
    double eval_t, double x_shift, bool need_density, bool need_cdf,
    double &base_density_out, double &base_cdf_out);

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

bool exact_eval_simple_acc_event_batch(
    const uuber::NativeContext &ctx, const uuber::IrEvent &event,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need,
    uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask = nullptr) {
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
  const uuber::TrialParamsSoA *trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);

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
  std::vector<std::uint8_t> direct_shifted_eval(point_count, 0u);
  bool need_pdf_vec = false;
  bool need_cdf_vec = false;

  for (std::size_t i = 0; i < point_count; ++i) {
    if (active_mask &&
        (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    const ExactScenarioPoint &point = points[i];
    const bool forced_complete = exact_forced_contains_id(
        ctx, point.forced_complete_bits, point.forced_complete_bits_valid,
        label_id);
    const bool forced_survive = exact_forced_contains_id(
        ctx, point.forced_survive_bits, point.forced_survive_bits_valid,
        label_id);
    const SourceTimeConstraint *constraint =
        time_constraints_find(&point.time_constraints, label_id);
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
    if (!needs_relative_onset) {
      direct_shifted_eval[i] = 1u;
      shifted_times[i] = point.t - onset_eff;
    } else {
      const SourceTimeConstraint *source_constraint =
          time_constraints_find(&point.time_constraints, onset_source_label_id);
      if (source_constraint != nullptr && source_constraint->has_exact &&
          std::isfinite(source_constraint->exact_time) &&
          !exact_forced_contains_id(ctx, point.forced_survive_bits,
                                    point.forced_survive_bits_valid,
                                    onset_source_label_id)) {
        direct_shifted_eval[i] = 1u;
        shifted_times[i] =
            point.t - source_constraint->exact_time - onset_lag - onset_eff;
      }
    }
    if (direct_shifted_eval[i] != 0u && needs_density(need) &&
        (kinds[i] == RelativeSelfConstraintKind::kNone ||
         kinds[i] == RelativeSelfConstraintKind::kBounded)) {
      need_pdf_vec = true;
    }
    if (direct_shifted_eval[i] != 0u &&
        ((needs_cdf(need) || needs_survival(need)) ||
         kinds[i] == RelativeSelfConstraintKind::kBounded) &&
        (kinds[i] == RelativeSelfConstraintKind::kNone ||
         kinds[i] == RelativeSelfConstraintKind::kBounded)) {
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

  for (std::size_t i = 0; i < point_count; ++i) {
    if (active_mask &&
        (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    const bool need_base_density =
        needs_density(need) &&
        (kinds[i] == RelativeSelfConstraintKind::kNone ||
         kinds[i] == RelativeSelfConstraintKind::kBounded);
    const bool need_base_cdf =
        (needs_cdf(need) || needs_survival(need) ||
         kinds[i] == RelativeSelfConstraintKind::kBounded) &&
        (kinds[i] == RelativeSelfConstraintKind::kNone ||
         kinds[i] == RelativeSelfConstraintKind::kBounded);
    double base_density = 0.0;
    double base_cdf = 0.0;
    if (direct_shifted_eval[i] != 0u) {
      if (need_base_density) {
        base_density = safe_density(success_prob * pdf_values[i]);
      }
      if (need_base_cdf) {
        base_cdf = clamp_probability(success_prob * cdf_values[i]);
      }
    } else if (needs_relative_onset &&
               !exact_eval_relative_onset_base_point(
                   ctx, points[i], onset_source_ref, component_idx,
                   trial_params, trial_type_key, trial_params_soa, cfg,
                   lower_bound, success_prob, points[i].t,
                   onset_lag + onset_eff,
                   need_base_density, need_base_cdf, base_density, base_cdf)) {
      return false;
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
          (!(points[i].t > self_bound_lower[i]) ||
           points[i].t > self_bound_upper[i])) {
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
                     ? clamp_probability(success_prob * eval_cdf_single_with_lower_bound(
                                               cfg, lower - (points[i].t - shifted_times[i]),
                                               lower_bound))
                     : [&]() {
                         double cdf_at_lower = 0.0;
                         double density_unused = 0.0;
                         const bool ok = exact_eval_relative_onset_base_point(
                             ctx, points[i], onset_source_ref, component_idx,
                             trial_params, trial_type_key, trial_params_soa, cfg,
                             lower_bound, success_prob, lower,
                             onset_lag + onset_eff, false, true,
                             density_unused, cdf_at_lower);
                         return ok ? clamp_probability(cdf_at_lower) : 0.0;
                       }())
              : 0.0;
      const double upper_cdf =
          std::isfinite(upper)
              ? (direct_shifted_eval[i] != 0u
                     ? clamp_probability(success_prob * eval_cdf_single_with_lower_bound(
                                               cfg, upper - (points[i].t - shifted_times[i]),
                                               lower_bound))
                     : [&]() {
                         double cdf_at_upper = 0.0;
                         double density_unused = 0.0;
                         const bool ok = exact_eval_relative_onset_base_point(
                             ctx, points[i], onset_source_ref, component_idx,
                             trial_params, trial_type_key, trial_params_soa, cfg,
                             lower_bound, success_prob, upper,
                             onset_lag + onset_eff, false, true,
                             density_unused, cdf_at_upper);
                         return ok ? clamp_probability(cdf_at_upper)
                                   : clamp_probability(success_prob);
                       }())
              : clamp_probability(success_prob);
      const double condition_mass = clamp_probability(upper_cdf - lower_cdf);
      if (!(upper > lower) || !(condition_mass > 0.0)) {
        break;
      }
      if (needs_density(need)) {
        if (std::isfinite(points[i].t) && points[i].t > lower &&
            points[i].t <= upper) {
          out_values.density[i] =
              safe_density(base_density / condition_mass);
        } else {
          out_values.density[i] = 0.0;
        }
      }
      if (needs_cdf(need) || needs_survival(need)) {
        double cdf_cond = 0.0;
        if (!std::isfinite(points[i].t) || points[i].t >= upper) {
          cdf_cond = 1.0;
        } else if (!(points[i].t > lower)) {
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

inline uuber::KernelNodeValues kernel_batch_values_at(
    const uuber::KernelNodeBatchValues &values, std::size_t idx) {
  uuber::KernelNodeValues out;
  if (idx < values.density.size()) {
    out.density = values.density[idx];
  }
  if (idx < values.survival.size()) {
    out.survival = values.survival[idx];
  }
  if (idx < values.cdf.size()) {
    out.cdf = values.cdf[idx];
  }
  return out;
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

inline bool time_constraint_map_equal(const TimeConstraintMap &lhs,
                                      const TimeConstraintMap &rhs) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  auto lhs_it = lhs.begin();
  auto rhs_it = rhs.begin();
  while (lhs_it != lhs.end() && rhs_it != rhs.end()) {
    if (lhs_it->first != rhs_it->first ||
        !time_constraint_equal(lhs_it->second, rhs_it->second)) {
      return false;
    }
    ++lhs_it;
    ++rhs_it;
  }
  return lhs_it == lhs.end() && rhs_it == rhs.end();
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

bool exact_points_share_relevant_state(
    const std::vector<ExactScenarioPoint> &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, const uuber::NativeContext &ctx,
    std::size_t &anchor_idx);

inline bool exact_prepare_shared_state_storage(
    const uuber::NativeContext &ctx,
    const std::vector<ExactScenarioPoint> &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, std::size_t &anchor_idx,
    ExactSharedStateStorage &storage) {
  anchor_idx = points.size();
  if (!exact_points_share_relevant_state(points, active_mask,
                                         relevant_source_ids, ctx,
                                         anchor_idx)) {
    return false;
  }
  const ExactScenarioPoint &anchor = points[anchor_idx];
  build_filtered_time_constraints(anchor.time_constraints, relevant_source_ids,
                                  storage.time_constraints);
  build_filtered_forced_bits(
      ctx, anchor.forced_complete_bits, anchor.forced_complete_bits_valid,
      relevant_source_ids, storage.forced_complete_bits,
      storage.forced_complete_bits_valid);
  build_filtered_forced_bits(
      ctx, anchor.forced_survive_bits, anchor.forced_survive_bits_valid,
      relevant_source_ids, storage.forced_survive_bits,
      storage.forced_survive_bits_valid);
  return true;
}

template <typename SharedFn, typename FallbackFn>
bool exact_eval_with_shared_relevant_state(
    const uuber::NativeContext &ctx,
    const std::vector<ExactScenarioPoint> &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    SharedFn &&shared_eval, FallbackFn &&fallback_eval) {
  std::size_t anchor_idx = points.size();
  ExactSharedStateStorage storage;
  if (!exact_prepare_shared_state_storage(ctx, points, active_mask,
                                          relevant_source_ids, anchor_idx,
                                          storage)) {
    return fallback_eval();
  }

  const ExactScenarioPoint &anchor = points[anchor_idx];
  NodeEvalState shared_state(
      ctx, anchor.t, component_idx, trial_params, trial_type_key, false, -1,
      &storage.time_constraints, nullptr,
      storage.forced_complete_bits_valid ? &storage.forced_complete_bits
                                         : nullptr,
      storage.forced_complete_bits_valid,
      storage.forced_survive_bits_valid ? &storage.forced_survive_bits
                                        : nullptr,
      storage.forced_survive_bits_valid, nullptr);
  return shared_eval(shared_state);
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

inline uuber::LabelRef exact_event_label_ref(const uuber::IrEvent &event) {
  uuber::LabelRef ref;
  ref.label_id = event.label_id;
  ref.acc_idx = event.acc_idx;
  ref.pool_idx = event.pool_idx;
  ref.outcome_idx = event.outcome_idx;
  return ref;
}

inline std::uint32_t exact_event_node_flags(const uuber::NativeContext &ctx,
                                            const uuber::IrEvent &event) {
  if (event.node_idx >= 0 &&
      event.node_idx < static_cast<int>(ctx.ir.nodes.size())) {
    return ctx.ir.nodes[static_cast<std::size_t>(event.node_idx)].flags;
  }
  return 0u;
}

inline ForcedStateView exact_make_point_forced_state(
    const uuber::NativeContext &ctx, const ExactScenarioPoint &point) {
  return make_forced_state_view(
      nullptr,
      point.forced_complete_bits_valid ? &point.forced_complete_bits : nullptr,
      &point.forced_complete_bits_valid,
      point.forced_survive_bits_valid ? &point.forced_survive_bits : nullptr,
      &point.forced_survive_bits_valid,
      ctx.ir.label_id_to_bit_idx.empty() ? nullptr : &ctx.ir.label_id_to_bit_idx);
}

inline NodeEvalState exact_make_local_state(
    const uuber::NativeContext &ctx, double t, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const ExactScenarioPoint &point) {
  return NodeEvalState(
      ctx, t, component_idx, trial_params, trial_type_key, false, -1,
      &point.time_constraints, nullptr,
      point.forced_complete_bits_valid ? &point.forced_complete_bits : nullptr,
      point.forced_complete_bits_valid,
      point.forced_survive_bits_valid ? &point.forced_survive_bits : nullptr,
      point.forced_survive_bits_valid, nullptr);
}

bool exact_eval_relative_onset_base_point(
    const uuber::NativeContext &ctx, const ExactScenarioPoint &point,
    const uuber::LabelRef &source_ref, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const uuber::TrialParamsSoA *trial_params_soa,
    const uuber::AccDistParams &cfg,
    const LowerBoundTransform &lower_bound, double success_prob,
    double eval_t, double x_shift, bool need_density, bool need_cdf,
    double &base_density_out, double &base_cdf_out) {
  base_density_out = 0.0;
  base_cdf_out = 0.0;
  if (!need_density && !need_cdf) {
    return true;
  }
  if (source_ref.acc_idx < 0 && source_ref.pool_idx < 0) {
    return true;
  }

  const int source_label_id = source_ref.label_id;
  const ForcedStateView forced_state = exact_make_point_forced_state(ctx, point);
  if (source_label_id >= 0 && source_label_id != NA_INTEGER &&
      evaluator_state_contains_survive_at(
          forced_state, &point.time_constraints, source_label_id, eval_t)) {
    return true;
  }

  double bound_lower = 0.0;
  double bound_upper = std::numeric_limits<double>::infinity();
  bool has_conditioning_bounds = false;
  bool has_exact = false;
  double exact_time = std::numeric_limits<double>::quiet_NaN();

  if (source_label_id >= 0 && source_label_id != NA_INTEGER) {
    if (forced_state_contains_complete(forced_state, source_label_id)) {
      bound_upper = std::min(bound_upper, eval_t);
      has_conditioning_bounds = true;
    }
    bool source_has_bounds = false;
    double source_bound_lower = 0.0;
    double source_bound_upper = std::numeric_limits<double>::infinity();
    (void)evaluator_resolve_label_time_constraint(
        &point.time_constraints, source_label_id, has_exact, exact_time,
        source_has_bounds, source_bound_lower, source_bound_upper);
    if (source_has_bounds) {
      bound_lower = std::max(bound_lower, source_bound_lower);
      bound_upper = std::min(bound_upper, source_bound_upper);
      has_conditioning_bounds = true;
    }
  }

  if (has_exact) {
    const double shifted_time = point.t - exact_time - x_shift;
    if (need_density) {
      base_density_out = safe_density(
          success_prob *
          eval_pdf_single_with_lower_bound(cfg, shifted_time, lower_bound));
    }
    if (need_cdf) {
      base_cdf_out = clamp_probability(
          success_prob *
          eval_cdf_single_with_lower_bound(cfg, shifted_time, lower_bound));
    }
    return true;
  }

  if (!(bound_upper > bound_lower)) {
    return true;
  }
  if (std::isfinite(eval_t) && (eval_t - x_shift <= bound_lower)) {
    return true;
  }

  auto source_eval = [&](double u, EvalNeed source_need) -> NodeEvalResult {
    return evaluator_eval_event_ref_idx(
        ctx, source_ref, 0u, u, component_idx, source_need, trial_params,
        trial_type_key, false, -1, &point.time_constraints, nullptr, nullptr,
        nullptr, nullptr, trial_params_soa, &forced_state);
  };
  auto source_cdf_at = [&](double u) -> double {
    if (u <= 0.0) {
      return 0.0;
    }
    const NodeEvalResult source = source_eval(u, EvalNeed::kCDF);
    return clamp_probability(source.cdf);
  };

  double source_mass = 1.0;
  if (has_conditioning_bounds) {
    const double lower_cdf =
        std::isfinite(bound_lower) ? source_cdf_at(bound_lower) : 0.0;
    const double upper_cdf =
        std::isfinite(bound_upper)
            ? source_cdf_at(bound_upper)
            : source_cdf_at(std::numeric_limits<double>::infinity());
    source_mass = upper_cdf - lower_cdf;
    if (!std::isfinite(source_mass) || source_mass <= 0.0) {
      return true;
    }
  }

  const double upper_int =
      std::min(bound_upper, std::max(bound_lower, eval_t - x_shift));
  if (!(upper_int > bound_lower)) {
    return true;
  }

  if (!std::isfinite(eval_t)) {
    if (need_cdf) {
      const double cdf_total =
          has_conditioning_bounds
              ? clamp_probability(success_prob)
              : clamp_probability(
                    success_prob *
                    source_cdf_at(std::numeric_limits<double>::infinity()));
      base_cdf_out = cdf_total;
    }
    return true;
  }

  auto source_pdf = [&](double u) -> double {
    const NodeEvalResult source = source_eval(u, EvalNeed::kDensity);
    return safe_density(source.density);
  };
  auto cond_source_pdf = [&](double u) -> double {
    if (!std::isfinite(u) || u <= bound_lower || u > bound_upper) {
      return 0.0;
    }
    const double dens = source_pdf(u);
    if (dens <= 0.0) {
      return 0.0;
    }
    if (!has_conditioning_bounds) {
      return dens;
    }
    return dens / source_mass;
  };

  const FusedIntegralResult fused = integrate_fused_onset_terms(
      cfg, point.t, x_shift, bound_lower, upper_int, need_density, need_cdf,
      &lower_bound, cond_source_pdf);

  if (need_density) {
    if (fused.ok) {
      base_density_out = safe_density(success_prob * fused.density);
    } else {
      auto integrand_density = [&](double u) -> double {
        const double fs = cond_source_pdf(u);
        if (fs <= 0.0) {
          return 0.0;
        }
        const double shifted_time = eval_t - u - x_shift;
        const double fx =
            eval_pdf_single_with_lower_bound(cfg, shifted_time, lower_bound);
        if (!std::isfinite(fx) || fx <= 0.0) {
          return 0.0;
        }
        return fs * fx;
      };
      base_density_out = safe_density(
          success_prob *
          uuber::integrate_boost_fn(integrand_density, bound_lower, upper_int,
                                    kDefaultRelTol, kDefaultAbsTol,
                                    kDefaultMaxDepth));
    }
  }

  if (need_cdf) {
    double cdf_cond = 0.0;
    if (fused.ok) {
      cdf_cond = fused.cdf;
    } else {
      auto integrand_cdf = [&](double u) -> double {
        const double fs = cond_source_pdf(u);
        if (fs <= 0.0) {
          return 0.0;
        }
        const double shifted_time = eval_t - u - x_shift;
        const double Fx =
            eval_cdf_single_with_lower_bound(cfg, shifted_time, lower_bound);
        if (!std::isfinite(Fx) || Fx <= 0.0) {
          return 0.0;
        }
        return fs * clamp_probability(Fx);
      };
      cdf_cond = uuber::integrate_boost_fn(integrand_cdf, bound_lower,
                                           upper_int, kDefaultRelTol,
                                           kDefaultAbsTol, kDefaultMaxDepth);
    }
    base_cdf_out = clamp_probability(success_prob * clamp_probability(cdf_cond));
  }

  return true;
}

template <typename PointVec>
bool exact_eval_event_pointwise_batch(
    const uuber::NativeContext &ctx, const PointVec &points,
    const std::vector<double> &batch_times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const uuber::IrEvent &event, EvalNeed event_need,
    const uuber::TrialParamsSoA *trial_params_soa,
    uuber::KernelNodeBatchValues &event_out,
    const std::vector<std::uint8_t> *active_mask) {
  const std::size_t point_count = points.size();
  if (batch_times.size() != point_count) {
    return false;
  }
  const uuber::LabelRef ref = exact_event_label_ref(event);
  const std::uint32_t node_flags = exact_event_node_flags(ctx, event);
  for (std::size_t i = 0; i < point_count; ++i) {
    if (!exact_point_active(active_mask, i)) {
      continue;
    }
    const ExactScenarioPoint &point = points[i];
    const ForcedStateView forced_state = exact_make_point_forced_state(ctx, point);
    const NodeEvalResult event_eval = evaluator_eval_event_ref_idx(
        ctx, ref, node_flags, batch_times[i], component_idx, event_need,
        trial_params, trial_type_key, false, -1, &point.time_constraints, nullptr,
        point.forced_complete_bits_valid ? &point.forced_complete_bits : nullptr,
        point.forced_survive_bits_valid ? &point.forced_survive_bits : nullptr,
        ctx.ir.label_id_to_bit_idx.empty() ? nullptr : &ctx.ir.label_id_to_bit_idx,
        trial_params_soa, &forced_state);
    event_out.density[i] = safe_density(event_eval.density);
    event_out.survival[i] = clamp_probability(event_eval.survival);
    event_out.cdf[i] = clamp_probability(event_eval.cdf);
  }
  return true;
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

template <typename PointVec>
bool exact_eval_guard_pointwise_batch(
    const uuber::NativeContext &ctx, const PointVec &points,
    const std::vector<double> &batch_times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const uuber::KernelOp &op,
    const uuber::KernelNodeBatchValues &reference_values,
    const uuber::KernelNodeBatchValues &blocker_values,
    const uuber::KernelEvalNeed &kneed,
    uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask) {
  const std::size_t point_count = points.size();
  if (!exact_guard_batch_inputs_valid(batch_times, point_count, reference_values,
                                      blocker_values)) {
    return false;
  }
  for (std::size_t i = 0; i < point_count; ++i) {
    if (!exact_point_active(active_mask, i)) {
      continue;
    }
    const ExactScenarioPoint &point = points[i];
    NodeEvalState local_state = exact_make_local_state(
        ctx, batch_times[i], component_idx, trial_params, trial_type_key, point);
    uuber::KernelGuardEvalFn guard_eval =
        evaluator_make_kernel_guard_eval(local_state);
    const uuber::KernelNodeValues values = guard_eval(
        op, kernel_batch_values_at(reference_values, i),
        kernel_batch_values_at(blocker_values, i), kneed);
    if (kneed.density) {
      out_values.density[i] = safe_density(values.density);
    }
    out_values.survival[i] = clamp_probability(values.survival);
    out_values.cdf[i] = clamp_probability(values.cdf);
  }
  return true;
}

template <typename PointVec, typename DirectEvalFn>
bool exact_eval_event_batch_common(
    const uuber::NativeContext &ctx, const PointVec &points,
    const std::vector<double> &batch_times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const uuber::TrialParamsSoA *trial_params_soa, int event_idx,
    const uuber::KernelEvalNeed &kneed,
    uuber::KernelNodeBatchValues &event_out,
    const std::vector<std::uint8_t> *active_mask, DirectEvalFn &&direct_eval) {
  const std::size_t point_count = points.size();
  if (batch_times.size() != point_count ||
      event_idx < 0 ||
      event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return false;
  }
  exact_init_batch_values(point_count, event_out);

  const uuber::IrEvent &event =
      ctx.ir.events[static_cast<std::size_t>(event_idx)];
  const EvalNeed event_need = exact_eval_need_from_kernel_need(kneed);

  if (direct_eval(event, event_need, event_out)) {
    return true;
  }
  if (event.node_idx >= 0) {
    const std::vector<int> event_source_ids =
        exact_relevant_source_ids_for_node(ctx, event.node_idx);
    return exact_eval_with_shared_relevant_state(
        ctx, points, active_mask, event_source_ids, component_idx,
        trial_params, trial_type_key,
        [&](NodeEvalState &shared_state) -> bool {
          uuber::KernelEventBatchEvalFn shared_event_eval =
              evaluator_make_kernel_event_eval_batch(shared_state);
          return shared_event_eval(event_idx, batch_times, kneed, event_out);
        },
        [&]() -> bool {
          return exact_eval_event_pointwise_batch(
              ctx, points, batch_times, component_idx, trial_params,
              trial_type_key, event, event_need, trial_params_soa, event_out,
              active_mask);
        });
  }
  return exact_eval_event_pointwise_batch(
      ctx, points, batch_times, component_idx, trial_params, trial_type_key,
      event, event_need, trial_params_soa, event_out, active_mask);
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
    bool allow_pointwise_without_node_idx) {
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
  if (op.node_idx >= 0) {
    const std::vector<int> guard_source_ids =
        exact_relevant_source_ids_for_node(ctx, op.node_idx);
    return exact_eval_with_shared_relevant_state(
        ctx, points, active_mask, guard_source_ids, component_idx,
        trial_params, trial_type_key,
        [&](NodeEvalState &shared_state) -> bool {
          uuber::KernelGuardBatchEvalFn shared_guard_eval =
              evaluator_make_kernel_guard_eval_batch(shared_state);
          return shared_guard_eval(op, batch_times, reference_values,
                                   blocker_values, kneed, guard_out);
        },
        [&]() -> bool {
          return exact_eval_guard_pointwise_batch(
              ctx, points, batch_times, component_idx, trial_params,
              trial_type_key, op, reference_values, blocker_values, kneed,
              guard_out, active_mask);
        });
  }
  if (!allow_pointwise_without_node_idx) {
    return true;
  }
  return exact_eval_guard_pointwise_batch(
      ctx, points, batch_times, component_idx, trial_params, trial_type_key,
      op, reference_values, blocker_values, kneed, guard_out, active_mask);
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask);

bool exact_points_share_relevant_state(
    const std::vector<ExactScenarioPoint> &points,
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
  const ExactScenarioPoint &anchor = points[anchor_idx];
  TimeConstraintMap anchor_constraints;
  build_filtered_time_constraints(anchor.time_constraints, relevant_source_ids,
                                  anchor_constraints);
  for (std::size_t i = anchor_idx + 1u; i < points.size(); ++i) {
    if (active_mask && (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    const ExactScenarioPoint &point = points[i];
    for (int source_id : relevant_source_ids) {
      if (filtered_forced_contains_id(anchor.forced_complete_bits,
                                      anchor.forced_complete_bits_valid,
                                      source_id, ctx) !=
              filtered_forced_contains_id(point.forced_complete_bits,
                                          point.forced_complete_bits_valid,
                                          source_id, ctx) ||
          filtered_forced_contains_id(anchor.forced_survive_bits,
                                      anchor.forced_survive_bits_valid,
                                      source_id, ctx) !=
              filtered_forced_contains_id(point.forced_survive_bits,
                                          point.forced_survive_bits_valid,
                                          source_id, ctx)) {
        return false;
      }
    }
    TimeConstraintMap point_constraints;
    build_filtered_time_constraints(point.time_constraints, relevant_source_ids,
                                    point_constraints);
    if (!time_constraint_map_equal(anchor_constraints, point_constraints)) {
      return false;
    }
  }
  return true;
}

bool exact_eval_node_batch_with_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const std::vector<std::uint8_t> *active_mask) {
  out_values = uuber::KernelNodeBatchValues{};
  if (points.empty()) {
    return true;
  }
  ExactNodeBatchDepthGuard depth_guard;
  const std::size_t active_point_count =
      exact_count_active_points(active_mask, points.size());
  record_unified_outcome_exact_node_batch_call(
      static_cast<std::uint64_t>(active_point_count),
      depth_guard.recursive_call());

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
              need, out_values, active_mask)) {
        return true;
      }
    }
  }

  std::vector<double> times(points.size(), 0.0);
  for (std::size_t i = 0; i < points.size(); ++i) {
    times[i] = points[i].t;
  }

  const std::vector<int> relevant_source_ids =
      exact_relevant_source_ids_for_node(ctx, node_idx);
  auto pointwise_fallback = [&]() -> bool {
    record_unified_outcome_exact_node_batch_pointwise_fallback(
        static_cast<std::uint64_t>(active_point_count));
    if (exact_batch_debug_enabled()) {
      const uuber::IrNode &node = ir_node_required(ctx, node_idx);
      Rcpp::Rcout << "exact_batch_fallback node_idx=" << node_idx
                  << " op=" << static_cast<int>(node.op)
                  << " event_idx=" << node.event_idx
                  << " active_points=" << active_point_count << "\n";
    }
    const uuber::TrialParamsSoA *trial_params_soa =
        resolve_trial_params_soa(ctx, trial_params);
    uuber::KernelBatchRuntimeState batch_runtime;
    uuber::KernelEventBatchEvalFn event_eval_batch =
        [&](int event_idx, const std::vector<double> &batch_times,
            const uuber::KernelEvalNeed &kneed,
            uuber::KernelNodeBatchValues &event_out) -> bool {
      return exact_eval_event_batch_common(
          ctx, points, batch_times, component_idx, trial_params,
          trial_type_key, trial_params_soa, event_idx, kneed, event_out,
          active_mask,
          [&](const uuber::IrEvent &event, EvalNeed event_need,
              uuber::KernelNodeBatchValues &out) -> bool {
            return exact_eval_simple_acc_event_batch(
                ctx, event, points, component_idx, trial_params, trial_type_key,
                event_need, out, active_mask);
          });
    };

    uuber::KernelGuardBatchEvalFn guard_eval_batch =
        [&](const uuber::KernelOp &op, const std::vector<double> &batch_times,
            const uuber::KernelNodeBatchValues &reference_values,
            const uuber::KernelNodeBatchValues &blocker_values,
            const uuber::KernelEvalNeed &kneed,
            uuber::KernelNodeBatchValues &guard_out) -> bool {
      return exact_eval_guard_batch_common(
          ctx, points, batch_times, component_idx, trial_params,
          trial_type_key, op, reference_values, blocker_values, kneed,
          guard_out, active_mask, true);
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
  };

  return exact_eval_with_shared_relevant_state(
      ctx, points, active_mask, relevant_source_ids, component_idx,
      trial_params, trial_type_key,
      [&](NodeEvalState &shared_state) -> bool {
        uuber::KernelBatchRuntimeState batch_runtime;
        return evaluator_eval_node_batch_with_state_dense(
            node_idx, times, shared_state, need, batch_runtime, out_values);
      },
      pointwise_fallback);
}

bool exact_collect_scenarios_batch_from_points_impl(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &scenario_points) {
  scenario_points.clear();
  const RankedNodeTransitionPlan &plan = compiler.plan_for_node(node_idx);
  if (!plan.valid || plan.transitions.empty()) {
    return false;
  }

  scenario_points.reserve(seed_points.size() * plan.transitions.size());
  for (const RankedTransitionTemplate &transition : plan.transitions) {
    std::vector<ExactScenarioPoint> transition_points = seed_points;
    std::vector<std::uint8_t> active(transition_points.size(), 0u);
    bool any_active = false;
    for (std::size_t i = 0; i < transition_points.size(); ++i) {
      const ExactScenarioPoint &point = transition_points[i];
      const bool is_active = std::isfinite(point.t) && point.t >= 0.0 &&
                             std::isfinite(point.weight) && point.weight > 0.0;
      active[i] = is_active ? 1u : 0u;
      any_active = any_active || is_active;
    }
    if (!any_active) {
      continue;
    }

    bool transition_valid = true;
    for (const RankedTransitionStep &step : transition.steps) {
      if (step.kind == RankedTransitionStepKind::EvalDensityNode ||
          step.kind == RankedTransitionStepKind::EvalCDFNode) {
        uuber::KernelNodeBatchValues values;
        const EvalNeed step_need =
            (step.kind == RankedTransitionStepKind::EvalDensityNode)
                ? EvalNeed::kDensity
                : EvalNeed::kCDF;
        if (!exact_eval_node_batch_with_points(
                ctx, step.node_idx, transition_points, component_idx,
                trial_params, trial_type_key, step_need, values, &active)) {
          transition_valid = false;
          break;
        }
        for (std::size_t i = 0; i < transition_points.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          const double factor =
              (step.kind == RankedTransitionStepKind::EvalDensityNode)
                  ? safe_density(values.density[i])
                  : clamp_probability(values.cdf[i]);
          if (!std::isfinite(factor) || factor <= 0.0) {
            active[i] = 0u;
            transition_points[i].weight = 0.0;
            continue;
          }
          transition_points[i].weight *= factor;
          if (!std::isfinite(transition_points[i].weight) ||
              transition_points[i].weight <= 0.0) {
            active[i] = 0u;
            transition_points[i].weight = 0.0;
          }
        }
      } else if (step.kind == RankedTransitionStepKind::EvalGuardEffective) {
        const uuber::IrNode &guard_node = ir_node_required(ctx, step.node_idx);
        if (guard_node.blocker_idx < 0) {
          transition_valid = false;
          break;
        }
        uuber::KernelNodeBatchValues blocker_values;
        if (!exact_eval_node_batch_with_points(
                ctx, guard_node.blocker_idx, transition_points, component_idx,
                trial_params, trial_type_key, EvalNeed::kSurvival,
                blocker_values, &active)) {
          transition_valid = false;
          break;
        }
        for (std::size_t i = 0; i < transition_points.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          const double factor = clamp_probability(blocker_values.survival[i]);
          if (!std::isfinite(factor) || factor <= 0.0) {
            active[i] = 0u;
            transition_points[i].weight = 0.0;
            continue;
          }
          transition_points[i].weight *= factor;
          if (!std::isfinite(transition_points[i].weight) ||
              transition_points[i].weight <= 0.0) {
            active[i] = 0u;
            transition_points[i].weight = 0.0;
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddCompleteSources) {
        for (std::size_t i = 0; i < transition_points.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          ExactScenarioPoint &point = transition_points[i];
          bool ok = true;
          for (int id : step.source_ids) {
            if (!time_constraints_mark_complete(
                    id, point.t, step.bind_exact_current_time,
                    point.time_constraints)) {
              ok = false;
              break;
            }
          }
          if (!ok) {
            active[i] = 0u;
            point.weight = 0.0;
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddSurviveSources) {
        for (std::size_t i = 0; i < transition_points.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          ExactScenarioPoint &point = transition_points[i];
          bool ok = true;
          for (int id : step.source_ids) {
            if (!time_constraints_mark_survive(id, point.t,
                                               point.time_constraints)) {
              ok = false;
              break;
            }
          }
          if (!ok) {
            active[i] = 0u;
            point.weight = 0.0;
          }
        }
      } else if (step.kind == RankedTransitionStepKind::AddOrWitnessFromSources) {
        for (std::size_t i = 0; i < transition_points.size(); ++i) {
          if (active[i] == 0u) {
            continue;
          }
          ExactScenarioPoint &point = transition_points[i];
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
                (point.forced_complete_bits_valid && bit_idx >= 0 &&
                 point.forced_complete_bits.test(bit_idx)) ||
                time_constraints_contains_complete_at(&point.time_constraints,
                                                     id, point.t);
            if (!is_forced) {
              witness = id;
              all_forced = false;
              break;
            }
          }
          if (all_forced || witness == NA_INTEGER ||
              !time_constraints_mark_survive(
                  witness, point.t, point.time_constraints)) {
            active[i] = 0u;
            point.weight = 0.0;
          }
        }
      }

      any_active = false;
      for (std::uint8_t bit : active) {
        if (bit != 0u) {
          any_active = true;
          break;
        }
      }
      if (!any_active) {
        transition_valid = false;
        break;
      }
    }

    if (!transition_valid) {
      continue;
    }
    for (std::size_t i = 0; i < transition_points.size(); ++i) {
      ExactScenarioPoint &point = transition_points[i];
      if (active[i] == 0u || !std::isfinite(point.weight) ||
          point.weight <= 0.0) {
        continue;
      }
      scenario_points.push_back(std::move(point));
    }
  }

  return !scenario_points.empty();
}

bool exact_collect_scenarios_batch_impl(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    std::vector<ExactScenarioPoint> &scenario_points) {
  std::vector<ExactScenarioPoint> seed_points(times.size());
  for (std::size_t i = 0; i < times.size(); ++i) {
    ExactScenarioPoint &point = seed_points[i];
    point.t = times[i];
    point.density_index = i;
    point.weight =
        (std::isfinite(times[i]) && times[i] >= 0.0) ? 1.0 : 0.0;
    point.forced_complete_bits = state.forced_complete_bits;
    point.forced_survive_bits = state.forced_survive_bits;
    point.forced_complete_bits_valid = state.forced_complete_bits_valid;
    point.forced_survive_bits_valid = state.forced_survive_bits_valid;
    point.time_constraints = state.time_constraints;
  }
  return exact_collect_scenarios_batch_from_points_impl(
      compiler, ctx, node_idx, seed_points, state.component_idx,
      state.trial_params, state.trial_type_key, scenario_points);
}

std::vector<int> bitset_to_label_ids(const uuber::NativeContext &ctx,
                                     const uuber::BitsetState &bits,
                                     bool bits_valid) {
  if (!bits_valid) {
    return {};
  }
  std::vector<int> out;
  out.reserve(ctx.ir.label_id_to_bit_idx.size());
  for (const auto &kv : ctx.ir.label_id_to_bit_idx) {
    if (bits.test(kv.second)) {
      out.push_back(kv.first);
    }
  }
  sort_unique(out);
  return out;
}

Rcpp::NumericVector exact_times_to_numeric(const TimeConstraintMap &constraints) {
  std::size_t count = 0;
  for (const auto &kv : constraints) {
    if (kv.second.has_exact) {
      ++count;
    }
  }
  Rcpp::NumericVector out(static_cast<R_xlen_t>(count));
  Rcpp::CharacterVector names(static_cast<R_xlen_t>(count));
  R_xlen_t idx = 0;
  for (const auto &kv : constraints) {
    if (!kv.second.has_exact) {
      continue;
    }
    out[idx] = kv.second.exact_time;
    names[idx] = std::to_string(kv.first);
    ++idx;
  }
  out.attr("names") = names;
  return out;
}

Rcpp::List bounds_to_list(const TimeConstraintMap &constraints) {
  std::size_t count = 0;
  for (const auto &kv : constraints) {
    if (kv.second.has_lower || kv.second.has_upper) {
      ++count;
    }
  }
  Rcpp::NumericVector lower(static_cast<R_xlen_t>(count));
  Rcpp::NumericVector upper(static_cast<R_xlen_t>(count));
  Rcpp::CharacterVector names(static_cast<R_xlen_t>(count));
  R_xlen_t idx = 0;
  for (const auto &kv : constraints) {
    if (!(kv.second.has_lower || kv.second.has_upper)) {
      continue;
    }
    lower[idx] = kv.second.has_lower ? kv.second.lower : 0.0;
    upper[idx] =
        kv.second.has_upper ? kv.second.upper : std::numeric_limits<double>::infinity();
    names[idx] = std::to_string(kv.first);
    ++idx;
  }
  lower.attr("names") = names;
  upper.attr("names") = names;
  return Rcpp::List::create(Rcpp::Named("lower") = lower,
                            Rcpp::Named("upper") = upper);
}

RankedTransitionCompiler &
exact_scenario_compiler_for_ctx(const uuber::NativeContext &ctx) {
  thread_local std::unordered_map<
      const uuber::NativeContext *, std::unique_ptr<RankedTransitionCompiler>>
      cache;
  std::unique_ptr<RankedTransitionCompiler> &entry = cache[&ctx];
  if (!entry) {
    entry = std::make_unique<RankedTransitionCompiler>(ctx);
  }
  return *entry;
}

Rcpp::List scenario_to_debug_list(const uuber::NativeContext &ctx,
                                  const ExactNodeScenario &scenario,
                                  double competitor_survival,
                                  double contribution) {
  return Rcpp::List::create(
      Rcpp::Named("weight") = scenario.weight,
      Rcpp::Named("competitor_survival") = competitor_survival,
      Rcpp::Named("contribution") = contribution,
      Rcpp::Named("forced_complete_ids") =
          bitset_to_label_ids(ctx, scenario.forced_complete_bits,
                              scenario.forced_complete_bits_valid),
      Rcpp::Named("forced_survive_ids") =
          bitset_to_label_ids(ctx, scenario.forced_survive_bits,
                              scenario.forced_survive_bits_valid),
      Rcpp::Named("exact_source_times") =
          exact_times_to_numeric(scenario.time_constraints),
      Rcpp::Named("source_time_bounds") =
          bounds_to_list(scenario.time_constraints));
}

bool exact_competitor_survival_batch_impl(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const std::vector<ExactScenarioPoint> &points,
    std::vector<double> &survival_out) {
  survival_out.assign(points.size(), 1.0);
  if (points.empty() || competitor_cache.compiled_ops.empty()) {
    return true;
  }
  record_unified_outcome_exact_competitor_batch_call(
      static_cast<std::uint64_t>(points.size()),
      static_cast<std::uint64_t>(competitor_cache.compiled_ops.size()));

  std::vector<ExactScenarioPoint> mutable_points = points;
  std::vector<double> times(points.size(), 0.0);
  for (std::size_t i = 0; i < points.size(); ++i) {
    times[i] = points[i].t;
  }

  uuber::KernelBatchRuntimeState batch_runtime;
  const uuber::TrialParamsSoA *trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);

  uuber::KernelEventBatchEvalFn event_eval_batch =
      [&](int event_idx, const std::vector<double> &batch_times,
          const uuber::KernelEvalNeed &kneed,
          uuber::KernelNodeBatchValues &out_values) -> bool {
    return exact_eval_event_batch_common(
        ctx, mutable_points, batch_times, component_idx, trial_params,
        trial_type_key, trial_params_soa, event_idx, kneed, out_values, nullptr,
        [&](const uuber::IrEvent &event, EvalNeed need,
            uuber::KernelNodeBatchValues &out) -> bool {
          return event.node_idx >= 0 &&
                 exact_eval_node_batch_with_points(
                     ctx, event.node_idx, mutable_points, component_idx,
                     trial_params, trial_type_key, need, out, nullptr);
        });
  };

  uuber::KernelGuardBatchEvalFn guard_eval_batch =
      [&](const uuber::KernelOp &op, const std::vector<double> &batch_times,
          const uuber::KernelNodeBatchValues &reference_values,
          const uuber::KernelNodeBatchValues &blocker_values,
          const uuber::KernelEvalNeed &kneed,
          uuber::KernelNodeBatchValues &out_values) -> bool {
    return exact_eval_guard_batch_common(
        ctx, mutable_points, batch_times, component_idx, trial_params,
        trial_type_key, op, reference_values, blocker_values, kneed,
        out_values, nullptr, false);
  };

  uuber::KernelBatchTransitionApplyFn apply_transition_batch =
      [&](const uuber::CompetitorCompiledOp &op,
          uuber::KernelBatchRuntimeState &runtime) {
    if (op.transition_mask_begin < 0 || op.transition_mask_count <= 0) {
      return;
    }
    for (ExactScenarioPoint &point : mutable_points) {
      apply_transition_mask_words(
          ctx, op.transition_mask_begin, op.transition_mask_count,
          op.transition_invalidate_slot, point.forced_survive_bits,
          point.forced_survive_bits_valid, nullptr, nullptr);
    }
    if (op.transition_invalidate_slot >= 0) {
      uuber::invalidate_kernel_batch_runtime_from_slot(
          runtime, op.transition_invalidate_slot);
    }
  };

  return uuber::eval_kernel_competitor_product_batch_incremental(
      ctx.kernel_program, batch_runtime, competitor_cache.compiled_ops, times,
      event_eval_batch, guard_eval_batch, apply_transition_batch,
      survival_out);
}

double exact_outcome_density_at_time(
    const uuber::NativeContext &ctx, RankedTransitionCompiler &compiler,
    int node_idx, const std::vector<int> &competitor_ids, double t,
    NodeEvalState &state,
    const uuber::CompetitorClusterCacheEntry *competitor_cache) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }

  std::vector<ExactNodeScenario> scenarios;
  if (!collect_exact_node_scenarios(
          compiler, ctx, node_idx, t, state.component_idx, state.trial_params,
          state.trial_type_key, &state.time_constraints,
          (state.kernel_runtime_ready ? state.kernel_runtime_ptr : nullptr),
          state.forced_complete_bits_valid ? &state.forced_complete_bits
                                           : nullptr,
          state.forced_complete_bits_valid,
          state.forced_survive_bits_valid ? &state.forced_survive_bits : nullptr,
          state.forced_survive_bits_valid, scenarios)) {
    return 0.0;
  }

  double total = 0.0;
  for (const ExactNodeScenario &scenario : scenarios) {
    if (!std::isfinite(scenario.weight) || scenario.weight <= 0.0) {
      continue;
    }
    double weight = scenario.weight;
    if (!competitor_ids.empty()) {
      NodeEvalState scenario_state(
          ctx, t, state.component_idx, state.trial_params, state.trial_type_key,
          false, -1, &scenario.time_constraints, nullptr,
          scenario.forced_complete_bits_valid ? &scenario.forced_complete_bits
                                              : nullptr,
          scenario.forced_complete_bits_valid,
          scenario.forced_survive_bits_valid ? &scenario.forced_survive_bits
                                             : nullptr,
          scenario.forced_survive_bits_valid, nullptr);
      const double surv = competitor_survival_from_state_compiled_ops(
          ctx, competitor_ids, scenario_state, competitor_cache);
      if (!std::isfinite(surv) || surv <= 0.0) {
        continue;
      }
      weight *= surv;
    }
    if (std::isfinite(weight) && weight > 0.0) {
      total += weight;
    }
  }

  return safe_density(total);
}

} // namespace

bool exact_eval_node_batch_from_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values) {
  return exact_eval_node_batch_with_points(ctx, node_idx, points, component_idx,
                                           trial_params, trial_type_key, need,
                                           out_values, nullptr);
}

bool exact_collect_scenarios_batch(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    std::vector<ExactScenarioPoint> &scenario_points) {
  return exact_collect_scenarios_batch_impl(compiler, ctx, node_idx, times,
                                            state, scenario_points);
}

bool exact_collect_scenarios_batch_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &scenario_points) {
  return exact_collect_scenarios_batch_from_points_impl(
      compiler, ctx, node_idx, seed_points, component_idx, trial_params,
      trial_type_key, scenario_points);
}

bool exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const std::vector<ExactScenarioPoint> &points,
    std::vector<double> &survival_out) {
  return exact_competitor_survival_batch_impl(
      ctx, competitor_cache, component_idx, trial_params, trial_type_key,
      points, survival_out);
}

double exact_outcome_density_from_state(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    const uuber::CompetitorClusterCacheEntry *competitor_cache) {
  if (!std::isfinite(state.t) || state.t < 0.0) {
    return 0.0;
  }
  record_unified_outcome_exact_scalar_density_call();
  RankedTransitionCompiler &compiler = exact_scenario_compiler_for_ctx(ctx);
  return exact_outcome_density_at_time(ctx, compiler, node_idx, competitor_ids,
                                       state.t, state, competitor_cache);
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
  std::vector<ExactScenarioPoint> scenario_points;
  if (!exact_collect_scenarios_batch(compiler, ctx, node_idx, times, state,
                                     scenario_points)) {
    state.t = saved_t;
    invalidate_kernel_runtime_root(state);
    return true;
  }

  if (competitor_ids.empty()) {
    for (const ExactScenarioPoint &point : scenario_points) {
      if (!std::isfinite(point.weight) || point.weight <= 0.0 ||
          point.density_index >= density_out.size()) {
        continue;
      }
      density_out[point.density_index] += point.weight;
    }
    for (double &value : density_out) {
      value = safe_density(value);
    }
  } else if (!scenario_points.empty()) {
    const uuber::CompetitorClusterCacheEntry &cache =
        competitor_cache ? *competitor_cache
                         : fetch_competitor_cluster_cache(ctx, competitor_ids);
    std::vector<double> survival_out;
    if (exact_competitor_survival_batch(
            ctx, cache, state.component_idx, state.trial_params,
            state.trial_type_key, scenario_points, survival_out) &&
        survival_out.size() == scenario_points.size()) {
      for (std::size_t i = 0; i < scenario_points.size(); ++i) {
        const ExactScenarioPoint &point = scenario_points[i];
        const double surv = clamp_probability(survival_out[i]);
        if (!std::isfinite(surv) || surv <= 0.0) {
          continue;
        }
        density_out[point.density_index] += point.weight * surv;
      }
      for (double &value : density_out) {
        value = safe_density(value);
      }
    } else {
      density_out.assign(times.size(), 0.0);
      for (const ExactScenarioPoint &point : scenario_points) {
        NodeEvalState scenario_state(
            ctx, point.t, state.component_idx, state.trial_params,
            state.trial_type_key, false, -1, &point.time_constraints, nullptr,
            point.forced_complete_bits_valid ? &point.forced_complete_bits
                                             : nullptr,
            point.forced_complete_bits_valid,
            point.forced_survive_bits_valid ? &point.forced_survive_bits
                                            : nullptr,
            point.forced_survive_bits_valid, nullptr);
        const double surv = competitor_survival_from_state_compiled_ops(
            ctx, competitor_ids, scenario_state, competitor_cache);
        if (!std::isfinite(surv) || surv <= 0.0) {
          continue;
        }
        density_out[point.density_index] += point.weight * surv;
      }
      for (double &value : density_out) {
        value = safe_density(value);
      }
    }
  }
  state.t = saved_t;
  invalidate_kernel_runtime_root(state);
  return true;
}

// [[Rcpp::export]]
Rcpp::List native_debug_exact_scenarios_cpp(
    SEXP ctxSEXP, int node_id, double t, int component_idx = -1,
    Rcpp::Nullable<Rcpp::DataFrame> trial_rows = R_NilValue) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  const int node_idx = resolve_dense_node_idx_required(*ctx, node_id);

  std::unique_ptr<TrialParamSet> trial_params;
  if (!trial_rows.isNull()) {
    trial_params = build_trial_params_from_df(*ctx, trial_rows);
  }

  RankedTransitionCompiler compiler(*ctx);
  std::vector<ExactNodeScenario> scenarios;
  collect_exact_node_scenarios(
      compiler, *ctx, node_idx, t, component_idx, trial_params.get(),
      std::string(), nullptr, nullptr, nullptr, false, nullptr, false,
      scenarios, false);

  Rcpp::List out(scenarios.size());
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(scenarios.size()); ++i) {
    const ExactNodeScenario &scenario = scenarios[static_cast<std::size_t>(i)];
    out[i] = Rcpp::List::create(
        Rcpp::Named("weight") = scenario.weight,
        Rcpp::Named("forced_complete_ids") =
            bitset_to_label_ids(*ctx, scenario.forced_complete_bits,
                                scenario.forced_complete_bits_valid),
        Rcpp::Named("forced_survive_ids") =
            bitset_to_label_ids(*ctx, scenario.forced_survive_bits,
                                scenario.forced_survive_bits_valid),
        Rcpp::Named("exact_source_times") =
            exact_times_to_numeric(scenario.time_constraints),
        Rcpp::Named("source_time_bounds") =
            bounds_to_list(scenario.time_constraints));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List native_debug_exact_outcome_contributions_cpp(
    SEXP ctxSEXP, int node_id, double t, Rcpp::IntegerVector competitor_ids,
    int component_idx = -1,
    Rcpp::Nullable<Rcpp::DataFrame> trial_rows = R_NilValue) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  const int node_idx = resolve_dense_node_idx_required(*ctx, node_id);

  std::unique_ptr<TrialParamSet> trial_params;
  if (!trial_rows.isNull()) {
    trial_params = build_trial_params_from_df(*ctx, trial_rows);
  }

  RankedTransitionCompiler compiler(*ctx);
  std::vector<ExactNodeScenario> scenarios;
  collect_exact_node_scenarios(
      compiler, *ctx, node_idx, t, component_idx, trial_params.get(),
      std::string(), nullptr, nullptr, nullptr, false, nullptr, false,
      scenarios, true);

  std::vector<int> competitor_vec = integer_vector_to_std(competitor_ids, false);
  const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr;
  if (!competitor_vec.empty()) {
    competitor_cache = &fetch_competitor_cluster_cache(*ctx, competitor_vec);
  }

  Rcpp::List out(scenarios.size());
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(scenarios.size()); ++i) {
    const ExactNodeScenario &scenario = scenarios[static_cast<std::size_t>(i)];
    double competitor_survival = 1.0;
    double contribution = scenario.weight;
    if (!competitor_vec.empty() && std::isfinite(scenario.weight) &&
        scenario.weight > 0.0) {
      NodeEvalState scenario_state(
          *ctx, t, component_idx, trial_params.get(), std::string(), false, -1,
          &scenario.time_constraints, nullptr,
          scenario.forced_complete_bits_valid ? &scenario.forced_complete_bits
                                              : nullptr,
          scenario.forced_complete_bits_valid,
          scenario.forced_survive_bits_valid ? &scenario.forced_survive_bits
                                             : nullptr,
          scenario.forced_survive_bits_valid, nullptr);
      competitor_survival = competitor_survival_from_state_compiled_ops(
          *ctx, competitor_vec, scenario_state, competitor_cache);
      if (!std::isfinite(competitor_survival) || competitor_survival <= 0.0) {
        contribution = 0.0;
      } else {
        contribution *= competitor_survival;
      }
    }
    out[i] = scenario_to_debug_list(*ctx, scenario, competitor_survival,
                                    contribution);
  }
  return out;
}
