#include "coupling_eval.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "evaluator_internal.h"
#include "pool_math.h"
#include "runtime_stats.h"

using uuber::AccDistParams;
using uuber::LabelRef;

namespace {

using CouplingEventPayload = uuber::VectorEventRefPayload;
using CouplingGenericPayload = uuber::VectorGenericNodePayload;

inline bool coupling_acc_specialization_allowed(const uuber::NativeContext &ctx,
                                                bool include_na_donors,
                                                int outcome_idx_context) {
  if (include_na_donors) {
    return false;
  }
  if (outcome_idx_context < 0) {
    return true;
  }
  if (!ctx.ir.valid ||
      outcome_idx_context >= static_cast<int>(ctx.ir.outcomes.size())) {
    return false;
  }
  const uuber::IrOutcome &outcome =
      ctx.ir.outcomes[static_cast<std::size_t>(outcome_idx_context)];
  return outcome.alias_count == 0 && outcome.guess_count == 0;
}

struct CouplingAccRef {
  double onset{0.0};
  double q{0.0};
  AccDistParams cfg{};
  LowerBoundTransform lower_bound{};
};

struct CouplingBatchScratch {
  std::vector<double> shifted;
  std::vector<double> pdf_values;
  std::vector<double> cdf_values;
  std::vector<double> temp_density;
  std::vector<double> temp_survival;
  std::vector<double> expanded_times;
  std::vector<const uuber::TrialParamsSoA *> expanded_trial_params_soa;
};

inline bool evaluate_labelref_batch(
    const uuber::NativeContext &ctx, const CouplingEventPayload &ref,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context, const std::vector<double> &times,
    bool need_density, bool need_survival, bool need_cdf,
    std::vector<double> *density_out, std::vector<double> *survival_out,
    std::vector<double> *cdf_out, CouplingBatchScratch *scratch = nullptr,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch =
        nullptr);

inline bool build_coupling_acc_ref(
    const uuber::NativeContext &ctx, const CouplingEventPayload &ref,
    int component_idx, const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *trial_params_soa, CouplingAccRef &out) {
  if (ref.ref.acc_idx < 0 ||
      ref.ref.acc_idx >= static_cast<int>(ctx.accumulators.size())) {
    return false;
  }
  if ((ref.node_flags & (uuber::IR_NODE_FLAG_SPECIAL_DEADLINE |
                         uuber::IR_NODE_FLAG_SPECIAL_GUESS)) != 0u) {
    return false;
  }
  const uuber::NativeAccumulator &acc =
      ctx.accumulators[static_cast<std::size_t>(ref.ref.acc_idx)];
  const TrialAccumulatorParams *override =
      evaluator_get_trial_param_entry(trial_params, ref.ref.acc_idx);
  if (!evaluator_component_active_idx(acc, component_idx, override)) {
    return false;
  }
  const int onset_kind = override ? override->onset_kind : acc.onset_kind;
  if (ctx.has_chained_onsets && onset_kind != uuber::ONSET_ABSOLUTE) {
    return false;
  }
  evaluator_resolve_event_numeric_params(
      acc, ref.ref.acc_idx, override, trial_params_soa, out.onset, out.q,
      out.cfg);
  out.lower_bound = default_lower_bound_transform(out.cfg);
  return true;
}

inline bool evaluate_labelref_pool_batch(
    const uuber::NativeContext &ctx, const CouplingEventPayload &ref,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const std::vector<double> &times,
    bool need_density, bool need_survival, bool need_cdf,
    std::vector<double> *density_out, std::vector<double> *survival_out,
    std::vector<double> *cdf_out, CouplingBatchScratch *scratch,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch) {
  const int pool_idx = ref.ref.pool_idx;
  if (pool_idx < 0 || pool_idx >= static_cast<int>(ctx.pools.size()) ||
      (need_density && !density_out) || (need_survival && !survival_out) ||
      (need_cdf && !cdf_out)) {
    return false;
  }

  const std::size_t point_count = times.size();
  if (need_density) {
    density_out->assign(point_count, 0.0);
  }
  if (need_survival) {
    survival_out->assign(point_count, 1.0);
  }
  if (need_cdf) {
    cdf_out->assign(point_count, 0.0);
  }
  if (point_count == 0u) {
    return true;
  }

  const uuber::NativePool &pool = ctx.pools[static_cast<std::size_t>(pool_idx)];
  if (pool.members.empty()) {
    return true;
  }

  std::vector<std::size_t> active_member_indices;
  active_member_indices.reserve(pool.member_refs.size());
  for (std::size_t member_idx = 0; member_idx < pool.member_refs.size();
       ++member_idx) {
    const LabelRef &member_ref = pool.member_refs[member_idx];
    const int member_acc_idx = member_ref.acc_idx;
    if (member_acc_idx >= 0 &&
        member_acc_idx < static_cast<int>(ctx.accumulators.size())) {
      const uuber::NativeAccumulator &acc =
          ctx.accumulators[static_cast<std::size_t>(member_acc_idx)];
      const TrialAccumulatorParams *override =
          evaluator_get_trial_param_entry(trial_params, member_acc_idx);
      if (!evaluator_component_active_idx(acc, component_idx, override)) {
        continue;
      }
    }
    active_member_indices.push_back(member_idx);
  }
  if (active_member_indices.empty()) {
    return true;
  }

  const bool need_member_survival = need_survival || need_cdf || need_density;
  if (!need_member_survival) {
    return true;
  }

  std::vector<double> density_local;
  std::vector<double> survival_local;
  std::vector<double> &temp_density =
      scratch ? scratch->temp_density : density_local;
  std::vector<double> &temp_survival =
      scratch ? scratch->temp_survival : survival_local;
  std::vector<double> member_density_matrix;
  std::vector<double> member_survival_matrix;
  std::vector<double> point_density;
  std::vector<double> point_survival;

  member_survival_matrix.assign(active_member_indices.size() * point_count, 1.0);
  if (need_density) {
    member_density_matrix.assign(active_member_indices.size() * point_count, 0.0);
  }
  point_survival.resize(active_member_indices.size());
  if (need_density) {
    point_density.resize(active_member_indices.size());
  }

  for (std::size_t member_pos = 0; member_pos < active_member_indices.size();
       ++member_pos) {
    CouplingEventPayload member_payload{};
    member_payload.ref =
        pool.member_refs[active_member_indices[member_pos]];
    temp_density.clear();
    temp_survival.clear();
    if (!evaluate_labelref_batch(
            ctx, member_payload, component_idx, trial_params, trial_type_key,
            false, -1, times, need_density, true, false,
            need_density ? &temp_density : nullptr, &temp_survival, nullptr,
            scratch, trial_params_soa_batch) ||
        temp_survival.size() != point_count ||
        (need_density && temp_density.size() != point_count)) {
      return false;
    }
    const std::size_t offset = member_pos * point_count;
    std::copy(temp_survival.begin(), temp_survival.end(),
              member_survival_matrix.begin() + static_cast<std::ptrdiff_t>(offset));
    if (need_density) {
      std::copy(temp_density.begin(), temp_density.end(),
                member_density_matrix.begin() +
                    static_cast<std::ptrdiff_t>(offset));
    }
  }

  for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
    for (std::size_t member_pos = 0; member_pos < active_member_indices.size();
         ++member_pos) {
      const std::size_t offset = member_pos * point_count + point_idx;
      point_survival[member_pos] =
          clamp_probability(member_survival_matrix[offset]);
      if (need_density) {
        point_density[member_pos] =
            safe_density(member_density_matrix[offset]);
      }
    }
    if (need_density) {
      const double density =
          pool_density_fast(point_density, point_survival, pool.k);
      (*density_out)[point_idx] =
          (std::isfinite(density) && density > 0.0) ? safe_density(density)
                                                    : 0.0;
    }
    if (need_survival || need_cdf) {
      const double survival = pool_survival_fast(point_survival, pool.k);
      const double clamped_survival = clamp_probability(survival);
      if (need_survival) {
        (*survival_out)[point_idx] = clamped_survival;
      }
      if (need_cdf) {
        (*cdf_out)[point_idx] = clamp_probability(1.0 - clamped_survival);
      }
    }
  }

  record_unified_outcome_generic_poolref_batch_fastpath_call();
  record_unified_outcome_generic_labelref_batch_fastpath_call();
  return true;
}

inline void expand_time_batch_by_trial_params(
    const std::vector<double> &times,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    std::vector<double> &expanded_times,
    std::vector<const uuber::TrialParamsSoA *> &expanded_trial_params_soa) {
  expanded_times.clear();
  expanded_trial_params_soa.clear();
  if (times.empty() || trial_params_soa_batch.empty()) {
    return;
  }
  expanded_times.reserve(times.size() * trial_params_soa_batch.size());
  expanded_trial_params_soa.reserve(times.size() * trial_params_soa_batch.size());
  for (const uuber::TrialParamsSoA *point_soa : trial_params_soa_batch) {
    for (double time_value : times) {
      expanded_times.push_back(time_value);
      expanded_trial_params_soa.push_back(point_soa);
    }
  }
}

inline bool coupling_competitor_node_allowed_idx(
    const uuber::NativeContext &ctx, int node_id, int component_idx) {
  if (component_idx < 0) {
    return true;
  }
  const int dense_node = resolve_dense_node_idx(ctx, node_id);
  if (dense_node < 0) {
    return true;
  }
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
  if (out_it == ctx.ir.node_idx_to_outcomes.end() || out_it->second.empty()) {
    return true;
  }
  for (int out_idx : out_it->second) {
    if (ir_outcome_allows_component(ctx, out_idx, component_idx)) {
      return true;
    }
  }
  return false;
}

inline const std::vector<int> &coupling_filter_competitor_ids(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    int component_idx, std::vector<int> &scratch) {
  if (competitor_ids.empty() || component_idx < 0) {
    return competitor_ids;
  }
  bool all_allowed = true;
  for (int node_id : competitor_ids) {
    if (!coupling_competitor_node_allowed_idx(ctx, node_id, component_idx)) {
      all_allowed = false;
      break;
    }
  }
  if (all_allowed) {
    return competitor_ids;
  }
  scratch.clear();
  scratch.reserve(competitor_ids.size());
  for (int node_id : competitor_ids) {
    if (coupling_competitor_node_allowed_idx(ctx, node_id, component_idx)) {
      scratch.push_back(node_id);
    }
  }
  return scratch;
}

inline bool coupling_evaluate_outcome_density_batch(
    const uuber::NativeContext &ctx, int outcome_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_label_context_idx,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out) {
  density_out.assign(times.size(), 0.0);
  if (outcome_idx < 0 ||
      outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return true;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
  if (info.node_id < 0) {
    return true;
  }
  std::vector<int> comp_filtered;
  const std::vector<int> &comp_use = coupling_filter_competitor_ids(
      ctx, info.competitor_ids, component_idx, comp_filtered);
  return evaluator_node_density_with_competitors_batch_internal(
      ctx, info.node_id, times, component_idx, nullptr, false, nullptr, false,
      comp_use, trial_params, trial_type_key, include_na_donors,
      outcome_label_context_idx, density_out, nullptr, nullptr,
      trial_params_soa_batch);
}

inline bool coupling_component_guess_matches_target(
    const uuber::NativeContext &ctx, int target_label_id, int target_outcome_idx,
    bool target_is_guess, int component_idx) {
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
    return target_is_guess;
  }
  if (guess.target_outcome_idx >= 0) {
    return target_outcome_idx == guess.target_outcome_idx;
  }
  if (guess.target_label_id != NA_INTEGER) {
    return target_label_id != NA_INTEGER &&
           target_label_id == guess.target_label_id;
  }
  return false;
}

inline bool coupling_outcome_has_alias_or_guess_donors(
    const uuber::NativeContext &ctx, int outcome_idx, bool include_na_donors) {
  if (outcome_idx < 0 ||
      outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return false;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
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

inline bool coupling_accumulate_component_guess_density_batch(
    const uuber::NativeContext &ctx, int target_label_id, int target_outcome_idx,
    bool target_is_guess, const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out) {
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

  std::vector<double> donor_density;
  for (const auto &entry : guess.keep_weights) {
    const int donor_outcome_idx = entry.first;
    const double release = 1.0 - entry.second;
    if (donor_outcome_idx < 0 || !(release > 0.0)) {
      continue;
    }
    if (!coupling_evaluate_outcome_density_batch(
            ctx, donor_outcome_idx, times, component_idx, trial_params,
            trial_type_key, false, donor_outcome_idx, trial_params_soa_batch,
            donor_density) ||
        donor_density.size() != density_out.size()) {
      return false;
    }
    for (std::size_t i = 0; i < density_out.size(); ++i) {
      density_out[i] += release * safe_density(donor_density[i]);
    }
  }
  return true;
}

inline bool coupling_accumulate_outcome_alias_density_batch(
    const uuber::NativeContext &ctx, int target_outcome_idx,
    const std::vector<double> &times, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch,
    std::vector<double> &density_out) {
  density_out.assign(times.size(), 0.0);
  if (target_outcome_idx < 0 ||
      target_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return true;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(target_outcome_idx)];
  std::vector<double> donor_density;
  for (int source_idx : info.alias_sources) {
    if (!coupling_evaluate_outcome_density_batch(
            ctx, source_idx, times, component_idx, trial_params,
            trial_type_key, include_na_donors, source_idx,
            trial_params_soa_batch, donor_density) ||
        donor_density.size() != density_out.size()) {
      return false;
    }
    for (std::size_t i = 0; i < density_out.size(); ++i) {
      density_out[i] += safe_density(donor_density[i]);
    }
  }
  for (const auto &donor : info.guess_donors) {
    if (donor.rt_policy == "na" && !include_na_donors) {
      continue;
    }
    if (donor.outcome_idx < 0) {
      continue;
    }
    if (!coupling_evaluate_outcome_density_batch(
            ctx, donor.outcome_idx, times, component_idx, trial_params,
            trial_type_key, include_na_donors, donor.outcome_idx,
            trial_params_soa_batch, donor_density) ||
        donor_density.size() != density_out.size()) {
      return false;
    }
    for (std::size_t i = 0; i < density_out.size(); ++i) {
      density_out[i] += donor.weight * safe_density(donor_density[i]);
    }
  }
  return true;
}

inline bool build_coupling_time_batch_mass_infinite_segment(
    double x_lo, double x_hi, uuber::TimeBatch &out) {
  constexpr double kXUpperCap = 1.0 - 1e-12;
  if (!(x_hi > x_lo) || x_lo < 0.0) {
    return false;
  }
  x_hi = std::min(x_hi, kXUpperCap);
  x_lo = std::min(x_lo, x_hi);
  if (!(x_hi > x_lo) || x_hi <= 0.0) {
    return false;
  }
  uuber::TimeBatch x_batch = uuber::build_time_batch(x_lo, x_hi);
  if (x_batch.nodes.empty() || x_batch.nodes.size() != x_batch.weights.size()) {
    return false;
  }
  out.lower = 0.0;
  out.upper = std::numeric_limits<double>::infinity();
  out.nodes.clear();
  out.weights.clear();
  out.nodes.reserve(x_batch.nodes.size());
  out.weights.reserve(x_batch.weights.size());
  for (std::size_t i = 0; i < x_batch.nodes.size(); ++i) {
    const double x = x_batch.nodes[i];
    const double wx = x_batch.weights[i];
    if (!std::isfinite(x) || !std::isfinite(wx) || wx <= 0.0 || x < 0.0 ||
        x >= 1.0) {
      continue;
    }
    const double one_minus_x = 1.0 - x;
    if (!(one_minus_x > 0.0)) {
      continue;
    }
    const double t = x / one_minus_x;
    const double wt = wx / (one_minus_x * one_minus_x);
    if (!std::isfinite(t) || !std::isfinite(wt) || wt <= 0.0) {
      continue;
    }
    out.nodes.push_back(t);
    out.weights.push_back(wt);
  }
  return !out.nodes.empty() && out.nodes.size() == out.weights.size();
}

inline bool evaluate_coupling_acc_batch(
    const CouplingAccRef &ref, const std::vector<double> &times,
    bool need_density, bool need_survival, bool need_cdf,
    std::vector<double> *density_out, std::vector<double> *survival_out,
    std::vector<double> *cdf_out, CouplingBatchScratch *scratch = nullptr) {
  const std::size_t n = times.size();
  if (n == 0) {
    if (need_density && density_out) {
      density_out->clear();
    }
    if (need_survival && survival_out) {
      survival_out->clear();
    }
    if (need_cdf && cdf_out) {
      cdf_out->clear();
    }
    return true;
  }
  if ((need_density && !density_out) || (need_survival && !survival_out) ||
      (need_cdf && !cdf_out)) {
    return false;
  }

  if (need_density) {
    density_out->resize(n);
  }
  if (need_survival) {
    survival_out->resize(n);
  }
  if (need_cdf) {
    cdf_out->resize(n);
  }

  std::vector<double> shifted_local;
  std::vector<double> pdf_values_local;
  std::vector<double> cdf_values_local;
  std::vector<double> &shifted = scratch ? scratch->shifted : shifted_local;
  shifted.resize(n);

  const double onset_eff = total_onset_with_t0(ref.onset, ref.cfg);
  const double success_prob = 1.0 - ref.q;
  for (std::size_t i = 0; i < n; ++i) {
    shifted[i] = times[i] - onset_eff;
  }

  const double *pdf_values_data = nullptr;
  const double *cdf_values_data = nullptr;
  if (need_density) {
    std::vector<double> &pdf_values =
        scratch ? scratch->pdf_values : pdf_values_local;
    pdf_values.resize(n);
    eval_pdf_vec_with_lower_bound(ref.cfg, ref.lower_bound, shifted.data(), n,
                                  pdf_values.data());
    pdf_values_data = pdf_values.data();
  }
  if (need_survival || need_cdf) {
    std::vector<double> &cdf_values =
        scratch ? scratch->cdf_values : cdf_values_local;
    cdf_values.resize(n);
    eval_cdf_vec_with_lower_bound(ref.cfg, ref.lower_bound, shifted.data(), n,
                                  cdf_values.data());
    cdf_values_data = cdf_values.data();
  }

  for (std::size_t i = 0; i < n; ++i) {
    const double t = times[i];

    if (need_density) {
      double dens = 0.0;
      if (std::isfinite(t) && t >= 0.0 && t >= onset_eff &&
          success_prob > 0.0) {
        dens = success_prob * safe_density(pdf_values_data[i]);
      }
      (*density_out)[i] = safe_density(dens);
    }

    if (need_survival) {
      double surv = 0.0;
      if (std::isfinite(t)) {
        if (t < 0.0 || t < onset_eff) {
          surv = 1.0;
        } else {
          const double cdf = clamp_probability(cdf_values_data[i]);
          const double surv_underlying = clamp(1.0 - cdf, 0.0, 1.0);
          surv = ref.q + success_prob * surv_underlying;
        }
      }
      (*survival_out)[i] = clamp_probability(surv);
    }

    if (need_cdf) {
      double cdf_success = 0.0;
      if (!std::isfinite(t)) {
        cdf_success = success_prob;
      } else if (t < 0.0 || t < onset_eff) {
        cdf_success = 0.0;
      } else {
        cdf_success = clamp_probability(success_prob * cdf_values_data[i]);
      }
      (*cdf_out)[i] = clamp_probability(cdf_success);
    }
  }
  return true;
}

inline bool evaluate_coupling_acc_batch_param_matrix(
    const CouplingAccRef &base_ref, int acc_idx,
    const uuber::TrialParamsSoA *base_trial_params_soa,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    const std::vector<double> &times, bool need_density, bool need_survival,
    bool need_cdf, std::vector<double> *density_out,
    std::vector<double> *survival_out, std::vector<double> *cdf_out,
    CouplingBatchScratch *scratch = nullptr) {
  const std::size_t n = times.size();
  if (n == 0) {
    if (need_density && density_out) {
      density_out->clear();
    }
    if (need_survival && survival_out) {
      survival_out->clear();
    }
    if (need_cdf && cdf_out) {
      cdf_out->clear();
    }
    return true;
  }
  if (trial_params_soa_batch.size() != n || !base_trial_params_soa ||
      !base_trial_params_soa->valid) {
    return false;
  }
  for (const uuber::TrialParamsSoA *point_soa : trial_params_soa_batch) {
    if (!point_soa ||
        !same_acc_batch_params_except_q(*base_trial_params_soa, *point_soa,
                                        acc_idx)) {
      return false;
    }
  }

  if ((need_density && !density_out) || (need_survival && !survival_out) ||
      (need_cdf && !cdf_out)) {
    return false;
  }
  if (need_density) {
    density_out->resize(n);
  }
  if (need_survival) {
    survival_out->resize(n);
  }
  if (need_cdf) {
    cdf_out->resize(n);
  }

  std::vector<double> shifted_local;
  std::vector<double> pdf_values_local;
  std::vector<double> cdf_values_local;
  std::vector<double> &shifted = scratch ? scratch->shifted : shifted_local;
  shifted.resize(n);

  const double onset_eff = total_onset_with_t0(base_ref.onset, base_ref.cfg);
  for (std::size_t i = 0; i < n; ++i) {
    shifted[i] = times[i] - onset_eff;
  }

  const double *pdf_values_data = nullptr;
  const double *cdf_values_data = nullptr;
  if (need_density) {
    std::vector<double> &pdf_values =
        scratch ? scratch->pdf_values : pdf_values_local;
    pdf_values.resize(n);
    eval_pdf_vec_with_lower_bound(base_ref.cfg, base_ref.lower_bound,
                                  shifted.data(), n, pdf_values.data());
    pdf_values_data = pdf_values.data();
  }
  if (need_survival || need_cdf) {
    std::vector<double> &cdf_values =
        scratch ? scratch->cdf_values : cdf_values_local;
    cdf_values.resize(n);
    eval_cdf_vec_with_lower_bound(base_ref.cfg, base_ref.lower_bound,
                                  shifted.data(), n, cdf_values.data());
    cdf_values_data = cdf_values.data();
  }

  for (std::size_t i = 0; i < n; ++i) {
    const uuber::TrialParamsSoA *point_soa = trial_params_soa_batch[i];
    double q = base_ref.q;
    if (point_soa && acc_idx >= 0 && acc_idx < point_soa->n_acc &&
        static_cast<std::size_t>(acc_idx) < point_soa->q.size()) {
      q = point_soa->q[static_cast<std::size_t>(acc_idx)];
    }
    const double success_prob = 1.0 - q;
    const double t = times[i];

    if (need_density) {
      double dens = 0.0;
      if (std::isfinite(t) && t >= 0.0 && t >= onset_eff &&
          success_prob > 0.0) {
        dens = success_prob * safe_density(pdf_values_data[i]);
      }
      (*density_out)[i] = safe_density(dens);
    }

    if (need_survival) {
      double surv = 0.0;
      if (std::isfinite(t)) {
        if (t < 0.0 || t < onset_eff) {
          surv = 1.0;
        } else {
          const double cdf = clamp_probability(cdf_values_data[i]);
          const double surv_underlying = clamp(1.0 - cdf, 0.0, 1.0);
          surv = q + success_prob * surv_underlying;
        }
      }
      (*survival_out)[i] = clamp_probability(surv);
    }

    if (need_cdf) {
      double cdf_success = 0.0;
      if (!std::isfinite(t)) {
        cdf_success = success_prob;
      } else if (t < 0.0 || t < onset_eff) {
        cdf_success = 0.0;
      } else {
        cdf_success = clamp_probability(success_prob * cdf_values_data[i]);
      }
      (*cdf_out)[i] = clamp_probability(cdf_success);
    }
  }
  return true;
}

inline bool evaluate_labelref_node_batch(
    const uuber::NativeContext &ctx,
    const CouplingEventPayload &ref, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    const std::vector<double> &times, bool need_density, bool need_survival,
    bool need_cdf, std::vector<double> *density_out,
    std::vector<double> *survival_out, std::vector<double> *cdf_out,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch) {
  if (ref.node_id < 0 || !ctx.kernel_program.valid) {
    return false;
  }
  EvalNeed need = static_cast<EvalNeed>(0u);
  if (need_density) {
    need = need | EvalNeed::kDensity;
  }
  if (need_survival) {
    need = need | EvalNeed::kSurvival;
  }
  if (need_cdf) {
    need = need | EvalNeed::kCDF;
  }
  if (need == static_cast<EvalNeed>(0u)) {
    return false;
  }

  const int node_idx = resolve_dense_node_idx_required(ctx, ref.node_id);
  NodeEvalState state(ctx, 0.0, component_idx, trial_params, trial_type_key,
                      include_na_donors, outcome_idx_context);
  state.trial_params_soa_batch = trial_params_soa_batch;
  uuber::KernelBatchRuntimeState batch_runtime;
  uuber::KernelNodeBatchValues values;
  if (!evaluator_eval_node_batch_with_state_dense(node_idx, times, state, need,
                                                  batch_runtime, values)) {
    return false;
  }

  if (need_density) {
    if (!density_out || values.density.size() != times.size()) {
      return false;
    }
    density_out->resize(times.size());
    for (std::size_t i = 0; i < times.size(); ++i) {
      (*density_out)[i] = safe_density(values.density[i]);
    }
  }
  if (need_survival) {
    if (!survival_out || values.survival.size() != times.size()) {
      return false;
    }
    survival_out->resize(times.size());
    for (std::size_t i = 0; i < times.size(); ++i) {
      (*survival_out)[i] = clamp_probability(values.survival[i]);
    }
  }
  if (need_cdf) {
    if (!cdf_out || values.cdf.size() != times.size()) {
      return false;
    }
    cdf_out->resize(times.size());
    for (std::size_t i = 0; i < times.size(); ++i) {
      (*cdf_out)[i] = clamp_probability(values.cdf[i]);
    }
  }
  record_unified_outcome_generic_labelref_batch_fastpath_call();
  return true;
}

inline bool evaluate_labelref_batch(
    const uuber::NativeContext &ctx, const CouplingEventPayload &ref,
    int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    const std::vector<double> &times, bool need_density, bool need_survival,
    bool need_cdf, std::vector<double> *density_out,
    std::vector<double> *survival_out, std::vector<double> *cdf_out,
    CouplingBatchScratch *scratch,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch) {
  if ((need_density && !density_out) || (need_survival && !survival_out) ||
      (need_cdf && !cdf_out)) {
    return false;
  }
  const std::size_t n = times.size();
  if (need_density) {
    density_out->assign(n, 0.0);
  }
  if (need_survival) {
    survival_out->assign(n, 1.0);
  }
  if (need_cdf) {
    cdf_out->assign(n, 0.0);
  }
  if (n == 0) {
    return true;
  }

  if (trial_params_soa_batch != nullptr &&
      trial_params_soa_batch->size() != n) {
    return false;
  }

  const uuber::TrialParamsSoA *trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);
  const bool is_special_deadline =
      (ref.node_flags & uuber::IR_NODE_FLAG_SPECIAL_DEADLINE) != 0u;
  const bool is_special_guess =
      (ref.node_flags & uuber::IR_NODE_FLAG_SPECIAL_GUESS) != 0u;
  if (is_special_deadline) {
    if (need_density) {
      std::fill(density_out->begin(), density_out->end(), 0.0);
    }
    if (need_survival) {
      for (std::size_t i = 0; i < n; ++i) {
        (*survival_out)[i] =
            (std::isfinite(times[i]) && times[i] >= 0.0) ? 0.0 : 1.0;
      }
    }
    if (need_cdf) {
      for (std::size_t i = 0; i < n; ++i) {
        (*cdf_out)[i] =
            (std::isfinite(times[i]) && times[i] >= 0.0) ? 1.0 : 0.0;
      }
    }
    return true;
  }

  std::vector<double> donor_density;
  const int target_outcome_idx =
      (outcome_idx_context >= 0) ? outcome_idx_context : ref.ref.outcome_idx;
  int target_label_id = ref.ref.label_id;
  if (outcome_idx_context >= 0 &&
      outcome_idx_context < static_cast<int>(ctx.outcome_label_ids.size())) {
    target_label_id =
        ctx.outcome_label_ids[static_cast<std::size_t>(outcome_idx_context)];
  }
  const bool need_donor_density =
      need_density &&
      (is_special_guess ||
       coupling_component_guess_matches_target(
           ctx, target_label_id, target_outcome_idx, is_special_guess,
           component_idx) ||
       coupling_outcome_has_alias_or_guess_donors(ctx, target_outcome_idx,
                                                  include_na_donors));
  if (need_donor_density) {
    std::vector<double> component_guess_density;
    std::vector<double> outcome_alias_density;
    if (!coupling_accumulate_component_guess_density_batch(
            ctx, target_label_id, target_outcome_idx, is_special_guess, times,
            component_idx, trial_params, trial_type_key, trial_params_soa_batch,
            component_guess_density) ||
        !coupling_accumulate_outcome_alias_density_batch(
            ctx, target_outcome_idx, times, component_idx, trial_params,
            trial_type_key, include_na_donors, trial_params_soa_batch,
            outcome_alias_density) ||
        component_guess_density.size() != n ||
        outcome_alias_density.size() != n) {
      return false;
    }
    donor_density.resize(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      donor_density[i] = safe_density(component_guess_density[i] +
                                      outcome_alias_density[i]);
    }
  }
  if (is_special_guess) {
    if (need_density) {
      *density_out = donor_density;
    }
    if (need_survival) {
      std::fill(survival_out->begin(), survival_out->end(), 0.0);
    }
    if (need_cdf) {
      std::fill(cdf_out->begin(), cdf_out->end(), 1.0);
    }
    return true;
  }

  if (ref.ref.pool_idx >= 0 &&
      evaluate_labelref_pool_batch(
          ctx, ref, component_idx, trial_params, trial_type_key, times,
          need_density, need_survival, need_cdf, density_out, survival_out,
          cdf_out, scratch, trial_params_soa_batch)) {
    if (need_density && donor_density.size() == n) {
      for (std::size_t i = 0; i < n; ++i) {
        (*density_out)[i] = safe_density((*density_out)[i] + donor_density[i]);
      }
    }
    return true;
  }

  CouplingAccRef ref_acc;
  if (build_coupling_acc_ref(ctx, ref, component_idx, trial_params,
                             trial_params_soa, ref_acc)) {
    bool acc_batch_ok = false;
    if (trial_params_soa_batch != nullptr &&
        evaluate_coupling_acc_batch_param_matrix(
            ref_acc, ref.ref.acc_idx, trial_params_soa, *trial_params_soa_batch,
            times, need_density, need_survival, need_cdf, density_out,
            survival_out, cdf_out, scratch)) {
      acc_batch_ok = true;
    } else if (trial_params_soa_batch == nullptr &&
               evaluate_coupling_acc_batch(ref_acc, times, need_density,
                                           need_survival, need_cdf, density_out,
                                           survival_out, cdf_out, scratch)) {
      acc_batch_ok = true;
    }
    if (acc_batch_ok) {
      if (need_density && donor_density.size() == n) {
        for (std::size_t i = 0; i < n; ++i) {
          (*density_out)[i] =
              safe_density((*density_out)[i] + donor_density[i]);
        }
      }
      record_unified_outcome_generic_labelref_batch_fastpath_call();
      return true;
    }
  }

  if (evaluate_labelref_node_batch(
          ctx, ref, component_idx, trial_params, trial_type_key,
          include_na_donors, outcome_idx_context, times, need_density,
          need_survival, need_cdf, density_out, survival_out, cdf_out,
          trial_params_soa_batch)) {
    return true;
  }
  // Valid coupling IR always resolves label refs back to an event node.
  // Missing all batched routes here means malformed metadata, not a reason to
  // silently re-enter the scalar evaluator.
  return false;
}

struct GenericCouplingProviderRuntime {
  CouplingEventPayload target_ref{};
  std::vector<CouplingEventPayload> competitor_refs;
  CouplingBatchScratch labelref_scratch;
  uuber::KernelBatchRuntimeState noderef_kernel_batch_runtime;
};

inline bool generic_coupling_runtime_has_labelref_fastpath(
    const GenericCouplingProviderRuntime &runtime) {
  const LabelRef &ref = runtime.target_ref.ref;
  return ref.acc_idx >= 0 || ref.pool_idx >= 0 || ref.outcome_idx >= 0;
}

inline void record_generic_coupling_provider_usage(
    const GenericCouplingProviderRuntime &runtime) {
  if (generic_coupling_runtime_has_labelref_fastpath(runtime)) {
    record_unified_outcome_generic_labelref_batch_fastpath_call();
  } else {
    record_unified_outcome_generic_noderef_batch_call();
  }
}

inline bool multiply_generic_survival_product(
    std::vector<double> &terms, const std::vector<double> &survival) {
  if (terms.size() != survival.size()) {
    return false;
  }
  for (std::size_t i = 0; i < terms.size(); ++i) {
    terms[i] *= clamp_probability(survival[i]);
  }
  return true;
}

inline bool finalize_generic_terms_from_density(
    const std::vector<double> &density, std::vector<double> &terms) {
  if (terms.size() != density.size()) {
    return false;
  }
  for (std::size_t i = 0; i < terms.size(); ++i) {
    const double fx = safe_density(density[i]);
    const double surv_prod = clamp_probability(terms[i]);
    if (!std::isfinite(fx) || fx <= 0.0 || surv_prod <= 0.0) {
      terms[i] = 0.0;
    } else {
      terms[i] = fx * surv_prod;
    }
  }
  return true;
}

inline bool generic_coupling_event_refs_resolve(
    const uuber::NativeContext &ctx,
    const CouplingGenericPayload &generic, CouplingEventPayload &target_ref,
    std::vector<CouplingEventPayload> &competitor_refs) {
  const int node_idx = resolve_dense_node_idx(ctx, generic.node_id);
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  const uuber::IrNode &target_node =
      ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  if (target_node.event_idx < 0) {
    return false;
  }
  target_ref = CouplingEventPayload{};
  target_ref.ref = node_label_ref(ctx, target_node);
  target_ref.node_id = generic.node_id;
  target_ref.node_flags = target_node.flags;
  if (!(target_ref.ref.acc_idx >= 0 || target_ref.ref.pool_idx >= 0 ||
        target_ref.ref.outcome_idx >= 0)) {
    return false;
  }
  competitor_refs.clear();
  competitor_refs.reserve(generic.competitor_node_ids.size());
  for (int comp_node_id : generic.competitor_node_ids) {
    const int comp_idx = resolve_dense_node_idx(ctx, comp_node_id);
    if (comp_idx < 0 || comp_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      competitor_refs.clear();
      return false;
    }
    const uuber::IrNode &comp_node =
        ctx.ir.nodes[static_cast<std::size_t>(comp_idx)];
    if (comp_node.event_idx < 0) {
      competitor_refs.clear();
      return false;
    }
    CouplingEventPayload comp_ref;
    comp_ref.ref = node_label_ref(ctx, comp_node);
    comp_ref.node_id = comp_node_id;
    comp_ref.node_flags = comp_node.flags;
    if (!(comp_ref.ref.acc_idx >= 0 || comp_ref.ref.pool_idx >= 0 ||
          comp_ref.ref.outcome_idx >= 0)) {
      competitor_refs.clear();
      return false;
    }
    competitor_refs.push_back(comp_ref);
  }
  return true;
}

inline bool evaluate_generic_terms_resolved_provider_batch(
    const uuber::NativeContext &ctx,
    const CouplingGenericPayload &generic,
    GenericCouplingProviderRuntime &runtime, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid, const std::vector<double> &nodes,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    std::vector<double> &terms) {
  if (trial_params_soa_batch.empty()) {
    terms.clear();
    return true;
  }
  expand_time_batch_by_trial_params(
      nodes, trial_params_soa_batch, runtime.labelref_scratch.expanded_times,
      runtime.labelref_scratch.expanded_trial_params_soa);
  const std::size_t point_count = trial_params_soa_batch.size();
  const std::size_t node_count = nodes.size();
  if (runtime.labelref_scratch.expanded_times.size() !=
      point_count * node_count) {
    return false;
  }

  if (!generic_coupling_runtime_has_labelref_fastpath(runtime)) {
    uuber::KernelBatchRuntimeState *kernel_batch_runtime =
        ctx.kernel_program.valid ? &runtime.noderef_kernel_batch_runtime
                                 : nullptr;
    return evaluator_node_density_with_competitors_batch_internal(
        ctx, generic.node_id, runtime.labelref_scratch.expanded_times,
        component_idx, forced_complete_bits, forced_complete_bits_valid,
        forced_survive_bits, forced_survive_bits_valid,
        generic.competitor_node_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, terms, nullptr,
        kernel_batch_runtime,
        &runtime.labelref_scratch.expanded_trial_params_soa);
  }

  const uuber::TrialParamsSoA *trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);
  const bool can_specialize = coupling_acc_specialization_allowed(
      ctx, include_na_donors, outcome_idx_context);
  if (can_specialize) {
    CouplingAccRef target_acc;
    if (build_coupling_acc_ref(ctx, runtime.target_ref, component_idx,
                               trial_params, trial_params_soa, target_acc)) {
      std::vector<CouplingAccRef> competitor_accs;
      competitor_accs.reserve(runtime.competitor_refs.size());
      bool all_competitors_specialized = true;
      for (const CouplingEventPayload &comp_ref :
           runtime.competitor_refs) {
        CouplingAccRef comp_acc;
        if (!build_coupling_acc_ref(ctx, comp_ref, component_idx, trial_params,
                                    trial_params_soa, comp_acc)) {
          all_competitors_specialized = false;
          break;
        }
        competitor_accs.push_back(comp_acc);
      }
      if (all_competitors_specialized) {
        if (!evaluate_coupling_acc_batch_param_matrix(
                target_acc, runtime.target_ref.ref.acc_idx, trial_params_soa,
                runtime.labelref_scratch.expanded_trial_params_soa,
                runtime.labelref_scratch.expanded_times, true, false, false,
                &runtime.labelref_scratch.pdf_values, nullptr, nullptr,
                &runtime.labelref_scratch)) {
          return false;
        }
        terms.assign(runtime.labelref_scratch.pdf_values.size(), 1.0);
        for (std::size_t comp_idx = 0; comp_idx < competitor_accs.size();
             ++comp_idx) {
          if (!evaluate_coupling_acc_batch_param_matrix(
                  competitor_accs[comp_idx],
                  runtime.competitor_refs[comp_idx].ref.acc_idx,
                  trial_params_soa,
                  runtime.labelref_scratch.expanded_trial_params_soa,
                  runtime.labelref_scratch.expanded_times, false, true, false,
                  nullptr, &runtime.labelref_scratch.temp_survival, nullptr,
                  &runtime.labelref_scratch) ||
              !multiply_generic_survival_product(
                  terms, runtime.labelref_scratch.temp_survival)) {
            return false;
          }
        }
        return finalize_generic_terms_from_density(
            runtime.labelref_scratch.pdf_values, terms);
      }
    }
  }

  if (!evaluate_labelref_batch(
          ctx, runtime.target_ref, component_idx, trial_params, trial_type_key,
          include_na_donors, outcome_idx_context,
          runtime.labelref_scratch.expanded_times, true, false, false,
          &runtime.labelref_scratch.pdf_values, nullptr, nullptr,
          &runtime.labelref_scratch,
          &runtime.labelref_scratch.expanded_trial_params_soa)) {
    return false;
  }
  terms.assign(runtime.labelref_scratch.pdf_values.size(), 1.0);
  for (const CouplingEventPayload &comp_ref :
       runtime.competitor_refs) {
    if (!evaluate_labelref_batch(
            ctx, comp_ref, component_idx, trial_params, trial_type_key,
            include_na_donors, outcome_idx_context,
            runtime.labelref_scratch.expanded_times, false, true, false,
            nullptr, &runtime.labelref_scratch.temp_survival, nullptr,
            &runtime.labelref_scratch,
            &runtime.labelref_scratch.expanded_trial_params_soa) ||
        !multiply_generic_survival_product(
            terms, runtime.labelref_scratch.temp_survival)) {
      return false;
    }
  }
  return finalize_generic_terms_from_density(runtime.labelref_scratch.pdf_values,
                                             terms);
}

inline GenericCouplingProviderRuntime build_generic_coupling_provider_runtime(
    const uuber::NativeContext &ctx,
    const CouplingGenericPayload &generic,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid) {
  if (generic.node_id < 0) {
    Rcpp::stop("IR generic coupling provider resolve failed: missing node_id");
  }
  const int target_idx = resolve_dense_node_idx(ctx, generic.node_id);
  if (target_idx < 0 || target_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    Rcpp::stop("IR generic coupling provider resolve failed: unknown node_id=%d",
               generic.node_id);
  }

  GenericCouplingProviderRuntime runtime;
  const bool forced_empty =
      !forced_bits_any(forced_complete_bits, forced_complete_bits_valid) &&
      !forced_bits_any(forced_survive_bits, forced_survive_bits_valid);
  if (!forced_empty) {
    return runtime;
  }

  if (!generic_coupling_event_refs_resolve(ctx, generic, runtime.target_ref,
                                           runtime.competitor_refs)) {
    return runtime;
  }
  return runtime;
}

[[maybe_unused]] inline double integrate_node_density_batch(
    const uuber::TimeBatch &batch, const std::vector<double> &density,
    const std::vector<double> *survival_product = nullptr) {
  if (batch.nodes.empty() || batch.nodes.size() != batch.weights.size() ||
      density.size() != batch.nodes.size()) {
    return 0.0;
  }
  if (survival_product && survival_product->size() != density.size()) {
    return 0.0;
  }
  double total = 0.0;
  for (std::size_t i = 0; i < density.size(); ++i) {
    const double w = batch.weights[i];
    if (!std::isfinite(w) || w <= 0.0) {
      continue;
    }
    double term = safe_density(density[i]);
    if (term <= 0.0) {
      continue;
    }
    if (survival_product) {
      const double surv = clamp_probability((*survival_product)[i]);
      if (surv <= 0.0) {
        continue;
      }
      term *= surv;
    }
    if (std::isfinite(term) && term > 0.0) {
      total += w * term;
    }
  }
  if (!std::isfinite(total) || total <= 0.0) {
    return 0.0;
  }
  return clamp_probability(total);
}

inline bool integrate_node_density_batch_per_point(
    const uuber::TimeBatch &batch, const std::vector<double> &density,
    std::size_t point_count, std::vector<double> &out) {
  if (batch.nodes.empty() || batch.nodes.size() != batch.weights.size()) {
    return false;
  }
  const std::size_t node_count = batch.nodes.size();
  if (density.size() != point_count * node_count) {
    return false;
  }
  out.assign(point_count, 0.0);
  for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
    double total = 0.0;
    const std::size_t offset = point_idx * node_count;
    for (std::size_t node_idx = 0; node_idx < node_count; ++node_idx) {
      const double w = batch.weights[node_idx];
      if (!std::isfinite(w) || w <= 0.0) {
        continue;
      }
      const double term = safe_density(density[offset + node_idx]);
      if (!std::isfinite(term) || term <= 0.0) {
        continue;
      }
      total += w * term;
    }
    out[point_idx] =
        (std::isfinite(total) && total > 0.0) ? clamp_probability(total) : 0.0;
  }
  return true;
}

struct CouplingAdaptiveSegmentInterval {
  double lo{0.0};
  double hi{0.0};
  int depth{0};
};

inline int coupling_adaptive_base_segments(double upper) {
  if (!std::isfinite(upper) || upper <= 0.0) {
    return 4;
  }
  if (upper <= 0.75) {
    return 1;
  }
  if (upper <= 3.0) {
    return 2;
  }
  if (upper <= 8.0) {
    return 3;
  }
  return 4;
}

inline int coupling_adaptive_max_segments(int base_segments, int depth_cap) {
  const int shift = std::max(0, std::min(2, depth_cap));
  const int scaled = base_segments * (1 << shift);
  return std::max(base_segments, std::min(24, scaled));
}

template <typename EvalBatchTermsFn>
double integrate_coupling_mass_batch_adaptive(
    double upper, double rel_tol, double abs_tol, int max_depth,
    EvalBatchTermsFn &&eval_batch_terms) {
  if (!(upper > 0.0)) {
    return 0.0;
  }

  std::vector<double> terms;
  const bool finite_upper = std::isfinite(upper);
  const int depth_cap = std::max(1, std::min(2, max_depth));
  const int base_segments =
      finite_upper ? std::max(1, coupling_adaptive_base_segments(upper)) : 4;
  const int max_segments =
      coupling_adaptive_max_segments(base_segments, depth_cap);

  std::vector<CouplingAdaptiveSegmentInterval> frontier;
  frontier.reserve(static_cast<std::size_t>(base_segments));
  if (finite_upper) {
    for (int seg = 0; seg < base_segments; ++seg) {
      const double lo = upper * static_cast<double>(seg) /
                        static_cast<double>(base_segments);
      const double hi = upper * static_cast<double>(seg + 1) /
                        static_cast<double>(base_segments);
      if (hi > lo) {
        frontier.push_back(CouplingAdaptiveSegmentInterval{lo, hi, 0});
      }
    }
  } else {
    constexpr double kXUpperCap = 1.0 - 1e-12;
    constexpr std::array<double, 5> kInitialXEdges{
        0.0, 0.80, 0.94, 0.985, kXUpperCap};
    for (std::size_t i = 0; i + 1 < kInitialXEdges.size(); ++i) {
      const double lo = kInitialXEdges[i];
      const double hi = kInitialXEdges[i + 1];
      if (hi > lo) {
        frontier.push_back(CouplingAdaptiveSegmentInterval{lo, hi, 0});
      }
    }
  }
  if (frontier.empty()) {
    return 0.0;
  }

  double accepted_total = 0.0;
  std::vector<CouplingAdaptiveSegmentInterval> next_frontier;
  while (!frontier.empty()) {
    next_frontier.clear();
    next_frontier.reserve(frontier.size() * 2u);
    const int current_leaf_count = static_cast<int>(frontier.size());
    int planned_splits = 0;

    for (const CouplingAdaptiveSegmentInterval &seg : frontier) {
      uuber::TimeBatch coarse_batch;
      const bool built_batch =
          finite_upper
              ? [&]() {
                  coarse_batch = uuber::build_time_batch(seg.lo, seg.hi);
                  return true;
                }()
              : build_coupling_time_batch_mass_infinite_segment(seg.lo, seg.hi,
                                                                coarse_batch);
      if (!built_batch) {
        continue;
      }
      if (coarse_batch.nodes.empty() ||
          coarse_batch.nodes.size() != coarse_batch.weights.size()) {
        continue;
      }
      if (!eval_batch_terms(coarse_batch.nodes, terms)) {
        return 0.0;
      }
      const double coarse_val = integrate_node_density_batch(coarse_batch, terms);
      double min_term = std::numeric_limits<double>::infinity();
      double max_term = 0.0;
      for (double term : terms) {
        const double safe = safe_density(term);
        if (safe <= 0.0) {
          continue;
        }
        min_term = std::min(min_term, safe);
        max_term = std::max(max_term, safe);
      }
      if (!std::isfinite(min_term)) {
        min_term = 0.0;
      }
      const double span = std::max(0.0, max_term - min_term);
      double weight_sum = 0.0;
      for (double w : coarse_batch.weights) {
        if (std::isfinite(w) && w > 0.0) {
          weight_sum += w;
        }
      }
      const double err = 0.2 * span * std::max(0.0, weight_sum);
      const double width_scale =
          finite_upper
              ? ((upper > 0.0) ? (std::max(0.0, seg.hi - seg.lo) / upper) : 1.0)
              : std::max(1e-6, std::max(0.0, seg.hi - seg.lo));
      const double tol = std::max(
          1e-10,
          8.0 * (std::max(0.0, abs_tol) * width_scale +
                 std::max(0.0, rel_tol) * std::fabs(coarse_val)));
      const double mid = 0.5 * (seg.lo + seg.hi);
      const bool can_split =
          err > tol && seg.depth < depth_cap && mid > seg.lo && seg.hi > mid &&
          (current_leaf_count + planned_splits + 1 <= max_segments);
      if (can_split) {
        record_unified_outcome_adaptive_segment_split();
        next_frontier.push_back(
            CouplingAdaptiveSegmentInterval{seg.lo, mid, seg.depth + 1});
        next_frontier.push_back(
            CouplingAdaptiveSegmentInterval{mid, seg.hi, seg.depth + 1});
        ++planned_splits;
      } else {
        record_unified_outcome_adaptive_segment_accept();
        accepted_total += coarse_val;
      }
    }

    frontier.swap(next_frontier);
  }

  if (!std::isfinite(accepted_total) || accepted_total <= 0.0) {
    return 0.0;
  }
  return clamp_probability(accepted_total);
}

template <typename EvalBatchTermsFn>
bool integrate_coupling_mass_batch_adaptive(
    double upper, double rel_tol, double abs_tol, int max_depth,
    std::size_t point_count, EvalBatchTermsFn &&eval_batch_terms,
    std::vector<double> &out) {
  out.assign(point_count, 0.0);
  if (point_count == 0u) {
    return true;
  }
  if (!(upper > 0.0)) {
    return true;
  }

  std::vector<double> terms;
  std::vector<double> segment_values;
  const bool finite_upper = std::isfinite(upper);
  const int depth_cap = std::max(1, std::min(2, max_depth));
  const int base_segments =
      finite_upper ? std::max(1, coupling_adaptive_base_segments(upper)) : 4;
  const int max_segments =
      coupling_adaptive_max_segments(base_segments, depth_cap);

  std::vector<CouplingAdaptiveSegmentInterval> frontier;
  frontier.reserve(static_cast<std::size_t>(base_segments));
  if (finite_upper) {
    for (int seg = 0; seg < base_segments; ++seg) {
      const double lo = upper * static_cast<double>(seg) /
                        static_cast<double>(base_segments);
      const double hi = upper * static_cast<double>(seg + 1) /
                        static_cast<double>(base_segments);
      if (hi > lo) {
        frontier.push_back(CouplingAdaptiveSegmentInterval{lo, hi, 0});
      }
    }
  } else {
    constexpr double kXUpperCap = 1.0 - 1e-12;
    constexpr std::array<double, 5> kInitialXEdges{
        0.0, 0.80, 0.94, 0.985, kXUpperCap};
    for (std::size_t i = 0; i + 1 < kInitialXEdges.size(); ++i) {
      const double lo = kInitialXEdges[i];
      const double hi = kInitialXEdges[i + 1];
      if (hi > lo) {
        frontier.push_back(CouplingAdaptiveSegmentInterval{lo, hi, 0});
      }
    }
  }
  if (frontier.empty()) {
    return true;
  }

  std::vector<CouplingAdaptiveSegmentInterval> next_frontier;
  while (!frontier.empty()) {
    next_frontier.clear();
    next_frontier.reserve(frontier.size() * 2u);
    const int current_leaf_count = static_cast<int>(frontier.size());
    int planned_splits = 0;

    for (const CouplingAdaptiveSegmentInterval &seg : frontier) {
      uuber::TimeBatch coarse_batch;
      const bool built_batch =
          finite_upper
              ? [&]() {
                  coarse_batch = uuber::build_time_batch(seg.lo, seg.hi);
                  return true;
                }()
              : build_coupling_time_batch_mass_infinite_segment(seg.lo, seg.hi,
                                                                coarse_batch);
      if (!built_batch) {
        continue;
      }
      if (coarse_batch.nodes.empty() ||
          coarse_batch.nodes.size() != coarse_batch.weights.size()) {
        continue;
      }
      if (!eval_batch_terms(coarse_batch.nodes, terms) ||
          !integrate_node_density_batch_per_point(coarse_batch, terms,
                                                 point_count,
                                                 segment_values)) {
        return false;
      }
      const std::size_t node_count = coarse_batch.nodes.size();
      double weight_sum = 0.0;
      for (double w : coarse_batch.weights) {
        if (std::isfinite(w) && w > 0.0) {
          weight_sum += w;
        }
      }
      const double width_scale =
          finite_upper
              ? ((upper > 0.0) ? (std::max(0.0, seg.hi - seg.lo) / upper) : 1.0)
              : std::max(1e-6, std::max(0.0, seg.hi - seg.lo));
      const double mid = 0.5 * (seg.lo + seg.hi);
      bool should_split = false;
      for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
        double min_term = std::numeric_limits<double>::infinity();
        double max_term = 0.0;
        const std::size_t offset = point_idx * node_count;
        for (std::size_t node_idx = 0; node_idx < node_count; ++node_idx) {
          const double safe = safe_density(terms[offset + node_idx]);
          if (safe <= 0.0) {
            continue;
          }
          min_term = std::min(min_term, safe);
          max_term = std::max(max_term, safe);
        }
        if (!std::isfinite(min_term)) {
          min_term = 0.0;
        }
        const double span = std::max(0.0, max_term - min_term);
        const double err = 0.2 * span * std::max(0.0, weight_sum);
        const double tol = std::max(
            1e-10,
            8.0 * (std::max(0.0, abs_tol) * width_scale +
                   std::max(0.0, rel_tol) *
                       std::fabs(segment_values[point_idx])));
        if (err > tol) {
          should_split = true;
          break;
        }
      }

      const bool can_split =
          should_split && seg.depth < depth_cap && mid > seg.lo &&
          seg.hi > mid &&
          (current_leaf_count + planned_splits + 1 <= max_segments);
      if (can_split) {
        record_unified_outcome_adaptive_segment_split();
        next_frontier.push_back(
            CouplingAdaptiveSegmentInterval{seg.lo, mid, seg.depth + 1});
        next_frontier.push_back(
            CouplingAdaptiveSegmentInterval{mid, seg.hi, seg.depth + 1});
        ++planned_splits;
      } else {
        record_unified_outcome_adaptive_segment_accept();
        for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
          out[point_idx] += segment_values[point_idx];
        }
      }
    }

    frontier.swap(next_frontier);
  }

  for (double &value : out) {
    value = (std::isfinite(value) && value > 0.0) ? clamp_probability(value)
                                                  : 0.0;
  }
  return true;
}

} // namespace

namespace {

struct CouplingProgramInterpreter {
  const uuber::NativeContext &ctx;
  const OutcomeCouplingProgram &program;
  double upper;
  int component_idx;
  const TrialParamSet *trial_params;
  const std::string &trial_type_key;
  double rel_tol;
  double abs_tol;
  int max_depth;
  bool include_na_donors;
  int outcome_idx_context;
  const uuber::BitsetState *forced_complete_bits;
  bool forced_complete_bits_valid;
  const uuber::BitsetState *forced_survive_bits;
  bool forced_survive_bits_valid;
  const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch;
  std::vector<GenericCouplingProviderRuntime> generic_runtimes;
  CouplingBatchScratch scratch;

  explicit CouplingProgramInterpreter(
      const uuber::NativeContext &ctx_, const OutcomeCouplingProgram &program_,
      double upper_, int component_idx_, const TrialParamSet *trial_params_,
      const std::string &trial_type_key_, double rel_tol_, double abs_tol_,
      int max_depth_, bool include_na_donors_, int outcome_idx_context_,
      const uuber::BitsetState *forced_complete_bits_,
      bool forced_complete_bits_valid_,
      const uuber::BitsetState *forced_survive_bits_,
      bool forced_survive_bits_valid_,
      const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch_)
      : ctx(ctx_), program(program_), upper(upper_), component_idx(component_idx_),
        trial_params(trial_params_), trial_type_key(trial_type_key_),
        rel_tol(rel_tol_), abs_tol(abs_tol_), max_depth(max_depth_),
        include_na_donors(include_na_donors_),
        outcome_idx_context(outcome_idx_context_),
        forced_complete_bits(forced_complete_bits_),
        forced_complete_bits_valid(forced_complete_bits_valid_),
        forced_survive_bits(forced_survive_bits_),
        forced_survive_bits_valid(forced_survive_bits_valid_),
        trial_params_soa_batch(trial_params_soa_batch_) {
    generic_runtimes.reserve(program.generic_payloads.size());
    for (const auto &payload : program.generic_payloads) {
      generic_runtimes.push_back(build_generic_coupling_provider_runtime(
          ctx, payload, forced_complete_bits, forced_complete_bits_valid,
          forced_survive_bits, forced_survive_bits_valid));
    }
  }

  std::size_t point_count() const noexcept { return trial_params_soa_batch.size(); }

  bool valid_program() const noexcept {
    return program.valid &&
           program.domain == uuber::VectorProgramDomain::OutcomeCoupling &&
           program.outputs.result_slot >= 0 &&
           program.outputs.result_slot <
               static_cast<int>(program.ops.size());
  }

  bool evaluate_result(std::vector<double> &out) {
    out.assign(point_count(), 0.0);
    if (point_count() == 0u) {
      return true;
    }
    std::vector<double> upper_nodes{upper};
    std::vector<double> raw;
    if (!evaluate_slot_terms(program.outputs.result_slot, upper_nodes, raw) ||
        raw.size() != point_count()) {
      return false;
    }
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] = (std::isfinite(raw[i]) && raw[i] > 0.0)
                   ? clamp_probability(raw[i])
                   : 0.0;
    }
    return true;
  }

private:
  const uuber::VectorOp *op_at(int slot) const noexcept {
    if (slot < 0 || slot >= static_cast<int>(program.ops.size())) {
      return nullptr;
    }
    return &program.ops[static_cast<std::size_t>(slot)];
  }

  std::size_t flattened_size(std::size_t node_count) const noexcept {
    return point_count() * node_count;
  }

  bool child_slots(const uuber::VectorOp &op,
                   std::vector<int> &child_slots_out) const {
    child_slots_out.clear();
    if (op.child_count <= 0) {
      return true;
    }
    if (op.child_begin < 0 ||
        op.child_begin + op.child_count >
            static_cast<int>(program.children.size())) {
      return false;
    }
    child_slots_out.insert(
        child_slots_out.end(),
        program.children.begin() + static_cast<std::ptrdiff_t>(op.child_begin),
        program.children.begin() +
            static_cast<std::ptrdiff_t>(op.child_begin + op.child_count));
    return true;
  }

  bool expand_nodes(const std::vector<double> &nodes) {
    expand_time_batch_by_trial_params(nodes, trial_params_soa_batch,
                                      scratch.expanded_times,
                                      scratch.expanded_trial_params_soa);
    return scratch.expanded_times.size() == flattened_size(nodes.size()) &&
           scratch.expanded_trial_params_soa.size() == flattened_size(nodes.size());
  }

  bool evaluate_event_terms(const CouplingEventPayload &payload,
                            uuber::VectorOpCode code,
                            const std::vector<double> &nodes,
                            std::vector<double> &out) {
    if (!expand_nodes(nodes)) {
      return false;
    }
    const bool need_density = code == uuber::VectorOpCode::EventDensity;
    const bool need_survival = code == uuber::VectorOpCode::EventSurvival;
    const bool need_cdf = code == uuber::VectorOpCode::EventCDF;
    std::vector<double> *density_out = need_density ? &out : nullptr;
    std::vector<double> *survival_out = need_survival ? &out : nullptr;
    std::vector<double> *cdf_out = need_cdf ? &out : nullptr;
    return evaluate_labelref_batch(
               ctx, payload, component_idx, trial_params, trial_type_key,
               include_na_donors, outcome_idx_context, scratch.expanded_times,
               need_density, need_survival, need_cdf, density_out, survival_out,
               cdf_out, &scratch, &scratch.expanded_trial_params_soa) &&
           out.size() == flattened_size(nodes.size());
  }

  bool evaluate_generic_direct_cdf_terms(std::size_t payload_idx,
                                         const std::vector<double> &nodes,
                                         std::vector<double> &out) {
    if (payload_idx >= program.generic_payloads.size() ||
        payload_idx >= generic_runtimes.size() || !expand_nodes(nodes)) {
      return false;
    }
    const CouplingGenericPayload &generic = program.generic_payloads[payload_idx];
    GenericCouplingProviderRuntime &runtime = generic_runtimes[payload_idx];
    record_generic_coupling_provider_usage(runtime);
    if (generic_coupling_runtime_has_labelref_fastpath(runtime)) {
      return evaluate_labelref_batch(
                 ctx, runtime.target_ref, component_idx, trial_params,
                 trial_type_key, include_na_donors, outcome_idx_context,
                 scratch.expanded_times, false, false, true, nullptr, nullptr,
                 &out, &runtime.labelref_scratch,
                 &scratch.expanded_trial_params_soa) &&
             out.size() == flattened_size(nodes.size());
    }

    NodeEvalState state(ctx, 0.0, component_idx, trial_params, trial_type_key,
                        false, -1, nullptr, nullptr,
                        forced_complete_bits_valid ? forced_complete_bits
                                                  : nullptr,
                        forced_complete_bits_valid,
                        forced_survive_bits_valid ? forced_survive_bits
                                                  : nullptr,
                        forced_survive_bits_valid);
    state.trial_params_soa_batch = &scratch.expanded_trial_params_soa;
    uuber::KernelNodeBatchValues batch_values;
    const int node_idx = resolve_dense_node_idx_required(ctx, generic.node_id);
    if (!evaluator_eval_node_batch_with_state_dense(
            node_idx, scratch.expanded_times, state, EvalNeed::kCDF,
            runtime.noderef_kernel_batch_runtime, batch_values) ||
        batch_values.cdf.size() != scratch.expanded_times.size()) {
      return false;
    }
    out.resize(batch_values.cdf.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] = clamp_probability(batch_values.cdf[i]);
    }
    return true;
  }

  bool evaluate_generic_integrand_terms(std::size_t payload_idx,
                                        const std::vector<double> &nodes,
                                        std::vector<double> &out) {
    if (payload_idx >= program.generic_payloads.size() ||
        payload_idx >= generic_runtimes.size()) {
      return false;
    }
    const CouplingGenericPayload &generic = program.generic_payloads[payload_idx];
    GenericCouplingProviderRuntime &runtime = generic_runtimes[payload_idx];
    record_generic_coupling_provider_usage(runtime);
    return evaluate_generic_terms_resolved_provider_batch(
               ctx, generic, runtime, component_idx, trial_params,
               trial_type_key, include_na_donors, outcome_idx_context,
               forced_complete_bits, forced_complete_bits_valid,
               forced_survive_bits, forced_survive_bits_valid, nodes,
               trial_params_soa_batch, out) &&
           out.size() == flattened_size(nodes.size());
  }

  bool evaluate_integral_terms(const uuber::VectorOp &op,
                               const std::vector<double> &nodes,
                               std::vector<double> &out) {
    std::vector<int> children;
    if (!child_slots(op, children) || children.size() != 1u) {
      return false;
    }
    std::vector<double> integral_values;
    if (!integrate_coupling_mass_batch_adaptive(
            upper, rel_tol, abs_tol, max_depth, point_count(),
            [&](const std::vector<double> &quad_nodes,
                std::vector<double> &terms) -> bool {
              return evaluate_slot_terms(children[0], quad_nodes, terms);
            },
            integral_values) ||
        integral_values.size() != point_count()) {
      return false;
    }
    out.assign(flattened_size(nodes.size()), 0.0);
    for (std::size_t point_idx = 0; point_idx < point_count(); ++point_idx) {
      const double value = integral_values[point_idx];
      const std::size_t offset = point_idx * nodes.size();
      for (std::size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {
        out[offset + node_idx] = value;
      }
    }
    return true;
  }

public:
  bool evaluate_slot_terms(int slot, const std::vector<double> &nodes,
                           std::vector<double> &out) {
    const uuber::VectorOp *op = op_at(slot);
    if (op == nullptr) {
      return false;
    }
    switch (op->code) {
    case uuber::VectorOpCode::EventDensity:
    case uuber::VectorOpCode::EventCDF:
    case uuber::VectorOpCode::EventSurvival:
      if (op->payload_idx < 0 ||
          op->payload_idx >= static_cast<int>(program.event_payloads.size())) {
        return false;
      }
      return evaluate_event_terms(
          program.event_payloads[static_cast<std::size_t>(op->payload_idx)],
          op->code, nodes, out);
    case uuber::VectorOpCode::GenericDirectCDF:
      if (op->payload_idx < 0) {
        return false;
      }
      return evaluate_generic_direct_cdf_terms(
          static_cast<std::size_t>(op->payload_idx), nodes, out);
    case uuber::VectorOpCode::GenericIntegrand:
      if (op->payload_idx < 0) {
        return false;
      }
      return evaluate_generic_integrand_terms(
          static_cast<std::size_t>(op->payload_idx), nodes, out);
    case uuber::VectorOpCode::Complement: {
      std::vector<int> children;
      if (!child_slots(*op, children) || children.size() != 1u ||
          !evaluate_slot_terms(children[0], nodes, out)) {
        return false;
      }
      for (double &value : out) {
        value = clamp_probability(1.0 - clamp_probability(value));
      }
      return true;
    }
    case uuber::VectorOpCode::Multiply: {
      std::vector<int> children;
      if (!child_slots(*op, children) || children.empty() ||
          !evaluate_slot_terms(children[0], nodes, out)) {
        return false;
      }
      std::vector<double> child_values;
      for (std::size_t i = 1; i < children.size(); ++i) {
        if (!evaluate_slot_terms(children[i], nodes, child_values) ||
            child_values.size() != out.size()) {
          return false;
        }
        for (std::size_t j = 0; j < out.size(); ++j) {
          const double lhs = out[j];
          const double rhs = child_values[j];
          const double product =
              (std::isfinite(lhs) && std::isfinite(rhs)) ? lhs * rhs : 0.0;
          out[j] = (std::isfinite(product) && product > 0.0) ? product : 0.0;
        }
      }
      return true;
    }
    case uuber::VectorOpCode::Add: {
      std::vector<int> children;
      if (!child_slots(*op, children) || children.empty()) {
        return false;
      }
      out.assign(flattened_size(nodes.size()), 0.0);
      std::vector<double> child_values;
      for (int child_slot : children) {
        if (!evaluate_slot_terms(child_slot, nodes, child_values) ||
            child_values.size() != out.size()) {
          return false;
        }
        for (std::size_t j = 0; j < out.size(); ++j) {
          const double lhs = out[j];
          const double rhs = child_values[j];
          out[j] = (std::isfinite(lhs) ? lhs : 0.0) +
                   (std::isfinite(rhs) ? rhs : 0.0);
        }
      }
      return true;
    }
    case uuber::VectorOpCode::Subtract: {
      std::vector<int> children;
      if (!child_slots(*op, children) || children.empty() ||
          !evaluate_slot_terms(children[0], nodes, out)) {
        return false;
      }
      std::vector<double> child_values;
      for (std::size_t i = 1; i < children.size(); ++i) {
        if (!evaluate_slot_terms(children[i], nodes, child_values) ||
            child_values.size() != out.size()) {
          return false;
        }
        for (std::size_t j = 0; j < out.size(); ++j) {
          const double lhs = std::isfinite(out[j]) ? out[j] : 0.0;
          const double rhs = std::isfinite(child_values[j]) ? child_values[j] : 0.0;
          const double value = lhs - rhs;
          out[j] = std::isfinite(value) ? value : 0.0;
        }
      }
      return true;
    }
    case uuber::VectorOpCode::Integral:
      return evaluate_integral_terms(*op, nodes, out);
    }
    return false;
  }
};

inline bool evaluate_outcome_coupling_vector_program_batch_internal(
    const uuber::NativeContext &ctx, const OutcomeCouplingProgram &program,
    double upper, int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, double rel_tol, double abs_tol,
    int max_depth, bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    std::vector<double> &out) {
  out.assign(trial_params_soa_batch.size(), 0.0);
  if (trial_params_soa_batch.empty()) {
    return true;
  }
  CouplingProgramInterpreter interpreter(
      ctx, program, upper, component_idx, trial_params, trial_type_key, rel_tol,
      abs_tol, max_depth, include_na_donors, outcome_idx_context,
      forced_complete_bits, forced_complete_bits_valid, forced_survive_bits,
      forced_survive_bits_valid, trial_params_soa_batch);
  if (!interpreter.valid_program()) {
    return true;
  }
  return interpreter.evaluate_result(out);
}

} // namespace

double evaluate_outcome_coupling_unified(
    const uuber::NativeContext &ctx, const OutcomeCouplingProgram &program,
    double upper, int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, double rel_tol, double abs_tol,
    int max_depth, bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid) {
  record_unified_outcome_coupling_eval_call();
  std::vector<const uuber::TrialParamsSoA *> singleton_trial_params_soa{
      resolve_trial_params_soa(ctx, trial_params)};
  std::vector<double> out;
  if (!evaluate_outcome_coupling_vector_program_batch_internal(
          ctx, program, upper, component_idx, trial_params, trial_type_key,
          rel_tol, abs_tol, max_depth, include_na_donors, outcome_idx_context,
          forced_complete_bits, forced_complete_bits_valid, forced_survive_bits,
          forced_survive_bits_valid, singleton_trial_params_soa, out) ||
      out.empty()) {
    return 0.0;
  }
  return out[0];
}

bool evaluate_outcome_coupling_unified_batch(
    const uuber::NativeContext &ctx, const OutcomeCouplingProgram &program,
    double upper, int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, double rel_tol, double abs_tol,
    int max_depth, bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    std::vector<double> &out) {
  record_unified_outcome_coupling_eval_call();
  return evaluate_outcome_coupling_vector_program_batch_internal(
      ctx, program, upper, component_idx, trial_params, trial_type_key,
      rel_tol, abs_tol, max_depth, include_na_donors, outcome_idx_context,
      forced_complete_bits, forced_complete_bits_valid, forced_survive_bits,
      forced_survive_bits_valid, trial_params_soa_batch, out);
}
