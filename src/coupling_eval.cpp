#include "coupling_eval.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "evaluator_internal.h"
#include "runtime_stats.h"

using uuber::AccDistParams;
using uuber::LabelRef;

namespace {

NodeEvalResult eval_event_unforced(
    const uuber::NativeContext &ctx, const LabelRef &label_ref, double t,
    int component_idx, EvalNeed need, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const uuber::TrialParamsSoA *trial_params_soa = nullptr) {
  if (!trial_params_soa) {
    trial_params_soa = resolve_trial_params_soa(ctx, trial_params);
  }
  return evaluator_eval_event_ref_idx(
      ctx, label_ref, 0u, t, component_idx, need, trial_params, trial_type_key,
      include_na_donors, outcome_idx_context, nullptr, nullptr, nullptr,
      nullptr, nullptr, trial_params_soa, nullptr);
}

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
  std::vector<double> temp_survival;
};

inline bool build_coupling_acc_ref(
    const uuber::NativeContext &ctx, const LabelRef &ref, int component_idx,
    const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *trial_params_soa, CouplingAccRef &out) {
  if (ref.acc_idx < 0 ||
      ref.acc_idx >= static_cast<int>(ctx.accumulators.size())) {
    return false;
  }
  const uuber::NativeAccumulator &acc =
      ctx.accumulators[static_cast<std::size_t>(ref.acc_idx)];
  const TrialAccumulatorParams *override =
      evaluator_get_trial_param_entry(trial_params, ref.acc_idx);
  if (!evaluator_component_active_idx(acc, component_idx, override)) {
    return false;
  }
  const int onset_kind = override ? override->onset_kind : acc.onset_kind;
  if (ctx.has_chained_onsets && onset_kind != uuber::ONSET_ABSOLUTE) {
    return false;
  }
  evaluator_resolve_event_numeric_params(
      acc, ref.acc_idx, override, trial_params_soa, out.onset, out.q, out.cfg);
  out.lower_bound = default_lower_bound_transform(out.cfg);
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

inline bool evaluate_labelref_batch(
    const uuber::NativeContext &ctx, const LabelRef &ref, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    const std::vector<double> &times, bool need_density, bool need_survival,
    bool need_cdf, std::vector<double> *density_out,
    std::vector<double> *survival_out, std::vector<double> *cdf_out,
    CouplingBatchScratch *scratch = nullptr) {
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

  const uuber::TrialParamsSoA *trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);
  const bool can_specialize = coupling_acc_specialization_allowed(
      ctx, include_na_donors, outcome_idx_context);
  CouplingAccRef ref_acc;
  if (can_specialize &&
      build_coupling_acc_ref(ctx, ref, component_idx, trial_params,
                             trial_params_soa, ref_acc)) {
    return evaluate_coupling_acc_batch(ref_acc, times, need_density,
                                       need_survival, need_cdf, density_out,
                                       survival_out, cdf_out, scratch);
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
  for (std::size_t i = 0; i < n; ++i) {
    NodeEvalResult ev = eval_event_unforced(
        ctx, ref, times[i], component_idx, need, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, trial_params_soa);
    if (need_density) {
      (*density_out)[i] = safe_density(ev.density);
    }
    if (need_survival) {
      (*survival_out)[i] = clamp_probability(ev.survival);
    }
    if (need_cdf) {
      (*cdf_out)[i] = clamp_probability(ev.cdf);
    }
  }
  return true;
}

struct GenericCouplingProviderRuntime {
  LabelRef target_ref{};
  std::vector<LabelRef> competitor_refs;
  CouplingBatchScratch labelref_scratch;
  uuber::KernelBatchRuntimeState noderef_kernel_batch_runtime;
};

inline bool generic_coupling_runtime_has_labelref_fastpath(
    const GenericCouplingProviderRuntime &runtime) {
  const LabelRef &ref = runtime.target_ref;
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
    const OutcomeCouplingGenericNodeIntegralPayload &generic,
    LabelRef &target_ref, std::vector<LabelRef> &competitor_refs) {
  const int node_idx = resolve_dense_node_idx(ctx, generic.node_id);
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  const uuber::IrNode &target_node =
      ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  if (target_node.event_idx < 0) {
    return false;
  }
  target_ref = node_label_ref(ctx, target_node);
  if (!(target_ref.acc_idx >= 0 || target_ref.pool_idx >= 0 ||
        target_ref.outcome_idx >= 0)) {
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
    LabelRef comp_ref = node_label_ref(ctx, comp_node);
    if (!(comp_ref.acc_idx >= 0 || comp_ref.pool_idx >= 0 ||
          comp_ref.outcome_idx >= 0)) {
      competitor_refs.clear();
      return false;
    }
    competitor_refs.push_back(comp_ref);
  }
  return true;
}

inline double evaluate_generic_direct_cdf_resolved_provider(
    const uuber::NativeContext &ctx,
    const OutcomeCouplingGenericNodeIntegralPayload &generic,
    GenericCouplingProviderRuntime &runtime, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid, double upper) {
  std::vector<double> upper_vec{upper};
  if (generic_coupling_runtime_has_labelref_fastpath(runtime)) {
    if (!evaluate_labelref_batch(
            ctx, runtime.target_ref, component_idx, trial_params,
            trial_type_key, include_na_donors, outcome_idx_context, upper_vec,
            false, false, true, nullptr, nullptr,
            &runtime.labelref_scratch.cdf_values, &runtime.labelref_scratch) ||
        runtime.labelref_scratch.cdf_values.empty()) {
      return 0.0;
    }
    return clamp_probability(runtime.labelref_scratch.cdf_values[0]);
  }

  const ForcedStateView forced_state = make_forced_state_view(
      nullptr, forced_complete_bits, &forced_complete_bits_valid,
      forced_survive_bits, &forced_survive_bits_valid,
      ctx.ir.label_id_to_bit_idx.empty() ? nullptr
                                         : &ctx.ir.label_id_to_bit_idx);
  uuber::KernelNodeBatchValues batch_values;
  if (!evaluator_eval_node_with_forced_state_view_batch(
          ctx, generic.node_id, upper_vec, component_idx, EvalNeed::kCDF,
          trial_params, trial_type_key, nullptr, forced_state,
          batch_values) ||
      batch_values.cdf.empty()) {
    return 0.0;
  }
  return clamp_probability(batch_values.cdf[0]);
}

inline bool evaluate_generic_terms_resolved_provider(
    const uuber::NativeContext &ctx,
    const OutcomeCouplingGenericNodeIntegralPayload &generic,
    GenericCouplingProviderRuntime &runtime, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid, const std::vector<double> &nodes,
    std::vector<double> &terms) {
  if (!generic_coupling_runtime_has_labelref_fastpath(runtime)) {
    uuber::KernelBatchRuntimeState *kernel_batch_runtime =
        ctx.kernel_program.valid ? &runtime.noderef_kernel_batch_runtime
                                 : nullptr;
    return evaluator_node_density_with_competitors_batch_internal(
        ctx, generic.node_id, nodes, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, generic.competitor_node_ids, trial_params,
        trial_type_key, include_na_donors, outcome_idx_context, terms, nullptr,
        kernel_batch_runtime);
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
      for (const LabelRef &comp_ref : runtime.competitor_refs) {
        CouplingAccRef comp_acc;
        if (!build_coupling_acc_ref(ctx, comp_ref, component_idx, trial_params,
                                    trial_params_soa, comp_acc)) {
          all_competitors_specialized = false;
          break;
        }
        competitor_accs.push_back(comp_acc);
      }
      if (all_competitors_specialized) {
        if (!evaluate_coupling_acc_batch(
                target_acc, nodes, true, false, false,
                &runtime.labelref_scratch.pdf_values, nullptr, nullptr,
                &runtime.labelref_scratch)) {
          return false;
        }
        terms.assign(nodes.size(), 1.0);
        for (const CouplingAccRef &comp_acc : competitor_accs) {
          if (!evaluate_coupling_acc_batch(
                  comp_acc, nodes, false, true, false, nullptr,
                  &runtime.labelref_scratch.temp_survival, nullptr,
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
          include_na_donors, outcome_idx_context, nodes, true, false, false,
          &runtime.labelref_scratch.pdf_values, nullptr, nullptr,
          &runtime.labelref_scratch)) {
    return false;
  }
  terms.assign(nodes.size(), 1.0);
  for (const LabelRef &comp_ref : runtime.competitor_refs) {
    if (!evaluate_labelref_batch(
            ctx, comp_ref, component_idx, trial_params, trial_type_key,
            include_na_donors, outcome_idx_context, nodes, false, true, false,
            nullptr, &runtime.labelref_scratch.temp_survival, nullptr,
            &runtime.labelref_scratch) ||
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
    const OutcomeCouplingGenericNodeIntegralPayload &generic,
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

inline double integrate_node_density_batch(
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
  if (!program.valid || program.kind == OutcomeCouplingOpKind::None) {
    return 0.0;
  }
  if (program.kind == OutcomeCouplingOpKind::GenericNodeIntegral) {
    const OutcomeCouplingGenericNodeIntegralPayload &generic = program.generic;
    GenericCouplingProviderRuntime provider_runtime =
        build_generic_coupling_provider_runtime(
            ctx, generic, forced_complete_bits, forced_complete_bits_valid,
            forced_survive_bits, forced_survive_bits_valid);
    if (generic.competitor_node_ids.empty()) {
      record_generic_coupling_provider_usage(provider_runtime);
      return evaluate_generic_direct_cdf_resolved_provider(
          ctx, generic, provider_runtime, component_idx, trial_params,
          trial_type_key, include_na_donors, outcome_idx_context,
          forced_complete_bits, forced_complete_bits_valid, forced_survive_bits,
          forced_survive_bits_valid, upper);
    }
    record_generic_coupling_provider_usage(provider_runtime);
    auto eval_generic_terms =
        [&](const std::vector<double> &nodes, std::vector<double> &terms)
        -> bool {
      return evaluate_generic_terms_resolved_provider(
          ctx, generic, provider_runtime, component_idx, trial_params,
          trial_type_key, include_na_donors, outcome_idx_context,
          forced_complete_bits, forced_complete_bits_valid, forced_survive_bits,
          forced_survive_bits_valid, nodes, terms);
    };
    return integrate_coupling_mass_batch_adaptive(
        upper, rel_tol, abs_tol, max_depth, eval_generic_terms);
  }
  if (program.kind == OutcomeCouplingOpKind::NWay) {
    const OutcomeCouplingNWayPayload &gate = program.nway;
    if (gate.competitor_refs.empty() || !(upper > 0.0)) {
      return 0.0;
    }
    CouplingBatchScratch eval_scratch;
    std::vector<double> gate_cdf_vec;
    std::vector<double> gate_upper_times{upper};
    if (!evaluate_labelref_batch(
            ctx, gate.gate_ref, component_idx, trial_params, trial_type_key,
            include_na_donors, outcome_idx_context, gate_upper_times, false,
            false, true, nullptr, nullptr, &gate_cdf_vec, &eval_scratch) ||
        gate_cdf_vec.empty()) {
      return 0.0;
    }
    const double gate_cdf = clamp_probability(gate_cdf_vec[0]);
    if (!std::isfinite(gate_cdf) || gate_cdf <= 0.0) {
      return 0.0;
    }

    auto eval_nway_terms =
        [&](const std::vector<double> &nodes,
            std::vector<double> &terms) -> bool {
      std::vector<double> target_density;
      if (!evaluate_labelref_batch(ctx, gate.target_ref, component_idx,
                                   trial_params, trial_type_key,
                                   include_na_donors, outcome_idx_context,
                                   nodes, true, false, false, &target_density,
                                   nullptr, nullptr, &eval_scratch)) {
        return false;
      }
      terms.assign(nodes.size(), 1.0);
      for (const LabelRef &comp_ref : gate.competitor_refs) {
        std::vector<double> comp_survival;
        if (!evaluate_labelref_batch(
                ctx, comp_ref, component_idx, trial_params, trial_type_key,
                include_na_donors, outcome_idx_context, nodes, false, true,
                false, nullptr, &comp_survival, nullptr, &eval_scratch)) {
          return false;
        }
        for (std::size_t i = 0; i < terms.size(); ++i) {
          terms[i] *= clamp_probability(comp_survival[i]);
        }
      }
      for (std::size_t i = 0; i < terms.size(); ++i) {
        const double fx = safe_density(target_density[i]);
        const double surv_prod = clamp_probability(terms[i]);
        if (!std::isfinite(fx) || fx <= 0.0 || surv_prod <= 0.0) {
          terms[i] = 0.0;
        } else {
          terms[i] = fx * surv_prod;
        }
      }
      return true;
    };

    const double order_mass = integrate_coupling_mass_batch_adaptive(
        upper, rel_tol, abs_tol, max_depth, eval_nway_terms);
    const double total = gate_cdf * order_mass;
    if (!std::isfinite(total) || total <= 0.0) {
      return 0.0;
    }
    return clamp_probability(total);
  }
  if (program.kind == OutcomeCouplingOpKind::Pair) {
    const OutcomeCouplingPairPayload &pair = program.pair;
    if (!(upper > 0.0)) {
      return 0.0;
    }

    CouplingBatchScratch eval_scratch;
    double coeff = 1.0;
    if (std::isfinite(upper)) {
      std::vector<double> c_cdf_vec;
      const std::vector<double> upper_vec{upper};
      if (!evaluate_labelref_batch(
              ctx, pair.c_ref, component_idx, trial_params, trial_type_key,
              include_na_donors, outcome_idx_context, upper_vec, false, false,
              true, nullptr, nullptr, &c_cdf_vec, &eval_scratch) ||
          c_cdf_vec.empty()) {
        return 0.0;
      }
      coeff = clamp_probability(c_cdf_vec[0]);
    }

    auto integrate_pair_term = [&](int term_kind) -> double {
      auto eval_terms =
          [&](const std::vector<double> &nodes,
              std::vector<double> &terms) -> bool {
        std::vector<double> x_density;
        std::vector<double> x_survival;
        std::vector<double> y_survival;
        std::vector<double> c_density;
        std::vector<double> c_survival;
        if (!evaluate_labelref_batch(
                ctx, pair.x_ref, component_idx, trial_params, trial_type_key,
                include_na_donors, outcome_idx_context, nodes, true, true,
                false, &x_density, &x_survival, nullptr, &eval_scratch) ||
            !evaluate_labelref_batch(
                ctx, pair.y_ref, component_idx, trial_params, trial_type_key,
                include_na_donors, outcome_idx_context, nodes, false, true,
                false, nullptr, &y_survival, nullptr, &eval_scratch) ||
            !evaluate_labelref_batch(
                ctx, pair.c_ref, component_idx, trial_params, trial_type_key,
                include_na_donors, outcome_idx_context, nodes, true, true,
                false, &c_density, &c_survival, nullptr, &eval_scratch)) {
          return false;
        }
        terms.assign(nodes.size(), 0.0);
        for (std::size_t i = 0; i < nodes.size(); ++i) {
          const double fx = safe_density(x_density[i]);
          const double fc = safe_density(c_density[i]);
          const double Fx =
              clamp_probability(1.0 - clamp_probability(x_survival[i]));
          const double Fy =
              clamp_probability(1.0 - clamp_probability(y_survival[i]));
          const double Fc =
              clamp_probability(1.0 - clamp_probability(c_survival[i]));
          double val = 0.0;
          if (term_kind == 0) {
            if (fc > 0.0 && Fx > 0.0) {
              val = fc * Fx;
            }
          } else if (term_kind == 1) {
            if (fx > 0.0 && Fc > 0.0) {
              val = fx * Fc;
            }
          } else {
            if (fx > 0.0 && Fy > 0.0) {
              val = fx * Fy;
            }
          }
          terms[i] = (std::isfinite(val) && val > 0.0) ? val : 0.0;
        }
        return true;
      };
      return integrate_coupling_mass_batch_adaptive(upper, rel_tol, abs_tol,
                                                    max_depth, eval_terms);
    };

    const double i_fc_fx = integrate_pair_term(0);
    const double i_fx_fc = integrate_pair_term(1);
    const double i_fx_fy = integrate_pair_term(2);
    const double total = i_fc_fx + i_fx_fc - coeff * i_fx_fy;
    if (!std::isfinite(total) || total <= 0.0) {
      return 0.0;
    }
    return clamp_probability(total);
  }
  return 0.0;
}
