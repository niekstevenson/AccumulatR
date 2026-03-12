#pragma once

#include <Rcpp.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "accumulator.h"
#include "context.h"
#include "kernel_executor.h"
#include "native_utils.h"

struct TrialAccumulatorParams {
  double onset{0.0};
  int onset_kind{uuber::ONSET_ABSOLUTE};
  int onset_source_acc_idx{-1};
  int onset_source_pool_idx{-1};
  double onset_lag{0.0};
  double q{0.0};
  double shared_q{std::numeric_limits<double>::quiet_NaN()};
  uuber::AccDistParams dist_cfg{};
  bool has_components{false};
  std::vector<std::string> components;
  std::vector<int> component_indices;
  std::string shared_trigger_id;
  bool has_override{false};
};

struct TrialParamSet {
  std::vector<TrialAccumulatorParams> acc_params;
  bool shared_trigger_layout_matches_context{true};
  bool has_any_override{false};
  std::uint64_t shared_trigger_source_fingerprint{0ULL};
  mutable uuber::TrialParamsSoA soa_cache;
  mutable bool soa_cache_valid{false};
  const TrialParamSet *shared_trigger_source_params{nullptr};
  bool shared_trigger_scratch_ready{false};
  const void *shared_trigger_plan_identity{nullptr};
  std::uint64_t shared_trigger_current_mask{0ULL};
  bool shared_trigger_mask_valid{false};
};

struct SharedTriggerPlan {
  std::vector<int> kernel_transition_begin;
  std::vector<int> kernel_transition_count;
  std::vector<int> kernel_invalidate_slot;
  std::vector<double> kernel_q;
  std::vector<double> mask_weights;
};

std::uint64_t compute_trial_param_fingerprint(const TrialParamSet &params);
void refresh_trial_param_fingerprint(TrialParamSet &params);
bool trial_paramsets_equivalent(const TrialParamSet &a,
                                const TrialParamSet &b);
uuber::NAMapCacheKey
build_na_map_cache_key_idx(const uuber::OutcomeContextInfo &info,
                           const TrialParamSet *params_ptr,
                           const std::vector<int> &component_indices,
                           const std::vector<double> &comp_weights,
                           int outcome_idx_context);
TrialAccumulatorParams base_params(const uuber::NativeAccumulator &base);
TrialParamSet build_base_paramset(const uuber::NativeContext &ctx);
void materialize_trial_params_soa(const uuber::NativeContext &ctx,
                                  const TrialParamSet *trial_params,
                                  uuber::TrialParamsSoA &out);
std::size_t shared_trigger_count(const SharedTriggerPlan &plan);
double shared_trigger_mask_weight(const SharedTriggerPlan &plan,
                                  std::uint64_t mask);
SharedTriggerPlan build_shared_trigger_plan(const uuber::NativeContext &ctx,
                                            const TrialParamSet *params_ptr);
void apply_trigger_state_inplace(TrialParamSet &params, int acc_idx, bool fail);
void ensure_trial_params_soa(const uuber::NativeContext &ctx,
                             TrialParamSet &params);
const uuber::TrialParamsSoA *
resolve_trial_params_soa(const uuber::NativeContext &ctx,
                         const TrialParamSet *trial_params);
void prepare_trigger_scratch(const uuber::NativeContext &ctx,
                             const TrialParamSet *base_params,
                             const SharedTriggerPlan &plan,
                             TrialParamSet &scratch);
void apply_trigger_q_soa_inplace(TrialParamSet &params, int acc_idx, bool fail);
void apply_trigger_q_soa_inplace_kernel(const uuber::NativeContext &ctx,
                                        TrialParamSet &params,
                                        const SharedTriggerPlan &plan,
                                        int trigger_bit_idx, bool fail);
std::unique_ptr<TrialParamSet>
build_trial_params_from_df(const uuber::NativeContext &ctx,
                           const Rcpp::Nullable<Rcpp::DataFrame> &rows_opt);

template <typename EvalFn>
inline double evaluate_shared_trigger_states(
    const uuber::NativeContext &ctx, const SharedTriggerPlan &plan,
    TrialParamSet &scratch, const char *too_many_triggers_msg, EvalFn &&eval_fn,
    uuber::KernelRuntimeState *kernel_runtime = nullptr,
    const std::vector<uuber::KernelRuntimeState *> *kernel_runtimes = nullptr) {
  const std::size_t trigger_count = shared_trigger_count(plan);
  if (trigger_count == 0) {
    return eval_fn();
  }
  if (trigger_count >= 63) {
    Rcpp::stop(too_many_triggers_msg);
  }
  const std::uint64_t n_states = 1ULL << trigger_count;
  double total = 0.0;
  std::uint64_t mask = 0ULL;
  for (std::uint64_t idx = 0; idx < n_states; ++idx) {
    const double w = plan.mask_weights.empty()
                         ? shared_trigger_mask_weight(plan, mask)
                         : plan.mask_weights[static_cast<std::size_t>(mask)];
    if (w > 0.0) {
      const double contrib = eval_fn();
      if (std::isfinite(contrib) && contrib > 0.0) {
        total += w * contrib;
      }
    }
    if (idx + 1 == n_states) {
      break;
    }
    const std::uint64_t next = idx + 1;
    const std::uint64_t next_gray = next ^ (next >> 1);
    const std::uint64_t diff = next_gray ^ mask;
    if (diff != 0ULL) {
      const int bit = static_cast<int>(__builtin_ctzll(diff));
      const bool fail = ((next_gray >> bit) & 1ULL) != 0ULL;
      apply_trigger_q_soa_inplace_kernel(ctx, scratch, plan, bit, fail);
      int invalidate_slot = 0;
      if (bit >= 0 && bit < static_cast<int>(plan.kernel_invalidate_slot.size())) {
        invalidate_slot =
            plan.kernel_invalidate_slot[static_cast<std::size_t>(bit)];
      }
      if (kernel_runtime) {
        uuber::invalidate_kernel_runtime_from_slot(*kernel_runtime,
                                                   invalidate_slot);
      }
      if (kernel_runtimes) {
        for (uuber::KernelRuntimeState *runtime_ptr : *kernel_runtimes) {
          if (!runtime_ptr || runtime_ptr == kernel_runtime) {
            continue;
          }
          uuber::invalidate_kernel_runtime_from_slot(*runtime_ptr,
                                                     invalidate_slot);
        }
      }
    }
    mask = next_gray;
  }
  scratch.shared_trigger_plan_identity = static_cast<const void *>(&plan);
  scratch.shared_trigger_current_mask = mask;
  scratch.shared_trigger_mask_valid = true;
  return total;
}
