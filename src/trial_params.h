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
  std::uint64_t value_fingerprint{0ULL};
  bool value_fingerprint_valid{false};
  mutable uuber::TrialParamsSoA soa_cache;
  mutable bool soa_cache_valid{false};
  const TrialParamSet *shared_trigger_source_params{nullptr};
  bool shared_trigger_scratch_ready{false};
  const void *shared_trigger_plan_identity{nullptr};
  std::uint64_t shared_trigger_current_mask{0ULL};
  bool shared_trigger_mask_valid{false};
};

struct SharedTriggerPlan {
  std::vector<int> trigger_acc_begin;
  std::vector<int> trigger_acc_count;
  std::vector<int> trigger_acc_indices;
  std::vector<double> trigger_q;
  std::vector<double> mask_weights;
};

struct SharedTriggerMaskSoABatch {
  std::vector<uuber::TrialParamsSoA> mask_params;
  std::vector<const uuber::TrialParamsSoA *> mask_param_ptrs;
  std::vector<double> mask_weights;
};

std::uint64_t compute_trial_param_fingerprint(const TrialParamSet &params);
void refresh_trial_param_fingerprint(TrialParamSet &params);
bool trial_paramsets_equivalent(const TrialParamSet &a,
                                const TrialParamSet &b);
std::uint64_t compute_trial_param_value_fingerprint(const TrialParamSet &params);
bool trial_paramsets_value_equivalent(const TrialParamSet &a,
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
bool build_shared_trigger_mask_soa_batch(const uuber::NativeContext &ctx,
                                         const TrialParamSet *base_params,
                                         const SharedTriggerPlan &plan,
                                         SharedTriggerMaskSoABatch &out);
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
void apply_trigger_q_soa_inplace_plan(const uuber::NativeContext &ctx,
                                      TrialParamSet &params,
                                      const SharedTriggerPlan &plan,
                                      int trigger_bit_idx, bool fail);
std::unique_ptr<TrialParamSet>
build_trial_params_from_df(const uuber::NativeContext &ctx,
                           const Rcpp::Nullable<Rcpp::DataFrame> &rows_opt);
