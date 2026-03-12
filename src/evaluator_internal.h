#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "distribution_core.h"
#include "forced_state.h"
#include "native_utils.h"
#include "trial_params.h"

using ExactSourceTimeMap = std::map<int, double>;
using SourceTimeBoundsMap = std::map<int, std::pair<double, double>>;

inline int accumulator_label_id_of(const uuber::NativeContext &ctx,
                                   int acc_idx) {
  if (acc_idx < 0 ||
      acc_idx >= static_cast<int>(ctx.accumulator_label_ids.size())) {
    return NA_INTEGER;
  }
  return ctx.accumulator_label_ids[static_cast<std::size_t>(acc_idx)];
}

inline int pool_label_id_of(const uuber::NativeContext &ctx, int pool_idx) {
  if (pool_idx < 0 || pool_idx >= static_cast<int>(ctx.pool_label_ids.size())) {
    return NA_INTEGER;
  }
  return ctx.pool_label_ids[static_cast<std::size_t>(pool_idx)];
}

inline int resolve_dense_node_idx(const uuber::NativeContext &ctx, int node_id) {
  auto it = ctx.ir.id_to_node_idx.find(node_id);
  if (it != ctx.ir.id_to_node_idx.end()) {
    return it->second;
  }
  if (node_id >= 0 && node_id < static_cast<int>(ctx.ir.nodes.size())) {
    return node_id;
  }
  return -1;
}

inline int resolve_dense_node_idx_required(const uuber::NativeContext &ctx,
                                           int node_id) {
  const int node_idx = resolve_dense_node_idx(ctx, node_id);
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    Rcpp::stop("IR node id %d not found in context", node_id);
  }
  return node_idx;
}

inline const uuber::IrNode &ir_node_required(const uuber::NativeContext &ctx,
                                             int node_idx) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    Rcpp::stop("IR node index %d out of bounds", node_idx);
  }
  return ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
}

inline uuber::LabelRef node_label_ref(const uuber::NativeContext &ctx,
                                      const uuber::IrNode &node) {
  uuber::LabelRef ref;
  const int event_idx = node.event_idx;
  if (event_idx < 0 || event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return ref;
  }
  const uuber::IrEvent &event =
      ctx.ir.events[static_cast<std::size_t>(event_idx)];
  ref.label_id = event.label_id;
  ref.acc_idx = event.acc_idx;
  ref.pool_idx = event.pool_idx;
  ref.outcome_idx = event.outcome_idx;
  return ref;
}

inline std::vector<int> ensure_source_ids(const uuber::NativeContext &ctx,
                                          const uuber::IrNode &node) {
  std::vector<int> ids;
  if (node.source_id_begin >= 0 && node.source_id_count > 0) {
    const int begin = node.source_id_begin;
    const int end = begin + node.source_id_count;
    if (end <= static_cast<int>(ctx.ir.node_source_label_ids.size())) {
      ids.assign(ctx.ir.node_source_label_ids.begin() + begin,
                 ctx.ir.node_source_label_ids.begin() + end);
    }
  }
  sort_unique(ids);
  return ids;
}

inline std::vector<int> ensure_source_ids(const uuber::NativeContext &ctx,
                                          int node_id_or_idx) {
  const int node_idx = resolve_dense_node_idx_required(ctx, node_id_or_idx);
  const uuber::IrNode &node = ir_node_required(ctx, node_idx);
  return ensure_source_ids(ctx, node);
}

std::string evaluator_component_cache_key(const uuber::NativeContext &ctx,
                                          int component_idx,
                                          const std::string &trial_key);

const TrialAccumulatorParams *evaluator_get_trial_param_entry(
    const TrialParamSet *trial_params, int acc_index);
void evaluator_resolve_event_numeric_params(
    const uuber::NativeAccumulator &acc, int acc_index,
    const TrialAccumulatorParams *override,
    const uuber::TrialParamsSoA *trial_params_soa, double &onset_out,
    double &q_out, uuber::AccDistParams &cfg_out);
bool evaluator_component_active_idx(
    const uuber::NativeAccumulator &acc, int component_idx,
    const TrialAccumulatorParams *override);

NodeEvalResult evaluator_eval_event_ref_idx(
    const uuber::NativeContext &ctx, const uuber::LabelRef &label_ref,
    std::uint32_t node_flags, double t, int component_idx, EvalNeed need,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_idx_context,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    const ForcedScopeFilter *forced_scope_filter,
    const uuber::BitsetState *forced_complete_bits,
    const uuber::BitsetState *forced_survive_bits,
    const std::unordered_map<int, int> *forced_label_id_to_bit_idx,
    const uuber::TrialParamsSoA *trial_params_soa,
    const ForcedStateView *forced_state_view);

struct NodeEvalState {
  NodeEvalState(const uuber::NativeContext &ctx_, double time,
                int component_idx_val = -1,
                const TrialParamSet *params_ptr = nullptr,
                const std::string &trial_key = std::string(),
                bool include_na = false, int outcome_idx_val = -1,
                const ExactSourceTimeMap *exact_source_times_ptr = nullptr,
                const SourceTimeBoundsMap *source_time_bounds_ptr = nullptr,
                const ForcedScopeFilter *forced_scope_filter_ptr = nullptr,
                const uuber::BitsetState *forced_complete_bits_ptr = nullptr,
                bool forced_complete_bits_ptr_valid = false,
                const uuber::BitsetState *forced_survive_bits_ptr = nullptr,
                bool forced_survive_bits_ptr_valid = false,
                uuber::KernelRuntimeState *external_kernel_runtime = nullptr,
                const ForcedStateView *forced_state_view = nullptr)
      : ctx(ctx_), t(time), trial_params(params_ptr),
        trial_type_key(
            evaluator_component_cache_key(ctx_, component_idx_val, trial_key)),
        include_na_donors(include_na), component_idx(component_idx_val),
        outcome_idx(outcome_idx_val),
        exact_source_times(exact_source_times_ptr),
        source_time_bounds(source_time_bounds_ptr),
        forced_scope_filter(nullptr),
        trial_params_soa(nullptr),
        forced_label_id_to_bit_idx(nullptr) {
    trial_params_soa = resolve_trial_params_soa(ctx, trial_params);
    bool forced_complete_valid_resolved = forced_complete_bits_ptr_valid;
    bool forced_survive_valid_resolved = forced_survive_bits_ptr_valid;
    const uuber::BitsetState *forced_complete_bits_resolved =
        forced_complete_bits_ptr;
    const uuber::BitsetState *forced_survive_bits_resolved =
        forced_survive_bits_ptr;
    forced_scope_filter = forced_scope_filter_ptr;
    forced_label_id_to_bit_idx = ctx.ir.label_id_to_bit_idx.empty()
                                     ? nullptr
                                     : &ctx.ir.label_id_to_bit_idx;
    if (forced_state_view) {
      if (forced_state_view->scope_filter) {
        forced_scope_filter = forced_state_view->scope_filter;
      }
      if (forced_state_view->label_id_to_bit_idx) {
        forced_label_id_to_bit_idx = forced_state_view->label_id_to_bit_idx;
      }
      if (forced_state_view->forced_complete_bits) {
        forced_complete_bits_resolved = forced_state_view->forced_complete_bits;
      }
      if (forced_state_view->forced_survive_bits) {
        forced_survive_bits_resolved = forced_state_view->forced_survive_bits;
      }
      if (forced_state_view->forced_complete_bits_valid) {
        forced_complete_valid_resolved =
            *forced_state_view->forced_complete_bits_valid;
      }
      if (forced_state_view->forced_survive_bits_valid) {
        forced_survive_valid_resolved =
            *forced_state_view->forced_survive_bits_valid;
      }
    }

    if (forced_complete_bits_resolved && forced_complete_valid_resolved) {
      forced_complete_bits = *forced_complete_bits_resolved;
      forced_complete_bits_valid = true;
    } else {
      (void)ensure_forced_bitset_capacity(ctx, forced_complete_bits,
                                          forced_complete_bits_valid);
    }
    if (forced_survive_bits_resolved && forced_survive_valid_resolved) {
      forced_survive_bits = *forced_survive_bits_resolved;
      forced_survive_bits_valid = true;
    } else {
      (void)ensure_forced_bitset_capacity(ctx, forced_survive_bits,
                                          forced_survive_bits_valid);
    }
    forced_state = make_forced_state_view(
        forced_scope_filter, &forced_complete_bits, &forced_complete_bits_valid,
        &forced_survive_bits, &forced_survive_bits_valid,
        forced_label_id_to_bit_idx);
    if (ctx.kernel_program.valid) {
      if (external_kernel_runtime) {
        kernel_runtime_ptr = external_kernel_runtime;
        if (!kernel_runtime_ptr->initialized ||
            kernel_runtime_ptr->program != &ctx.kernel_program ||
            kernel_runtime_ptr->slots.size() != ctx.kernel_program.ops.size()) {
          uuber::reset_kernel_runtime(ctx.kernel_program, *kernel_runtime_ptr);
        }
      } else {
        uuber::reset_kernel_runtime(ctx.kernel_program, kernel_runtime_storage);
        kernel_runtime_ptr = &kernel_runtime_storage;
      }
      kernel_runtime_ready = (kernel_runtime_ptr != nullptr);
    }
  }

  const uuber::NativeContext &ctx;
  double t;
  const TrialParamSet *trial_params;
  std::string trial_type_key;
  bool include_na_donors;
  int component_idx;
  int outcome_idx;
  const ExactSourceTimeMap *exact_source_times;
  const SourceTimeBoundsMap *source_time_bounds;
  const ForcedScopeFilter *forced_scope_filter;
  const uuber::TrialParamsSoA *trial_params_soa;
  const std::unordered_map<int, int> *forced_label_id_to_bit_idx;
  ForcedStateView forced_state{};
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  uuber::KernelRuntimeState kernel_runtime_storage;
  uuber::KernelRuntimeState *kernel_runtime_ptr{nullptr};
  bool kernel_runtime_ready{false};
};

inline bool kernel_runtime_cache_safe(const NodeEvalState &state) {
  return state.forced_scope_filter == nullptr &&
         !(state.forced_complete_bits_valid && state.forced_complete_bits.any()) &&
         !(state.forced_survive_bits_valid && state.forced_survive_bits.any()) &&
         state.exact_source_times == nullptr &&
         state.source_time_bounds == nullptr;
}

inline void apply_transition_mask_words(
    const uuber::NativeContext &ctx, int source_mask_begin,
    int source_mask_count, int invalidate_slot, uuber::BitsetState &bits_out,
    bool &bits_valid, uuber::KernelRuntimeState *kernel_runtime,
    uuber::KernelBatchRuntimeState *kernel_batch_runtime = nullptr) {
  if (source_mask_begin < 0 || source_mask_count <= 0 ||
      source_mask_begin + source_mask_count >
          static_cast<int>(ctx.ir.node_source_masks.size())) {
    Rcpp::stop("IR guard transition mask invalid begin=%d count=%d",
               source_mask_begin, source_mask_count);
  }
  if (!ensure_forced_bitset_capacity(ctx, bits_out, bits_valid)) {
    Rcpp::stop("IR guard transition bitset unavailable for mask update");
  }
  const std::uint64_t *mask_words =
      &ctx.ir.node_source_masks[static_cast<std::size_t>(source_mask_begin)];
  bits_out.or_words(mask_words, source_mask_count);
  if (kernel_runtime && invalidate_slot >= 0) {
    uuber::invalidate_kernel_runtime_from_slot(*kernel_runtime, invalidate_slot);
  }
  if (kernel_batch_runtime && invalidate_slot >= 0) {
    uuber::invalidate_kernel_batch_runtime_from_slot(*kernel_batch_runtime,
                                                     invalidate_slot);
  }
}

struct IntegrationSettings {
  double rel_tol{1e-5};
  double abs_tol{1e-6};
  int max_depth{12};
};

struct CanonicalIntegrationSettings {
  double rel_tol{1e-5};
  double abs_tol{1e-6};
  int max_depth{12};
  double effective_tol{1e-5};
};

inline CanonicalIntegrationSettings
canonicalize_integration_settings(const IntegrationSettings &settings) {
  CanonicalIntegrationSettings canonical;
  canonical.rel_tol =
      (settings.rel_tol > 0.0 && std::isfinite(settings.rel_tol))
          ? settings.rel_tol
          : 1e-5;
  canonical.abs_tol =
      (settings.abs_tol > 0.0 && std::isfinite(settings.abs_tol))
          ? settings.abs_tol
          : 1e-6;
  canonical.max_depth = (settings.max_depth > 0) ? settings.max_depth : 12;
  canonical.effective_tol =
      std::max(std::max(canonical.rel_tol, canonical.abs_tol), 1e-10);
  return canonical;
}

struct GuardEvalInput {
  const uuber::NativeContext &ctx;
  int node_idx{-1};
  int guard_transition_idx{-1};
  const uuber::KernelGuardTransition *guard_transition{nullptr};
  int component_idx{-1};
  const std::string *trial_type_key{nullptr};
  ForcedStateView forced_state{};
  const uuber::BitsetState *forced_complete_bits{nullptr};
  bool forced_complete_bits_valid{false};
  const uuber::BitsetState *forced_survive_bits{nullptr};
  bool forced_survive_bits_valid{false};
  const std::unordered_map<int, int> *forced_label_id_to_bit_idx{nullptr};
  const ForcedScopeFilter *forced_scope_filter{nullptr};
  ForcedScopeFilter local_scope_filter{};
  bool has_scoped_forced{false};
  const TrialParamSet *trial_params;
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  const ExactSourceTimeMap *exact_source_times{nullptr};
  const SourceTimeBoundsMap *source_time_bounds{nullptr};
};

bool evaluator_eval_node_with_forced_state_view_batch(
    const uuber::NativeContext &ctx, int node_id_or_idx,
    const std::vector<double> &times, int component_idx, EvalNeed need,
    const TrialParamSet *trial_params, const std::string &trial_key,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    const ForcedStateView &forced_state,
    uuber::KernelNodeBatchValues &out_values);
GuardEvalInput evaluator_make_guard_input(
    const uuber::NativeContext &ctx, int node_idx, int component_idx,
    const std::string *trial_type_key, const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *trial_params_soa,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    const ForcedScopeFilter *forced_scope_filter,
    const uuber::BitsetState *forced_complete_bits,
    const uuber::BitsetState *forced_survive_bits,
    const std::unordered_map<int, int> *forced_label_id_to_bit_idx,
    const ForcedStateView *forced_state_view);

inline GuardEvalInput make_guard_input_forced_state(
    const uuber::NativeContext &ctx, int node_idx, int component_idx,
    const std::string *trial_type_key, const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *trial_params_soa,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    const ForcedStateView &forced_state_view) {
  return evaluator_make_guard_input(
      ctx, node_idx, component_idx, trial_type_key, trial_params,
      trial_params_soa, exact_source_times, source_time_bounds, nullptr,
      nullptr, nullptr, nullptr, &forced_state_view);
}

std::vector<int> evaluator_gather_blocker_sources(
    const uuber::NativeContext &ctx, int guard_node_idx);
double evaluator_guard_effective_survival_internal(
    const GuardEvalInput &input, double t,
    const IntegrationSettings &settings);
NodeEvalResult evaluator_eval_node_recursive_dense(int node_idx,
                                                   NodeEvalState &state,
                                                   EvalNeed need);
uuber::KernelEventEvalFn evaluator_make_kernel_event_eval(NodeEvalState &state);
uuber::KernelEventBatchEvalFn evaluator_make_kernel_event_eval_batch(
    NodeEvalState &state);
uuber::KernelGuardEvalFn evaluator_make_kernel_guard_eval(NodeEvalState &state);
uuber::KernelGuardBatchEvalFn evaluator_make_kernel_guard_eval_batch(
    NodeEvalState &state);
bool evaluator_eval_node_batch_with_state_dense(
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    EvalNeed need, uuber::KernelBatchRuntimeState &batch_runtime,
    uuber::KernelNodeBatchValues &out_values);
double evaluator_evaluate_survival_with_forced(
    int node_id, const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid, int component_idx, double t,
    const uuber::NativeContext &ctx, const std::string &trial_key,
    const TrialParamSet *trial_params,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    uuber::KernelRuntimeState *kernel_runtime);
double competitor_survival_internal(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits_in,
    bool forced_complete_bits_in_valid,
    const uuber::BitsetState *forced_survive_bits_in,
    bool forced_survive_bits_in_valid, const std::string &trial_type_key,
    const TrialParamSet *trial_params,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    uuber::KernelRuntimeState *kernel_runtime);
double competitor_survival_from_state_compiled_ops(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    NodeEvalState &state,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr);
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
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    uuber::KernelBatchRuntimeState *kernel_batch_runtime);
void invalidate_kernel_runtime_root(NodeEvalState &state);
