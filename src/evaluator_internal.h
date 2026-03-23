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
#include "tree_eval.h"
#include "time_constraints.h"
#include "trial_params.h"
#include "vector_lane_utils.h"

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

inline uuber::LabelRef evaluator_make_onset_source_ref(
    const uuber::NativeContext &ctx, int onset_kind, int onset_source_acc_idx,
    int onset_source_pool_idx) {
  uuber::LabelRef ref;
  if (onset_kind == uuber::ONSET_AFTER_ACCUMULATOR) {
    if (onset_source_acc_idx >= 0 &&
        onset_source_acc_idx < static_cast<int>(ctx.accumulators.size())) {
      ref.acc_idx = onset_source_acc_idx;
      ref.label_id = accumulator_label_id_of(ctx, onset_source_acc_idx);
    }
  } else if (onset_kind == uuber::ONSET_AFTER_POOL) {
    if (onset_source_pool_idx >= 0 &&
        onset_source_pool_idx < static_cast<int>(ctx.pools.size())) {
      ref.pool_idx = onset_source_pool_idx;
      ref.label_id = pool_label_id_of(ctx, onset_source_pool_idx);
    }
  }
  return ref;
}

inline bool same_acc_batch_params_except_q(const uuber::TrialParamsSoA &lhs,
                                           const uuber::TrialParamsSoA &rhs,
                                           int acc_idx) {
  if (acc_idx < 0 || acc_idx >= lhs.n_acc || acc_idx >= rhs.n_acc) {
    return false;
  }
  const std::size_t idx = static_cast<std::size_t>(acc_idx);
  return lhs.valid && rhs.valid &&
         idx < lhs.dist_code.size() && idx < rhs.dist_code.size() &&
         idx < lhs.onset.size() && idx < rhs.onset.size() &&
         idx < lhs.t0.size() && idx < rhs.t0.size() &&
         idx < lhs.p1.size() && idx < rhs.p1.size() &&
         idx < lhs.p2.size() && idx < rhs.p2.size() &&
         idx < lhs.p3.size() && idx < rhs.p3.size() &&
         idx < lhs.p4.size() && idx < rhs.p4.size() &&
         idx < lhs.p5.size() && idx < rhs.p5.size() &&
         idx < lhs.p6.size() && idx < rhs.p6.size() &&
         idx < lhs.p7.size() && idx < rhs.p7.size() &&
         idx < lhs.p8.size() && idx < rhs.p8.size() &&
         lhs.dist_code[idx] == rhs.dist_code[idx] &&
         canonical_double_bits(lhs.onset[idx]) ==
             canonical_double_bits(rhs.onset[idx]) &&
         canonical_double_bits(lhs.t0[idx]) ==
             canonical_double_bits(rhs.t0[idx]) &&
         canonical_double_bits(lhs.p1[idx]) ==
             canonical_double_bits(rhs.p1[idx]) &&
         canonical_double_bits(lhs.p2[idx]) ==
             canonical_double_bits(rhs.p2[idx]) &&
         canonical_double_bits(lhs.p3[idx]) ==
             canonical_double_bits(rhs.p3[idx]) &&
         canonical_double_bits(lhs.p4[idx]) ==
             canonical_double_bits(rhs.p4[idx]) &&
         canonical_double_bits(lhs.p5[idx]) ==
             canonical_double_bits(rhs.p5[idx]) &&
         canonical_double_bits(lhs.p6[idx]) ==
             canonical_double_bits(rhs.p6[idx]) &&
         canonical_double_bits(lhs.p7[idx]) ==
             canonical_double_bits(rhs.p7[idx]) &&
         canonical_double_bits(lhs.p8[idx]) ==
             canonical_double_bits(rhs.p8[idx]);
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

inline bool ir_mask_has_component(const uuber::NativeContext &ctx,
                                  int mask_offset, int component_idx) {
  if (component_idx < 0) {
    return true;
  }
  if (!ctx.ir.valid || ctx.ir.component_mask_words <= 0 || mask_offset < 0) {
    return true;
  }
  const int word_idx = component_idx / 64;
  const int bit = component_idx % 64;
  if (word_idx < 0 || word_idx >= ctx.ir.component_mask_words) {
    return false;
  }
  const std::size_t idx = static_cast<std::size_t>(mask_offset + word_idx);
  if (idx >= ctx.ir.component_masks.size()) {
    return false;
  }
  const std::uint64_t word = ctx.ir.component_masks[idx];
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

std::vector<int> evaluator_gather_blocker_sources(
    const uuber::NativeContext &ctx, int guard_node_idx);

inline std::vector<int>
relevant_source_ids_for_batch_node(const uuber::NativeContext &ctx,
                                   int node_idx) {
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

inline bool evaluator_resolve_label_time_constraint(
    const TimeConstraintMap *time_constraints, int label_id, bool &has_exact,
    double &exact_time, bool &has_bounds, double &bound_lower,
    double &bound_upper) {
  has_exact = false;
  exact_time = std::numeric_limits<double>::quiet_NaN();
  has_bounds = false;
  bound_lower = 0.0;
  bound_upper = std::numeric_limits<double>::infinity();
  const SourceTimeConstraint *constraint =
      time_constraints_find(time_constraints, label_id);
  if (constraint == nullptr) {
    return false;
  }
  if (constraint->has_exact) {
    has_exact = true;
    exact_time = constraint->exact_time;
  }
  if (constraint->has_lower) {
    bound_lower = std::max(bound_lower, constraint->lower);
    has_bounds = true;
  }
  if (constraint->has_upper) {
    bound_upper = std::min(bound_upper, constraint->upper);
    has_bounds = true;
  }
  return true;
}

struct CompactTimeConstraintLookup {
  std::vector<int> source_ids;
  std::vector<SourceTimeConstraint> constraints;
};

inline const SourceTimeConstraint *compact_time_constraint_find(
    const CompactTimeConstraintLookup *lookup, int source_id) {
  if (lookup == nullptr || lookup->source_ids.empty() ||
      lookup->source_ids.size() != lookup->constraints.size()) {
    return nullptr;
  }
  auto it = std::lower_bound(lookup->source_ids.begin(),
                             lookup->source_ids.end(), source_id);
  if (it == lookup->source_ids.end() || *it != source_id) {
    return nullptr;
  }
  const std::size_t idx = static_cast<std::size_t>(
      std::distance(lookup->source_ids.begin(), it));
  return &lookup->constraints[idx];
}

inline void build_compact_time_constraint_lookup(
    const TimeConstraintMap *time_constraints,
    const std::vector<int> &relevant_source_ids,
    CompactTimeConstraintLookup &out) {
  out.source_ids = relevant_source_ids;
  sort_unique(out.source_ids);
  out.constraints.assign(out.source_ids.size(), SourceTimeConstraint{});
  if (time_constraints == nullptr || time_constraints->empty() ||
      out.source_ids.empty()) {
    return;
  }
  for (std::size_t i = 0; i < out.source_ids.size(); ++i) {
    const SourceTimeConstraint *constraint =
        time_constraints_find(time_constraints, out.source_ids[i]);
    if (constraint != nullptr) {
      out.constraints[i] = *constraint;
    }
  }
}

inline bool evaluator_state_contains_survive_at(
    const ForcedStateView &forced_state,
    const TimeConstraintMap *time_constraints, int label_id, double t) {
  return forced_state_contains_survive(forced_state, label_id) ||
         time_constraints_contains_survive_at(time_constraints, label_id, t);
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
    const TimeConstraintMap *time_constraints,
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
                const TimeConstraintMap *time_constraints_ptr = nullptr,
                const ForcedScopeFilter *forced_scope_filter_ptr = nullptr,
                const uuber::BitsetState *forced_complete_bits_ptr = nullptr,
                bool forced_complete_bits_ptr_valid = false,
                const uuber::BitsetState *forced_survive_bits_ptr = nullptr,
                bool forced_survive_bits_ptr_valid = false,
                const ForcedStateView *forced_state_view = nullptr)
      : ctx(ctx_), t(time), trial_params(params_ptr),
        trial_type_key(
            evaluator_component_cache_key(ctx_, component_idx_val, trial_key)),
        include_na_donors(include_na), component_idx(component_idx_val),
        outcome_idx(outcome_idx_val),
        forced_scope_filter(nullptr),
        trial_params_soa(nullptr),
        forced_label_id_to_bit_idx(nullptr) {
    if (time_constraints_ptr != nullptr && !time_constraints_ptr->empty()) {
      time_constraints = *time_constraints_ptr;
    }
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
  }

  const uuber::NativeContext &ctx;
  double t;
  const TrialParamSet *trial_params;
  std::string trial_type_key;
  bool include_na_donors;
  int component_idx;
  int outcome_idx;
  TimeConstraintMap time_constraints;
  const ForcedScopeFilter *forced_scope_filter;
  const uuber::TrialParamsSoA *trial_params_soa;
  const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch{
      nullptr};
  const std::unordered_map<int, int> *forced_label_id_to_bit_idx;
  ForcedStateView forced_state{};
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
};

template <typename PointVec, typename SharedFn>
inline bool eval_vector_lanes_with_shared_node_state_group(
    const uuber::NativeContext &ctx, const PointVec &points,
    const std::vector<std::uint8_t> *active_mask, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::size_t anchor_idx, const uuber::VectorRelevantStateStorage &storage,
    SharedFn &&shared_eval,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  if (anchor_idx >= points.size()) {
    return false;
  }
  const uuber::VectorLaneRef anchor = uuber::vector_lane_ref(points, anchor_idx);
  NodeEvalState shared_state(
      ctx, anchor.t, component_idx, trial_params, trial_type_key, false, -1,
      &storage.time_constraints, nullptr,
      storage.forced_complete_bits_valid ? &storage.forced_complete_bits
                                         : nullptr,
      storage.forced_complete_bits_valid,
      storage.forced_survive_bits_valid ? &storage.forced_survive_bits
                                        : nullptr,
      storage.forced_survive_bits_valid);
  shared_state.trial_params_soa =
      uniform_trial_params_soa != nullptr
          ? uniform_trial_params_soa
          : uuber::resolve_vector_lane_trial_params_soa(
                ctx, anchor, trial_params, shared_state.trial_params_soa);
  if (!shared_state.trial_params_soa) {
    return false;
  }

  std::vector<const uuber::TrialParamsSoA *> point_trial_params_soa;
  if (uniform_trial_params_soa == nullptr) {
    std::size_t first_trial_param_mismatch = points.size();
    const uuber::TrialParamsSoA *mismatch_trial_params_soa = nullptr;
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (!uuber::vector_lane_active(active_mask, i)) {
        continue;
      }
      const uuber::TrialParamsSoA *point_soa =
          uuber::resolve_vector_lane_trial_params_soa(
              ctx, uuber::vector_lane_ref(points, i), trial_params,
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
        if (!uuber::vector_lane_active(active_mask, i)) {
          continue;
        }
        point_trial_params_soa[i] =
            uuber::resolve_vector_lane_trial_params_soa(
                ctx, uuber::vector_lane_ref(points, i), trial_params,
                shared_state.trial_params_soa);
      }
      shared_state.trial_params_soa_batch = &point_trial_params_soa;
    }
  }

  return shared_eval(shared_state, active_mask);
}

template <typename PointVec, typename SharedFn, typename FallbackFn,
          typename RecordPartitionFn>
inline bool eval_vector_lanes_with_shared_node_state(
    const uuber::NativeContext &ctx, const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    SharedFn &&shared_eval, FallbackFn &&fallback_eval,
    RecordPartitionFn &&record_partition,
    bool fallback_on_singleton_partition = false,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr) {
  auto &shared_eval_ref = shared_eval;
  auto &fallback_eval_ref = fallback_eval;
  auto &record_partition_ref = record_partition;
  return uuber::eval_vector_lanes_with_shared_state(
      ctx, points, active_mask, relevant_source_ids,
      [&](std::size_t group_anchor_idx,
          const uuber::VectorRelevantStateStorage &group_storage,
          const std::vector<std::uint8_t> *group_mask) -> bool {
        return eval_vector_lanes_with_shared_node_state_group(
            ctx, points, group_mask, component_idx, trial_params,
            trial_type_key, group_anchor_idx, group_storage, shared_eval_ref,
            uniform_trial_params_soa);
      },
      [&]() -> bool { return fallback_eval_ref(); },
      [&](std::uint64_t active_count, std::uint64_t group_count) {
        record_partition_ref(active_count, group_count);
      },
      fallback_on_singleton_partition);
}

inline void apply_transition_mask_words(
    const uuber::NativeContext &ctx, int source_mask_begin,
    int source_mask_count, int invalidate_slot, uuber::BitsetState &bits_out,
    bool &bits_valid) {
  (void)invalidate_slot;
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
  const TimeConstraintMap *time_constraints{nullptr};
  const CompactTimeConstraintLookup *time_constraint_lookup{nullptr};
  CompactTimeConstraintLookup local_time_constraint_lookup{};
};

bool evaluator_eval_node_with_forced_state_view_batch(
    const uuber::NativeContext &ctx, int node_id_or_idx,
    const std::vector<double> &times, int component_idx, EvalNeed need,
    const TrialParamSet *trial_params, const std::string &trial_key,
    const TimeConstraintMap *time_constraints,
    const ForcedStateView &forced_state,
    uuber::KernelNodeBatchValues &out_values);
GuardEvalInput evaluator_make_guard_input(
    const uuber::NativeContext &ctx, int node_idx, int component_idx,
    const std::string *trial_type_key, const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *trial_params_soa,
    const TimeConstraintMap *time_constraints,
    const ForcedScopeFilter *forced_scope_filter,
    const uuber::BitsetState *forced_complete_bits,
    const uuber::BitsetState *forced_survive_bits,
    const std::unordered_map<int, int> *forced_label_id_to_bit_idx,
    const ForcedStateView *forced_state_view);

inline GuardEvalInput make_guard_input_forced_state(
    const uuber::NativeContext &ctx, int node_idx, int component_idx,
    const std::string *trial_type_key, const TrialParamSet *trial_params,
    const uuber::TrialParamsSoA *trial_params_soa,
    const TimeConstraintMap *time_constraints,
    const ForcedStateView &forced_state_view) {
  return evaluator_make_guard_input(
      ctx, node_idx, component_idx, trial_type_key, trial_params,
      trial_params_soa, time_constraints, nullptr,
      nullptr, nullptr, nullptr, &forced_state_view);
}

std::vector<int> evaluator_gather_blocker_sources(
    const uuber::NativeContext &ctx, int guard_node_idx);
double evaluator_guard_effective_survival_internal(
    const GuardEvalInput &input, double t,
    const IntegrationSettings &settings);
bool evaluator_guard_cdf_batch_prepared(const GuardEvalInput &input,
                                        const std::vector<double> &times,
                                        std::vector<double> &cdf_out);
NodeEvalResult evaluator_eval_node_recursive_dense(int node_idx,
                                                   NodeEvalState &state,
                                                   EvalNeed need);
uuber::KernelEventBatchEvalFn evaluator_make_kernel_event_eval_batch(
    NodeEvalState &state);
uuber::KernelGuardBatchEvalFn evaluator_make_kernel_guard_eval_batch(
    NodeEvalState &state);
bool evaluator_eval_node_batch_with_state_dense(
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values);
double evaluator_evaluate_survival_with_forced(
    int node_id, const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid, int component_idx, double t,
    const uuber::NativeContext &ctx, const std::string &trial_key,
    const TrialParamSet *trial_params,
    const TimeConstraintMap *time_constraints);
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
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_batch =
        nullptr);
