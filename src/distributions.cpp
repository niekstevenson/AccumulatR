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
#include "component_plan.h"
#include "competitor_cache.h"
#include "context.h"
#include "coupling_eval.h"
#include "coupling_program.h"
#include "distribution_core.h"
#include "dist_vector.h"
#include "exact_outcome_density.h"
#include "evaluator_internal.h"
#include "forced_state.h"
#include "integrate.h"
#include "kernel_executor.h"
#include "native_utils.h"
#include "pool_math.h"
#include "proto.h"
#include "quadrature_batch.h"
#include "ranked_transitions.h"
#include "runtime_stats.h"
#include "trial_params.h"

using uuber::AccDistParams;
using uuber::ComponentMap;
using uuber::LabelRef;

namespace {

constexpr double kDefaultRelTol = 1e-5;
constexpr double kDefaultAbsTol = 1e-6;
constexpr int kDefaultMaxDepth = 12;
constexpr int kOutcomeIdxNA = -2;

inline int resolve_observed_label_id(const uuber::NativeContext &ctx,
                                     const std::string &label) {
  if (label.empty())
    return -1;
  auto it = ctx.label_to_id.find(label);
  if (it == ctx.label_to_id.end())
    return -1;
  return it->second;
}

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

inline int outcome_index_of(const uuber::NativeContext &ctx,
                            const std::string &label) {
  if (label.empty())
    return -1;
  int observed_label_id = resolve_observed_label_id(ctx, label);
  if (observed_label_id < 0)
    return -1;
  auto it_out = ctx.ir.label_id_to_outcomes.find(observed_label_id);
  if (it_out == ctx.ir.label_id_to_outcomes.end() || it_out->second.empty()) {
    return -1;
  }
  return it_out->second.front();
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

inline bool state_allows_exact_scalar_dispatch(const uuber::NativeContext &ctx,
                                               const NodeEvalState &state) {
  if (state.include_na_donors || state.outcome_idx < 0 ||
      state.outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return !state.include_na_donors;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(state.outcome_idx)];
  return info.alias_sources.empty() && info.guess_donors.empty();
}

inline bool coupling_program_requires_exact_scenario_eval(
    const OutcomeCouplingProgram *program) {
  return program != nullptr && program->valid &&
         program->kind == OutcomeCouplingOpKind::GenericNodeIntegral &&
         program->generic.requires_exact_scenario_eval;
}

inline OutcomeCouplingProgram resolve_density_coupling_program(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return OutcomeCouplingProgram{};
  }
  const int node_id =
      ctx.ir.nodes[static_cast<std::size_t>(node_idx)].node_id;
  return resolve_outcome_coupling_program_with_generic(
      ctx, node_id >= 0 ? node_id : node_idx, competitor_ids, false);
}

inline bool should_use_exact_density_program(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, const NodeEvalState &state,
    const OutcomeCouplingProgram *program) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size()) ||
      !state_allows_exact_scalar_dispatch(ctx, state)) {
    return false;
  }
  OutcomeCouplingProgram local_program;
  const OutcomeCouplingProgram *program_ptr = program;
  if (program_ptr == nullptr) {
    local_program = resolve_density_coupling_program(ctx, node_idx, competitor_ids);
    program_ptr = local_program.valid ? &local_program : nullptr;
  }
  return coupling_program_requires_exact_scenario_eval(program_ptr);
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
    int max_depth, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan, TrialParamSet *trigger_scratch);

double node_density_with_competitors_internal(
    const uuber::NativeContext &ctx, int node_id, double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const TimeConstraintMap *time_constraints = nullptr,
    uuber::KernelRuntimeState *kernel_runtime = nullptr);

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
    uuber::KernelBatchRuntimeState *kernel_batch_runtime = nullptr);

double kernel_node_density_entry_idx(
    const uuber::NativeContext &ctx, int node_id, double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const TimeConstraintMap *time_constraints = nullptr,
    uuber::KernelRuntimeState *kernel_runtime = nullptr);

double node_density_entry_idx(
    const uuber::NativeContext &ctx, int node_id, double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan = nullptr,
    TrialParamSet *scratch_params = nullptr,
    bool use_shared_trigger_eval = false,
    const TimeConstraintMap *time_constraints = nullptr,
    uuber::KernelRuntimeState *kernel_runtime = nullptr);

double evaluate_outcome_density_idx(
    const uuber::NativeContext &ctx, int outcome_idx, double t, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_label_context_idx = -1,
    const SharedTriggerPlan *trigger_plan = nullptr,
    TrialParamSet *trigger_scratch = nullptr) {
  if (outcome_idx < 0 ||
      outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return 0.0;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
  if (info.node_id < 0)
    return 0.0;
  std::vector<int> comp_filtered;
  const std::vector<int> &comp_use =
      filter_competitor_ids(ctx, info.competitor_ids, component_idx,
                            comp_filtered);
  return node_density_entry_idx(
      ctx, info.node_id, t, component_idx, nullptr, false, nullptr, false,
      comp_use, trial_params, trial_type_key, include_na_donors,
      outcome_label_context_idx, trigger_plan, trigger_scratch,
      trigger_plan != nullptr);
}

double accumulate_component_guess_density_idx(
    const uuber::NativeContext &ctx, int target_label_id,
    int target_outcome_idx, bool target_is_guess, double t, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const SharedTriggerPlan *trigger_plan = nullptr,
    TrialParamSet *trigger_scratch = nullptr) {
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.component_info.size())) {
    return 0.0;
  }
  const uuber::ComponentGuessPolicy &guess =
      ctx.component_info[static_cast<std::size_t>(component_idx)].guess;
  if (!guess.valid)
    return 0.0;
  if (guess.target_is_guess) {
    if (!target_is_guess)
      return 0.0;
  } else if (guess.target_outcome_idx >= 0) {
    if (target_outcome_idx != guess.target_outcome_idx)
      return 0.0;
  } else if (guess.target_label_id != NA_INTEGER) {
    if (target_label_id == NA_INTEGER || target_label_id != guess.target_label_id)
      return 0.0;
  } else if (!guess.target.empty()) {
    // Unknown label target that could not be resolved at context build time.
    return 0.0;
  }
  double total = 0.0;
  for (const auto &entry : guess.keep_weights) {
    int donor_outcome_idx = entry.first;
    if (donor_outcome_idx < 0)
      continue;
    double keep_prob = entry.second;
    double release = 1.0 - keep_prob;
    if (release <= 0.0)
      continue;
    total += release * evaluate_outcome_density_idx(
                          ctx, donor_outcome_idx, t, component_idx, trial_params,
                          trial_type_key, false, donor_outcome_idx,
                          trigger_plan, trigger_scratch);
  }
  return total;
}

double accumulate_outcome_alias_density(const uuber::NativeContext &ctx,
                                        int target_outcome_idx, double t,
                                        int component_idx,
                                        const TrialParamSet *trial_params,
                                        const std::string &trial_type_key,
                                        bool include_na_donors) {
  if (target_outcome_idx < 0 ||
      target_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return 0.0;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(target_outcome_idx)];
  double total = 0.0;
  for (int source_idx : info.alias_sources) {
    total += evaluate_outcome_density_idx(
        ctx, source_idx, t, component_idx, trial_params,
        trial_type_key, include_na_donors, source_idx);
  }
  for (const auto &donor : info.guess_donors) {
    if (donor.rt_policy == "na" && !include_na_donors)
      continue;
    if (donor.outcome_idx < 0)
      continue;
    total += donor.weight * evaluate_outcome_density_idx(
                                ctx, donor.outcome_idx, t, component_idx,
                                trial_params, trial_type_key,
                                include_na_donors, donor.outcome_idx);
  }
  return total;
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

inline bool resolve_label_time_constraint(const TimeConstraintMap *time_constraints,
                                          int label_id,
                                          bool &has_exact,
                                          double &exact_time,
                                          bool &has_bounds,
                                          double &bound_lower,
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

inline bool state_contains_survive_at(const ForcedStateView &forced_state,
                                      const TimeConstraintMap *time_constraints,
                                      int label_id, double t) {
  return forced_state_contains_survive(forced_state, label_id) ||
         time_constraints_contains_survive_at(time_constraints, label_id, t);
}

NodeEvalResult
eval_event_ref_idx(const uuber::NativeContext &ctx, const LabelRef &label_ref,
                   std::uint32_t node_flags, double t,
                   int component_idx, EvalNeed need,
                   const TrialParamSet *trial_params = nullptr,
                   const std::string &trial_type_key = std::string(),
                   bool include_na_donors = false,
                   int outcome_idx_context = -1,
                   const TimeConstraintMap *time_constraints = nullptr,
                   const ForcedScopeFilter *forced_scope_filter = nullptr,
                   const uuber::BitsetState *forced_complete_bits = nullptr,
                   const uuber::BitsetState *forced_survive_bits = nullptr,
                   const std::unordered_map<int, int>
                       *forced_label_id_to_bit_idx = nullptr,
                   const uuber::TrialParamsSoA *trial_params_soa = nullptr,
                   const ForcedStateView *forced_state_view = nullptr) {
  bool forced_complete_bits_valid_fallback = (forced_complete_bits != nullptr);
  bool forced_survive_bits_valid_fallback = (forced_survive_bits != nullptr);
  const ForcedScopeFilter *resolved_scope_filter = forced_scope_filter;
  const uuber::BitsetState *resolved_complete_bits = forced_complete_bits;
  const uuber::BitsetState *resolved_survive_bits = forced_survive_bits;
  const std::unordered_map<int, int> *resolved_label_id_to_bit_idx =
      forced_label_id_to_bit_idx;
  const bool *resolved_complete_bits_valid =
      &forced_complete_bits_valid_fallback;
  const bool *resolved_survive_bits_valid = &forced_survive_bits_valid_fallback;
  if (forced_state_view) {
    if (forced_state_view->scope_filter) {
      resolved_scope_filter = forced_state_view->scope_filter;
    }
    if (forced_state_view->forced_complete_bits) {
      resolved_complete_bits = forced_state_view->forced_complete_bits;
    }
    if (forced_state_view->forced_survive_bits) {
      resolved_survive_bits = forced_state_view->forced_survive_bits;
    }
    if (forced_state_view->label_id_to_bit_idx) {
      resolved_label_id_to_bit_idx = forced_state_view->label_id_to_bit_idx;
    }
    if (forced_state_view->forced_complete_bits_valid) {
      resolved_complete_bits_valid =
          forced_state_view->forced_complete_bits_valid;
    }
    if (forced_state_view->forced_survive_bits_valid) {
      resolved_survive_bits_valid = forced_state_view->forced_survive_bits_valid;
    }
  }
  const ForcedStateView forced_state =
      make_forced_state_view(resolved_scope_filter, resolved_complete_bits,
                             resolved_complete_bits_valid,
                             resolved_survive_bits, resolved_survive_bits_valid,
                             resolved_label_id_to_bit_idx);
  const bool is_special_deadline =
      (node_flags & uuber::IR_NODE_FLAG_SPECIAL_DEADLINE) != 0u;
  const bool is_special_guess =
      (node_flags & uuber::IR_NODE_FLAG_SPECIAL_GUESS) != 0u;
  if (is_special_deadline) {
    double survival = std::numeric_limits<double>::quiet_NaN();
    double cdf = std::numeric_limits<double>::quiet_NaN();
    if (needs_survival(need) || needs_cdf(need)) {
      survival = 0.0;
      cdf = 1.0;
      if (std::isfinite(t)) {
        if (t < 0.0) {
          survival = 1.0;
          cdf = 0.0;
        } else {
          survival = 0.0;
          cdf = 1.0;
        }
      }
    }
    return make_node_result(need, 0.0, survival, cdf);
  }
  double donor_density = 0.0;
  if (needs_density(need)) {
    int target_outcome_idx =
        (outcome_idx_context >= 0) ? outcome_idx_context : label_ref.outcome_idx;
    int target_label_id = label_ref.label_id;
    if (outcome_idx_context >= 0 &&
        outcome_idx_context < static_cast<int>(ctx.outcome_label_ids.size())) {
      target_label_id =
          ctx.outcome_label_ids[static_cast<std::size_t>(outcome_idx_context)];
    }
    donor_density += accumulate_component_guess_density_idx(
        ctx, target_label_id, target_outcome_idx, is_special_guess, t,
        component_idx, trial_params, trial_type_key);
    donor_density += accumulate_outcome_alias_density(
        ctx, target_outcome_idx, t, component_idx, trial_params, trial_type_key,
        include_na_donors);
  }
  if (is_special_guess) {
    return make_node_result(need, donor_density, 0.0, 1.0);
  }
  int label_idx = label_ref.label_id;
  bool self_has_exact_source_time = false;
  double self_exact_source_time = std::numeric_limits<double>::quiet_NaN();
  bool self_has_source_bounds = false;
  double self_bound_lower = 0.0;
  double self_bound_upper = std::numeric_limits<double>::infinity();
  if (label_idx >= 0 && label_idx != NA_INTEGER) {
    (void)resolve_label_time_constraint(
        time_constraints, label_idx, self_has_exact_source_time,
        self_exact_source_time, self_has_source_bounds, self_bound_lower,
        self_bound_upper);
  }
  const bool self_is_time_conditioned =
      self_has_exact_source_time || self_has_source_bounds;
  if (label_idx >= 0 && label_idx != NA_INTEGER) {
    if (!self_is_time_conditioned &&
        forced_state_contains_complete(forced_state, label_idx)) {
      return make_node_result(need, 0.0, 0.0, 1.0);
    }
    if (forced_state_contains_survive(forced_state, label_idx)) {
      if (self_is_time_conditioned) {
        return make_node_result(need, 0.0, 0.0, 0.0);
      }
      return make_node_result(need, 0.0, 1.0, 0.0);
    }
  }

  int acc_idx = label_ref.acc_idx;
  if (acc_idx >= 0 && acc_idx < static_cast<int>(ctx.accumulators.size())) {
    const uuber::NativeAccumulator &acc = ctx.accumulators[acc_idx];
    const TrialAccumulatorParams *override =
        get_trial_param_entry(trial_params, acc_idx);
    if (!component_active_idx(acc, component_idx, override)) {
      return make_node_result(need, 0.0, 1.0, 0.0);
    }
    int onset_kind = override ? override->onset_kind : acc.onset_kind;
    int onset_source_acc_idx =
        override ? override->onset_source_acc_idx : acc.onset_source_acc_idx;
    int onset_source_pool_idx =
        override ? override->onset_source_pool_idx : acc.onset_source_pool_idx;
    double onset_lag = override ? override->onset_lag : acc.onset_lag;
    double onset = 0.0;
    double q = 0.0;
    AccDistParams cfg;
    resolve_event_numeric_params(acc, acc_idx, override, trial_params_soa, onset,
                                 q, cfg);
    const LowerBoundTransform lower_bound = default_lower_bound_transform(cfg);
    double density = std::numeric_limits<double>::quiet_NaN();
    double survival = std::numeric_limits<double>::quiet_NaN();
    double cdf = std::numeric_limits<double>::quiet_NaN();
    auto eval_accumulator_base = [&](double eval_t,
                                     EvalNeed inner_need) -> NodeEvalResult {
      double inner_density = std::numeric_limits<double>::quiet_NaN();
      double inner_survival = std::numeric_limits<double>::quiet_NaN();
      double inner_cdf = std::numeric_limits<double>::quiet_NaN();

      if (!ctx.has_chained_onsets || onset_kind == uuber::ONSET_ABSOLUTE) {
        if (needs_density(inner_need)) {
          inner_density = acc_density_from_cfg(eval_t, onset, q, cfg);
        }
        if (needs_survival(inner_need)) {
          inner_survival = acc_survival_from_cfg(eval_t, onset, q, cfg);
        }
        if (needs_cdf(inner_need)) {
          inner_cdf = acc_cdf_success_from_cfg(eval_t, onset, q, cfg);
        }
        return make_node_result(inner_need, inner_density, inner_survival,
                                inner_cdf);
      }

      LabelRef source_ref;
      if (onset_kind == uuber::ONSET_AFTER_ACCUMULATOR) {
        if (onset_source_acc_idx >= 0 &&
            onset_source_acc_idx < static_cast<int>(ctx.accumulators.size())) {
          source_ref.acc_idx = onset_source_acc_idx;
          source_ref.label_id =
              accumulator_label_id_of(ctx, onset_source_acc_idx);
        }
      } else if (onset_kind == uuber::ONSET_AFTER_POOL) {
        if (onset_source_pool_idx >= 0 &&
            onset_source_pool_idx < static_cast<int>(ctx.pools.size())) {
          source_ref.pool_idx = onset_source_pool_idx;
          source_ref.label_id = pool_label_id_of(ctx, onset_source_pool_idx);
        }
      }
      if (source_ref.acc_idx < 0 && source_ref.pool_idx < 0) {
        return make_node_result(inner_need, 0.0, 1.0, 0.0);
      }

      const double lag_total = onset_lag + onset;
      const double x_shift = lag_total + cfg.t0;
      const int source_label_id = source_ref.label_id;

      if (source_label_id >= 0 && source_label_id != NA_INTEGER &&
          state_contains_survive_at(forced_state, time_constraints,
                                    source_label_id, eval_t)) {
        return make_node_result(inner_need, 0.0, 1.0, 0.0);
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
        (void)resolve_label_time_constraint(
            time_constraints, source_label_id, has_exact, exact_time,
            source_has_bounds, source_bound_lower, source_bound_upper);
        if (source_has_bounds) {
          bound_lower = std::max(bound_lower, source_bound_lower);
          bound_upper = std::min(bound_upper, source_bound_upper);
          has_conditioning_bounds = true;
        }
      }

      if (has_exact) {
        double x = eval_t - exact_time - x_shift;
        double cdf_success = clamp_probability(
            eval_cdf_single_with_lower_bound(cfg, x, lower_bound));
        if (needs_density(inner_need)) {
          inner_density = (1.0 - q) *
                          eval_pdf_single_with_lower_bound(cfg, x, lower_bound);
          inner_density = safe_density(inner_density);
        }
        if (needs_cdf(inner_need)) {
          inner_cdf = clamp_probability((1.0 - q) * cdf_success);
        }
        if (needs_survival(inner_need)) {
          inner_survival =
              clamp_probability(q + (1.0 - q) * (1.0 - cdf_success));
        }
        return make_node_result(inner_need, inner_density, inner_survival,
                                inner_cdf);
      }

      if (bound_upper <= bound_lower) {
        return make_node_result(inner_need, 0.0, 1.0, 0.0);
      }
      if (std::isfinite(eval_t) && (eval_t - x_shift <= bound_lower)) {
        return make_node_result(inner_need, 0.0, 1.0, 0.0);
      }

      auto source_eval = [&](double u, EvalNeed source_need) -> NodeEvalResult {
        return eval_event_ref_idx(
            ctx, source_ref, 0u, u, component_idx, source_need, trial_params,
            trial_type_key, false, -1, time_constraints, nullptr, nullptr,
            nullptr, nullptr, trial_params_soa, &forced_state);
      };
      auto source_cdf_at = [&](double u) -> double {
        if (u <= 0.0) {
          return 0.0;
        }
        NodeEvalResult source = source_eval(u, EvalNeed::kCDF);
        return clamp_probability(source.cdf);
      };

      double source_mass = 1.0;
      if (has_conditioning_bounds) {
        double lower_cdf =
            std::isfinite(bound_lower) ? source_cdf_at(bound_lower) : 0.0;
        double upper_cdf = std::isfinite(bound_upper)
                               ? source_cdf_at(bound_upper)
                               : source_cdf_at(
                                     std::numeric_limits<double>::infinity());
        source_mass = upper_cdf - lower_cdf;
        if (!std::isfinite(source_mass) || source_mass <= 0.0) {
          return make_node_result(inner_need, 0.0, 1.0, 0.0);
        }
      }

      double upper_int =
          std::min(bound_upper, std::max(bound_lower, eval_t - x_shift));
      if (upper_int <= bound_lower) {
        return make_node_result(inner_need, 0.0, 1.0, 0.0);
      }

      auto source_pdf = [&](double u) -> double {
        NodeEvalResult source = source_eval(u, EvalNeed::kDensity);
        return safe_density(source.density);
      };

      if (!std::isfinite(eval_t)) {
        double cdf_total = 0.0;
        if (has_conditioning_bounds) {
          cdf_total = clamp_probability(1.0 - q);
        } else {
          NodeEvalResult source_inf =
              source_eval(std::numeric_limits<double>::infinity(),
                          EvalNeed::kCDF);
          double source_finish = clamp_probability(source_inf.cdf);
          cdf_total = clamp_probability((1.0 - q) * source_finish);
        }
        if (needs_density(inner_need)) {
          inner_density = 0.0;
        }
        if (needs_cdf(inner_need)) {
          inner_cdf = cdf_total;
        }
        if (needs_survival(inner_need)) {
          inner_survival = clamp_probability(1.0 - cdf_total);
        }
        return make_node_result(inner_need, inner_density, inner_survival,
                                inner_cdf);
      }

      auto cond_source_pdf = [&](double u) -> double {
        if (!std::isfinite(u) || u <= bound_lower || u > bound_upper) {
          return 0.0;
        }
        double dens = source_pdf(u);
        if (dens <= 0.0) {
          return 0.0;
        }
        if (!has_conditioning_bounds) {
          return dens;
        }
        return dens / source_mass;
      };

      const bool need_density_term = needs_density(inner_need);
      const bool need_cdf_term = needs_cdf(inner_need) || needs_survival(inner_need);
      const FusedIntegralResult fused = integrate_fused_onset_terms(
          cfg, eval_t, x_shift, bound_lower, upper_int, need_density_term,
          need_cdf_term, &lower_bound, cond_source_pdf);

      if (need_density_term) {
        if (fused.ok) {
          inner_density = safe_density((1.0 - q) * fused.density);
        } else {
          auto integrand_density = [&](double u) -> double {
            double fs = cond_source_pdf(u);
            if (fs <= 0.0) {
              return 0.0;
            }
            double x = eval_t - u - x_shift;
            double fx = eval_pdf_single_with_lower_bound(cfg, x, lower_bound);
            if (!std::isfinite(fx) || fx <= 0.0) {
              return 0.0;
            }
            return fs * fx;
          };
          inner_density = (1.0 - q) * uuber::integrate_boost_fn(
                                            integrand_density, bound_lower,
                                            upper_int, kDefaultRelTol,
                                            kDefaultAbsTol, kDefaultMaxDepth);
          inner_density = safe_density(inner_density);
        }
      }

      if (need_cdf_term) {
        double cdf_cond = 0.0;
        if (fused.ok) {
          cdf_cond = fused.cdf;
        } else {
          auto integrand_cdf = [&](double u) -> double {
            double fs = cond_source_pdf(u);
            if (fs <= 0.0) {
              return 0.0;
            }
            double x = eval_t - u - x_shift;
            double Fx = eval_cdf_single_with_lower_bound(cfg, x, lower_bound);
            if (!std::isfinite(Fx) || Fx <= 0.0) {
              return 0.0;
            }
            return fs * clamp_probability(Fx);
          };
          cdf_cond = uuber::integrate_boost_fn(
              integrand_cdf, bound_lower, upper_int, kDefaultRelTol,
              kDefaultAbsTol, kDefaultMaxDepth);
        }
        cdf_cond = clamp_probability(cdf_cond);
        double cdf_total = clamp_probability((1.0 - q) * cdf_cond);
        if (needs_cdf(inner_need)) {
          inner_cdf = cdf_total;
        }
        if (needs_survival(inner_need)) {
          inner_survival = clamp_probability(1.0 - cdf_total);
        }
      }

      return make_node_result(inner_need, inner_density, inner_survival,
                              inner_cdf);
    };

    if (self_is_time_conditioned) {
      auto impossible_self_condition = [&]() -> NodeEvalResult {
        return make_node_result(need, 0.0, 0.0, 0.0);
      };
      auto base_cdf_at = [&](double u) -> double {
        return clamp_probability(
            eval_accumulator_base(u, EvalNeed::kCDF).cdf);
      };
      auto base_density_at = [&](double u) -> double {
        if (!std::isfinite(u)) {
          return 0.0;
        }
        return safe_density(
            eval_accumulator_base(u, EvalNeed::kDensity).density);
      };

      if (self_has_exact_source_time) {
        if (self_has_source_bounds &&
            (!(self_exact_source_time > self_bound_lower) ||
             self_exact_source_time > self_bound_upper)) {
          return impossible_self_condition();
        }
        if (needs_density(need)) {
          density = 0.0;
        }
        if (needs_cdf(need) || needs_survival(need)) {
          const bool before_exact =
              std::isfinite(t) && t < self_exact_source_time;
          const double cdf_exact = before_exact ? 0.0 : 1.0;
          if (needs_cdf(need)) {
            cdf = cdf_exact;
          }
          if (needs_survival(need)) {
            survival = 1.0 - cdf_exact;
          }
        }
      } else {
        const double lower_cdf = std::isfinite(self_bound_lower)
                                     ? base_cdf_at(self_bound_lower)
                                     : 0.0;
        const double upper_cdf =
            std::isfinite(self_bound_upper)
                ? base_cdf_at(self_bound_upper)
                : base_cdf_at(std::numeric_limits<double>::infinity());
        const double condition_mass = upper_cdf - lower_cdf;
        if (self_bound_upper <= self_bound_lower ||
            !std::isfinite(condition_mass) || condition_mass <= 0.0) {
          return impossible_self_condition();
        }

        if (needs_density(need)) {
          if (std::isfinite(t) && t > self_bound_lower && t <= self_bound_upper) {
            density = safe_density(base_density_at(t) / condition_mass);
          } else {
            density = 0.0;
          }
        }
        if (needs_cdf(need) || needs_survival(need)) {
          double cdf_cond = 0.0;
          if (!std::isfinite(t) || t >= self_bound_upper) {
            cdf_cond = 1.0;
          } else if (!(t > self_bound_lower)) {
            cdf_cond = 0.0;
          } else {
            cdf_cond =
                clamp_probability((base_cdf_at(t) - lower_cdf) / condition_mass);
          }
          if (needs_cdf(need)) {
            cdf = cdf_cond;
          }
          if (needs_survival(need)) {
            survival = 1.0 - cdf_cond;
          }
        }
      }
    } else {
      NodeEvalResult base_res = eval_accumulator_base(t, need);
      density = base_res.density;
      survival = base_res.survival;
      cdf = base_res.cdf;
    }

    NodeEvalResult res = make_node_result(need, density, survival, cdf);
    if (needs_density(need) && donor_density > 0.0) {
      res.density = safe_density(res.density + donor_density);
    }
    return res;
  }

  int pool_idx = label_ref.pool_idx;
  if (pool_idx >= 0 && pool_idx < static_cast<int>(ctx.pools.size())) {
    const uuber::NativePool &pool = ctx.pools[pool_idx];
    if (pool.members.empty()) {
      NodeEvalResult res = make_node_result(need, 0.0, 1.0, 0.0);
      if (needs_density(need) && donor_density > 0.0) {
        res.density = safe_density(res.density + donor_density);
      }
      return res;
    }
    PoolEvalScratchGuard scratch_guard;
    PoolEvalScratch &scratch = scratch_guard.scratch();
    scratch.active_indices.clear();
    scratch.active_indices.reserve(pool.members.size());
    for (std::size_t member_idx = 0; member_idx < pool.members.size();
         ++member_idx) {
      const LabelRef &member_ref = pool.member_refs[member_idx];
      int member_acc_idx = member_ref.acc_idx;
      if (member_acc_idx >= 0 &&
          member_acc_idx < static_cast<int>(ctx.accumulators.size())) {
        const uuber::NativeAccumulator &acc = ctx.accumulators[member_acc_idx];
        const TrialAccumulatorParams *override =
            get_trial_param_entry(trial_params, member_acc_idx);
        if (!component_active_idx(acc, component_idx, override))
          continue;
      }
      scratch.active_indices.push_back(member_idx);
    }
    if (scratch.active_indices.empty()) {
      NodeEvalResult res = make_node_result(need, 0.0, 1.0, 0.0);
      if (needs_density(need) && donor_density > 0.0) {
        res.density = safe_density(res.density + donor_density);
      }
      return res;
    }
    const bool need_density = needs_density(need);
    const bool need_survival =
        needs_survival(need) || needs_cdf(need) || need_density;
    EvalNeed child_need = need_density
                              ? (EvalNeed::kDensity | EvalNeed::kSurvival)
                              : EvalNeed::kSurvival;
    if (need_density) {
      scratch.density.resize(scratch.active_indices.size());
    }
    if (need_survival) {
      scratch.survival.resize(scratch.active_indices.size());
    }
    for (std::size_t i = 0; i < scratch.active_indices.size(); ++i) {
      std::size_t member_idx = scratch.active_indices[i];
      const LabelRef &member_ref = pool.member_refs[member_idx];
      NodeEvalResult child = eval_event_ref_idx(
          ctx, member_ref, 0u, t, component_idx, child_need, trial_params,
          std::string(), false, -1, time_constraints, nullptr, nullptr,
          nullptr, nullptr, trial_params_soa, &forced_state);
      if (need_density) {
        scratch.density[i] = child.density;
      }
      if (need_survival) {
        scratch.survival[i] = child.survival;
      }
    }
    double density = std::numeric_limits<double>::quiet_NaN();
    double survival = std::numeric_limits<double>::quiet_NaN();
    double cdf = std::numeric_limits<double>::quiet_NaN();
    if (need_density) {
      density = pool_density_fast(scratch.density, scratch.survival, pool.k);
      if (!std::isfinite(density) || density < 0.0)
        density = 0.0;
    }
    if (needs_survival(need) || needs_cdf(need)) {
      survival = pool_survival_fast(scratch.survival, pool.k);
      if (!std::isfinite(survival))
        survival = 0.0;
    }
    if (needs_cdf(need)) {
      cdf = clamp_probability(1.0 - survival);
    }
    NodeEvalResult res = make_node_result(need, density, survival, cdf);
    if (needs_density(need) && donor_density > 0.0) {
      res.density = safe_density(res.density + donor_density);
    }
    return res;
  }

  return make_node_result(need, 0.0, 1.0, 0.0);
}

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

double accumulate_plan_guess_density(const uuber::NativeContext &ctx,
                                     const Rcpp::List &donors, double t,
                                     int component_idx,
                                     const std::string &trial_type_key,
                                     const TrialParamSet *trial_params,
                                     bool include_na_donors) {
  double total = 0.0;
  for (R_xlen_t i = 0; i < donors.size(); ++i) {
    Rcpp::RObject donor_obj = donors[i];
    if (donor_obj.isNULL())
      continue;
    Rcpp::List donor(donor_obj);
    std::string rt_policy = donor.containsElementNamed("rt_policy")
                                ? Rcpp::as<std::string>(donor["rt_policy"])
                                : std::string("keep");
    bool donor_na = rt_policy == "na";
    if (donor_na != include_na_donors)
      continue;
    double release = donor.containsElementNamed("release")
                         ? Rcpp::as<double>(donor["release"])
                         : 0.0;
    if (!std::isfinite(release))
      release = 0.0;
    if (release <= 0.0)
      continue;
    int donor_node = donor.containsElementNamed("node_id")
                         ? Rcpp::as<int>(donor["node_id"])
                         : NA_INTEGER;
    if (donor_node == NA_INTEGER)
      continue;
    Rcpp::IntegerVector comp_ids =
        donor.containsElementNamed("competitor_ids")
            ? Rcpp::IntegerVector(donor["competitor_ids"])
            : Rcpp::IntegerVector();
    std::vector<int> comp_vec = integer_vector_to_std(comp_ids, false);
    std::vector<int> comp_filtered;
    const std::vector<int> &comp_use =
        filter_competitor_ids(ctx, comp_vec, component_idx, comp_filtered);
    std::string donor_label;
    if (donor.containsElementNamed("source_label") &&
        !Rf_isNull(donor["source_label"])) {
      donor_label = Rcpp::as<std::string>(donor["source_label"]);
    }
    int outcome_idx_context =
        donor_label.empty() ? -1 : outcome_index_of(ctx, donor_label);
    double density = kernel_node_density_entry_idx(
        ctx, donor_node, t, component_idx, nullptr, false, nullptr, false,
        comp_use, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context);
    if (!std::isfinite(density) || density <= 0.0)
      continue;
    total += release * density;
  }
  return total;
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

double guard_cdf_internal(const GuardEvalInput &input, double t,
                          const IntegrationSettings &settings);

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState &state,
                                   EvalNeed need);

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
    if (!ctx.kernel_state_graph.valid ||
        node_idx < 0 ||
        node_idx >= static_cast<int>(
                        ctx.kernel_state_graph.node_guard_transition_idx.size())) {
      Rcpp::stop("IR guard transition metadata missing for guard node_idx=%d",
                 node_idx);
    }
    const int tr_idx = ctx.kernel_state_graph.node_guard_transition_idx
        [static_cast<std::size_t>(node_idx)];
    if (tr_idx < 0 ||
        tr_idx >=
            static_cast<int>(ctx.kernel_state_graph.guard_transitions.size())) {
      Rcpp::stop("IR guard transition mapping missing for guard node_idx=%d",
                 node_idx);
    }
    const uuber::KernelGuardTransition &tr =
        ctx.kernel_state_graph.guard_transitions[static_cast<std::size_t>(
            tr_idx)];
    if (tr.node_idx != node_idx) {
      Rcpp::stop("IR guard transition node mismatch for guard node_idx=%d",
                 node_idx);
    }
    input.guard_transition_idx = tr_idx;
    input.guard_transition = &tr;
  }
  return input;
}

uuber::KernelEventEvalFn make_kernel_event_eval(NodeEvalState &state) {
  return [&](int event_idx) -> uuber::KernelNodeValues {
    uuber::KernelNodeValues out{};
    if (event_idx < 0 || event_idx >= static_cast<int>(state.ctx.ir.events.size())) {
      return out;
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
      node_flags = state.ctx.ir.nodes[static_cast<std::size_t>(event.node_idx)].flags;
    }
    NodeEvalResult event_eval = eval_event_ref_idx(
        state.ctx, ref, node_flags, state.t, state.component_idx, kEvalAll,
        state.trial_params, state.trial_type_key, state.include_na_donors,
        state.outcome_idx, &state.time_constraints,
        nullptr, nullptr, nullptr, nullptr, state.trial_params_soa,
        &state.forced_state);
    out.density = event_eval.density;
    out.survival = event_eval.survival;
    out.cdf = event_eval.cdf;
    return out;
  };
}

uuber::KernelEventBatchEvalFn make_kernel_event_eval_batch(NodeEvalState &state) {
  return [&](int event_idx, const std::vector<double> &times,
             const uuber::KernelEvalNeed &kneed,
             uuber::KernelNodeBatchValues &out) -> bool {
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
    const double saved_t = state.t;
    for (std::size_t i = 0; i < times.size(); ++i) {
      state.t = times[i];
      NodeEvalResult event_eval = eval_event_ref_idx(
          state.ctx, ref, node_flags, state.t, state.component_idx, need,
          state.trial_params, state.trial_type_key, state.include_na_donors,
          state.outcome_idx, &state.time_constraints,
          nullptr, nullptr, nullptr, nullptr, state.trial_params_soa,
          &state.forced_state);
      out.density[i] = safe_density(event_eval.density);
      out.survival[i] = clamp_probability(event_eval.survival);
      out.cdf[i] = clamp_probability(event_eval.cdf);
    }
    state.t = saved_t;
    return true;
  };
}

uuber::KernelGuardEvalFn make_kernel_guard_eval(NodeEvalState &state) {
  return [&](const uuber::KernelOp &op,
             const uuber::KernelNodeValues &reference_value,
             const uuber::KernelNodeValues &blocker_value,
             const uuber::KernelEvalNeed &kneed) -> uuber::KernelNodeValues {
    uuber::KernelNodeValues out{};
    GuardEvalInput guard_input = make_guard_input_forced_state(
        state.ctx, op.node_idx, state.component_idx, &state.trial_type_key,
        state.trial_params, state.trial_params_soa, &state.time_constraints,
        state.forced_state);
    IntegrationSettings settings;
    if (kneed.density) {
      const double ref_density = safe_density(reference_value.density);
      const double blocker_survival = clamp_probability(blocker_value.survival);
      out.density = safe_density(ref_density * blocker_survival);
    }
    if (kneed.cdf || kneed.survival) {
      out.cdf = guard_cdf_internal(guard_input, state.t, settings);
      out.survival = clamp_probability(1.0 - out.cdf);
    }
    return out;
  };
}

uuber::KernelGuardBatchEvalFn make_kernel_guard_eval_batch(NodeEvalState &state) {
  return [&](const uuber::KernelOp &op, const std::vector<double> &times,
             const uuber::KernelNodeBatchValues &reference_values,
             const uuber::KernelNodeBatchValues &blocker_values,
             const uuber::KernelEvalNeed &kneed,
             uuber::KernelNodeBatchValues &out) -> bool {
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
    GuardEvalInput guard_input = make_guard_input_forced_state(
        state.ctx, op.node_idx, state.component_idx, &state.trial_type_key,
        state.trial_params, state.trial_params_soa, &state.time_constraints,
        state.forced_state);
    IntegrationSettings settings;
    for (std::size_t i = 0; i < point_count; ++i) {
      if (kneed.density) {
        const double ref_density = safe_density(reference_values.density[i]);
        const double blocker_survival =
            clamp_probability(blocker_values.survival[i]);
        out.density[i] = safe_density(ref_density * blocker_survival);
      }
      if (kneed.cdf || kneed.survival) {
        out.cdf[i] = guard_cdf_internal(guard_input, times[i], settings);
        out.survival[i] = clamp_probability(1.0 - out.cdf[i]);
      }
    }
    return true;
  };
}

bool eval_node_batch_with_state_dense(
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    EvalNeed need, uuber::KernelBatchRuntimeState &batch_runtime,
    uuber::KernelNodeBatchValues &out_values) {
  out_values = uuber::KernelNodeBatchValues{};
  if (times.empty()) {
    return true;
  }
  uuber::KernelEvalNeed kernel_need;
  kernel_need.density = needs_density(need);
  kernel_need.survival = needs_survival(need);
  kernel_need.cdf = needs_cdf(need);
  if (!kernel_need.density && !kernel_need.survival && !kernel_need.cdf) {
    kernel_need.density = true;
    kernel_need.survival = true;
    kernel_need.cdf = true;
  }
  uuber::KernelEventBatchEvalFn event_eval_batch_cb =
      make_kernel_event_eval_batch(state);
  uuber::KernelGuardBatchEvalFn guard_eval_batch_cb =
      make_kernel_guard_eval_batch(state);
  return uuber::eval_kernel_node_batch_incremental(
      state.ctx.kernel_program, batch_runtime, node_idx, times, kernel_need,
      event_eval_batch_cb, guard_eval_batch_cb, out_values);
}

NodeEvalResult eval_node_recursive_dense(int node_idx, NodeEvalState &state,
                                         EvalNeed need) {
  if (!state.ctx.kernel_program.valid) {
    Rcpp::stop("IR kernel execution unavailable for node evaluation");
  }
  uuber::KernelEvalNeed kernel_need;
  // Evaluate all channels once; composition ops require mixed channels internally.
  kernel_need.density = true;
  kernel_need.survival = true;
  kernel_need.cdf = true;
  uuber::KernelNodeValues kernel_values;
  uuber::KernelEventEvalFn event_eval_cb = make_kernel_event_eval(state);
  uuber::KernelGuardEvalFn guard_eval_cb = make_kernel_guard_eval(state);
  bool kernel_ok = false;
  if (state.kernel_runtime_ready && state.kernel_runtime_ptr) {
    if (!kernel_runtime_cache_safe(state)) {
      uuber::invalidate_kernel_runtime_from_slot(*state.kernel_runtime_ptr, 0);
    }
    kernel_ok = uuber::eval_kernel_node_incremental(
        state.ctx.kernel_program, *state.kernel_runtime_ptr, node_idx, kernel_need,
        event_eval_cb, guard_eval_cb, kernel_values);
  } else {
    kernel_ok = uuber::eval_kernel_node(state.ctx.kernel_program, node_idx,
                                        kernel_need, event_eval_cb, guard_eval_cb,
                                        kernel_values);
  }
  if (!kernel_ok) {
    Rcpp::stop("IR kernel execution failed for node evaluation");
  }
  return make_node_result(need, kernel_values.density, kernel_values.survival,
                          kernel_values.cdf);
}

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState &state,
                                   EvalNeed need) {
  int node_idx = resolve_dense_node_idx_required(state.ctx, node_id);
  return eval_node_recursive_dense(node_idx, state, need);
}

bool eval_node_with_forced_state_view_batch(
    const uuber::NativeContext &ctx, int node_id_or_idx,
    const std::vector<double> &times, int component_idx, EvalNeed need,
    const TrialParamSet *trial_params, const std::string &trial_key,
    const TimeConstraintMap *time_constraints,
    const ForcedStateView &forced_state,
    uuber::KernelNodeBatchValues &out_values) {
  out_values = uuber::KernelNodeBatchValues{};
  if (times.empty()) {
    return true;
  }
  if (!ctx.kernel_program.valid) {
    Rcpp::stop("IR kernel execution unavailable for batch node evaluation");
  }
  const int node_idx = resolve_dense_node_idx_required(ctx, node_id_or_idx);
  NodeEvalState local(ctx, 0.0, component_idx, trial_params, trial_key, false, -1,
                      time_constraints, nullptr, nullptr,
                      false, nullptr, false, nullptr, &forced_state);
  uuber::KernelBatchRuntimeState batch_runtime;
  const bool ok = eval_node_batch_with_state_dense(node_idx, times, local, need,
                                                   batch_runtime, out_values);
  if (!ok) {
    Rcpp::stop("IR batch kernel execution failed for node evaluation");
  }
  return true;
}

inline const uuber::KernelGuardTransition &
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
  if (!input.ctx.kernel_state_graph.valid ||
      input.node_idx < 0 ||
      input.node_idx >= static_cast<int>(
                            input.ctx.kernel_state_graph.node_guard_transition_idx
                                .size())) {
    Rcpp::stop("IR guard transition metadata missing for guard node_idx=%d",
               input.node_idx);
  }
  const int tr_idx = input.ctx.kernel_state_graph.node_guard_transition_idx
      [static_cast<std::size_t>(input.node_idx)];
  if (tr_idx < 0 ||
      tr_idx >=
          static_cast<int>(input.ctx.kernel_state_graph.guard_transitions.size())) {
    Rcpp::stop("IR guard transition mapping missing for guard node_idx=%d",
               input.node_idx);
  }
  const uuber::KernelGuardTransition &tr =
      input.ctx.kernel_state_graph.guard_transitions[static_cast<std::size_t>(
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

double guard_effective_survival_internal(const GuardEvalInput &input, double t,
                                         const IntegrationSettings &settings) {
  if (!std::isfinite(t)) {
    return 1.0;
  }
  if (t <= 0.0) {
    return 1.0;
  }
  (void)settings;
  const uuber::KernelGuardTransition &tr =
      guard_transition_required(input, true);

  NodeEvalState local(input.ctx, t, input.component_idx, input.trial_params,
                      guard_trial_type_key(input), false, -1,
                      input.time_constraints,
                      nullptr, nullptr, false, nullptr, false, nullptr,
                      &input.forced_state);
  NodeEvalResult block =
      eval_node_recursive_dense(tr.blocker_node_idx, local, EvalNeed::kSurvival);
  return clamp_probability(block.survival);
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

bool fast_event_info_dense_idx(const GuardEvalInput &input, int node_idx,
                               FastEventInfo &info) {
  if (time_constraints_any(input.time_constraints)) {
    return false;
  }
  const uuber::IrNode &node = ir_node_required(input.ctx, node_idx);
  if (!(node.op == uuber::IrNodeOp::EventAcc || node.op == uuber::IrNodeOp::EventPool)) {
    return false;
  }
  if ((node.flags & uuber::IR_NODE_FLAG_SPECIAL_DEADLINE) != 0u ||
      (node.flags & uuber::IR_NODE_FLAG_SPECIAL_GUESS) != 0u) {
    return false;
  }
  LabelRef ref = node_label_ref(input.ctx, node);

  const bool has_forced_complete = forced_state_complete_valid(input.forced_state) &&
                                   input.forced_state.forced_complete_bits->any();
  const bool has_forced_survive = forced_state_survive_valid(input.forced_state) &&
                                  input.forced_state.forced_survive_bits->any();
  if (has_forced_complete || has_forced_survive) {
    if (!input.forced_state.label_id_to_bit_idx || ref.label_id == NA_INTEGER ||
        ref.label_id < 0) {
      return false;
    }
    const bool forced_complete_label =
        forced_state_contains_complete(input.forced_state, ref.label_id);
    const bool forced_survive_label =
        forced_state_contains_survive(input.forced_state, ref.label_id);
    if (forced_complete_label || forced_survive_label) {
      return false;
    }
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

inline NodeEvalResult eval_node_from_guard_input(const GuardEvalInput &input,
                                                 int node_idx, double t,
                                                 EvalNeed need) {
  NodeEvalState local(input.ctx, t, input.component_idx, input.trial_params,
                      guard_trial_type_key(input), false, -1,
                      input.time_constraints,
                      nullptr, nullptr, false, nullptr, false, nullptr,
                      &input.forced_state);
  return eval_node_recursive_dense(node_idx, local, need);
}

LinearGuardChain linear_guard_chain_from_transition(
    const uuber::NativeContext &ctx, const uuber::KernelGuardTransition &tr) {
  LinearGuardChain chain;
  if (tr.eval_mode != uuber::KernelGuardEvalMode::LinearChainODE) {
    return chain;
  }
  if (tr.linear_chain_begin < 0 || tr.linear_chain_count <= 0 ||
      tr.linear_chain_begin + tr.linear_chain_count >
          static_cast<int>(ctx.kernel_state_graph.guard_linear_chain_nodes.size()) ||
      tr.linear_chain_leaf_idx < 0 ||
      tr.linear_chain_leaf_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return {};
  }
  const int begin = tr.linear_chain_begin;
  chain.reference_indices =
      &ctx.kernel_state_graph
           .guard_linear_chain_nodes[static_cast<std::size_t>(begin)];
  chain.reference_count = static_cast<std::size_t>(tr.linear_chain_count);
  chain.leaf_blocker_idx = tr.linear_chain_leaf_idx;
  if (chain.reference_count == 0u) {
    return {};
  }
  chain.valid = true;
  return chain;
}

constexpr std::size_t kLinearGuardChainMaxDepth = 8u;
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

template <std::size_t MaxDepth>
class LinearChainIntegrator final {
public:
  using ChainState = std::array<double, MaxDepth>;
  static constexpr std::size_t kMaxGridPoints =
      static_cast<std::size_t>(2 * kLinearChainMaxSteps + 1);
  static constexpr std::size_t kMaxValueWidth = MaxDepth + 1u;

  LinearChainIntegrator(const GuardEvalInput &input,
                        const LinearGuardChain &chain,
                        double t,
                        const CanonicalIntegrationSettings &canonical)
      : input_(input), chain_(chain), t_(t), canonical_(canonical),
        depth_(chain.reference_count), value_width_(depth_ + 1u),
        ref_infos_(depth_) {
    ref_info_ptrs_.fill(nullptr);
    ref_density_ok_.fill(false);
    for (std::size_t i = 0; i < depth_; ++i) {
      if (fast_event_info_dense_idx(input_, chain_.reference_indices[i],
                                    ref_infos_[i])) {
        ref_info_ptrs_[i] = &ref_infos_[i];
        ref_density_ok_[i] = fast_event_density_supported(input_, ref_infos_[i]);
      }
    }
    leaf_info_ptr_ =
        fast_event_info_dense_idx(input_, chain_.leaf_blocker_idx, leaf_info_)
            ? &leaf_info_
            : nullptr;
  }

  bool solve(double &out_cdf) {
    if (!is_valid()) {
      return false;
    }
    const int n_steps = linear_chain_step_budget(depth_, canonical_);
    return solve_fixed_rk4(n_steps, out_cdf);
  }

private:
  bool is_valid() const {
    return chain_.valid && depth_ >= 1u && depth_ <= MaxDepth &&
           chain_.leaf_blocker_idx >= 0 && std::isfinite(t_) && t_ > 0.0;
  }

  inline double grid_x_value(std::size_t grid_idx, double h) const {
    return 0.5 * h * static_cast<double>(grid_idx);
  }

  inline double *grid_slot(std::size_t grid_idx) {
    return grid_vals_storage_.data() + (grid_idx * value_width_);
  }

  bool precompute_uniform_grid_values(int n_steps, double h) {
    if (n_steps < 1 || !std::isfinite(h) || h <= 0.0) {
      return false;
    }
    const std::size_t grid_points = static_cast<std::size_t>(2 * n_steps + 1);
    if (grid_points < 3u) {
      return false;
    }
    if (grid_points > kMaxGridPoints || value_width_ > kMaxValueWidth) {
      return false;
    }
    grid_points_ = grid_points;
    std::fill_n(grid_vals_storage_.data(), grid_points_ * value_width_, 0.0);
    double *batch_shifted = batch_shifted_storage_.data();
    double *batch_values = batch_values_storage_.data();

    for (std::size_t i = 0; i < depth_; ++i) {
      const FastEventInfo *info_ptr = ref_info_ptrs_[i];
      if (info_ptr && ref_density_ok_[i]) {
        const FastEventInfo &info = *info_ptr;
        const double onset_eff = total_onset_with_t0(info.onset, info.cfg);
        const double success_prob = clamp(1.0 - info.q, 0.0, 1.0);
        for (std::size_t g = 0; g < grid_points; ++g) {
          batch_shifted[g] = grid_x_value(g, h) - onset_eff;
        }
        eval_pdf_vec_with_lower_bound(info.cfg, info.lower_bound, batch_shifted,
                                      grid_points, batch_values);
        for (std::size_t g = 0; g < grid_points; ++g) {
          double dens = 0.0;
          const double x = grid_x_value(g, h);
          if (info.component_ok && success_prob > 0.0 && x >= onset_eff) {
            dens = success_prob * safe_density(batch_values[g]);
          }
          double *vals = grid_slot(g);
          vals[i] = dens;
        }
      } else {
        for (std::size_t g = 0; g < grid_points; ++g) {
          const double x = grid_x_value(g, h);
          NodeEvalResult ref_eval =
              eval_node_from_guard_input(input_, chain_.reference_indices[i], x,
                                         EvalNeed::kDensity);
          double dens = safe_density(ref_eval.density);
          if (!std::isfinite(dens) || dens <= 0.0) {
            dens = 0.0;
          }
          double *vals = grid_slot(g);
          vals[i] = dens;
        }
      }
    }

    if (leaf_info_ptr_) {
      const FastEventInfo &info = *leaf_info_ptr_;
      const double onset_eff = total_onset_with_t0(info.onset, info.cfg);
      const double success_prob = clamp(1.0 - info.q, 0.0, 1.0);
      for (std::size_t g = 0; g < grid_points; ++g) {
        batch_shifted[g] = grid_x_value(g, h) - onset_eff;
      }
      eval_cdf_vec_with_lower_bound(info.cfg, info.lower_bound, batch_shifted,
                                    grid_points, batch_values);
      for (std::size_t g = 0; g < grid_points; ++g) {
        double surv = 1.0;
        const double x = grid_x_value(g, h);
        if (!info.component_ok) {
          surv = 1.0;
        } else if (x >= onset_eff) {
          const double cdf = clamp_probability(batch_values[g]);
          surv = info.q + success_prob * (1.0 - cdf);
        }
        double *vals = grid_slot(g);
        vals[depth_] = clamp_probability(surv);
      }
    } else {
      for (std::size_t g = 0; g < grid_points; ++g) {
        const double x = grid_x_value(g, h);
        NodeEvalResult leaf_eval =
            eval_node_from_guard_input(input_, chain_.leaf_blocker_idx, x,
                                       EvalNeed::kSurvival);
        const double leaf_surv = clamp_probability(leaf_eval.survival);
        double *vals = grid_slot(g);
        vals[depth_] = clamp_probability(leaf_surv);
      }
    }
    return true;
  }

  void deriv_from_vals(const double *vals, const ChainState &state, ChainState &dst) {
    dst.fill(0.0);
    for (std::size_t i = depth_; i-- > 0;) {
      const double f_ref = vals[i];
      if (!std::isfinite(f_ref) || f_ref <= 0.0) {
        continue;
      }
      const double s_down = (i + 1 < depth_)
                                ? clamp_probability(1.0 - state[i + 1])
                                : vals[depth_];
      if (!std::isfinite(s_down) || s_down <= 0.0) {
        continue;
      }
      const double val = f_ref * s_down;
      if (std::isfinite(val) && val > 0.0) {
        dst[i] = val;
      }
    }
  }

  void rk4_step(const ChainState &state, double h, const double *vals_x,
                const double *vals_half, const double *vals_next,
                ChainState &out) {
    out.fill(0.0);
    deriv_from_vals(vals_x, state, k1_);
    for (std::size_t i = 0; i < depth_; ++i) {
      tmp_[i] = state[i] + 0.5 * h * k1_[i];
    }
    deriv_from_vals(vals_half, tmp_, k2_);
    for (std::size_t i = 0; i < depth_; ++i) {
      tmp_[i] = state[i] + 0.5 * h * k2_[i];
    }
    deriv_from_vals(vals_half, tmp_, k3_);
    for (std::size_t i = 0; i < depth_; ++i) {
      tmp_[i] = state[i] + h * k3_[i];
    }
    deriv_from_vals(vals_next, tmp_, k4_);
    for (std::size_t i = 0; i < depth_; ++i) {
      out[i] = state[i] + (h / 6.0) *
                            (k1_[i] + 2.0 * k2_[i] + 2.0 * k3_[i] + k4_[i]);
    }
  }

  bool finalize_state(const ChainState &state, double &out_val) const {
    const double cdf = state.front();
    if (!std::isfinite(cdf)) {
      return false;
    }
    out_val = clamp_probability(cdf);
    return std::isfinite(out_val);
  }

  bool solve_fixed_rk4(int n_steps, double &out_val) {
    if (n_steps < 1 || !std::isfinite(t_) || t_ <= 0.0) {
      return false;
    }
    const double h = t_ / static_cast<double>(n_steps);
    if (!std::isfinite(h) || h <= 0.0) {
      return false;
    }
    if (!precompute_uniform_grid_values(n_steps, h)) {
      return false;
    }
    ChainState state{};
    ChainState next{};
    for (int i = 0; i < n_steps; ++i) {
      const std::size_t g0 = static_cast<std::size_t>(2 * i);
      const double *vals_x = grid_slot(g0);
      const double *vals_half = grid_slot(g0 + 1u);
      const double *vals_next = grid_slot(g0 + 2u);
      rk4_step(state, h, vals_x, vals_half, vals_next, next);
      state = next;
      for (std::size_t j = 0; j < depth_; ++j) {
        if (!std::isfinite(state[j])) {
          return false;
        }
      }
    }
    return finalize_state(state, out_val);
  }

  const GuardEvalInput &input_;
  const LinearGuardChain &chain_;
  double t_{0.0};
  CanonicalIntegrationSettings canonical_;
  std::size_t depth_{0u};
  std::size_t value_width_{0u};

  std::vector<FastEventInfo> ref_infos_;
  std::array<const FastEventInfo *, MaxDepth> ref_info_ptrs_{};
  std::array<bool, MaxDepth> ref_density_ok_{};
  FastEventInfo leaf_info_{};
  const FastEventInfo *leaf_info_ptr_{nullptr};

  std::size_t grid_points_{0u};
  std::array<double, kMaxGridPoints * kMaxValueWidth> grid_vals_storage_;
  std::array<double, kMaxGridPoints> batch_shifted_storage_;
  std::array<double, kMaxGridPoints> batch_values_storage_;

  ChainState k1_{};
  ChainState k2_{};
  ChainState k3_{};
  ChainState k4_{};
  ChainState tmp_{};
};

bool eval_optimized_linear_guard_chain_ode(
    const GuardEvalInput &input, const LinearGuardChain &chain, double t,
    const CanonicalIntegrationSettings &canonical, double &out_cdf) {
  const std::size_t depth = chain.reference_count;
  if (!chain.valid || depth < 1u || chain.leaf_blocker_idx < 0 ||
      !std::isfinite(t) || t <= 0.0) {
    return false;
  }

  if (depth <= kLinearGuardChainMaxDepth) {
    LinearChainIntegrator<kLinearGuardChainMaxDepth> integrator(input, chain, t,
                                                                canonical);
    return integrator.solve(out_cdf);
  }

  class DynamicLinearChainIntegrator final {
  public:
    DynamicLinearChainIntegrator(const GuardEvalInput &input,
                                 const LinearGuardChain &chain, double t,
                                 const CanonicalIntegrationSettings &canonical)
        : input_(input), chain_(chain), t_(t), canonical_(canonical),
          depth_(chain.reference_count), value_width_(depth_ + 1u),
          ref_infos_(depth_), ref_info_ptrs_(depth_, nullptr),
          ref_density_ok_(depth_, false), state_(depth_, 0.0),
          next_(depth_, 0.0), k1_(depth_, 0.0), k2_(depth_, 0.0),
          k3_(depth_, 0.0), k4_(depth_, 0.0), tmp_(depth_, 0.0) {
      for (std::size_t i = 0; i < depth_; ++i) {
        if (fast_event_info_dense_idx(input_, chain_.reference_indices[i],
                                      ref_infos_[i])) {
          ref_info_ptrs_[i] = &ref_infos_[i];
          ref_density_ok_[i] =
              fast_event_density_supported(input_, ref_infos_[i]);
        }
      }
      leaf_info_ptr_ =
          fast_event_info_dense_idx(input_, chain_.leaf_blocker_idx, leaf_info_)
              ? &leaf_info_
              : nullptr;
    }

    bool solve(double &out_cdf) {
      if (!is_valid()) {
        return false;
      }
      const int n_steps = linear_chain_step_budget(depth_, canonical_);
      return solve_fixed_rk4(n_steps, out_cdf);
    }

  private:
    bool is_valid() const {
      return chain_.valid && depth_ >= 1u && chain_.leaf_blocker_idx >= 0 &&
             std::isfinite(t_) && t_ > 0.0;
    }

    inline double grid_x_value(std::size_t grid_idx, double h) const {
      return 0.5 * h * static_cast<double>(grid_idx);
    }

    inline double *grid_slot(std::size_t grid_idx) {
      return grid_vals_.data() + (grid_idx * value_width_);
    }

    bool precompute_uniform_grid_values(int n_steps, double h) {
      if (n_steps < 1 || !std::isfinite(h) || h <= 0.0) {
        return false;
      }
      const std::size_t grid_points = static_cast<std::size_t>(2 * n_steps + 1);
      if (grid_points < 3u) {
        return false;
      }
      grid_vals_.assign(grid_points * value_width_, 0.0);
      batch_shifted_.assign(grid_points, 0.0);
      batch_values_.assign(grid_points, 0.0);

      for (std::size_t i = 0; i < depth_; ++i) {
        const FastEventInfo *info_ptr = ref_info_ptrs_[i];
        if (info_ptr && ref_density_ok_[i]) {
          const FastEventInfo &info = *info_ptr;
          const double onset_eff = total_onset_with_t0(info.onset, info.cfg);
          const double success_prob = clamp(1.0 - info.q, 0.0, 1.0);
          for (std::size_t g = 0; g < grid_points; ++g) {
            batch_shifted_[g] = grid_x_value(g, h) - onset_eff;
          }
          eval_pdf_vec_with_lower_bound(
              info.cfg, info.lower_bound, batch_shifted_.data(), grid_points,
              batch_values_.data());
          for (std::size_t g = 0; g < grid_points; ++g) {
            double dens = 0.0;
            const double x = grid_x_value(g, h);
            if (info.component_ok && success_prob > 0.0 && x >= onset_eff) {
              dens = success_prob * safe_density(batch_values_[g]);
            }
            double *vals = grid_slot(g);
            vals[i] = dens;
          }
        } else {
          for (std::size_t g = 0; g < grid_points; ++g) {
            const double x = grid_x_value(g, h);
            NodeEvalResult ref_eval = eval_node_from_guard_input(
                input_, chain_.reference_indices[i], x, EvalNeed::kDensity);
            double dens = safe_density(ref_eval.density);
            if (!std::isfinite(dens) || dens <= 0.0) {
              dens = 0.0;
            }
            double *vals = grid_slot(g);
            vals[i] = dens;
          }
        }
      }

      if (leaf_info_ptr_) {
        const FastEventInfo &info = *leaf_info_ptr_;
        const double onset_eff = total_onset_with_t0(info.onset, info.cfg);
        const double success_prob = clamp(1.0 - info.q, 0.0, 1.0);
        for (std::size_t g = 0; g < grid_points; ++g) {
          batch_shifted_[g] = grid_x_value(g, h) - onset_eff;
        }
        eval_cdf_vec_with_lower_bound(
            info.cfg, info.lower_bound, batch_shifted_.data(), grid_points,
            batch_values_.data());
        for (std::size_t g = 0; g < grid_points; ++g) {
          double surv = 1.0;
          const double x = grid_x_value(g, h);
          if (!info.component_ok) {
            surv = 1.0;
          } else if (x >= onset_eff) {
            const double cdf = clamp_probability(batch_values_[g]);
            surv = info.q + success_prob * (1.0 - cdf);
          }
          double *vals = grid_slot(g);
          vals[depth_] = clamp_probability(surv);
        }
      } else {
        for (std::size_t g = 0; g < grid_points; ++g) {
          const double x = grid_x_value(g, h);
          NodeEvalResult leaf_eval = eval_node_from_guard_input(
              input_, chain_.leaf_blocker_idx, x, EvalNeed::kSurvival);
          const double leaf_surv = clamp_probability(leaf_eval.survival);
          double *vals = grid_slot(g);
          vals[depth_] = clamp_probability(leaf_surv);
        }
      }
      return true;
    }

    void deriv_from_vals(const double *vals, const std::vector<double> &state,
                         std::vector<double> &dst) {
      std::fill(dst.begin(), dst.end(), 0.0);
      for (std::size_t i = depth_; i-- > 0;) {
        const double f_ref = vals[i];
        if (!std::isfinite(f_ref) || f_ref <= 0.0) {
          continue;
        }
        const double s_down = (i + 1 < depth_)
                                  ? clamp_probability(1.0 - state[i + 1])
                                  : vals[depth_];
        if (!std::isfinite(s_down) || s_down <= 0.0) {
          continue;
        }
        const double val = f_ref * s_down;
        if (std::isfinite(val) && val > 0.0) {
          dst[i] = val;
        }
      }
    }

    void rk4_step(const std::vector<double> &state, double h,
                  const double *vals_x, const double *vals_half,
                  const double *vals_next, std::vector<double> &out) {
      std::fill(out.begin(), out.end(), 0.0);
      deriv_from_vals(vals_x, state, k1_);
      for (std::size_t i = 0; i < depth_; ++i) {
        tmp_[i] = state[i] + 0.5 * h * k1_[i];
      }
      deriv_from_vals(vals_half, tmp_, k2_);
      for (std::size_t i = 0; i < depth_; ++i) {
        tmp_[i] = state[i] + 0.5 * h * k2_[i];
      }
      deriv_from_vals(vals_half, tmp_, k3_);
      for (std::size_t i = 0; i < depth_; ++i) {
        tmp_[i] = state[i] + h * k3_[i];
      }
      deriv_from_vals(vals_next, tmp_, k4_);
      for (std::size_t i = 0; i < depth_; ++i) {
        out[i] = state[i] + (h / 6.0) *
                              (k1_[i] + 2.0 * k2_[i] + 2.0 * k3_[i] + k4_[i]);
      }
    }

    bool solve_fixed_rk4(int n_steps, double &out_val) {
      if (n_steps < 1 || !std::isfinite(t_) || t_ <= 0.0) {
        return false;
      }
      const double h = t_ / static_cast<double>(n_steps);
      if (!std::isfinite(h) || h <= 0.0) {
        return false;
      }
      if (!precompute_uniform_grid_values(n_steps, h)) {
        return false;
      }
      std::fill(state_.begin(), state_.end(), 0.0);
      std::fill(next_.begin(), next_.end(), 0.0);
      for (int i = 0; i < n_steps; ++i) {
        const std::size_t g0 = static_cast<std::size_t>(2 * i);
        const double *vals_x = grid_slot(g0);
        const double *vals_half = grid_slot(g0 + 1u);
        const double *vals_next = grid_slot(g0 + 2u);
        rk4_step(state_, h, vals_x, vals_half, vals_next, next_);
        state_.swap(next_);
        for (std::size_t j = 0; j < depth_; ++j) {
          if (!std::isfinite(state_[j])) {
            return false;
          }
        }
      }
      const double cdf = state_.front();
      if (!std::isfinite(cdf)) {
        return false;
      }
      out_val = clamp_probability(cdf);
      return std::isfinite(out_val);
    }

    const GuardEvalInput &input_;
    const LinearGuardChain &chain_;
    double t_{0.0};
    CanonicalIntegrationSettings canonical_;
    std::size_t depth_{0u};
    std::size_t value_width_{0u};

    std::vector<FastEventInfo> ref_infos_;
    std::vector<const FastEventInfo *> ref_info_ptrs_;
    std::vector<bool> ref_density_ok_;
    FastEventInfo leaf_info_{};
    const FastEventInfo *leaf_info_ptr_{nullptr};

    std::vector<double> grid_vals_;
    std::vector<double> batch_shifted_;
    std::vector<double> batch_values_;

    std::vector<double> state_;
    std::vector<double> next_;
    std::vector<double> k1_;
    std::vector<double> k2_;
    std::vector<double> k3_;
    std::vector<double> k4_;
    std::vector<double> tmp_;
  };

  DynamicLinearChainIntegrator integrator(input, chain, t, canonical);
  return integrator.solve(out_cdf);
}

double guard_cdf_internal(const GuardEvalInput &input, double t,
                          const IntegrationSettings &settings) {
  if (!std::isfinite(t)) {
    return 1.0;
  }
  if (t <= 0.0) {
    return 0.0;
  }
  const CanonicalIntegrationSettings canonical =
      canonicalize_integration_settings(settings);

  const uuber::KernelGuardTransition &tr = guard_transition_required(input);
  if (tr.eval_mode != uuber::KernelGuardEvalMode::LinearChainODE) {
    Rcpp::stop("IR guard transition mode must be LinearChainODE for node_idx=%d",
               input.node_idx);
  }
  LinearGuardChain chain = linear_guard_chain_from_transition(input.ctx, tr);
  if (!chain.valid) {
    Rcpp::stop("IR guard linear-chain metadata invalid for guard node_idx=%d",
               input.node_idx);
  }
  double cdf_chain = 0.0;
  if (!eval_optimized_linear_guard_chain_ode(input, chain, t, canonical,
                                             cdf_chain)) {
    Rcpp::stop("IR guard linear-chain ODE evaluation failed for node_idx=%d",
               input.node_idx);
  }
  return cdf_chain;
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

double
evaluate_survival_with_forced(int node_id,
                              const uuber::BitsetState *forced_complete_bits,
                              bool forced_complete_bits_valid,
                              const uuber::BitsetState *forced_survive_bits,
                              bool forced_survive_bits_valid,
                              int component_idx, double t,
                              const uuber::NativeContext &ctx,
                              const std::string &trial_key = std::string(),
                              const TrialParamSet *trial_params = nullptr,
                              const TimeConstraintMap *time_constraints =
                                  nullptr,
                              uuber::KernelRuntimeState *kernel_runtime = nullptr) {
  const bool kernel_runtime_usable =
      kernel_runtime &&
      !forced_bits_any(forced_complete_bits, forced_complete_bits_valid) &&
      !forced_bits_any(forced_survive_bits, forced_survive_bits_valid) &&
      !time_constraints_any(time_constraints);
  NodeEvalState state(ctx, t, component_idx, trial_params, trial_key, false,
                      -1,
                      time_constraints, nullptr, forced_complete_bits,
                      forced_complete_bits_valid, forced_survive_bits,
                      forced_survive_bits_valid,
                      kernel_runtime_usable ? kernel_runtime : nullptr);
  NodeEvalResult res = eval_node_recursive(node_id, state, EvalNeed::kSurvival);
  return clamp_probability(res.survival);
}

inline double node_density_with_competitors_from_state_dense(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr,
    const OutcomeCouplingProgram *coupling_program = nullptr) {
  if (should_use_exact_density_program(ctx, node_idx, competitor_ids, state,
                                       coupling_program)) {
    return exact_outcome_density_from_state(ctx, node_idx, competitor_ids, state,
                                            competitor_cache);
  }
  NodeEvalResult base =
      eval_node_recursive_dense(node_idx, state, EvalNeed::kDensity);
  double density = base.density;
  if (!std::isfinite(density) || density <= 0.0) {
    return 0.0;
  }
  if (!competitor_ids.empty()) {
    const double surv =
        competitor_survival_from_state_compiled_ops(ctx, competitor_ids, state,
                                                    competitor_cache);
    if (!std::isfinite(surv) || surv <= 0.0) {
      return 0.0;
    }
    density *= surv;
  }
  return density;
}

inline double node_density_with_competitors_from_state(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    uuber::KernelRuntimeState *kernel_runtime,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr,
    const OutcomeCouplingProgram *coupling_program = nullptr) {
  (void)kernel_runtime;
  const int node_idx = resolve_dense_node_idx_required(ctx, node_id);
  return node_density_with_competitors_from_state_dense(
      ctx, node_idx, competitor_ids, state, competitor_cache, coupling_program);
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
    uuber::KernelBatchRuntimeState *kernel_batch_runtime) {
  density_out.assign(times.size(), 0.0);
  if (times.empty()) {
    return true;
  }

  const bool kernel_runtime_usable =
      kernel_batch_runtime &&
      !forced_bits_any(forced_complete_bits, forced_complete_bits_valid) &&
      !forced_bits_any(forced_survive_bits, forced_survive_bits_valid) &&
      !time_constraints_any(time_constraints);
  NodeEvalState state(
      ctx, 0.0, component_idx, trial_params, trial_type_key, include_na_donors,
      outcome_idx_context, time_constraints, nullptr,
      forced_complete_bits_valid ? forced_complete_bits : nullptr,
      forced_complete_bits_valid,
      forced_survive_bits_valid ? forced_survive_bits : nullptr,
      forced_survive_bits_valid, nullptr);

  const uuber::BitsetState forced_complete_seed = state.forced_complete_bits;
  const uuber::BitsetState forced_survive_seed = state.forced_survive_bits;
  const bool forced_complete_seed_valid = state.forced_complete_bits_valid;
  const bool forced_survive_seed_valid = state.forced_survive_bits_valid;

  const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr;
  if (!competitor_ids.empty()) {
    competitor_cache = &fetch_competitor_cluster_cache(ctx, competitor_ids);
  }
  const int node_idx = resolve_dense_node_idx_required(ctx, node_id);
  OutcomeCouplingProgram coupling_program =
      resolve_density_coupling_program(ctx, node_idx, competitor_ids);
  if (should_use_exact_density_program(ctx, node_idx, competitor_ids, state,
                                       &coupling_program)) {
    return exact_outcome_density_batch_from_state(
        ctx, node_idx, competitor_ids, state, times, density_out,
        competitor_cache);
  }
  const bool has_competitors =
      competitor_cache && !competitor_cache->compiled_ops.empty();
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
  uuber::KernelBatchRuntimeState local_batch_runtime;
  uuber::KernelBatchRuntimeState *batch_runtime_ptr =
      kernel_runtime_usable ? kernel_batch_runtime : &local_batch_runtime;

  state.include_na_donors = include_na_donors;
  state.outcome_idx = outcome_idx_context;
  state.forced_complete_bits = forced_complete_seed;
  state.forced_survive_bits = forced_survive_seed;
  state.forced_complete_bits_valid = forced_complete_seed_valid;
  state.forced_survive_bits_valid = forced_survive_seed_valid;
  uuber::invalidate_kernel_batch_runtime_from_slot(*batch_runtime_ptr, 0);
  uuber::KernelNodeBatchValues base_values;
  if (!eval_node_batch_with_state_dense(node_idx, eval_times, state,
                                        EvalNeed::kDensity, *batch_runtime_ptr,
                                        base_values) ||
      base_values.density.size() != times.size()) {
    Rcpp::stop("IR batch kernel execution failed for node density");
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
    uuber::invalidate_kernel_batch_runtime_from_slot(*batch_runtime_ptr, 0);
    uuber::KernelEventBatchEvalFn event_eval_batch_cb =
        make_kernel_event_eval_batch(state);
    uuber::KernelGuardBatchEvalFn guard_eval_batch_cb =
        make_kernel_guard_eval_batch(state);
    uuber::KernelBatchTransitionApplyFn apply_transition_batch =
        [&](const uuber::CompetitorCompiledOp &op,
            uuber::KernelBatchRuntimeState &runtime) {
          if (op.transition_mask_begin < 0 || op.transition_mask_count <= 0) {
            return;
          }
          apply_transition_mask_words(
              state.ctx, op.transition_mask_begin, op.transition_mask_count,
              op.transition_invalidate_slot, state.forced_survive_bits,
              state.forced_survive_bits_valid,
              (state.kernel_runtime_ready && state.kernel_runtime_ptr)
                  ? state.kernel_runtime_ptr
                  : nullptr,
              &runtime);
        };
    std::vector<double> survival_product;
    if (!uuber::eval_kernel_competitor_product_batch_incremental(
            state.ctx.kernel_program, *batch_runtime_ptr,
            competitor_cache->compiled_ops, eval_times, event_eval_batch_cb,
            guard_eval_batch_cb, apply_transition_batch, survival_product) ||
        survival_product.size() != times.size()) {
      Rcpp::stop("IR batch kernel execution failed for competitor product");
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
    uuber::KernelBatchRuntimeState *kernel_batch_runtime) {
  if (!use_shared_trigger_eval || !trial_params) {
    return node_density_with_competitors_batch_internal(
        ctx, node_id, times, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, density_out, time_constraints,
        kernel_batch_runtime);
  }

  SharedTriggerPlan local_plan;
  const SharedTriggerPlan *plan_ptr = trigger_plan;
  if (!plan_ptr) {
    local_plan = build_shared_trigger_plan(ctx, trial_params);
    plan_ptr = &local_plan;
  }
  if (shared_trigger_count(*plan_ptr) == 0u) {
    return node_density_with_competitors_batch_internal(
        ctx, node_id, times, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, density_out, time_constraints,
        kernel_batch_runtime);
  }

  return false;
}

double node_density_with_competitors_internal(
    const uuber::NativeContext &ctx, int node_id, double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const TimeConstraintMap *time_constraints,
    uuber::KernelRuntimeState *kernel_runtime) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  NodeEvalState state(
      ctx, t, component_idx, trial_params, trial_type_key, include_na_donors,
      outcome_idx_context, time_constraints, nullptr,
      forced_complete_bits_valid ? forced_complete_bits : nullptr,
      forced_complete_bits_valid,
      forced_survive_bits_valid ? forced_survive_bits : nullptr,
      forced_survive_bits_valid, kernel_runtime);
  const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr;
  if (!competitor_ids.empty()) {
    competitor_cache = &fetch_competitor_cluster_cache(ctx, competitor_ids);
  }
  OutcomeCouplingProgram coupling_program =
      resolve_density_coupling_program(ctx, node_id, competitor_ids);
  return node_density_with_competitors_from_state(
      ctx, node_id, competitor_ids, state, kernel_runtime, competitor_cache,
      &coupling_program);
}

double kernel_node_density_entry_idx(
    const uuber::NativeContext &ctx, int node_id, double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const TimeConstraintMap *time_constraints,
    uuber::KernelRuntimeState *kernel_runtime) {
  record_unified_outcome_kernel_density_call();
  return node_density_with_competitors_internal(
      ctx, node_id, t, component_idx, forced_complete_bits,
      forced_complete_bits_valid, forced_survive_bits,
      forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
      include_na_donors, outcome_idx_context, time_constraints, kernel_runtime);
}

inline double node_density_entry_idx(
    const uuber::NativeContext &ctx, int node_id, double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan, TrialParamSet *scratch_params,
    bool use_shared_trigger_eval,
    const TimeConstraintMap *time_constraints,
    uuber::KernelRuntimeState *kernel_runtime) {
  if (!use_shared_trigger_eval || !trial_params) {
    return kernel_node_density_entry_idx(
        ctx, node_id, t, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, time_constraints,
        kernel_runtime);
  }
  SharedTriggerPlan local_plan;
  const SharedTriggerPlan *plan_ptr = trigger_plan;
  if (!plan_ptr) {
    local_plan = build_shared_trigger_plan(ctx, trial_params);
    plan_ptr = &local_plan;
  }
  if (shared_trigger_count(*plan_ptr) == 0) {
    return kernel_node_density_entry_idx(
        ctx, node_id, t, component_idx, forced_complete_bits,
        forced_complete_bits_valid, forced_survive_bits,
        forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_idx_context, time_constraints,
        kernel_runtime);
  }
  TrialParamSet local_scratch;
  TrialParamSet &scratch = scratch_params ? *scratch_params : local_scratch;
  prepare_trigger_scratch(ctx, trial_params, *plan_ptr, scratch);
  const bool forced_empty =
      !forced_bits_any(forced_complete_bits, forced_complete_bits_valid) &&
      !forced_bits_any(forced_survive_bits, forced_survive_bits_valid);
  uuber::KernelRuntimeState shared_kernel_runtime;
  uuber::KernelRuntimeState *shared_kernel_runtime_ptr = nullptr;
  if (ctx.kernel_program.valid && shared_trigger_count(*plan_ptr) > 0 &&
      forced_empty) {
    uuber::reset_kernel_runtime(ctx.kernel_program, shared_kernel_runtime);
    shared_kernel_runtime_ptr = &shared_kernel_runtime;
  }
  NodeEvalState prebuilt_state(
      ctx, t, component_idx, &scratch, trial_type_key, include_na_donors,
      outcome_idx_context, time_constraints, nullptr,
      forced_complete_bits_valid ? forced_complete_bits : nullptr,
      forced_complete_bits_valid,
      forced_survive_bits_valid ? forced_survive_bits : nullptr,
      forced_survive_bits_valid, shared_kernel_runtime_ptr);
  const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr;
  if (!competitor_ids.empty()) {
    competitor_cache = &fetch_competitor_cluster_cache(ctx, competitor_ids);
  }
  OutcomeCouplingProgram coupling_program =
      resolve_density_coupling_program(ctx, node_id, competitor_ids);
  double total = evaluate_shared_trigger_states(
      ctx, *plan_ptr, scratch, "Too many shared triggers for density evaluation",
      [&]() -> double {
        return node_density_with_competitors_from_state(
            ctx, node_id, competitor_ids, prebuilt_state,
            shared_kernel_runtime_ptr, competitor_cache, &coupling_program);
      },
      shared_kernel_runtime_ptr);
  if (!std::isfinite(total) || total <= 0.0)
    return 0.0;
  return total;
}

double native_outcome_probability_bits_impl_idx(
    const uuber::NativeContext &ctx, int node_id, double upper, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<int> &competitor_ids_raw, double rel_tol, double abs_tol,
    int max_depth, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan,
    TrialParamSet *trigger_scratch,
    const OutcomeCouplingProgram *precompiled_coupling_program) {
  if (upper <= 0.0) {
    return 0.0;
  }
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.components.ids.size())) {
    component_idx = -1;
  }
  const bool forced_empty =
      !forced_bits_any(forced_complete_bits, forced_complete_bits_valid) &&
      !forced_bits_any(forced_survive_bits, forced_survive_bits_valid);
  std::vector<int> comp_vec_filtered;
  const std::vector<int> &comp_vec = filter_competitor_ids(
      ctx, competitor_ids_raw, component_idx, comp_vec_filtered);
  const OutcomeCouplingProgram *coupling_program_ptr = nullptr;
  OutcomeCouplingProgram resolved_coupling_program;
  if (forced_empty && precompiled_coupling_program &&
      precompiled_coupling_program->valid) {
    coupling_program_ptr = precompiled_coupling_program;
  } else {
    resolved_coupling_program = resolve_outcome_coupling_program_with_generic(
        ctx, node_id, comp_vec, forced_empty);
    if (resolved_coupling_program.valid) {
      coupling_program_ptr = &resolved_coupling_program;
    }
  }
  if (!coupling_program_ptr || !coupling_program_ptr->valid) {
    return 0.0;
  }
  auto integrate_with_params = [&](const TrialParamSet *params_ptr) -> double {
    return evaluate_outcome_coupling_unified(
        ctx, *coupling_program_ptr, upper, component_idx, params_ptr,
        trial_type_key, rel_tol, abs_tol, max_depth, include_na_donors,
        outcome_idx_context, forced_complete_bits, forced_complete_bits_valid,
        forced_survive_bits, forced_survive_bits_valid);
  };

  // Ensure we have a concrete parameter set to modify for shared triggers.
  TrialParamSet base_params_holder;
  const TrialParamSet *params_ptr = trial_params;
  if (!params_ptr) {
    base_params_holder = build_base_paramset(ctx);
    params_ptr = &base_params_holder;
  }

  // Enumerate shared-trigger states to preserve correlation across linked
  // accumulators. Reuse precomputed plans/scratch when available.
  SharedTriggerPlan local_trigger_plan;
  const SharedTriggerPlan *plan_ptr = trigger_plan;
  if (!plan_ptr) {
    local_trigger_plan = build_shared_trigger_plan(ctx, params_ptr);
    plan_ptr = &local_trigger_plan;
  }
  const std::size_t trigger_count = shared_trigger_count(*plan_ptr);
  if (trigger_count == 0) {
    return integrate_with_params(params_ptr);
  }

  if (trigger_count >= 63) {
    Rcpp::stop("Too many shared triggers for outcome probability evaluation");
  }
  TrialParamSet local_scratch;
  TrialParamSet &scratch = trigger_scratch ? *trigger_scratch : local_scratch;
  prepare_trigger_scratch(ctx, params_ptr, *plan_ptr, scratch);
  double total = evaluate_shared_trigger_states(
      ctx, *plan_ptr, scratch,
      "Too many shared triggers for outcome probability evaluation",
      [&]() -> double { return integrate_with_params(&scratch); });
  return clamp_probability(total);
}

double native_outcome_probability_bits_impl_idx(
    const uuber::NativeContext &ctx, int node_id, double upper, int component_idx,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const Rcpp::IntegerVector &competitor_ids, double rel_tol, double abs_tol,
    int max_depth, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan,
    TrialParamSet *trigger_scratch,
    const OutcomeCouplingProgram *precompiled_coupling_program) {
  std::vector<int> competitor_ids_raw = integer_vector_to_std(competitor_ids, false);
  return native_outcome_probability_bits_impl_idx(
      ctx, node_id, upper, component_idx,
      forced_complete_bits, forced_complete_bits_valid,
      forced_survive_bits, forced_survive_bits_valid,
      competitor_ids_raw, rel_tol, abs_tol, max_depth, trial_params,
      trial_type_key, include_na_donors, outcome_idx_context,
      trigger_plan, trigger_scratch, precompiled_coupling_program);
}

double native_outcome_probability_impl_idx(
    const uuber::NativeContext &ctx, int node_id, double upper, int component_idx,
    SEXP forced_complete, SEXP forced_survive,
    const Rcpp::IntegerVector &competitor_ids, double rel_tol, double abs_tol,
    int max_depth, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan,
    TrialParamSet *trigger_scratch) {
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
      forced_survive_bits_valid, competitor_ids, rel_tol, abs_tol, max_depth,
      trial_params, trial_type_key, include_na_donors, outcome_idx_context,
      trigger_plan, trigger_scratch, nullptr);
}

double native_outcome_probability_impl_idx(
    SEXP ctxSEXP, int node_id, double upper, int component_idx,
    SEXP forced_complete, SEXP forced_survive,
    const Rcpp::IntegerVector &competitor_ids, double rel_tol, double abs_tol,
    int max_depth, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan, TrialParamSet *trigger_scratch) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  return native_outcome_probability_impl_idx(
      *ctx, node_id, upper, component_idx, forced_complete, forced_survive,
      competitor_ids, rel_tol, abs_tol, max_depth, trial_params, trial_type_key,
      include_na_donors, outcome_idx_context, trigger_plan, trigger_scratch);
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
    const ForcedStateView *forced_state_view) {
  return eval_event_ref_idx(
      ctx, label_ref, node_flags, t, component_idx, need, trial_params,
      trial_type_key, include_na_donors, outcome_idx_context, time_constraints,
      forced_scope_filter, forced_complete_bits,
      forced_survive_bits, forced_label_id_to_bit_idx, trial_params_soa,
      forced_state_view);
}

bool evaluator_eval_node_with_forced_state_view_batch(
    const uuber::NativeContext &ctx, int node_id_or_idx,
    const std::vector<double> &times, int component_idx, EvalNeed need,
    const TrialParamSet *trial_params, const std::string &trial_key,
    const TimeConstraintMap *time_constraints,
    const ForcedStateView &forced_state,
    uuber::KernelNodeBatchValues &out_values) {
  return eval_node_with_forced_state_view_batch(
      ctx, node_id_or_idx, times, component_idx, need, trial_params, trial_key,
      time_constraints, forced_state, out_values);
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

double evaluator_guard_effective_survival_internal(
    const GuardEvalInput &input, double t,
    const IntegrationSettings &settings) {
  return guard_effective_survival_internal(input, t, settings);
}

NodeEvalResult evaluator_eval_node_recursive_dense(int node_idx,
                                                   NodeEvalState &state,
                                                   EvalNeed need) {
  return eval_node_recursive_dense(node_idx, state, need);
}

uuber::KernelEventEvalFn
evaluator_make_kernel_event_eval(NodeEvalState &state) {
  return make_kernel_event_eval(state);
}

uuber::KernelEventBatchEvalFn
evaluator_make_kernel_event_eval_batch(NodeEvalState &state) {
  return make_kernel_event_eval_batch(state);
}

uuber::KernelGuardEvalFn
evaluator_make_kernel_guard_eval(NodeEvalState &state) {
  return make_kernel_guard_eval(state);
}

uuber::KernelGuardBatchEvalFn
evaluator_make_kernel_guard_eval_batch(NodeEvalState &state) {
  return make_kernel_guard_eval_batch(state);
}

bool evaluator_eval_node_batch_with_state_dense(
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    EvalNeed need, uuber::KernelBatchRuntimeState &batch_runtime,
    uuber::KernelNodeBatchValues &out_values) {
  return eval_node_batch_with_state_dense(node_idx, times, state, need,
                                          batch_runtime, out_values);
}

double evaluator_evaluate_survival_with_forced(
    int node_id, const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid, int component_idx, double t,
    const uuber::NativeContext &ctx, const std::string &trial_key,
    const TrialParamSet *trial_params,
    const TimeConstraintMap *time_constraints,
    uuber::KernelRuntimeState *kernel_runtime) {
  return evaluate_survival_with_forced(
      node_id, forced_complete_bits, forced_complete_bits_valid,
      forced_survive_bits, forced_survive_bits_valid, component_idx, t, ctx,
      trial_key, trial_params, time_constraints,
      kernel_runtime);
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
    uuber::KernelBatchRuntimeState *kernel_batch_runtime) {
  return node_density_with_competitors_batch_internal(
      ctx, node_id, times, component_idx, forced_complete_bits,
      forced_complete_bits_valid, forced_survive_bits,
      forced_survive_bits_valid, competitor_ids, trial_params, trial_type_key,
      include_na_donors, outcome_idx_context, density_out, time_constraints,
      kernel_batch_runtime);
}

double native_trial_mixture_internal_idx(
    uuber::NativeContext &ctx, int node_id, double t,
    const std::vector<int> &component_indices, const std::vector<double> &weights,
    const std::vector<int> &competitor_ids, const TrialParamSet *params_ptr,
    const std::vector<ComponentCacheEntry> &cache_entries,
    int keep_outcome_idx_override = std::numeric_limits<int>::min(),
    int keep_label_id_override = NA_INTEGER, bool keep_is_guess = false,
    const Rcpp::List &guess_donors = Rcpp::List(),
    const SharedTriggerPlan *shared_trigger_plan = nullptr,
    TrialParamSet *shared_trigger_scratch = nullptr) {
  if (!std::isfinite(t) || t < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double total = 0.0;
  const bool has_keep =
      (keep_outcome_idx_override != std::numeric_limits<int>::min()) ||
      (keep_label_id_override != NA_INTEGER) || keep_is_guess;
  int keep_outcome_idx = (keep_outcome_idx_override !=
                          std::numeric_limits<int>::min())
                             ? keep_outcome_idx_override
                             : -1;
  int keep_outcome_idx_context = (keep_outcome_idx >= 0) ? keep_outcome_idx : -1;
  int keep_label_id = keep_label_id_override;
  if (keep_label_id == NA_INTEGER && keep_outcome_idx >= 0 &&
      keep_outcome_idx < static_cast<int>(ctx.outcome_label_ids.size())) {
    keep_label_id =
        ctx.outcome_label_ids[static_cast<std::size_t>(keep_outcome_idx)];
  }
  SharedTriggerPlan local_trigger_plan;
  const SharedTriggerPlan *trigger_plan_ptr = shared_trigger_plan;
  if (!trigger_plan_ptr) {
    local_trigger_plan = build_shared_trigger_plan(ctx, params_ptr);
    trigger_plan_ptr = &local_trigger_plan;
  }
  TrialParamSet local_trigger_scratch;
  TrialParamSet *trigger_scratch_ptr =
      shared_trigger_scratch ? shared_trigger_scratch : &local_trigger_scratch;
  for (std::size_t i = 0; i < component_indices.size(); ++i) {
    int comp_idx = component_indices[i];
    const ComponentCacheEntry *cache_entry_ptr = nullptr;
    ComponentCacheEntry fallback;
    if (i < cache_entries.size()) {
      cache_entry_ptr = &cache_entries[i];
    } else {
      fallback = default_component_cache_entry_idx(ctx, comp_idx);
      cache_entry_ptr = &fallback;
    }
    double contrib = 0.0;
    bool used_guess_shortcut = false;
    if (has_keep) {
      if (comp_idx >= 0 &&
          comp_idx < static_cast<int>(ctx.component_info.size())) {
        const uuber::ComponentGuessPolicy &guess =
            ctx.component_info[static_cast<std::size_t>(comp_idx)].guess;
        bool guess_applies = false;
        if (guess.valid) {
          if (guess.target_is_guess) {
            guess_applies = keep_is_guess;
          } else if (guess.target_outcome_idx >= 0) {
            guess_applies = (keep_outcome_idx == guess.target_outcome_idx);
          } else if (guess.target_label_id != NA_INTEGER) {
            guess_applies =
                (keep_label_id != NA_INTEGER &&
                 keep_label_id == guess.target_label_id);
          }
        }
        if (guess_applies) {
          contrib = accumulate_component_guess_density_idx(
              ctx, keep_label_id, keep_outcome_idx, keep_is_guess, t, comp_idx, params_ptr,
              cache_entry_ptr->trial_type_key);
          used_guess_shortcut = true;
        }
      }
    }
    if (!used_guess_shortcut) {
      std::vector<int> comp_competitors;
      const std::vector<int> &comp_use = filter_competitor_ids(
          ctx, competitor_ids, comp_idx, comp_competitors);
      contrib = node_density_entry_idx(
          ctx, node_id, t, comp_idx, nullptr, false, nullptr, false, comp_use,
          params_ptr,
          cache_entry_ptr->trial_type_key, false,
          keep_outcome_idx_context, trigger_plan_ptr, trigger_scratch_ptr, true);
      if (!std::isfinite(contrib) || contrib < 0.0)
        contrib = 0.0;
      if (has_keep) {
        double keep_w = component_keep_weight(ctx, comp_idx, keep_outcome_idx);
        contrib *= keep_w;
      }
      if (guess_donors.size() > 0) {
        double donor_contrib = accumulate_plan_guess_density(
            ctx, guess_donors, t, comp_idx,
            cache_entry_ptr->trial_type_key, params_ptr, false);
        contrib += donor_contrib;
      }
    }
    const double w = (i < weights.size()) ? weights[i] : 0.0;
    total += w * contrib;
  }
  if (!std::isfinite(total) || total <= 0.0) {
    return 0.0;
  }
  return total;
}

// [[Rcpp::export(name = "native_trial_mixture_cpp")]]
double native_trial_mixture_driver(
    SEXP ctxSEXP, int node_id, double t, Rcpp::CharacterVector component_ids,
    Rcpp::NumericVector weights, Rcpp::Nullable<Rcpp::String> forced_component,
    Rcpp::IntegerVector competitor_ids,
    Rcpp::Nullable<Rcpp::DataFrame> trial_rows,
    Rcpp::Nullable<Rcpp::List> guess_donors) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::vector<std::string> components;
  std::vector<double> mix_weights;
  if (!build_component_mix(component_ids, weights, forced_component, components,
                           mix_weights)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::vector<int> component_indices;
  component_indices.reserve(components.size());
  for (const auto &comp_id : components) {
    component_indices.push_back(component_index_of(*ctx, comp_id));
  }
  std::vector<int> comp_ids = integer_vector_to_std(competitor_ids, false);
  std::unique_ptr<TrialParamSet> param_holder;
  const TrialParamSet *params_ptr = nullptr;
  if (!trial_rows.isNull()) {
    param_holder = build_trial_params_from_df(*ctx, trial_rows);
    params_ptr = param_holder.get();
  }
  std::vector<ComponentCacheEntry> cache_entries =
      build_component_cache_entries_from_indices(*ctx, component_indices);
  return native_trial_mixture_internal_idx(
      *ctx, node_id, t, component_indices, mix_weights, comp_ids,
      params_ptr, cache_entries, std::numeric_limits<int>::min(),
      NA_INTEGER, false,
      guess_donors.isNotNull() ? Rcpp::List(guess_donors.get()) : Rcpp::List());
}

double mix_outcome_mass_idx(uuber::NativeContext &ctx,
                            const uuber::OutcomeContextInfo &info,
                            const std::vector<int> &component_indices,
                            const std::vector<double> &comp_weights,
                            const TrialParamSet *params_ptr, double rel_tol,
                            double abs_tol, int max_depth,
                            int outcome_idx_context,
                            const SharedTriggerPlan *trigger_plan = nullptr,
                            TrialParamSet *trigger_scratch = nullptr,
                            bool use_na_cache = true,
                            const std::vector<OutcomeCouplingProgram>
                                *precompiled_coupling_programs = nullptr) {
  if (info.node_id < 0)
    return 0.0;
  // Cache NA-map integrals when no RT is provided (maps_to_na cases).
  uuber::NAMapCacheKey cache_key;
  bool use_cache = use_na_cache && params_ptr != nullptr;
  if (use_cache) {
    cache_key = build_na_map_cache_key_idx(info, params_ptr, component_indices,
                                           comp_weights, outcome_idx_context);
    auto hit = ctx.na_map_cache.find(cache_key);
    if (hit != ctx.na_map_cache.end()) {
      return hit->second;
    }
  }
  uuber::BitsetState empty_forced_complete_bits;
  uuber::BitsetState empty_forced_survive_bits;
  bool empty_forced_complete_bits_valid = false;
  bool empty_forced_survive_bits_valid = false;
  (void)ensure_forced_bitset_capacity(ctx, empty_forced_complete_bits,
                                      empty_forced_complete_bits_valid);
  (void)ensure_forced_bitset_capacity(ctx, empty_forced_survive_bits,
                                      empty_forced_survive_bits_valid);
  double total = 0.0;
  for (std::size_t i = 0; i < component_indices.size(); ++i) {
    double w = (i < comp_weights.size()) ? comp_weights[i] : 0.0;
    if (w <= 0.0)
      continue;
    const int comp_idx = component_indices[i];
    if (!outcome_allows_component_idx(ctx, info, outcome_idx_context, comp_idx))
      continue;
    const OutcomeCouplingProgram *coupling_program_ptr = nullptr;
    if (precompiled_coupling_programs &&
        i < precompiled_coupling_programs->size()) {
      const OutcomeCouplingProgram &candidate =
          (*precompiled_coupling_programs)[i];
      if (candidate.valid) {
        coupling_program_ptr = &candidate;
      }
    }
    double p = native_outcome_probability_bits_impl_idx(
        ctx, info.node_id, std::numeric_limits<double>::infinity(), comp_idx,
        empty_forced_complete_bits_valid ? &empty_forced_complete_bits : nullptr,
        empty_forced_complete_bits_valid,
        empty_forced_survive_bits_valid ? &empty_forced_survive_bits : nullptr,
        empty_forced_survive_bits_valid, info.competitor_ids, rel_tol, abs_tol,
        max_depth, params_ptr, std::string(), false, outcome_idx_context,
        trigger_plan, trigger_scratch, coupling_program_ptr);
    if (std::isfinite(p) && p > 0.0) {
      double keep_w = component_keep_weight(ctx, comp_idx, outcome_idx_context);
      total += w * p * keep_w;
    }
  }
  double out = clamp_probability(total);
  if (use_cache) {
    if (static_cast<int>(ctx.na_map_cache.size()) >= ctx.na_cache_limit) {
      ctx.na_map_cache.clear();
      ctx.na_cache_order.clear();
    }
    ctx.na_map_cache.emplace(std::move(cache_key), out);
  }
  return out;
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
  int mass_outcome_idx{-1};
  const std::vector<OutcomeCouplingProgram> *mass_coupling_programs_ptr{nullptr};
  std::vector<OutcomeCouplingProgram> mass_coupling_programs;
  int runtime_slot{-1};
};

inline const std::string &empty_trial_type_key_ref() {
  static const std::string kEmptyTrialTypeKey;
  return kEmptyTrialTypeKey;
}

inline const std::vector<int> &empty_competitor_ids_ref() {
  static const std::vector<int> kEmptyCompetitorIds;
  return kEmptyCompetitorIds;
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
  const std::vector<std::vector<OutcomeCouplingProgram>>
      *nonresponse_coupling_programs{nullptr};
};

struct TrialEvalScratch {
  std::vector<TrialContributionSpec> contributions;
  std::unique_ptr<RankedTransitionCompiler> sequence_transition_compiler;
  std::vector<uuber::KernelRuntimeState> contribution_runtimes;
  std::vector<uuber::KernelRuntimeState *> kernel_runtime_ptrs;
};

bool build_trial_contributions_unified(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    std::vector<TrialContributionSpec> &contributions,
    std::unique_ptr<RankedTransitionCompiler> &sequence_transition_compiler);

struct SingleStepDirectBatchSpec {
  int node_id{-1};
  int component_idx{-1};
  int outcome_idx_context{-1};
  int param_group_idx{-1};
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
  std::string trial_type_key;
  std::vector<int> competitor_ids;

  bool operator==(const SingleStepDirectBatchKey &other) const {
    return node_id == other.node_id && component_idx == other.component_idx &&
           outcome_idx_context == other.outcome_idx_context &&
           param_group_idx == other.param_group_idx &&
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
  std::unique_ptr<RankedTransitionCompiler> sequence_transition_compiler;
  if (!build_trial_contributions_unified(
          ctx, eval_input, component_indices, component_weights, cache_entries,
          contributions, sequence_transition_compiler) ||
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
  eval_input.nonresponse_coupling_programs = nullptr;

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

inline void prepare_ranked_trial_eval_input(
    TrialEvalInput &eval_input, const std::vector<int> &trial_ranked_ids,
    const std::vector<double> &trial_ranked_rt,
    const std::vector<int> *nonresponse_outcome_indices) {
  if (trial_ranked_ids.size() != trial_ranked_rt.size()) {
    reset_trial_eval_input(eval_input, nonresponse_outcome_indices);
    return;
  }
  initialize_ranked_trial_eval_input_stage(eval_input, trial_ranked_ids.size());
  eval_input.sequence_label_storage = trial_ranked_ids;
  eval_input.sequence_time_storage = trial_ranked_rt;
  finalize_ranked_trial_eval_input_inline(eval_input,
                                          nonresponse_outcome_indices);
}

inline void build_nonresponse_coupling_program_cache(
    const uuber::NativeContext &ctx,
    const std::vector<int> &nonresponse_outcome_indices,
    const std::vector<int> &component_indices,
    std::vector<std::vector<OutcomeCouplingProgram>> &out) {
  out.clear();
  out.resize(nonresponse_outcome_indices.size());
  if (component_indices.empty()) {
    return;
  }
  for (std::size_t i = 0; i < nonresponse_outcome_indices.size(); ++i) {
    const int oi = nonresponse_outcome_indices[i];
    if (oi < 0 || oi >= static_cast<int>(ctx.outcome_info.size())) {
      continue;
    }
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(oi)];
    if (info.node_id < 0) {
      continue;
    }
    std::vector<OutcomeCouplingProgram> &row = out[i];
    row.reserve(component_indices.size());
    for (int comp_idx : component_indices) {
      OutcomeCouplingProgram coupling_program;
      if (outcome_allows_component_idx(ctx, info, oi, comp_idx)) {
        std::vector<int> filtered_competitors;
        const std::vector<int> &competitors = filter_competitor_ids(
            ctx, info.competitor_ids, comp_idx, filtered_competitors);
        coupling_program = resolve_outcome_coupling_program_with_generic(
            ctx, info.node_id, competitors, true);
      }
      row.push_back(std::move(coupling_program));
    }
  }
}

bool build_trial_contributions_unified(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    std::vector<TrialContributionSpec> &contributions,
    std::unique_ptr<RankedTransitionCompiler> &sequence_transition_compiler) {
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
      spec.mass_outcome_idx = oi;
      if (info.node_id >= 0 && !component_indices.empty()) {
        const std::vector<OutcomeCouplingProgram> *precomputed_programs = nullptr;
        if (use_precomputed && eval_input.nonresponse_coupling_programs &&
            i < eval_input.nonresponse_coupling_programs->size()) {
          const std::vector<OutcomeCouplingProgram> &candidate =
              (*eval_input.nonresponse_coupling_programs)[i];
          if (candidate.size() == component_indices.size()) {
            precomputed_programs = &candidate;
          }
        }
        if (precomputed_programs) {
          spec.mass_coupling_programs_ptr = precomputed_programs;
        } else {
          spec.mass_coupling_programs.reserve(component_indices.size());
          for (int comp_idx : component_indices) {
            OutcomeCouplingProgram coupling_program;
            if (outcome_allows_component_idx(ctx, info, oi, comp_idx)) {
              std::vector<int> filtered_competitors;
              const std::vector<int> &competitors = filter_competitor_ids(
                  ctx, info.competitor_ids, comp_idx, filtered_competitors);
              coupling_program = resolve_outcome_coupling_program_with_generic(
                  ctx, info.node_id, competitors, true);
            }
            spec.mass_coupling_programs.push_back(std::move(coupling_program));
          }
          spec.mass_coupling_programs_ptr = &spec.mass_coupling_programs;
        }
      }
      const bool internal_program_storage =
          (spec.mass_coupling_programs_ptr == &spec.mass_coupling_programs);
      contributions.push_back(std::move(spec));
      if (internal_program_storage) {
        TrialContributionSpec &stored = contributions.back();
        stored.mass_coupling_programs_ptr = &stored.mass_coupling_programs;
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

  if (uses_sequence_state && !contributions.empty() &&
      !sequence_transition_compiler) {
    sequence_transition_compiler =
        std::make_unique<RankedTransitionCompiler>(ctx);
  }
  return true;
}

inline double evaluate_sequence_density_single_step_lower_layer(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, double t, int component_idx,
    int outcome_idx_context, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    uuber::KernelRuntimeState *kernel_runtime) {
  if (node_idx < 0 || !std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr;
  if (!competitor_ids.empty()) {
    competitor_cache = &fetch_competitor_cluster_cache(ctx, competitor_ids);
  }
  NodeEvalState state(ctx, t, component_idx, trial_params, trial_type_key, false,
                      outcome_idx_context, nullptr, nullptr, nullptr, false,
                      nullptr, false, kernel_runtime);
  return node_density_with_competitors_from_state_dense(
      ctx, node_idx, competitor_ids, state, competitor_cache);
}

double evaluate_sequence_density_kernel_idx(
    uuber::NativeContext &ctx, const TrialContributionSpec &spec,
    const TrialParamSet *active_params,
    RankedTransitionCompiler *sequence_transition_compiler,
    SharedTriggerPlan *single_state_plan, TrialParamSet *prob_scratch_ptr,
    uuber::KernelRuntimeState *runtime_ptr) {
  const std::string &trial_type_key =
      spec.trial_type_key_ptr ? *spec.trial_type_key_ptr
                              : empty_trial_type_key_ref();
  const std::size_t sequence_length = spec.sequence_outcome_indices.size();
  if (spec.intent != TrialKernelIntent::SequenceDensity || sequence_length == 0u ||
      spec.sequence_node_indices.size() != sequence_length ||
      spec.sequence_time_data == nullptr) {
    return 0.0;
  }
  switch (spec.sequence_execution) {
  case TrialSequenceExecutionKind::LowerLayerDirect: {
    if (sequence_length != 1u) {
      return 0.0;
    }
    const int sequence_outcome_idx = spec.sequence_outcome_indices[0];
    const int sequence_node_idx = spec.sequence_node_indices[0];
    const double sequence_time = spec.sequence_time_data[0];
    if (sequence_outcome_idx < 0 || sequence_node_idx < 0 ||
        !std::isfinite(sequence_time) || sequence_time < 0.0) {
      return 0.0;
    }
    const std::vector<int> &competitors =
        (!spec.step_competitor_ids_ptrs.empty() &&
         spec.step_competitor_ids_ptrs[0] != nullptr)
            ? *spec.step_competitor_ids_ptrs[0]
            : empty_competitor_ids_ref();
    return evaluate_sequence_density_single_step_lower_layer(
        ctx, sequence_node_idx, competitors, sequence_time, spec.component_idx,
        sequence_outcome_idx, active_params, trial_type_key, runtime_ptr);
  }
  case TrialSequenceExecutionKind::SequenceState:
    if (!sequence_transition_compiler) {
      return 0.0;
    }
    return sequence_prefix_density_resolved(
        ctx, spec.sequence_outcome_indices, spec.sequence_node_indices,
        spec.sequence_time_data, spec.component_idx, active_params,
        trial_type_key, runtime_ptr, sequence_transition_compiler,
        &spec.step_competitor_ids_ptrs, &spec.step_persistent_sources);
  }
  return 0.0;
}

double evaluate_trial_contribution_kernel_idx(
    uuber::NativeContext &ctx, const TrialContributionSpec &spec,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const TrialParamSet *active_params, double rel_tol, double abs_tol,
    int max_depth, RankedTransitionCompiler *sequence_transition_compiler,
    SharedTriggerPlan *single_state_plan, TrialParamSet *prob_scratch_ptr,
    uuber::KernelRuntimeState *runtime_ptr) {
  if (spec.intent == TrialKernelIntent::OutcomeMass) {
    const int oi = spec.mass_outcome_idx;
    if (oi < 0 || oi >= static_cast<int>(ctx.outcome_info.size())) {
      return 0.0;
    }
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(oi)];
    const std::vector<OutcomeCouplingProgram> *coupling_programs_ptr =
        spec.mass_coupling_programs_ptr;
    if (!coupling_programs_ptr && !spec.mass_coupling_programs.empty()) {
      coupling_programs_ptr = &spec.mass_coupling_programs;
    }
    return mix_outcome_mass_idx(
        ctx, info, component_indices, component_weights, active_params, rel_tol,
        abs_tol, max_depth, oi, single_state_plan, prob_scratch_ptr, true,
        coupling_programs_ptr);
  }
  return evaluate_sequence_density_kernel_idx(
      ctx, spec, active_params, sequence_transition_compiler, single_state_plan,
      prob_scratch_ptr, runtime_ptr);
}

double evaluate_trial_probability_kernel_idx(
    uuber::NativeContext &ctx, const TrialEvalInput &eval_input,
    const std::vector<int> &component_indices,
    const std::vector<double> &component_weights,
    const std::vector<ComponentCacheEntry> &cache_entries,
    const TrialParamSet *params_ptr, double rel_tol, double abs_tol,
    int max_depth, const SharedTriggerPlan *trigger_plan,
    TrialParamSet *prob_trigger_scratch, TrialEvalScratch *eval_scratch) {
  if (!eval_input.valid) {
    return eval_input.probability_transform ==
                   TrialProbabilityTransform::Complement
               ? 1.0
               : 0.0;
  }

  TrialEvalScratch local_eval_scratch;
  TrialEvalScratch *eval_scratch_ptr =
      eval_scratch ? eval_scratch : &local_eval_scratch;
  std::vector<TrialContributionSpec> &contributions =
      eval_scratch_ptr->contributions;
  std::unique_ptr<RankedTransitionCompiler> &sequence_transition_compiler =
      eval_scratch_ptr->sequence_transition_compiler;

  TrialParamSet local_prob_scratch;
  TrialParamSet *prob_scratch_ptr =
      prob_trigger_scratch ? prob_trigger_scratch : &local_prob_scratch;
  auto finalize_probability = [&](double total) -> double {
    if (!std::isfinite(total)) {
      return 0.0;
    }
    if (eval_input.probability_transform ==
        TrialProbabilityTransform::Complement) {
      return clamp_probability(1.0 - total);
    }
    return (total > 0.0) ? total : 0.0;
  };

  if (!build_trial_contributions_unified(
          ctx, eval_input, component_indices, component_weights, cache_entries,
          contributions, sequence_transition_compiler)) {
    return eval_input.probability_transform ==
                   TrialProbabilityTransform::Complement
               ? 1.0
               : 0.0;
  }

  if (contributions.empty()) {
    return eval_input.probability_transform ==
                   TrialProbabilityTransform::Complement
               ? 1.0
               : 0.0;
  }

  int runtime_count = 0;
  for (TrialContributionSpec &spec : contributions) {
    spec.runtime_slot = -1;
    if (!ctx.kernel_program.valid ||
        spec.intent != TrialKernelIntent::SequenceDensity) {
      continue;
    }
    spec.runtime_slot = runtime_count;
    ++runtime_count;
  }
  std::vector<uuber::KernelRuntimeState> &contribution_runtimes =
      eval_scratch_ptr->contribution_runtimes;
  contribution_runtimes.clear();
  contribution_runtimes.resize(
      static_cast<std::size_t>(std::max(0, runtime_count)));

  if (!sequence_transition_compiler) {
    bool needs_sequence_compiler = false;
    for (const TrialContributionSpec &spec : contributions) {
      if (spec.intent == TrialKernelIntent::SequenceDensity &&
          spec.sequence_execution == TrialSequenceExecutionKind::SequenceState) {
        needs_sequence_compiler = true;
        break;
      }
    }
    if (needs_sequence_compiler) {
      sequence_transition_compiler =
          std::make_unique<RankedTransitionCompiler>(ctx);
    }
  }

  SharedTriggerPlan single_state_plan;
  auto evaluate_contribution = [&](const TrialContributionSpec &spec,
                                   const TrialParamSet *active_params) -> double {
    uuber::KernelRuntimeState *runtime_ptr = nullptr;
    if (spec.runtime_slot >= 0 &&
        spec.runtime_slot < static_cast<int>(contribution_runtimes.size())) {
      runtime_ptr =
          &contribution_runtimes[static_cast<std::size_t>(spec.runtime_slot)];
    }
    return evaluate_trial_contribution_kernel_idx(
        ctx, spec, component_indices, component_weights, active_params, rel_tol,
        abs_tol, max_depth, sequence_transition_compiler.get(),
        &single_state_plan, prob_scratch_ptr, runtime_ptr);
  };

  auto eval_single_state = [&](const TrialParamSet *active_params) -> double {
    double total = 0.0;
    for (const TrialContributionSpec &spec : contributions) {
      const double contrib = evaluate_contribution(spec, active_params);
      if (std::isfinite(contrib) && contrib > 0.0) {
        total += spec.scaled_weight * contrib;
      }
    }
    return total;
  };

  if (!params_ptr) {
    const double total = eval_single_state(nullptr);
    return finalize_probability(total);
  }

  SharedTriggerPlan local_plan;
  const SharedTriggerPlan *plan_ptr = trigger_plan;
  if (!plan_ptr) {
    local_plan = build_shared_trigger_plan(ctx, params_ptr);
    plan_ptr = &local_plan;
  }
  const std::size_t trigger_count = shared_trigger_count(*plan_ptr);
  if (trigger_count == 0) {
    const double total = eval_single_state(params_ptr);
    return finalize_probability(total);
  }
  if (trigger_count >= 63) {
    Rcpp::stop("Too many shared triggers for trial likelihood evaluation");
  }

  TrialParamSet local_outer_scratch;
  TrialParamSet *outer_scratch_ptr =
      prob_scratch_ptr ? prob_scratch_ptr : &local_outer_scratch;
  if (!outer_scratch_ptr) {
    outer_scratch_ptr = &local_outer_scratch;
  }
  prepare_trigger_scratch(ctx, params_ptr, *plan_ptr, *outer_scratch_ptr);

  std::vector<uuber::KernelRuntimeState *> &kernel_runtimes =
      eval_scratch_ptr->kernel_runtime_ptrs;
  kernel_runtimes.clear();
  kernel_runtimes.reserve(contribution_runtimes.size());
  for (std::size_t i = 0; i < contribution_runtimes.size(); ++i) {
    kernel_runtimes.push_back(&contribution_runtimes[i]);
  }
  const double total = evaluate_shared_trigger_states(
      ctx, *plan_ptr, *outer_scratch_ptr,
      "Too many shared triggers for trial likelihood evaluation",
      [&]() -> double { return eval_single_state(outer_scratch_ptr); }, nullptr,
      kernel_runtimes.empty() ? nullptr : &kernel_runtimes);
  return finalize_probability(total);
}

// [[Rcpp::export]]
double cpp_loglik(SEXP ctxSEXP, Rcpp::NumericMatrix params_mat,
                  Rcpp::DataFrame data_df, Rcpp::LogicalVector ok,
                  Rcpp::IntegerVector expand, double min_ll, double rel_tol,
                  double abs_tol, int max_depth) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  if (!ctx->ir.valid) {
    Rcpp::stop("loglik requires a compiled IR context");
  }
  if (!ctx->kernel_program.valid) {
    Rcpp::stop("loglik requires a compiled kernel program");
  }
  ctx->na_map_cache.clear();
  ctx->na_cache_order.clear();
  const int q_col = 0;
  const int w_col = 1;
  const int t0_col = 2;
  std::array<int, 8> p_cols{};
  p_cols.fill(-1);
  const int n_param_cols = std::max(0, params_mat.ncol() - 3);
  const int n_slot_cols = std::min(8, n_param_cols);
  for (int i = 0; i < n_slot_cols; ++i) {
    p_cols[static_cast<std::size_t>(i)] = 3 + i;
  }

  const ComponentMap &comp_map = ctx->components;
  const std::vector<std::string> &comp_ids = comp_map.ids;
  const std::vector<int> &outcome_label_ids = ctx->outcome_label_ids;
  if (outcome_label_ids.empty()) {
    Rcpp::stop("IR loglik requires native outcome labels in context");
  }

  auto decode_outcome_factor = [&](SEXP col_sexp) -> Rcpp::IntegerVector {
    const R_xlen_t n = Rf_xlength(col_sexp);
    Rcpp::IntegerVector out(n, NA_INTEGER);
    Rcpp::IntegerVector codes(col_sexp);
    for (R_xlen_t i = 0; i < n; ++i) {
      const int code = codes[i];
      if (code == NA_INTEGER) {
        continue;
      }
      out[i] = outcome_label_ids[static_cast<std::size_t>(code - 1)];
    }
    return out;
  };

  auto decode_component_factor = [&](SEXP col_sexp) -> Rcpp::IntegerVector {
    const R_xlen_t n = Rf_xlength(col_sexp);
    Rcpp::IntegerVector out(n, NA_INTEGER);
    Rcpp::IntegerVector codes(col_sexp);
    for (R_xlen_t i = 0; i < n; ++i) {
      const int code = codes[i];
      if (code == NA_INTEGER) {
        continue;
      }
      out[i] = code - 1;
    }
    return out;
  };

  Rcpp::NumericVector trial_col = data_df["trial"];
  const bool has_r = data_df.containsElementNamed("R");
  Rcpp::IntegerVector outcome_id_col =
      has_r ? decode_outcome_factor(data_df["R"]) : Rcpp::IntegerVector();
  Rcpp::NumericVector rt_col = data_df["rt"];
  std::vector<Rcpp::IntegerVector> rank_outcome_id_cols;
  std::vector<Rcpp::NumericVector> rank_rt_cols;
  rank_outcome_id_cols.push_back(outcome_id_col);
  rank_rt_cols.push_back(rt_col);
  for (int rank = 2;; ++rank) {
    std::string r_name = "R" + std::to_string(rank);
    std::string rt_name = "rt" + std::to_string(rank);
    bool has_r = data_df.containsElementNamed(r_name.c_str());
    bool has_rt = data_df.containsElementNamed(rt_name.c_str());
    if (!(has_r && has_rt)) {
      break;
    }
    rank_outcome_id_cols.push_back(decode_outcome_factor(data_df[r_name]));
    rank_rt_cols.push_back(Rcpp::NumericVector(data_df[rt_name]));
  }
  const int rank_width = static_cast<int>(rank_outcome_id_cols.size());
  const bool ranked_mode = rank_width > 1;
  if (outcome_id_col.size() == 0) {
    Rcpp::stop("IR loglik requires factor column 'R'");
  }
  const bool has_component_label = data_df.containsElementNamed("component");
  bool has_component = has_component_label;
  Rcpp::IntegerVector component_idx_col =
      has_component_label
          ? decode_component_factor(data_df["component"])
          : Rcpp::IntegerVector();
  bool has_onset = data_df.containsElementNamed("onset");
  Rcpp::NumericVector onset_col =
      has_onset ? Rcpp::NumericVector(data_df["onset"]) : Rcpp::NumericVector();

  const std::vector<TrialAccumulatorParams> base_acc_params =
      build_base_paramset(*ctx).acc_params;
  std::vector<int> all_acc_indices(ctx->accumulators.size());
  for (std::size_t i = 0; i < all_acc_indices.size(); ++i)
    all_acc_indices[i] = static_cast<int>(i);
  std::vector<std::vector<int>> comp_acc_indices(comp_ids.size());
  for (std::size_t c = 0; c < comp_ids.size(); ++c) {
    const std::string &cid = comp_ids[c];
    for (std::size_t a = 0; a < ctx->accumulators.size(); ++a) {
      const auto &comps = ctx->accumulators[a].components;
      if (std::find(comps.begin(), comps.end(), cid) != comps.end()) {
        comp_acc_indices[c].push_back(static_cast<int>(a));
      }
    }
    if (comp_acc_indices[c].empty()) {
      comp_acc_indices[c] = all_acc_indices;
    }
  }

  std::vector<TrialParamSet> param_sets;
  std::vector<TrialEvalInput> trial_eval_inputs;
  std::vector<int> outcome_label_id_by_trial;
  std::vector<double> rt_by_trial;
  std::vector<int> comp_idx_by_trial;
  std::vector<std::vector<double>> weight_override_by_trial;

  const std::size_t comp_count = comp_map.ids.size();

  const R_xlen_t n_rows = params_mat.nrow();
  if (trial_col.size() != n_rows) {
    Rcpp::stop(
        "Parameter matrix rows (%lld) must match nested likelihood rows (%lld)",
        static_cast<long long>(n_rows),
        static_cast<long long>(trial_col.size()));
  }
  int current_trial_label = std::numeric_limits<int>::min();
  int trial_idx = -1;
  std::size_t acc_cursor = 0;
  const std::vector<int> *current_acc_order = &all_acc_indices;
  for (R_xlen_t r = 0; r < n_rows; ++r) {
    int trial_label = static_cast<int>(trial_col[r]);
    if (r == 0 || trial_label != current_trial_label) {
      current_trial_label = trial_label;
      ++trial_idx;
      TrialParamSet ps;
      ps.acc_params = base_acc_params;
      param_sets.push_back(std::move(ps));
      trial_eval_inputs.push_back(TrialEvalInput{});
      outcome_label_id_by_trial.push_back(-1);
      rt_by_trial.push_back(std::numeric_limits<double>::quiet_NaN());
      comp_idx_by_trial.push_back(-1);
      weight_override_by_trial.push_back(std::vector<double>(
          comp_count, std::numeric_limits<double>::quiet_NaN()));
      if (ranked_mode) {
        initialize_ranked_trial_eval_input_stage(
            trial_eval_inputs.back(), static_cast<std::size_t>(rank_width));
      }
      acc_cursor = 0;
      current_acc_order = &all_acc_indices;
    }
    if (has_component && comp_idx_by_trial[trial_idx] < 0 &&
        r < component_idx_col.size()) {
      int comp_idx_val = component_idx_col[r];
      if (comp_idx_val != NA_INTEGER) {
        comp_idx_by_trial[trial_idx] = comp_idx_val;
        current_acc_order =
            &comp_acc_indices[static_cast<std::size_t>(comp_idx_val)];
      }
    }
    int acc_idx = (acc_cursor < current_acc_order->size())
                      ? (*current_acc_order)[acc_cursor]
                      : all_acc_indices[acc_cursor % all_acc_indices.size()];
    ++acc_cursor;
    TrialAccumulatorParams &tap =
        param_sets[static_cast<std::size_t>(trial_idx)]
            .acc_params[static_cast<std::size_t>(acc_idx)];
    double q_val = clamp_probability(params_mat(r, q_col));
    if (!tap.shared_trigger_id.empty()) {
      tap.shared_q = q_val;
      tap.q = 0.0;
    } else {
      tap.q = q_val;
      tap.shared_q = q_val;
    }
    tap.dist_cfg.t0 = params_mat(r, t0_col);
    if (has_onset && r < onset_col.size()) {
      double onset_val = static_cast<double>(onset_col[r]);
      if (std::isfinite(onset_val)) {
        tap.onset = onset_val;
      }
    }
    tap.has_override = true;
    param_sets[static_cast<std::size_t>(trial_idx)].has_any_override = true;
    param_sets[static_cast<std::size_t>(trial_idx)].soa_cache_valid = false;
    const int p1_col = p_cols[0];
    const int p2_col = p_cols[1];
    const int p3_col = p_cols[2];
    const int p4_col = p_cols[3];
    const int p5_col = p_cols[4];
    const int p6_col = p_cols[5];
    const int p7_col = p_cols[6];
    const int p8_col = p_cols[7];
    if (p1_col >= 0)
      tap.dist_cfg.p1 = params_mat(r, p1_col);
    if (p2_col >= 0)
      tap.dist_cfg.p2 = params_mat(r, p2_col);
    if (p3_col >= 0)
      tap.dist_cfg.p3 = params_mat(r, p3_col);
    if (p4_col >= 0)
      tap.dist_cfg.p4 = params_mat(r, p4_col);
    if (p5_col >= 0)
      tap.dist_cfg.p5 = params_mat(r, p5_col);
    if (p6_col >= 0)
      tap.dist_cfg.p6 = params_mat(r, p6_col);
    if (p7_col >= 0)
      tap.dist_cfg.p7 = params_mat(r, p7_col);
    if (p8_col >= 0)
      tap.dist_cfg.p8 = params_mat(r, p8_col);

    for (std::size_t c = 0; c < comp_map.leader_idx.size(); ++c) {
      if (comp_map.leader_idx[c] == acc_idx) {
        weight_override_by_trial[static_cast<std::size_t>(trial_idx)][c] =
            params_mat(r, w_col);
      }
    }

    if (outcome_label_id_by_trial[trial_idx] < 0 && r < outcome_id_col.size()) {
      int id_val = outcome_id_col[r];
      if (id_val != NA_INTEGER) {
        outcome_label_id_by_trial[trial_idx] = id_val;
      }
    }
    if (!std::isfinite(rt_by_trial[trial_idx]) && r < rt_col.size()) {
      rt_by_trial[trial_idx] = static_cast<double>(rt_col[r]);
    }
    if (ranked_mode) {
      TrialEvalInput &trial_eval_input =
          trial_eval_inputs[static_cast<std::size_t>(trial_idx)];
      std::vector<int> &trial_ranked_ids = trial_eval_input.sequence_label_storage;
      std::vector<double> &trial_ranked_rt =
          trial_eval_input.sequence_time_storage;
      for (int rank_idx = 0; rank_idx < rank_width; ++rank_idx) {
        if (trial_ranked_ids[static_cast<std::size_t>(rank_idx)] < 0 &&
            r < rank_outcome_id_cols[static_cast<std::size_t>(rank_idx)].size()) {
          int id_val = rank_outcome_id_cols[static_cast<std::size_t>(rank_idx)][r];
          if (id_val != NA_INTEGER) {
            trial_ranked_ids[static_cast<std::size_t>(rank_idx)] = id_val;
          }
        }
        if (!std::isfinite(
                trial_ranked_rt[static_cast<std::size_t>(rank_idx)]) &&
            r < rank_rt_cols[static_cast<std::size_t>(rank_idx)].size()) {
          double cand_rt =
              rank_rt_cols[static_cast<std::size_t>(rank_idx)][r];
          if (std::isfinite(cand_rt)) {
            trial_ranked_rt[static_cast<std::size_t>(rank_idx)] = cand_rt;
          }
        }
      }
    }
  }

  int n_trials = static_cast<int>(param_sets.size());
  std::vector<SharedTriggerPlan> trigger_plans;
  trigger_plans.reserve(static_cast<std::size_t>(n_trials));
  std::vector<int> trigger_plan_index_by_trial(static_cast<std::size_t>(n_trials),
                                               -1);
  if (ctx->shared_trigger_groups.empty()) {
    for (TrialParamSet &ps : param_sets) {
      ps.shared_trigger_source_fingerprint = 1ULL;
    }
    if (n_trials > 0) {
      trigger_plans.emplace_back();
      std::fill(trigger_plan_index_by_trial.begin(),
                trigger_plan_index_by_trial.end(), 0);
    }
  } else {
    for (TrialParamSet &ps : param_sets) {
      refresh_trial_param_fingerprint(ps);
    }
    std::vector<int> trigger_plan_representative_trial;
    trigger_plan_representative_trial.reserve(static_cast<std::size_t>(n_trials));
    std::unordered_map<std::uint64_t, std::vector<int>> trigger_plan_candidates;
    trigger_plan_candidates.reserve(static_cast<std::size_t>(n_trials));
    for (int t = 0; t < n_trials; ++t) {
      const TrialParamSet &params = param_sets[static_cast<std::size_t>(t)];
      const std::uint64_t fp = params.shared_trigger_source_fingerprint;
      std::vector<int> &candidates = trigger_plan_candidates[fp];
      int plan_idx = -1;
      for (int candidate_idx : candidates) {
        if (candidate_idx < 0 ||
            candidate_idx >=
                static_cast<int>(trigger_plan_representative_trial.size())) {
          continue;
        }
        const int rep_trial_idx =
            trigger_plan_representative_trial[static_cast<std::size_t>(
                candidate_idx)];
        if (rep_trial_idx < 0 || rep_trial_idx >= n_trials) {
          continue;
        }
        if (trial_paramsets_equivalent(
                params, param_sets[static_cast<std::size_t>(rep_trial_idx)])) {
          plan_idx = candidate_idx;
          break;
        }
      }
      if (plan_idx < 0) {
        plan_idx = static_cast<int>(trigger_plans.size());
        trigger_plans.push_back(
            build_shared_trigger_plan(*ctx,
                                      &param_sets[static_cast<std::size_t>(t)]));
        trigger_plan_representative_trial.push_back(t);
        candidates.push_back(plan_idx);
      }
      trigger_plan_index_by_trial[static_cast<std::size_t>(t)] = plan_idx;
    }
  }

  std::vector<int> param_value_group_by_trial(static_cast<std::size_t>(n_trials),
                                              -1);
  std::vector<int> param_value_group_representative_trial;
  param_value_group_representative_trial.reserve(
      static_cast<std::size_t>(n_trials));
  std::unordered_map<std::uint64_t, std::vector<int>> param_value_group_candidates;
  param_value_group_candidates.reserve(static_cast<std::size_t>(n_trials));
  for (int t = 0; t < n_trials; ++t) {
    const TrialParamSet &params = param_sets[static_cast<std::size_t>(t)];
    const std::uint64_t fp = compute_trial_param_value_fingerprint(params);
    std::vector<int> &candidates = param_value_group_candidates[fp];
    int group_idx = -1;
    for (int candidate_idx : candidates) {
      if (candidate_idx < 0 ||
          candidate_idx >=
              static_cast<int>(param_value_group_representative_trial.size())) {
        continue;
      }
      const int rep_trial_idx =
          param_value_group_representative_trial[static_cast<std::size_t>(
              candidate_idx)];
      if (rep_trial_idx < 0 || rep_trial_idx >= n_trials) {
        continue;
      }
      if (trial_paramsets_value_equivalent(
              params, param_sets[static_cast<std::size_t>(rep_trial_idx)])) {
        group_idx = candidate_idx;
        break;
      }
    }
    if (group_idx < 0) {
      group_idx = static_cast<int>(param_value_group_representative_trial.size());
      param_value_group_representative_trial.push_back(t);
      candidates.push_back(group_idx);
    }
    param_value_group_by_trial[static_cast<std::size_t>(t)] = group_idx;
  }

  std::vector<int> default_component_indices;
  default_component_indices.reserve(comp_map.ids.size());
  for (const auto &comp_id : comp_map.ids) {
    default_component_indices.push_back(component_index_of(*ctx, comp_id));
  }
  const std::size_t n_components = default_component_indices.size();
  const std::vector<ComponentCacheEntry> default_cache_entries =
      build_component_cache_entries_from_indices(*ctx, default_component_indices);
  std::vector<int> nonresponse_outcome_indices;
  nonresponse_outcome_indices.reserve(ctx->outcome_info.size());
  for (std::size_t oi = 0; oi < ctx->outcome_info.size(); ++oi) {
    bool is_deadline = false;
    const bool maps_to_na =
        ctx->outcome_info[static_cast<std::size_t>(oi)].maps_to_na;
    if (oi < ctx->ir.outcomes.size()) {
      int node_idx = ctx->ir.outcomes[oi].node_idx;
      if (node_idx >= 0 && node_idx < static_cast<int>(ctx->ir.nodes.size())) {
        const auto flags = ctx->ir.nodes[static_cast<std::size_t>(node_idx)].flags;
        is_deadline = (flags & uuber::IR_NODE_FLAG_SPECIAL_DEADLINE) != 0u;
      }
    }
    if (!is_deadline && !maps_to_na) {
      nonresponse_outcome_indices.push_back(static_cast<int>(oi));
    }
  }
  std::vector<std::vector<OutcomeCouplingProgram>>
      default_nonresponse_coupling_programs;
  build_nonresponse_coupling_program_cache(
      *ctx, nonresponse_outcome_indices, default_component_indices,
      default_nonresponse_coupling_programs);

  std::vector<std::vector<double>> weights_by_trial(
      static_cast<std::size_t>(n_trials),
      std::vector<double>(n_components, 0.0));
  for (int t = 0; t < n_trials; ++t) {
    std::vector<double> w = comp_map.base_weights;
    const auto &overrides =
        weight_override_by_trial[static_cast<std::size_t>(t)];
    for (std::size_t c = 0; c < w.size(); ++c) {
      double cand = overrides[c];
      if (std::isfinite(cand))
        w[c] = cand;
    }
    double total = 0.0;
    for (double v : w)
      total += v;
    if (!std::isfinite(total) || total <= 0.0) {
      double uni = 1.0 / static_cast<double>(w.size());
      std::fill(w.begin(), w.end(), uni);
    } else {
      for (double &v : w)
        v /= total;
    }
    weights_by_trial[static_cast<std::size_t>(t)] = std::move(w);
  }

  if (ok.size() == 0) {
    ok = Rcpp::LogicalVector(n_trials, true);
  } else if (ok.size() != n_trials) {
    Rcpp::stop("Length of ok must equal number of trials");
  }

  R_xlen_t n_expand = expand.size();
  if (n_expand == 0) {
    expand = Rcpp::seq(1, n_trials);
    n_expand = expand.size();
  }
  for (R_xlen_t i = 0; i < n_expand; ++i) {
    int comp_idx = expand[i];
    if (comp_idx == NA_INTEGER || comp_idx < 1 || comp_idx > n_trials) {
      Rcpp::stop("expand indices must reference trials");
    }
  }

  for (int t = 0; t < n_trials; ++t) {
    if (t >= static_cast<int>(trial_eval_inputs.size())) {
      trial_eval_inputs.resize(static_cast<std::size_t>(t + 1));
    }
    TrialEvalInput &eval_input = trial_eval_inputs[static_cast<std::size_t>(t)];
    if (ranked_mode) {
      finalize_ranked_trial_eval_input_inline(eval_input,
                                              &nonresponse_outcome_indices);
      continue;
    }
    const int outcome_label_id =
        (t < static_cast<int>(outcome_label_id_by_trial.size()))
            ? outcome_label_id_by_trial[static_cast<std::size_t>(t)]
            : -1;
    const double rt = (t < static_cast<int>(rt_by_trial.size()))
                          ? rt_by_trial[static_cast<std::size_t>(t)]
                          : std::numeric_limits<double>::quiet_NaN();
    prepare_observed_trial_eval_input(eval_input, outcome_label_id, rt,
                                      &nonresponse_outcome_indices);
  }

  struct ForcedComponentBundle {
    std::vector<int> component_indices;
    std::vector<double> component_weights;
    std::vector<ComponentCacheEntry> cache_entries;
    std::vector<std::vector<OutcomeCouplingProgram>>
        nonresponse_coupling_programs;
  };
  std::unordered_map<int, ForcedComponentBundle> forced_component_bundle_cache;
  forced_component_bundle_cache.reserve(static_cast<std::size_t>(n_trials));
  auto forced_bundle_for = [&](int forced_component_idx)
      -> const ForcedComponentBundle & {
    auto it = forced_component_bundle_cache.find(forced_component_idx);
    if (it != forced_component_bundle_cache.end()) {
      return it->second;
    }
    ForcedComponentBundle bundle;
    bundle.component_indices.push_back(forced_component_idx);
    bundle.component_weights.push_back(1.0);
    bundle.cache_entries.push_back(
        default_component_cache_entry_idx(*ctx, forced_component_idx));
    build_nonresponse_coupling_program_cache(
        *ctx, nonresponse_outcome_indices, bundle.component_indices,
        bundle.nonresponse_coupling_programs);
    auto inserted = forced_component_bundle_cache.emplace(forced_component_idx,
                                                          std::move(bundle));
    return inserted.first->second;
  };

  std::vector<double> trial_loglik(static_cast<std::size_t>(n_trials), min_ll);
  std::vector<std::uint8_t> trial_loglik_batched(
      static_cast<std::size_t>(n_trials), 0u);

  std::vector<SingleStepDirectBatchSpec> single_step_batch_specs(
      static_cast<std::size_t>(n_trials));
  std::unordered_map<SingleStepDirectBatchKey, std::vector<int>,
                     SingleStepDirectBatchKeyHash>
      single_step_batch_groups;
  single_step_batch_groups.reserve(static_cast<std::size_t>(n_trials));

  for (int t = 0; t < n_trials; ++t) {
    if (t >= ok.size() || !ok[t]) {
      continue;
    }

    const TrialEvalInput &eval_input =
        trial_eval_inputs[static_cast<std::size_t>(t)];
    const std::vector<int> *comp_indices_ptr = &default_component_indices;
    const std::vector<ComponentCacheEntry> *cache_entries_ptr =
        &default_cache_entries;
    const std::vector<double> *comp_weights_ptr =
        &weights_by_trial[static_cast<std::size_t>(t)];

    const int forced_component_idx =
        (t < static_cast<int>(comp_idx_by_trial.size()))
            ? comp_idx_by_trial[static_cast<std::size_t>(t)]
            : -1;
    if (forced_component_idx >= 0) {
      const ForcedComponentBundle &forced_bundle =
          forced_bundle_for(forced_component_idx);
      comp_indices_ptr = &forced_bundle.component_indices;
      comp_weights_ptr = &forced_bundle.component_weights;
      cache_entries_ptr = &forced_bundle.cache_entries;
    }

    const int trigger_plan_idx =
        (t >= 0 && t < static_cast<int>(trigger_plan_index_by_trial.size()))
            ? trigger_plan_index_by_trial[static_cast<std::size_t>(t)]
            : -1;
    const SharedTriggerPlan *trigger_plan_ptr =
        (trigger_plan_idx >= 0 &&
         trigger_plan_idx < static_cast<int>(trigger_plans.size()))
            ? &trigger_plans[static_cast<std::size_t>(trigger_plan_idx)]
            : nullptr;
    if (trigger_plan_ptr != nullptr &&
        shared_trigger_count(*trigger_plan_ptr) > 0u) {
      continue;
    }

    SingleStepDirectBatchSpec batch_spec;
    if (!build_single_step_direct_batch_spec(
            *ctx, eval_input, *comp_indices_ptr, *comp_weights_ptr,
            *cache_entries_ptr, batch_spec)) {
      continue;
    }

    const int param_group_idx =
        (t < static_cast<int>(param_value_group_by_trial.size()))
            ? param_value_group_by_trial[static_cast<std::size_t>(t)]
            : -1;
    if (param_group_idx < 0) {
      continue;
    }
    batch_spec.param_group_idx = param_group_idx;
    single_step_batch_specs[static_cast<std::size_t>(t)] = batch_spec;

    SingleStepDirectBatchKey key;
    key.node_id = batch_spec.node_id;
    key.component_idx = batch_spec.component_idx;
    key.outcome_idx_context = batch_spec.outcome_idx_context;
    key.param_group_idx = batch_spec.param_group_idx;
    key.trial_type_key = batch_spec.trial_type_key;
    key.competitor_ids = batch_spec.competitor_ids;
    single_step_batch_groups[std::move(key)].push_back(t);
  }

  for (auto &kv : single_step_batch_groups) {
    const std::vector<int> &trial_indices = kv.second;
    if (trial_indices.size() < 2u) {
      continue;
    }
    const int first_trial = trial_indices.front();
    if (first_trial < 0 || first_trial >= n_trials) {
      continue;
    }
    const SingleStepDirectBatchSpec &group_spec =
        single_step_batch_specs[static_cast<std::size_t>(first_trial)];
    if (group_spec.param_group_idx < 0 ||
        group_spec.param_group_idx >=
            static_cast<int>(param_value_group_representative_trial.size())) {
      continue;
    }
    const int rep_trial_idx =
        param_value_group_representative_trial[static_cast<std::size_t>(
            group_spec.param_group_idx)];
    if (rep_trial_idx < 0 || rep_trial_idx >= n_trials) {
      continue;
    }

    std::vector<double> times;
    times.reserve(trial_indices.size());
    for (int trial_idx : trial_indices) {
      const SingleStepDirectBatchSpec &trial_spec =
          single_step_batch_specs[static_cast<std::size_t>(trial_idx)];
      times.push_back(trial_spec.time);
    }

    std::vector<double> density_out;
    uuber::KernelBatchRuntimeState batch_runtime;
    const bool batched = node_density_entry_batch_idx(
        *ctx, group_spec.node_id, times, group_spec.component_idx, nullptr,
        false, nullptr, false, group_spec.competitor_ids,
        &param_sets[static_cast<std::size_t>(rep_trial_idx)],
        group_spec.trial_type_key, false, group_spec.outcome_idx_context,
        nullptr, false, density_out, nullptr, &batch_runtime);
    if (!batched || density_out.size() != trial_indices.size()) {
      continue;
    }

    for (std::size_t i = 0; i < trial_indices.size(); ++i) {
      const int trial_idx = trial_indices[i];
      const SingleStepDirectBatchSpec &trial_spec =
          single_step_batch_specs[static_cast<std::size_t>(trial_idx)];
      const double prob = trial_spec.scaled_weight * density_out[i];
      double ll_val = min_ll;
      if (std::isfinite(prob) && prob > 0.0) {
        const double lp = std::log(prob);
        if (std::isfinite(lp) && lp > min_ll) {
          ll_val = lp;
        }
      }
      trial_loglik[static_cast<std::size_t>(trial_idx)] = ll_val;
      trial_loglik_batched[static_cast<std::size_t>(trial_idx)] = 1u;
    }
  }

  TrialParamSet prob_trigger_scratch;
  TrialEvalScratch trial_eval_scratch;

  for (int t = 0; t < n_trials; ++t) {
    if (t >= ok.size() || !ok[t]) {
      trial_loglik[static_cast<std::size_t>(t)] = min_ll;
      continue;
    }
    if (trial_loglik_batched[static_cast<std::size_t>(t)] != 0u) {
      continue;
    }

    TrialParamSet *params_ptr = &param_sets[static_cast<std::size_t>(t)];
    TrialEvalInput eval_input = trial_eval_inputs[static_cast<std::size_t>(t)];
    const std::vector<int> *comp_indices_ptr = &default_component_indices;
    const std::vector<ComponentCacheEntry> *cache_entries_ptr =
        &default_cache_entries;
    const std::vector<double> *comp_weights_ptr =
        &weights_by_trial[static_cast<std::size_t>(t)];
    const std::vector<std::vector<OutcomeCouplingProgram>>
        *nonresponse_coupling_programs_ptr =
            &default_nonresponse_coupling_programs;

    int forced_component_idx =
        (t < static_cast<int>(comp_idx_by_trial.size()))
            ? comp_idx_by_trial[static_cast<std::size_t>(t)]
            : -1;
    if (forced_component_idx >= 0) {
      const ForcedComponentBundle &forced_bundle =
          forced_bundle_for(forced_component_idx);
      comp_indices_ptr = &forced_bundle.component_indices;
      comp_weights_ptr = &forced_bundle.component_weights;
      cache_entries_ptr = &forced_bundle.cache_entries;
      nonresponse_coupling_programs_ptr =
          &forced_bundle.nonresponse_coupling_programs;
    }

    eval_input.nonresponse_coupling_programs = nonresponse_coupling_programs_ptr;

    const int trigger_plan_idx =
        (t >= 0 && t < static_cast<int>(trigger_plan_index_by_trial.size()))
            ? trigger_plan_index_by_trial[static_cast<std::size_t>(t)]
            : -1;
    const SharedTriggerPlan *trigger_plan_ptr =
        (trigger_plan_idx >= 0 &&
         trigger_plan_idx < static_cast<int>(trigger_plans.size()))
            ? &trigger_plans[static_cast<std::size_t>(trigger_plan_idx)]
            : nullptr;

    const double prob = evaluate_trial_probability_kernel_idx(
        *ctx, eval_input, *comp_indices_ptr, *comp_weights_ptr,
        *cache_entries_ptr, params_ptr, rel_tol, abs_tol, max_depth,
        trigger_plan_ptr, &prob_trigger_scratch, &trial_eval_scratch);

    double ll_val = min_ll;
    if (!std::isfinite(prob) || prob <= 0.0) {
      ll_val = min_ll;
    } else {
      double lp = std::log(prob);
      ll_val = (std::isfinite(lp) && lp > min_ll) ? lp : min_ll;
    }
    trial_loglik[static_cast<std::size_t>(t)] = ll_val;
  }

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < n_expand; ++i) {
    int comp_idx = expand[i];
    total_loglik += trial_loglik[static_cast<std::size_t>(comp_idx - 1)];
  }

  return total_loglik;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_loglik_multiple(SEXP ctxSEXP, Rcpp::List params_list,
                                        Rcpp::DataFrame data_df,
                                        Rcpp::LogicalVector ok,
                                        Rcpp::IntegerVector expand,
                                        double min_ll, double rel_tol,
                                        double abs_tol, int max_depth) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  Rcpp::NumericVector out(params_list.size());
  for (R_xlen_t i = 0; i < params_list.size(); ++i) {
    if (params_list[i] == R_NilValue) {
      out[i] = NA_REAL;
      continue;
    }
    Rcpp::NumericMatrix pm(params_list[i]);
    out[i] = cpp_loglik(ctxSEXP, pm, data_df, ok, expand, min_ll, rel_tol,
                        abs_tol, max_depth);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List native_component_plan_exported(SEXP structureSEXP,
                                          SEXP trial_rowsSEXP,
                                          SEXP forced_componentSEXP) {
  Rcpp::List structure = Rcpp::as<Rcpp::List>(structureSEXP);
  Rcpp::DataFrame trial_rows_df;
  const Rcpp::DataFrame *trial_rows_ptr = nullptr;
  if (!Rf_isNull(trial_rowsSEXP)) {
    trial_rows_df = Rcpp::DataFrame(trial_rowsSEXP);
    trial_rows_ptr = &trial_rows_df;
  }
  double trial_id = extract_trial_id(trial_rows_ptr);
  std::string forced_component_value;
  const std::string *forced_component_ptr = nullptr;
  if (!Rf_isNull(forced_componentSEXP)) {
    forced_component_value = Rcpp::as<std::string>(forced_componentSEXP);
    forced_component_ptr = &forced_component_value;
  }
  return native_component_plan_impl(structure, trial_rows_ptr, trial_id,
                                    forced_component_ptr);
}

// [[Rcpp::export]]
double native_outcome_probability_params_cpp_idx(
    SEXP ctxSEXP, int node_id, double upper, int component_idx,
    SEXP forced_complete, SEXP forced_survive,
    Rcpp::IntegerVector competitor_ids, double rel_tol, double abs_tol,
    int max_depth, Rcpp::Nullable<Rcpp::DataFrame> trial_rows) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::unique_ptr<TrialParamSet> params_holder =
      build_trial_params_from_df(*ctx, trial_rows);
  return native_outcome_probability_impl_idx(
      ctxSEXP, node_id, upper, component_idx, forced_complete, forced_survive,
      competitor_ids, rel_tol, abs_tol, max_depth,
      params_holder ? params_holder.get() : nullptr, std::string(), false, -1,
      nullptr, nullptr);
}
