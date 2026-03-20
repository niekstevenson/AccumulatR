#include "ranked_transitions.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "competitor_cache.h"
#include "evaluator_internal.h"
#include "exact_outcome_density.h"

namespace {

struct SequenceStateKey {
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  int forced_complete_bit_count{0};
  int forced_survive_bit_count{0};
  std::uint64_t forced_complete_hash{0};
  std::uint64_t forced_survive_hash{0};
  std::uint64_t time_constraint_hash{0};
  std::size_t time_constraint_size{0};

  bool operator==(const SequenceStateKey &other) const noexcept {
    return forced_complete_bits_valid == other.forced_complete_bits_valid &&
           forced_survive_bits_valid == other.forced_survive_bits_valid &&
           forced_complete_bit_count == other.forced_complete_bit_count &&
           forced_survive_bit_count == other.forced_survive_bit_count &&
           forced_complete_hash == other.forced_complete_hash &&
           forced_survive_hash == other.forced_survive_hash &&
           time_constraint_hash == other.time_constraint_hash &&
           time_constraint_size == other.time_constraint_size;
  }
};

struct SequenceStateKeyHash {
  std::size_t operator()(const SequenceStateKey &key) const noexcept {
    std::size_t seed = static_cast<std::size_t>(key.forced_complete_hash);
    auto combine = [](std::size_t &seed_ref, std::size_t value) {
      seed_ref ^= value + 0x9e3779b97f4a7c15ULL + (seed_ref << 6) +
                  (seed_ref >> 2);
    };
    combine(seed, static_cast<std::size_t>(key.forced_survive_hash));
    combine(seed, static_cast<std::size_t>(key.time_constraint_hash));
    combine(seed, static_cast<std::size_t>(key.forced_complete_bits_valid));
    combine(seed, static_cast<std::size_t>(key.forced_survive_bits_valid));
    combine(seed, static_cast<std::size_t>(key.forced_complete_bit_count));
    combine(seed, static_cast<std::size_t>(key.forced_survive_bit_count));
    combine(seed, key.time_constraint_size);
    return seed;
  }
};

struct TrialSequenceState {
  int trial_slot{-1};
  double weight{0.0};
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  TimeConstraintMap time_constraints;
};

struct RankedFrontierStateRef {
  std::size_t state_index{0};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  SequenceStateKey state_key;
};

struct RankedStateCollapseIndex {
  std::size_t candidate_index{0};
  int trial_slot{-1};
  SequenceStateKey state_key;
};

struct RankedFrontierGroup {
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  std::size_t begin{0};
  std::size_t end{0};
};

inline std::uint64_t ranked_hash_bitset(const uuber::BitsetState &bits,
                                        bool bits_valid) {
  if (!bits_valid) {
    return 0ULL;
  }
  std::uint64_t hash = kFNV64Offset;
  const std::int32_t bit_count = static_cast<std::int32_t>(bits.bit_count());
  hash_append_bytes(hash, &bit_count, sizeof(bit_count));
  if (bits.words().empty()) {
    const std::uint64_t word = bits.small_word();
    hash_append_bytes(hash, &word, sizeof(word));
  } else {
    for (std::uint64_t word : bits.words()) {
      hash_append_bytes(hash, &word, sizeof(word));
    }
  }
  return hash;
}

inline std::uint64_t ranked_hash_time_constraints(
    const TimeConstraintMap &time_constraints) {
  std::uint64_t hash = kFNV64Offset;
  for (const auto &kv : time_constraints) {
    const std::int32_t id = static_cast<std::int32_t>(kv.first);
    hash_append_bytes(hash, &id, sizeof(id));
    const std::uint8_t has_exact = kv.second.has_exact ? 1u : 0u;
    const std::uint8_t has_lower = kv.second.has_lower ? 1u : 0u;
    const std::uint8_t has_upper = kv.second.has_upper ? 1u : 0u;
    hash_append_bytes(hash, &has_exact, sizeof(has_exact));
    hash_append_bytes(hash, &has_lower, sizeof(has_lower));
    hash_append_bytes(hash, &has_upper, sizeof(has_upper));
    if (kv.second.has_exact) {
      const std::uint64_t exact_bits = canonical_double_bits(kv.second.exact_time);
      hash_append_bytes(hash, &exact_bits, sizeof(exact_bits));
    }
    if (kv.second.has_lower) {
      const std::uint64_t lower_bits = canonical_double_bits(kv.second.lower);
      hash_append_bytes(hash, &lower_bits, sizeof(lower_bits));
    }
    if (kv.second.has_upper) {
      const std::uint64_t upper_bits = canonical_double_bits(kv.second.upper);
      hash_append_bytes(hash, &upper_bits, sizeof(upper_bits));
    }
  }
  return hash;
}

inline SequenceStateKey sequence_state_key(
    const uuber::BitsetState &forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState &forced_survive_bits,
    bool forced_survive_bits_valid,
    const TimeConstraintMap &time_constraints) {
  SequenceStateKey key;
  key.forced_complete_bits_valid = forced_complete_bits_valid;
  key.forced_survive_bits_valid = forced_survive_bits_valid;
  key.forced_complete_bit_count =
      forced_complete_bits_valid ? forced_complete_bits.bit_count() : 0;
  key.forced_survive_bit_count =
      forced_survive_bits_valid ? forced_survive_bits.bit_count() : 0;
  key.forced_complete_hash =
      ranked_hash_bitset(forced_complete_bits, forced_complete_bits_valid);
  key.forced_survive_hash =
      ranked_hash_bitset(forced_survive_bits, forced_survive_bits_valid);
  key.time_constraint_hash = ranked_hash_time_constraints(time_constraints);
  key.time_constraint_size = time_constraints.size();
  return key;
}

inline bool sequence_state_key_less(const SequenceStateKey &lhs,
                                    const SequenceStateKey &rhs) noexcept {
  if (lhs.forced_complete_bits_valid != rhs.forced_complete_bits_valid) {
    return lhs.forced_complete_bits_valid < rhs.forced_complete_bits_valid;
  }
  if (lhs.forced_survive_bits_valid != rhs.forced_survive_bits_valid) {
    return lhs.forced_survive_bits_valid < rhs.forced_survive_bits_valid;
  }
  if (lhs.forced_complete_bit_count != rhs.forced_complete_bit_count) {
    return lhs.forced_complete_bit_count < rhs.forced_complete_bit_count;
  }
  if (lhs.forced_survive_bit_count != rhs.forced_survive_bit_count) {
    return lhs.forced_survive_bit_count < rhs.forced_survive_bit_count;
  }
  if (lhs.forced_complete_hash != rhs.forced_complete_hash) {
    return lhs.forced_complete_hash < rhs.forced_complete_hash;
  }
  if (lhs.forced_survive_hash != rhs.forced_survive_hash) {
    return lhs.forced_survive_hash < rhs.forced_survive_hash;
  }
  if (lhs.time_constraint_hash != rhs.time_constraint_hash) {
    return lhs.time_constraint_hash < rhs.time_constraint_hash;
  }
  return lhs.time_constraint_size < rhs.time_constraint_size;
}

inline bool ranked_frontier_state_ref_less(const RankedFrontierStateRef &lhs,
                                           const RankedFrontierStateRef &rhs) noexcept {
  if (lhs.trial_params_soa != rhs.trial_params_soa) {
    return std::less<const uuber::TrialParamsSoA *>{}(lhs.trial_params_soa,
                                                      rhs.trial_params_soa);
  }
  return sequence_state_key_less(lhs.state_key, rhs.state_key);
}

inline bool ranked_frontier_state_ref_same_group(
    const RankedFrontierStateRef &lhs,
    const RankedFrontierStateRef &rhs) noexcept {
  return lhs.trial_params_soa == rhs.trial_params_soa &&
         lhs.state_key == rhs.state_key;
}

inline bool ranked_state_collapse_index_less(
    const RankedStateCollapseIndex &lhs,
    const RankedStateCollapseIndex &rhs) noexcept {
  if (lhs.trial_slot != rhs.trial_slot) {
    return lhs.trial_slot < rhs.trial_slot;
  }
  return sequence_state_key_less(lhs.state_key, rhs.state_key);
}

inline bool ranked_state_collapse_index_same_key(
    const RankedStateCollapseIndex &lhs,
    const RankedStateCollapseIndex &rhs) noexcept {
  return lhs.trial_slot == rhs.trial_slot && lhs.state_key == rhs.state_key;
}

inline void collapse_ranked_state_candidates(
    std::vector<TrialSequenceState> &states,
    std::vector<RankedStateCollapseIndex> &collapse_order,
    double branch_eps) {
  if (states.empty() || collapse_order.empty()) {
    states.clear();
    collapse_order.clear();
    return;
  }

  std::sort(collapse_order.begin(), collapse_order.end(),
            ranked_state_collapse_index_less);

  std::vector<TrialSequenceState> collapsed_states;
  std::vector<RankedStateCollapseIndex> collapsed_order;
  collapsed_states.reserve(collapse_order.size());
  collapsed_order.reserve(collapse_order.size());

  bool have_prev = false;
  RankedStateCollapseIndex prev_key;
  for (const RankedStateCollapseIndex &entry : collapse_order) {
    if (entry.candidate_index >= states.size()) {
      continue;
    }
    TrialSequenceState &candidate = states[entry.candidate_index];
    if (!std::isfinite(candidate.weight) || candidate.weight <= branch_eps) {
      continue;
    }
    if (have_prev && ranked_state_collapse_index_same_key(prev_key, entry)) {
      collapsed_states.back().weight += candidate.weight;
      continue;
    }
    const std::size_t collapsed_idx = collapsed_states.size();
    collapsed_states.push_back(std::move(candidate));
    RankedStateCollapseIndex collapsed_entry = entry;
    collapsed_entry.candidate_index = collapsed_idx;
    collapsed_order.push_back(std::move(collapsed_entry));
    prev_key = entry;
    have_prev = true;
  }

  states.swap(collapsed_states);
  collapse_order.swap(collapsed_order);
}

inline std::vector<int>
ranked_plan_child_ids(const uuber::NativeContext &ctx,
                      const uuber::IrNode &node) {
  std::vector<int> out;
  if (node.child_begin >= 0 && node.child_count > 0 &&
      node.child_begin + node.child_count <=
          static_cast<int>(ctx.ir.node_children.size())) {
    out.assign(ctx.ir.node_children.begin() + node.child_begin,
               ctx.ir.node_children.begin() + node.child_begin +
                   node.child_count);
  }
  return out;
}

inline std::vector<int>
ranked_source_bits_for_ids(const uuber::NativeContext &ctx,
                           const std::vector<int> &source_ids) {
  std::vector<int> out;
  out.reserve(source_ids.size());
  for (int id : source_ids) {
    int bit_idx = -1;
    auto bit_it = ctx.ir.label_id_to_bit_idx.find(id);
    if (bit_it != ctx.ir.label_id_to_bit_idx.end()) {
      bit_idx = bit_it->second;
    }
    out.push_back(bit_idx);
  }
  return out;
}

inline bool ranked_mask_covers_source_ids(const uuber::NativeContext &ctx,
                                          const RankedTransitionStep &step) {
  if (step.source_mask_begin < 0 || step.source_mask_count <= 0 ||
      step.source_ids.empty()) {
    return false;
  }
  if (step.source_mask_begin + step.source_mask_count >
      static_cast<int>(ctx.ir.node_source_masks.size())) {
    return false;
  }
  const std::uint64_t *mask_words =
      &ctx.ir.node_source_masks[static_cast<std::size_t>(step.source_mask_begin)];
  for (std::size_t i = 0; i < step.source_ids.size(); ++i) {
    const int source_id = step.source_ids[i];
    if (source_id == NA_INTEGER || source_id < 0) {
      continue;
    }
    int source_bit_idx = (i < step.source_bits.size()) ? step.source_bits[i] : -1;
    if (source_bit_idx < 0) {
      auto bit_it = ctx.ir.label_id_to_bit_idx.find(source_id);
      if (bit_it == ctx.ir.label_id_to_bit_idx.end()) {
        return false;
      }
      source_bit_idx = bit_it->second;
    }
    if (source_bit_idx < 0) {
      return false;
    }
    const int mask_word_idx = source_bit_idx / 64;
    const int mask_bit_offset = source_bit_idx % 64;
    if (mask_word_idx < 0 || mask_word_idx >= step.source_mask_count) {
      return false;
    }
    const std::uint64_t word =
        mask_words[static_cast<std::size_t>(mask_word_idx)];
    const std::uint64_t bit = (std::uint64_t{1} << mask_bit_offset);
    if ((word & bit) == 0u) {
      return false;
    }
  }
  return true;
}

inline int ranked_invalidate_slot_for_node(const uuber::NativeContext &ctx,
                                           int node_idx) {
  if (!ctx.kernel_program.valid) {
    return 0;
  }
  if (node_idx < 0 || node_idx >=
                          static_cast<int>(
                              ctx.kernel_program.outputs.node_idx_to_slot.size())) {
    return 0;
  }
  const int out_slot = ctx.kernel_program.outputs.node_idx_to_slot
      [static_cast<std::size_t>(node_idx)];
  if (out_slot < 0 ||
      out_slot >= static_cast<int>(ctx.kernel_program.ops.size())) {
    return 0;
  }
  return out_slot;
}

inline void ranked_attach_source_mask_for_node(const uuber::NativeContext &ctx,
                                               int node_idx,
                                               RankedTransitionStep &step) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return;
  }
  const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  if (node.source_mask_begin < 0 || node.source_mask_count <= 0) {
    return;
  }
  if (node.source_mask_begin + node.source_mask_count >
      static_cast<int>(ctx.ir.node_source_masks.size())) {
    return;
  }
  step.source_mask_begin = node.source_mask_begin;
  step.source_mask_count = node.source_mask_count;
  step.invalidate_slot = ranked_invalidate_slot_for_node(ctx, node_idx);
  step.source_mask_covers_ids = ranked_mask_covers_source_ids(ctx, step);
}

} // namespace

RankedTransitionCompiler::RankedTransitionCompiler(
    const uuber::NativeContext &ctx_)
    : ctx(ctx_), plans(ctx_.ir.nodes.size()) {}

const RankedNodeTransitionPlan &
RankedTransitionCompiler::plan_for_node(int node_idx) {
  static const RankedNodeTransitionPlan kInvalidPlan;
  if (node_idx < 0 || node_idx >= static_cast<int>(plans.size())) {
    return kInvalidPlan;
  }
  RankedNodeTransitionPlan &plan = plans[static_cast<std::size_t>(node_idx)];
  if (!plan.compiled) {
    compile_node(node_idx);
  }
  return plan;
}

void RankedTransitionCompiler::compile_node(int node_idx) {
  RankedNodeTransitionPlan &plan = plans[static_cast<std::size_t>(node_idx)];
  if (plan.compiled || plan.compiling) {
    return;
  }
  plan.compiling = true;
  plan.transitions.clear();
  plan.valid = false;

  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    plan.compiled = true;
    plan.compiling = false;
    return;
  }

  const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  switch (node.op) {
  case uuber::IrNodeOp::EventAcc:
  case uuber::IrNodeOp::EventPool: {
    RankedTransitionTemplate tr;
    RankedTransitionStep eval_step;
    eval_step.kind = RankedTransitionStepKind::EvalDensityNode;
    eval_step.node_idx = node_idx;
    tr.steps.push_back(std::move(eval_step));
    std::vector<int> source_ids = ensure_source_ids(ctx, node);
    if (!source_ids.empty()) {
      RankedTransitionStep add_step;
      add_step.kind = RankedTransitionStepKind::AddCompleteSources;
      add_step.bind_exact_current_time =
          (node.op == uuber::IrNodeOp::EventAcc);
      add_step.source_ids = std::move(source_ids);
      add_step.source_bits = ranked_source_bits_for_ids(ctx, add_step.source_ids);
      ranked_attach_source_mask_for_node(ctx, node_idx, add_step);
      tr.steps.push_back(std::move(add_step));
    }
    plan.transitions.push_back(std::move(tr));
    break;
  }
  case uuber::IrNodeOp::And: {
    std::vector<int> child_ids = ranked_plan_child_ids(ctx, node);
    if (!child_ids.empty()) {
      std::vector<std::vector<int>> child_sources;
      child_sources.reserve(child_ids.size());
      for (int child_idx : child_ids) {
        child_sources.push_back(
            ensure_source_ids(ctx, ir_node_required(ctx, child_idx)));
      }
      for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
        const RankedNodeTransitionPlan &child_plan =
            plan_for_node(child_ids[idx]);
        if (!child_plan.valid || child_plan.transitions.empty()) {
          continue;
        }
        for (const RankedTransitionTemplate &child_tr : child_plan.transitions) {
          RankedTransitionTemplate tr = child_tr;
          for (std::size_t other_idx = 0; other_idx < child_ids.size();
               ++other_idx) {
            if (other_idx == idx) {
              continue;
            }
            RankedTransitionStep cdf_step;
            cdf_step.kind = RankedTransitionStepKind::EvalCDFNode;
            cdf_step.node_idx = child_ids[other_idx];
            tr.steps.push_back(std::move(cdf_step));
            if (!child_sources[other_idx].empty()) {
              RankedTransitionStep add_step;
              add_step.kind = RankedTransitionStepKind::AddCompleteSources;
              add_step.source_ids = child_sources[other_idx];
              add_step.source_bits =
                  ranked_source_bits_for_ids(ctx, add_step.source_ids);
              ranked_attach_source_mask_for_node(ctx, child_ids[other_idx],
                                                 add_step);
              tr.steps.push_back(std::move(add_step));
            }
          }
          plan.transitions.push_back(std::move(tr));
        }
      }
    }
    break;
  }
  case uuber::IrNodeOp::Or: {
    std::vector<int> child_ids = ranked_plan_child_ids(ctx, node);
    if (!child_ids.empty()) {
      std::vector<std::vector<int>> child_sources;
      child_sources.reserve(child_ids.size());
      for (int child_idx : child_ids) {
        child_sources.push_back(
            ensure_source_ids(ctx, ir_node_required(ctx, child_idx)));
      }
      for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
        const RankedNodeTransitionPlan &child_plan =
            plan_for_node(child_ids[idx]);
        if (!child_plan.valid || child_plan.transitions.empty()) {
          continue;
        }
        for (const RankedTransitionTemplate &child_tr : child_plan.transitions) {
          RankedTransitionTemplate tr = child_tr;
          for (std::size_t other_idx = 0; other_idx < child_ids.size();
               ++other_idx) {
            if (other_idx == idx || child_sources[other_idx].empty()) {
              continue;
            }
            RankedTransitionStep witness_step;
            witness_step.kind =
                RankedTransitionStepKind::AddOrWitnessFromSources;
            witness_step.source_ids = child_sources[other_idx];
            witness_step.source_bits =
                ranked_source_bits_for_ids(ctx, witness_step.source_ids);
            tr.steps.push_back(std::move(witness_step));
          }
          plan.transitions.push_back(std::move(tr));
        }
      }
    }
    break;
  }
  case uuber::IrNodeOp::Guard: {
    if (node.reference_idx >= 0) {
      const RankedNodeTransitionPlan &ref_plan = plan_for_node(node.reference_idx);
      if (ref_plan.valid && !ref_plan.transitions.empty()) {
        const std::vector<int> blocker_sources =
            evaluator_gather_blocker_sources(ctx, node_idx);
        for (const RankedTransitionTemplate &ref_tr : ref_plan.transitions) {
          RankedTransitionTemplate tr = ref_tr;
          RankedTransitionStep eff_step;
          eff_step.kind = RankedTransitionStepKind::EvalGuardEffective;
          eff_step.node_idx = node_idx;
          tr.steps.push_back(std::move(eff_step));
          if (!blocker_sources.empty()) {
            RankedTransitionStep add_step;
            add_step.kind = RankedTransitionStepKind::AddSurviveSources;
            add_step.source_ids = blocker_sources;
            add_step.source_bits =
                ranked_source_bits_for_ids(ctx, add_step.source_ids);
            ranked_attach_source_mask_for_node(ctx, node.blocker_idx, add_step);
            tr.steps.push_back(std::move(add_step));
          }
          plan.transitions.push_back(std::move(tr));
        }
      }
    }
    break;
  }
  default:
    break;
  }

  plan.valid = !plan.transitions.empty();
  plan.compiled = true;
  plan.compiling = false;
}

std::vector<int>
collect_competitor_sources(const uuber::NativeContext &ctx,
                           const std::vector<int> &competitor_ids) {
  if (competitor_ids.empty()) {
    return {};
  }
  std::vector<int> out;
  for (int node_id : competitor_ids) {
    if (node_id == NA_INTEGER) {
      continue;
    }
    std::vector<int> source_ids = ensure_source_ids(ctx, node_id);
    out.insert(out.end(), source_ids.begin(), source_ids.end());
  }
  sort_unique(out);
  return out;
}

bool sequence_prefix_density_batch_resolved(
    const uuber::NativeContext &ctx, const std::vector<int> &outcome_indices,
    const std::vector<int> &node_indices,
    const std::vector<const double *> &times_by_trial, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::vector<double> &density_out,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_by_trial,
    RankedTransitionCompiler *transition_compiler,
    const std::vector<const std::vector<int> *> *step_competitor_ids_ptrs,
    const std::vector<std::vector<int>> *step_persistent_sources) {
  density_out.assign(times_by_trial.size(), 0.0);
  const std::size_t trial_count = times_by_trial.size();
  const std::size_t rank_count = outcome_indices.size();
  if (trial_count == 0u) {
    return true;
  }
  if (outcome_indices.empty() || node_indices.size() != rank_count ||
      step_competitor_ids_ptrs == nullptr ||
      step_persistent_sources == nullptr ||
      step_competitor_ids_ptrs->size() != rank_count ||
      step_persistent_sources->size() != rank_count) {
    return false;
  }
  for (const double *times_ptr : times_by_trial) {
    if (times_ptr == nullptr) {
      return false;
    }
  }

  RankedTransitionCompiler local_transition_compiler(ctx);
  RankedTransitionCompiler *compiler =
      transition_compiler ? transition_compiler : &local_transition_compiler;
  const uuber::TrialParamsSoA *default_trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);

  constexpr double kBranchEps = 1e-18;
  std::vector<TrialSequenceState> states;
  states.reserve(trial_count);
  for (std::size_t trial_slot = 0; trial_slot < trial_count; ++trial_slot) {
    TrialSequenceState init_state;
    init_state.trial_slot = static_cast<int>(trial_slot);
    init_state.weight = 1.0;
    (void)ensure_forced_bitset_capacity(ctx, init_state.forced_complete_bits,
                                        init_state.forced_complete_bits_valid);
    (void)ensure_forced_bitset_capacity(ctx, init_state.forced_survive_bits,
                                        init_state.forced_survive_bits_valid);
    states.push_back(std::move(init_state));
  }

  for (std::size_t rank_idx = 0; rank_idx < rank_count; ++rank_idx) {
    const int outcome_idx = outcome_indices[rank_idx];
    if (outcome_idx < 0 ||
        outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
      return false;
    }
    const int outcome_node_idx = node_indices[rank_idx];
    if (outcome_node_idx < 0 ||
        outcome_node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      return false;
    }
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
    const int info_node_idx = resolve_dense_node_idx_required(ctx, info.node_id);
    static const std::vector<int> kEmptyCompetitorIds;
    const std::vector<int> &competitors =
        ((*step_competitor_ids_ptrs)[rank_idx] != nullptr)
            ? *(*step_competitor_ids_ptrs)[rank_idx]
            : kEmptyCompetitorIds;
    const std::vector<int> &persistent_sources =
        (*step_persistent_sources)[rank_idx];

    std::vector<TrialSequenceState> next_states_collapsed;
    next_states_collapsed.reserve(states.size() * 2u);
    std::vector<RankedStateCollapseIndex> next_state_collapse_order;
    next_state_collapse_order.reserve(states.size() * 2u);

    const bool use_frontier_groups = (trial_params_soa_by_trial != nullptr);
    std::vector<RankedFrontierStateRef> active_state_refs;
    active_state_refs.reserve(states.size());
    std::vector<RankedFrontierGroup> frontier_groups;
    if (use_frontier_groups) {
      frontier_groups.reserve(states.size());
    }

    bool valid_rank = true;
    const uuber::TrialParamsSoA *frontier_uniform_trial_params_soa =
        default_trial_params_soa;
    bool frontier_uniform_trial_params_known = false;
    for (std::size_t state_idx = 0; state_idx < states.size(); ++state_idx) {
      const TrialSequenceState &state = states[state_idx];
      if (!std::isfinite(state.weight) || state.weight <= kBranchEps) {
        continue;
      }
      const int trial_slot = state.trial_slot;
      if (trial_slot < 0 ||
          static_cast<std::size_t>(trial_slot) >= times_by_trial.size()) {
        valid_rank = false;
        break;
      }
      const double t =
          times_by_trial[static_cast<std::size_t>(trial_slot)][rank_idx];
      if (!std::isfinite(t) || t < 0.0) {
        valid_rank = false;
        break;
      }
      if (rank_idx > 0u) {
        const double lower_t =
            times_by_trial[static_cast<std::size_t>(trial_slot)][rank_idx - 1u];
        if (!std::isfinite(lower_t) || lower_t < 0.0) {
          valid_rank = false;
          break;
        }
      }

      RankedFrontierStateRef state_ref;
      state_ref.state_index = state_idx;
      state_ref.trial_params_soa = default_trial_params_soa;
      if (trial_params_soa_by_trial != nullptr &&
          static_cast<std::size_t>(trial_slot) < trial_params_soa_by_trial->size()) {
        state_ref.trial_params_soa =
            (*trial_params_soa_by_trial)[static_cast<std::size_t>(trial_slot)];
      }
      state_ref.state_key = sequence_state_key(
          state.forced_complete_bits, state.forced_complete_bits_valid,
          state.forced_survive_bits, state.forced_survive_bits_valid,
          state.time_constraints);
      active_state_refs.push_back(std::move(state_ref));

      const uuber::TrialParamsSoA *state_trial_params_soa =
          active_state_refs.back().trial_params_soa;
      if (!frontier_uniform_trial_params_known) {
        frontier_uniform_trial_params_soa = state_trial_params_soa;
        frontier_uniform_trial_params_known = true;
      } else if (frontier_uniform_trial_params_soa != state_trial_params_soa) {
        frontier_uniform_trial_params_soa = nullptr;
      }
    }
    if (!valid_rank) {
      return false;
    }
    if (active_state_refs.empty()) {
      return true;
    }
    if (use_frontier_groups && active_state_refs.size() > 1u) {
      std::sort(active_state_refs.begin(), active_state_refs.end(),
                ranked_frontier_state_ref_less);
      std::size_t group_begin = 0u;
      while (group_begin < active_state_refs.size()) {
        std::size_t group_end = group_begin + 1u;
        while (group_end < active_state_refs.size() &&
               ranked_frontier_state_ref_same_group(
                   active_state_refs[group_begin],
                   active_state_refs[group_end])) {
          ++group_end;
        }
        RankedFrontierGroup group;
        group.trial_params_soa = active_state_refs[group_begin].trial_params_soa;
        group.begin = group_begin;
        group.end = group_end;
        frontier_groups.push_back(group);
        group_begin = group_end;
      }
    }

    std::vector<ExactScenarioPoint> rank_seed_points;
    rank_seed_points.reserve(active_state_refs.size());
    std::vector<ExactScenarioPoint> denom_points;
    denom_points.reserve(active_state_refs.size());
    std::vector<double> denom;
    denom.reserve(active_state_refs.size());
    std::vector<ExactScenarioPoint> scenario_points;
    scenario_points.reserve(active_state_refs.size() * 2u);
    std::vector<double> scenario_survival;
    scenario_survival.reserve(active_state_refs.size() * 2u);
    std::vector<ExactScenarioPoint> deterministic_scenario_points;
    deterministic_scenario_points.reserve(active_state_refs.size());
    std::vector<std::uint8_t> deterministic_scenario_active;
    deterministic_scenario_active.reserve(active_state_refs.size());
    std::vector<ExactScenarioPoint> deterministic_competitor_points;
    deterministic_competitor_points.reserve(active_state_refs.size());
    std::vector<double> deterministic_competitor_survival;
    deterministic_competitor_survival.reserve(active_state_refs.size());
    const uuber::CompetitorClusterCacheEntry *competitor_cache_ptr = nullptr;
    if (!competitors.empty()) {
      competitor_cache_ptr = &fetch_competitor_cluster_cache(ctx, competitors);
    }

    const bool split_frontier_groups =
        use_frontier_groups && frontier_groups.size() > 1u &&
        frontier_groups.size() * 2u < active_state_refs.size();
    const std::size_t frontier_begin = 0u;
    const std::size_t frontier_end = active_state_refs.size();

    auto process_state_group =
        [&](std::size_t group_begin, std::size_t group_end,
            const uuber::TrialParamsSoA *uniform_trial_params_soa) -> bool {
      rank_seed_points.clear();
      denom_points.clear();
      denom.clear();

      std::size_t local_idx = 0u;
      for (std::size_t ref_idx = group_begin; ref_idx < group_end;
           ++ref_idx, ++local_idx) {
        const RankedFrontierStateRef &state_ref = active_state_refs[ref_idx];
        const TrialSequenceState &state = states[state_ref.state_index];
        const std::size_t trial_slot = static_cast<std::size_t>(state.trial_slot);
        ExactScenarioPoint rank_point;
        rank_point.t = times_by_trial[trial_slot][rank_idx];
        rank_point.density_index = local_idx;
        rank_point.weight = 1.0;
        rank_point.trial_params_soa = uniform_trial_params_soa;
        if (rank_point.trial_params_soa == nullptr) {
          rank_point.trial_params_soa = state_ref.trial_params_soa;
        }
        rank_point.forced_complete_bits = state.forced_complete_bits;
        rank_point.forced_survive_bits = state.forced_survive_bits;
        rank_point.forced_complete_bits_valid = state.forced_complete_bits_valid;
        rank_point.forced_survive_bits_valid = state.forced_survive_bits_valid;
        rank_point.time_constraints = state.time_constraints;
        rank_seed_points.push_back(rank_point);

        if (rank_idx > 0u) {
          ExactScenarioPoint denom_point = rank_point;
          denom_point.t = times_by_trial[trial_slot][rank_idx - 1u];
          denom_points.push_back(std::move(denom_point));
        }
      }

      const std::size_t group_size = group_end - group_begin;
      denom.assign(group_size, 1.0);
      if (rank_idx > 0u) {
        uuber::KernelNodeBatchValues denom_values;
        if (!exact_eval_node_batch_from_points(
                ctx, info_node_idx, denom_points, component_idx, trial_params,
                trial_type_key, EvalNeed::kSurvival, denom_values,
                uniform_trial_params_soa) ||
            denom_values.survival.size() != denom.size()) {
          return false;
        }
        for (std::size_t i = 0; i < denom.size(); ++i) {
          denom[i] = denom_values.survival[i];
        }
      }

      scenario_points.clear();
      std::vector<TrialSequenceState> group_next_states;
      std::vector<RankedStateCollapseIndex> group_next_state_order;

      deterministic_scenario_points.clear();
      deterministic_scenario_active.clear();
      if (exact_collect_deterministic_scenarios_batch_from_points(
              *compiler, ctx, outcome_node_idx, rank_seed_points, component_idx,
              trial_params, trial_type_key, deterministic_scenario_points,
              deterministic_scenario_active, uniform_trial_params_soa)) {
        group_next_states.reserve(group_size);
        group_next_state_order.reserve(group_size);

        deterministic_competitor_points.clear();
        deterministic_competitor_survival.clear();
        if (competitor_cache_ptr != nullptr) {
          deterministic_competitor_points.reserve(group_size);
          for (std::size_t i = 0; i < deterministic_scenario_points.size(); ++i) {
            if (i < deterministic_scenario_active.size() &&
                deterministic_scenario_active[i] != 0u) {
              deterministic_competitor_points.push_back(
                  deterministic_scenario_points[i]);
            }
          }
          deterministic_competitor_survival.assign(
              deterministic_competitor_points.size(), 1.0);
          if (!deterministic_competitor_points.empty()) {
            exact_competitor_survival_batch(
                ctx, *competitor_cache_ptr, component_idx, trial_params,
                trial_type_key, deterministic_competitor_points,
                deterministic_competitor_survival, uniform_trial_params_soa);
          }
        }

        std::size_t competitor_survival_idx = 0u;
        for (std::size_t local_idx = 0; local_idx < group_size; ++local_idx) {
          if (local_idx >= deterministic_scenario_active.size() ||
              deterministic_scenario_active[local_idx] == 0u) {
            continue;
          }
          const ExactScenarioPoint &point = deterministic_scenario_points[local_idx];
          const RankedFrontierStateRef &source_state_ref =
              active_state_refs[group_begin + local_idx];
          const TrialSequenceState &source_state =
              states[source_state_ref.state_index];
          const double denom_val =
              (local_idx < denom.size()) ? denom[local_idx] : 0.0;
          if (!std::isfinite(denom_val) || denom_val <= kBranchEps) {
            continue;
          }
          double weight = source_state.weight * (point.weight / denom_val);
          if (!std::isfinite(weight) || weight <= kBranchEps) {
            continue;
          }
          if (competitor_cache_ptr != nullptr) {
            if (competitor_survival_idx >= deterministic_competitor_survival.size()) {
              return false;
            }
            const double surv =
                clamp_probability(deterministic_competitor_survival[
                    competitor_survival_idx++]);
            if (!std::isfinite(surv) || surv <= kBranchEps) {
              continue;
            }
            weight *= surv;
            if (!std::isfinite(weight) || weight <= kBranchEps) {
              continue;
            }
          }

          TrialSequenceState next_state;
          next_state.trial_slot = source_state.trial_slot;
          next_state.weight = weight;
          next_state.forced_complete_bits = point.forced_complete_bits;
          next_state.forced_survive_bits = point.forced_survive_bits;
          next_state.forced_complete_bits_valid = point.forced_complete_bits_valid;
          next_state.forced_survive_bits_valid = point.forced_survive_bits_valid;
          next_state.time_constraints = point.time_constraints;

          if (competitor_cache_ptr != nullptr && !persistent_sources.empty()) {
            bool keep = true;
            for (int src_id : persistent_sources) {
              set_forced_id_bit_strict(ctx, src_id, next_state.forced_survive_bits,
                                       next_state.forced_survive_bits_valid);
              if (!time_constraints_mark_survive(src_id, point.t,
                                                 next_state.time_constraints)) {
                keep = false;
                break;
              }
            }
            if (!keep) {
              continue;
            }
          }

          const std::size_t candidate_index = next_states_collapsed.size();
          next_states_collapsed.push_back(std::move(next_state));
          RankedStateCollapseIndex collapse_index;
          collapse_index.candidate_index = candidate_index;
          collapse_index.trial_slot = next_states_collapsed.back().trial_slot;
          collapse_index.state_key = sequence_state_key(
              next_states_collapsed.back().forced_complete_bits,
              next_states_collapsed.back().forced_complete_bits_valid,
              next_states_collapsed.back().forced_survive_bits,
              next_states_collapsed.back().forced_survive_bits_valid,
              next_states_collapsed.back().time_constraints);
          next_state_collapse_order.push_back(std::move(collapse_index));
        }
        return true;
      }

      if (!exact_collect_scenarios_batch_from_points(
              *compiler, ctx, outcome_node_idx, rank_seed_points, component_idx,
              trial_params, trial_type_key, scenario_points,
              uniform_trial_params_soa)) {
        return true;
      }

      scenario_survival.assign(scenario_points.size(), 1.0);
      if (competitor_cache_ptr != nullptr) {
        exact_competitor_survival_batch(
            ctx, *competitor_cache_ptr, component_idx, trial_params,
            trial_type_key, scenario_points, scenario_survival,
            uniform_trial_params_soa);
      }

      group_next_states.reserve(scenario_points.size());
      group_next_state_order.reserve(scenario_points.size());

      for (std::size_t point_idx = 0; point_idx < scenario_points.size();
           ++point_idx) {
        const ExactScenarioPoint &point = scenario_points[point_idx];
        if (point.density_index >= group_size) {
          return false;
        }
        const RankedFrontierStateRef &source_state_ref =
            active_state_refs[group_begin + point.density_index];
        const TrialSequenceState &source_state =
            states[source_state_ref.state_index];
        const double denom_val =
            (point.density_index < denom.size()) ? denom[point.density_index]
                                                : 0.0;
        if (!std::isfinite(denom_val) || denom_val <= kBranchEps) {
          continue;
        }
        double weight = source_state.weight * (point.weight / denom_val);
        if (!std::isfinite(weight) || weight <= kBranchEps) {
          continue;
        }

        if (competitor_cache_ptr != nullptr) {
          const double surv = clamp_probability(scenario_survival[point_idx]);
          if (!std::isfinite(surv) || surv <= kBranchEps) {
            continue;
          }
          weight *= surv;
          if (!std::isfinite(weight) || weight <= kBranchEps) {
            continue;
          }
        }

        TrialSequenceState next_state;
        next_state.trial_slot = source_state.trial_slot;
        next_state.weight = weight;
        next_state.forced_complete_bits = point.forced_complete_bits;
        next_state.forced_survive_bits = point.forced_survive_bits;
        next_state.forced_complete_bits_valid = point.forced_complete_bits_valid;
        next_state.forced_survive_bits_valid = point.forced_survive_bits_valid;
        next_state.time_constraints = point.time_constraints;

        if (competitor_cache_ptr != nullptr && !persistent_sources.empty()) {
          bool keep = true;
          for (int src_id : persistent_sources) {
            set_forced_id_bit_strict(ctx, src_id, next_state.forced_survive_bits,
                                     next_state.forced_survive_bits_valid);
            if (!time_constraints_mark_survive(src_id, point.t,
                                               next_state.time_constraints)) {
              keep = false;
              break;
            }
          }
          if (!keep) {
            continue;
          }
        }

        const std::size_t candidate_index = group_next_states.size();
        group_next_states.push_back(std::move(next_state));
        RankedStateCollapseIndex collapse_index;
        collapse_index.candidate_index = candidate_index;
        collapse_index.trial_slot = group_next_states.back().trial_slot;
        collapse_index.state_key = sequence_state_key(
            group_next_states.back().forced_complete_bits,
            group_next_states.back().forced_complete_bits_valid,
            group_next_states.back().forced_survive_bits,
            group_next_states.back().forced_survive_bits_valid,
            group_next_states.back().time_constraints);
        group_next_state_order.push_back(std::move(collapse_index));
      }

      collapse_ranked_state_candidates(group_next_states, group_next_state_order,
                                       kBranchEps);
      const std::size_t global_base = next_states_collapsed.size();
      next_states_collapsed.reserve(global_base + group_next_states.size());
      next_state_collapse_order.reserve(global_base + group_next_state_order.size());
      for (std::size_t i = 0; i < group_next_states.size(); ++i) {
        next_states_collapsed.push_back(std::move(group_next_states[i]));
      }
      for (RankedStateCollapseIndex &entry : group_next_state_order) {
        entry.candidate_index += global_base;
        next_state_collapse_order.push_back(std::move(entry));
      }
      return true;
    };

    if (split_frontier_groups) {
      for (const RankedFrontierGroup &group : frontier_groups) {
        if (!process_state_group(group.begin, group.end, group.trial_params_soa)) {
          return false;
        }
      }
    } else if (!process_state_group(frontier_begin, frontier_end,
                                    frontier_uniform_trial_params_soa)) {
      return false;
    }

    collapse_ranked_state_candidates(next_states_collapsed, next_state_collapse_order,
                                     kBranchEps);
    states.swap(next_states_collapsed);
    if (states.empty()) {
      return true;
    }
  }

  for (const TrialSequenceState &state : states) {
    if (!std::isfinite(state.weight) || state.weight <= 0.0 ||
        state.trial_slot < 0 ||
        static_cast<std::size_t>(state.trial_slot) >= density_out.size()) {
      continue;
    }
    density_out[static_cast<std::size_t>(state.trial_slot)] += state.weight;
  }
  for (double &value : density_out) {
    if (!std::isfinite(value) || value <= 0.0) {
      value = 0.0;
    }
  }

  return true;
}
