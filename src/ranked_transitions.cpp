#include "ranked_transitions.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "competitor_cache.h"
#include "evaluator_internal.h"
#include "exact_outcome_density.h"

namespace {

struct SequenceState {
  double weight{0.0};
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  TimeConstraintMap time_constraints;
};

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

struct TrialSequenceStateKey {
  int trial_slot{-1};
  SequenceStateKey state_key;

  bool operator==(const TrialSequenceStateKey &other) const noexcept {
    return trial_slot == other.trial_slot && state_key == other.state_key;
  }
};

struct TrialSequenceStateKeyHash {
  std::size_t operator()(const TrialSequenceStateKey &key) const noexcept {
    std::size_t seed =
        SequenceStateKeyHash{}(key.state_key) ^ static_cast<std::size_t>(key.trial_slot);
    seed ^= 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

inline bool ranked_batch_debug_enabled() {
  static const bool enabled = []() {
    const char *raw = std::getenv("ACCUMULATR_DEBUG_RANKED_BATCH");
    return raw != nullptr && raw[0] != '\0' && raw[0] != '0';
  }();
  return enabled;
}

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

template <typename EmitFn>
bool for_each_sequence_node_transition(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, double t, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const TimeConstraintMap *time_constraints,
    uuber::KernelRuntimeState *kernel_runtime,
    const uuber::BitsetState *base_forced_complete_bits_in,
    bool base_forced_complete_bits_in_valid,
    const uuber::BitsetState *base_forced_survive_bits_in,
    bool base_forced_survive_bits_in_valid, EmitFn &&emit,
    bool evaluate_step_weights = true) {
  const RankedNodeTransitionPlan &plan = compiler.plan_for_node(node_idx);
  if (!plan.valid || plan.transitions.empty()) {
    return false;
  }

  const uuber::TrialParamsSoA *trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);
  IntegrationSettings settings;
  bool emitted_any = false;
  uuber::BitsetState base_forced_complete_bits;
  uuber::BitsetState base_forced_survive_bits;
  bool base_forced_complete_bits_valid = false;
  bool base_forced_survive_bits_valid = false;
  if (base_forced_complete_bits_in && base_forced_complete_bits_in_valid) {
    base_forced_complete_bits = *base_forced_complete_bits_in;
    base_forced_complete_bits_valid = true;
  } else {
    (void)ensure_forced_bitset_capacity(ctx, base_forced_complete_bits,
                                        base_forced_complete_bits_valid);
  }
  if (base_forced_survive_bits_in && base_forced_survive_bits_in_valid) {
    base_forced_survive_bits = *base_forced_survive_bits_in;
    base_forced_survive_bits_valid = true;
  } else {
    (void)ensure_forced_bitset_capacity(ctx, base_forced_survive_bits,
                                        base_forced_survive_bits_valid);
  }

  for (const RankedTransitionTemplate &transition : plan.transitions) {
    double weight = 1.0;
    TimeConstraintMap transition_time_constraints =
        time_constraints ? *time_constraints : TimeConstraintMap{};
    NodeEvalState transition_state(
        ctx, t, component_idx, trial_params, trial_type_key, false, -1,
        &transition_time_constraints, nullptr,
        base_forced_complete_bits_valid ? &base_forced_complete_bits : nullptr,
        base_forced_complete_bits_valid,
        base_forced_survive_bits_valid ? &base_forced_survive_bits : nullptr,
        base_forced_survive_bits_valid, kernel_runtime);
    invalidate_kernel_runtime_root(transition_state);

    bool valid = true;
    for (const RankedTransitionStep &step : transition.steps) {
      if (step.kind == RankedTransitionStepKind::EvalDensityNode) {
        if (!evaluate_step_weights) {
          continue;
        }
        NodeEvalResult eval =
            evaluator_eval_node_recursive_dense(step.node_idx, transition_state,
                                                EvalNeed::kDensity);
        const double d = eval.density;
        if (!std::isfinite(d) || d <= 0.0) {
          valid = false;
          break;
        }
        weight *= d;
      } else if (step.kind == RankedTransitionStepKind::EvalCDFNode) {
        if (!evaluate_step_weights) {
          continue;
        }
        NodeEvalResult eval = evaluator_eval_node_recursive_dense(
            step.node_idx, transition_state, EvalNeed::kCDF);
        const double Fj = clamp_probability(eval.cdf);
        if (!std::isfinite(Fj) || Fj <= 0.0) {
          valid = false;
          break;
        }
        weight *= Fj;
      } else if (step.kind == RankedTransitionStepKind::EvalGuardEffective) {
        if (!evaluate_step_weights) {
          continue;
        }
        GuardEvalInput guard_input = make_guard_input_forced_state(
            ctx, step.node_idx, component_idx, &trial_type_key, trial_params,
            trial_params_soa, &transition_state.time_constraints,
            transition_state.forced_state);
        const double eff = evaluator_guard_effective_survival_internal(
            guard_input, t, settings);
        if (!std::isfinite(eff) || eff <= 0.0) {
          valid = false;
          break;
        }
        weight *= eff;
      } else if (step.kind == RankedTransitionStepKind::AddCompleteSources) {
        bool mutated_state = false;
        for (std::size_t i = 0; i < step.source_ids.size(); ++i) {
          const int id = step.source_ids[i];
          if (!time_constraints_mark_complete(
                  id, t, step.bind_exact_current_time,
                  transition_state.time_constraints)) {
            valid = false;
            break;
          }
          mutated_state = true;
        }
        if (!valid) {
          break;
        }
        if (mutated_state) {
          invalidate_kernel_runtime_root(transition_state);
        }
      } else if (step.kind == RankedTransitionStepKind::AddSurviveSources) {
        bool mutated_state = false;
        for (std::size_t i = 0; i < step.source_ids.size(); ++i) {
          const int id = step.source_ids[i];
          if (!time_constraints_mark_survive(id, t,
                                             transition_state.time_constraints)) {
            valid = false;
            break;
          }
          mutated_state = true;
        }
        if (!valid) {
          break;
        }
        if (mutated_state) {
          invalidate_kernel_runtime_root(transition_state);
        }
      } else if (step.kind == RankedTransitionStepKind::AddOrWitnessFromSources) {
        int witness = NA_INTEGER;
        int witness_bit = -1;
        bool all_forced = true;
        for (std::size_t i = 0; i < step.source_ids.size(); ++i) {
          const int id = step.source_ids[i];
          if (id == NA_INTEGER || id < 0) {
            continue;
          }
          const int bit_idx =
              (i < step.source_bits.size()) ? step.source_bits[i] : -1;
          int bit_to_check = bit_idx;
          if (bit_to_check < 0) {
            auto bit_it = ctx.ir.label_id_to_bit_idx.find(id);
            if (bit_it == ctx.ir.label_id_to_bit_idx.end()) {
              Rcpp::stop("IR ranked source id %d missing bit index", id);
            }
            bit_to_check = bit_it->second;
          }
          const bool is_forced =
              (transition_state.forced_complete_bits_valid &&
               transition_state.forced_complete_bits.test(bit_to_check)) ||
              time_constraints_contains_complete_at(
                  &transition_state.time_constraints, id, t);
          if (!is_forced) {
            witness = id;
            witness_bit = bit_to_check;
            all_forced = false;
            break;
          }
        }
        if (all_forced || witness == NA_INTEGER) {
          valid = false;
          break;
        }
        (void)witness_bit;
        if (!time_constraints_mark_survive(witness, t,
                                           transition_state.time_constraints)) {
          valid = false;
          break;
        }
        invalidate_kernel_runtime_root(transition_state);
      }
      if (!std::isfinite(weight) || weight <= 0.0) {
        valid = false;
        break;
      }
    }
    if (!valid || !std::isfinite(weight) || weight <= 0.0) {
      continue;
    }
    emit(weight, std::move(transition_state.forced_complete_bits),
         transition_state.forced_complete_bits_valid,
         std::move(transition_state.forced_survive_bits),
         transition_state.forced_survive_bits_valid,
         std::move(transition_state.time_constraints));
    emitted_any = true;
  }
  return emitted_any;
}

} // namespace

bool collect_exact_node_scenarios(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, double t, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const TimeConstraintMap *time_constraints,
    uuber::KernelRuntimeState *kernel_runtime,
    const uuber::BitsetState *base_forced_complete_bits_in,
    bool base_forced_complete_bits_in_valid,
    const uuber::BitsetState *base_forced_survive_bits_in,
    bool base_forced_survive_bits_in_valid,
    std::vector<ExactNodeScenario> &out,
    bool evaluate_step_weights) {
  out.clear();
  const bool emitted = for_each_sequence_node_transition(
      compiler, ctx, node_idx, t, component_idx, trial_params, trial_type_key,
      time_constraints, kernel_runtime,
      base_forced_complete_bits_in, base_forced_complete_bits_in_valid,
      base_forced_survive_bits_in, base_forced_survive_bits_in_valid,
      [&](double weight, uuber::BitsetState &&forced_complete_bits,
          bool forced_complete_bits_valid,
          uuber::BitsetState &&forced_survive_bits,
          bool forced_survive_bits_valid,
          TimeConstraintMap &&scenario_time_constraints) {
        ExactNodeScenario scenario;
        scenario.weight = weight;
        scenario.forced_complete_bits = std::move(forced_complete_bits);
        scenario.forced_survive_bits = std::move(forced_survive_bits);
        scenario.forced_complete_bits_valid = forced_complete_bits_valid;
        scenario.forced_survive_bits_valid = forced_survive_bits_valid;
        scenario.time_constraints = std::move(scenario_time_constraints);
        out.push_back(std::move(scenario));
      },
      evaluate_step_weights);
  return emitted;
}

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

double sequence_prefix_density_resolved(
    const uuber::NativeContext &ctx, const std::vector<int> &outcome_indices,
    const std::vector<int> &node_indices, const double *times,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, uuber::KernelRuntimeState *kernel_runtime,
    RankedTransitionCompiler *transition_compiler,
    const std::vector<const std::vector<int> *> *step_competitor_ids_ptrs,
    const std::vector<std::vector<int>> *step_persistent_sources) {
  const std::size_t rank_count = outcome_indices.size();
  if (outcome_indices.empty() || times == nullptr ||
      node_indices.size() != rank_count ||
      step_competitor_ids_ptrs == nullptr ||
      step_persistent_sources == nullptr ||
      step_competitor_ids_ptrs->size() != rank_count ||
      step_persistent_sources->size() != rank_count) {
    return 0.0;
  }
  RankedTransitionCompiler local_transition_compiler(ctx);
  RankedTransitionCompiler *compiler =
      transition_compiler ? transition_compiler : &local_transition_compiler;

  constexpr double kBranchEps = 1e-18;
  std::vector<SequenceState> states;
  SequenceState init_state;
  init_state.weight = 1.0;
  (void)ensure_forced_bitset_capacity(ctx, init_state.forced_complete_bits,
                                      init_state.forced_complete_bits_valid);
  (void)ensure_forced_bitset_capacity(ctx, init_state.forced_survive_bits,
                                      init_state.forced_survive_bits_valid);
  states.push_back(std::move(init_state));

  for (std::size_t rank_idx = 0; rank_idx < rank_count; ++rank_idx) {
    if (kernel_runtime) {
      uuber::invalidate_kernel_runtime_from_slot(*kernel_runtime, 0);
    }
    double t = times[rank_idx];
    if (!std::isfinite(t) || t < 0.0) {
      return 0.0;
    }

    int outcome_idx = outcome_indices[rank_idx];
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
    static const std::vector<int> kEmptyCompetitorIds;
    const std::vector<int> &competitors =
        ((*step_competitor_ids_ptrs)[rank_idx] != nullptr)
            ? *(*step_competitor_ids_ptrs)[rank_idx]
            : kEmptyCompetitorIds;
    const std::vector<int> &persistent_sources =
        (*step_persistent_sources)[rank_idx];

    std::unordered_map<SequenceStateKey, std::size_t, SequenceStateKeyHash>
        next_state_index;
    std::vector<SequenceState> next_states_collapsed;
    next_states_collapsed.reserve(states.size() * 2);

    auto accumulate_next_state = [&](SequenceState &&candidate) {
      if (!std::isfinite(candidate.weight) || candidate.weight <= kBranchEps) {
        return;
      }
      SequenceStateKey key = sequence_state_key(
          candidate.forced_complete_bits, candidate.forced_complete_bits_valid,
          candidate.forced_survive_bits, candidate.forced_survive_bits_valid,
          candidate.time_constraints);
      auto it = next_state_index.find(key);
      if (it == next_state_index.end()) {
        std::size_t idx = next_states_collapsed.size();
        next_state_index.emplace(std::move(key), idx);
        next_states_collapsed.push_back(std::move(candidate));
      } else {
        next_states_collapsed[it->second].weight += candidate.weight;
      }
    };

    for (const SequenceState &state : states) {
      if (!std::isfinite(state.weight) || state.weight <= kBranchEps) {
        continue;
      }
      TimeConstraintMap merged_time_constraints = state.time_constraints;
      uuber::BitsetState forced_complete_bits = state.forced_complete_bits;
      uuber::BitsetState forced_survive_bits = state.forced_survive_bits;
      bool forced_complete_bits_valid = state.forced_complete_bits_valid;
      bool forced_survive_bits_valid = state.forced_survive_bits_valid;
      double denom = 1.0;
      if (rank_idx > 0) {
        double lower_t = times[rank_idx - 1];
        denom = evaluator_evaluate_survival_with_forced(
            info.node_id,
            forced_complete_bits_valid ? &forced_complete_bits : nullptr,
            forced_complete_bits_valid,
            forced_survive_bits_valid ? &forced_survive_bits : nullptr,
            forced_survive_bits_valid, component_idx, lower_t, ctx,
            trial_type_key, trial_params, &merged_time_constraints,
            kernel_runtime);
        if (!std::isfinite(denom) || denom <= kBranchEps) {
          continue;
        }
      }
      int outcome_node_idx = node_indices[rank_idx];
      const bool has_transition = for_each_sequence_node_transition(
          *compiler, ctx, outcome_node_idx, t, component_idx, trial_params,
          trial_type_key, &merged_time_constraints, kernel_runtime,
          forced_complete_bits_valid ? &forced_complete_bits : nullptr,
          forced_complete_bits_valid,
          forced_survive_bits_valid ? &forced_survive_bits : nullptr,
          forced_survive_bits_valid,
          [&](double scenario_weight,
              uuber::BitsetState &&scenario_complete_bits,
              bool scenario_complete_bits_valid,
              uuber::BitsetState &&scenario_survive_bits,
              bool scenario_survive_bits_valid,
              TimeConstraintMap &&scenario_time_constraints) {
            if (!std::isfinite(scenario_weight) ||
                scenario_weight <= kBranchEps) {
              return;
            }
            double weight = state.weight * (scenario_weight / denom);
            if (!std::isfinite(weight) || weight <= kBranchEps) {
              return;
            }
            uuber::BitsetState next_complete_bits =
                std::move(scenario_complete_bits);
            uuber::BitsetState next_survive_bits =
                std::move(scenario_survive_bits);
            bool next_complete_bits_valid = scenario_complete_bits_valid;
            bool next_survive_bits_valid = scenario_survive_bits_valid;
            if (!competitors.empty()) {
              double surv = competitor_survival_internal(
                  ctx, competitors, t, component_idx,
                  next_complete_bits_valid ? &next_complete_bits : nullptr,
                  next_complete_bits_valid,
                  next_survive_bits_valid ? &next_survive_bits : nullptr,
                  next_survive_bits_valid, trial_type_key, trial_params,
                  &scenario_time_constraints,
                  kernel_runtime);
              if (!std::isfinite(surv) || surv <= kBranchEps) {
                return;
              }
              weight *= surv;
              if (!std::isfinite(weight) || weight <= kBranchEps) {
                return;
              }
              if (!persistent_sources.empty()) {
                for (int src_id : persistent_sources) {
                  set_forced_id_bit_strict(ctx, src_id, next_survive_bits,
                                           next_survive_bits_valid);
                  if (!time_constraints_mark_survive(
                          src_id, t, scenario_time_constraints)) {
                    return;
                  }
                }
              }
            }
            SequenceState next_state;
            next_state.weight = weight;
            next_state.forced_complete_bits = std::move(next_complete_bits);
            next_state.forced_survive_bits = std::move(next_survive_bits);
            next_state.forced_complete_bits_valid = next_complete_bits_valid;
            next_state.forced_survive_bits_valid = next_survive_bits_valid;
            next_state.time_constraints = std::move(scenario_time_constraints);
            accumulate_next_state(std::move(next_state));
          });
      if (!has_transition) {
        continue;
      }
    }

    states.clear();
    states.reserve(next_states_collapsed.size());
    for (SequenceState &next_state : next_states_collapsed) {
      if (std::isfinite(next_state.weight) && next_state.weight > kBranchEps) {
        states.push_back(std::move(next_state));
      }
    }
    if (states.empty()) {
      return 0.0;
    }
  }

  double total = 0.0;
  for (const SequenceState &state : states) {
    if (std::isfinite(state.weight) && state.weight > 0.0) {
      total += state.weight;
    }
  }
  if (!std::isfinite(total) || total <= 0.0) {
    return 0.0;
  }
  return total;
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

    std::unordered_map<TrialSequenceStateKey, std::size_t,
                       TrialSequenceStateKeyHash>
        next_state_index;
    std::vector<TrialSequenceState> next_states_collapsed;
    next_states_collapsed.reserve(states.size() * 2u);

    auto accumulate_next_state = [&](TrialSequenceState &&candidate) {
      if (!std::isfinite(candidate.weight) || candidate.weight <= kBranchEps) {
        return;
      }
      TrialSequenceStateKey key;
      key.trial_slot = candidate.trial_slot;
      key.state_key = sequence_state_key(
          candidate.forced_complete_bits, candidate.forced_complete_bits_valid,
          candidate.forced_survive_bits, candidate.forced_survive_bits_valid,
          candidate.time_constraints);
      auto it = next_state_index.find(key);
      if (it == next_state_index.end()) {
        const std::size_t idx = next_states_collapsed.size();
        next_state_index.emplace(std::move(key), idx);
        next_states_collapsed.push_back(std::move(candidate));
      } else {
        next_states_collapsed[it->second].weight += candidate.weight;
      }
    };

    std::vector<std::size_t> active_state_indices;
    active_state_indices.reserve(states.size());
    std::vector<ExactScenarioPoint> rank_seed_points;
    rank_seed_points.reserve(states.size());
    std::vector<ExactScenarioPoint> denom_points;
    if (rank_idx > 0u) {
      denom_points.reserve(states.size());
    }

    bool valid_rank = true;
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

      ExactScenarioPoint rank_point;
      rank_point.t = t;
      rank_point.density_index = active_state_indices.size();
      rank_point.weight = 1.0;
      if (trial_params_soa_by_trial != nullptr &&
          static_cast<std::size_t>(trial_slot) < trial_params_soa_by_trial->size()) {
        rank_point.trial_params_soa =
            (*trial_params_soa_by_trial)[static_cast<std::size_t>(trial_slot)];
      }
      rank_point.forced_complete_bits = state.forced_complete_bits;
      rank_point.forced_survive_bits = state.forced_survive_bits;
      rank_point.forced_complete_bits_valid = state.forced_complete_bits_valid;
      rank_point.forced_survive_bits_valid = state.forced_survive_bits_valid;
      rank_point.time_constraints = state.time_constraints;

      if (rank_idx > 0u) {
        const double lower_t =
            times_by_trial[static_cast<std::size_t>(trial_slot)][rank_idx - 1u];
        if (!std::isfinite(lower_t) || lower_t < 0.0) {
          valid_rank = false;
          break;
        }
        ExactScenarioPoint denom_point = rank_point;
        denom_point.t = lower_t;
        denom_points.push_back(std::move(denom_point));
      }

      active_state_indices.push_back(state_idx);
      rank_seed_points.push_back(std::move(rank_point));
    }
    if (!valid_rank) {
      return false;
    }
    if (active_state_indices.empty()) {
      return true;
    }

    std::vector<double> denom(active_state_indices.size(), 1.0);
    if (rank_idx > 0u) {
      uuber::KernelNodeBatchValues denom_values;
      if (!exact_eval_node_batch_from_points(
              ctx, info_node_idx, denom_points, component_idx, trial_params,
              trial_type_key, EvalNeed::kSurvival, denom_values) ||
          denom_values.survival.size() != denom.size()) {
        return false;
      }
      for (std::size_t i = 0; i < denom.size(); ++i) {
        denom[i] = denom_values.survival[i];
      }
    }

    if (ranked_batch_debug_enabled() && rank_idx > 0u &&
        trial_params_soa_by_trial == nullptr) {
      static int denom_mismatch_budget = 10;
      constexpr double kDebugTol = 1e-10;
      for (std::size_t i = 0; i < active_state_indices.size(); ++i) {
        const TrialSequenceState &source_state = states[active_state_indices[i]];
        const double lower_t =
            times_by_trial[static_cast<std::size_t>(source_state.trial_slot)][rank_idx - 1u];
        const double scalar_denom = evaluator_evaluate_survival_with_forced(
            info.node_id,
            source_state.forced_complete_bits_valid
                ? &source_state.forced_complete_bits
                : nullptr,
            source_state.forced_complete_bits_valid,
            source_state.forced_survive_bits_valid
                ? &source_state.forced_survive_bits
                : nullptr,
            source_state.forced_survive_bits_valid, component_idx, lower_t, ctx,
            trial_type_key, trial_params, &source_state.time_constraints,
            nullptr);
        if (std::fabs(scalar_denom - denom[i]) > kDebugTol &&
            denom_mismatch_budget > 0) {
          --denom_mismatch_budget;
          Rcpp::Rcout << "ranked_batch_denom_mismatch rank=" << rank_idx
                      << " trial_slot=" << source_state.trial_slot
                      << " batch=" << denom[i]
                      << " scalar=" << scalar_denom << "\n";
        }
      }
    }

    std::vector<ExactScenarioPoint> scenario_points;
    if (!exact_collect_scenarios_batch_from_points(
            *compiler, ctx, outcome_node_idx, rank_seed_points, component_idx,
            trial_params, trial_type_key, scenario_points)) {
      return true;
    }

    if (ranked_batch_debug_enabled() && trial_params_soa_by_trial == nullptr) {
      static int scenario_mismatch_budget = 10;
      constexpr double kDebugTol = 1e-10;
      std::vector<double> batch_transition_weight(active_state_indices.size(),
                                                  0.0);
      for (const ExactScenarioPoint &point : scenario_points) {
        if (point.density_index < batch_transition_weight.size()) {
          batch_transition_weight[point.density_index] += point.weight;
        }
      }
      for (std::size_t i = 0; i < active_state_indices.size(); ++i) {
        const TrialSequenceState &source_state = states[active_state_indices[i]];
        const double t =
            times_by_trial[static_cast<std::size_t>(source_state.trial_slot)][rank_idx];
        double scalar_transition_weight = 0.0;
        (void)for_each_sequence_node_transition(
            *compiler, ctx, outcome_node_idx, t, component_idx, trial_params,
            trial_type_key, &source_state.time_constraints, nullptr,
            source_state.forced_complete_bits_valid
                ? &source_state.forced_complete_bits
                : nullptr,
            source_state.forced_complete_bits_valid,
            source_state.forced_survive_bits_valid
                ? &source_state.forced_survive_bits
                : nullptr,
            source_state.forced_survive_bits_valid,
            [&](double scenario_weight, uuber::BitsetState &&,
                bool, uuber::BitsetState &&, bool, TimeConstraintMap &&) {
              if (std::isfinite(scenario_weight) && scenario_weight > 0.0) {
                scalar_transition_weight += scenario_weight;
              }
            });
        if (std::fabs(scalar_transition_weight - batch_transition_weight[i]) >
                kDebugTol &&
            scenario_mismatch_budget > 0) {
          --scenario_mismatch_budget;
          Rcpp::Rcout << "ranked_batch_transition_mismatch rank=" << rank_idx
                      << " trial_slot=" << source_state.trial_slot
                      << " batch=" << batch_transition_weight[i]
                      << " scalar=" << scalar_transition_weight << "\n";
        }
      }
    }

    std::vector<double> scenario_survival(scenario_points.size(), 1.0);
    if (!competitors.empty()) {
      const uuber::CompetitorClusterCacheEntry &competitor_cache =
          fetch_competitor_cluster_cache(ctx, competitors);
      if (!exact_competitor_survival_batch(
              ctx, competitor_cache, component_idx, trial_params,
              trial_type_key, scenario_points, scenario_survival) ||
          scenario_survival.size() != scenario_points.size()) {
        return false;
      }
    }

    for (std::size_t point_idx = 0; point_idx < scenario_points.size();
         ++point_idx) {
      const ExactScenarioPoint &point = scenario_points[point_idx];
      if (point.density_index >= active_state_indices.size()) {
        return false;
      }
      const TrialSequenceState &source_state =
          states[active_state_indices[point.density_index]];
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

      if (!competitors.empty()) {
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

      if (!competitors.empty() && !persistent_sources.empty()) {
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

      accumulate_next_state(std::move(next_state));
    }

    states.clear();
    states.reserve(next_states_collapsed.size());
    for (TrialSequenceState &next_state : next_states_collapsed) {
      if (std::isfinite(next_state.weight) && next_state.weight > kBranchEps) {
        states.push_back(std::move(next_state));
      }
    }
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

  if (ranked_batch_debug_enabled() && trial_params_soa_by_trial == nullptr) {
    static int mismatch_budget = 10;
    constexpr double kTol = 1e-10;
    for (std::size_t trial_slot = 0; trial_slot < trial_count; ++trial_slot) {
      const double scalar = sequence_prefix_density_resolved(
          ctx, outcome_indices, node_indices, times_by_trial[trial_slot],
          component_idx, trial_params, trial_type_key, nullptr, compiler,
          step_competitor_ids_ptrs, step_persistent_sources);
      if (std::fabs(scalar - density_out[trial_slot]) > kTol &&
          mismatch_budget > 0) {
        --mismatch_budget;
        Rcpp::Rcout << "ranked_batch_mismatch trial_slot=" << trial_slot
                    << " batch=" << density_out[trial_slot]
                    << " scalar=" << scalar << " times=[";
        for (std::size_t rank_i = 0; rank_i < rank_count; ++rank_i) {
          if (rank_i > 0u) {
            Rcpp::Rcout << ",";
          }
          Rcpp::Rcout << times_by_trial[trial_slot][rank_i];
        }
        Rcpp::Rcout << "]\n";
      }
    }
  }

  return true;
}
