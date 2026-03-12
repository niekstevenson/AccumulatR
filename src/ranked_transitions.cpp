#include "ranked_transitions.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "evaluator_internal.h"

namespace {

struct SequenceState {
  double weight{0.0};
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  ExactSourceTimeMap exact_source_times;
  SourceTimeBoundsMap source_time_bounds;
};

struct SequenceStateKey {
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  int forced_complete_bit_count{0};
  int forced_survive_bit_count{0};
  std::uint64_t forced_complete_hash{0};
  std::uint64_t forced_survive_hash{0};
  std::uint64_t exact_hash{0};
  std::uint64_t bounds_hash{0};
  std::size_t exact_size{0};
  std::size_t bounds_size{0};

  bool operator==(const SequenceStateKey &other) const noexcept {
    return forced_complete_bits_valid == other.forced_complete_bits_valid &&
           forced_survive_bits_valid == other.forced_survive_bits_valid &&
           forced_complete_bit_count == other.forced_complete_bit_count &&
           forced_survive_bit_count == other.forced_survive_bit_count &&
           forced_complete_hash == other.forced_complete_hash &&
           forced_survive_hash == other.forced_survive_hash &&
           exact_hash == other.exact_hash && bounds_hash == other.bounds_hash &&
           exact_size == other.exact_size && bounds_size == other.bounds_size;
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
    combine(seed, static_cast<std::size_t>(key.exact_hash));
    combine(seed, static_cast<std::size_t>(key.bounds_hash));
    combine(seed, static_cast<std::size_t>(key.forced_complete_bits_valid));
    combine(seed, static_cast<std::size_t>(key.forced_survive_bits_valid));
    combine(seed, static_cast<std::size_t>(key.forced_complete_bit_count));
    combine(seed, static_cast<std::size_t>(key.forced_survive_bit_count));
    combine(seed, key.exact_size);
    combine(seed, key.bounds_size);
    return seed;
  }
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

inline std::uint64_t ranked_hash_exact_times(
    const ExactSourceTimeMap &exact_source_times) {
  std::uint64_t hash = kFNV64Offset;
  for (const auto &kv : exact_source_times) {
    const std::int32_t id = static_cast<std::int32_t>(kv.first);
    hash_append_bytes(hash, &id, sizeof(id));
    const std::uint64_t bits = canonical_double_bits(kv.second);
    hash_append_bytes(hash, &bits, sizeof(bits));
  }
  return hash;
}

inline std::uint64_t ranked_hash_bounds(
    const SourceTimeBoundsMap &source_time_bounds) {
  std::uint64_t hash = kFNV64Offset;
  for (const auto &kv : source_time_bounds) {
    const std::int32_t id = static_cast<std::int32_t>(kv.first);
    hash_append_bytes(hash, &id, sizeof(id));
    const std::uint64_t lower_bits = canonical_double_bits(kv.second.first);
    const std::uint64_t upper_bits = canonical_double_bits(kv.second.second);
    hash_append_bytes(hash, &lower_bits, sizeof(lower_bits));
    hash_append_bytes(hash, &upper_bits, sizeof(upper_bits));
  }
  return hash;
}

inline SequenceStateKey sequence_state_key(
    const uuber::BitsetState &forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState &forced_survive_bits,
    bool forced_survive_bits_valid,
    const ExactSourceTimeMap &exact_source_times,
    const SourceTimeBoundsMap &source_time_bounds) {
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
  key.exact_hash = ranked_hash_exact_times(exact_source_times);
  key.bounds_hash = ranked_hash_bounds(source_time_bounds);
  key.exact_size = exact_source_times.size();
  key.bounds_size = source_time_bounds.size();
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

inline void ranked_add_source_id(const uuber::NativeContext &ctx, int source_id,
                                 int source_bit_idx,
                                 uuber::BitsetState &bits_out,
                                 bool &bits_valid) {
  if (source_id == NA_INTEGER || source_id < 0) {
    return;
  }
  if (source_bit_idx < 0) {
    set_forced_id_bit_strict(ctx, source_id, bits_out, bits_valid);
    return;
  }
  ensure_forced_bitset_capacity(ctx, bits_out, bits_valid);
  if (!bits_valid) {
    Rcpp::stop("IR forced-state bitset unavailable for ranked source id %d",
               source_id);
  }
  bits_out.set(source_bit_idx);
}

template <typename EmitFn>
bool for_each_sequence_node_transition(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, double t, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    uuber::KernelRuntimeState *kernel_runtime,
    const uuber::BitsetState *base_forced_complete_bits_in,
    bool base_forced_complete_bits_in_valid,
    const uuber::BitsetState *base_forced_survive_bits_in,
    bool base_forced_survive_bits_in_valid, EmitFn &&emit) {
  const RankedNodeTransitionPlan &plan = compiler.plan_for_node(node_idx);
  if (!plan.valid || plan.transitions.empty()) {
    return false;
  }

  const uuber::TrialParamsSoA *trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);
  IntegrationSettings settings;
  auto apply_ranked_transition_mask =
      [&](const RankedTransitionStep &step, uuber::BitsetState &bits,
          bool &bits_valid) -> bool {
    if (step.source_mask_begin < 0 || step.source_mask_count <= 0) {
      return false;
    }
    apply_transition_mask_words(ctx, step.source_mask_begin,
                                step.source_mask_count, step.invalidate_slot,
                                bits, bits_valid, kernel_runtime);
    return true;
  };
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
    NodeEvalState transition_state(
        ctx, t, component_idx, trial_params, trial_type_key, false, -1,
        exact_source_times, source_time_bounds, nullptr,
        base_forced_complete_bits_valid ? &base_forced_complete_bits : nullptr,
        base_forced_complete_bits_valid,
        base_forced_survive_bits_valid ? &base_forced_survive_bits : nullptr,
        base_forced_survive_bits_valid, kernel_runtime);
    invalidate_kernel_runtime_root(transition_state);

    bool valid = true;
    for (const RankedTransitionStep &step : transition.steps) {
      if (step.kind == RankedTransitionStepKind::EvalDensityNode) {
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
        NodeEvalResult eval = evaluator_eval_node_recursive_dense(
            step.node_idx, transition_state, EvalNeed::kCDF);
        const double Fj = clamp_probability(eval.cdf);
        if (!std::isfinite(Fj) || Fj <= 0.0) {
          valid = false;
          break;
        }
        weight *= Fj;
      } else if (step.kind == RankedTransitionStepKind::EvalGuardEffective) {
        GuardEvalInput guard_input = make_guard_input_forced_state(
            ctx, step.node_idx, component_idx, &trial_type_key, trial_params,
            trial_params_soa, exact_source_times, source_time_bounds,
            transition_state.forced_state);
        const double eff = evaluator_guard_effective_survival_internal(
            guard_input, t, settings);
        if (!std::isfinite(eff) || eff <= 0.0) {
          valid = false;
          break;
        }
        weight *= eff;
      } else if (step.kind == RankedTransitionStepKind::AddCompleteSources) {
        const bool applied_mask = apply_ranked_transition_mask(
            step, transition_state.forced_complete_bits,
            transition_state.forced_complete_bits_valid);
        if (applied_mask && step.source_mask_covers_ids) {
          continue;
        }
        bool mutated_bits = false;
        for (std::size_t i = 0; i < step.source_ids.size(); ++i) {
          const int id = step.source_ids[i];
          const int bit_idx =
              (i < step.source_bits.size()) ? step.source_bits[i] : -1;
          ranked_add_source_id(ctx, id, bit_idx,
                               transition_state.forced_complete_bits,
                               transition_state.forced_complete_bits_valid);
          mutated_bits = true;
        }
        if (mutated_bits) {
          invalidate_kernel_runtime_root(transition_state);
        }
      } else if (step.kind == RankedTransitionStepKind::AddSurviveSources) {
        const bool applied_mask = apply_ranked_transition_mask(
            step, transition_state.forced_survive_bits,
            transition_state.forced_survive_bits_valid);
        if (applied_mask && step.source_mask_covers_ids) {
          continue;
        }
        bool mutated_bits = false;
        for (std::size_t i = 0; i < step.source_ids.size(); ++i) {
          const int id = step.source_ids[i];
          const int bit_idx =
              (i < step.source_bits.size()) ? step.source_bits[i] : -1;
          ranked_add_source_id(ctx, id, bit_idx,
                               transition_state.forced_survive_bits,
                               transition_state.forced_survive_bits_valid);
          mutated_bits = true;
        }
        if (mutated_bits) {
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
              transition_state.forced_complete_bits_valid &&
              transition_state.forced_complete_bits.test(bit_to_check);
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
        ranked_add_source_id(ctx, witness, witness_bit,
                             transition_state.forced_survive_bits,
                             transition_state.forced_survive_bits_valid);
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
         transition_state.forced_survive_bits_valid);
    emitted_any = true;
  }
  return emitted_any;
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
        child_sources.push_back(ensure_source_ids(ctx, child_idx));
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
        child_sources.push_back(ensure_source_ids(ctx, child_idx));
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
  std::vector<int> onset_source_ids;
  onset_source_ids.reserve(ctx.accumulators.size());
  const std::vector<TrialAccumulatorParams> *acc_params =
      trial_params ? &trial_params->acc_params : nullptr;
  for (std::size_t acc_i = 0; acc_i < ctx.accumulators.size(); ++acc_i) {
    const uuber::NativeAccumulator &acc = ctx.accumulators[acc_i];
    const TrialAccumulatorParams *override =
        (acc_params && acc_i < acc_params->size()) ? &((*acc_params)[acc_i])
                                                   : nullptr;
    int onset_kind = override ? override->onset_kind : acc.onset_kind;
    if (onset_kind == uuber::ONSET_AFTER_ACCUMULATOR) {
      int src_idx =
          override ? override->onset_source_acc_idx : acc.onset_source_acc_idx;
      if (src_idx >= 0 && src_idx < static_cast<int>(ctx.accumulators.size())) {
        int src_id = accumulator_label_id_of(ctx, src_idx);
        if (src_id >= 0 && src_id != NA_INTEGER) {
          onset_source_ids.push_back(src_id);
        }
      }
    } else if (onset_kind == uuber::ONSET_AFTER_POOL) {
      int src_idx = override ? override->onset_source_pool_idx
                             : acc.onset_source_pool_idx;
      if (src_idx >= 0 && src_idx < static_cast<int>(ctx.pools.size())) {
        int src_id = pool_label_id_of(ctx, src_idx);
        if (src_id >= 0 && src_id != NA_INTEGER) {
          onset_source_ids.push_back(src_id);
        }
      }
    }
  }
  sort_unique(onset_source_ids);
  auto onset_source_contains = [&](int src_id) -> bool {
    return std::binary_search(onset_source_ids.begin(), onset_source_ids.end(),
                              src_id);
  };

  std::vector<SequenceState> states;
  SequenceState init_state;
  init_state.weight = 1.0;
  (void)ensure_forced_bitset_capacity(ctx, init_state.forced_complete_bits,
                                      init_state.forced_complete_bits_valid);
  (void)ensure_forced_bitset_capacity(ctx, init_state.forced_survive_bits,
                                      init_state.forced_survive_bits_valid);
  states.push_back(std::move(init_state));
  ExactSourceTimeMap prefix_exact;

  for (std::size_t rank_idx = 0; rank_idx < rank_count; ++rank_idx) {
    if (kernel_runtime) {
      uuber::invalidate_kernel_runtime_from_slot(*kernel_runtime, 0);
    }
    double t = times[rank_idx];
    if (!std::isfinite(t) || t < 0.0) {
      return 0.0;
    }

    int outcome_idx = outcome_indices[rank_idx];
    std::vector<int> rank_onset_sources;
    if (!onset_source_ids.empty()) {
      std::vector<int> src_ids = ensure_source_ids(ctx, node_indices[rank_idx]);
      rank_onset_sources.reserve(src_ids.size());
      for (int src_id : src_ids) {
        if (onset_source_contains(src_id)) {
          rank_onset_sources.push_back(src_id);
        }
      }
      sort_unique(rank_onset_sources);
    }
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
          candidate.exact_source_times, candidate.source_time_bounds);
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
      ExactSourceTimeMap merged_exact = state.exact_source_times;
      if (!prefix_exact.empty()) {
        for (const auto &kv : prefix_exact) {
          if (std::isfinite(kv.second)) {
            merged_exact[kv.first] = kv.second;
          }
        }
      }
      SourceTimeBoundsMap merged_bounds = state.source_time_bounds;
      if (!merged_exact.empty() && !merged_bounds.empty()) {
        for (const auto &kv : merged_exact) {
          merged_bounds.erase(kv.first);
        }
      }
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
            trial_type_key, trial_params, &merged_exact, &merged_bounds,
            kernel_runtime);
        if (!std::isfinite(denom) || denom <= kBranchEps) {
          continue;
        }
      }
      int outcome_node_idx = node_indices[rank_idx];
      const bool has_transition = for_each_sequence_node_transition(
          *compiler, ctx, outcome_node_idx, t, component_idx, trial_params,
          trial_type_key, &merged_exact, &merged_bounds, kernel_runtime,
          forced_complete_bits_valid ? &forced_complete_bits : nullptr,
          forced_complete_bits_valid,
          forced_survive_bits_valid ? &forced_survive_bits : nullptr,
          forced_survive_bits_valid,
          [&](double scenario_weight,
              uuber::BitsetState &&scenario_complete_bits,
              bool scenario_complete_bits_valid,
              uuber::BitsetState &&scenario_survive_bits,
              bool scenario_survive_bits_valid) {
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
                  &merged_exact, &merged_bounds, kernel_runtime);
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
                }
              }
            }
            SequenceState next_state;
            next_state.weight = weight;
            next_state.forced_complete_bits = std::move(next_complete_bits);
            next_state.forced_survive_bits = std::move(next_survive_bits);
            next_state.forced_complete_bits_valid = next_complete_bits_valid;
            next_state.forced_survive_bits_valid = next_survive_bits_valid;
            next_state.exact_source_times = state.exact_source_times;
            next_state.source_time_bounds = state.source_time_bounds;
            if (!onset_source_ids.empty()) {
              bool bounds_ok = true;
              for (int src_id : onset_source_ids) {
                if (!forced_bits_contains_label_id_strict(
                        ctx, src_id, next_state.forced_survive_bits,
                        next_state.forced_survive_bits_valid)) {
                  continue;
                }
                if (next_state.exact_source_times.find(src_id) !=
                    next_state.exact_source_times.end()) {
                  continue;
                }
                if (prefix_exact.find(src_id) != prefix_exact.end()) {
                  continue;
                }
                auto bound_it = next_state.source_time_bounds.find(src_id);
                if (bound_it == next_state.source_time_bounds.end()) {
                  bound_it = next_state.source_time_bounds
                                 .emplace(src_id,
                                          std::make_pair(
                                              0.0,
                                              std::numeric_limits<double>::infinity()))
                                 .first;
                }
                auto &bound = bound_it->second;
                bound.first = std::max(bound.first, t);
                if (!(bound.second > bound.first)) {
                  bounds_ok = false;
                  break;
                }
              }
              if (!bounds_ok) {
                return;
              }
              for (int src_id : onset_source_ids) {
                if (!forced_bits_contains_label_id_strict(
                        ctx, src_id, next_state.forced_complete_bits,
                        next_state.forced_complete_bits_valid)) {
                  continue;
                }
                if (next_state.exact_source_times.find(src_id) !=
                    next_state.exact_source_times.end()) {
                  continue;
                }
                if (prefix_exact.find(src_id) != prefix_exact.end()) {
                  continue;
                }
                auto bound_it = next_state.source_time_bounds.find(src_id);
                if (bound_it == next_state.source_time_bounds.end()) {
                  bound_it = next_state.source_time_bounds
                                 .emplace(src_id,
                                          std::make_pair(
                                              0.0,
                                              std::numeric_limits<double>::infinity()))
                                 .first;
                }
                auto &bound = bound_it->second;
                bound.second = std::min(bound.second, t);
                if (!(bound.second > bound.first)) {
                  bounds_ok = false;
                  break;
                }
              }
              if (!bounds_ok) {
                return;
              }
            }
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

    if (!rank_onset_sources.empty()) {
      for (int src_id : rank_onset_sources) {
        prefix_exact[src_id] = t;
      }
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
