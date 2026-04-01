#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "context.h"
#include "evaluator_internal.h"
#include "trial_params.h"

enum class RankedProgramEvalKind : std::uint8_t {
  Density = 0,
  CDF = 1,
  Survival = 2
};

struct RankedProgramEvalRef {
  RankedProgramEvalKind kind{RankedProgramEvalKind::Density};
  int node_idx{-1};
  int slot{-1};
};

struct RankedProgramSliceRef {
  std::vector<RankedProgramEvalRef> evals;
  std::vector<int> node_indices;
  std::vector<int> relevant_source_ids;

  bool empty() const noexcept { return evals.empty(); }
};

inline bool ranked_program_eval_ref_same(
    const RankedProgramEvalRef &lhs,
    const RankedProgramEvalRef &rhs) noexcept {
  return lhs.kind == rhs.kind && lhs.node_idx == rhs.node_idx &&
         lhs.slot == rhs.slot;
}

inline std::uint64_t ranked_program_eval_ref_hash(
    const RankedProgramEvalRef &ref) noexcept {
  std::uint64_t hash = kFNV64Offset;
  const std::uint8_t kind = static_cast<std::uint8_t>(ref.kind);
  hash_append_bytes(hash, &kind, sizeof(kind));
  hash_append_bytes(hash, &ref.node_idx, sizeof(ref.node_idx));
  hash_append_bytes(hash, &ref.slot, sizeof(ref.slot));
  return hash;
}

inline bool ranked_program_slice_same(const RankedProgramSliceRef &lhs,
                                      const RankedProgramSliceRef &rhs) noexcept {
  if (lhs.evals.size() != rhs.evals.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.evals.size(); ++i) {
    if (!ranked_program_eval_ref_same(lhs.evals[i], rhs.evals[i])) {
      return false;
    }
  }
  return true;
}

inline std::uint64_t
ranked_program_slice_hash(const RankedProgramSliceRef &slice) noexcept {
  std::uint64_t hash = kFNV64Offset;
  const std::size_t eval_count = slice.evals.size();
  hash_append_bytes(hash, &eval_count, sizeof(eval_count));
  for (const RankedProgramEvalRef &eval_ref : slice.evals) {
    const std::uint64_t eval_hash = ranked_program_eval_ref_hash(eval_ref);
    hash_append_bytes(hash, &eval_hash, sizeof(eval_hash));
  }
  return hash;
}

enum class RankedStateDeltaKind : std::uint8_t {
  CompleteSources = 0,
  SurviveSources = 1,
  OrWitnessSurvive = 2
};

struct RankedStateDelta {
  RankedStateDeltaKind kind{RankedStateDeltaKind::CompleteSources};
  int trigger_bit{-1};
  int source_mask_begin{-1};
  int source_mask_count{0};
  int invalidate_slot{0};
  bool bind_exact_current_time{false};
  std::vector<int> source_ids;
  std::vector<int> source_bits;
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

inline std::uint64_t sequence_state_hash_bitset(
    const uuber::BitsetState *bits, bool bits_valid) {
  if (!bits_valid || bits == nullptr) {
    return 0ULL;
  }
  std::uint64_t hash = kFNV64Offset;
  const std::int32_t bit_count = static_cast<std::int32_t>(bits->bit_count());
  hash_append_bytes(hash, &bit_count, sizeof(bit_count));
  if (bits->words().empty()) {
    const std::uint64_t word = bits->small_word();
    hash_append_bytes(hash, &word, sizeof(word));
  } else {
    for (std::uint64_t word : bits->words()) {
      hash_append_bytes(hash, &word, sizeof(word));
    }
  }
  return hash;
}

inline std::uint64_t
sequence_state_hash_time_constraints(const TimeConstraintMap *time_constraints) {
  if (time_constraints == nullptr) {
    return 0ULL;
  }
  std::uint64_t hash = kFNV64Offset;
  for (const auto &kv : *time_constraints) {
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
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const TimeConstraintMap *time_constraints) {
  SequenceStateKey key;
  key.forced_complete_bits_valid = forced_complete_bits_valid;
  key.forced_survive_bits_valid = forced_survive_bits_valid;
  key.forced_complete_bit_count =
      (forced_complete_bits_valid && forced_complete_bits != nullptr)
          ? forced_complete_bits->bit_count()
          : 0;
  key.forced_survive_bit_count =
      (forced_survive_bits_valid && forced_survive_bits != nullptr)
          ? forced_survive_bits->bit_count()
          : 0;
  key.forced_complete_hash =
      sequence_state_hash_bitset(forced_complete_bits,
                                 forced_complete_bits_valid);
  key.forced_survive_hash =
      sequence_state_hash_bitset(forced_survive_bits,
                                 forced_survive_bits_valid);
  key.time_constraint_hash =
      sequence_state_hash_time_constraints(time_constraints);
  key.time_constraint_size =
      time_constraints != nullptr ? time_constraints->size() : 0u;
  return key;
}

inline SequenceStateKey sequence_state_key(
    const uuber::BitsetState &forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState &forced_survive_bits,
    bool forced_survive_bits_valid,
    const TimeConstraintMap &time_constraints) {
  return sequence_state_key(&forced_complete_bits, forced_complete_bits_valid,
                            &forced_survive_bits, forced_survive_bits_valid,
                            &time_constraints);
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

inline bool apply_ranked_state_delta_raw(
    double t, const RankedStateDelta &delta,
    const uuber::BitsetState &forced_complete_bits,
    bool forced_complete_bits_valid, uuber::BitsetState &forced_survive_bits,
    bool &forced_survive_bits_valid, TimeConstraintMap &time_constraints) {
  switch (delta.kind) {
  case RankedStateDeltaKind::CompleteSources:
    for (int id : delta.source_ids) {
      if (!time_constraints_mark_complete(id, t, delta.bind_exact_current_time,
                                          time_constraints)) {
        return false;
      }
    }
    return true;
  case RankedStateDeltaKind::SurviveSources:
    for (int id : delta.source_ids) {
      if (!time_constraints_mark_survive(id, t, time_constraints)) {
        return false;
      }
    }
    return true;
  case RankedStateDeltaKind::OrWitnessSurvive: {
    int witness = NA_INTEGER;
    bool all_forced = true;
    for (std::size_t src_i = 0; src_i < delta.source_ids.size(); ++src_i) {
      const int id = delta.source_ids[src_i];
      if (id == NA_INTEGER || id < 0) {
        continue;
      }
      const int bit_idx =
          (src_i < delta.source_bits.size()) ? delta.source_bits[src_i] : -1;
      const bool is_forced =
          (forced_complete_bits_valid && bit_idx >= 0 &&
           forced_complete_bits.test(bit_idx)) ||
          time_constraints_contains_complete_at(&time_constraints, id, t);
      if (!is_forced) {
        witness = id;
        all_forced = false;
        break;
      }
    }
    return !all_forced && witness != NA_INTEGER &&
           time_constraints_mark_survive(witness, t, time_constraints);
  }
  }
  return false;
}

struct RankedBatchPlan {
  std::vector<RankedStateDelta> deltas;
  RankedProgramSliceRef slice;
  bool deterministic{false};
  bool valid{false};
};

struct RankedBatchPlanGroup {
  std::size_t leader_index{0};
  std::vector<std::size_t> plan_indices;
  bool has_deltas{false};
};

struct RankedNodeBatchPlan {
  bool compiled{false};
  bool compiling{false};
  bool valid{false};
  bool has_shared_slice_groups{false};
  std::vector<RankedBatchPlan> batch_plans;
  std::vector<RankedBatchPlanGroup> batch_plan_groups;
};

class RankedBatchPlanner {
public:
  explicit RankedBatchPlanner(const uuber::NativeContext &ctx);

  const RankedNodeBatchPlan &plan_for_node(int node_idx);

private:
  void compile_node(int node_idx);

  const uuber::NativeContext &ctx;
  std::vector<RankedNodeBatchPlan> plans;
};

RankedBatchPlanner &ranked_transition_planner_for_ctx(
    const uuber::NativeContext &ctx);
void clear_ranked_transition_runtime_caches(
    std::uint64_t runtime_cache_instance_id) noexcept;

bool ranked_batch_plan_supports_deterministic_fastpath(
    const RankedBatchPlan &plan);

void compile_ranked_node_batch_plan(
    const uuber::NativeContext &ctx, int node_idx,
    const std::function<const RankedNodeBatchPlan &(int)> &plan_lookup,
    RankedNodeBatchPlan &plan);

std::vector<int>
collect_competitor_sources(const uuber::NativeContext &ctx,
                           const std::vector<int> &competitor_ids);

bool sequence_prefix_density_batch_resolved(
    const uuber::NativeContext &ctx, const std::vector<int> &outcome_indices,
    const std::vector<int> &node_indices,
    const std::vector<const double *> &times_by_trial, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::vector<double> &density_out,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_by_trial =
        nullptr,
    RankedBatchPlanner *transition_planner = nullptr,
    const std::vector<const std::vector<int> *> *step_competitor_ids_ptrs =
        nullptr,
    const std::vector<std::vector<int>> *step_persistent_sources = nullptr);
