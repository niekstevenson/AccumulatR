#pragma once

#include <Rcpp.h>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "accumulator.h"

namespace uuber {

enum OnsetKind : int {
  ONSET_ABSOLUTE = 0,
  ONSET_AFTER_ACCUMULATOR = 1,
  ONSET_AFTER_POOL = 2
};

enum NodeKind : int {
  NODE_UNKNOWN = 0,
  NODE_EVENT = 1,
  NODE_AND = 2,
  NODE_OR = 3,
  NODE_NOT = 4,
  NODE_GUARD = 5
};

enum class IrNodeOp : std::uint8_t {
  EventAcc = 0,
  EventPool = 1,
  And = 2,
  Or = 3,
  Not = 4,
  Guard = 5
};

enum class IrOutcomeCouplingKind : std::uint8_t {
  None = 0,
  Pair = 1,
  NWay = 2,
  GenericNodeIntegral = 3
};

enum IrNodeFlags : std::uint32_t {
  IR_NODE_FLAG_NONE = 0u,
  IR_NODE_FLAG_NEEDS_FORCED = 1u << 0,
  IR_NODE_FLAG_SCENARIO_SENSITIVE = 1u << 1,
  IR_NODE_FLAG_SPECIAL_DEADLINE = 1u << 2,
  IR_NODE_FLAG_SPECIAL_GUESS = 1u << 3
};

struct NativeAccumulator {
  std::string id;
  std::string dist;
  double onset{};
  int onset_kind{ONSET_ABSOLUTE};
  int onset_source_acc_idx{-1};
  int onset_source_pool_idx{-1};
  double onset_lag{0.0};
  double q{};
  std::vector<std::string> components;
  std::vector<int> component_indices;
  std::string shared_trigger_id;
  AccDistParams dist_cfg;
};

struct LabelRef {
  int label_id{-1};
  int acc_idx{-1};
  int pool_idx{-1};
  int outcome_idx{-1};
};

struct NativePool {
  std::string id;
  int k{};
  std::vector<std::string> members;
  std::vector<LabelRef> member_refs;
};

struct IrNode {
  IrNodeOp op{IrNodeOp::And};
  int child_begin{-1};
  int child_count{0};
  int source_id_begin{-1};
  int source_id_count{0};
  int source_mask_begin{-1};
  int source_mask_count{0};
  int component_mask_offset{-1};
  int reference_idx{-1};
  int blocker_idx{-1};
  std::uint32_t flags{IR_NODE_FLAG_NONE};
  int event_idx{-1};
  int node_id{-1};
};

struct IrEvent {
  int node_idx{-1};
  int acc_idx{-1};
  int pool_idx{-1};
  int label_id{-1};
  int outcome_idx{-1};
  int component_mask_offset{-1};
};

struct IrOutcomeGuessDonor {
  int outcome_idx{-1};
  double weight{0.0};
  int rt_policy_code{0}; // 0=keep, 1=na, 2=drop, 3=unknown
};

struct IrOutcomeCouplingOp {
  IrOutcomeCouplingKind kind{IrOutcomeCouplingKind::None};
  int node_idx{-1};
  int competitor_begin{-1};
  int competitor_count{0};
  bool requires_exact_scenario_eval{false};
  int gate_event_idx{-1};
  int target_event_idx{-1};
  int aux_event_begin{-1};
  int aux_event_count{0};
};

struct IrOutcome {
  int label_id{-1};
  int node_idx{-1};
  int competitor_begin{-1};
  int competitor_count{0};
  int allowed_component_mask_offset{-1};
  bool maps_to_na{false};
  int alias_begin{-1};
  int alias_count{0};
  int guess_begin{-1};
  int guess_count{0};
  int coupling_op_idx{-1};
};

struct IrContext {
  std::vector<IrNode> nodes;
  std::vector<IrEvent> events;
  std::vector<IrOutcome> outcomes;
  std::vector<int> node_children;
  std::vector<int> node_source_label_ids;
  std::vector<std::uint64_t> node_source_masks;
  std::vector<int> outcome_competitors;
  std::vector<std::uint64_t> component_masks;
  std::vector<int> outcome_alias_sources;
  std::vector<IrOutcomeGuessDonor> outcome_guess_donors;
  std::vector<IrOutcomeCouplingOp> outcome_coupling_ops;
  std::vector<int> outcome_coupling_event_indices;
  std::unordered_map<std::uint64_t, int> outcome_coupling_lookup;
  std::unordered_map<int, std::vector<int>> label_id_to_outcomes;
  std::unordered_map<int, std::vector<int>> node_idx_to_outcomes;
  std::unordered_map<int, int> id_to_node_idx;
  std::unordered_map<int, int> label_id_to_bit_idx;
  int source_mask_words{0};
  int component_mask_words{0};
  bool valid{false};
};

enum class KernelOpCode : std::uint8_t {
  Event = 0,
  And = 1,
  Or = 2,
  Not = 3,
  Guard = 4
};

struct KernelOp {
  KernelOpCode code{KernelOpCode::Event};
  int node_idx{-1};
  int out_slot{-1};
  int event_idx{-1};
  int child_begin{-1};
  int child_count{0};
  int ref_slot{-1};
  int blocker_slot{-1};
  std::uint32_t flags{IR_NODE_FLAG_NONE};
};

struct KernelOutputMap {
  std::vector<int> node_idx_to_slot;
  std::vector<int> slot_to_node_idx;
  std::vector<int> outcome_idx_to_slot;
};

struct KernelProgram {
  std::vector<KernelOp> ops;
  std::vector<int> children;
  KernelOutputMap outputs;
  int max_child_count{0};
  bool has_guard{false};
  bool valid{false};
};

struct KernelStateTransition {
  int trigger_bit{-1};
  int acc_idx{-1};
  int op_begin{-1};
  int op_count{0};
};

enum class KernelGuardEvalMode : std::uint8_t {
  LinearChainODE = 0
};

struct KernelGuardTransition {
  int node_idx{-1};
  int source_mask_begin{-1};
  int source_mask_count{0};
  int invalidate_slot{0};
  int reference_node_idx{-1};
  int blocker_node_idx{-1};
  int reference_slot{-1};
  int blocker_slot{-1};
  KernelGuardEvalMode eval_mode{KernelGuardEvalMode::LinearChainODE};
  int linear_chain_begin{-1};
  int linear_chain_count{0};
  int linear_chain_leaf_idx{-1};
};

struct KernelStateGraph {
  int forced_bit_count{0};
  std::vector<int> bit_idx_to_label_id;
  std::vector<KernelStateTransition> trigger_transitions;
  std::vector<int> trigger_transition_begin;
  std::vector<int> trigger_transition_count;
  std::vector<int> trigger_op_indices;
  std::vector<KernelGuardTransition> guard_transitions;
  std::vector<int> guard_linear_chain_nodes;
  std::vector<int> node_guard_transition_idx;
  std::vector<std::uint8_t> node_contains_guard;
  std::vector<int> node_competitor_guard_transition_idx;
  std::vector<int> node_competitor_transition_mask_begin;
  std::vector<int> node_competitor_transition_mask_count;
  std::vector<int> node_competitor_transition_invalidate_slot;
  bool valid{false};
};

struct CompetitorCompiledOp {
  std::vector<int> target_node_indices;
  int transition_guard_idx{-1};
  int transition_mask_begin{-1};
  int transition_mask_count{0};
  int transition_invalidate_slot{0};
};

struct CompetitorClusterCacheEntry {
  std::vector<CompetitorCompiledOp> compiled_ops;
  bool mutates_forced_survive{false};
};

struct CompetitorCacheRecord {
  std::vector<int> competitor_ids;
  CompetitorClusterCacheEntry entry;
};

using CompetitorCacheMap =
    std::unordered_map<std::uint64_t, std::deque<CompetitorCacheRecord>>;

struct TrialParamsSoA {
  int n_acc{0};
  std::vector<int> dist_code;
  std::vector<double> onset;
  std::vector<double> q;
  std::vector<double> t0;
  std::vector<double> p1;
  std::vector<double> p2;
  std::vector<double> p3;
  std::vector<double> p4;
  std::vector<double> p5;
  std::vector<double> p6;
  std::vector<double> p7;
  std::vector<double> p8;
  bool valid{false};
};

struct ComponentGuessPolicy {
  std::string target;
  int target_outcome_idx{-1};
  int target_label_id{NA_INTEGER};
  bool target_is_guess{false};
  std::vector<std::pair<int, double>> keep_weights;
  double keep_weight_na{1.0};
  bool has_keep_weight_na{false};
  bool valid{false};
};

struct ComponentContextInfo {
  ComponentGuessPolicy guess;
};

struct OutcomeGuessDonor {
  int outcome_idx{-1};
  double weight{0.0};
  std::string rt_policy;
};

struct OutcomeContextInfo {
  int node_id{-1};
  std::vector<int> competitor_ids;
  std::vector<OutcomeGuessDonor> guess_donors;
  std::vector<int> alias_sources;
  bool maps_to_na{false};
  std::vector<std::string> allowed_components;
};

struct ComponentMap {
  std::vector<std::string> ids;
  std::vector<int> leader_idx;
  std::vector<double> base_weights;
};

struct SharedTriggerGroup {
  std::vector<int> acc_indices;
  int q_acc_idx{-1};
};

struct NAMapCacheKey {
  std::uint64_t hash1{0ULL};
  std::uint64_t hash2{0ULL};

  bool operator==(const NAMapCacheKey &other) const noexcept {
    return hash1 == other.hash1 && hash2 == other.hash2;
  }
};

struct NAMapCacheKeyHash {
  std::size_t operator()(const NAMapCacheKey &key) const noexcept {
    std::uint64_t x = key.hash1 ^ (key.hash2 + 0x9e3779b97f4a7c15ULL +
                                   (key.hash1 << 6) + (key.hash1 >> 2));
    return static_cast<std::size_t>(x);
  }
};

struct NativeContext {
  std::vector<NativeAccumulator> accumulators;
  std::vector<NativePool> pools;
  std::unordered_map<std::string, int> accumulator_index;
  std::unordered_map<std::string, int> pool_index;
  std::unordered_map<std::string, int> label_to_id;
  mutable std::unordered_map<NAMapCacheKey, double, NAMapCacheKeyHash>
      na_map_cache;
  mutable std::unordered_map<std::string, std::deque<std::string>>
      na_cache_order;
  mutable CompetitorCacheMap competitor_cache;
  int na_cache_limit{128};
  std::unordered_map<std::string, std::vector<int>> shared_trigger_map;
  std::vector<SharedTriggerGroup> shared_trigger_groups;
  std::unordered_map<std::string, int> component_index;
  std::vector<ComponentContextInfo> component_info;
  std::unordered_map<std::string, std::vector<int>> outcome_index;
  std::vector<std::string> outcome_labels;
  std::vector<int> outcome_label_ids;
  std::vector<OutcomeContextInfo> outcome_info;
  std::vector<int> accumulator_label_ids;
  std::vector<int> pool_label_ids;
  ComponentMap components;
  IrContext ir;
  KernelProgram kernel_program;
  KernelStateGraph kernel_state_graph;
  TrialParamsSoA base_params_soa;
  bool has_chained_onsets{false};
};

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep);

} // namespace uuber
