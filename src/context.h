#pragma once

#include <Rcpp.h>
#include <cstdint>
#include <deque>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "accumulator.h"

namespace uuber {

struct NativeAccumulator {
  std::string id;
  std::string dist;
  double onset{};
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

struct PoolTemplateEntry {
  int finisher_idx{0};
  std::vector<int> complete_idx;
  std::vector<int> survivor_idx;
  std::vector<int> forced_complete_ids;
  std::vector<int> forced_survive_ids;
};

struct PoolTemplateCacheEntry {
  std::vector<PoolTemplateEntry> templates;
  std::vector<std::vector<int>> finisher_map;
  std::vector<int> shared_index;
};

struct NativeNode {
  int id{};
  std::string kind;
  std::string source;
  LabelRef source_ref;
  std::vector<int> args;
  std::vector<int> source_ids;
  int reference_id{-1};
  int blocker_id{-1};

  int arg_id{-1};
  bool needs_forced{false};
  bool scenario_sensitive{false};
};

struct ComponentGuessPolicy {
  std::string target;
  int target_outcome_idx{-1};
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

struct NativeContext {
  std::vector<NativeAccumulator> accumulators;
  std::vector<NativePool> pools;
  std::vector<NativeNode> nodes;
  std::unordered_map<std::string, int> accumulator_index;
  std::unordered_map<std::string, int> pool_index;
  std::unordered_map<int, int> node_index;
  std::unordered_map<std::string, int> label_to_id;
  mutable std::unordered_map<std::string, PoolTemplateCacheEntry>
      pool_template_cache;
  mutable std::unordered_map<std::string, double> na_map_cache;
  mutable std::unordered_map<std::string, std::deque<std::string>>
      na_cache_order;
  int na_cache_limit{128};
  std::unordered_map<std::string, std::vector<int>> shared_trigger_map;
  std::unordered_map<std::string, int> component_index;
  std::vector<ComponentContextInfo> component_info;
  std::unordered_map<std::string, std::vector<int>> outcome_index;
  std::vector<std::string> outcome_labels;
  std::vector<OutcomeContextInfo> outcome_info;
  ComponentMap components;
};

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep);

} // namespace uuber
