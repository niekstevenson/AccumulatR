#pragma once

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include <cstdint>
#include <deque>
#include <deque>

#include "native_accumulator.hpp"

namespace uuber {

struct NativeAccumulator {
  std::string id;
  std::string dist;
  double onset{};
  double q{};
  std::vector<std::string> components;
  std::string shared_trigger_id;
  AccDistParams dist_cfg;
};

struct NativePool {
  std::string id;
  int k{};
  std::vector<std::string> members;
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
  std::vector<int> args;
  std::vector<int> source_ids;
  int reference_id{-1};
  int blocker_id{-1};
  std::vector<int> unless_ids;
  int arg_id{-1};
  bool needs_forced{false};
  bool scenario_sensitive{false};
};

struct ComponentGuessPolicy {
  std::string target;
  std::unordered_map<std::string, double> keep_weights;
  bool valid{false};
};

struct ComponentContextInfo {
  double deadline{std::numeric_limits<double>::quiet_NaN()};
  ComponentGuessPolicy guess;
};

struct OutcomeGuessDonor {
  std::string label;
  double weight{0.0};
  std::string rt_policy;
};

struct OutcomeContextInfo {
  int node_id{-1};
  std::vector<int> competitor_ids;
  std::vector<OutcomeGuessDonor> guess_donors;
  std::vector<std::string> map_sources;
  bool maps_to_na{false};
};

struct CacheMetrics {
  std::uint64_t na_hits{0};
  std::uint64_t na_misses{0};
};

struct NativeContext {
  std::vector<NativeAccumulator> accumulators;
  std::vector<NativePool> pools;
  std::vector<NativeNode> nodes;
  std::unordered_map<std::string, int> accumulator_index;
  std::unordered_map<std::string, int> pool_index;
  std::unordered_map<int, int> node_index;
  std::unordered_map<std::string, int> label_to_id;
  mutable std::unordered_map<std::string, PoolTemplateCacheEntry> pool_template_cache;
  mutable std::unordered_map<std::string, double> na_map_cache;
  mutable std::unordered_map<std::string, std::deque<std::string>> na_cache_order;
  int na_cache_limit{128};
  mutable CacheMetrics cache_metrics;
  mutable std::uint64_t context_builds{0};
  mutable std::uint64_t context_reuses{0};
  std::unordered_map<std::string, std::vector<int>> shared_trigger_map;
  std::unordered_map<std::string, ComponentContextInfo> component_info;
  std::unordered_map<std::string, OutcomeContextInfo> outcome_info;
  std::unordered_map<std::string, std::vector<std::string>> alias_sources;
  std::unordered_map<std::string, std::vector<OutcomeGuessDonor>> guess_target_map;
};

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep);

} // namespace uuber
