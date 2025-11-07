#pragma once

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>

namespace uuber {

struct NativeAccumulator {
  std::string id;
  std::string dist;
  double onset{};
  double q{};
  Rcpp::List params;
  std::vector<std::string> components;
  std::string shared_trigger_id;
};

struct NativePool {
  std::string id;
  int k{};
  std::vector<std::string> members;
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

struct NativeContext {
  std::vector<NativeAccumulator> accumulators;
  std::vector<NativePool> pools;
  std::vector<NativeNode> nodes;
  std::unordered_map<std::string, int> accumulator_index;
  std::unordered_map<std::string, int> pool_index;
  std::unordered_map<int, int> node_index;
  std::unordered_map<std::string, int> label_to_id;
};

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep);

} // namespace uuber
