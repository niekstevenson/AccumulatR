#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace uuber {

enum class ParamValueTag : std::uint8_t {
  NumericScalar = 0,
  NumericVector = 1,
  LogicalScalar = 2,
  LogicalVector = 3
};

struct ProtoParamEntry {
  std::string name;
  ParamValueTag tag{ParamValueTag::NumericScalar};
  double numeric_scalar{};
  std::vector<double> numeric_values;
  int logical_scalar{std::numeric_limits<int>::min()};
  std::vector<int> logical_values;
};

struct ProtoAccumulator {
  std::string id;
  std::string dist;
  double onset{};
  int onset_kind{};
  std::string onset_source;
  std::string onset_source_kind;
  double onset_lag{};
  double q{};
  std::vector<std::string> components;
  std::string shared_trigger_id;
  std::vector<ProtoParamEntry> params;
};

struct ProtoPool {
  std::string id;
  int k{};
  std::vector<std::string> members;
};

struct ProtoNode {
  int id{};
  std::uint8_t op{2}; // IrNodeOp::And
  std::vector<int> children;
  std::vector<int> source_label_ids;
  int reference_id{-1};
  int blocker_id{-1};
  int event_acc_idx{-1};
  int event_pool_idx{-1};
  int event_label_id{-1};
  int event_outcome_idx{-1};
  std::uint32_t flags{0};
};

struct ProtoLabelEntry {
  std::string label;
  int id{};
};

struct ProtoOutcomeGuessDonor {
  std::string label;
  double weight{0.0};
  std::string rt_policy;
};

struct ProtoOutcome {
  std::string label;
  int node_id{-1};
  std::vector<int> competitor_ids;
  std::vector<std::string> allowed_components;
  bool maps_to_na{false};
  std::string map_target;
  std::vector<ProtoOutcomeGuessDonor> guess_donors;
};

struct ProtoComponentGuess {
  std::string target;
  std::vector<std::string> weight_labels;
  std::vector<double> weights;
};

struct ProtoComponent {
  std::string id;
  double weight{1.0};
  bool has_guess{false};
  ProtoComponentGuess guess;
};

struct NativePrepProto {
  std::vector<ProtoAccumulator> accumulators;
  std::vector<ProtoPool> pools;
  std::vector<ProtoNode> nodes;
  std::vector<ProtoLabelEntry> label_index;
  std::vector<ProtoOutcome> outcomes;
  std::vector<ProtoComponent> components;
};

std::vector<std::uint8_t> serialize_native_prep(const NativePrepProto& proto);
NativePrepProto deserialize_native_prep(const std::uint8_t* data, std::size_t size);

} // namespace uuber
