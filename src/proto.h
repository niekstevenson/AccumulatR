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
  std::string kind;
  std::string source;
  std::vector<int> args;
  std::vector<int> source_ids;
  int reference_id{-1};
  int blocker_id{-1};

  int arg_id{-1};
  bool needs_forced{false};
  bool scenario_sensitive{false};
};

struct ProtoLabelEntry {
  std::string label;
  int id{};
};

struct NativePrepProto {
  std::vector<ProtoAccumulator> accumulators;
  std::vector<ProtoPool> pools;
  std::vector<ProtoNode> nodes;
  std::vector<ProtoLabelEntry> label_index;
};

std::vector<std::uint8_t> serialize_native_prep(const NativePrepProto& proto);
NativePrepProto deserialize_native_prep(const std::uint8_t* data, std::size_t size);

} // namespace uuber
