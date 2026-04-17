#include "proto.h"

#include <cstring>
#include <stdexcept>

namespace uuber {
namespace {

constexpr std::uint32_t kPrepMagic = 0x55425052; // 'UBPR'
constexpr std::uint32_t kPrepVersion = 3;

class Serializer {
public:
  void write_u8(std::uint8_t value) { buffer.push_back(value); }

  void write_u32(std::uint32_t value) {
    std::uint8_t bytes[4];
    std::memcpy(bytes, &value, sizeof(value));
    buffer.insert(buffer.end(), bytes, bytes + sizeof(bytes));
  }

  void write_i32(std::int32_t value) {
    std::uint8_t bytes[4];
    std::memcpy(bytes, &value, sizeof(value));
    buffer.insert(buffer.end(), bytes, bytes + sizeof(bytes));
  }

  void write_double(double value) {
    std::uint8_t bytes[sizeof(double)];
    std::memcpy(bytes, &value, sizeof(value));
    buffer.insert(buffer.end(), bytes, bytes + sizeof(bytes));
  }

  void write_bool(bool value) {
    buffer.push_back(static_cast<std::uint8_t>(value ? 1 : 0));
  }

  void write_string(const std::string &str) {
    write_u32(static_cast<std::uint32_t>(str.size()));
    buffer.insert(buffer.end(), str.begin(), str.end());
  }

  void write_string_vec(const std::vector<std::string> &vec) {
    write_u32(static_cast<std::uint32_t>(vec.size()));
    for (const auto &entry : vec) {
      write_string(entry);
    }
  }

  void write_int_vec(const std::vector<int> &vec) {
    write_u32(static_cast<std::uint32_t>(vec.size()));
    for (int value : vec) {
      write_i32(static_cast<std::int32_t>(value));
    }
  }

  void write_double_vec(const std::vector<double> &vec) {
    write_u32(static_cast<std::uint32_t>(vec.size()));
    for (double value : vec) {
      write_double(value);
    }
  }

  std::vector<std::uint8_t> buffer;
};

class Deserializer {
public:
  Deserializer(const std::uint8_t *data, std::size_t size)
      : ptr(data), end(data + size) {}

  std::uint8_t read_u8() {
    ensure(1);
    return *ptr++;
  }

  std::uint32_t read_u32() {
    ensure(sizeof(std::uint32_t));
    std::uint32_t value;
    std::memcpy(&value, ptr, sizeof(value));
    ptr += sizeof(value);
    return value;
  }

  std::int32_t read_i32() {
    ensure(sizeof(std::int32_t));
    std::int32_t value;
    std::memcpy(&value, ptr, sizeof(value));
    ptr += sizeof(value);
    return value;
  }

  double read_double() {
    ensure(sizeof(double));
    double value;
    std::memcpy(&value, ptr, sizeof(value));
    ptr += sizeof(value);
    return value;
  }

  bool read_bool() {
    ensure(1);
    bool value = (*ptr != 0);
    ++ptr;
    return value;
  }

  std::string read_string() {
    std::uint32_t size = read_u32();
    ensure(size);
    std::string out(reinterpret_cast<const char *>(ptr), size);
    ptr += size;
    return out;
  }

  std::vector<std::string> read_string_vec() {
    std::uint32_t size = read_u32();
    std::vector<std::string> out;
    out.reserve(size);
    for (std::uint32_t i = 0; i < size; ++i) {
      out.push_back(read_string());
    }
    return out;
  }

  std::vector<int> read_int_vec() {
    std::uint32_t size = read_u32();
    std::vector<int> out;
    out.reserve(size);
    for (std::uint32_t i = 0; i < size; ++i) {
      out.push_back(static_cast<int>(read_i32()));
    }
    return out;
  }

  std::vector<double> read_double_vec() {
    std::uint32_t size = read_u32();
    std::vector<double> out;
    out.reserve(size);
    for (std::uint32_t i = 0; i < size; ++i) {
      out.push_back(read_double());
    }
    return out;
  }

private:
  void ensure(std::size_t bytes) {
    if (static_cast<std::size_t>(end - ptr) < bytes) {
      throw std::runtime_error(
          "native prep deserialization: unexpected end of buffer");
    }
  }

  const std::uint8_t *ptr;
  const std::uint8_t *end;
};

} // namespace

std::vector<std::uint8_t> serialize_native_prep(const NativePrepProto &proto) {
  Serializer writer;
  writer.write_u32(kPrepMagic);
  writer.write_u32(kPrepVersion);

  writer.write_u32(static_cast<std::uint32_t>(proto.accumulators.size()));
  for (const auto &acc : proto.accumulators) {
    writer.write_string(acc.id);
    writer.write_string(acc.dist);
    writer.write_double(acc.onset);
    writer.write_i32(acc.onset_kind);
    writer.write_string(acc.onset_source);
    writer.write_string(acc.onset_source_kind);
    writer.write_double(acc.onset_lag);
    writer.write_double(acc.q);
    writer.write_string(acc.shared_trigger_id);
    writer.write_string_vec(acc.components);
    writer.write_u32(static_cast<std::uint32_t>(acc.params.size()));
    for (const auto &param : acc.params) {
      writer.write_string(param.name);
      writer.write_u8(static_cast<std::uint8_t>(param.tag));
      switch (param.tag) {
      case ParamValueTag::NumericScalar:
        writer.write_double(param.numeric_scalar);
        break;
      case ParamValueTag::NumericVector:
        writer.write_u32(
            static_cast<std::uint32_t>(param.numeric_values.size()));
        for (double val : param.numeric_values) {
          writer.write_double(val);
        }
        break;
      case ParamValueTag::LogicalScalar:
        writer.write_i32(param.logical_scalar);
        break;
      case ParamValueTag::LogicalVector:
        writer.write_u32(
            static_cast<std::uint32_t>(param.logical_values.size()));
        for (int val : param.logical_values) {
          writer.write_i32(val);
        }
        break;
      default:
        throw std::runtime_error(
            "native prep serialization: unknown param tag");
      }
    }
  }

  writer.write_u32(static_cast<std::uint32_t>(proto.pools.size()));
  for (const auto &pool : proto.pools) {
    writer.write_string(pool.id);
    writer.write_i32(pool.k);
    writer.write_string_vec(pool.members);
  }

  writer.write_u32(static_cast<std::uint32_t>(proto.nodes.size()));
  for (const auto &node : proto.nodes) {
    writer.write_i32(node.id);
    writer.write_u8(node.op);
    writer.write_int_vec(node.children);
    writer.write_int_vec(node.source_label_ids);
    writer.write_i32(node.reference_id);
    writer.write_i32(node.blocker_id);
    writer.write_i32(node.event_acc_idx);
    writer.write_i32(node.event_pool_idx);
    writer.write_i32(node.event_label_id);
    writer.write_i32(node.event_outcome_idx);
    writer.write_u32(node.flags);
  }

  writer.write_u32(static_cast<std::uint32_t>(proto.label_index.size()));
  for (const auto &label : proto.label_index) {
    writer.write_string(label.label);
    writer.write_i32(label.id);
  }

  writer.write_u32(static_cast<std::uint32_t>(proto.outcomes.size()));
  for (const auto &outcome : proto.outcomes) {
    writer.write_string(outcome.label);
    writer.write_i32(outcome.node_id);
    writer.write_int_vec(outcome.competitor_ids);
    writer.write_string_vec(outcome.allowed_components);
    writer.write_bool(outcome.maps_to_na);
    writer.write_string(outcome.map_target);
    writer.write_u32(static_cast<std::uint32_t>(outcome.guess_donors.size()));
    for (const auto &donor : outcome.guess_donors) {
      writer.write_string(donor.label);
      writer.write_double(donor.weight);
      writer.write_string(donor.rt_policy);
    }
  }

  writer.write_u32(static_cast<std::uint32_t>(proto.components.size()));
  for (const auto &component : proto.components) {
    writer.write_string(component.id);
    writer.write_double(component.weight);
    writer.write_bool(component.has_guess);
    if (component.has_guess) {
      writer.write_string(component.guess.target);
      writer.write_string_vec(component.guess.weight_labels);
      writer.write_double_vec(component.guess.weights);
    }
  }

  return std::move(writer.buffer);
}

NativePrepProto deserialize_native_prep(const std::uint8_t *data,
                                        std::size_t size) {
  Deserializer reader(data, size);
  std::uint32_t magic = reader.read_u32();
  if (magic != kPrepMagic) {
    throw std::runtime_error("native prep deserialization: invalid header");
  }
  std::uint32_t version = reader.read_u32();
  if (version != kPrepVersion) {
    throw std::runtime_error(
        "native prep deserialization: unsupported version");
  }

  NativePrepProto proto;

  std::uint32_t acc_count = reader.read_u32();
  proto.accumulators.reserve(acc_count);
  for (std::uint32_t i = 0; i < acc_count; ++i) {
    ProtoAccumulator acc;
    acc.id = reader.read_string();
    acc.dist = reader.read_string();
    acc.onset = reader.read_double();
    acc.onset_kind = reader.read_i32();
    acc.onset_source = reader.read_string();
    acc.onset_source_kind = reader.read_string();
    acc.onset_lag = reader.read_double();
    acc.q = reader.read_double();
    acc.shared_trigger_id = reader.read_string();
    acc.components = reader.read_string_vec();
    std::uint32_t param_count = reader.read_u32();
    acc.params.reserve(param_count);
    for (std::uint32_t j = 0; j < param_count; ++j) {
      ProtoParamEntry param;
      param.name = reader.read_string();
      param.tag = static_cast<ParamValueTag>(reader.read_u8());
      switch (param.tag) {
      case ParamValueTag::NumericScalar:
        param.numeric_scalar = reader.read_double();
        break;
      case ParamValueTag::NumericVector: {
        std::uint32_t n = reader.read_u32();
        param.numeric_values.reserve(n);
        for (std::uint32_t k = 0; k < n; ++k) {
          param.numeric_values.push_back(reader.read_double());
        }
        break;
      }
      case ParamValueTag::LogicalScalar:
        param.logical_scalar = reader.read_i32();
        break;
      case ParamValueTag::LogicalVector: {
        std::uint32_t n = reader.read_u32();
        param.logical_values.reserve(n);
        for (std::uint32_t k = 0; k < n; ++k) {
          param.logical_values.push_back(reader.read_i32());
        }
        break;
      }
      default:
        throw std::runtime_error(
            "native prep deserialization: unknown param tag");
      }
      acc.params.push_back(std::move(param));
    }
    proto.accumulators.push_back(std::move(acc));
  }

  std::uint32_t pool_count = reader.read_u32();
  proto.pools.reserve(pool_count);
  for (std::uint32_t i = 0; i < pool_count; ++i) {
    ProtoPool pool;
    pool.id = reader.read_string();
    pool.k = reader.read_i32();
    pool.members = reader.read_string_vec();
    proto.pools.push_back(std::move(pool));
  }

  std::uint32_t node_count = reader.read_u32();
  proto.nodes.reserve(node_count);
  for (std::uint32_t i = 0; i < node_count; ++i) {
    ProtoNode node;
    node.id = reader.read_i32();
    node.op = reader.read_u8();
    node.children = reader.read_int_vec();
    node.source_label_ids = reader.read_int_vec();
    node.reference_id = reader.read_i32();
    node.blocker_id = reader.read_i32();
    node.event_acc_idx = reader.read_i32();
    node.event_pool_idx = reader.read_i32();
    node.event_label_id = reader.read_i32();
    node.event_outcome_idx = reader.read_i32();
    node.flags = reader.read_u32();
    proto.nodes.push_back(std::move(node));
  }

  std::uint32_t label_count = reader.read_u32();
  proto.label_index.reserve(label_count);
  for (std::uint32_t i = 0; i < label_count; ++i) {
    ProtoLabelEntry label;
    label.label = reader.read_string();
    label.id = reader.read_i32();
    proto.label_index.push_back(std::move(label));
  }

  std::uint32_t outcome_count = reader.read_u32();
  proto.outcomes.reserve(outcome_count);
  for (std::uint32_t i = 0; i < outcome_count; ++i) {
    ProtoOutcome outcome;
    outcome.label = reader.read_string();
    outcome.node_id = reader.read_i32();
    outcome.competitor_ids = reader.read_int_vec();
    outcome.allowed_components = reader.read_string_vec();
    outcome.maps_to_na = reader.read_bool();
    outcome.map_target = reader.read_string();
    std::uint32_t donor_count = reader.read_u32();
    outcome.guess_donors.reserve(donor_count);
    for (std::uint32_t j = 0; j < donor_count; ++j) {
      ProtoOutcomeGuessDonor donor;
      donor.label = reader.read_string();
      donor.weight = reader.read_double();
      donor.rt_policy = reader.read_string();
      outcome.guess_donors.push_back(std::move(donor));
    }
    proto.outcomes.push_back(std::move(outcome));
  }

  std::uint32_t component_count = reader.read_u32();
  proto.components.reserve(component_count);
  for (std::uint32_t i = 0; i < component_count; ++i) {
    ProtoComponent component;
    component.id = reader.read_string();
    component.weight = reader.read_double();
    component.has_guess = reader.read_bool();
    if (component.has_guess) {
      component.guess.target = reader.read_string();
      component.guess.weight_labels = reader.read_string_vec();
      component.guess.weights = reader.read_double_vec();
    }
    proto.components.push_back(std::move(component));
  }

  return proto;
}

} // namespace uuber
