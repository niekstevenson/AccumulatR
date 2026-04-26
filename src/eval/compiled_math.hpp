#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "../semantic/model.hpp"

namespace accumulatr::eval {
namespace detail {

enum class CompiledMathValueKind : std::uint8_t {
  Scalar = 0,
  Pdf = 1,
  Cdf = 2,
  Survival = 3,
  Density = 4
};

enum class CompiledMathNodeKind : std::uint8_t {
  Constant = 0,
  SourcePdf = 1,
  SourceCdf = 2,
  SourceSurvival = 3,
  ExprDensity = 4,
  ExprCdf = 5,
  ExprSurvival = 6,
  Product = 7,
  Sum = 8,
  CleanSignedSum = 9,
  ClampProbability = 10,
  Complement = 11,
  Negate = 12,
  SimpleGuardDensity = 13,
  SimpleGuardCdf = 14,
  IntegralZeroToCurrent = 17,
  ExprUpperBoundDensity = 18,
  ExprUpperBoundCdf = 19,
  UnionKernelDensity = 20,
  UnionKernelCdf = 21,
  UnionKernelMultiSubsetDensity = 22,
  UnionKernelMultiSubsetCdf = 23,
  OutcomeSubsetUnused = 24,
  RootValue = 25,
  IntegralZeroToCurrentRaw = 26,
  LazyProduct = 27,
  SourceChannelLoad = 28,
  TimeGate = 29
};

enum class CompiledMathIntegralKernelKind : std::uint8_t {
  Schedule = 0,
  SourceProduct = 1,
  SourceProductSum = 2
};

enum class CompiledSourceChannelKernelKind : std::uint8_t {
  LeafAbsolute = 0,
  LeafOnsetConvolution = 1,
  PoolKOfN = 2,
  Invalid = 255
};

struct CompiledMathIndexSpan {
  semantic::Index offset{0};
  semantic::Index size{0};

  [[nodiscard]] bool empty() const noexcept {
    return size == 0;
  }
};

struct CompiledMathPairFactEntry {
  semantic::Index first{semantic::kInvalidIndex};
  semantic::Index second{semantic::kInvalidIndex};
  CompiledMathIndexSpan facts{};
};

struct CompiledMathPairFactLookup {
  std::vector<CompiledMathPairFactEntry> entries;
  std::vector<semantic::Index> fact_indices;
};

struct CompiledMathNode {
  CompiledMathNodeKind kind{CompiledMathNodeKind::Constant};
  CompiledMathValueKind value_kind{CompiledMathValueKind::Scalar};
  semantic::Index subject_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index time_id{0};
  semantic::Index aux_id{semantic::kInvalidIndex};
  semantic::Index aux2_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  CompiledMathIndexSpan children{};
  semantic::Index cache_slot{semantic::kInvalidIndex};
  semantic::Index source_channel_slot{semantic::kInvalidIndex};
  semantic::Index integral_kernel_slot{semantic::kInvalidIndex};
  double constant{0.0};
};

struct CompiledMathRoot {
  semantic::Index node_id{semantic::kInvalidIndex};
  CompiledMathIndexSpan schedule{};
};

struct CompiledMathIntegralSourceProductTerm {
  CompiledMathIndexSpan source_value_nodes{};
  double sign{1.0};
};

struct CompiledMathIntegralKernel {
  CompiledMathIntegralKernelKind kind{CompiledMathIntegralKernelKind::Schedule};
  semantic::Index root_id{semantic::kInvalidIndex};
  semantic::Index result_node_id{semantic::kInvalidIndex};
  CompiledMathIndexSpan schedule{};
  semantic::Index bind_time_id{semantic::kInvalidIndex};
  CompiledMathIndexSpan source_value_nodes{};
  CompiledMathIndexSpan source_product_terms{};
  bool clean_signed_source_sum{false};
};

struct CompiledMathNodeKey {
  CompiledMathNodeKind kind{CompiledMathNodeKind::Constant};
  CompiledMathValueKind value_kind{CompiledMathValueKind::Scalar};
  semantic::Index subject_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index time_id{0};
  semantic::Index aux_id{semantic::kInvalidIndex};
  semantic::Index aux2_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  double constant{0.0};
  std::vector<semantic::Index> children;

  bool operator==(const CompiledMathNodeKey &other) const noexcept {
    return kind == other.kind &&
           value_kind == other.value_kind &&
           subject_id == other.subject_id &&
           condition_id == other.condition_id &&
           time_id == other.time_id &&
           aux_id == other.aux_id &&
           aux2_id == other.aux2_id &&
           source_view_id == other.source_view_id &&
           constant == other.constant &&
           children == other.children;
  }
};

struct CompiledMathNodeKeyHash {
  std::size_t operator()(const CompiledMathNodeKey &key) const noexcept {
    std::size_t seed = static_cast<std::size_t>(key.kind);
    hash_combine(&seed, static_cast<std::size_t>(key.value_kind));
    hash_combine(&seed, static_cast<std::size_t>(key.subject_id));
    hash_combine(&seed, static_cast<std::size_t>(key.condition_id));
    hash_combine(&seed, static_cast<std::size_t>(key.time_id));
    hash_combine(&seed, static_cast<std::size_t>(key.aux_id));
    hash_combine(&seed, static_cast<std::size_t>(key.aux2_id));
    hash_combine(&seed, static_cast<std::size_t>(key.source_view_id));
    hash_combine(&seed, std::hash<double>{}(key.constant));
    for (const auto child : key.children) {
      hash_combine(&seed, static_cast<std::size_t>(child));
    }
    return seed;
  }

private:
  static void hash_combine(std::size_t *seed, const std::size_t value) noexcept {
    *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
  }
};

struct CompiledMathSourceChannelKey {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};

  bool operator==(const CompiledMathSourceChannelKey &other) const noexcept {
    return source_id == other.source_id &&
           condition_id == other.condition_id &&
           time_id == other.time_id &&
           time_cap_id == other.time_cap_id &&
           source_view_id == other.source_view_id;
  }
};

struct CompiledSourceChannelPlan {
  CompiledMathSourceChannelKey request{};
  CompiledSourceChannelKernelKind kernel{
      CompiledSourceChannelKernelKind::Invalid};
  bool has_source_condition_overlay{false};
};

struct CompiledMathSourceChannelKeyHash {
  std::size_t operator()(const CompiledMathSourceChannelKey &key) const
      noexcept {
    std::size_t seed = static_cast<std::size_t>(key.source_id);
    hash_combine(&seed, static_cast<std::size_t>(key.condition_id));
    hash_combine(&seed, static_cast<std::size_t>(key.time_id));
    hash_combine(&seed, static_cast<std::size_t>(key.time_cap_id));
    hash_combine(&seed, static_cast<std::size_t>(key.source_view_id));
    return seed;
  }

private:
  static void hash_combine(std::size_t *seed, const std::size_t value) noexcept {
    *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
  }
};

struct CompiledMathConditionKey {
  bool impossible{false};
  std::vector<semantic::Index> source_ids;
  std::vector<std::uint8_t> relations;
  std::vector<std::uint8_t> fact_kinds;
  std::vector<semantic::Index> fact_subject_ids;
  std::vector<semantic::Index> fact_aux_ids;
  std::vector<semantic::Index> fact_aux2_ids;
  std::vector<semantic::Index> fact_time_ids;
  std::vector<semantic::Index> fact_normalizer_node_ids;
  std::vector<semantic::Index> fact_normalizer_root_ids;

  bool operator==(const CompiledMathConditionKey &other) const noexcept {
    return impossible == other.impossible &&
           source_ids == other.source_ids &&
           relations == other.relations &&
           fact_kinds == other.fact_kinds &&
           fact_subject_ids == other.fact_subject_ids &&
           fact_aux_ids == other.fact_aux_ids &&
           fact_aux2_ids == other.fact_aux2_ids &&
           fact_time_ids == other.fact_time_ids &&
           fact_normalizer_node_ids == other.fact_normalizer_node_ids &&
           fact_normalizer_root_ids == other.fact_normalizer_root_ids;
  }
};

struct CompiledMathConditionKeyHash {
  std::size_t operator()(const CompiledMathConditionKey &key) const noexcept {
    std::size_t seed = 0;
    hash_combine(&seed, key.impossible ? 1U : 0U);
    for (const auto source_id : key.source_ids) {
      hash_combine(&seed, static_cast<std::size_t>(source_id));
    }
    for (const auto relation : key.relations) {
      hash_combine(&seed, static_cast<std::size_t>(relation));
    }
    for (const auto kind : key.fact_kinds) {
      hash_combine(&seed, static_cast<std::size_t>(kind));
    }
    for (const auto subject_id : key.fact_subject_ids) {
      hash_combine(&seed, static_cast<std::size_t>(subject_id));
    }
    for (const auto aux_id : key.fact_aux_ids) {
      hash_combine(&seed, static_cast<std::size_t>(aux_id));
    }
    for (const auto aux2_id : key.fact_aux2_ids) {
      hash_combine(&seed, static_cast<std::size_t>(aux2_id));
    }
    for (const auto time_id : key.fact_time_ids) {
      hash_combine(&seed, static_cast<std::size_t>(time_id));
    }
    for (const auto node_id : key.fact_normalizer_node_ids) {
      hash_combine(&seed, static_cast<std::size_t>(node_id));
    }
    for (const auto root_id : key.fact_normalizer_root_ids) {
      hash_combine(&seed, static_cast<std::size_t>(root_id));
    }
    return seed;
  }

private:
  static void hash_combine(std::size_t *seed, const std::size_t value) noexcept {
    *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
  }
};

enum class CompiledMathConditionFactKind : std::uint8_t {
  SourceExact = 0,
  SourceUpperBound = 1,
  SourceLowerBound = 2,
  ExprUpperBound = 3,
  SourceOrder = 4,
  GuardUpperBound = 5
};

enum class CompiledMathTimeSlot : semantic::Index {
  Observed = 0,
  Readiness = 1,
  Active = 2,
  Zero = 3
};

struct CompiledMathCondition {
  bool impossible{false};
  std::vector<semantic::Index> source_ids;
  std::vector<std::uint8_t> relations;
  std::vector<std::uint8_t> fact_kinds;
  std::vector<semantic::Index> fact_subject_ids;
  std::vector<semantic::Index> fact_aux_ids;
  std::vector<semantic::Index> fact_aux2_ids;
  std::vector<semantic::Index> fact_time_ids;
  std::vector<semantic::Index> fact_normalizer_node_ids;
  std::vector<semantic::Index> fact_normalizer_root_ids;
  std::vector<CompiledMathIndexSpan> source_exact_fact_spans;
  std::vector<semantic::Index> source_exact_fact_indices;
  std::vector<CompiledMathIndexSpan> source_lower_fact_spans;
  std::vector<semantic::Index> source_lower_fact_indices;
  std::vector<CompiledMathIndexSpan> source_upper_fact_spans;
  std::vector<semantic::Index> source_upper_fact_indices;
  std::vector<CompiledMathIndexSpan> expr_upper_fact_spans;
  std::vector<semantic::Index> expr_upper_fact_indices;
  CompiledMathPairFactLookup source_order_fact_lookup;
  CompiledMathPairFactLookup guard_upper_fact_lookup;
};

struct CompiledMathProgram {
  std::vector<CompiledMathNode> nodes;
  std::vector<semantic::Index> child_nodes;
  std::vector<CompiledMathRoot> roots;
  std::vector<semantic::Index> root_schedule_nodes;
  std::vector<CompiledMathCondition> conditions;
  std::vector<CompiledMathIntegralKernel> integral_kernels;
  std::vector<CompiledMathIntegralSourceProductTerm>
      integral_kernel_source_product_terms;
  std::vector<semantic::Index> integral_kernel_source_value_nodes;
  std::vector<CompiledMathSourceChannelKey> source_channel_keys;
  std::vector<CompiledSourceChannelPlan> source_channel_plans;
  std::vector<semantic::Index> source_channel_use_counts;
  semantic::Index source_channel_count{0};
  std::unordered_map<
      CompiledMathNodeKey,
      semantic::Index,
      CompiledMathNodeKeyHash>
      node_index;
  std::unordered_map<
      CompiledMathSourceChannelKey,
      semantic::Index,
      CompiledMathSourceChannelKeyHash>
      source_channel_index;
  std::unordered_map<
      CompiledMathConditionKey,
      semantic::Index,
      CompiledMathConditionKeyHash>
      condition_index;
};

struct CompiledMathWorkspace {
  CompiledMathWorkspace() = default;

  explicit CompiledMathWorkspace(const CompiledMathProgram &program) {
    resize(program);
  }

  void resize(const CompiledMathProgram &program) {
    values.assign(program.nodes.size(), 0.0);
    cache_valid.assign(program.nodes.size(), 0U);
    cache_condition_ids.assign(program.nodes.size(), 0);
    cache_times.assign(program.nodes.size(), 0.0);
    cache_evaluators.assign(program.nodes.size(), nullptr);
    cache_values.assign(program.nodes.size(), 0.0);
    source_channel_valid.assign(program.source_channel_count, 0U);
    source_channel_condition_ids.assign(program.source_channel_count, 0);
    source_channel_times.assign(program.source_channel_count, 0.0);
    source_channel_evaluators.assign(program.source_channel_count, nullptr);
    source_channel_pdf.assign(program.source_channel_count, 0.0);
    source_channel_cdf.assign(program.source_channel_count, 0.0);
    source_channel_survival.assign(program.source_channel_count, 1.0);
    time_values.assign(4U, 0.0);
    time_valid.assign(4U, 0U);
  }

  void ensure_size(const CompiledMathProgram &program) {
    if (values.size() < program.nodes.size()) {
      const auto size = program.nodes.size();
      values.resize(size, 0.0);
      cache_valid.resize(size, 0U);
      cache_condition_ids.resize(size, 0);
      cache_times.resize(size, 0.0);
      cache_evaluators.resize(size, nullptr);
      cache_values.resize(size, 0.0);
    }
    if (source_channel_valid.size() <
        static_cast<std::size_t>(program.source_channel_count)) {
      const auto size =
          static_cast<std::size_t>(program.source_channel_count);
      source_channel_valid.resize(size, 0U);
      source_channel_condition_ids.resize(size, 0);
      source_channel_times.resize(size, 0.0);
      source_channel_evaluators.resize(size, nullptr);
      source_channel_pdf.resize(size, 0.0);
      source_channel_cdf.resize(size, 0.0);
      source_channel_survival.resize(size, 1.0);
    }
  }

  void reset_cache() {
    std::fill(cache_valid.begin(), cache_valid.end(), 0U);
    std::fill(source_channel_valid.begin(), source_channel_valid.end(), 0U);
  }

  void set_time(const semantic::Index time_id, const double value) {
    const auto pos = static_cast<std::size_t>(time_id);
    if (time_values.size() <= pos) {
      time_values.resize(pos + 1U, 0.0);
      time_valid.resize(pos + 1U, 0U);
    }
    time_values[pos] = value;
    time_valid[pos] = 1U;
  }

  bool has_time(const semantic::Index time_id) const {
    const auto pos = static_cast<std::size_t>(time_id);
    return pos < time_valid.size() && time_valid[pos] != 0U;
  }

  double time(const semantic::Index time_id) const {
    const auto pos = static_cast<std::size_t>(time_id);
    if (pos >= time_valid.size() || time_valid[pos] == 0U) {
      throw std::runtime_error("compiled math time slot is unbound");
    }
    return time_values[pos];
  }

  void copy_times_from(const CompiledMathWorkspace &other) {
    time_values = other.time_values;
    time_valid = other.time_valid;
    used_outcomes = other.used_outcomes;
  }

  class TimeBinding {
  public:
    TimeBinding(CompiledMathWorkspace *workspace,
                const semantic::Index time_id,
                const double value)
        : workspace_(workspace), time_id_(time_id) {
      if (workspace_ == nullptr) {
        return;
      }
      const auto pos = static_cast<std::size_t>(time_id_);
      if (workspace_->time_values.size() <= pos) {
        workspace_->time_values.resize(pos + 1U, 0.0);
        workspace_->time_valid.resize(pos + 1U, 0U);
      }
      had_previous_ = workspace_->time_valid[pos] != 0U;
      previous_ = workspace_->time_values[pos];
      workspace_->time_values[pos] = value;
      workspace_->time_valid[pos] = 1U;
    }

    TimeBinding(const TimeBinding &) = delete;
    TimeBinding &operator=(const TimeBinding &) = delete;

    ~TimeBinding() {
      if (workspace_ == nullptr) {
        return;
      }
      const auto pos = static_cast<std::size_t>(time_id_);
      if (had_previous_) {
        workspace_->time_values[pos] = previous_;
        workspace_->time_valid[pos] = 1U;
      } else {
        workspace_->time_values[pos] = 0.0;
        workspace_->time_valid[pos] = 0U;
      }
    }

  private:
    CompiledMathWorkspace *workspace_{nullptr};
    semantic::Index time_id_{0};
    double previous_{0.0};
    bool had_previous_{false};
  };

  CompiledMathWorkspace &integral_workspace_for(
      const CompiledMathProgram &program) {
    if (!integral_workspace) {
      integral_workspace = std::make_unique<CompiledMathWorkspace>(program);
    } else {
      integral_workspace->ensure_size(program);
      integral_workspace->reset_cache();
    }
    integral_workspace->copy_times_from(*this);
    return *integral_workspace;
  }

  std::vector<double> values;
  std::vector<std::uint8_t> cache_valid;
  std::vector<semantic::Index> cache_condition_ids;
  std::vector<double> cache_times;
  std::vector<const void *> cache_evaluators;
  std::vector<double> cache_values;
  std::vector<std::uint8_t> source_channel_valid;
  std::vector<semantic::Index> source_channel_condition_ids;
  std::vector<double> source_channel_times;
  std::vector<const void *> source_channel_evaluators;
  std::vector<double> source_channel_pdf;
  std::vector<double> source_channel_cdf;
  std::vector<double> source_channel_survival;
  std::vector<double> time_values;
  std::vector<std::uint8_t> time_valid;
  const std::vector<std::uint8_t> *used_outcomes{nullptr};
  std::unique_ptr<CompiledMathWorkspace> integral_workspace;
};

inline bool compiled_math_is_source_value_node(
    const CompiledMathNodeKind kind) noexcept {
  return kind == CompiledMathNodeKind::SourcePdf ||
         kind == CompiledMathNodeKind::SourceCdf ||
         kind == CompiledMathNodeKind::SourceSurvival;
}

inline bool compiled_math_is_source_channel_node(
    const CompiledMathNodeKind kind) noexcept {
  return compiled_math_is_source_value_node(kind) ||
         kind == CompiledMathNodeKind::SourceChannelLoad;
}

inline bool compiled_math_is_integral_node(
    const CompiledMathNodeKind kind) noexcept {
  return kind == CompiledMathNodeKind::IntegralZeroToCurrent ||
         kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw ||
         kind == CompiledMathNodeKind::UnionKernelMultiSubsetCdf;
}

inline semantic::Index compiled_math_integral_root_id(
    const CompiledMathNode &node) noexcept {
  if (node.kind == CompiledMathNodeKind::UnionKernelMultiSubsetCdf) {
    return node.aux_id;
  }
  return node.subject_id;
}

inline semantic::Index compiled_math_integral_bind_time_id(
    const CompiledMathNode &node) noexcept {
  return node.aux2_id == semantic::kInvalidIndex ? node.time_id : node.aux2_id;
}

inline semantic::Index compiled_math_source_channel_slot(
    CompiledMathProgram *program,
    const CompiledMathNode &node) {
  if (!compiled_math_is_source_channel_node(node.kind)) {
    return semantic::kInvalidIndex;
  }
  CompiledMathSourceChannelKey key;
  key.source_id = node.subject_id;
  key.condition_id = node.condition_id;
  key.time_id = node.time_id;
  key.time_cap_id = node.aux_id;
  key.source_view_id = node.source_view_id;
  const auto found = program->source_channel_index.find(key);
  if (found != program->source_channel_index.end()) {
    if (compiled_math_is_source_value_node(node.kind)) {
      const auto slot_pos = static_cast<std::size_t>(found->second);
      if (slot_pos >= program->source_channel_use_counts.size()) {
        program->source_channel_use_counts.resize(slot_pos + 1U, 0);
      }
      ++program->source_channel_use_counts[slot_pos];
    }
    return found->second;
  }
  const auto slot = program->source_channel_count++;
  program->source_channel_keys.push_back(key);
  program->source_channel_use_counts.push_back(
      compiled_math_is_source_value_node(node.kind) ? 1 : 0);
  program->source_channel_index.emplace(std::move(key), slot);
  return slot;
}

inline bool compiled_math_source_product_kernel_nodes(
    const CompiledMathProgram &program,
    const semantic::Index node_id,
    std::vector<semantic::Index> *source_value_nodes) {
  if (node_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
  if (compiled_math_is_source_value_node(node.kind)) {
    source_value_nodes->push_back(node_id);
    return true;
  }
  if (node.kind != CompiledMathNodeKind::Product) {
    return false;
  }
  if (node.children.empty()) {
    return false;
  }
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + i)];
    if (!compiled_math_source_product_kernel_nodes(
            program, child_id, source_value_nodes)) {
      return false;
    }
  }
  return !source_value_nodes->empty();
}

inline bool compiled_math_collect_source_product_terms(
    CompiledMathProgram *program,
    const semantic::Index node_id,
    const double sign,
    std::vector<CompiledMathIntegralSourceProductTerm> *terms,
    bool *clean_signed) {
  if (node_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &node = program->nodes[static_cast<std::size_t>(node_id)];
  if (node.kind == CompiledMathNodeKind::Sum ||
      node.kind == CompiledMathNodeKind::CleanSignedSum) {
    if (node.kind == CompiledMathNodeKind::CleanSignedSum &&
        clean_signed != nullptr) {
      *clean_signed = true;
    }
    if (node.children.empty()) {
      return false;
    }
    for (semantic::Index i = 0; i < node.children.size; ++i) {
      const auto child_id = program->child_nodes[
          static_cast<std::size_t>(node.children.offset + i)];
      if (!compiled_math_collect_source_product_terms(
              program, child_id, sign, terms, clean_signed)) {
        return false;
      }
    }
    return true;
  }
  if (node.kind == CompiledMathNodeKind::Negate && node.children.size == 1U) {
    const auto child_id = program->child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    return compiled_math_collect_source_product_terms(
        program, child_id, -sign, terms, clean_signed);
  }

  std::vector<semantic::Index> source_value_nodes;
  if (!compiled_math_source_product_kernel_nodes(
          *program, node_id, &source_value_nodes)) {
    return false;
  }
  const auto offset = static_cast<semantic::Index>(
      program->integral_kernel_source_value_nodes.size());
  program->integral_kernel_source_value_nodes.insert(
      program->integral_kernel_source_value_nodes.end(),
      source_value_nodes.begin(),
      source_value_nodes.end());
  terms->push_back(
      CompiledMathIntegralSourceProductTerm{
          CompiledMathIndexSpan{
              offset,
              static_cast<semantic::Index>(source_value_nodes.size())},
          sign});
  return true;
}

inline semantic::Index compiled_math_integral_kernel_slot(
    CompiledMathProgram *program,
    const CompiledMathNode &node) {
  if (!compiled_math_is_integral_node(node.kind)) {
    return semantic::kInvalidIndex;
  }
  const auto root_id = compiled_math_integral_root_id(node);
  if (root_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto root_pos = static_cast<std::size_t>(root_id);
  if (root_pos >= program->roots.size()) {
    throw std::runtime_error(
        "compiled integral node references an unplanned integrand root");
  }
  const auto &root = program->roots[root_pos];
  const auto slot =
      static_cast<semantic::Index>(program->integral_kernels.size());
  CompiledMathIntegralKernel kernel;
  kernel.root_id = root_id;
  kernel.result_node_id = root.node_id;
  kernel.schedule = root.schedule;
  kernel.bind_time_id = compiled_math_integral_bind_time_id(node);

  std::vector<CompiledMathIntegralSourceProductTerm> source_product_terms;
  bool clean_signed_source_sum = false;
  const auto source_value_node_mark =
      program->integral_kernel_source_value_nodes.size();
  if (compiled_math_collect_source_product_terms(
          program,
          root.node_id,
          1.0,
          &source_product_terms,
          &clean_signed_source_sum)) {
    if (source_product_terms.size() == 1U &&
        source_product_terms.front().sign == 1.0 &&
        !clean_signed_source_sum) {
      kernel.kind = CompiledMathIntegralKernelKind::SourceProduct;
      kernel.source_value_nodes =
          source_product_terms.front().source_value_nodes;
    } else {
      kernel.kind = CompiledMathIntegralKernelKind::SourceProductSum;
      kernel.source_product_terms = CompiledMathIndexSpan{
          static_cast<semantic::Index>(
              program->integral_kernel_source_product_terms.size()),
          static_cast<semantic::Index>(source_product_terms.size())};
      kernel.clean_signed_source_sum = clean_signed_source_sum;
      program->integral_kernel_source_product_terms.insert(
          program->integral_kernel_source_product_terms.end(),
          source_product_terms.begin(),
          source_product_terms.end());
    }
  } else {
    program->integral_kernel_source_value_nodes.resize(
        source_value_node_mark);
  }
  program->integral_kernels.push_back(kernel);
  return slot;
}

inline semantic::Index compiled_math_intern_node(
    CompiledMathProgram *program,
    CompiledMathNodeKey key) {
  const auto found = program->node_index.find(key);
  if (found != program->node_index.end()) {
    return found->second;
  }

  const auto node_id = static_cast<semantic::Index>(program->nodes.size());
  const auto child_offset =
      static_cast<semantic::Index>(program->child_nodes.size());
  program->child_nodes.insert(
      program->child_nodes.end(), key.children.begin(), key.children.end());

  CompiledMathNode node;
  node.kind = key.kind;
  node.value_kind = key.value_kind;
  node.subject_id = key.subject_id;
  node.condition_id = key.condition_id;
  node.time_id = key.time_id;
  node.aux_id = key.aux_id;
  node.aux2_id = key.aux2_id;
  node.source_view_id = key.source_view_id;
  node.children = CompiledMathIndexSpan{
      child_offset,
      static_cast<semantic::Index>(key.children.size())};
  node.cache_slot = node_id;
  node.source_channel_slot = compiled_math_source_channel_slot(program, node);
  node.integral_kernel_slot = compiled_math_integral_kernel_slot(program, node);
  node.constant = key.constant;
  program->nodes.push_back(node);
  program->node_index.emplace(std::move(key), node_id);
  return node_id;
}

inline void compiled_math_append_fact_index(
    std::vector<std::vector<semantic::Index>> *by_key,
    const semantic::Index key,
    const semantic::Index fact_index) {
  if (key == semantic::kInvalidIndex) {
    return;
  }
  const auto pos = static_cast<std::size_t>(key);
  if (by_key->size() <= pos) {
    by_key->resize(pos + 1U);
  }
  (*by_key)[pos].push_back(fact_index);
}

inline void compiled_math_finish_fact_spans(
    const std::vector<std::vector<semantic::Index>> &by_key,
    std::vector<CompiledMathIndexSpan> *spans,
    std::vector<semantic::Index> *indices) {
  spans->assign(by_key.size(), CompiledMathIndexSpan{});
  indices->clear();
  for (std::size_t i = 0; i < by_key.size(); ++i) {
    const auto offset = static_cast<semantic::Index>(indices->size());
    indices->insert(indices->end(), by_key[i].begin(), by_key[i].end());
    (*spans)[i] = CompiledMathIndexSpan{
        offset,
        static_cast<semantic::Index>(by_key[i].size())};
  }
}

inline void compiled_math_append_pair_fact_index(
    std::vector<CompiledMathPairFactEntry> *entries,
    std::vector<semantic::Index> *indices,
    const semantic::Index first,
    const semantic::Index second,
    const semantic::Index fact_index) {
  if (first == semantic::kInvalidIndex ||
      second == semantic::kInvalidIndex) {
    return;
  }
  entries->push_back(CompiledMathPairFactEntry{
      first,
      second,
      CompiledMathIndexSpan{
          static_cast<semantic::Index>(indices->size()),
          1}});
  indices->push_back(fact_index);
}

inline void compiled_math_finish_pair_lookup(
    CompiledMathPairFactLookup *lookup) {
  if (lookup->entries.empty()) {
    return;
  }
  std::vector<CompiledMathPairFactEntry> old_entries =
      std::move(lookup->entries);
  std::vector<semantic::Index> old_indices = std::move(lookup->fact_indices);
  std::vector<std::size_t> order(old_entries.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }
  std::sort(
      order.begin(),
      order.end(),
      [&](const auto lhs, const auto rhs) {
        const auto &a = old_entries[lhs];
        const auto &b = old_entries[rhs];
        return a.first < b.first ||
               (a.first == b.first && a.second < b.second);
      });
  lookup->entries.clear();
  lookup->fact_indices.clear();
  std::size_t order_pos = 0;
  while (order_pos < order.size()) {
    const auto first = old_entries[order[order_pos]].first;
    const auto second = old_entries[order[order_pos]].second;
    const auto offset =
        static_cast<semantic::Index>(lookup->fact_indices.size());
    while (order_pos < order.size() &&
           old_entries[order[order_pos]].first == first &&
           old_entries[order[order_pos]].second == second) {
      const auto old_offset = static_cast<std::size_t>(
          old_entries[order[order_pos]].facts.offset);
      lookup->fact_indices.push_back(old_indices[old_offset]);
      ++order_pos;
    }
    lookup->entries.push_back(CompiledMathPairFactEntry{
        first,
        second,
        CompiledMathIndexSpan{
            offset,
            static_cast<semantic::Index>(
                lookup->fact_indices.size() -
                static_cast<std::size_t>(offset))}});
  }
}

inline void compiled_math_build_condition_access(
    CompiledMathCondition *condition) {
  std::vector<std::vector<semantic::Index>> source_exact;
  std::vector<std::vector<semantic::Index>> source_lower;
  std::vector<std::vector<semantic::Index>> source_upper;
  std::vector<std::vector<semantic::Index>> expr_upper;
  condition->source_order_fact_lookup.entries.clear();
  condition->source_order_fact_lookup.fact_indices.clear();
  condition->guard_upper_fact_lookup.entries.clear();
  condition->guard_upper_fact_lookup.fact_indices.clear();

  for (std::size_t i = 0; i < condition->fact_kinds.size(); ++i) {
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(condition->fact_kinds[i]);
    const auto fact_index = static_cast<semantic::Index>(i);
    switch (kind) {
    case CompiledMathConditionFactKind::SourceExact:
      compiled_math_append_fact_index(
          &source_exact, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::SourceLowerBound:
      compiled_math_append_fact_index(
          &source_lower, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::SourceUpperBound:
      compiled_math_append_fact_index(
          &source_upper, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::ExprUpperBound:
      compiled_math_append_fact_index(
          &expr_upper, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::SourceOrder:
      compiled_math_append_pair_fact_index(
          &condition->source_order_fact_lookup.entries,
          &condition->source_order_fact_lookup.fact_indices,
          condition->fact_subject_ids[i],
          condition->fact_aux_ids[i],
          fact_index);
      break;
    case CompiledMathConditionFactKind::GuardUpperBound:
      compiled_math_append_pair_fact_index(
          &condition->guard_upper_fact_lookup.entries,
          &condition->guard_upper_fact_lookup.fact_indices,
          condition->fact_aux_ids[i],
          condition->fact_aux2_ids[i],
          fact_index);
      break;
    }
  }

  compiled_math_finish_fact_spans(
      source_exact,
      &condition->source_exact_fact_spans,
      &condition->source_exact_fact_indices);
  compiled_math_finish_fact_spans(
      source_lower,
      &condition->source_lower_fact_spans,
      &condition->source_lower_fact_indices);
  compiled_math_finish_fact_spans(
      source_upper,
      &condition->source_upper_fact_spans,
      &condition->source_upper_fact_indices);
  compiled_math_finish_fact_spans(
      expr_upper,
      &condition->expr_upper_fact_spans,
      &condition->expr_upper_fact_indices);
  compiled_math_finish_pair_lookup(&condition->source_order_fact_lookup);
  compiled_math_finish_pair_lookup(&condition->guard_upper_fact_lookup);
}

inline semantic::Index compiled_math_constant(CompiledMathProgram *program,
                                             const double value) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::Constant;
  key.value_kind = CompiledMathValueKind::Scalar;
  key.constant = value;
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_intern_condition(
    CompiledMathProgram *program,
    CompiledMathConditionKey key) {
  if (key.source_ids.empty() && key.fact_kinds.empty()) {
    return 0;
  }
  const auto found = program->condition_index.find(key);
  if (found != program->condition_index.end()) {
    return found->second;
  }
  const auto condition_id =
      static_cast<semantic::Index>(program->conditions.size() + 1U);
  program->conditions.push_back(
      CompiledMathCondition{
          key.impossible,
          key.source_ids,
          key.relations,
          key.fact_kinds,
          key.fact_subject_ids,
          key.fact_aux_ids,
          key.fact_aux2_ids,
          key.fact_time_ids,
          key.fact_normalizer_node_ids,
          key.fact_normalizer_root_ids});
  compiled_math_build_condition_access(&program->conditions.back());
  program->condition_index.emplace(std::move(key), condition_id);
  return condition_id;
}

inline bool compiled_condition_impossible(
    const CompiledMathProgram &program,
    const semantic::Index condition_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto pos = static_cast<std::size_t>(condition_id - 1U);
  return pos >= program.conditions.size() || program.conditions[pos].impossible;
}

inline void compiled_math_append_condition_to_key(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    CompiledMathConditionKey *key) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return;
  }
  const auto pos = static_cast<std::size_t>(condition_id - 1U);
  if (pos >= program.conditions.size()) {
    key->impossible = true;
    return;
  }
  const auto &condition = program.conditions[pos];
  key->impossible = key->impossible || condition.impossible;
  key->source_ids.insert(
      key->source_ids.end(),
      condition.source_ids.begin(),
      condition.source_ids.end());
  key->relations.insert(
      key->relations.end(),
      condition.relations.begin(),
      condition.relations.end());
  key->fact_kinds.insert(
      key->fact_kinds.end(),
      condition.fact_kinds.begin(),
      condition.fact_kinds.end());
  key->fact_subject_ids.insert(
      key->fact_subject_ids.end(),
      condition.fact_subject_ids.begin(),
      condition.fact_subject_ids.end());
  key->fact_aux_ids.insert(
      key->fact_aux_ids.end(),
      condition.fact_aux_ids.begin(),
      condition.fact_aux_ids.end());
  key->fact_aux2_ids.insert(
      key->fact_aux2_ids.end(),
      condition.fact_aux2_ids.begin(),
      condition.fact_aux2_ids.end());
  key->fact_time_ids.insert(
      key->fact_time_ids.end(),
      condition.fact_time_ids.begin(),
      condition.fact_time_ids.end());
  key->fact_normalizer_node_ids.insert(
      key->fact_normalizer_node_ids.end(),
      condition.fact_normalizer_node_ids.begin(),
      condition.fact_normalizer_node_ids.end());
  key->fact_normalizer_root_ids.insert(
      key->fact_normalizer_root_ids.end(),
      condition.fact_normalizer_root_ids.begin(),
      condition.fact_normalizer_root_ids.end());
}

inline semantic::Index compiled_math_merge_conditions(
    CompiledMathProgram *program,
    const std::initializer_list<semantic::Index> condition_ids) {
  CompiledMathConditionKey key;
  for (const auto condition_id : condition_ids) {
    compiled_math_append_condition_to_key(*program, condition_id, &key);
  }
  return compiled_math_intern_condition(program, std::move(key));
}

inline semantic::Index compiled_math_source_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    const semantic::Index source_id,
    const semantic::Index condition_id = 0,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index time_cap_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey load_key;
  load_key.kind = CompiledMathNodeKind::SourceChannelLoad;
  load_key.subject_id = source_id;
  load_key.condition_id = condition_id;
  load_key.time_id = time_id;
  load_key.aux_id = time_cap_id;
  load_key.source_view_id = source_view_id;
  const auto load_node_id =
      compiled_math_intern_node(program, std::move(load_key));

  CompiledMathNodeKey key;
  key.kind = kind;
  key.subject_id = source_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.aux_id = time_cap_id;
  key.source_view_id = source_view_id;
  key.children.push_back(load_node_id);
  if (kind == CompiledMathNodeKind::SourcePdf) {
    key.value_kind = CompiledMathValueKind::Pdf;
  } else if (kind == CompiledMathNodeKind::SourceCdf) {
    key.value_kind = CompiledMathValueKind::Cdf;
  } else {
    key.value_kind = CompiledMathValueKind::Survival;
  }
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_algebra_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    std::vector<semantic::Index> children,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  if (children.empty()) {
    if (kind == CompiledMathNodeKind::Product) {
      return compiled_math_constant(program, 1.0);
    }
    return compiled_math_constant(program, 0.0);
  }
  if (children.size() == 1U &&
      (kind == CompiledMathNodeKind::Product ||
       kind == CompiledMathNodeKind::Sum ||
       kind == CompiledMathNodeKind::CleanSignedSum)) {
    return children.front();
  }
  CompiledMathNodeKey key;
  key.kind = kind;
  key.value_kind = value_kind;
  key.children = std::move(children);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_unary_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    const semantic::Index child,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  CompiledMathNodeKey key;
  key.kind = kind;
  key.value_kind = value_kind;
  key.children.push_back(child);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_lazy_product_node(
    CompiledMathProgram *program,
    std::vector<semantic::Index> children,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  if (children.empty()) {
    return compiled_math_constant(program, 1.0);
  }
  if (children.size() == 1U) {
    return children.front();
  }
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::LazyProduct;
  key.value_kind = value_kind;
  key.children = std::move(children);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_time_gate_node(
    CompiledMathProgram *program,
    const semantic::Index child,
    const semantic::Index current_time_id,
    const semantic::Index gate_time_id,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::TimeGate;
  key.value_kind = value_kind;
  key.time_id = current_time_id;
  key.aux_id = gate_time_id;
  key.children.push_back(child);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_integral_zero_to_current_node(
    CompiledMathProgram *program,
    const semantic::Index integrand_root_id,
    const semantic::Index condition_id = 0,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::IntegralZeroToCurrent;
  key.value_kind = CompiledMathValueKind::Cdf;
  key.subject_id = integrand_root_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.aux2_id = bind_time_id;
  key.source_view_id = source_view_id;
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_raw_integral_zero_to_current_node(
    CompiledMathProgram *program,
    const semantic::Index integrand_root_id,
    const semantic::Index condition_id = 0,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::IntegralZeroToCurrentRaw;
  key.value_kind = CompiledMathValueKind::Scalar;
  key.subject_id = integrand_root_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.aux2_id = bind_time_id;
  key.source_view_id = source_view_id;
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_root_value_node(
    CompiledMathProgram *program,
    const semantic::Index root_id,
    const semantic::Index source_view_id = 0,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::RootValue;
  key.value_kind = value_kind;
  key.subject_id = root_id;
  key.source_view_id = source_view_id;
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_union_kernel_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    const semantic::Index expr_id,
    const semantic::Index condition_id,
    std::vector<semantic::Index> children,
    const semantic::Index aux_id = semantic::kInvalidIndex,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey key;
  key.kind = kind;
  key.subject_id = expr_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.aux_id = aux_id;
  key.aux2_id = bind_time_id;
  key.source_view_id = source_view_id;
  key.children = std::move(children);
  if (kind == CompiledMathNodeKind::UnionKernelDensity ||
      kind == CompiledMathNodeKind::UnionKernelMultiSubsetDensity) {
    key.value_kind = CompiledMathValueKind::Density;
  } else if (kind == CompiledMathNodeKind::UnionKernelCdf) {
    key.value_kind = CompiledMathValueKind::Cdf;
  } else {
    key.value_kind = CompiledMathValueKind::Scalar;
  }
  return compiled_math_intern_node(program, std::move(key));
}

inline void compiled_math_append_schedule_node(
    const CompiledMathProgram &program,
    const semantic::Index node_id,
    std::vector<std::uint8_t> *visited,
    std::vector<semantic::Index> *schedule) {
  const auto pos = static_cast<std::size_t>(node_id);
  if ((*visited)[pos] != 0U) {
    return;
  }
  (*visited)[pos] = 1U;
  const auto &node = program.nodes[pos];
  if (node.kind == CompiledMathNodeKind::LazyProduct) {
    schedule->push_back(node_id);
    return;
  }
  if (node.condition_id != 0 &&
      node.condition_id != semantic::kInvalidIndex) {
    const auto condition_pos = static_cast<std::size_t>(node.condition_id - 1U);
    if (condition_pos < program.conditions.size()) {
      const auto &condition = program.conditions[condition_pos];
      for (const auto normalizer_node_id :
           condition.fact_normalizer_node_ids) {
        if (normalizer_node_id != semantic::kInvalidIndex) {
          compiled_math_append_schedule_node(
              program, normalizer_node_id, visited, schedule);
        }
      }
    }
  }
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + i)];
    compiled_math_append_schedule_node(program, child_id, visited, schedule);
  }
  schedule->push_back(node_id);
}

inline semantic::Index compiled_math_make_root(CompiledMathProgram *program,
                                              const semantic::Index node_id) {
  const auto root_id = static_cast<semantic::Index>(program->roots.size());
  const auto offset =
      static_cast<semantic::Index>(program->root_schedule_nodes.size());
  std::vector<std::uint8_t> visited(program->nodes.size(), 0U);
  std::vector<semantic::Index> schedule;
  schedule.reserve(program->nodes.size());
  compiled_math_append_schedule_node(*program, node_id, &visited, &schedule);
  program->root_schedule_nodes.insert(
      program->root_schedule_nodes.end(), schedule.begin(), schedule.end());
  program->roots.push_back(
      CompiledMathRoot{
          node_id,
          CompiledMathIndexSpan{
              offset,
              static_cast<semantic::Index>(schedule.size())}});
  return root_id;
}

inline void compiled_math_release_planning_fields(
    CompiledMathProgram *program) {
  decltype(program->node_index)().swap(program->node_index);
  decltype(program->source_channel_index)().swap(
      program->source_channel_index);
  decltype(program->condition_index)().swap(program->condition_index);
}

} // namespace detail
} // namespace accumulatr::eval
