#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iterator>
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
  IntegralZeroToCurrentRaw = 26,
  SourceChannelLoad = 28,
  TimeGate = 29
};

enum class CompiledMathIntegralKernelKind : std::uint8_t {
  SourceProduct = 0,
  SourceProductSum = 1
};

enum class CompiledMathSourceProductOpKind : std::uint8_t {
  ConstantZero = 0,
  ConstantOne = 1,
  LeafLognormalPdf = 2,
  LeafLognormalCdf = 3,
  LeafLognormalSurvival = 4,
  LeafGammaPdf = 5,
  LeafGammaCdf = 6,
  LeafGammaSurvival = 7,
  LeafExgaussPdf = 8,
  LeafExgaussCdf = 9,
  LeafExgaussSurvival = 10,
  LeafLbaPdf = 11,
  LeafLbaCdf = 12,
  LeafLbaSurvival = 13,
  LeafRdmPdf = 14,
  LeafRdmCdf = 15,
  LeafRdmSurvival = 16,
  GenericPdf = 17,
  GenericCdf = 18,
  GenericSurvival = 19
};

inline std::uint8_t compiled_math_source_product_op_channel_mask(
    const CompiledMathSourceProductOpKind kind) noexcept {
  switch (kind) {
  case CompiledMathSourceProductOpKind::LeafLognormalPdf:
  case CompiledMathSourceProductOpKind::LeafGammaPdf:
  case CompiledMathSourceProductOpKind::LeafExgaussPdf:
  case CompiledMathSourceProductOpKind::LeafLbaPdf:
  case CompiledMathSourceProductOpKind::LeafRdmPdf:
  case CompiledMathSourceProductOpKind::GenericPdf:
    return 1U;
  case CompiledMathSourceProductOpKind::LeafLognormalCdf:
  case CompiledMathSourceProductOpKind::LeafGammaCdf:
  case CompiledMathSourceProductOpKind::LeafExgaussCdf:
  case CompiledMathSourceProductOpKind::LeafLbaCdf:
  case CompiledMathSourceProductOpKind::LeafRdmCdf:
  case CompiledMathSourceProductOpKind::GenericCdf:
    return 2U;
  case CompiledMathSourceProductOpKind::LeafLognormalSurvival:
  case CompiledMathSourceProductOpKind::LeafGammaSurvival:
  case CompiledMathSourceProductOpKind::LeafExgaussSurvival:
  case CompiledMathSourceProductOpKind::LeafLbaSurvival:
  case CompiledMathSourceProductOpKind::LeafRdmSurvival:
  case CompiledMathSourceProductOpKind::GenericSurvival:
    return 4U;
  case CompiledMathSourceProductOpKind::ConstantZero:
  case CompiledMathSourceProductOpKind::ConstantOne:
    break;
  }
  return 0U;
}

inline bool compiled_math_source_product_op_is_pdf(
    const CompiledMathSourceProductOpKind kind) noexcept {
  return compiled_math_source_product_op_channel_mask(kind) == 1U;
}

inline bool compiled_math_source_product_op_is_cdf(
    const CompiledMathSourceProductOpKind kind) noexcept {
  return compiled_math_source_product_op_channel_mask(kind) == 2U;
}

inline bool compiled_math_source_product_op_is_survival(
    const CompiledMathSourceProductOpKind kind) noexcept {
  return compiled_math_source_product_op_channel_mask(kind) == 4U;
}

inline bool compiled_math_source_product_op_is_direct_leaf(
    const CompiledMathSourceProductOpKind kind) noexcept {
  return kind >= CompiledMathSourceProductOpKind::LeafLognormalPdf &&
         kind <= CompiledMathSourceProductOpKind::LeafRdmSurvival;
}

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
  CompiledMathIndexSpan source_value_factors{};
  CompiledMathIndexSpan outcome_gate_nodes{};
  CompiledMathIndexSpan time_gate_nodes{};
  CompiledMathIndexSpan integral_factor_nodes{};
  CompiledMathIndexSpan expr_upper_factors{};
  double sign{1.0};
  CompiledMathIndexSpan source_product_ops{};
};

struct CompiledSourceBoundPlan {
  CompiledMathIndexSpan exact{};
  CompiledMathIndexSpan lower{};
  CompiledMathIndexSpan upper{};
  bool has_condition_exact{false};
  bool has_condition_lower{false};
  bool has_condition_upper{false};
  bool use_sequence_exact{true};
  bool use_sequence_lower{true};
  bool use_sequence_upper{true};
};

struct CompiledMathSourceProductChannel {
  semantic::Index source_channel_slot{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  std::uint8_t required_channels{0U};
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index source_kernel_slot{semantic::kInvalidIndex};
  semantic::Index leaf_index{semantic::kInvalidIndex};
  CompiledSourceChannelKernelKind kernel{
      CompiledSourceChannelKernelKind::Invalid};
  std::uint8_t static_source_view_relation{0U};
  std::uint8_t leaf_dist_kind{0U};
  int leaf_param_count{0};
  double leaf_onset_abs_value{0.0};
  CompiledSourceBoundPlan bounds{};
  CompiledMathIndexSpan condition_time_dependencies{};
  semantic::Index condition_static_cache_id{0};
  semantic::Index scalar_op_count{0};
  bool has_static_source_view_relation{false};
  bool direct_leaf_absolute_candidate{false};
  bool has_source_condition_overlay{false};
  bool condition_cache_dynamic{false};
};

enum class CompiledMathIntegralExprUpperMode : std::uint8_t {
  BeforeScale = 0,
  AfterOne = 1
};

struct CompiledMathIntegralExprUpperFactor {
  semantic::Index node_id{semantic::kInvalidIndex};
  CompiledMathIntegralExprUpperMode mode{
      CompiledMathIntegralExprUpperMode::BeforeScale};
};

struct CompiledMathTimedUpperBoundTerm {
  semantic::Index time_id{semantic::kInvalidIndex};
  semantic::Index normalizer_node_id{semantic::kInvalidIndex};
};

struct CompiledMathIntegralKernel {
  CompiledMathIntegralKernelKind kind{
      CompiledMathIntegralKernelKind::SourceProductSum};
  semantic::Index root_id{semantic::kInvalidIndex};
  semantic::Index bind_time_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  CompiledMathIndexSpan source_value_nodes{};
  CompiledMathIndexSpan source_value_factors{};
  CompiledMathIndexSpan source_product_ops{};
  CompiledMathIndexSpan source_product_channels{};
  CompiledMathIndexSpan source_product_terms{};
  bool clean_signed_source_sum{false};
};

struct CompiledMathSourceValueFactor {
  semantic::Index source_channel_slot{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  CompiledMathNodeKind kind{CompiledMathNodeKind::SourceSurvival};
  semantic::Index source_product_channel_id{semantic::kInvalidIndex};
};

struct CompiledMathSourceProductOp {
  CompiledMathSourceProductOpKind kind{
      CompiledMathSourceProductOpKind::GenericSurvival};
  semantic::Index source_product_channel_id{semantic::kInvalidIndex};
  std::uint8_t value_channel_mask{0U};
  std::uint8_t fill_channel_mask{0U};
  double constant_value{0.0};
  bool cache_result{true};
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

struct CompiledSourceBoundTerm {
  semantic::Index time_id{0};
};

struct CompiledConditionCachePlan {
  CompiledMathIndexSpan time_dependencies{};
  semantic::Index static_cache_id{0};
  bool dynamic{false};
};

struct CompiledSourceChannelPlan {
  CompiledMathSourceChannelKey request{};
  CompiledSourceChannelKernelKind kernel{
      CompiledSourceChannelKernelKind::Invalid};
  semantic::Index source_kernel_slot{semantic::kInvalidIndex};
  semantic::Index bound_plan_slot{semantic::kInvalidIndex};
  std::uint8_t required_channels{7U};
  CompiledSourceBoundPlan bounds{};
  CompiledMathIndexSpan condition_time_dependencies{};
  semantic::Index condition_static_cache_id{0};
  semantic::Index source_view_condition_static_cache_id{0};
  bool has_source_condition_overlay{false};
  bool direct_leaf_absolute_candidate{false};
  bool condition_cache_dynamic{false};
  bool source_view_condition_cache_dynamic{false};
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

inline void compiled_math_condition_cache_hash_part(
    std::size_t *seed,
    const std::size_t value) noexcept {
  *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
}

inline semantic::Index compiled_math_condition_cache_index(
    const std::size_t seed) noexcept {
  return static_cast<semantic::Index>(
      (seed & static_cast<std::size_t>(0x3fffffffU)) + 1U);
}

inline semantic::Index compiled_math_static_condition_cache_id(
    const semantic::Index condition_id) noexcept {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return 0;
  }
  std::size_t seed = 0x517cc1b727220a95ULL;
  compiled_math_condition_cache_hash_part(
      &seed, static_cast<std::size_t>(condition_id));
  return compiled_math_condition_cache_index(seed);
}

inline semantic::Index compiled_math_static_source_view_condition_cache_id(
    const semantic::Index source_view_id,
    const semantic::Index condition_id) noexcept {
  const auto evaluator_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return evaluator_id;
  }
  std::size_t seed = 0x9e3779b97f4a7c15ULL;
  compiled_math_condition_cache_hash_part(
      &seed, static_cast<std::size_t>(evaluator_id));
  compiled_math_condition_cache_hash_part(
      &seed, static_cast<std::size_t>(condition_id));
  return compiled_math_condition_cache_index(seed);
}

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
  std::vector<semantic::Index> time_dependency_ids;
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
  std::vector<CompiledMathSourceValueFactor>
      integral_kernel_source_value_factors;
  std::vector<CompiledMathSourceProductChannel>
      integral_kernel_source_product_channels;
  std::vector<CompiledMathSourceProductOp>
      integral_kernel_source_product_ops;
  std::vector<semantic::Index> integral_kernel_outcome_gate_nodes;
  std::vector<semantic::Index> integral_kernel_time_gate_nodes;
  std::vector<semantic::Index> integral_kernel_integral_factor_nodes;
  std::vector<CompiledMathIntegralExprUpperFactor>
      integral_kernel_expr_upper_factors;
  std::vector<CompiledMathTimedUpperBoundTerm> timed_upper_bound_terms;
  std::vector<CompiledMathSourceChannelKey> source_channel_keys;
  std::vector<std::uint8_t> source_channel_required_channels;
  std::vector<CompiledSourceChannelPlan> source_channel_plans;
  std::vector<CompiledSourceBoundPlan> source_condition_bound_plans;
  std::vector<CompiledSourceBoundTerm> source_condition_bound_terms;
  std::vector<CompiledConditionCachePlan> condition_cache_plans;
  std::vector<semantic::Index> condition_cache_time_dependencies;
  std::vector<std::uint8_t> condition_source_relations;
  semantic::Index source_condition_bound_source_count{0};
  semantic::Index condition_source_relation_source_count{0};
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
    const auto source_product_channel_count =
        program.integral_kernel_source_product_channels.size();
    source_product_channel_epoch.assign(source_product_channel_count, 0U);
    source_product_channel_valid_mask.assign(source_product_channel_count, 0U);
    source_product_channel_pdf.assign(source_product_channel_count, 0.0);
    source_product_channel_cdf.assign(source_product_channel_count, 0.0);
    source_product_channel_survival.assign(source_product_channel_count, 1.0);
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
    if (source_product_channel_epoch.size() <
        program.integral_kernel_source_product_channels.size()) {
      const auto size = program.integral_kernel_source_product_channels.size();
      source_product_channel_epoch.resize(size, 0U);
      source_product_channel_valid_mask.resize(size, 0U);
      source_product_channel_pdf.resize(size, 0.0);
      source_product_channel_cdf.resize(size, 0.0);
      source_product_channel_survival.resize(size, 1.0);
    }
  }

  void reset_cache() {
    std::fill(cache_valid.begin(), cache_valid.end(), 0U);
    std::fill(source_channel_valid.begin(), source_channel_valid.end(), 0U);
  }

  void reset_source_product_channel_cache() {
    ++source_product_channel_current_epoch;
    if (source_product_channel_current_epoch == 0U) {
      source_product_channel_current_epoch = 1U;
      std::fill(
          source_product_channel_epoch.begin(),
          source_product_channel_epoch.end(),
          0U);
      std::fill(
          source_product_channel_valid_mask.begin(),
          source_product_channel_valid_mask.end(),
          0U);
    }
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

  class RebindableTimeBinding {
  public:
    RebindableTimeBinding(CompiledMathWorkspace *workspace,
                          const semantic::Index time_id)
        : workspace_(workspace), time_id_(time_id) {
      if (workspace_ == nullptr) {
        return;
      }
      pos_ = static_cast<std::size_t>(time_id_);
      if (workspace_->time_values.size() <= pos_) {
        workspace_->time_values.resize(pos_ + 1U, 0.0);
        workspace_->time_valid.resize(pos_ + 1U, 0U);
      }
      had_previous_ = workspace_->time_valid[pos_] != 0U;
      previous_ = workspace_->time_values[pos_];
    }

    RebindableTimeBinding(const RebindableTimeBinding &) = delete;
    RebindableTimeBinding &operator=(const RebindableTimeBinding &) = delete;

    ~RebindableTimeBinding() {
      if (workspace_ == nullptr) {
        return;
      }
      if (had_previous_) {
        workspace_->time_values[pos_] = previous_;
        workspace_->time_valid[pos_] = 1U;
      } else {
        workspace_->time_values[pos_] = 0.0;
        workspace_->time_valid[pos_] = 0U;
      }
    }

    void set(const double value) {
      if (workspace_ == nullptr) {
        return;
      }
      workspace_->time_values[pos_] = value;
      workspace_->time_valid[pos_] = 1U;
    }

  private:
    CompiledMathWorkspace *workspace_{nullptr};
    semantic::Index time_id_{0};
    std::size_t pos_{0};
    double previous_{0.0};
    bool had_previous_{false};
  };

  CompiledMathWorkspace &integral_workspace_for(
      const CompiledMathProgram &program) {
    if (!integral_workspace) {
      integral_workspace = std::make_unique<CompiledMathWorkspace>(program);
    } else {
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
  std::vector<std::uint32_t> source_product_channel_epoch;
  std::vector<std::uint8_t> source_product_channel_valid_mask;
  std::vector<double> source_product_channel_pdf;
  std::vector<double> source_product_channel_cdf;
  std::vector<double> source_product_channel_survival;
  std::uint32_t source_product_channel_current_epoch{1U};
  std::vector<double> time_values;
  std::vector<std::uint8_t> time_valid;
  std::vector<std::uint8_t> integral_term_open;
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

inline std::uint8_t compiled_math_source_factor_channel_mask(
    const CompiledMathNodeKind kind) noexcept;

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
  const auto required_channels =
      compiled_math_source_factor_channel_mask(node.kind);
  if (found != program->source_channel_index.end()) {
    if (compiled_math_is_source_value_node(node.kind)) {
      const auto slot_pos = static_cast<std::size_t>(found->second);
      if (slot_pos >= program->source_channel_required_channels.size()) {
        program->source_channel_required_channels.resize(slot_pos + 1U, 0U);
      }
      program->source_channel_required_channels[slot_pos] |= required_channels;
    }
    return found->second;
  }
  const auto slot = program->source_channel_count++;
  program->source_channel_keys.push_back(key);
  program->source_channel_required_channels.push_back(required_channels);
  program->source_channel_index.emplace(std::move(key), slot);
  return slot;
}

struct CompiledMathSourceProductTermBuild {
  std::vector<semantic::Index> source_value_nodes;
  std::vector<semantic::Index> outcome_gate_nodes;
  std::vector<semantic::Index> time_gate_nodes;
  std::vector<semantic::Index> integral_factor_nodes;
  std::vector<CompiledMathIntegralExprUpperFactor> expr_upper_factors;
  double sign{1.0};
};

inline bool compiled_math_expand_source_product_terms(
    const CompiledMathProgram &program,
    const semantic::Index node_id,
    const double sign,
    std::vector<CompiledMathSourceProductTermBuild> *terms,
    bool *clean_signed) {
  if (node_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
  if (compiled_math_is_source_value_node(node.kind)) {
    terms->push_back(
        CompiledMathSourceProductTermBuild{{node_id}, {}, {}, {}, {}, sign});
    return true;
  }
  if (node.kind == CompiledMathNodeKind::Constant) {
    const double term_sign = sign * node.constant;
    if (term_sign != 0.0) {
      terms->push_back(
          CompiledMathSourceProductTermBuild{{}, {}, {}, {}, {}, term_sign});
    }
    return true;
  }
  if (node.kind == CompiledMathNodeKind::OutcomeSubsetUnused) {
    terms->push_back(
        CompiledMathSourceProductTermBuild{{}, {node_id}, {}, {}, {}, sign});
    return true;
  }
  if (compiled_math_is_integral_node(node.kind)) {
    if (node.integral_kernel_slot == semantic::kInvalidIndex) {
      return false;
    }
    terms->push_back(
        CompiledMathSourceProductTermBuild{{}, {}, {}, {node_id}, {}, sign});
    return true;
  }
  if ((node.kind == CompiledMathNodeKind::ExprUpperBoundDensity ||
       node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) &&
      node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    std::vector<CompiledMathSourceProductTermBuild> bounded_terms;
    if (!compiled_math_expand_source_product_terms(
            program, child_id, sign, &bounded_terms, clean_signed)) {
      return false;
    }
    for (auto &term : bounded_terms) {
      term.expr_upper_factors.push_back(
          CompiledMathIntegralExprUpperFactor{
              node_id,
              CompiledMathIntegralExprUpperMode::BeforeScale});
    }
    terms->insert(
        terms->end(),
        std::make_move_iterator(bounded_terms.begin()),
        std::make_move_iterator(bounded_terms.end()));
    if (node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
      terms->push_back(
          CompiledMathSourceProductTermBuild{
              {},
              {},
              {},
              {},
              {CompiledMathIntegralExprUpperFactor{
                  node_id,
                  CompiledMathIntegralExprUpperMode::AfterOne}},
              sign});
    }
    return true;
  }
  if (node.kind == CompiledMathNodeKind::TimeGate &&
      node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    std::vector<CompiledMathSourceProductTermBuild> gated_terms;
    if (!compiled_math_expand_source_product_terms(
            program, child_id, sign, &gated_terms, clean_signed)) {
      return false;
    }
    for (auto &term : gated_terms) {
      term.time_gate_nodes.push_back(node_id);
    }
    terms->insert(
        terms->end(),
        std::make_move_iterator(gated_terms.begin()),
        std::make_move_iterator(gated_terms.end()));
    return true;
  }
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
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset + i)];
      if (!compiled_math_expand_source_product_terms(
              program, child_id, sign, terms, clean_signed)) {
        return false;
      }
    }
    return true;
  }
  if (node.kind == CompiledMathNodeKind::Negate && node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    return compiled_math_expand_source_product_terms(
        program, child_id, -sign, terms, clean_signed);
  }
  if (node.kind == CompiledMathNodeKind::ClampProbability &&
      node.children.size == 1U) {
    if (clean_signed != nullptr) {
      *clean_signed = true;
    }
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    return compiled_math_expand_source_product_terms(
        program, child_id, sign, terms, clean_signed);
  }
  if (node.kind == CompiledMathNodeKind::Complement &&
      node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    terms->push_back(CompiledMathSourceProductTermBuild{{}, {}, {}, {}, {}, sign});
    return compiled_math_expand_source_product_terms(
        program, child_id, -sign, terms, clean_signed);
  }
  if (node.kind != CompiledMathNodeKind::Product) {
    throw std::runtime_error(
        "source-product integral collector reached unsupported node kind " +
        std::to_string(static_cast<int>(node.kind)));
  }

  std::vector<CompiledMathSourceProductTermBuild> product_terms(1);
  product_terms.front().sign = sign;
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + i)];
    std::vector<CompiledMathSourceProductTermBuild> child_terms;
    if (!compiled_math_expand_source_product_terms(
            program, child_id, 1.0, &child_terms, clean_signed)) {
      return false;
    }
    std::vector<CompiledMathSourceProductTermBuild> next_terms;
    next_terms.reserve(product_terms.size() * child_terms.size());
    for (const auto &lhs : product_terms) {
      for (const auto &rhs : child_terms) {
        CompiledMathSourceProductTermBuild combined;
        combined.sign = lhs.sign * rhs.sign;
        combined.source_value_nodes.reserve(
            lhs.source_value_nodes.size() + rhs.source_value_nodes.size());
        combined.source_value_nodes.insert(
            combined.source_value_nodes.end(),
            lhs.source_value_nodes.begin(),
            lhs.source_value_nodes.end());
        combined.source_value_nodes.insert(
            combined.source_value_nodes.end(),
            rhs.source_value_nodes.begin(),
            rhs.source_value_nodes.end());
        combined.outcome_gate_nodes.reserve(
            lhs.outcome_gate_nodes.size() + rhs.outcome_gate_nodes.size());
        combined.outcome_gate_nodes.insert(
            combined.outcome_gate_nodes.end(),
            lhs.outcome_gate_nodes.begin(),
            lhs.outcome_gate_nodes.end());
        combined.outcome_gate_nodes.insert(
            combined.outcome_gate_nodes.end(),
            rhs.outcome_gate_nodes.begin(),
            rhs.outcome_gate_nodes.end());
        combined.time_gate_nodes.reserve(
            lhs.time_gate_nodes.size() + rhs.time_gate_nodes.size());
        combined.time_gate_nodes.insert(
            combined.time_gate_nodes.end(),
            lhs.time_gate_nodes.begin(),
            lhs.time_gate_nodes.end());
        combined.time_gate_nodes.insert(
            combined.time_gate_nodes.end(),
            rhs.time_gate_nodes.begin(),
            rhs.time_gate_nodes.end());
        combined.integral_factor_nodes.reserve(
            lhs.integral_factor_nodes.size() +
            rhs.integral_factor_nodes.size());
        combined.integral_factor_nodes.insert(
            combined.integral_factor_nodes.end(),
            lhs.integral_factor_nodes.begin(),
            lhs.integral_factor_nodes.end());
        combined.integral_factor_nodes.insert(
            combined.integral_factor_nodes.end(),
            rhs.integral_factor_nodes.begin(),
            rhs.integral_factor_nodes.end());
        combined.expr_upper_factors.reserve(
            lhs.expr_upper_factors.size() + rhs.expr_upper_factors.size());
        combined.expr_upper_factors.insert(
            combined.expr_upper_factors.end(),
            lhs.expr_upper_factors.begin(),
            lhs.expr_upper_factors.end());
        combined.expr_upper_factors.insert(
            combined.expr_upper_factors.end(),
            rhs.expr_upper_factors.begin(),
            rhs.expr_upper_factors.end());
        if (combined.sign != 0.0) {
          next_terms.push_back(std::move(combined));
        }
      }
    }
    product_terms = std::move(next_terms);
  }
  terms->insert(terms->end(), product_terms.begin(), product_terms.end());
  return true;
}

inline bool compiled_math_collect_source_product_terms(
    CompiledMathProgram *program,
    const semantic::Index node_id,
    const double sign,
    std::vector<CompiledMathIntegralSourceProductTerm> *terms,
    bool *clean_signed) {
  std::vector<CompiledMathSourceProductTermBuild> built_terms;
  if (!compiled_math_expand_source_product_terms(
          *program, node_id, sign, &built_terms, clean_signed)) {
    return false;
  }
  if (built_terms.empty()) {
    terms->push_back(
        CompiledMathIntegralSourceProductTerm{
            CompiledMathIndexSpan{},
            CompiledMathIndexSpan{},
            CompiledMathIndexSpan{},
            CompiledMathIndexSpan{},
            CompiledMathIndexSpan{},
            CompiledMathIndexSpan{},
            0.0});
    return true;
  }
  for (const auto &built : built_terms) {
    const auto offset = static_cast<semantic::Index>(
        program->integral_kernel_source_value_nodes.size());
    program->integral_kernel_source_value_nodes.insert(
        program->integral_kernel_source_value_nodes.end(),
        built.source_value_nodes.begin(),
        built.source_value_nodes.end());
    const auto factor_offset = static_cast<semantic::Index>(
        program->integral_kernel_source_value_factors.size());
    for (const auto source_node_id : built.source_value_nodes) {
      const auto &source_node =
          program->nodes[static_cast<std::size_t>(source_node_id)];
      if (!compiled_math_is_source_value_node(source_node.kind) ||
          source_node.source_channel_slot == semantic::kInvalidIndex) {
        throw std::runtime_error(
            "source-product integral factor is not a planned source value");
      }
      program->integral_kernel_source_value_factors.push_back(
          CompiledMathSourceValueFactor{
              source_node.source_channel_slot,
              source_node.source_view_id,
              source_node.time_id,
              source_node.aux_id,
              source_node.kind});
    }
    const auto gate_offset = static_cast<semantic::Index>(
        program->integral_kernel_outcome_gate_nodes.size());
    program->integral_kernel_outcome_gate_nodes.insert(
        program->integral_kernel_outcome_gate_nodes.end(),
        built.outcome_gate_nodes.begin(),
        built.outcome_gate_nodes.end());
    const auto time_gate_offset = static_cast<semantic::Index>(
        program->integral_kernel_time_gate_nodes.size());
    program->integral_kernel_time_gate_nodes.insert(
        program->integral_kernel_time_gate_nodes.end(),
        built.time_gate_nodes.begin(),
        built.time_gate_nodes.end());
    const auto integral_factor_offset = static_cast<semantic::Index>(
        program->integral_kernel_integral_factor_nodes.size());
    program->integral_kernel_integral_factor_nodes.insert(
        program->integral_kernel_integral_factor_nodes.end(),
        built.integral_factor_nodes.begin(),
        built.integral_factor_nodes.end());
    const auto expr_upper_offset = static_cast<semantic::Index>(
        program->integral_kernel_expr_upper_factors.size());
    program->integral_kernel_expr_upper_factors.insert(
        program->integral_kernel_expr_upper_factors.end(),
        built.expr_upper_factors.begin(),
        built.expr_upper_factors.end());
    terms->push_back(
        CompiledMathIntegralSourceProductTerm{
            CompiledMathIndexSpan{
                offset,
                static_cast<semantic::Index>(
                    built.source_value_nodes.size())},
            CompiledMathIndexSpan{
                factor_offset,
                static_cast<semantic::Index>(
                    built.source_value_nodes.size())},
            CompiledMathIndexSpan{
                gate_offset,
                static_cast<semantic::Index>(
                    built.outcome_gate_nodes.size())},
            CompiledMathIndexSpan{
                time_gate_offset,
                static_cast<semantic::Index>(
                    built.time_gate_nodes.size())},
            CompiledMathIndexSpan{
                integral_factor_offset,
                static_cast<semantic::Index>(
                    built.integral_factor_nodes.size())},
            CompiledMathIndexSpan{
                expr_upper_offset,
                static_cast<semantic::Index>(
                    built.expr_upper_factors.size())},
            built.sign});
  }
  return true;
}

inline bool compiled_math_same_source_product_channel(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathSourceValueFactor &factor,
    const semantic::Index effective_source_view_id) noexcept {
  return channel.source_channel_slot == factor.source_channel_slot &&
         channel.source_view_id == effective_source_view_id &&
         channel.time_id == factor.time_id &&
         channel.time_cap_id == factor.time_cap_id;
}

inline std::uint8_t compiled_math_source_factor_channel_mask(
    const CompiledMathNodeKind kind) noexcept {
  switch (kind) {
  case CompiledMathNodeKind::SourcePdf:
    return 1U;
  case CompiledMathNodeKind::SourceCdf:
    return 2U;
  case CompiledMathNodeKind::SourceSurvival:
    return 4U;
  default:
    break;
  }
  return 0U;
}

inline void compiled_math_compile_source_product_channels(
    CompiledMathProgram *program,
    CompiledMathIntegralKernel *kernel,
    const std::size_t source_value_factor_mark) {
  const auto channel_offset = static_cast<semantic::Index>(
      program->integral_kernel_source_product_channels.size());
  for (std::size_t factor_pos = source_value_factor_mark;
       factor_pos < program->integral_kernel_source_value_factors.size();
       ++factor_pos) {
    auto &factor = program->integral_kernel_source_value_factors[factor_pos];
    const auto required_channels =
        compiled_math_source_factor_channel_mask(factor.kind);
    const auto factor_source_view_id =
        factor.source_view_id == semantic::kInvalidIndex ? 0 : factor.source_view_id;
    const auto effective_source_view_id =
        factor_source_view_id == 0 ? kernel->source_view_id : factor_source_view_id;
    semantic::Index channel_id = semantic::kInvalidIndex;
    for (std::size_t channel_pos = static_cast<std::size_t>(channel_offset);
         channel_pos < program->integral_kernel_source_product_channels.size();
         ++channel_pos) {
      if (compiled_math_same_source_product_channel(
              program->integral_kernel_source_product_channels[channel_pos],
              factor,
              effective_source_view_id)) {
        channel_id = static_cast<semantic::Index>(channel_pos);
        break;
      }
    }
    if (channel_id == semantic::kInvalidIndex) {
      channel_id = static_cast<semantic::Index>(
          program->integral_kernel_source_product_channels.size());
      program->integral_kernel_source_product_channels.push_back(
          CompiledMathSourceProductChannel{
              factor.source_channel_slot,
              effective_source_view_id,
              factor.time_id,
              factor.time_cap_id,
              required_channels});
    } else {
      program->integral_kernel_source_product_channels[
          static_cast<std::size_t>(channel_id)].required_channels |=
          required_channels;
    }
    factor.source_product_channel_id = channel_id;
  }
  kernel->source_product_channels = CompiledMathIndexSpan{
      channel_offset,
      static_cast<semantic::Index>(
          program->integral_kernel_source_product_channels.size() -
          static_cast<std::size_t>(channel_offset))};
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
  kernel.bind_time_id = compiled_math_integral_bind_time_id(node);
  kernel.source_view_id =
      node.source_view_id == semantic::kInvalidIndex ? 0 : node.source_view_id;

  std::vector<CompiledMathIntegralSourceProductTerm> source_product_terms;
  bool clean_signed_source_sum = false;
  const auto source_value_node_mark =
      program->integral_kernel_source_value_nodes.size();
  const auto source_value_factor_mark =
      program->integral_kernel_source_value_factors.size();
  const auto source_product_channel_mark =
      program->integral_kernel_source_product_channels.size();
  const auto outcome_gate_node_mark =
      program->integral_kernel_outcome_gate_nodes.size();
  const auto time_gate_node_mark =
      program->integral_kernel_time_gate_nodes.size();
  const auto integral_factor_node_mark =
      program->integral_kernel_integral_factor_nodes.size();
  const auto expr_upper_factor_mark =
      program->integral_kernel_expr_upper_factors.size();
  if (compiled_math_collect_source_product_terms(
          program,
          root.node_id,
          1.0,
          &source_product_terms,
          &clean_signed_source_sum)) {
    if (source_product_terms.size() == 1U &&
        source_product_terms.front().sign == 1.0 &&
        source_product_terms.front().outcome_gate_nodes.empty() &&
        source_product_terms.front().time_gate_nodes.empty() &&
        source_product_terms.front().integral_factor_nodes.empty() &&
        source_product_terms.front().expr_upper_factors.empty() &&
        !clean_signed_source_sum) {
      kernel.kind = CompiledMathIntegralKernelKind::SourceProduct;
      kernel.source_value_nodes =
          source_product_terms.front().source_value_nodes;
      kernel.source_value_factors =
          source_product_terms.front().source_value_factors;
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
    compiled_math_compile_source_product_channels(
        program,
        &kernel,
        source_value_factor_mark);
  } else {
    program->integral_kernel_source_value_nodes.resize(
        source_value_node_mark);
    program->integral_kernel_source_value_factors.resize(
        source_value_factor_mark);
    program->integral_kernel_source_product_channels.resize(
        source_product_channel_mark);
    program->integral_kernel_outcome_gate_nodes.resize(
        outcome_gate_node_mark);
    program->integral_kernel_time_gate_nodes.resize(time_gate_node_mark);
    program->integral_kernel_integral_factor_nodes.resize(
        integral_factor_node_mark);
    program->integral_kernel_expr_upper_factors.resize(
        expr_upper_factor_mark);
    throw std::runtime_error(
        "compiled integral root is not a planned source-product kernel");
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
  condition->time_dependency_ids.clear();
  condition->source_order_fact_lookup.entries.clear();
  condition->source_order_fact_lookup.fact_indices.clear();
  condition->guard_upper_fact_lookup.entries.clear();
  condition->guard_upper_fact_lookup.fact_indices.clear();

  for (std::size_t i = 0; i < condition->fact_kinds.size(); ++i) {
    const auto time_id = condition->fact_time_ids[i];
    if (time_id != semantic::kInvalidIndex &&
        std::find(
            condition->time_dependency_ids.begin(),
            condition->time_dependency_ids.end(),
            time_id) == condition->time_dependency_ids.end()) {
      condition->time_dependency_ids.push_back(time_id);
    }
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
