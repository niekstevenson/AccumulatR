#pragma once

#include <cstdint>
#include <functional>
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
  IntegralZeroToCurrent = 17,
  ExprUpperBoundDensity = 18,
  ExprUpperBoundCdf = 19,
  OutcomeSubsetUnused = 20,
  OutcomeSubsetUsed = 21,
  IntegralZeroToCurrentRaw = 22,
  TimeGate = 29,
  StrictTimeGate = 30
};

enum class CompiledMathIntegralKernelKind : std::uint8_t {
  SourceProduct = 0,
  SourceProductSum = 1,
  Generic = 2
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

enum class CompiledMathSourceProductProgramKind : std::uint8_t {
  ConstantZero = 0,
  ConstantOne = 1,
  LeafAbsolute = 2,
  ExactGate = 3,
  Conditioned = 4,
  OnsetConvolution = 5,
  PoolKOfN = 6
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
  semantic::Index source_program_id{semantic::kInvalidIndex};
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
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index source_view_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  std::uint8_t required_channels{0U};
  semantic::Index source_product_program_id{semantic::kInvalidIndex};
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

struct CompiledMathSourceProductProgram {
  CompiledMathSourceProductProgramKind kind{
      CompiledMathSourceProductProgramKind::ConstantZero};
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index time_id{semantic::kInvalidIndex};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  semantic::Index child_program_id{semantic::kInvalidIndex};
  semantic::Index onset_source_program_id{semantic::kInvalidIndex};
  CompiledMathIndexSpan member_programs{};
  CompiledSourceBoundPlan bounds{};
  CompiledSourceBoundPlan onset_bounds{};
  semantic::Index leaf_index{semantic::kInvalidIndex};
  std::uint8_t leaf_dist_kind{0U};
  double leaf_onset_abs_value{0.0};
  double leaf_onset_lag{0.0};
  semantic::Index pool_k{0};
  semantic::Index source_product_scratch_offset{0};
  semantic::Index source_product_scratch_size{0};
  std::uint8_t static_source_view_relation{0U};
  bool has_static_source_view_relation{false};
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
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
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
  semantic::Index source_product_program_id{semantic::kInvalidIndex};
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

struct CompiledSourceBoundTerm {
  semantic::Index time_id{0};
};

struct CompiledConditionCachePlan {
  CompiledMathIndexSpan time_dependencies{};
  semantic::Index static_cache_id{0};
  bool dynamic{false};
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
  SourceOrder = 4
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
  std::vector<CompiledMathSourceProductProgram>
      integral_kernel_source_product_programs;
  std::vector<semantic::Index>
      integral_kernel_source_product_program_members;
  semantic::Index integral_kernel_source_product_scratch_size{0};
  std::vector<semantic::Index> integral_kernel_outcome_gate_nodes;
  std::vector<semantic::Index> integral_kernel_time_gate_nodes;
  std::vector<semantic::Index> integral_kernel_integral_factor_nodes;
  std::vector<CompiledMathIntegralExprUpperFactor>
      integral_kernel_expr_upper_factors;
  std::vector<CompiledMathTimedUpperBoundTerm> timed_upper_bound_terms;
  std::vector<CompiledSourceBoundPlan> source_condition_bound_plans;
  std::vector<CompiledSourceBoundTerm> source_condition_bound_terms;
  std::vector<CompiledConditionCachePlan> condition_cache_plans;
  std::vector<semantic::Index> condition_cache_time_dependencies;
  std::vector<std::uint8_t> condition_source_relations;
  semantic::Index source_condition_bound_source_count{0};
  semantic::Index condition_source_relation_source_count{0};
  std::unordered_map<
      CompiledMathNodeKey,
      semantic::Index,
      CompiledMathNodeKeyHash>
      node_index;
  std::unordered_map<
      CompiledMathConditionKey,
      semantic::Index,
      CompiledMathConditionKeyHash>
      condition_index;
};

} // namespace detail
} // namespace accumulatr::eval
