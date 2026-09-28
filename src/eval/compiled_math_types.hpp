#pragma once

#include <cstdint>
#include <array>
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
  Constant,
  SourcePdf,
  SourceCdf,
  SourceSurvival,
  Product,
  Sum,
  CleanSignedSum,
  ClampProbability,
  Complement,
  Negate,
  IntegralZeroToCurrent,
  ExprUpperBoundDensity,
  ExprUpperBoundCdf,
  OutcomeSelect,
  IntegralZeroToCurrentRaw,
  TimeGate,
  StrictTimeGate
};

enum class CompiledMathExecutionKind : std::uint8_t {
  SourceProduct = 0,
  SourceProductSum = 1,
  Schedule = 2
};

enum class CompiledMathSourceProductProgramKind : std::uint8_t {
  ConstantZero = 0,
  LeafAbsolute = 1,
  ExactGate = 2,
  Conditioned = 3,
  OnsetConvolution = 4,
  PoolKOfN = 5
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

struct CompiledMathSourceProductOps {
  semantic::Index offset{0};
  semantic::Index size{0};
  bool can_overflow{false};

  [[nodiscard]] bool empty() const noexcept {
    return size == 0;
  }
};

struct CompiledMathNode {
  CompiledMathNodeKind kind{CompiledMathNodeKind::Constant};
  semantic::Index subject_id{semantic::kInvalidIndex};
  semantic::Index time_id{0};
  semantic::Index aux_id{semantic::kInvalidIndex};
  semantic::Index aux2_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  CompiledMathIndexSpan children{};
  std::array<semantic::Index, 2> branch_roots{
      semantic::kInvalidIndex, semantic::kInvalidIndex};
  semantic::Index source_program_id{semantic::kInvalidIndex};
  semantic::Index integral_kernel_slot{semantic::kInvalidIndex};
  bool cache_source_program{false};
  double constant{0.0};
};

struct CompiledMathSourceProductTerm {
  CompiledMathIndexSpan source_value_factors{};
  CompiledMathIndexSpan time_gate_nodes{};
  CompiledMathIndexSpan integral_factor_nodes{};
  CompiledMathIndexSpan expr_upper_factors{};
  double sign{1.0};
  CompiledMathSourceProductOps source_product_ops{};
};

struct CompiledMathSourceProductChannel {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  std::uint8_t required_channels{0U};
  semantic::Index source_product_program_id{semantic::kInvalidIndex};
  std::uint8_t static_source_view_relation{0U};
};

constexpr semantic::Index kInitialCertainSourceProgramId{-2};

struct CompiledMathSourceProductProgram {
  CompiledMathSourceProductProgramKind kind{
      CompiledMathSourceProductProgramKind::ConstantZero};
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index child_program_id{semantic::kInvalidIndex};
  semantic::Index onset_source_program_id{semantic::kInvalidIndex};
  CompiledMathIndexSpan member_programs{};
  semantic::Index leaf_index{semantic::kInvalidIndex};
  std::uint8_t leaf_dist_kind{0U};
  double leaf_onset_lag{0.0};
  semantic::Index pool_k{0};
  std::uint8_t static_source_view_relation{0U};
  semantic::Index initial_without_pdf_program_id{semantic::kInvalidIndex};
  semantic::Index initial_with_pdf_program_id{semantic::kInvalidIndex};
};

enum class CompiledMathIntegralExprUpperMode : std::uint8_t {
  BeforeScale = 0,
  AfterOne = 1
};

struct CompiledMathExprUpperFactor {
  semantic::Index node_id{semantic::kInvalidIndex};
  CompiledMathIntegralExprUpperMode mode{
      CompiledMathIntegralExprUpperMode::BeforeScale};
};

struct CompiledMathExecutionPlan {
  CompiledMathExecutionKind kind{CompiledMathExecutionKind::Schedule};
  CompiledMathIndexSpan source_value_factors{};
  CompiledMathSourceProductOps source_product_ops{};
  CompiledMathIndexSpan source_product_terms{};
  bool clean_signed_source_sum{false};
};

struct CompiledMathRoot {
  semantic::Index node_id{semantic::kInvalidIndex};
  CompiledMathIndexSpan schedule{};
  CompiledMathExecutionPlan execution{};
  CompiledMathExecutionPlan initial_execution{};
};

struct CompiledMathIntegralKernel {
  semantic::Index root_id{semantic::kInvalidIndex};
  semantic::Index bind_time_id{semantic::kInvalidIndex};
  CompiledMathExecutionPlan execution{};
  CompiledMathExecutionPlan initial_execution{};
  // Leaves determining a univariate cumulative integrand (empty if ineligible).
  std::vector<semantic::Index> cumulative_leaves;
  std::vector<semantic::Index> support_leaves;
};

struct CompiledMathSourceValueFactor {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  CompiledMathNodeKind kind{CompiledMathNodeKind::SourceSurvival};
  semantic::Index source_product_channel_id{semantic::kInvalidIndex};
};

struct CompiledMathSourceProductOp {
  semantic::Index source_product_program_id{semantic::kInvalidIndex};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};
  std::uint8_t value_channel_mask{0U};
  std::uint8_t fill_channel_mask{0U};
  bool cache_result{false};
};

struct CompiledMathNodeKey {
  CompiledMathNodeKind kind{CompiledMathNodeKind::Constant};
  CompiledMathValueKind value_kind{CompiledMathValueKind::Scalar};
  semantic::Index subject_id{semantic::kInvalidIndex};
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

enum class CompiledMathTimeSlot : semantic::Index {
  Observed = 0,
  Active = 1,
  Zero = 2
};

struct CompiledMathProgram {
  std::vector<CompiledMathNode> nodes;
  std::vector<semantic::Index> child_nodes;
  std::vector<CompiledMathRoot> roots;
  std::vector<semantic::Index> root_schedule_nodes;
  std::vector<CompiledMathIntegralKernel> integral_kernels;
  std::vector<CompiledMathSourceProductTerm> source_product_terms;
  std::vector<CompiledMathSourceValueFactor>
      source_value_factors;
  std::vector<CompiledMathSourceProductChannel>
      source_product_channels;
  std::vector<CompiledMathSourceProductOp>
      source_product_ops;
  std::vector<CompiledMathSourceProductProgram>
      source_programs;
  std::vector<semantic::Index>
      source_program_members;
  std::vector<semantic::Index> source_program_cache_slots;
  semantic::Index source_program_cache_count{0};
  std::vector<semantic::Index> time_gate_nodes;
  std::vector<semantic::Index> integral_factor_nodes;
  std::vector<CompiledMathExprUpperFactor> expr_upper_factors;
  semantic::Index time_slot_count{
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) + 1U};
  std::unordered_map<
      CompiledMathNodeKey,
      semantic::Index,
      CompiledMathNodeKeyHash>
      node_index;
};

} // namespace detail
} // namespace accumulatr::eval
