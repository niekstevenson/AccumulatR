// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include <array>
#include <string>
#include <vector>

#include "eval/likelihood_context.hpp"
#include "eval/observed_kernel.hpp"
#include "leaf/dist_kind.hpp"

namespace {

using accumulatr::eval::detail::CompiledMathIndexSpan;
using accumulatr::eval::detail::CompiledMathIntegralKernel;
using accumulatr::eval::detail::CompiledMathIntegralKernelKind;
using accumulatr::eval::detail::CompiledMathNode;
using accumulatr::eval::detail::CompiledMathNodeKind;
using accumulatr::eval::detail::CompiledMathProgram;
using accumulatr::eval::detail::CompiledMathSourceProductOpKind;
using accumulatr::eval::detail::CompiledMathSourceProductProgramKind;
using accumulatr::eval::detail::CompiledSourceBoundPlan;
using accumulatr::eval::detail::ExactProbabilityOp;
using accumulatr::eval::detail::ExactProbabilityOpKind;
using accumulatr::eval::detail::ExactProbabilityProgram;
using accumulatr::eval::detail::ExactVariantPlan;
using accumulatr::eval::detail::NativeLikelihoodContext;
using accumulatr::semantic::Index;

constexpr std::array<const char *, 6> kBatchCoverageCategories{{
    "BatchComplete",
    "BatchGroupedButScalarLeafMath",
    "BatchGroupedButScalarBounds",
    "BatchGroupedButScalarExprUpper",
    "ScalarIdentityShortcut",
    "Unsupported",
}};

struct BatchCoverageRows {
  std::vector<std::string> scope;
  std::vector<int> variant_index;
  std::vector<int> program_index;
  std::vector<int> root_id;
  std::vector<std::string> category;
  std::vector<std::string> reason;

  void add(const char *scope_value,
           const int variant_value,
           const int program_value,
           const int root_value,
           const char *category_value,
           std::string reason_value) {
    for (std::size_t i = 0; i < category.size(); ++i) {
      if (scope[i] == scope_value &&
          variant_index[i] == variant_value &&
          program_index[i] == program_value &&
          root_id[i] == root_value &&
          category[i] == category_value) {
        return;
      }
    }
    scope.emplace_back(scope_value);
    variant_index.push_back(variant_value);
    program_index.push_back(program_value);
    root_id.push_back(root_value);
    category.emplace_back(category_value);
    reason.emplace_back(std::move(reason_value));
  }

  std::size_t size() const noexcept {
    return category.size();
  }
};

inline int r_na_int() noexcept {
  return Rcpp::IntegerVector::get_na();
}

inline int coverage_index(const Index value) noexcept {
  return value == accumulatr::semantic::kInvalidIndex
             ? r_na_int()
             : static_cast<int>(value);
}

inline bool coverage_bound_plan_has_work(
    const CompiledSourceBoundPlan &bounds) noexcept {
  return !bounds.exact.empty() || !bounds.lower.empty() ||
         !bounds.upper.empty() || bounds.has_condition_exact ||
         bounds.has_condition_lower || bounds.has_condition_upper ||
         !bounds.use_sequence_exact || !bounds.use_sequence_lower ||
         !bounds.use_sequence_upper;
}

inline const char *coverage_node_kind_name(
    const CompiledMathNodeKind kind) noexcept {
  switch (kind) {
  case CompiledMathNodeKind::Constant:
    return "Constant";
  case CompiledMathNodeKind::SourcePdf:
    return "SourcePdf";
  case CompiledMathNodeKind::SourceCdf:
    return "SourceCdf";
  case CompiledMathNodeKind::SourceSurvival:
    return "SourceSurvival";
  case CompiledMathNodeKind::ExprDensity:
    return "ExprDensity";
  case CompiledMathNodeKind::ExprCdf:
    return "ExprCdf";
  case CompiledMathNodeKind::ExprSurvival:
    return "ExprSurvival";
  case CompiledMathNodeKind::Product:
    return "Product";
  case CompiledMathNodeKind::Sum:
    return "Sum";
  case CompiledMathNodeKind::CleanSignedSum:
    return "CleanSignedSum";
  case CompiledMathNodeKind::ClampProbability:
    return "ClampProbability";
  case CompiledMathNodeKind::Complement:
    return "Complement";
  case CompiledMathNodeKind::Negate:
    return "Negate";
  case CompiledMathNodeKind::SimpleGuardDensity:
    return "SimpleGuardDensity";
  case CompiledMathNodeKind::SimpleGuardCdf:
    return "SimpleGuardCdf";
  case CompiledMathNodeKind::IntegralZeroToCurrent:
    return "IntegralZeroToCurrent";
  case CompiledMathNodeKind::ExprUpperBoundDensity:
    return "ExprUpperBoundDensity";
  case CompiledMathNodeKind::ExprUpperBoundCdf:
    return "ExprUpperBoundCdf";
  case CompiledMathNodeKind::UnionKernelDensity:
    return "UnionKernelDensity";
  case CompiledMathNodeKind::UnionKernelCdf:
    return "UnionKernelCdf";
  case CompiledMathNodeKind::UnionKernelMultiSubsetDensity:
    return "UnionKernelMultiSubsetDensity";
  case CompiledMathNodeKind::UnionKernelMultiSubsetCdf:
    return "UnionKernelMultiSubsetCdf";
  case CompiledMathNodeKind::OutcomeSubsetUnused:
    return "OutcomeSubsetUnused";
  case CompiledMathNodeKind::IntegralZeroToCurrentRaw:
    return "IntegralZeroToCurrentRaw";
  case CompiledMathNodeKind::TimeGate:
    return "TimeGate";
  }
  return "UnknownNode";
}

inline const char *coverage_program_kind_name(
    const CompiledMathSourceProductProgramKind kind) noexcept {
  switch (kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
    return "ConstantZero";
  case CompiledMathSourceProductProgramKind::ConstantOne:
    return "ConstantOne";
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    return "LeafAbsolute";
  case CompiledMathSourceProductProgramKind::ExactGate:
    return "ExactGate";
  case CompiledMathSourceProductProgramKind::Conditioned:
    return "Conditioned";
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
    return "OnsetConvolution";
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    return "PoolKOfN";
  }
  return "UnknownSourceProgram";
}

inline const char *coverage_probability_op_kind_name(
    const ExactProbabilityOpKind kind) noexcept {
  switch (kind) {
  case ExactProbabilityOpKind::Constant:
    return "Constant";
  case ExactProbabilityOpKind::Top1LeafRaceDensity:
    return "Top1LeafRaceDensity";
  case ExactProbabilityOpKind::TerminalNoResponseProbability:
    return "TerminalNoResponseProbability";
  case ExactProbabilityOpKind::GenericTransitionDensity:
    return "GenericTransitionDensity";
  case ExactProbabilityOpKind::GenericTransitionProbability:
    return "GenericTransitionProbability";
  case ExactProbabilityOpKind::RankedTransitionSequence:
    return "RankedTransitionSequence";
  case ExactProbabilityOpKind::Integral:
    return "Integral";
  case ExactProbabilityOpKind::WeightedTriggerSum:
    return "WeightedTriggerSum";
  case ExactProbabilityOpKind::Log:
    return "Log";
  }
  return "UnknownProbabilityOp";
}

inline const char *coverage_dist_kind_name(const std::uint8_t kind) noexcept {
  switch (static_cast<accumulatr::leaf::DistKind>(kind)) {
  case accumulatr::leaf::DistKind::Lognormal:
    return "lognormal";
  case accumulatr::leaf::DistKind::Gamma:
    return "gamma";
  case accumulatr::leaf::DistKind::Exgauss:
    return "exgauss";
  case accumulatr::leaf::DistKind::LBA:
    return "LBA";
  case accumulatr::leaf::DistKind::RDM:
    return "RDM";
  }
  return "unknown";
}

inline bool coverage_leaf_math_is_scalar(
    const std::uint8_t leaf_dist_kind,
    const std::uint8_t channel_mask) noexcept {
#if !defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  (void)leaf_dist_kind;
  (void)channel_mask;
  return true;
#else
  switch (static_cast<accumulatr::leaf::DistKind>(leaf_dist_kind)) {
  case accumulatr::leaf::DistKind::Lognormal:
  case accumulatr::leaf::DistKind::Exgauss:
    return false;
  case accumulatr::leaf::DistKind::Gamma:
    return (channel_mask & ~1U) != 0U;
  case accumulatr::leaf::DistKind::LBA:
  case accumulatr::leaf::DistKind::RDM:
    return true;
  }
  return true;
#endif
}

inline bool coverage_source_product_op_is_generic(
    const CompiledMathSourceProductOpKind kind) noexcept {
  return kind == CompiledMathSourceProductOpKind::GenericPdf ||
         kind == CompiledMathSourceProductOpKind::GenericCdf ||
         kind == CompiledMathSourceProductOpKind::GenericSurvival;
}

void scan_source_product_program_coverage(const CompiledMathProgram &program,
                                          Index source_program_id,
                                          int variant_index,
                                          int probability_program_index,
                                          int root_id,
                                          std::size_t depth,
                                          BatchCoverageRows *rows);

void scan_source_product_ops_coverage(const CompiledMathProgram &program,
                                      const CompiledMathIndexSpan ops,
                                      const int variant_index,
                                      const int probability_program_index,
                                      const int root_id,
                                      const std::size_t depth,
                                      BatchCoverageRows *rows) {
  if (depth > 32U) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        root_id,
        "Unsupported",
        "source-product op recursion exceeded the coverage depth limit");
    return;
  }
  for (Index i = 0; i < ops.size; ++i) {
    const auto op_index = static_cast<std::size_t>(ops.offset + i);
    if (op_index >= program.integral_kernel_source_product_ops.size()) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          root_id,
          "Unsupported",
          "source-product op span points outside the compiled op table");
      continue;
    }
    const auto &op = program.integral_kernel_source_product_ops[op_index];
    if (op.value_channel_mask == 0U) {
      continue;
    }
    if (coverage_source_product_op_is_generic(op.kind)) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          root_id,
          "BatchGroupedButScalarLeafMath",
          "generic source-product op has no distribution-specific vector kernel");
    }
    if (op.source_product_channel_id != accumulatr::semantic::kInvalidIndex) {
      const auto channel_index =
          static_cast<std::size_t>(op.source_product_channel_id);
      if (channel_index < program.integral_kernel_source_product_channels.size()) {
        const auto &channel =
            program.integral_kernel_source_product_channels[channel_index];
        if (coverage_bound_plan_has_work(channel.bounds)) {
          rows->add(
              "program",
              variant_index,
              probability_program_index,
              root_id,
              "BatchGroupedButScalarBounds",
              "source-product channel resolves bounds per lane");
        }
        if (channel.leaf_index != accumulatr::semantic::kInvalidIndex &&
            coverage_leaf_math_is_scalar(
                channel.leaf_dist_kind,
                op.value_channel_mask)) {
          rows->add(
              "program",
              variant_index,
              probability_program_index,
              root_id,
              "BatchGroupedButScalarLeafMath",
              std::string("direct leaf math is not vector-native for ") +
                  coverage_dist_kind_name(channel.leaf_dist_kind));
        }
      } else {
        rows->add(
            "program",
            variant_index,
            probability_program_index,
            root_id,
            "Unsupported",
            "source-product op points outside the compiled channel table");
      }
    }
    scan_source_product_program_coverage(
        program,
        op.source_product_program_id,
        variant_index,
        probability_program_index,
        root_id,
        depth + 1U,
        rows);
  }
}

void scan_source_product_program_coverage(const CompiledMathProgram &program,
                                          const Index source_program_id,
                                          const int variant_index,
                                          const int probability_program_index,
                                          const int root_id,
                                          const std::size_t depth,
                                          BatchCoverageRows *rows) {
  if (source_program_id == accumulatr::semantic::kInvalidIndex) {
    return;
  }
  if (depth > 32U) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        root_id,
        "Unsupported",
        "source-product program recursion exceeded the coverage depth limit");
    return;
  }
  const auto source_program_index = static_cast<std::size_t>(source_program_id);
  if (source_program_index >=
      program.integral_kernel_source_product_programs.size()) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        root_id,
        "Unsupported",
        "source-product op points outside the compiled source-program table");
    return;
  }
  const auto &source_program =
      program.integral_kernel_source_product_programs[source_program_index];
  if (coverage_bound_plan_has_work(source_program.bounds) ||
      coverage_bound_plan_has_work(source_program.onset_bounds)) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        root_id,
        "BatchGroupedButScalarBounds",
        std::string("source-product program ") +
            coverage_program_kind_name(source_program.kind) +
            " resolves bounds per lane");
  }
  switch (source_program.kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
  case CompiledMathSourceProductProgramKind::ConstantOne:
    return;
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    if (coverage_leaf_math_is_scalar(
            source_program.leaf_dist_kind,
            1U | 2U | 4U)) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          root_id,
          "BatchGroupedButScalarLeafMath",
          std::string("source-product leaf program is not vector-native for ") +
              coverage_dist_kind_name(source_program.leaf_dist_kind));
    }
    return;
  case CompiledMathSourceProductProgramKind::ExactGate:
  case CompiledMathSourceProductProgramKind::Conditioned:
    scan_source_product_program_coverage(
        program,
        source_program.child_program_id,
        variant_index,
        probability_program_index,
        root_id,
        depth + 1U,
        rows);
    return;
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        root_id,
        "Unsupported",
        std::string("source-product program kind is not lane-native: ") +
            coverage_program_kind_name(source_program.kind));
    return;
  }
}

void scan_integral_kernel_coverage(const CompiledMathProgram &program,
                                   const CompiledMathIntegralKernel &kernel,
                                   const int variant_index,
                                   const int probability_program_index,
                                   const int root_id,
                                   const std::size_t depth,
                                   BatchCoverageRows *rows) {
  const bool lane_native =
      accumulatr::eval::detail::compiled_math_batch_kernel_supported(
          program,
          kernel,
          0U);
  const bool direct_leaf =
      accumulatr::eval::detail::batch_finite_integral_kernel_direct_leaf_supported(
          program,
          kernel);
  if (!lane_native && !direct_leaf) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        root_id,
        "Unsupported",
        "compiled integral kernel would fall back to scalar evaluation");
  }
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
    scan_source_product_ops_coverage(
        program,
        kernel.source_product_ops,
        variant_index,
        probability_program_index,
        root_id,
        depth + 1U,
        rows);
    return;
  }
  if (kernel.kind != CompiledMathIntegralKernelKind::SourceProductSum) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        root_id,
        "Unsupported",
        "unknown compiled integral kernel kind");
    return;
  }
  for (Index term_idx = 0; term_idx < kernel.source_product_terms.size;
       ++term_idx) {
    const auto term_index =
        static_cast<std::size_t>(kernel.source_product_terms.offset + term_idx);
    if (term_index >= program.integral_kernel_source_product_terms.size()) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          root_id,
          "Unsupported",
          "source-product term span points outside the compiled term table");
      continue;
    }
    const auto &term =
        program.integral_kernel_source_product_terms[term_index];
    if (!term.expr_upper_factors.empty()) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          root_id,
          "BatchGroupedButScalarExprUpper",
          "source-product term applies expr-upper factors per lane");
    }
    scan_source_product_ops_coverage(
        program,
        term.source_product_ops,
        variant_index,
        probability_program_index,
        root_id,
        depth + 1U,
        rows);
    for (Index i = 0; i < term.integral_factor_nodes.size; ++i) {
      const auto factor_index =
          static_cast<std::size_t>(term.integral_factor_nodes.offset + i);
      if (factor_index >= program.integral_kernel_integral_factor_nodes.size()) {
        rows->add(
            "program",
            variant_index,
            probability_program_index,
            root_id,
            "Unsupported",
            "integral-factor span points outside the compiled factor table");
        continue;
      }
      const auto node_id =
          program.integral_kernel_integral_factor_nodes[factor_index];
      const auto node_index = static_cast<std::size_t>(node_id);
      if (node_id == accumulatr::semantic::kInvalidIndex ||
          node_index >= program.nodes.size()) {
        rows->add(
            "program",
            variant_index,
            probability_program_index,
            root_id,
            "Unsupported",
            "integral-factor node points outside the compiled node table");
        continue;
      }
      const auto &node = program.nodes[node_index];
      if (!accumulatr::eval::detail::compiled_math_is_integral_node(
              node.kind) ||
          node.integral_kernel_slot == accumulatr::semantic::kInvalidIndex ||
          static_cast<std::size_t>(node.integral_kernel_slot) >=
              program.integral_kernels.size()) {
        rows->add(
            "program",
            variant_index,
            probability_program_index,
            root_id,
            "Unsupported",
            "integral-factor node is not a compiled integral node");
        continue;
      }
      scan_integral_kernel_coverage(
          program,
          program.integral_kernels[
              static_cast<std::size_t>(node.integral_kernel_slot)],
          variant_index,
          probability_program_index,
          root_id,
          depth + 1U,
          rows);
    }
  }
}

void scan_compiled_root_coverage(const ExactVariantPlan &variant,
                                 const Index root_id,
                                 const int variant_index,
                                 const int probability_program_index,
                                 BatchCoverageRows *rows) {
  if (root_id == accumulatr::semantic::kInvalidIndex ||
      static_cast<std::size_t>(root_id) >= variant.compiled_math.roots.size()) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        coverage_index(root_id),
        "Unsupported",
        "probability program points to an invalid compiled math root");
    return;
  }
  const auto &program = variant.compiled_math;
  const auto &root = program.roots[static_cast<std::size_t>(root_id)];
  for (Index schedule_idx = 0; schedule_idx < root.schedule.size;
       ++schedule_idx) {
    const auto schedule_index =
        static_cast<std::size_t>(root.schedule.offset + schedule_idx);
    if (schedule_index >= program.root_schedule_nodes.size()) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          static_cast<int>(root_id),
          "Unsupported",
          "compiled root schedule points outside the schedule-node table");
      continue;
    }
    const auto node_id = program.root_schedule_nodes[schedule_index];
    const auto node_index = static_cast<std::size_t>(node_id);
    if (node_id == accumulatr::semantic::kInvalidIndex ||
        node_index >= program.nodes.size()) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          static_cast<int>(root_id),
          "Unsupported",
          "compiled root schedule points outside the node table");
      continue;
    }
    const CompiledMathNode &node = program.nodes[node_index];
    if (accumulatr::eval::detail::compiled_math_is_integral_node(node.kind)) {
      if (node.integral_kernel_slot == accumulatr::semantic::kInvalidIndex ||
          static_cast<std::size_t>(node.integral_kernel_slot) >=
              program.integral_kernels.size()) {
        rows->add(
            "program",
            variant_index,
            probability_program_index,
            static_cast<int>(root_id),
            "Unsupported",
            "compiled integral node points outside the kernel table");
        continue;
      }
      scan_integral_kernel_coverage(
          program,
          program.integral_kernels[
              static_cast<std::size_t>(node.integral_kernel_slot)],
          variant_index,
          probability_program_index,
          static_cast<int>(root_id),
          0U,
          rows);
      continue;
    }
    switch (node.kind) {
    case CompiledMathNodeKind::SourcePdf:
    case CompiledMathNodeKind::SourceCdf:
    case CompiledMathNodeKind::SourceSurvival:
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          static_cast<int>(root_id),
          "BatchGroupedButScalarLeafMath",
          std::string("compiled source node is still evaluated per lane: ") +
              coverage_node_kind_name(node.kind));
      break;
    case CompiledMathNodeKind::ExprUpperBoundDensity:
    case CompiledMathNodeKind::ExprUpperBoundCdf:
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          static_cast<int>(root_id),
          "BatchGroupedButScalarExprUpper",
          std::string("compiled expr-upper node is still evaluated per lane: ") +
              coverage_node_kind_name(node.kind));
      break;
    case CompiledMathNodeKind::ExprDensity:
    case CompiledMathNodeKind::ExprCdf:
    case CompiledMathNodeKind::ExprSurvival:
    case CompiledMathNodeKind::UnionKernelMultiSubsetCdf:
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          static_cast<int>(root_id),
          "Unsupported",
          std::string("compiled batch node path rejects node kind: ") +
              coverage_node_kind_name(node.kind));
      break;
    default:
      break;
    }
  }
}

void scan_probability_op_coverage(const ExactVariantPlan &variant,
                                  const ExactProbabilityProgram &program,
                                  const ExactProbabilityOp &op,
                                  const int variant_index,
                                  const int probability_program_index,
                                  const std::size_t depth,
                                  BatchCoverageRows *rows) {
  if (depth > 32U) {
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        r_na_int(),
        "Unsupported",
        "probability op recursion exceeded the coverage depth limit");
    return;
  }
  switch (op.kind) {
  case ExactProbabilityOpKind::Constant:
    return;
  case ExactProbabilityOpKind::Top1LeafRaceDensity:
  case ExactProbabilityOpKind::TerminalNoResponseProbability:
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        r_na_int(),
        "BatchGroupedButScalarLeafMath",
        std::string("probability op uses scalar simple-race leaf math: ") +
            coverage_probability_op_kind_name(op.kind));
    return;
  case ExactProbabilityOpKind::GenericTransitionDensity:
  case ExactProbabilityOpKind::GenericTransitionProbability: {
    if (op.outcome_index == accumulatr::semantic::kInvalidIndex ||
        static_cast<std::size_t>(op.outcome_index) >=
            variant.runtime.outcomes.size()) {
      rows->add(
          "program",
          variant_index,
          probability_program_index,
          r_na_int(),
          "Unsupported",
          "generic transition op points outside runtime outcomes");
      return;
    }
    const auto root_id =
        variant.runtime.outcomes[static_cast<std::size_t>(op.outcome_index)]
            .successor_distribution.total_probability_root_id;
    scan_compiled_root_coverage(
        variant,
        root_id,
        variant_index,
        probability_program_index,
        rows);
    return;
  }
  case ExactProbabilityOpKind::RankedTransitionSequence:
    rows->add(
        "program",
        variant_index,
        probability_program_index,
        r_na_int(),
        "Unsupported",
        "ranked transition sequence is outside the current batch contract");
    return;
  case ExactProbabilityOpKind::Integral:
  case ExactProbabilityOpKind::WeightedTriggerSum:
  case ExactProbabilityOpKind::Log:
    for (Index i = 0; i < op.children.size; ++i) {
      const auto child_index =
          static_cast<std::size_t>(op.children.offset + i);
      if (child_index >= program.child_ops.size()) {
        rows->add(
            "program",
            variant_index,
            probability_program_index,
            r_na_int(),
            "Unsupported",
            "probability op child span points outside child table");
        continue;
      }
      const auto child_op_id = program.child_ops[child_index];
      if (child_op_id == accumulatr::semantic::kInvalidIndex ||
          static_cast<std::size_t>(child_op_id) >= program.ops.size()) {
        rows->add(
            "program",
            variant_index,
            probability_program_index,
            r_na_int(),
            "Unsupported",
            "probability op child points outside op table");
        continue;
      }
      scan_probability_op_coverage(
          variant,
          program,
          program.ops[static_cast<std::size_t>(child_op_id)],
          variant_index,
          probability_program_index,
          depth + 1U,
          rows);
    }
    return;
  }
}

Rcpp::List make_batch_coverage_report(const NativeLikelihoodContext &ctx) {
  BatchCoverageRows rows;

  for (std::size_t variant_pos = 0; variant_pos < ctx.exact_plans.size();
       ++variant_pos) {
    const auto &variant = ctx.exact_plans[variant_pos];
    for (std::size_t program_pos = 0;
         program_pos < variant.probability_programs.programs.size();
         ++program_pos) {
      const auto before = rows.size();
      const auto &program =
          variant.probability_programs.programs[program_pos];
      if (program.empty() ||
          static_cast<std::size_t>(program.root) >= program.ops.size()) {
        rows.add(
            "program",
            static_cast<int>(variant_pos),
            static_cast<int>(program_pos),
            r_na_int(),
            "Unsupported",
            "probability program is empty or has an invalid root op");
      } else {
        scan_probability_op_coverage(
            variant,
            program,
            program.ops[static_cast<std::size_t>(program.root)],
            static_cast<int>(variant_pos),
            static_cast<int>(program_pos),
            0U,
            &rows);
      }
      if (rows.size() == before) {
        rows.add(
            "program",
            static_cast<int>(variant_pos),
            static_cast<int>(program_pos),
            r_na_int(),
            "BatchComplete",
            "no unsupported, scalar leaf-math, scalar-bound, or scalar expr-upper marker found");
      }
    }
  }

  std::array<int, kBatchCoverageCategories.size()> counts{};
  for (const auto &category : rows.category) {
    for (std::size_t i = 0; i < kBatchCoverageCategories.size(); ++i) {
      if (category == kBatchCoverageCategories[i]) {
        ++counts[i];
        break;
      }
    }
  }

  std::vector<std::string> summary_category;
  std::vector<int> summary_count;
  summary_category.reserve(kBatchCoverageCategories.size());
  summary_count.reserve(kBatchCoverageCategories.size());
  for (std::size_t i = 0; i < kBatchCoverageCategories.size(); ++i) {
    summary_category.emplace_back(kBatchCoverageCategories[i]);
    summary_count.push_back(counts[i]);
  }

  return Rcpp::List::create(
      Rcpp::Named("summary") =
          Rcpp::DataFrame::create(
              Rcpp::Named("category") = summary_category,
              Rcpp::Named("count") = summary_count,
              Rcpp::Named("stringsAsFactors") = false),
      Rcpp::Named("programs") =
          Rcpp::DataFrame::create(
              Rcpp::Named("scope") = rows.scope,
              Rcpp::Named("variant_index") = rows.variant_index,
              Rcpp::Named("program_index") = rows.program_index,
              Rcpp::Named("root_id") = rows.root_id,
              Rcpp::Named("category") = rows.category,
              Rcpp::Named("reason") = rows.reason,
              Rcpp::Named("stringsAsFactors") = false));
}

} // namespace

// [[Rcpp::export]]
SEXP semantic_make_likelihood_context_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  auto ctx = accumulatr::eval::detail::build_native_likelihood_context(prep);
  auto ptr = Rcpp::XPtr<accumulatr::eval::detail::NativeLikelihoodContext>(
      new accumulatr::eval::detail::NativeLikelihoodContext(std::move(ctx)),
      true);
  return Rcpp::List::create(
      Rcpp::Named("native") = ptr,
      Rcpp::Named("observed_identity") = ptr->observed_identity,
      Rcpp::Named("identity_backend") = "exact",
      Rcpp::Named("ranked_supported") = ptr->ranked_supported,
      Rcpp::Named("batch_coverage") = make_batch_coverage_report(*ptr));
}

// [[Rcpp::export]]
SEXP semantic_loglik_context_cpp(SEXP contextSEXP,
                                 SEXP paramsSEXP,
                                 SEXP dataSEXP,
                                 SEXP okSEXP,
                                 SEXP expandSEXP,
                                 SEXP minLLSEXP) {
  auto &ctx =
      accumulatr::eval::detail::mutable_likelihood_context_from_xptr(contextSEXP);
  if (!ctx.observation_batch_workspace) {
    ctx.observation_batch_workspace.reset(
        new accumulatr::eval::detail::ObservationBatchWorkspace());
  }
  const auto layout =
      accumulatr::eval::detail::read_prepared_trial_layout(dataSEXP);
  const int *ok = Rf_isNull(okSEXP) ? nullptr : LOGICAL(okSEXP);
  const double min_ll = REAL(minLLSEXP)[0];
  Rcpp::List observed = accumulatr::eval::detail::evaluate_observed_trials_cached(
      ctx.observed_plans_by_component_code,
      ctx.observed_identity,
      ctx.model,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      layout,
      paramsSEXP,
      dataSEXP,
      min_ll,
      expandSEXP,
      ctx.observation_batch_workspace.get(),
      ok);
  return observed;
}

extern "C" {

double accumulatr_cpp_loglik_ccallable(SEXP contextSEXP,
                                       SEXP paramsSEXP,
                                       SEXP dataSEXP,
                                       SEXP okSEXP,
                                       SEXP expandSEXP,
                                       double min_ll) {
  try {
    Rcpp::List observed = semantic_loglik_context_cpp(
        contextSEXP,
        paramsSEXP,
        dataSEXP,
        okSEXP,
        expandSEXP,
        Rcpp::wrap(min_ll));
    return Rcpp::as<double>(observed["total_loglik"]);
  } catch (const std::exception &e) {
    ::Rf_error("%s", e.what());
  } catch (...) {
    ::Rf_error("Unknown C++ exception in AccumulatR::cpp_loglik");
  }
  return NA_REAL;
}

} // extern "C"

// [[Rcpp::init]]
void accumulatr_register_ccallables(DllInfo *dll) {
  (void)dll;
  R_RegisterCCallable(
      "AccumulatR",
      "cpp_loglik",
      reinterpret_cast<DL_FUNC>(accumulatr_cpp_loglik_ccallable));
}

// [[Rcpp::export]]
SEXP semantic_probability_context_cpp(SEXP contextSEXP,
                                      SEXP paramsSEXP,
                                      SEXP dataSEXP) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  const auto layout =
      accumulatr::eval::detail::read_prepared_trial_layout(dataSEXP);
  return accumulatr::eval::detail::evaluate_outcome_queries_cached(
      ctx.observed_plans_by_component_code,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      layout,
      paramsSEXP,
      dataSEXP);
}
