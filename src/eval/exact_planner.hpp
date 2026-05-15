#pragma once

#include <utility>

#include "exact_types.hpp"
#include "exact_compiled_math_lowering.hpp"
#include "exact_transition_lowering.hpp"
#include "exact_outcome_planning.hpp"
#include "exact_expr_canonicalization.hpp"
#include "exact_kernel_planning.hpp"
#include "exact_compiled_math_finalize.hpp"

namespace accumulatr::eval {
namespace detail {

inline ExactVariantPlan make_exact_variant_plan(
    runtime::ExactEvaluationProgram program,
    const std::size_t n_outcome_codes,
    ExactComplexityMetrics *complexity_metrics = nullptr) {
  ExactVariantBuildState build;
  build.program = std::move(program);
  build.complexity_metrics = complexity_metrics;
  canonicalize_exact_evaluation_program_expressions(&build.program);
  compile_exact_support_context(&build);
  compile_program_source_runtime_fields(&build);
  compile_source_kernels(&build);
  compile_exact_expr_kernels(&build);
  compile_shared_trigger_state_table(&build);
  compile_exact_outcome_transition_scenarios(&build, n_outcome_codes);
  build.no_response = compile_terminal_no_response_plan(build);
  auto competitor_plans = compile_target_competitor_plans(&build);
  compile_sequence_plan(&build, competitor_plans);
  build.compiled_outcomes =
      compile_exact_outcome_plans(&build, competitor_plans);
  compile_source_arithmetic_dependencies(&build);
  compile_source_view_relation_tables(&build);
  compile_source_product_channel_programs(&build);
  compile_source_product_scalar_ops(&build);
  compile_source_node_programs(&build);
  validate_source_product_relations_materialized(build.compiled_math);
  validate_compiled_math_has_no_interpreter_expr_nodes(build.compiled_math);
  exact_complexity_finalize(&build);
  compiled_math_release_planning_fields(&build.compiled_math);

  return finalize_exact_variant_plan(std::move(build));
}

} // namespace detail
} // namespace accumulatr::eval
