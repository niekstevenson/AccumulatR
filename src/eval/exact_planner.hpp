#pragma once

#include "exact_types.hpp"
#include "exact_compiled_math_lowering.hpp"
#include "exact_transition_lowering.hpp"
#include "exact_runtime_lowering.hpp"
#include "exact_kernel_planning.hpp"
#include "exact_compiled_math_finalize.hpp"

namespace accumulatr::eval {
namespace detail {

inline ExactVariantPlan make_exact_variant_plan(
    const runtime::LoweredExactVariant &lowered,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_outcome_codes) {
  ExactVariantPlan plan;
  plan.lowered = lowered;
  compile_exact_support_context(&plan);
  compile_program_source_runtime_fields(&plan);
  compile_source_kernels(&plan);
  compile_exact_expr_kernels(&plan);
  compile_shared_trigger_state_table(&plan);
  compile_exact_outcome_transition_scenarios(
      &plan,
      outcome_code_by_label,
      n_outcome_codes);
  auto competitor_plans = compile_target_competitor_plans(&plan);
  compile_sequence_expr_upper_bound_roots(&plan, competitor_plans);
  plan.runtime = compile_exact_runtime_plan(&plan, competitor_plans);
  compile_source_arithmetic_dependencies(&plan);
  compile_source_view_relation_tables(&plan);
  compile_source_product_channel_programs(&plan);
  compile_source_product_scalar_ops(&plan);
  compile_source_node_programs(&plan);
  validate_source_product_relations_materialized(plan.compiled_math);
  validate_compiled_math_has_no_interpreter_expr_nodes(plan.compiled_math);
  exact_complexity_finalize(&plan);
  compiled_math_release_planning_fields(&plan.compiled_math);

  return plan;
}

} // namespace detail
} // namespace accumulatr::eval
