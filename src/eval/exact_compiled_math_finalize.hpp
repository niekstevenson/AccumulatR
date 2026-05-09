#pragma once

#include "exact_types.hpp"
#include "exact_compiled_math_lowering.hpp"

namespace accumulatr::eval {
namespace detail {

inline void validate_compiled_math_has_no_interpreter_expr_nodes(
    const CompiledMathProgram &program) {
  for (std::size_t i = 0; i < program.nodes.size(); ++i) {
    const auto kind = program.nodes[i].kind;
    if (kind == CompiledMathNodeKind::ExprDensity ||
        kind == CompiledMathNodeKind::ExprCdf ||
        kind == CompiledMathNodeKind::ExprSurvival) {
      throw std::runtime_error(
          "exact compiled math contains an unlowered semantic node at " +
          std::to_string(i));
    }
  }
}

inline semantic::Index exact_complexity_integral_depth_for_node(
    const CompiledMathProgram &program,
    const semantic::Index node_id,
    std::vector<semantic::Index> *memo,
    std::vector<std::uint8_t> *visiting) {
  if (node_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(node_id) >= program.nodes.size()) {
    return 0;
  }
  const auto node_pos = static_cast<std::size_t>(node_id);
  if ((*memo)[node_pos] != semantic::kInvalidIndex) {
    return (*memo)[node_pos];
  }
  if ((*visiting)[node_pos] != 0U) {
    return 0;
  }
  (*visiting)[node_pos] = 1U;
  const auto &node = program.nodes[node_pos];
  semantic::Index depth{0};
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_pos =
        static_cast<std::size_t>(node.children.offset + i);
    if (child_pos >= program.child_nodes.size()) {
      continue;
    }
    depth = std::max(
        depth,
        exact_complexity_integral_depth_for_node(
            program,
            program.child_nodes[child_pos],
            memo,
            visiting));
  }
  if (compiled_math_is_integral_node(node.kind)) {
    const auto root_id = compiled_math_integral_root_id(node);
    semantic::Index child_depth{0};
    if (root_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(root_id) < program.roots.size()) {
      child_depth =
          exact_complexity_integral_depth_for_node(
              program,
              program.roots[static_cast<std::size_t>(root_id)].node_id,
              memo,
              visiting);
    }
    depth = std::max(depth, static_cast<semantic::Index>(child_depth + 1));
  }
  (*visiting)[node_pos] = 0U;
  (*memo)[node_pos] = depth;
  return depth;
}

inline semantic::Index exact_complexity_max_integral_depth(
    const CompiledMathProgram &program) {
  std::vector<semantic::Index> memo(
      program.nodes.size(), semantic::kInvalidIndex);
  std::vector<std::uint8_t> visiting(program.nodes.size(), 0U);
  semantic::Index out{0};
  for (const auto &root : program.roots) {
    out =
        std::max(
            out,
            exact_complexity_integral_depth_for_node(
                program, root.node_id, &memo, &visiting));
  }
  return out;
}

inline void exact_complexity_finalize(ExactVariantPlan *plan) {
  auto &metrics = plan->complexity;
  const auto &program = plan->compiled_math;
  metrics.compiled_root_count =
      static_cast<semantic::Index>(program.roots.size());
  metrics.compiled_node_count =
      static_cast<semantic::Index>(program.nodes.size());
  metrics.integral_node_count = 0;
  for (const auto &node : program.nodes) {
    if (compiled_math_is_integral_node(node.kind)) {
      ++metrics.integral_node_count;
    }
  }
  metrics.integral_kernel_count =
      static_cast<semantic::Index>(program.integral_kernels.size());
  metrics.source_product_integral_kernel_count = 0;
  metrics.generic_integral_kernel_count = 0;
  for (const auto &kernel : program.integral_kernels) {
    if (kernel.kind == CompiledMathIntegralKernelKind::Generic) {
      ++metrics.generic_integral_kernel_count;
    } else {
      ++metrics.source_product_integral_kernel_count;
    }
  }
  metrics.max_integral_depth =
      exact_complexity_max_integral_depth(program);
}

inline semantic::Index compiled_source_bound_plan_slot(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  if (source_id == semantic::kInvalidIndex ||
      program.source_condition_bound_source_count <= 0) {
    return semantic::kInvalidIndex;
  }
  const auto source_count =
      static_cast<std::size_t>(program.source_condition_bound_source_count);
  const auto source_pos = static_cast<std::size_t>(source_id);
  if (source_pos >= source_count) {
    return semantic::kInvalidIndex;
  }
  const auto condition_slot =
      condition_id == semantic::kInvalidIndex ? 0 : condition_id;
  const auto condition_pos = static_cast<std::size_t>(condition_slot);
  if (condition_pos > program.conditions.size()) {
    return semantic::kInvalidIndex;
  }
  const auto slot = condition_pos * source_count + source_pos;
  if (slot >= program.source_condition_bound_plans.size()) {
    return semantic::kInvalidIndex;
  }
  return static_cast<semantic::Index>(slot);
}

inline void compile_source_condition_bound_plans(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  const auto source_count = static_cast<std::size_t>(plan->source_count);
  program.source_condition_bound_source_count = plan->source_count;
  program.condition_source_relation_source_count = plan->source_count;
  program.source_condition_bound_plans.assign(
      (program.conditions.size() + 1U) * source_count,
      CompiledSourceBoundPlan{});
  program.source_condition_bound_terms.clear();
  program.condition_source_relations.assign(
      (program.conditions.size() + 1U) * source_count,
      static_cast<std::uint8_t>(ExactRelation::Unknown));
  if (source_count == 0U) {
    return;
  }
  auto append_bound_terms = [&program](
                                const CompiledMathCondition &condition,
                                const std::vector<semantic::Index> &fact_indices,
                                const CompiledMathIndexSpan fact_span) {
    const auto offset =
        static_cast<semantic::Index>(
            program.source_condition_bound_terms.size());
    for (semantic::Index i = 0; i < fact_span.size; ++i) {
      const auto fact_pos = static_cast<std::size_t>(
          fact_indices[
              static_cast<std::size_t>(fact_span.offset + i)]);
      const auto time_id =
          fact_pos < condition.fact_time_ids.size()
              ? condition.fact_time_ids[fact_pos]
              : static_cast<semantic::Index>(
                    CompiledMathTimeSlot::Observed);
      program.source_condition_bound_terms.push_back(
          CompiledSourceBoundTerm{time_id});
    }
    return CompiledMathIndexSpan{
        offset,
        static_cast<semantic::Index>(
            program.source_condition_bound_terms.size() -
            static_cast<std::size_t>(offset))};
  };
  for (std::size_t condition_pos = 0;
       condition_pos < program.conditions.size();
       ++condition_pos) {
    const auto &condition = program.conditions[condition_pos];
    const auto condition_id =
        static_cast<semantic::Index>(condition_pos + 1U);
    const auto relation_offset = (condition_pos + 1U) * source_count;
    for (std::size_t i = 0; i < condition.source_ids.size(); ++i) {
      const auto source_id = condition.source_ids[i];
      if (source_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(source_id) >= source_count) {
        continue;
      }
      program.condition_source_relations[
          relation_offset + static_cast<std::size_t>(source_id)] =
          condition.relations[i];
    }
    for (const auto &entry : condition.source_order_fact_lookup.entries) {
      const auto before_pos = static_cast<std::size_t>(entry.first);
      const auto after_pos = static_cast<std::size_t>(entry.second);
      if (entry.first == semantic::kInvalidIndex ||
          entry.second == semantic::kInvalidIndex ||
          before_pos >= source_count ||
          after_pos >= source_count) {
        continue;
      }
      program.condition_source_relations[relation_offset + before_pos] =
          static_cast<std::uint8_t>(ExactRelation::Before);
      program.condition_source_relations[relation_offset + after_pos] =
          static_cast<std::uint8_t>(ExactRelation::After);
    }
    for (std::size_t source_pos = 0; source_pos < source_count; ++source_pos) {
      const auto slot_id =
          compiled_source_bound_plan_slot(
              program,
              condition_id,
              static_cast<semantic::Index>(source_pos));
      if (slot_id == semantic::kInvalidIndex) {
        continue;
      }
      const auto slot = static_cast<std::size_t>(slot_id);
      auto &bounds = program.source_condition_bound_plans[slot];
      if (source_pos < condition.source_exact_fact_spans.size()) {
        bounds.exact =
            append_bound_terms(
                condition,
                condition.source_exact_fact_indices,
                condition.source_exact_fact_spans[source_pos]);
        bounds.has_condition_exact = !bounds.exact.empty();
      }
      if (source_pos < condition.source_lower_fact_spans.size()) {
        bounds.lower =
            append_bound_terms(
                condition,
                condition.source_lower_fact_indices,
                condition.source_lower_fact_spans[source_pos]);
        bounds.has_condition_lower = !bounds.lower.empty();
      }
      if (source_pos < condition.source_upper_fact_spans.size()) {
        bounds.upper =
            append_bound_terms(
                condition,
                condition.source_upper_fact_indices,
                condition.source_upper_fact_spans[source_pos]);
        bounds.has_condition_upper = !bounds.upper.empty();
      }
    }
  }
}

inline void compile_condition_cache_plans(CompiledMathProgram *program) {
  program->condition_cache_plans.clear();
  program->condition_cache_time_dependencies.clear();
  program->condition_cache_plans.reserve(program->conditions.size() + 1U);
  program->condition_cache_plans.push_back(CompiledConditionCachePlan{});
  for (std::size_t condition_pos = 0;
       condition_pos < program->conditions.size();
       ++condition_pos) {
    const auto condition_id =
        static_cast<semantic::Index>(condition_pos + 1U);
    const auto &condition = program->conditions[condition_pos];
    const auto offset =
        static_cast<semantic::Index>(
            program->condition_cache_time_dependencies.size());
    program->condition_cache_time_dependencies.insert(
        program->condition_cache_time_dependencies.end(),
        condition.time_dependency_ids.begin(),
        condition.time_dependency_ids.end());
    const auto size =
        static_cast<semantic::Index>(
            program->condition_cache_time_dependencies.size() -
            static_cast<std::size_t>(offset));
    program->condition_cache_plans.push_back(
        CompiledConditionCachePlan{
            CompiledMathIndexSpan{offset, size},
            compiled_math_static_condition_cache_id(condition_id),
            size > 0});
  }
}

inline void compile_source_arithmetic_dependencies(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  compile_condition_cache_plans(&program);
  compile_source_condition_bound_plans(plan);
}

inline const CompiledSourceBoundPlan &source_product_bound_plan_for(
    const ExactVariantPlan &plan,
    const semantic::Index condition_id,
    const semantic::Index source_id);

inline void compile_source_product_channel_fields(
    ExactVariantPlan *plan,
    CompiledMathSourceProductChannel *channel) {
  if (channel == nullptr) {
    return;
  }
  if (channel->source_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(channel->source_id) <
          plan->source_kernels.size()) {
    channel->source_kernel_slot = channel->source_id;
    channel->kernel =
        plan->source_kernels[static_cast<std::size_t>(channel->source_id)]
            .kind;
  }
  if (channel->source_id != semantic::kInvalidIndex) {
    channel->static_source_view_relation = static_cast<std::uint8_t>(
        exact_compiled_source_view_relation(
            *plan, channel->source_view_id, channel->source_id));
    channel->has_static_source_view_relation = true;
  }
  if (channel->source_kernel_slot != semantic::kInvalidIndex &&
      static_cast<std::size_t>(channel->source_kernel_slot) <
          plan->source_kernels.size()) {
    const auto &source_kernel =
        plan->source_kernels[
            static_cast<std::size_t>(channel->source_kernel_slot)];
    channel->leaf_index = source_kernel.leaf_index;
    if (channel->leaf_index != semantic::kInvalidIndex &&
        static_cast<std::size_t>(channel->leaf_index) <
            plan->lowered.program.leaf_descriptors.size()) {
      const auto &leaf =
          plan->lowered.program.leaf_descriptors[
              static_cast<std::size_t>(channel->leaf_index)];
      channel->leaf_dist_kind = leaf.dist_kind;
      channel->leaf_param_count = leaf.param_count;
      channel->leaf_onset_abs_value = leaf.onset_abs_value;
    }
  }
  channel->bounds =
      source_product_bound_plan_for(
          *plan, channel->condition_id, channel->source_id);
  channel->has_source_condition_overlay =
      channel->bounds.has_condition_exact ||
      channel->bounds.has_condition_lower ||
      channel->bounds.has_condition_upper;
  channel->direct_leaf_absolute_candidate =
      channel->source_id != semantic::kInvalidIndex &&
      channel->source_kernel_slot != semantic::kInvalidIndex &&
      channel->kernel == CompiledSourceChannelKernelKind::LeafAbsolute &&
      !channel->has_source_condition_overlay;
}

inline void compile_source_product_channel_programs(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  for (auto &channel : program.integral_kernel_source_product_channels) {
    compile_source_product_channel_fields(plan, &channel);
  }
}

inline const CompiledSourceBoundPlan &source_product_bound_plan_for(
    const ExactVariantPlan &plan,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  static const CompiledSourceBoundPlan empty{};
  const auto slot = compiled_source_bound_plan_slot(
      plan.compiled_math, condition_id, source_id);
  if (slot == semantic::kInvalidIndex ||
      static_cast<std::size_t>(slot) >=
          plan.compiled_math.source_condition_bound_plans.size()) {
    return empty;
  }
  return plan.compiled_math.source_condition_bound_plans[
      static_cast<std::size_t>(slot)];
}

inline semantic::Index push_source_product_program(
    CompiledMathProgram *program,
    CompiledMathSourceProductProgram source_program) {
  const auto id = static_cast<semantic::Index>(
      program->integral_kernel_source_product_programs.size());
  program->integral_kernel_source_product_programs.push_back(
      source_program);
  return id;
}

inline semantic::Index compile_source_product_base_program(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index source_view_id);

inline semantic::Index compile_source_product_exact_gate_program(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index source_view_id,
    const semantic::Index child_program_id) {
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::ExactGate;
  source_program.source_id = source_id;
  source_program.condition_id = condition_id;
  source_program.source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  source_program.child_program_id = child_program_id;
  source_program.bounds =
      source_product_bound_plan_for(*plan, condition_id, source_id);
  if (source_id != semantic::kInvalidIndex) {
    source_program.static_source_view_relation = static_cast<std::uint8_t>(
        exact_compiled_source_view_relation(
            *plan, source_program.source_view_id, source_id));
    source_program.has_static_source_view_relation = true;
  }
  return push_source_product_program(&plan->compiled_math, source_program);
}

inline semantic::Index compile_source_product_leaf_program(
    ExactVariantPlan *plan,
    const ExactSourceKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::LeafAbsolute;
  source_program.source_id = kernel.source_id;
  source_program.condition_id = condition_id;
  source_program.source_view_id = source_view_id;
  source_program.leaf_index = kernel.leaf_index;
  if (kernel.leaf_index != semantic::kInvalidIndex &&
      static_cast<std::size_t>(kernel.leaf_index) <
          plan->lowered.program.leaf_descriptors.size()) {
    const auto &leaf =
        plan->lowered.program.leaf_descriptors[
            static_cast<std::size_t>(kernel.leaf_index)];
    source_program.leaf_dist_kind = leaf.dist_kind;
    source_program.leaf_onset_abs_value = leaf.onset_abs_value;
  }
  const auto leaf_program_id =
      push_source_product_program(&plan->compiled_math, source_program);
  return compile_source_product_exact_gate_program(
      plan, kernel.source_id, condition_id, source_view_id, leaf_program_id);
}

inline semantic::Index compile_source_product_onset_program(
    ExactVariantPlan *plan,
    const ExactSourceKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::OnsetConvolution;
  source_program.source_id = kernel.source_id;
  source_program.condition_id = condition_id;
  source_program.source_view_id = source_view_id;
  source_program.leaf_index = kernel.leaf_index;
  source_program.onset_source_program_id =
      compile_source_product_base_program(
          plan, kernel.onset_source_id, condition_id, source_view_id);
  source_program.onset_bounds =
      source_product_bound_plan_for(
          *plan, condition_id, kernel.onset_source_id);
  if (kernel.leaf_index != semantic::kInvalidIndex &&
      static_cast<std::size_t>(kernel.leaf_index) <
          plan->lowered.program.leaf_descriptors.size()) {
    const auto &leaf =
        plan->lowered.program.leaf_descriptors[
            static_cast<std::size_t>(kernel.leaf_index)];
    source_program.leaf_dist_kind = leaf.dist_kind;
    source_program.leaf_onset_lag = leaf.onset_lag;
  }
  const auto onset_program_id =
      push_source_product_program(&plan->compiled_math, source_program);
  return compile_source_product_exact_gate_program(
      plan, kernel.source_id, condition_id, source_view_id, onset_program_id);
}

inline semantic::Index compile_source_product_pool_program(
    ExactVariantPlan *plan,
    const ExactSourceKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  auto &program = plan->compiled_math;
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::PoolKOfN;
  source_program.source_id = kernel.source_id;
  source_program.condition_id = condition_id;
  source_program.source_view_id = source_view_id;
  source_program.pool_k = kernel.pool_k;
  const auto member_offset = static_cast<semantic::Index>(
      program.integral_kernel_source_product_program_members.size());
  const auto member_end = kernel.pool_member_offset + kernel.pool_member_count;
  for (semantic::Index i = kernel.pool_member_offset; i < member_end; ++i) {
    const auto member_source =
        plan->lowered.program.pool_member_source_ids[
            static_cast<std::size_t>(i)];
    program.integral_kernel_source_product_program_members.push_back(
        compile_source_product_base_program(
            plan, member_source, condition_id, source_view_id));
  }
  source_program.member_programs = CompiledMathIndexSpan{
      member_offset,
      kernel.pool_member_count};
  const auto member_count =
      static_cast<semantic::Index>(kernel.pool_member_count);
  const auto width = member_count + 1;
  const auto table_size = width * width;
  source_program.source_product_scratch_offset =
      program.integral_kernel_source_product_scratch_size;
  source_program.source_product_scratch_size =
      3 * member_count + 2 * table_size;
  program.integral_kernel_source_product_scratch_size +=
      source_program.source_product_scratch_size;
  const auto pool_program_id =
      push_source_product_program(&program, source_program);
  return compile_source_product_exact_gate_program(
      plan, kernel.source_id, condition_id, source_view_id, pool_program_id);
}

inline semantic::Index compile_source_product_base_program(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  if (source_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(source_id) >= plan->source_kernels.size()) {
    CompiledMathSourceProductProgram source_program;
    source_program.kind = CompiledMathSourceProductProgramKind::ConstantZero;
    return push_source_product_program(&plan->compiled_math, source_program);
  }
  const auto &kernel =
      plan->source_kernels[static_cast<std::size_t>(source_id)];
  switch (kernel.kind) {
  case CompiledSourceChannelKernelKind::LeafAbsolute:
    return compile_source_product_leaf_program(
        plan, kernel, condition_id, source_view_id);
  case CompiledSourceChannelKernelKind::LeafOnsetConvolution:
    return compile_source_product_onset_program(
        plan, kernel, condition_id, source_view_id);
  case CompiledSourceChannelKernelKind::PoolKOfN:
    return compile_source_product_pool_program(
        plan, kernel, condition_id, source_view_id);
  case CompiledSourceChannelKernelKind::Invalid:
    break;
  }
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::ConstantZero;
  return push_source_product_program(&plan->compiled_math, source_program);
}

inline semantic::Index compile_source_product_channel_program(
    ExactVariantPlan *plan,
    CompiledMathSourceProductChannel *channel) {
  if (channel->source_product_program_id != semantic::kInvalidIndex) {
    return channel->source_product_program_id;
  }
  const auto child_program_id =
      compile_source_product_base_program(
          plan,
          channel->source_id,
          channel->condition_id,
          channel->source_view_id);
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::Conditioned;
  source_program.source_id = channel->source_id;
  source_program.condition_id = channel->condition_id;
  source_program.time_id = channel->time_id;
  source_program.time_cap_id = channel->time_cap_id;
  source_program.source_view_id = channel->source_view_id;
  source_program.child_program_id = child_program_id;
  source_program.bounds = channel->bounds;
  source_program.static_source_view_relation =
      channel->static_source_view_relation;
  source_program.has_static_source_view_relation =
      channel->has_static_source_view_relation;
  channel->source_product_program_id =
      push_source_product_program(&plan->compiled_math, source_program);
  return channel->source_product_program_id;
}

inline CompiledMathSourceProductOpKind source_product_generic_op_kind(
    const CompiledMathNodeKind factor_kind) noexcept {
  switch (factor_kind) {
  case CompiledMathNodeKind::SourcePdf:
    return CompiledMathSourceProductOpKind::GenericPdf;
  case CompiledMathNodeKind::SourceCdf:
    return CompiledMathSourceProductOpKind::GenericCdf;
  case CompiledMathNodeKind::SourceSurvival:
    return CompiledMathSourceProductOpKind::GenericSurvival;
  default:
    break;
  }
  return CompiledMathSourceProductOpKind::ConstantZero;
}

inline CompiledMathSourceProductOpKind source_product_forced_relation_op_kind(
    const ExactRelation relation,
    const CompiledMathNodeKind factor_kind) noexcept {
  if (relation == ExactRelation::Before) {
    return factor_kind == CompiledMathNodeKind::SourceCdf
               ? CompiledMathSourceProductOpKind::ConstantOne
               : CompiledMathSourceProductOpKind::ConstantZero;
  }
  if (relation == ExactRelation::At) {
    if (factor_kind == CompiledMathNodeKind::SourcePdf) {
      return source_product_generic_op_kind(factor_kind);
    }
    return factor_kind == CompiledMathNodeKind::SourceCdf
               ? CompiledMathSourceProductOpKind::ConstantOne
               : CompiledMathSourceProductOpKind::ConstantZero;
  }
  if (relation == ExactRelation::After) {
    return factor_kind == CompiledMathNodeKind::SourceSurvival
               ? CompiledMathSourceProductOpKind::ConstantOne
               : CompiledMathSourceProductOpKind::ConstantZero;
  }
  return source_product_generic_op_kind(factor_kind);
}

inline CompiledMathSourceProductOpKind source_product_leaf_op_kind(
    const std::uint8_t leaf_dist_kind,
    const CompiledMathNodeKind factor_kind) noexcept {
  const auto dist_kind = static_cast<leaf::DistKind>(leaf_dist_kind);
  switch (dist_kind) {
  case leaf::DistKind::Lognormal:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafLognormalPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafLognormalCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafLognormalSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::Gamma:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafGammaPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafGammaCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafGammaSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::Exgauss:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafExgaussPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafExgaussCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafExgaussSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::LBA:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafLbaPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafLbaCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafLbaSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::RDM:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafRdmPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafRdmCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafRdmSurvival;
    default:
      break;
    }
    break;
  }
  return source_product_generic_op_kind(factor_kind);
}

inline CompiledMathSourceProductOpKind source_product_op_kind_for_factor(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathNodeKind factor_kind) noexcept {
  if (channel.has_static_source_view_relation &&
      !channel.has_source_condition_overlay) {
    const auto relation =
        static_cast<ExactRelation>(channel.static_source_view_relation);
    if (relation != ExactRelation::Unknown) {
      if (relation != ExactRelation::At ||
          factor_kind != CompiledMathNodeKind::SourcePdf) {
        return source_product_forced_relation_op_kind(relation, factor_kind);
      }
    }
  }
  if (channel.direct_leaf_absolute_candidate) {
    return source_product_leaf_op_kind(channel.leaf_dist_kind, factor_kind);
  }
  return source_product_generic_op_kind(factor_kind);
}

inline std::uint8_t source_product_op_fill_mask(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathSourceProductOpKind op_kind) noexcept {
  const auto requested =
      compiled_math_source_product_op_channel_mask(op_kind);
  if ((requested & (2U | 4U)) != 0U) {
    const auto paired =
        channel.required_channels & (2U | 4U);
    return paired == 0U ? requested : paired;
  }
  return requested;
}

inline CompiledMathIndexSpan compile_source_product_ops_for_factor_span(
    ExactVariantPlan *plan,
    const CompiledMathIndexSpan factors) {
  auto *program = &plan->compiled_math;
  const auto offset = static_cast<semantic::Index>(
      program->integral_kernel_source_product_ops.size());
  for (semantic::Index i = 0; i < factors.size; ++i) {
    const auto factor_pos =
        static_cast<std::size_t>(factors.offset + i);
    const auto &factor =
        program->integral_kernel_source_value_factors[factor_pos];
    const auto channel_pos =
        static_cast<std::size_t>(factor.source_product_channel_id);
    const auto &channel =
        program->integral_kernel_source_product_channels[channel_pos];
    const auto kind = source_product_op_kind_for_factor(channel, factor.kind);
    if (kind == CompiledMathSourceProductOpKind::ConstantZero) {
      program->integral_kernel_source_product_ops.resize(
          static_cast<std::size_t>(offset));
      program->integral_kernel_source_product_ops.push_back(
          CompiledMathSourceProductOp{
              kind,
              semantic::kInvalidIndex,
              semantic::kInvalidIndex,
              0U,
              0U,
              0.0});
      return CompiledMathIndexSpan{offset, 1};
    }
    if (kind == CompiledMathSourceProductOpKind::ConstantOne) {
      continue;
    }
    const auto value_mask =
        compiled_math_source_product_op_channel_mask(kind);
    const auto fill_mask =
        static_cast<std::uint8_t>(
            value_mask == 0U
                ? 0U
                : source_product_op_fill_mask(channel, kind));
    program->integral_kernel_source_product_channels[channel_pos]
        .scalar_op_count++;
    const auto program_id =
        compile_source_product_channel_program(
            plan,
            &program->integral_kernel_source_product_channels[channel_pos]);
    program->integral_kernel_source_product_ops.push_back(
        CompiledMathSourceProductOp{
            kind,
            factor.source_product_channel_id,
            program_id,
            value_mask,
            fill_mask,
            kind == CompiledMathSourceProductOpKind::ConstantOne ? 1.0 : 0.0});
  }
  return CompiledMathIndexSpan{
      offset,
      static_cast<semantic::Index>(
          program->integral_kernel_source_product_ops.size() -
          static_cast<std::size_t>(offset))};
}

inline void compile_source_product_scalar_ops(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  program.integral_kernel_source_product_ops.clear();
  program.integral_kernel_source_product_programs.clear();
  program.integral_kernel_source_product_program_members.clear();
  program.integral_kernel_source_product_scratch_size = 0;
  for (auto &channel : program.integral_kernel_source_product_channels) {
    channel.scalar_op_count = 0;
    channel.source_product_program_id = semantic::kInvalidIndex;
  }
  for (auto &kernel : program.integral_kernels) {
    if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
      kernel.source_product_ops =
          compile_source_product_ops_for_factor_span(
              plan, kernel.source_value_factors);
      continue;
    }
    for (semantic::Index i = 0; i < kernel.source_product_terms.size; ++i) {
      auto &term =
          program.integral_kernel_source_product_terms[
              static_cast<std::size_t>(
                  kernel.source_product_terms.offset + i)];
      term.source_product_ops =
          compile_source_product_ops_for_factor_span(
              plan, term.source_value_factors);
    }
  }
  for (auto &op : program.integral_kernel_source_product_ops) {
    if (op.value_channel_mask == 0U) {
      op.cache_result = false;
      continue;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];
    op.cache_result = op.source_product_program_id != semantic::kInvalidIndex &&
                      (channel.scalar_op_count > 1 ||
                      op.fill_channel_mask != op.value_channel_mask);
  }
}

struct SourceArithmeticNodeProgramKey {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index source_view_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};

  bool operator==(const SourceArithmeticNodeProgramKey &other) const noexcept {
    return source_id == other.source_id &&
           condition_id == other.condition_id &&
           source_view_id == other.source_view_id &&
           time_id == other.time_id &&
           time_cap_id == other.time_cap_id;
  }
};

struct SourceArithmeticNodeProgramKeyHash {
  std::size_t operator()(const SourceArithmeticNodeProgramKey &key) const
      noexcept {
    std::size_t seed = static_cast<std::size_t>(key.source_id);
    hash_combine(&seed, static_cast<std::size_t>(key.condition_id));
    hash_combine(&seed, static_cast<std::size_t>(key.source_view_id));
    hash_combine(&seed, static_cast<std::size_t>(key.time_id));
    hash_combine(&seed, static_cast<std::size_t>(key.time_cap_id));
    return seed;
  }

private:
  static void hash_combine(std::size_t *seed, const std::size_t value) noexcept {
    *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
  }
};

inline semantic::Index compile_source_node_program(
    ExactVariantPlan *plan,
    const CompiledMathNode &node) {
  CompiledMathSourceProductChannel channel;
  channel.source_id = node.subject_id;
  channel.condition_id = node.condition_id;
  channel.source_view_id =
      node.source_view_id == semantic::kInvalidIndex ? 0 : node.source_view_id;
  channel.time_id = node.time_id;
  channel.time_cap_id = node.aux_id;
  channel.required_channels =
      compiled_math_source_factor_channel_mask(node.kind);
  compile_source_product_channel_fields(plan, &channel);
  return compile_source_product_channel_program(plan, &channel);
}

inline void compile_source_node_programs(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  std::unordered_map<
      SourceArithmeticNodeProgramKey,
      semantic::Index,
      SourceArithmeticNodeProgramKeyHash>
      program_index;
  for (auto &node : program.nodes) {
    if (!compiled_math_is_source_value_node(node.kind)) {
      continue;
    }
    SourceArithmeticNodeProgramKey key;
    key.source_id = node.subject_id;
    key.condition_id = node.condition_id;
    key.source_view_id =
        node.source_view_id == semantic::kInvalidIndex ? 0 : node.source_view_id;
    key.time_id = node.time_id;
    key.time_cap_id = node.aux_id;
    const auto found = program_index.find(key);
    if (found != program_index.end()) {
      node.source_program_id = found->second;
      continue;
    }
    const auto program_id = compile_source_node_program(plan, node);
    program_index.emplace(key, program_id);
    node.source_program_id = program_id;
  }
}

inline void validate_source_product_relations_materialized(
    const CompiledMathProgram &program) {
  for (std::size_t i = 0;
       i < program.integral_kernel_source_product_channels.size();
       ++i) {
    const auto &channel = program.integral_kernel_source_product_channels[i];
    if (channel.source_id != semantic::kInvalidIndex &&
        !channel.has_static_source_view_relation) {
      throw std::runtime_error(
          "source-product channel " + std::to_string(i) +
          " has no compiled source-view relation");
    }
  }
  for (std::size_t i = 0;
       i < program.integral_kernel_source_product_programs.size();
       ++i) {
    const auto &source_program =
        program.integral_kernel_source_product_programs[i];
    const bool relation_sensitive =
        source_program.kind ==
            CompiledMathSourceProductProgramKind::Conditioned ||
        source_program.kind ==
            CompiledMathSourceProductProgramKind::ExactGate;
    if (relation_sensitive &&
        source_program.source_id != semantic::kInvalidIndex &&
        !source_program.has_static_source_view_relation) {
      throw std::runtime_error(
          "source-product program " + std::to_string(i) +
          " has no compiled source-view relation");
    }
  }
  for (std::size_t i = 0;
       i < program.integral_kernel_source_product_ops.size();
       ++i) {
    const auto &op = program.integral_kernel_source_product_ops[i];
    if (op.value_channel_mask != 0U &&
        (op.source_product_program_id == semantic::kInvalidIndex ||
         static_cast<std::size_t>(op.source_product_program_id) >=
             program.integral_kernel_source_product_programs.size())) {
      throw std::runtime_error(
          "source-product op " + std::to_string(i) +
          " has no compiled source-product program");
    }
  }
  for (std::size_t i = 0; i < program.nodes.size(); ++i) {
    const auto &node = program.nodes[i];
    if (!compiled_math_is_source_value_node(node.kind)) {
      continue;
    }
    if (node.source_program_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(node.source_program_id) >=
            program.integral_kernel_source_product_programs.size()) {
      throw std::runtime_error(
          "source node " + std::to_string(i) +
          " has no compiled source arithmetic program");
    }
  }
}

inline void compile_source_view_relation_tables(ExactVariantPlan *plan) {
  const auto source_count = static_cast<std::size_t>(plan->source_count);
  plan->compiled_source_view_source_count = plan->source_count;
  plan->compiled_source_view_relations.assign(
      plan->compiled_source_views.size() * source_count,
      static_cast<std::uint8_t>(ExactRelation::Unknown));
  if (source_count == 0U) {
    return;
  }
  for (std::size_t view_pos = 0;
       view_pos < plan->compiled_source_views.size();
       ++view_pos) {
    const auto &view = plan->compiled_source_views[view_pos];
    const auto view_offset = view_pos * source_count;
    for (std::size_t i = 0; i < view.source_ids.size(); ++i) {
      const auto source_id = view.source_ids[i];
      if (source_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(source_id) >= source_count) {
        continue;
      }
      plan->compiled_source_view_relations[
          view_offset + static_cast<std::size_t>(source_id)] =
          static_cast<std::uint8_t>(view.relations[i]);
    }
  }
}
} // namespace detail
} // namespace accumulatr::eval
