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

inline void exact_complexity_finalize(ExactVariantBuildState *plan) {
  if (plan->complexity_metrics == nullptr) {
    return;
  }
  auto &metrics = *plan->complexity_metrics;
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
    if (kernel.execution.kind == CompiledMathExecutionKind::Schedule) {
      ++metrics.generic_integral_kernel_count;
    } else {
      ++metrics.source_product_integral_kernel_count;
    }
  }
  metrics.max_integral_depth =
      exact_complexity_max_integral_depth(program);
}

inline void compile_source_product_channel_fields(
    ExactVariantBuildState *plan,
    CompiledMathSourceProductChannel *channel) {
  if (channel->source_id != semantic::kInvalidIndex) {
    channel->static_source_view_relation = static_cast<std::uint8_t>(
        exact_compiled_source_view_relation(
            *plan, channel->source_view_id, channel->source_id));
    channel->has_static_source_view_relation = true;
  }
}

inline void compile_source_product_channel_programs(ExactVariantBuildState *plan) {
  auto &program = plan->compiled_math;
  for (auto &channel : program.source_product_channels) {
    compile_source_product_channel_fields(plan, &channel);
  }
}

inline semantic::Index push_source_product_program(
    CompiledMathProgram *program,
    CompiledMathSourceProductProgram source_program) {
  const auto id = static_cast<semantic::Index>(
      program->source_programs.size());
  program->source_programs.push_back(
      source_program);
  return id;
}

inline semantic::Index compile_source_product_base_program(
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index source_view_id);

inline semantic::Index compile_source_product_exact_gate_program(
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const semantic::Index source_view_id,
    const semantic::Index child_program_id) {
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::ExactGate;
  source_program.source_id = source_id;
  source_program.child_program_id = child_program_id;
  if (source_id != semantic::kInvalidIndex) {
    source_program.static_source_view_relation = static_cast<std::uint8_t>(
        exact_compiled_source_view_relation(
            *plan,
            source_view_id == semantic::kInvalidIndex ? 0 : source_view_id,
            source_id));
    source_program.has_static_source_view_relation = true;
  }
  return push_source_product_program(&plan->compiled_math, source_program);
}

inline semantic::Index compile_source_product_leaf_program(
    ExactVariantBuildState *plan,
    const ExactSourceKernel &kernel) {
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::LeafAbsolute;
  source_program.source_id = kernel.source_id;
  source_program.leaf_index = kernel.leaf_index;
  if (kernel.leaf_index != semantic::kInvalidIndex &&
      static_cast<std::size_t>(kernel.leaf_index) <
          plan->program.leaf_descriptors.size()) {
    const auto &leaf =
        plan->program.leaf_descriptors[
            static_cast<std::size_t>(kernel.leaf_index)];
    source_program.leaf_dist_kind = leaf.dist_kind;
  }
  return push_source_product_program(&plan->compiled_math, source_program);
}

inline semantic::Index compile_source_product_onset_program(
    ExactVariantBuildState *plan,
    const ExactSourceKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::OnsetConvolution;
  source_program.source_id = kernel.source_id;
  source_program.leaf_index = kernel.leaf_index;
  source_program.onset_source_program_id =
      compile_source_product_base_program(
          plan, kernel.onset_source_id, condition_id, source_view_id);
  if (kernel.leaf_index != semantic::kInvalidIndex &&
      static_cast<std::size_t>(kernel.leaf_index) <
          plan->program.leaf_descriptors.size()) {
    const auto &leaf =
        plan->program.leaf_descriptors[
            static_cast<std::size_t>(kernel.leaf_index)];
    source_program.leaf_dist_kind = leaf.dist_kind;
    source_program.leaf_onset_lag = leaf.onset_lag;
  }
  return push_source_product_program(&plan->compiled_math, source_program);
}

inline semantic::Index compile_source_product_pool_program(
    ExactVariantBuildState *plan,
    const ExactSourceKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  auto &program = plan->compiled_math;
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::PoolKOfN;
  source_program.source_id = kernel.source_id;
  source_program.pool_k = kernel.pool_k;
  std::vector<semantic::Index> member_programs;
  member_programs.reserve(
      static_cast<std::size_t>(kernel.pool_member_count));
  const auto member_end = kernel.pool_member_offset + kernel.pool_member_count;
  for (semantic::Index i = kernel.pool_member_offset; i < member_end; ++i) {
    const auto member_source =
        plan->program.pool_member_source_ids[
            static_cast<std::size_t>(i)];
    member_programs.push_back(
        compile_source_product_base_program(
            plan, member_source, condition_id, source_view_id));
  }
  const auto member_offset = static_cast<semantic::Index>(
      program.source_program_members.size());
  program.source_program_members.insert(
      program.source_program_members.end(),
      member_programs.begin(),
      member_programs.end());
  source_program.member_programs = CompiledMathIndexSpan{
      member_offset,
      kernel.pool_member_count};
  return push_source_product_program(&program, source_program);
}

inline semantic::Index compile_source_product_kernel_program(
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  const ExactSourceProgramCompileKey key{
      source_id,
      condition_id,
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id};
  const auto existing = plan->source_kernel_program_index.find(key);
  if (existing != plan->source_kernel_program_index.end()) {
    return existing->second;
  }
  semantic::Index program_id = semantic::kInvalidIndex;
  if (source_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(source_id) >= plan->source_kernels.size()) {
    CompiledMathSourceProductProgram source_program;
    source_program.kind = CompiledMathSourceProductProgramKind::ConstantZero;
    program_id =
        push_source_product_program(&plan->compiled_math, source_program);
  } else {
    const auto &kernel =
        plan->source_kernels[static_cast<std::size_t>(source_id)];
    switch (kernel.kind) {
    case CompiledSourceChannelKernelKind::LeafAbsolute:
      program_id = compile_source_product_leaf_program(
          plan, kernel);
      break;
    case CompiledSourceChannelKernelKind::LeafOnsetConvolution:
      program_id = compile_source_product_onset_program(
          plan, kernel, condition_id, key.source_view_id);
      break;
    case CompiledSourceChannelKernelKind::PoolKOfN:
      program_id = compile_source_product_pool_program(
          plan, kernel, condition_id, key.source_view_id);
      break;
    case CompiledSourceChannelKernelKind::Invalid:
      break;
    }
    if (program_id == semantic::kInvalidIndex) {
      CompiledMathSourceProductProgram source_program;
      source_program.kind = CompiledMathSourceProductProgramKind::ConstantZero;
      program_id =
          push_source_product_program(&plan->compiled_math, source_program);
    }
  }
  plan->source_kernel_program_index.emplace(key, program_id);
  return program_id;
}

inline semantic::Index compile_source_product_base_program(
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index source_view_id) {
  const ExactSourceProgramCompileKey key{
      source_id,
      condition_id,
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id};
  const auto existing = plan->source_base_program_index.find(key);
  if (existing != plan->source_base_program_index.end()) {
    return existing->second;
  }
  const auto kernel_program_id = compile_source_product_kernel_program(
      plan, source_id, condition_id, key.source_view_id);
  semantic::Index program_id = kernel_program_id;
  if (source_id == semantic::kInvalidIndex) {
    plan->source_base_program_index.emplace(key, program_id);
    return program_id;
  }
  program_id = compile_source_product_exact_gate_program(
      plan,
      source_id,
      key.source_view_id,
      kernel_program_id);
  plan->source_base_program_index.emplace(key, program_id);
  return program_id;
}

inline semantic::Index compile_source_product_channel_program(
    ExactVariantBuildState *plan,
    CompiledMathSourceProductChannel *channel) {
  if (channel->source_product_program_id != semantic::kInvalidIndex) {
    return channel->source_product_program_id;
  }
  const ExactConditionedSourceProgramCompileKey key{
      channel->source_id,
      channel->condition_id,
      channel->source_view_id == semantic::kInvalidIndex
          ? 0
          : channel->source_view_id,
      channel->time_id,
      channel->time_cap_id};
  const auto existing = plan->source_conditioned_program_index.find(key);
  if (existing != plan->source_conditioned_program_index.end()) {
    channel->source_product_program_id = existing->second;
    return channel->source_product_program_id;
  }
  const auto child_program_id =
      compile_source_product_kernel_program(
          plan,
          key.source_id,
          key.condition_id,
          key.source_view_id);
  CompiledMathSourceProductProgram source_program;
  source_program.kind = CompiledMathSourceProductProgramKind::Conditioned;
  source_program.source_id = key.source_id;
  source_program.child_program_id = child_program_id;
  source_program.static_source_view_relation =
      channel->static_source_view_relation;
  source_program.has_static_source_view_relation =
      channel->has_static_source_view_relation;
  channel->source_product_program_id =
      push_source_product_program(&plan->compiled_math, source_program);
  plan->source_conditioned_program_index.emplace(
      key, channel->source_product_program_id);
  return channel->source_product_program_id;
}

inline int source_product_forced_value(
    const CompiledMathSourceProductChannel &channel,
    const std::uint8_t value_mask) noexcept {
  if (!channel.has_static_source_view_relation) {
    return -1;
  }
  const auto relation =
      static_cast<ExactRelation>(channel.static_source_view_relation);
  if (relation == ExactRelation::Unknown ||
      (relation == ExactRelation::At && value_mask == 1U)) {
    return -1;
  }
  if (relation == ExactRelation::Before || relation == ExactRelation::At) {
    return value_mask == 2U ? 1 : 0;
  }
  if (relation == ExactRelation::After) {
    return value_mask == 4U ? 1 : 0;
  }
  return -1;
}

inline CompiledMathSourceProductOps compile_source_product_ops_for_factor_span(
    ExactVariantBuildState *plan,
    const CompiledMathIndexSpan factors) {
  auto *program = &plan->compiled_math;
  const auto offset = static_cast<semantic::Index>(
      program->source_product_ops.size());
  std::size_t pdf_count = 0U;
  for (semantic::Index i = 0; i < factors.size; ++i) {
    const auto factor_pos =
        static_cast<std::size_t>(factors.offset + i);
    const auto &factor =
        program->source_value_factors[factor_pos];
    const auto channel_pos =
        static_cast<std::size_t>(factor.source_product_channel_id);
    const auto &channel =
        program->source_product_channels[channel_pos];
    const auto value_mask =
        compiled_math_source_factor_channel_mask(factor.kind);
    const int forced_value =
        value_mask == 0U ? 0 : source_product_forced_value(channel, value_mask);
    if (forced_value == 0) {
      program->source_product_ops.resize(
          static_cast<std::size_t>(offset));
      program->source_product_ops.push_back(
          CompiledMathSourceProductOp{});
      return CompiledMathSourceProductOps{offset, 1, false};
    }
    if (forced_value == 1) {
      continue;
    }
    const auto program_id =
        compile_source_product_channel_program(
            plan,
            &program->source_product_channels[channel_pos]);
    CompiledMathSourceProductOp op;
    op.source_product_program_id = program_id;
    op.time_id = channel.time_id;
    op.time_cap_id = channel.time_cap_id;
    op.value_channel_mask = value_mask;
    op.fill_channel_mask = value_mask;
    program->source_product_ops.push_back(op);
    pdf_count += factor.kind == CompiledMathNodeKind::SourcePdf;
  }
  return CompiledMathSourceProductOps{
      offset,
      static_cast<semantic::Index>(
          program->source_product_ops.size() -
          static_cast<std::size_t>(offset)),
      pdf_count > 1U};
}

inline void compile_source_product_execution_ops(
    ExactVariantBuildState *plan,
    CompiledMathExecutionPlan *execution) {
  auto &program = plan->compiled_math;
  if (execution->kind == CompiledMathExecutionKind::SourceProduct) {
    execution->source_product_ops =
        compile_source_product_ops_for_factor_span(
            plan, execution->source_value_factors);
    return;
  }
  if (execution->kind != CompiledMathExecutionKind::SourceProductSum) {
    return;
  }
  execution->source_product_ops =
      compile_source_product_ops_for_factor_span(
          plan, execution->source_value_factors);
  for (semantic::Index i = 0; i < execution->source_product_terms.size; ++i) {
    auto &term = program.source_product_terms[
        static_cast<std::size_t>(execution->source_product_terms.offset + i)];
    term.source_product_ops =
        compile_source_product_ops_for_factor_span(
            plan, term.source_value_factors);
  }
}

inline void compile_source_product_execution_programs(
    ExactVariantBuildState *plan) {
  auto &program = plan->compiled_math;
  plan->source_kernel_program_index.clear();
  plan->source_base_program_index.clear();
  plan->source_conditioned_program_index.clear();
  program.source_product_ops.clear();
  program.source_programs.clear();
  program.source_program_members.clear();
  for (auto &channel : program.source_product_channels) {
    channel.source_product_program_id = semantic::kInvalidIndex;
  }
  for (auto &kernel : program.integral_kernels) {
    compile_source_product_execution_ops(plan, &kernel.execution);
    compile_source_product_execution_ops(plan, &kernel.initial_execution);
  }
  for (auto &root : program.roots) {
    compile_source_product_execution_ops(plan, &root.execution);
    compile_source_product_execution_ops(plan, &root.initial_execution);
  }
}

inline semantic::Index compile_source_node_program(
    ExactVariantBuildState *plan,
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

inline void compile_source_node_programs(ExactVariantBuildState *plan) {
  auto &program = plan->compiled_math;
  for (auto &node : program.nodes) {
    if (!compiled_math_is_source_value_node(node.kind)) {
      continue;
    }
    node.source_program_id = compile_source_node_program(plan, node);
  }
}

inline void finalize_source_program_initial_resolutions(
    CompiledMathProgram *program) {
  for (std::size_t i = 0; i < program->source_programs.size(); ++i) {
    auto &source_program = program->source_programs[i];
    const auto program_id = static_cast<semantic::Index>(i);
    source_program.initial_without_pdf_program_id = program_id;
    source_program.initial_with_pdf_program_id = program_id;
    if (source_program.kind ==
        CompiledMathSourceProductProgramKind::ConstantZero) {
      source_program.initial_without_pdf_program_id = semantic::kInvalidIndex;
      source_program.initial_with_pdf_program_id = semantic::kInvalidIndex;
      continue;
    }
    const bool wrapper =
        source_program.kind ==
            CompiledMathSourceProductProgramKind::ExactGate ||
        source_program.kind ==
            CompiledMathSourceProductProgramKind::Conditioned;
    if (!wrapper) {
      continue;
    }
    const auto relation = source_program.has_static_source_view_relation
                              ? static_cast<ExactRelation>(
                                    source_program.static_source_view_relation)
                              : ExactRelation::Unknown;
    const auto &child = program->source_programs[
        static_cast<std::size_t>(source_program.child_program_id)];
    if (relation == ExactRelation::Unknown) {
      source_program.initial_with_pdf_program_id =
          child.initial_with_pdf_program_id;
      source_program.initial_without_pdf_program_id =
          child.initial_without_pdf_program_id;
    } else if (relation == ExactRelation::Before) {
      source_program.initial_with_pdf_program_id =
          kInitialCertainSourceProgramId;
      source_program.initial_without_pdf_program_id =
          kInitialCertainSourceProgramId;
    } else if (relation == ExactRelation::At) {
      source_program.initial_without_pdf_program_id =
          kInitialCertainSourceProgramId;
    } else {
      source_program.initial_with_pdf_program_id = semantic::kInvalidIndex;
      source_program.initial_without_pdf_program_id = semantic::kInvalidIndex;
    }
  }
}

template <typename Visitor>
inline void visit_source_product_execution_ops(
    const CompiledMathProgram &program,
    const CompiledMathExecutionPlan &execution,
    Visitor &&visit) {
  if (execution.kind == CompiledMathExecutionKind::SourceProduct) {
    visit(execution.source_product_ops);
    return;
  }
  if (execution.kind != CompiledMathExecutionKind::SourceProductSum) {
    return;
  }
  visit(execution.source_product_ops);
  for (semantic::Index i = 0; i < execution.source_product_terms.size; ++i) {
    const auto &term = program.source_product_terms[
        static_cast<std::size_t>(execution.source_product_terms.offset + i)];
    visit(term.source_product_ops);
  }
}

inline void compile_source_program_cache_slots(CompiledMathProgram *program) {
  const auto program_count = program->source_programs.size();
  program->source_program_cache_slots.assign(
      program_count, semantic::kInvalidIndex);
  std::vector<std::uint8_t> reused(program_count, 0U);
  std::vector<std::size_t> schedule_count(program_count, 0U);
  std::vector<std::size_t> op_count(program_count, 0U);
  std::vector<std::uint8_t> op_fill_mask(program_count, 0U);
  for (auto &node : program->nodes) {
    node.cache_source_program = false;
  }
  for (auto &op : program->source_product_ops) {
    op.cache_result = false;
  }
  for (const auto &root : program->roots) {
    if (root.execution.kind != CompiledMathExecutionKind::Schedule) {
      continue;
    }
    std::fill(schedule_count.begin(), schedule_count.end(), 0U);
    for (semantic::Index i = 0; i < root.schedule.size; ++i) {
      const auto node_id = program->root_schedule_nodes[
          static_cast<std::size_t>(root.schedule.offset + i)];
      const auto &node = program->nodes[static_cast<std::size_t>(node_id)];
      if (compiled_math_is_source_value_node(node.kind) &&
          node.source_program_id != semantic::kInvalidIndex) {
        const auto source_program =
            static_cast<std::size_t>(node.source_program_id);
        if (++schedule_count[source_program] > 1U) {
          reused[source_program] = 1U;
        }
      }
    }
    for (semantic::Index i = 0; i < root.schedule.size; ++i) {
      const auto node_id = program->root_schedule_nodes[
          static_cast<std::size_t>(root.schedule.offset + i)];
      auto &node = program->nodes[static_cast<std::size_t>(node_id)];
      if (!compiled_math_is_source_value_node(node.kind) ||
          node.source_program_id == semantic::kInvalidIndex) {
        continue;
      }
      const auto source_program =
          static_cast<std::size_t>(node.source_program_id);
      if (schedule_count[source_program] > 1U) {
        node.cache_source_program = true;
        reused[source_program] = 1U;
      }
    }
  }
  const auto mark_execution = [&](const CompiledMathExecutionPlan &execution) {
    std::fill(op_count.begin(), op_count.end(), 0U);
    std::fill(op_fill_mask.begin(), op_fill_mask.end(), 0U);
    visit_source_product_execution_ops(
        *program,
        execution,
        [&](const CompiledMathSourceProductOps ops) {
          for (semantic::Index i = 0; i < ops.size; ++i) {
            const auto &op = program->source_product_ops[
                static_cast<std::size_t>(ops.offset + i)];
            if (op.value_channel_mask != 0U &&
                op.source_product_program_id != semantic::kInvalidIndex) {
              const auto program_id = static_cast<std::size_t>(
                  op.source_product_program_id);
              ++op_count[program_id];
              op_fill_mask[program_id] |= op.value_channel_mask;
            }
          }
        });
    visit_source_product_execution_ops(
        *program,
        execution,
        [&](const CompiledMathSourceProductOps ops) {
          for (semantic::Index i = 0; i < ops.size; ++i) {
            auto &op = program->source_product_ops[
                static_cast<std::size_t>(ops.offset + i)];
            if (op.value_channel_mask == 0U ||
                op.source_product_program_id == semantic::kInvalidIndex) {
              continue;
            }
            const auto program_id = static_cast<std::size_t>(
                op.source_product_program_id);
            if (op_count[program_id] > 1U) {
              op.cache_result = true;
              op.fill_channel_mask = op_fill_mask[program_id];
              reused[program_id] = 1U;
            }
          }
        });
  };
  for (const auto &kernel : program->integral_kernels) {
    mark_execution(kernel.execution);
    mark_execution(kernel.initial_execution);
  }
  for (const auto &root : program->roots) {
    mark_execution(root.execution);
    mark_execution(root.initial_execution);
  }
  semantic::Index next_slot = 0;
  for (std::size_t i = 0; i < program_count; ++i) {
    if (reused[i] != 0U) {
      program->source_program_cache_slots[i] = next_slot++;
    }
  }
  program->source_program_cache_count = next_slot;
}

inline void validate_source_product_relations_materialized(
    const CompiledMathProgram &program) {
  for (std::size_t i = 0;
       i < program.source_product_channels.size();
       ++i) {
    const auto &channel = program.source_product_channels[i];
    if (channel.source_id != semantic::kInvalidIndex &&
        !channel.has_static_source_view_relation) {
      throw std::runtime_error(
          "source-product channel " + std::to_string(i) +
          " has no compiled source-view relation");
    }
  }
  for (std::size_t i = 0;
       i < program.source_programs.size();
       ++i) {
    const auto &source_program =
        program.source_programs[i];
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
       i < program.source_product_ops.size();
       ++i) {
    const auto &op = program.source_product_ops[i];
    if (op.value_channel_mask != 0U &&
        (op.source_product_program_id == semantic::kInvalidIndex ||
         static_cast<std::size_t>(op.source_product_program_id) >=
             program.source_programs.size())) {
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
            program.source_programs.size()) {
      throw std::runtime_error(
          "source node " + std::to_string(i) +
          " has no compiled source arithmetic program");
    }
  }
}

inline void finalize_compiled_math_time_slots(CompiledMathProgram *program) {
  semantic::Index count =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) + 1U;
  const auto include = [&](const semantic::Index time_id) {
    if (time_id != semantic::kInvalidIndex && time_id >= count) {
      count = time_id + 1U;
    }
  };
  for (const auto &node : program->nodes) {
    include(node.time_id);
    if (compiled_math_is_source_value_node(node.kind) ||
        node.kind == CompiledMathNodeKind::TimeGate ||
        node.kind == CompiledMathNodeKind::StrictTimeGate) {
      include(node.aux_id);
    }
  }
  for (const auto &kernel : program->integral_kernels) {
    include(kernel.bind_time_id);
  }
  for (const auto &op : program->source_product_ops) {
    if (op.value_channel_mask == 0U) {
      continue;
    }
    include(op.time_id);
    include(op.time_cap_id);
  }
  program->time_slot_count = count;
}

inline void compile_source_view_relation_tables(ExactVariantBuildState *plan) {
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
