#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline void compile_exact_support_context(ExactVariantPlan *plan) {
  ExactSupportBuilder builder(plan->lowered);
  plan->leaf_supports = builder.build_leaf_supports();
  plan->pool_supports = builder.build_pool_supports();
  plan->expr_supports = builder.build_expr_supports();
  plan->source_count = static_cast<semantic::Index>(
      plan->lowered.program.layout.n_leaves +
      plan->lowered.program.layout.n_pools);
  plan->leaf_source_ids.resize(
      static_cast<std::size_t>(plan->lowered.program.layout.n_leaves));
  plan->pool_source_ids.resize(
      static_cast<std::size_t>(plan->lowered.program.layout.n_pools));
  for (semantic::Index i = 0; i < plan->lowered.program.layout.n_leaves; ++i) {
    plan->leaf_source_ids[static_cast<std::size_t>(i)] = i;
  }
  for (semantic::Index i = 0; i < plan->lowered.program.layout.n_pools; ++i) {
    plan->pool_source_ids[static_cast<std::size_t>(i)] =
        static_cast<semantic::Index>(plan->lowered.program.layout.n_leaves + i);
  }
}

inline void compile_program_source_runtime_fields(ExactVariantPlan *plan) {
  auto &program = plan->lowered.program;

  for (semantic::Index i = 0; i < program.layout.n_leaves; ++i) {
    const auto pos = static_cast<std::size_t>(i);
    const auto source_id = source_ordinal(
        *plan,
        static_cast<semantic::SourceKind>(program.onset_source_kind[pos]),
        program.onset_source_index[pos]);
    program.onset_source_ids[pos] = source_id;
    program.leaf_descriptors[pos].onset_source_id = source_id;
  }

  for (std::size_t i = 0; i < program.pool_member_indices.size(); ++i) {
    program.pool_member_source_ids[i] = source_ordinal(
        *plan,
        static_cast<semantic::SourceKind>(program.pool_member_kind[i]),
        program.pool_member_indices[i]);
  }

  for (std::size_t i = 0; i < program.expr_source_index.size(); ++i) {
    program.expr_source_ids[i] = source_ordinal(
        *plan,
        static_cast<semantic::SourceKind>(program.expr_source_kind[i]),
        program.expr_source_index[i]);
  }
}

inline void compile_source_kernels(ExactVariantPlan *plan) {
  const auto &program = plan->lowered.program;
  plan->source_kernels.assign(
      static_cast<std::size_t>(plan->source_count), ExactSourceKernel{});

  for (semantic::Index i = 0; i < program.layout.n_leaves; ++i) {
    const auto pos = static_cast<std::size_t>(i);
    auto &kernel = plan->source_kernels[pos];
    kernel.source_id = i;
    kernel.leaf_index = i;
    kernel.onset_source_id = program.onset_source_ids[pos];
    kernel.kind = static_cast<semantic::OnsetKind>(
                      program.onset_kind[pos]) ==
                          semantic::OnsetKind::Absolute
                      ? CompiledSourceChannelKernelKind::LeafAbsolute
                      : CompiledSourceChannelKernelKind::LeafOnsetConvolution;
  }

  for (semantic::Index i = 0; i < program.layout.n_pools; ++i) {
    const auto source_id =
        static_cast<semantic::Index>(program.layout.n_leaves + i);
    const auto pos = static_cast<std::size_t>(i);
    auto &kernel =
        plan->source_kernels[static_cast<std::size_t>(source_id)];
    kernel.kind = CompiledSourceChannelKernelKind::PoolKOfN;
    kernel.source_id = source_id;
    kernel.pool_index = i;
    kernel.pool_member_offset = program.pool_member_offsets[pos];
    kernel.pool_member_count =
        program.pool_member_offsets[pos + 1U] -
        program.pool_member_offsets[pos];
    kernel.pool_k = program.pool_k[pos];
  }
}

inline void append_expr_union_subset(
    ExactVariantPlan *plan,
    const std::vector<semantic::Index> &children,
    const std::vector<semantic::Index> &child_positions) {
  const auto offset = static_cast<semantic::Index>(
      plan->expr_union_subset_children.size());
  plan->expr_union_subset_children.insert(
      plan->expr_union_subset_children.end(),
      children.begin(),
      children.end());
  const auto position_offset = static_cast<semantic::Index>(
      plan->expr_union_subset_child_positions.size());
  plan->expr_union_subset_child_positions.insert(
      plan->expr_union_subset_child_positions.end(),
      child_positions.begin(),
      child_positions.end());
  plan->expr_union_subsets.push_back(
      ExactExprUnionSubset{
          ExactIndexSpan{
              offset,
              static_cast<semantic::Index>(children.size())},
          ExactIndexSpan{
              position_offset,
              static_cast<semantic::Index>(child_positions.size())},
          children.size() % 2U == 1U ? 1 : -1});
}

inline void enumerate_expr_union_subsets(
    ExactVariantPlan *plan,
    const std::vector<semantic::Index> &children,
    const std::size_t start,
    std::vector<semantic::Index> *current,
    std::vector<semantic::Index> *current_positions) {
  for (std::size_t i = start; i < children.size(); ++i) {
    current->push_back(children[i]);
    current_positions->push_back(static_cast<semantic::Index>(i));
    append_expr_union_subset(plan, *current, *current_positions);
    enumerate_expr_union_subsets(
        plan, children, i + 1U, current, current_positions);
    current_positions->pop_back();
    current->pop_back();
  }
}

inline ExactIndexSpan compile_expr_union_subsets(
    ExactVariantPlan *plan,
    const ExactIndexSpan children_span) {
  const auto offset =
      static_cast<semantic::Index>(plan->expr_union_subsets.size());
  std::vector<semantic::Index> children;
  children.reserve(static_cast<std::size_t>(children_span.size));
  const auto &program = plan->lowered.program;
  for (semantic::Index i = 0; i < children_span.size; ++i) {
    children.push_back(
        program.expr_args[static_cast<std::size_t>(children_span.offset + i)]);
  }
  std::vector<semantic::Index> current;
  current.reserve(children.size());
  std::vector<semantic::Index> current_positions;
  current_positions.reserve(children.size());
  enumerate_expr_union_subsets(
      plan, children, 0U, &current, &current_positions);
  return ExactIndexSpan{
      offset,
      static_cast<semantic::Index>(
          plan->expr_union_subsets.size() - static_cast<std::size_t>(offset))};
}

inline bool expr_children_supports_overlap(
    const ExactVariantPlan &plan,
    const ExactIndexSpan children_span) {
  const auto &program = plan.lowered.program;
  for (semantic::Index i = 0; i < children_span.size; ++i) {
    const auto lhs_idx =
        program.expr_args[static_cast<std::size_t>(children_span.offset + i)];
    const auto &lhs = plan.expr_supports[static_cast<std::size_t>(lhs_idx)];
    for (semantic::Index j = i + 1; j < children_span.size; ++j) {
      const auto rhs_idx =
          program.expr_args[static_cast<std::size_t>(children_span.offset + j)];
      const auto &rhs = plan.expr_supports[static_cast<std::size_t>(rhs_idx)];
      if (supports_overlap(lhs, rhs)) {
        return true;
      }
    }
  }
  return false;
}

inline bool expr_static_subset_of_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id);

inline bool expr_static_can_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id);

inline bool expr_span_static_has_event_source_subset(
    const ExactVariantPlan &plan,
    const ExactIndexSpan span,
    const semantic::Index source_id) {
  const auto &program = plan.lowered.program;
  for (semantic::Index i = span.offset; i < span.offset + span.size; ++i) {
    if (expr_static_subset_of_event_source(
            plan,
            program.expr_args[static_cast<std::size_t>(i)],
            source_id)) {
      return true;
    }
  }
  return false;
}

inline bool expr_static_subset_of_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id) {
  if (expr_idx == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_idx);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  switch (kind) {
  case semantic::ExprKind::Event:
    return program.expr_source_ids[pos] == source_id;
  case semantic::ExprKind::And:
    return expr_span_static_has_event_source_subset(
        plan,
        ExactIndexSpan{
            program.expr_arg_offsets[pos],
            static_cast<semantic::Index>(
                program.expr_arg_offsets[pos + 1U] -
                program.expr_arg_offsets[pos])},
        source_id);
  case semantic::ExprKind::Or: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    if (begin == end) {
      return false;
    }
    for (semantic::Index i = begin; i < end; ++i) {
      if (!expr_static_subset_of_event_source(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              source_id)) {
        return false;
      }
    }
    return true;
  }
  case semantic::ExprKind::Guard:
    return expr_static_subset_of_event_source(
        plan,
        program.expr_ref_child[pos],
        source_id);
  case semantic::ExprKind::Impossible:
    return true;
  case semantic::ExprKind::TrueExpr:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool expr_static_can_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id) {
  if (expr_idx == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_idx);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind == semantic::ExprKind::Event) {
    return program.expr_source_ids[pos] == source_id;
  }
  if (kind == semantic::ExprKind::Guard) {
    const auto &kernel = plan.expr_kernels[pos];
    return kernel.simple_event_guard &&
           !kernel.has_unless &&
           kernel.guard_ref_source_id == source_id;
  }
  if (kind != semantic::ExprKind::And) {
    return false;
  }
  const auto begin = program.expr_arg_offsets[pos];
  const auto end = program.expr_arg_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    if (expr_static_can_equiv_event_source(
            plan,
            program.expr_args[static_cast<std::size_t>(i)],
            source_id)) {
      return true;
    }
  }
  return false;
}

inline ExactIndexSpan compile_expr_union_absorption_candidates(
    ExactVariantPlan *plan,
    const semantic::Index expr_idx,
    const ExactExprKernel &kernel) {
  const auto offset = static_cast<semantic::Index>(
      plan->expr_union_absorption_sources.size());
  const auto &program = plan->lowered.program;
  const auto &support = plan->expr_supports[static_cast<std::size_t>(expr_idx)];
  for (const auto source_id : support) {
    bool all_subset = true;
    bool any_equiv = false;
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child = program.expr_args[
          static_cast<std::size_t>(kernel.children.offset + i)];
      if (!expr_static_subset_of_event_source(*plan, child, source_id)) {
        all_subset = false;
        break;
      }
      any_equiv =
          any_equiv ||
          expr_static_can_equiv_event_source(*plan, child, source_id);
    }
    if (all_subset && any_equiv) {
      plan->expr_union_absorption_sources.push_back(source_id);
    }
  }
  return ExactIndexSpan{
      offset,
      static_cast<semantic::Index>(
          plan->expr_union_absorption_sources.size() -
          static_cast<std::size_t>(offset))};
}

inline void compile_exact_expr_kernels(ExactVariantPlan *plan) {
  const auto &program = plan->lowered.program;
  plan->expr_kernels.assign(program.expr_kind.size(), ExactExprKernel{});
  plan->expr_union_subsets.clear();
  plan->expr_union_subset_children.clear();
  plan->expr_union_subset_child_positions.clear();
  plan->expr_union_absorption_sources.clear();
  plan->expr_union_kernel_cache_slot_count = 0;

  for (semantic::Index expr_idx = 0;
       expr_idx < static_cast<semantic::Index>(program.expr_kind.size());
       ++expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    auto &kernel = plan->expr_kernels[pos];
    kernel.kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
    kernel.children = ExactIndexSpan{
        program.expr_arg_offsets[pos],
        static_cast<semantic::Index>(
            program.expr_arg_offsets[pos + 1U] -
            program.expr_arg_offsets[pos])};

    if (kernel.kind == semantic::ExprKind::Event) {
      kernel.event_source_id = program.expr_source_ids[pos];
      continue;
    }

    if (kernel.kind == semantic::ExprKind::Or) {
      kernel.children_overlap =
          expr_children_supports_overlap(*plan, kernel.children);
      if (kernel.children_overlap) {
        kernel.union_subset_span =
            compile_expr_union_subsets(plan, kernel.children);
        kernel.union_absorption_candidate_span =
            compile_expr_union_absorption_candidates(
                plan, expr_idx, kernel);
        kernel.union_kernel_cache_slot =
            plan->expr_union_kernel_cache_slot_count++;
      }
      continue;
    }

    if (kernel.kind != semantic::ExprKind::Guard) {
      continue;
    }

    kernel.guard_ref_expr_id = program.expr_ref_child[pos];
    kernel.guard_blocker_expr_id = program.expr_blocker_child[pos];
    kernel.has_unless = !kernel.children.empty();
    const auto ref_pos =
        static_cast<std::size_t>(kernel.guard_ref_expr_id);
    const auto blocker_pos =
        static_cast<std::size_t>(kernel.guard_blocker_expr_id);
    kernel.guard_ref_kind =
        static_cast<semantic::ExprKind>(program.expr_kind[ref_pos]);
    kernel.guard_blocker_kind =
        static_cast<semantic::ExprKind>(program.expr_kind[blocker_pos]);
    kernel.simple_event_guard =
        kernel.guard_ref_kind == semantic::ExprKind::Event &&
        kernel.guard_blocker_kind == semantic::ExprKind::Event;
    if (kernel.simple_event_guard) {
      kernel.guard_ref_source_id = program.expr_source_ids[ref_pos];
      kernel.guard_blocker_source_id = program.expr_source_ids[blocker_pos];
    }
    if (kernel.guard_ref_kind == semantic::ExprKind::And) {
      kernel.guard_ref_child_span = ExactIndexSpan{
          program.expr_arg_offsets[ref_pos],
          static_cast<semantic::Index>(
              program.expr_arg_offsets[ref_pos + 1U] -
              program.expr_arg_offsets[ref_pos])};
    }
  }
}
inline void compile_trigger_state_table(ExactVariantPlan *plan) {
  struct TriggerStateBuilder {
    double fixed_weight{1.0};
    std::vector<ExactCompiledTriggerWeightTerm> weight_terms;
    std::vector<std::uint8_t> shared_started;
  };

  auto &table = plan->trigger_state_table;
  table.states.clear();
  table.weight_terms.clear();
  table.shared_started_values.clear();
  table.trigger_count = plan->lowered.program.layout.n_triggers;

  std::vector<TriggerStateBuilder> builders;
  builders.push_back(TriggerStateBuilder{});
  builders.front().shared_started.assign(
      static_cast<std::size_t>(table.trigger_count), 2U);

  const auto &program = plan->lowered.program;
  for (const auto trigger_index : plan->shared_trigger_indices) {
    const auto trigger_pos = static_cast<std::size_t>(trigger_index);
    const bool fixed_q = program.trigger_has_fixed_q[trigger_pos] != 0U;
    const double raw_fixed_q =
        fixed_q ? program.trigger_fixed_q[trigger_pos] : 0.0;
    const double q_fixed =
        fixed_q
            ? (!std::isfinite(raw_fixed_q)
                   ? 0.0
                   : (raw_fixed_q <= 0.0
                          ? 0.0
                          : (raw_fixed_q >= 1.0 ? 1.0 : raw_fixed_q)))
            : 0.0;
    semantic::Index q_leaf_index{semantic::kInvalidIndex};
    if (!fixed_q) {
      const auto member_begin = program.trigger_member_offsets[trigger_pos];
      const auto member_end = program.trigger_member_offsets[trigger_pos + 1U];
      if (member_begin != member_end) {
        q_leaf_index =
            program.trigger_member_indices[
                static_cast<std::size_t>(member_begin)];
      }
    }

    std::vector<TriggerStateBuilder> next;
    next.reserve(builders.size() * 2U);
    for (const auto &builder : builders) {
      auto append_state =
          [&](const std::uint8_t shared_started, const double fixed_weight) {
            auto out = builder;
            out.fixed_weight *= fixed_weight;
            out.shared_started[trigger_pos] = shared_started;
            next.push_back(std::move(out));
          };
      auto append_variable_state =
          [&](const std::uint8_t shared_started) {
            auto out = builder;
            out.shared_started[trigger_pos] = shared_started;
            out.weight_terms.push_back(
                ExactCompiledTriggerWeightTerm{
                    q_leaf_index, shared_started});
            next.push_back(std::move(out));
          };

      if (fixed_q) {
        if (q_fixed > 0.0) {
          append_state(0U, q_fixed);
        }
        if (q_fixed < 1.0) {
          append_state(1U, 1.0 - q_fixed);
        }
      } else {
        append_variable_state(0U);
        append_variable_state(1U);
      }
    }
    builders.swap(next);
  }

  table.states.reserve(builders.size());
  table.shared_started_values.reserve(
      builders.size() * static_cast<std::size_t>(table.trigger_count));
  for (const auto &builder : builders) {
    const auto shared_offset =
        static_cast<semantic::Index>(table.shared_started_values.size());
    table.shared_started_values.insert(
        table.shared_started_values.end(),
        builder.shared_started.begin(),
        builder.shared_started.end());
    const auto weight_offset =
        static_cast<semantic::Index>(table.weight_terms.size());
    table.weight_terms.insert(
        table.weight_terms.end(),
        builder.weight_terms.begin(),
        builder.weight_terms.end());
    table.states.push_back(
        ExactCompiledTriggerState{
            builder.fixed_weight,
            ExactIndexSpan{
                weight_offset,
                static_cast<semantic::Index>(
                    table.weight_terms.size() -
                    static_cast<std::size_t>(weight_offset))},
            shared_offset});
  }
}

inline void compile_shared_trigger_state_table(ExactVariantPlan *plan) {
  plan->shared_trigger_indices.clear();
  const auto &program = plan->lowered.program;
  for (int i = 0; i < program.layout.n_triggers; ++i) {
    const auto trigger_pos = static_cast<std::size_t>(i);
    if (static_cast<semantic::TriggerKind>(program.trigger_kind[trigger_pos]) ==
            semantic::TriggerKind::Shared &&
        program.trigger_member_offsets[trigger_pos + 1U] -
                program.trigger_member_offsets[trigger_pos] >
            1) {
      plan->shared_trigger_indices.push_back(i);
    }
  }
  compile_trigger_state_table(plan);
}

} // namespace detail
} // namespace accumulatr::eval
