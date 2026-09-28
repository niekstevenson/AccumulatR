#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline void compile_exact_support_context(ExactVariantBuildState *plan) {
  ExactSupportBuilder builder(plan->program);
  plan->leaf_supports = builder.build_leaf_supports();
  plan->pool_supports = builder.build_pool_supports();
  plan->expr_supports = builder.build_expr_supports();
  plan->source_count = static_cast<semantic::Index>(
      plan->program.layout.n_leaves +
      plan->program.layout.n_pools);
}

inline void compile_program_source_runtime_fields(ExactVariantBuildState *plan) {
  auto &program = plan->program;

  for (semantic::Index i = 0; i < program.layout.n_leaves; ++i) {
    const auto pos = static_cast<std::size_t>(i);
    auto &leaf = program.leaf_descriptors[pos];
    const auto source_id = source_ordinal(
        *plan,
        static_cast<semantic::SourceKind>(leaf.onset_source_kind),
        leaf.onset_source_index);
    leaf.onset_source_id = source_id;
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

inline void analyze_aggregate_pool_transition_safety(
    ExactVariantBuildState *plan) {
  const auto &program = plan->program;
  std::vector<std::uint8_t> used_event_sources(
      static_cast<std::size_t>(plan->source_count), 0U);
  for (std::size_t i = 0; i < program.expr_kind.size(); ++i) {
    if (static_cast<semantic::ExprKind>(program.expr_kind[i]) !=
        semantic::ExprKind::Event) {
      continue;
    }
    const auto source_id = program.expr_source_ids[i];
    if (source_id >= 0 && source_id < plan->source_count) {
      used_event_sources[static_cast<std::size_t>(source_id)] = 1U;
    }
  }

  const auto leaf_count = plan->program.layout.n_leaves;
  plan->pool_transition_can_stay_aggregate.assign(
      static_cast<std::size_t>(program.layout.n_pools), 1U);
  for (semantic::Index pool_index = 0;
       pool_index < program.layout.n_pools;
       ++pool_index) {
    const auto pool_source_id = leaf_count + pool_index;
    const auto &pool_support =
        plan->pool_supports[static_cast<std::size_t>(pool_index)];
    for (semantic::Index other_source_id = 0;
         other_source_id < plan->source_count;
         ++other_source_id) {
      if (other_source_id == pool_source_id ||
          used_event_sources[static_cast<std::size_t>(other_source_id)] == 0U) {
        continue;
      }
      const auto &other_support =
          other_source_id < leaf_count
              ? plan->leaf_supports[static_cast<std::size_t>(other_source_id)]
              : plan->pool_supports[
                    static_cast<std::size_t>(other_source_id - leaf_count)];
      if (supports_overlap(pool_support, other_support)) {
        plan->pool_transition_can_stay_aggregate[
            static_cast<std::size_t>(pool_index)] = 0U;
        break;
      }
    }
  }
}

inline void compile_source_kernels(ExactVariantBuildState *plan) {
  const auto &program = plan->program;
  plan->source_kernels.assign(
      static_cast<std::size_t>(plan->source_count), ExactSourceKernel{});

  for (semantic::Index i = 0; i < program.layout.n_leaves; ++i) {
    const auto pos = static_cast<std::size_t>(i);
    const auto &leaf = program.leaf_descriptors[pos];
    auto &kernel = plan->source_kernels[pos];
    kernel.source_id = i;
    kernel.leaf_index = i;
    kernel.onset_source_id = leaf.onset_source_id;
    kernel.kind = static_cast<semantic::OnsetKind>(leaf.onset_kind) ==
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
    kernel.pool_member_offset = program.pool_member_offsets[pos];
    kernel.pool_member_count =
        program.pool_member_offsets[pos + 1U] -
        program.pool_member_offsets[pos];
    kernel.pool_k = program.pool_k[pos];
  }
}

inline void compile_exact_expr_kernels(ExactVariantBuildState *plan) {
  const auto &program = plan->program;
  plan->expr_kernels.assign(program.expr_kind.size(), ExactExprKernel{});

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

    if (kernel.kind != semantic::ExprKind::Guard) {
      continue;
    }

    kernel.guard_ref_expr_id = program.expr_ref_child[pos];
    kernel.guard_blocker_expr_id = program.expr_blocker_child[pos];
  }
}
inline void compile_shared_trigger_state_table(ExactVariantBuildState *plan) {
  struct TriggerStateBuilder {
    std::vector<ExactCompiledTriggerWeightTerm> weight_terms;
    std::vector<std::uint8_t> shared_started;
  };

  auto &table = plan->trigger_state_table;
  table.states.clear();
  table.weight_terms.clear();
  table.shared_started_values.clear();
  table.trigger_count = plan->program.layout.n_triggers;

  std::vector<TriggerStateBuilder> builders;
  builders.push_back(TriggerStateBuilder{});
  builders.front().shared_started.assign(
      static_cast<std::size_t>(table.trigger_count), 2U);

  const auto &program = plan->program;
  for (semantic::Index trigger_index = 0;
       trigger_index < program.layout.n_triggers; ++trigger_index) {
    const auto trigger_pos = static_cast<std::size_t>(trigger_index);
    const auto member_begin = program.trigger_member_offsets[trigger_pos];
    const auto member_end = program.trigger_member_offsets[trigger_pos + 1U];
    if (member_end - member_begin <= 1) continue;
    const auto q_leaf_index = program.trigger_member_indices[member_begin];

    std::vector<TriggerStateBuilder> next;
    next.reserve(builders.size() * 2U);
    for (const auto &builder : builders) {
      auto append_variable_state =
          [&](const std::uint8_t shared_started) {
            auto out = builder;
            out.shared_started[trigger_pos] = shared_started;
            out.weight_terms.push_back(
                ExactCompiledTriggerWeightTerm{
                    q_leaf_index, shared_started});
            next.push_back(std::move(out));
          };

      append_variable_state(0U);
      append_variable_state(1U);
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
            ExactIndexSpan{
                weight_offset,
                static_cast<semantic::Index>(
                    table.weight_terms.size() -
                    static_cast<std::size_t>(weight_offset))},
            shared_offset});
  }
}

} // namespace detail
} // namespace accumulatr::eval
