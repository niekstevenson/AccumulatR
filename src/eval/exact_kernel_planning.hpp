#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline void compile_exact_support_context(ExactVariantBuildState *plan) {
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

inline void compile_program_source_runtime_fields(ExactVariantBuildState *plan) {
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

inline void compile_source_kernels(ExactVariantBuildState *plan) {
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

inline void compile_exact_expr_kernels(ExactVariantBuildState *plan) {
  const auto &program = plan->lowered.program;
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
    kernel.has_unless = !kernel.children.empty();
  }
}
inline void compile_trigger_state_table(ExactVariantBuildState *plan) {
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

inline void compile_shared_trigger_state_table(ExactVariantBuildState *plan) {
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
