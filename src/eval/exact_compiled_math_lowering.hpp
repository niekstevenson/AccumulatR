#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool relation_template_equal(const ExactRelationTemplate &lhs,
                                    const ExactRelationTemplate &rhs) {
  return lhs.source_ids == rhs.source_ids && lhs.relations == rhs.relations;
}

inline semantic::Index compile_source_view_id(
    ExactVariantBuildState *plan,
    const ExactRelationTemplate &relation_template) {
  if (relation_template.empty()) {
    return 0;
  }
  for (semantic::Index i = 0;
       i < static_cast<semantic::Index>(plan->compiled_source_views.size());
       ++i) {
    if (relation_template_equal(
            plan->compiled_source_views[static_cast<std::size_t>(i)],
            relation_template)) {
      return i + 1U;
    }
  }
  plan->compiled_source_views.push_back(relation_template);
  return static_cast<semantic::Index>(plan->compiled_source_views.size());
}

inline semantic::Index compile_expr_value_node(
    ExactVariantBuildState *plan,
    semantic::Index expr_id,
    CompiledMathValueKind value_kind,
    semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    semantic::Index source_view_id = 0);

inline semantic::Index compile_expr_distribution_node(
    ExactVariantBuildState *plan,
    semantic::Index expr_id,
    CompiledMathValueKind value_kind,
    semantic::Index time_id,
    semantic::Index source_view_id);

inline semantic::Index compile_expr_source_node(
    ExactVariantBuildState *plan,
    const CompiledMathNodeKind kind,
    const semantic::Index source_id,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0) {
  return compiled_math_source_node(
      &plan->compiled_math,
      kind,
      source_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_integral_zero_to_current_node(
    ExactVariantBuildState *plan,
    const semantic::Index integrand_node,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  const auto integrand_root =
      compiled_math_make_root(&plan->compiled_math, integrand_node);
  return compiled_math_integral_node(
      &plan->compiled_math,
      CompiledMathNodeKind::IntegralZeroToCurrent,
      integrand_root,
      time_id,
      source_view_id,
      bind_time_id);
}

inline semantic::Index compile_outcome_select_node(
    ExactVariantBuildState *plan,
    const std::vector<semantic::Index> &outcome_indices,
    const semantic::Index unused_node,
    const semantic::Index used_node) {
  if (outcome_indices.empty() || unused_node == used_node) return unused_node;
  auto inserted = plan->outcome_predicates.emplace(
      outcome_indices, static_cast<semantic::Index>(
          plan->compiled_outcome_gate_indices.size()));
  if (inserted.second) {
    plan->compiled_outcome_gate_indices.insert(
        plan->compiled_outcome_gate_indices.end(),
        outcome_indices.begin(), outcome_indices.end());
  }
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::OutcomeSelect;
  key.subject_id = inserted.first->second;
  key.aux_id = static_cast<semantic::Index>(outcome_indices.size());
  key.children = {unused_node, used_node};
  const auto branch_root = [&](const semantic::Index child) {
    return plan->compiled_math.nodes[child].kind == CompiledMathNodeKind::Constant
        ? semantic::kInvalidIndex
        : compiled_math_make_root(&plan->compiled_math, child);
  };
  const auto unused_root = branch_root(unused_node);
  const auto used_root = branch_root(used_node);
  const auto node = compiled_math_intern_node(&plan->compiled_math, std::move(key));
  plan->compiled_math.nodes[node].branch_roots = {unused_root, used_root};
  return node;
}

inline semantic::Index compile_expr_value_node_raw(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const CompiledMathValueKind value_kind,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  const auto &program = plan->program;
  const auto &kernel = plan->expr_kernels[static_cast<std::size_t>(expr_id)];
  const auto constant = [&](const double value) {
    return compiled_math_constant(&plan->compiled_math, value);
  };
  const auto complement = [&](const semantic::Index child) {
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Complement,
        child,
        CompiledMathValueKind::Cdf);
  };
  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
    if (value_kind == CompiledMathValueKind::Survival) {
      return constant(1.0);
    }
    return constant(0.0);

  case semantic::ExprKind::TrueExpr:
    if (value_kind == CompiledMathValueKind::Density) {
      return constant(0.0);
    }
    return constant(1.0);

  case semantic::ExprKind::Event:
    if (value_kind == CompiledMathValueKind::Density) {
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourcePdf,
          kernel.event_source_id,
          time_id,
          source_view_id);
    }
    if (value_kind == CompiledMathValueKind::Cdf) {
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourceCdf,
          kernel.event_source_id,
          time_id,
          source_view_id);
    }
    return compile_expr_source_node(
        plan,
        CompiledMathNodeKind::SourceSurvival,
        kernel.event_source_id,
        time_id,
        source_view_id);

  case semantic::ExprKind::And:
  case semantic::ExprKind::Or:
    return compile_expr_distribution_node(
        plan,
        expr_id,
        value_kind,
        time_id,
        source_view_id);

  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    if (value_kind == CompiledMathValueKind::Cdf) {
      return complement(
          compile_expr_value_node(
              plan,
              child,
              CompiledMathValueKind::Cdf,
              time_id,
              source_view_id));
    }
    if (value_kind == CompiledMathValueKind::Survival) {
      return compile_expr_value_node(
          plan,
          child,
          CompiledMathValueKind::Cdf,
          time_id,
          source_view_id);
    }
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Negate,
        compile_expr_value_node(
            plan,
            child,
            CompiledMathValueKind::Density,
            time_id,
            source_view_id),
        CompiledMathValueKind::Density);
  }

  case semantic::ExprKind::Guard:
    return compile_expr_distribution_node(
        plan,
        expr_id,
        value_kind,
        time_id,
        source_view_id);
  }

  throw std::runtime_error(
      "exact expression compilation reached unsupported expression " +
      std::to_string(expr_id));
}

inline semantic::Index compile_expr_upper_bound_node(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const semantic::Index child_node,
    const CompiledMathValueKind value_kind,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  CompiledMathNodeKey key;
  key.kind = value_kind == CompiledMathValueKind::Density
                 ? CompiledMathNodeKind::ExprUpperBoundDensity
                 : CompiledMathNodeKind::ExprUpperBoundCdf;
  key.value_kind = value_kind;
  key.subject_id = expr_id;
  key.time_id = time_id;
  key.source_view_id = source_view_id;
  key.children.push_back(child_node);
  return compiled_math_intern_node(&plan->compiled_math, std::move(key));
}

inline bool sequence_expr_upper_bound_used(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id) {
  return expr_id != semantic::kInvalidIndex &&
         static_cast<std::size_t>(expr_id) <
             plan.sequence.expr_upper_bound_used.size() &&
         plan.sequence.expr_upper_bound_used[
             static_cast<std::size_t>(expr_id)] != 0U;
}

inline semantic::Index compile_expr_value_node(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const CompiledMathValueKind value_kind,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  if (sequence_expr_upper_bound_used(*plan, expr_id)) {
    if (value_kind == CompiledMathValueKind::Survival) {
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Complement,
          compile_expr_value_node(
              plan,
              expr_id,
              CompiledMathValueKind::Cdf,
              time_id,
              source_view_id),
          CompiledMathValueKind::Survival);
    }
    const auto raw_node =
        compile_expr_value_node_raw(
            plan, expr_id, value_kind, time_id, source_view_id);
    return compile_expr_upper_bound_node(
        plan,
        expr_id,
        raw_node,
        value_kind,
        time_id,
        source_view_id);
  }
  return compile_expr_value_node_raw(
      plan, expr_id, value_kind, time_id, source_view_id);
}

inline semantic::Index compiled_math_root_node_id(
    const CompiledMathProgram &program,
    const semantic::Index root_id) {
  if (root_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  return program.roots[static_cast<std::size_t>(root_id)].node_id;
}

} // namespace detail
} // namespace accumulatr::eval
