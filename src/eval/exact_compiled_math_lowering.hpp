#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline semantic::Index compile_relation_condition_id(
    ExactVariantBuildState *plan,
    const ExactRelationTemplate &relation_template) {
  CompiledMathConditionKey key;
  key.source_ids = relation_template.source_ids;
  key.relations.reserve(relation_template.relations.size());
  for (const auto relation : relation_template.relations) {
    key.relations.push_back(static_cast<std::uint8_t>(relation));
  }
  return compiled_math_intern_condition(&plan->compiled_math, std::move(key));
}

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
    CompiledMathNodeKind value_kind,
    semantic::Index condition_id,
    semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    semantic::Index source_view_id = 0);

inline semantic::Index compile_expr_distribution_node(
    ExactVariantBuildState *plan,
    semantic::Index expr_id,
    CompiledMathNodeKind value_kind,
    semantic::Index condition_id,
    semantic::Index time_id,
    semantic::Index source_view_id);

inline semantic::Index compile_expr_source_node(
    ExactVariantBuildState *plan,
    const CompiledMathNodeKind kind,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0) {
  return compiled_math_source_node(
      &plan->compiled_math,
      kind,
      source_id,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_expr_child_value_node(
    ExactVariantBuildState *plan,
    const ExactIndexSpan children,
    const semantic::Index index,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0) {
  const auto &program = plan->lowered.program;
  return compile_expr_value_node(
      plan,
      program.expr_args[static_cast<std::size_t>(children.offset + index)],
      value_kind,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_expr_unsupported_node(
    ExactVariantBuildState *plan,
    const CompiledMathNodeKind kind,
    const semantic::Index expr_id,
    const semantic::Index condition_id) {
  (void)plan;
  (void)kind;
  (void)condition_id;
  throw std::runtime_error(
      "exact expression compilation reached unsupported expression " +
      std::to_string(expr_id) +
      "; runtime expression interpretation is disabled");
}

inline semantic::Index compile_guard_unless_allowed_node(
    ExactVariantBuildState *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  const auto &program = plan->lowered.program;
  std::vector<semantic::Index> blocked_factors;
  blocked_factors.push_back(
      compile_expr_value_node(
          plan,
          kernel.guard_blocker_expr_id,
          CompiledMathNodeKind::ExprCdf,
          condition_id,
          time_id,
          source_view_id));
  bool any_unless_true = false;
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    const auto child = program.expr_args[
        static_cast<std::size_t>(kernel.children.offset + i)];
    const auto child_kind = static_cast<semantic::ExprKind>(
        program.expr_kind[static_cast<std::size_t>(child)]);
    if (child_kind == semantic::ExprKind::TrueExpr) {
      any_unless_true = true;
      continue;
    }
    if (child_kind == semantic::ExprKind::Impossible) {
      continue;
    }
    blocked_factors.push_back(
        compile_expr_value_node(
            plan,
            child,
            CompiledMathNodeKind::ExprSurvival,
            condition_id,
            time_id,
            source_view_id));
  }
  if (any_unless_true) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  const auto blocked_node =
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Product,
          std::move(blocked_factors),
          CompiledMathValueKind::Cdf);
  return compiled_math_unary_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Complement,
      blocked_node,
      CompiledMathValueKind::Cdf);
}

inline semantic::Index compile_guard_unless_density_node(
    ExactVariantBuildState *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  const auto allowed_node =
      compile_guard_unless_allowed_node(
          plan,
          kernel,
          condition_id,
          time_id,
          source_view_id);
  const auto ref_density =
      compile_expr_value_node(
          plan,
          kernel.guard_ref_expr_id,
          CompiledMathNodeKind::ExprDensity,
          condition_id,
          time_id,
          source_view_id);
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Product,
      std::vector<semantic::Index>{ref_density, allowed_node},
      CompiledMathValueKind::Density);
}

inline bool compiled_condition_has_source_order(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      before_source_id == semantic::kInvalidIndex ||
      after_source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]);
    if (kind == CompiledMathConditionFactKind::SourceOrder &&
        condition.fact_subject_ids[i] == before_source_id &&
        condition.fact_aux_ids[i] == after_source_id) {
      return true;
    }
  }
  return false;
}

inline semantic::Index compiled_condition_source_exact_time_id(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return semantic::kInvalidIndex;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            CompiledMathConditionFactKind::SourceExact &&
        condition.fact_subject_ids[i] == source_id) {
      return condition.fact_time_ids[i];
    }
  }
  return semantic::kInvalidIndex;
}

inline bool compiled_condition_has_expr_upper_bound(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index expr_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]);
    if (kind == CompiledMathConditionFactKind::ExprUpperBound &&
        condition.fact_subject_ids[i] == expr_id) {
      return true;
    }
  }
  return false;
}

inline CompiledMathIndexSpan compiled_condition_expr_upper_bound_fact_span(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index expr_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      expr_id == semantic::kInvalidIndex) {
    return {};
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return {};
  }
  const auto &condition = program.conditions[condition_pos];
  const auto expr_pos = static_cast<std::size_t>(expr_id);
  if (expr_pos >= condition.expr_upper_fact_spans.size()) {
    return {};
  }
  return condition.expr_upper_fact_spans[expr_pos];
}

inline CompiledMathIndexSpan compile_timed_upper_bound_terms(
    CompiledMathProgram *program,
    const semantic::Index condition_id,
    const std::vector<semantic::Index> &fact_indices,
    const CompiledMathIndexSpan fact_span) {
  if (condition_id == 0 ||
      condition_id == semantic::kInvalidIndex ||
      fact_span.empty()) {
    return {};
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program->conditions.size()) {
    return {};
  }
  const auto &condition = program->conditions[condition_pos];
  const auto offset =
      static_cast<semantic::Index>(
          program->timed_upper_bound_terms.size());
  for (semantic::Index i = 0; i < fact_span.size; ++i) {
    const auto fact_pos = static_cast<std::size_t>(
        fact_indices[static_cast<std::size_t>(fact_span.offset + i)]);
    const auto time_id =
        fact_pos < condition.fact_time_ids.size()
            ? condition.fact_time_ids[fact_pos]
            : semantic::kInvalidIndex;
    const auto normalizer_node_id =
        fact_pos < condition.fact_normalizer_node_ids.size()
            ? condition.fact_normalizer_node_ids[fact_pos]
            : semantic::kInvalidIndex;
    program->timed_upper_bound_terms.push_back(
        CompiledMathTimedUpperBoundTerm{time_id, normalizer_node_id});
  }
  return CompiledMathIndexSpan{
      offset,
      static_cast<semantic::Index>(
          program->timed_upper_bound_terms.size() -
          static_cast<std::size_t>(offset))};
}

inline bool compiled_condition_has_source_relation(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id,
    const ExactRelation relation) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.source_ids.size(); ++i) {
    if (condition.source_ids[i] == source_id &&
        static_cast<ExactRelation>(condition.relations[i]) == relation) {
      return true;
    }
  }
  return false;
}

inline bool compiled_source_view_knows_before(
    const ExactVariantBuildState &plan,
    const semantic::Index source_view_id,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  if (source_view_id == 0 ||
      source_view_id == semantic::kInvalidIndex ||
      before_source_id == semantic::kInvalidIndex ||
      after_source_id == semantic::kInvalidIndex ||
      before_source_id == after_source_id) {
    return false;
  }
  const auto pos = static_cast<std::size_t>(source_view_id - 1U);
  if (pos >= plan.compiled_source_views.size()) {
    return false;
  }
  const auto &source_view = plan.compiled_source_views[pos];
  auto relation_for = [&](const semantic::Index source_id) {
    for (std::size_t i = 0; i < source_view.source_ids.size(); ++i) {
      if (source_view.source_ids[i] == source_id) {
        return source_view.relations[i];
      }
    }
    return ExactRelation::Unknown;
  };
  const auto before_relation = relation_for(before_source_id);
  const auto after_relation = relation_for(after_source_id);
  return before_relation != ExactRelation::Unknown &&
         after_relation != ExactRelation::Unknown &&
         before_relation < after_relation;
}

inline bool compiled_guard_order_blocks(
    const ExactVariantBuildState &plan,
    const semantic::Index condition_id,
    const semantic::Index source_view_id,
    const semantic::Index blocker_source_id,
    const semantic::Index ref_source_id) {
  return compiled_condition_has_source_order(
             plan.compiled_math,
             condition_id,
             blocker_source_id,
             ref_source_id) ||
         compiled_source_view_knows_before(
             plan,
             source_view_id,
             blocker_source_id,
             ref_source_id);
}

inline bool compiled_condition_has_source_fact(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const CompiledMathConditionFactKind fact_kind,
    const semantic::Index source_id,
    const CompiledMathTimeSlot time_slot) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  const auto target_time = static_cast<semantic::Index>(time_slot);
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            fact_kind &&
        condition.fact_subject_ids[i] == source_id &&
        condition.fact_time_ids[i] == target_time) {
      return true;
    }
  }
  return false;
}

inline bool compiled_condition_forces_source_after_observed(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  if (compiled_condition_has_source_relation(
          program, condition_id, source_id, ExactRelation::After)) {
    return true;
  }
  return compiled_condition_has_source_fact(
             program,
             condition_id,
             CompiledMathConditionFactKind::SourceLowerBound,
             source_id,
             CompiledMathTimeSlot::Observed) ||
         compiled_condition_has_source_fact(
             program,
             condition_id,
             CompiledMathConditionFactKind::SourceLowerBound,
             source_id,
             CompiledMathTimeSlot::Readiness) ||
         compiled_condition_has_source_fact(
             program,
             condition_id,
             CompiledMathConditionFactKind::SourceLowerBound,
             source_id,
             CompiledMathTimeSlot::Active);
}

inline bool compiled_condition_forces_source_certain(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  return compiled_condition_has_source_relation(
             program, condition_id, source_id, ExactRelation::Before) ||
         compiled_condition_has_source_relation(
             program, condition_id, source_id, ExactRelation::At);
}

inline semantic::Index compile_integral_zero_to_current_node(
    ExactVariantBuildState *plan,
    const semantic::Index integrand_node,
    const semantic::Index condition_id,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  const auto integrand_root =
      compiled_math_make_root(&plan->compiled_math, integrand_node);
  return compiled_math_integral_zero_to_current_node(
      &plan->compiled_math,
      integrand_root,
      condition_id,
      time_id,
      source_view_id,
      bind_time_id);
}

inline semantic::Index compile_signed_term_node(
    ExactVariantBuildState *plan,
    const semantic::Index node_id,
    const int sign) {
  if (sign >= 0) {
    return node_id;
  }
  return compiled_math_unary_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Negate,
      node_id);
}

inline semantic::Index compile_outcome_subset_unused_node(
    ExactVariantBuildState *plan,
    const std::vector<semantic::Index> &outcome_indices,
    const bool used = false) {
  if (outcome_indices.empty()) {
    return compiled_math_constant(&plan->compiled_math, used ? 0.0 : 1.0);
  }
  const auto offset =
      static_cast<semantic::Index>(
          plan->compiled_outcome_gate_indices.size());
  plan->compiled_outcome_gate_indices.insert(
      plan->compiled_outcome_gate_indices.end(),
      outcome_indices.begin(),
      outcome_indices.end());
  CompiledMathNodeKey key;
  key.kind = used ? CompiledMathNodeKind::OutcomeSubsetUsed
                  : CompiledMathNodeKind::OutcomeSubsetUnused;
  key.value_kind = CompiledMathValueKind::Scalar;
  key.subject_id = offset;
  key.aux_id = static_cast<semantic::Index>(outcome_indices.size());
  return compiled_math_intern_node(&plan->compiled_math, std::move(key));
}

inline semantic::Index compile_expr_value_node(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id);

inline semantic::Index compile_expr_value_node_raw(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  const auto &program = plan->lowered.program;
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
  const auto unsupported = [&]() {
    return compile_expr_unsupported_node(plan, value_kind, expr_id, condition_id);
  };

  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return constant(1.0);
    }
    return constant(0.0);

  case semantic::ExprKind::TrueExpr:
    if (value_kind == CompiledMathNodeKind::ExprDensity) {
      return constant(0.0);
    }
    return constant(1.0);

  case semantic::ExprKind::Event:
    if (value_kind == CompiledMathNodeKind::ExprDensity) {
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourcePdf,
          kernel.event_source_id,
          condition_id,
          time_id,
          source_view_id);
    }
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourceCdf,
          kernel.event_source_id,
          condition_id,
          time_id,
          source_view_id);
    }
    return compile_expr_source_node(
        plan,
        CompiledMathNodeKind::SourceSurvival,
        kernel.event_source_id,
        condition_id,
        time_id,
        source_view_id);

  case semantic::ExprKind::And: {
    return compile_expr_distribution_node(
        plan,
        expr_id,
        value_kind,
        condition_id,
        time_id,
        source_view_id);
  }
  case semantic::ExprKind::Or: {
    return compile_expr_distribution_node(
        plan,
        expr_id,
        value_kind,
        condition_id,
        time_id,
        source_view_id);
  }

  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return complement(
          compile_expr_value_node(
              plan,
              child,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return compile_expr_value_node(
          plan,
          child,
          CompiledMathNodeKind::ExprCdf,
          condition_id,
          time_id,
          source_view_id);
    }
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Negate,
        compile_expr_value_node(
            plan,
            child,
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id),
        CompiledMathValueKind::Density);
  }

  case semantic::ExprKind::Guard:
    if (kernel.has_unless) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        return compile_guard_unless_density_node(
            plan,
            kernel,
            condition_id,
            time_id,
            source_view_id);
      } else if (value_kind == CompiledMathNodeKind::ExprCdf) {
        const auto bind_time_id =
            time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
                ? time_id
                : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
        const auto density_node =
            compile_guard_unless_density_node(
                plan,
                kernel,
                condition_id,
                bind_time_id,
                source_view_id);
        return compile_integral_zero_to_current_node(
            plan,
            density_node,
            condition_id,
            time_id,
            source_view_id,
            bind_time_id);
      } else {
        return compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Complement,
            compile_expr_value_node_raw(
                plan,
                expr_id,
                CompiledMathNodeKind::ExprCdf,
                condition_id,
                time_id,
                source_view_id),
            CompiledMathValueKind::Survival);
      }
    }
    if (!kernel.has_unless) {
      return compile_expr_distribution_node(
          plan,
          expr_id,
          value_kind,
          condition_id,
          time_id,
          source_view_id);
    }
    return unsupported();
  }

  return unsupported();
}

inline semantic::Index compile_expr_upper_bound_node(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const semantic::Index child_node,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  CompiledMathNodeKey key;
  key.kind = value_kind == CompiledMathNodeKind::ExprDensity
                 ? CompiledMathNodeKind::ExprUpperBoundDensity
                 : CompiledMathNodeKind::ExprUpperBoundCdf;
  key.value_kind = value_kind == CompiledMathNodeKind::ExprDensity
                       ? CompiledMathValueKind::Density
                       : CompiledMathValueKind::Cdf;
  key.subject_id = expr_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.source_view_id = source_view_id;
  const auto upper_span =
      compiled_condition_expr_upper_bound_fact_span(
          plan->compiled_math, condition_id, expr_id);
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  const auto term_span =
      condition_id == 0 ||
              condition_id == semantic::kInvalidIndex ||
              condition_pos >= plan->compiled_math.conditions.size()
          ? CompiledMathIndexSpan{}
          : compile_timed_upper_bound_terms(
                &plan->compiled_math,
                condition_id,
                plan->compiled_math.conditions[condition_pos]
                    .expr_upper_fact_indices,
                upper_span);
  key.aux_id =
      term_span.empty() ? semantic::kInvalidIndex : term_span.offset;
  key.aux2_id = term_span.size;
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
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  if (compiled_condition_has_expr_upper_bound(
          plan->compiled_math,
          condition_id,
          expr_id) ||
      sequence_expr_upper_bound_used(*plan, expr_id)) {
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Complement,
          compile_expr_value_node(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id),
          CompiledMathValueKind::Survival);
    }
    const auto raw_node =
        compile_expr_value_node_raw(
            plan, expr_id, value_kind, condition_id, time_id, source_view_id);
    if (value_kind == CompiledMathNodeKind::ExprDensity ||
        value_kind == CompiledMathNodeKind::ExprCdf) {
      return compile_expr_upper_bound_node(
          plan,
          expr_id,
          raw_node,
          value_kind,
          condition_id,
          time_id,
          source_view_id);
    }
    return raw_node;
  }
  return compile_expr_value_node_raw(
      plan, expr_id, value_kind, condition_id, time_id, source_view_id);
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
