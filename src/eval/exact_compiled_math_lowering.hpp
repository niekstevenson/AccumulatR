#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline semantic::Index compile_relation_condition_id(
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
    semantic::Index expr_id,
    CompiledMathNodeKind value_kind,
    semantic::Index condition_id,
    semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    semantic::Index source_view_id = 0);

inline semantic::Index compile_expr_source_node(
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
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

inline bool compiled_condition_has_union_conditioning(
    const CompiledMathProgram &program,
    const semantic::Index condition_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (const auto kind_value : condition.fact_kinds) {
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(kind_value);
    if (kind == CompiledMathConditionFactKind::SourceExact ||
        kind == CompiledMathConditionFactKind::SourceUpperBound ||
        kind == CompiledMathConditionFactKind::SourceLowerBound ||
        kind == CompiledMathConditionFactKind::ExprUpperBound ||
        kind == CompiledMathConditionFactKind::SourceOrder ||
        kind == CompiledMathConditionFactKind::GuardUpperBound) {
      return true;
    }
  }
  return false;
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

inline CompiledMathIndexSpan compiled_condition_guard_upper_bound_fact_span(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      ref_source_id == semantic::kInvalidIndex ||
      blocker_source_id == semantic::kInvalidIndex) {
    return {};
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return {};
  }
  const auto &lookup =
      program.conditions[condition_pos].guard_upper_fact_lookup;
  const auto found = std::lower_bound(
      lookup.entries.begin(),
      lookup.entries.end(),
      std::pair<semantic::Index, semantic::Index>{
          ref_source_id,
          blocker_source_id},
      [](const CompiledMathPairFactEntry &entry, const auto &target) {
        return entry.first < target.first ||
               (entry.first == target.first && entry.second < target.second);
      });
  if (found == lookup.entries.end() ||
      found->first != ref_source_id ||
      found->second != blocker_source_id) {
    return {};
  }
  return found->facts;
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
    const ExactVariantPlan &plan,
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
    const ExactVariantPlan &plan,
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

inline bool expr_condition_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index source_id,
    const semantic::Index condition_id);

inline bool compiled_condition_has_guard_upper_bound_expr(
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
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            CompiledMathConditionFactKind::GuardUpperBound &&
        condition.fact_subject_ids[i] == expr_id) {
      return true;
    }
  }
  return false;
}

inline bool compiled_condition_has_guard_upper_bound_fact(
    const CompiledMathProgram &program,
    const semantic::Index condition_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (const auto kind_value : condition.fact_kinds) {
    if (static_cast<CompiledMathConditionFactKind>(kind_value) ==
        CompiledMathConditionFactKind::GuardUpperBound) {
      return true;
    }
  }
  return false;
}

inline bool compiled_condition_has_guard_upper_bound_signature(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      ref_source_id == semantic::kInvalidIndex ||
      blocker_source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            CompiledMathConditionFactKind::GuardUpperBound &&
        condition.fact_aux_ids[i] == ref_source_id &&
        condition.fact_aux2_ids[i] == blocker_source_id) {
      return true;
    }
  }
  return false;
}

inline bool expr_condition_certain(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  switch (kind) {
  case semantic::ExprKind::TrueExpr:
    return true;
  case semantic::ExprKind::Event:
    return compiled_condition_forces_source_certain(
        plan.compiled_math,
        condition_id,
        program.expr_source_ids[pos]);
  case semantic::ExprKind::And: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    if (begin == end) {
      return false;
    }
    for (semantic::Index i = begin; i < end; ++i) {
      if (!expr_condition_certain(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return false;
      }
    }
    return true;
  }
  case semantic::ExprKind::Or: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    for (semantic::Index i = begin; i < end; ++i) {
      if (expr_condition_certain(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return true;
      }
    }
    return false;
  }
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::Guard:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool expr_condition_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index source_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind == semantic::ExprKind::Event) {
    return program.expr_source_ids[pos] == source_id;
  }
  if (kind == semantic::ExprKind::Guard) {
    const auto &kernel = plan.expr_kernels[pos];
    return kernel.simple_event_guard &&
           !kernel.has_unless &&
           kernel.guard_ref_source_id == source_id &&
           compiled_condition_forces_source_after_observed(
               plan.compiled_math,
               condition_id,
               kernel.guard_blocker_source_id);
  }
  if (kind != semantic::ExprKind::And) {
    return false;
  }
  bool found_candidate = false;
  const auto begin = program.expr_arg_offsets[pos];
  const auto end = program.expr_arg_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    const auto child = program.expr_args[static_cast<std::size_t>(i)];
    if (expr_condition_equiv_event_source(
            plan, child, source_id, condition_id)) {
      if (found_candidate) {
        return false;
      }
      found_candidate = true;
      continue;
    }
    if (!expr_condition_certain(plan, child, condition_id)) {
      return false;
    }
  }
  return found_candidate;
}

inline bool expr_condition_certain_at_observed_cdf(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  switch (kind) {
  case semantic::ExprKind::TrueExpr:
    return true;
  case semantic::ExprKind::Event: {
    const auto source_id = program.expr_source_ids[pos];
    return compiled_condition_forces_source_certain(
               plan.compiled_math,
               condition_id,
               source_id) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceExact,
               source_id,
               CompiledMathTimeSlot::Observed) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceExact,
               source_id,
               CompiledMathTimeSlot::Readiness) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceExact,
               source_id,
               CompiledMathTimeSlot::Active) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceUpperBound,
               source_id,
               CompiledMathTimeSlot::Observed) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceUpperBound,
               source_id,
               CompiledMathTimeSlot::Readiness) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceUpperBound,
               source_id,
               CompiledMathTimeSlot::Active);
    }
  case semantic::ExprKind::Guard: {
    const auto &kernel = plan.expr_kernels[pos];
    return compiled_condition_has_guard_upper_bound_expr(
               plan.compiled_math,
               condition_id,
               expr_id) ||
           (kernel.simple_event_guard &&
            !kernel.has_unless &&
            compiled_condition_has_guard_upper_bound_signature(
                plan.compiled_math,
                condition_id,
                kernel.guard_ref_source_id,
                kernel.guard_blocker_source_id));
  }
  case semantic::ExprKind::And: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    if (begin == end) {
      return false;
    }
    for (semantic::Index i = begin; i < end; ++i) {
      if (!expr_condition_certain_at_observed_cdf(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return false;
      }
    }
    return true;
  }
  case semantic::ExprKind::Or: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    for (semantic::Index i = begin; i < end; ++i) {
      if (expr_condition_certain_at_observed_cdf(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return true;
      }
    }
    return false;
  }
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool expr_condition_equiv_event_source_at_observed_cdf(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index source_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind == semantic::ExprKind::Event) {
    return program.expr_source_ids[pos] == source_id;
  }
  if (kind == semantic::ExprKind::Guard) {
    return expr_condition_equiv_event_source(
        plan, expr_id, source_id, condition_id);
  }
  if (kind != semantic::ExprKind::And) {
    return false;
  }
  bool found_candidate = false;
  const auto begin = program.expr_arg_offsets[pos];
  const auto end = program.expr_arg_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    const auto child = program.expr_args[static_cast<std::size_t>(i)];
    if (expr_condition_equiv_event_source_at_observed_cdf(
            plan, child, source_id, condition_id)) {
      if (found_candidate) {
        return false;
      }
      found_candidate = true;
      continue;
    }
    if (!expr_condition_certain_at_observed_cdf(
            plan, child, condition_id)) {
      return false;
    }
  }
  return found_candidate;
}

inline bool compiled_condition_absorbs_union_to_source(
    const ExactVariantPlan &plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const bool observed_cdf_context,
    semantic::Index *source_id) {
  for (semantic::Index i = 0;
       i < kernel.union_absorption_candidate_span.size;
       ++i) {
    const auto candidate = plan.expr_union_absorption_sources[
        static_cast<std::size_t>(
            kernel.union_absorption_candidate_span.offset + i)];
    for (semantic::Index j = 0; j < kernel.children.size; ++j) {
      const auto child = plan.lowered.program.expr_args[
          static_cast<std::size_t>(kernel.children.offset + j)];
      const bool child_equiv =
          observed_cdf_context
              ? expr_condition_equiv_event_source_at_observed_cdf(
                    plan, child, candidate, condition_id)
              : expr_condition_equiv_event_source(
                    plan, child, candidate, condition_id);
      if (child_equiv) {
        if (source_id != nullptr) {
          *source_id = candidate;
        }
        return true;
      }
    }
  }
  return false;
}

inline semantic::Index compile_integral_zero_to_current_node(
    ExactVariantPlan *plan,
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

inline semantic::Index compile_simple_guard_raw_density_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  std::vector<semantic::Index> factors;
  factors.reserve(2U);
  factors.push_back(
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourcePdf,
          kernel.guard_ref_source_id,
          condition_id,
          time_id,
          source_view_id));
  factors.push_back(
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourceSurvival,
          kernel.guard_blocker_source_id,
          condition_id,
          time_id,
          source_view_id));
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Product,
      std::move(factors),
      CompiledMathValueKind::Density);
}

inline semantic::Index compile_simple_guard_upper_wrapper_node(
    ExactVariantPlan *plan,
    const CompiledMathNodeKind kind,
    const semantic::Index expr_id,
    const semantic::Index raw_node,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  CompiledMathNodeKey key;
  key.kind = kind;
  key.value_kind = kind == CompiledMathNodeKind::SimpleGuardDensity
                       ? CompiledMathValueKind::Density
                       : CompiledMathValueKind::Cdf;
  key.subject_id = expr_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.source_view_id = source_view_id;
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan->expr_kernels.size()) {
    throw std::runtime_error(
        "compiled simple guard upper wrapper has no expression kernel");
  }
  const auto &kernel = plan->expr_kernels[static_cast<std::size_t>(expr_id)];
  const auto upper_span =
      compiled_condition_guard_upper_bound_fact_span(
          plan->compiled_math,
          condition_id,
          kernel.guard_ref_source_id,
          kernel.guard_blocker_source_id);
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (upper_span.empty() ||
      condition_id == 0 ||
      condition_id == semantic::kInvalidIndex ||
      condition_pos >= plan->compiled_math.conditions.size()) {
    throw std::runtime_error(
        "compiled simple guard upper wrapper has no planned fact span");
  }
  const auto term_span =
      compile_timed_upper_bound_terms(
          &plan->compiled_math,
          condition_id,
          plan->compiled_math.conditions[condition_pos]
              .guard_upper_fact_lookup.fact_indices,
          upper_span);
  key.aux_id =
      term_span.empty() ? semantic::kInvalidIndex : term_span.offset;
  key.aux2_id = term_span.size;
  key.children.push_back(raw_node);
  return compiled_math_intern_node(&plan->compiled_math, std::move(key));
}

inline bool compile_simple_guard_needs_upper_wrapper(
    const ExactVariantPlan &plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id) {
  return compiled_condition_has_guard_upper_bound_signature(
             plan.compiled_math,
             condition_id,
             kernel.guard_ref_source_id,
             kernel.guard_blocker_source_id);
}

inline semantic::Index compile_simple_guard_density_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (compiled_guard_order_blocks(
          *plan,
          condition_id,
          source_view_id,
          kernel.guard_blocker_source_id,
          kernel.guard_ref_source_id)) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  const auto raw_node =
      compile_simple_guard_raw_density_node(
          plan, kernel, condition_id, time_id, source_view_id);
  if (!compile_simple_guard_needs_upper_wrapper(
          *plan, kernel, condition_id)) {
    return raw_node;
  }
  return compile_simple_guard_upper_wrapper_node(
      plan,
      CompiledMathNodeKind::SimpleGuardDensity,
      expr_id,
      raw_node,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_simple_guard_raw_cdf_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (compiled_guard_order_blocks(
          *plan,
          condition_id,
          source_view_id,
          kernel.guard_blocker_source_id,
          kernel.guard_ref_source_id)) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }

  const auto ref_exact_time_id =
      compiled_condition_source_exact_time_id(
          plan->compiled_math,
          condition_id,
          kernel.guard_ref_source_id);
  if (ref_exact_time_id != semantic::kInvalidIndex) {
    const auto blocker_survival_at_ref =
        compiled_math_source_node(
            &plan->compiled_math,
            CompiledMathNodeKind::SourceSurvival,
            kernel.guard_blocker_source_id,
            condition_id,
            ref_exact_time_id,
            source_view_id);
    return compiled_math_time_gate_node(
        &plan->compiled_math,
        blocker_survival_at_ref,
        time_id,
        ref_exact_time_id,
        CompiledMathValueKind::Cdf);
  }

  const auto blocker_exact_time_id =
      compiled_condition_source_exact_time_id(
          plan->compiled_math,
          condition_id,
          kernel.guard_blocker_source_id);
  if (blocker_exact_time_id != semantic::kInvalidIndex) {
    return compiled_math_source_node(
        &plan->compiled_math,
        CompiledMathNodeKind::SourceCdf,
        kernel.guard_ref_source_id,
        condition_id,
        time_id,
        source_view_id,
        blocker_exact_time_id);
  }

  (void)expr_id;
  const auto bind_time_id =
      time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
          ? time_id
          : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
  const auto raw_density_node =
      compile_simple_guard_raw_density_node(
          plan, kernel, condition_id, bind_time_id, source_view_id);
  return compile_integral_zero_to_current_node(
      plan,
      raw_density_node,
      condition_id,
      time_id,
      source_view_id,
      bind_time_id);
}

inline semantic::Index compile_simple_guard_cdf_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  const auto raw_node =
      compile_simple_guard_raw_cdf_node(
          plan, expr_id, kernel, condition_id, time_id, source_view_id);
  if (!compile_simple_guard_needs_upper_wrapper(
          *plan, kernel, condition_id)) {
    return raw_node;
  }
  return compile_simple_guard_upper_wrapper_node(
      plan,
      CompiledMathNodeKind::SimpleGuardCdf,
      expr_id,
      raw_node,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_signed_term_node(
    ExactVariantPlan *plan,
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
    ExactVariantPlan *plan,
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

inline semantic::Index compile_expr_span_conjunction_density_node(
    ExactVariantPlan *plan,
    const std::vector<semantic::Index> &arena,
    const ExactIndexSpan span,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (span.empty()) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  if (span.size == 1U) {
    return compile_expr_value_node(
        plan,
        arena[static_cast<std::size_t>(span.offset)],
        CompiledMathNodeKind::ExprDensity,
        condition_id,
        time_id,
        source_view_id);
  }

  std::vector<semantic::Index> terms;
  terms.reserve(static_cast<std::size_t>(span.size));
  for (semantic::Index i = 0; i < span.size; ++i) {
    std::vector<semantic::Index> factors;
    factors.reserve(static_cast<std::size_t>(span.size));
    factors.push_back(
        compile_expr_value_node(
            plan,
            arena[static_cast<std::size_t>(span.offset + i)],
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id));
    for (semantic::Index j = 0; j < span.size; ++j) {
      if (i == j) {
        continue;
      }
      factors.push_back(
          compile_expr_value_node(
              plan,
              arena[static_cast<std::size_t>(span.offset + j)],
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    terms.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(factors)));
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::CleanSignedSum,
      std::move(terms),
      CompiledMathValueKind::Density);
}

inline bool exact_union_has_multi_subset(const ExactVariantPlan &plan,
                                         const ExactExprKernel &kernel) {
  for (semantic::Index i = 0; i < kernel.union_subset_span.size; ++i) {
    const auto &subset = plan.expr_union_subsets[
        static_cast<std::size_t>(kernel.union_subset_span.offset + i)];
    if (subset.child_positions.size > 1U) {
      return true;
    }
  }
  return false;
}

inline semantic::Index compile_union_subset_density_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const ExactExprUnionSubset &subset,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (subset.child_positions.empty()) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  if (subset.child_positions.size == 1U) {
    const auto active_pos =
        plan->expr_union_subset_child_positions[
            static_cast<std::size_t>(subset.child_positions.offset)];
    return compile_expr_child_value_node(
        plan,
        kernel.children,
        active_pos,
        CompiledMathNodeKind::ExprDensity,
        condition_id,
        time_id,
        source_view_id);
  }

  std::vector<semantic::Index> terms;
  terms.reserve(static_cast<std::size_t>(subset.child_positions.size));
  for (semantic::Index i = 0; i < subset.child_positions.size; ++i) {
    const auto active_pos =
        plan->expr_union_subset_child_positions[
            static_cast<std::size_t>(subset.child_positions.offset + i)];
    std::vector<semantic::Index> factors;
    factors.reserve(static_cast<std::size_t>(subset.child_positions.size));
    factors.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            active_pos,
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id));
    for (semantic::Index j = 0; j < subset.child_positions.size; ++j) {
      if (i == j) {
        continue;
      }
      const auto other_pos =
          plan->expr_union_subset_child_positions[
              static_cast<std::size_t>(subset.child_positions.offset + j)];
      factors.push_back(
          compile_expr_child_value_node(
              plan,
              kernel.children,
              other_pos,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    terms.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(factors),
            CompiledMathValueKind::Density));
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::CleanSignedSum,
      std::move(terms),
      CompiledMathValueKind::Density);
}

inline semantic::Index compile_union_density_sum_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id,
    const bool multi_subset_only) {
  std::vector<semantic::Index> terms;
  for (semantic::Index i = 0; i < kernel.union_subset_span.size; ++i) {
    const auto &subset =
        plan->expr_union_subsets[
            static_cast<std::size_t>(kernel.union_subset_span.offset + i)];
    if (multi_subset_only && subset.child_positions.size <= 1U) {
      continue;
    }
    const auto subset_node =
        compile_union_subset_density_node(
            plan,
            kernel,
            subset,
            condition_id,
            time_id,
            source_view_id);
    terms.push_back(
        compile_signed_term_node(plan, subset_node, subset.sign));
  }
  if (terms.empty()) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::CleanSignedSum,
      std::move(terms),
      CompiledMathValueKind::Density);
}

inline std::vector<semantic::Index> compile_union_child_density_cdf_nodes(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  std::vector<semantic::Index> children;
  children.reserve(static_cast<std::size_t>(kernel.children.size) * 2U);
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    children.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            i,
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id));
  }
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    children.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            i,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id));
  }
  return children;
}

inline semantic::Index compile_union_multi_subset_density_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  const auto density_time_id =
      bind_time_id == semantic::kInvalidIndex ? time_id : bind_time_id;
  (void)expr_id;
  return compile_union_density_sum_node(
      plan,
      kernel,
      condition_id,
      density_time_id,
      source_view_id,
      true);
}

inline semantic::Index compile_overlapping_union_kernel_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (value_kind == CompiledMathNodeKind::ExprSurvival) {
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Complement,
        compile_overlapping_union_kernel_node(
            plan,
            expr_id,
            kernel,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id),
        CompiledMathValueKind::Survival);
  }

  auto *program = &plan->compiled_math;
  if (value_kind == CompiledMathNodeKind::ExprDensity) {
    (void)expr_id;
    return compile_union_density_sum_node(
        plan,
        kernel,
        condition_id,
        time_id,
        source_view_id,
        false);
  }

  std::vector<semantic::Index> children;
  children.reserve(static_cast<std::size_t>(kernel.children.size) + 1U);
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    children.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            i,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id));
  }
  if (exact_union_has_multi_subset(*plan, kernel)) {
    const auto bind_time_id =
        time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
            ? time_id
            : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
    const auto multi_density_node =
        compile_union_multi_subset_density_node(
            plan,
            expr_id,
            kernel,
            condition_id,
            time_id,
            source_view_id,
            bind_time_id);
    const auto multi_density_root =
        compiled_math_make_root(program, multi_density_node);
    children.push_back(
        compiled_math_union_kernel_node(
            program,
            CompiledMathNodeKind::UnionKernelMultiSubsetCdf,
            expr_id,
            condition_id,
            {},
            multi_density_root,
            time_id,
            source_view_id,
            bind_time_id));
  }
  return compiled_math_union_kernel_node(
      program,
      CompiledMathNodeKind::UnionKernelCdf,
      expr_id,
      condition_id,
      std::move(children),
      semantic::kInvalidIndex,
      time_id,
      source_view_id);
}

inline semantic::Index compile_overlapping_union_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (value_kind == CompiledMathNodeKind::ExprSurvival) {
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Complement,
        compile_overlapping_union_node(
            plan,
            expr_id,
            kernel,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id),
        CompiledMathValueKind::Survival);
  }

  if (!kernel.union_absorption_candidate_span.empty() &&
      compiled_condition_has_union_conditioning(plan->compiled_math, condition_id)) {
    semantic::Index absorbed_source_id{semantic::kInvalidIndex};
    const bool absorbed = compiled_condition_absorbs_union_to_source(
        *plan,
        kernel,
        condition_id,
        value_kind == CompiledMathNodeKind::ExprCdf,
        &absorbed_source_id);
    if (absorbed) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        return compile_expr_source_node(
            plan,
            CompiledMathNodeKind::SourcePdf,
            absorbed_source_id,
            condition_id,
            time_id,
            source_view_id);
      }
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourceCdf,
          absorbed_source_id,
          condition_id,
          time_id,
          source_view_id);
    }
  }

  return compile_overlapping_union_kernel_node(
      plan, expr_id, kernel, value_kind, condition_id, time_id, source_view_id);
}

inline semantic::Index compile_expr_value_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id);

inline semantic::Index compile_expr_value_node_raw(
    ExactVariantPlan *plan,
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
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      std::vector<semantic::Index> children;
      children.reserve(static_cast<std::size_t>(kernel.children.size));
      for (semantic::Index i = 0; i < kernel.children.size; ++i) {
        children.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                i,
                CompiledMathNodeKind::ExprCdf,
                condition_id,
                time_id,
                source_view_id));
      }
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::ClampProbability,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(children)),
          CompiledMathValueKind::Cdf);
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return complement(
          compile_expr_value_node_raw(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    std::vector<semantic::Index> terms;
    terms.reserve(static_cast<std::size_t>(kernel.children.size));
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      std::vector<semantic::Index> factors;
      factors.reserve(static_cast<std::size_t>(kernel.children.size));
      factors.push_back(
          compile_expr_child_value_node(
              plan,
              kernel.children,
              i,
              CompiledMathNodeKind::ExprDensity,
              condition_id,
              time_id,
              source_view_id));
      for (semantic::Index j = 0; j < kernel.children.size; ++j) {
        if (i == j) {
          continue;
        }
        factors.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                j,
                CompiledMathNodeKind::ExprCdf,
                condition_id,
                time_id,
                source_view_id));
      }
      terms.push_back(
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(factors)));
    }
    return compiled_math_algebra_node(
        &plan->compiled_math,
        CompiledMathNodeKind::CleanSignedSum,
        std::move(terms),
        CompiledMathValueKind::Density);
  }

  case semantic::ExprKind::Or: {
    if (kernel.children_overlap) {
	      return compile_overlapping_union_node(
	          plan,
	          expr_id,
	          kernel,
	          value_kind,
	          condition_id,
	          time_id,
	          source_view_id);
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      std::vector<semantic::Index> children;
      children.reserve(static_cast<std::size_t>(kernel.children.size));
      for (semantic::Index i = 0; i < kernel.children.size; ++i) {
        children.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                i,
                CompiledMathNodeKind::ExprSurvival,
                condition_id,
                time_id,
                source_view_id));
      }
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::ClampProbability,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(children)),
          CompiledMathValueKind::Survival);
    }
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return complement(
          compile_expr_value_node_raw(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprSurvival,
              condition_id,
              time_id,
              source_view_id));
    }
    std::vector<semantic::Index> terms;
    terms.reserve(static_cast<std::size_t>(kernel.children.size));
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      std::vector<semantic::Index> factors;
      factors.reserve(static_cast<std::size_t>(kernel.children.size));
      factors.push_back(
          compile_expr_child_value_node(
              plan,
              kernel.children,
              i,
              CompiledMathNodeKind::ExprDensity,
              condition_id,
              time_id,
              source_view_id));
      for (semantic::Index j = 0; j < kernel.children.size; ++j) {
        if (i == j) {
          continue;
        }
        factors.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                j,
                CompiledMathNodeKind::ExprSurvival,
                condition_id,
                time_id,
                source_view_id));
      }
      terms.push_back(
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(factors)));
    }
    return compiled_math_algebra_node(
        &plan->compiled_math,
        CompiledMathNodeKind::CleanSignedSum,
        std::move(terms),
        CompiledMathValueKind::Density);
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
    if (kernel.simple_event_guard && !kernel.has_unless) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        return compile_simple_guard_density_node(
            plan,
            expr_id,
            kernel,
            condition_id,
            time_id,
            source_view_id);
      }
      if (value_kind == CompiledMathNodeKind::ExprCdf) {
        return compile_simple_guard_cdf_node(
            plan,
            expr_id,
            kernel,
            condition_id,
            time_id,
            source_view_id);
      }
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
    if (!kernel.has_unless) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        std::vector<semantic::Index> factors;
        factors.reserve(2U);
        factors.push_back(
            compile_expr_value_node(
                plan,
                kernel.guard_ref_expr_id,
                CompiledMathNodeKind::ExprDensity,
                condition_id,
                time_id,
                source_view_id));
        factors.push_back(
            compile_expr_value_node(
                plan,
                kernel.guard_blocker_expr_id,
                CompiledMathNodeKind::ExprSurvival,
                condition_id,
                time_id,
                source_view_id));
        return compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(factors),
            CompiledMathValueKind::Density);
      }
      if (value_kind == CompiledMathNodeKind::ExprCdf) {
        const auto bind_time_id =
            time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
                ? time_id
                : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
        const auto density_node =
            compile_expr_value_node(
                plan,
                expr_id,
                CompiledMathNodeKind::ExprDensity,
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
      }
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
    return unsupported();
  }

  return unsupported();
}

inline semantic::Index compile_expr_upper_bound_node(
    ExactVariantPlan *plan,
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
    const ExactVariantPlan &plan,
    const semantic::Index expr_id) {
  return expr_id != semantic::kInvalidIndex &&
         static_cast<std::size_t>(expr_id) <
             plan.sequence_expr_upper_bound_used.size() &&
         plan.sequence_expr_upper_bound_used[
             static_cast<std::size_t>(expr_id)] != 0U;
}

inline semantic::Index compile_expr_value_node(
    ExactVariantPlan *plan,
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

struct CompiledMathSourceViewCloneMemo {
  struct Entry {
    semantic::Index original_id{semantic::kInvalidIndex};
    semantic::Index source_view_id{0};
    semantic::Index cloned_id{semantic::kInvalidIndex};
  };

  std::vector<Entry> nodes;
  std::vector<Entry> roots;
  std::vector<Entry> conditions;
};

inline bool compose_source_view_template(
    const ExactRelationTemplate &base,
    const ExactRelationTemplate &overlay,
    ExactRelationTemplate *out) {
  *out = base;
  for (std::size_t i = 0; i < overlay.source_ids.size(); ++i) {
    const auto source_id = overlay.source_ids[i];
    const auto relation = overlay.relations[i];
    const auto it =
        std::lower_bound(out->source_ids.begin(), out->source_ids.end(), source_id);
    const auto pos = static_cast<std::size_t>(it - out->source_ids.begin());
    if (it != out->source_ids.end() && *it == source_id) {
      if (out->relations[pos] != relation) {
        return false;
      }
      continue;
    }
    out->source_ids.insert(it, source_id);
    out->relations.insert(out->relations.begin() + static_cast<std::ptrdiff_t>(pos), relation);
  }
  return true;
}

inline bool composed_source_view_id(
    ExactVariantPlan *plan,
    const semantic::Index base_source_view_id,
    const semantic::Index overlay_source_view_id,
    semantic::Index *out_source_view_id) {
  if (overlay_source_view_id == 0 ||
      overlay_source_view_id == semantic::kInvalidIndex) {
    *out_source_view_id =
        base_source_view_id == semantic::kInvalidIndex ? 0 : base_source_view_id;
    return true;
  }
  if (base_source_view_id == 0 ||
      base_source_view_id == semantic::kInvalidIndex) {
    *out_source_view_id = overlay_source_view_id;
    return true;
  }
  const auto base_pos = static_cast<std::size_t>(base_source_view_id - 1U);
  const auto overlay_pos = static_cast<std::size_t>(overlay_source_view_id - 1U);
  if (base_pos >= plan->compiled_source_views.size() ||
      overlay_pos >= plan->compiled_source_views.size()) {
    return false;
  }
  ExactRelationTemplate composed;
  if (!compose_source_view_template(
          plan->compiled_source_views[base_pos],
          plan->compiled_source_views[overlay_pos],
          &composed)) {
    return false;
  }
  *out_source_view_id = compile_source_view_id(plan, composed);
  return true;
}

inline semantic::Index clone_compiled_math_node_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index node_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo);

inline semantic::Index clone_compiled_math_root_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index root_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo) {
  if (root_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto normalized_source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  for (const auto &entry : memo->roots) {
    if (entry.original_id == root_id &&
        entry.source_view_id == normalized_source_view_id) {
      return entry.cloned_id;
    }
  }
  const auto root_pos = static_cast<std::size_t>(root_id);
  if (root_pos >= plan->compiled_math.roots.size()) {
    return semantic::kInvalidIndex;
  }
  const auto cloned_node_id =
      clone_compiled_math_node_for_source_view(
          plan,
          plan->compiled_math.roots[root_pos].node_id,
          normalized_source_view_id,
          memo);
  const auto cloned_root_id =
      compiled_math_make_root(&plan->compiled_math, cloned_node_id);
  memo->roots.push_back(
      CompiledMathSourceViewCloneMemo::Entry{
          root_id, normalized_source_view_id, cloned_root_id});
  return cloned_root_id;
}

inline semantic::Index clone_compiled_math_condition_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index condition_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return 0;
  }
  const auto normalized_source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  for (const auto &entry : memo->conditions) {
    if (entry.original_id == condition_id &&
        entry.source_view_id == normalized_source_view_id) {
      return entry.cloned_id;
    }
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= plan->compiled_math.conditions.size()) {
    return semantic::kInvalidIndex;
  }
  const auto &condition = plan->compiled_math.conditions[condition_pos];
  CompiledMathConditionKey key;
  key.impossible = condition.impossible;
  key.source_ids = condition.source_ids;
  key.relations = condition.relations;
  key.fact_kinds = condition.fact_kinds;
  key.fact_subject_ids = condition.fact_subject_ids;
  key.fact_aux_ids = condition.fact_aux_ids;
  key.fact_aux2_ids = condition.fact_aux2_ids;
  key.fact_time_ids = condition.fact_time_ids;
  key.fact_normalizer_node_ids.reserve(
      condition.fact_normalizer_node_ids.size());
  for (const auto normalizer_node_id : condition.fact_normalizer_node_ids) {
    key.fact_normalizer_node_ids.push_back(
        normalizer_node_id == semantic::kInvalidIndex
            ? semantic::kInvalidIndex
            : clone_compiled_math_node_for_source_view(
                  plan,
                  normalizer_node_id,
                  normalized_source_view_id,
                  memo));
  }
  key.fact_normalizer_root_ids.reserve(
      condition.fact_normalizer_root_ids.size());
  for (const auto normalizer_root_id : condition.fact_normalizer_root_ids) {
    key.fact_normalizer_root_ids.push_back(
        normalizer_root_id == semantic::kInvalidIndex
            ? semantic::kInvalidIndex
            : clone_compiled_math_root_for_source_view(
                  plan,
                  normalizer_root_id,
                  normalized_source_view_id,
                  memo));
  }
  const auto cloned_condition_id =
      compiled_math_intern_condition(&plan->compiled_math, std::move(key));
  memo->conditions.push_back(
      CompiledMathSourceViewCloneMemo::Entry{
          condition_id, normalized_source_view_id, cloned_condition_id});
  return cloned_condition_id;
}

inline semantic::Index clone_compiled_math_node_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index node_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo) {
  if (node_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto normalized_source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  for (const auto &entry : memo->nodes) {
    if (entry.original_id == node_id &&
        entry.source_view_id == normalized_source_view_id) {
      return entry.cloned_id;
    }
  }
  const auto node_pos = static_cast<std::size_t>(node_id);
  if (node_pos >= plan->compiled_math.nodes.size()) {
    return semantic::kInvalidIndex;
  }
  const auto &node = plan->compiled_math.nodes[node_pos];
  semantic::Index composed_view_id{0};
  if (!composed_source_view_id(
          plan,
          normalized_source_view_id,
          node.source_view_id,
          &composed_view_id)) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }

  CompiledMathNodeKey key;
  key.kind = node.kind;
  key.value_kind = node.value_kind;
  key.subject_id = node.subject_id;
  key.condition_id =
      clone_compiled_math_condition_for_source_view(
          plan,
          node.condition_id,
          normalized_source_view_id,
          memo);
  key.time_id = node.time_id;
  key.aux_id = node.aux_id;
  key.aux2_id = node.aux2_id;
  key.source_view_id = composed_view_id;
  key.constant = node.constant;
  key.children.reserve(static_cast<std::size_t>(node.children.size));
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    key.children.push_back(
        clone_compiled_math_node_for_source_view(
            plan,
            plan->compiled_math.child_nodes[
                static_cast<std::size_t>(node.children.offset + i)],
            normalized_source_view_id,
            memo));
  }
  if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent ||
      node.kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw) {
    key.subject_id =
        clone_compiled_math_root_for_source_view(
            plan,
            node.subject_id,
            normalized_source_view_id,
            memo);
  }
  const auto cloned_node_id =
      compiled_math_intern_node(&plan->compiled_math, std::move(key));
  memo->nodes.push_back(
      CompiledMathSourceViewCloneMemo::Entry{
          node_id, normalized_source_view_id, cloned_node_id});
  return cloned_node_id;
}

inline semantic::Index compile_runtime_root_value_node(
    ExactVariantPlan *plan,
    const semantic::Index root_id,
    const semantic::Index source_view_id = 0,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  if (root_id == semantic::kInvalidIndex) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  if ((source_view_id == 0 || source_view_id == semantic::kInvalidIndex) &&
      static_cast<std::size_t>(root_id) < plan->compiled_math.roots.size()) {
    return plan->compiled_math.roots[static_cast<std::size_t>(root_id)].node_id;
  }
  (void)value_kind;
  CompiledMathSourceViewCloneMemo memo;
  const auto cloned_root_id =
      clone_compiled_math_root_for_source_view(
          plan,
          root_id,
          source_view_id,
          &memo);
  if (cloned_root_id == semantic::kInvalidIndex) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  return plan->compiled_math.roots[
      static_cast<std::size_t>(cloned_root_id)].node_id;
}


} // namespace detail
} // namespace accumulatr::eval
