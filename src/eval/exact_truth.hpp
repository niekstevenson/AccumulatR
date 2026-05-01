#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#if defined(__APPLE__)
#include <vecLib/vForce.h>
#define ACCUMULATR_USE_ACCELERATE_VFORCE 1
#endif

#include "batch_source_product.hpp"
#include "exact_source_channels.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
inline void compiled_math_batch_normal_cdf_from_z(
    double *z_by_lane,
    double *exp_scratch,
    const std::size_t count) {
  if (z_by_lane == nullptr || exp_scratch == nullptr || count == 0U) {
    return;
  }
  if (count < 16U) {
    for (std::size_t i = 0; i < count; ++i) {
      z_by_lane[i] = normal_cdf_rational_fast(z_by_lane[i]);
    }
    return;
  }
  for (std::size_t i = 0; i < count; ++i) {
    const double z = z_by_lane[i];
    if (!std::isfinite(z)) {
      exp_scratch[i] = 0.0;
      continue;
    }
    const double x = z * kInvSqrtTwo;
    z_by_lane[i] = x;
    exp_scratch[i] = -x * x;
  }
  const auto n = static_cast<int>(count);
  vvexp(exp_scratch, exp_scratch, &n);
  for (std::size_t i = 0; i < count; ++i) {
    const double x = z_by_lane[i];
    if (!std::isfinite(x)) {
      z_by_lane[i] = x > 0.0 ? 1.0 : 0.0;
      continue;
    }
    const double abs_x = std::abs(x);
    if (abs_x < kInvSqrtTwo) {
      z_by_lane[i] =
          clamp_probability(0.5 + 0.5 * normal_erf_small(x));
      continue;
    }
    double value =
        0.5 * normal_erfc_positive_from_exp(abs_x, exp_scratch[i]);
    if (x > 0.0) {
      value = 1.0 - value;
    }
    z_by_lane[i] = clamp_probability(value);
  }
}
#endif

class CompiledSourceView {
public:
  explicit CompiledSourceView(const ExactVariantPlan &plan)
      : plan_(plan) {}

  void reset(ExactSourceChannels *source_channels,
             const semantic::Index source_view_id = 0) {
    source_channels_ = source_channels;
    source_view_id_ =
        source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  }

  void invalidate_cache() noexcept {}

  ExactSourceChannels *source_channels() const {
    return source_channels_;
  }

  const ExactVariantPlan &plan() const {
    return plan_;
  }

  bool source_order_known_before(
      const semantic::Index before_source_id,
      const semantic::Index after_source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    if (relation_view_knows_before(before_source_id, after_source_id)) {
      return true;
    }
    return source_channels_->source_order_known_before(
        before_source_id, after_source_id, condition_id, workspace);
  }

  const double *exact_time_for_source(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    return source_channels_->exact_time_for_source(
        source_id, condition_id, workspace);
  }

  bool expr_upper_bound_for(
      const semantic::Index expr_id,
      ExactTimedExprUpperBound *out) const {
    return source_channels_ != nullptr &&
           source_channels_->expr_upper_bound_for(expr_id, out);
  }

private:
  bool relation_view_knows_before(const semantic::Index before_source_id,
                                  const semantic::Index after_source_id) const {
    if (before_source_id == semantic::kInvalidIndex ||
        after_source_id == semantic::kInvalidIndex ||
        before_source_id == after_source_id) {
      return false;
    }
    const auto before_relation =
        exact_compiled_source_view_relation(
            plan_, source_view_id_, before_source_id);
    const auto after_relation =
        exact_compiled_source_view_relation(
            plan_, source_view_id_, after_source_id);
    return before_relation != ExactRelation::Unknown &&
           after_relation != ExactRelation::Unknown &&
           before_relation < after_relation;
  }

public:
  semantic::Index condition_cache_id() const noexcept {
    return source_view_id_;
  }

  semantic::Index source_view_id() const noexcept {
    return source_view_id_;
  }

private:
  const ExactVariantPlan &plan_;
  ExactSourceChannels *source_channels_{nullptr};
  semantic::Index source_view_id_{0};
};

inline semantic::Index exact_condition_cache_id(
    const CompiledSourceView *evaluator) {
  return evaluator == nullptr ? 0 : evaluator->condition_cache_id();
}

inline void exact_hash_condition_part(std::size_t *seed,
                                      const std::size_t value) noexcept {
  *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
}

inline semantic::Index compiled_node_condition_cache_id(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace) {
  const auto evaluator_id = exact_condition_cache_id(evaluator);
  if (node.condition_id == 0 ||
      node.condition_id == semantic::kInvalidIndex) {
    return evaluator_id;
  }
  const auto condition_pos = static_cast<std::size_t>(node.condition_id);
  const auto &cache_plan = program.condition_cache_plans[condition_pos];
  if (!cache_plan.dynamic) {
    return compiled_math_static_source_view_condition_cache_id(
        evaluator_id,
        node.condition_id);
  }
  std::size_t seed = 0x9e3779b97f4a7c15ULL;
  exact_hash_condition_part(&seed, static_cast<std::size_t>(evaluator_id));
  exact_hash_condition_part(&seed, static_cast<std::size_t>(node.condition_id));
  for (semantic::Index i = 0; i < cache_plan.time_dependencies.size; ++i) {
    const auto time_id =
        program.condition_cache_time_dependencies[
            static_cast<std::size_t>(
                cache_plan.time_dependencies.offset + i)];
    const auto time_pos = static_cast<std::size_t>(time_id);
    exact_hash_condition_part(&seed, static_cast<std::size_t>(time_id) + 1U);
    exact_hash_condition_part(
        &seed, std::hash<double>{}(workspace.time_values[time_pos]));
  }
  return static_cast<semantic::Index>(
      (seed & static_cast<std::size_t>(0x3fffffffU)) + 1U);
}

struct CompiledEvalWorkspace {
  explicit CompiledEvalWorkspace(const ExactVariantPlan &plan)
      : evaluator(plan),
        compiled_math(plan.compiled_math) {
    const auto source_view_count = plan.compiled_source_views.size();
    source_view_evaluators.reserve(source_view_count);
    for (std::size_t i = 0; i < source_view_count; ++i) {
      source_view_evaluators.emplace_back(plan);
    }
  }

  void reset(ExactSourceChannels *source_channels,
             const semantic::Index source_view_id = 0) {
    evaluator.reset(source_channels, source_view_id);
    if (source_views_channels_ != source_channels) {
      for (std::size_t i = 0; i < source_view_evaluators.size(); ++i) {
        source_view_evaluators[i].reset(
            source_channels, static_cast<semantic::Index>(i + 1U));
      }
      source_views_channels_ = source_channels;
    }
    compiled_math.reset_cache();
  }

  void reset_planned_caches() {}

  CompiledSourceView *source_view_evaluator(
      const semantic::Index source_view_id,
      CompiledSourceView *parent) {
    if (source_view_id == 0 || source_view_id == semantic::kInvalidIndex) {
      return parent;
    }
    if (parent == nullptr) {
      return nullptr;
    }
    const auto pos = static_cast<std::size_t>(source_view_id - 1U);
    if (pos >= parent->plan().compiled_source_views.size() ||
        pos >= source_view_evaluators.size()) {
      throw std::runtime_error("compiled source view id is outside the plan");
    }
    return &source_view_evaluators[pos];
  }

  CompiledSourceView evaluator;
  CompiledMathWorkspace compiled_math;
  std::vector<CompiledSourceView> source_view_evaluators;
  ExactSourceChannels *source_views_channels_{nullptr};
};

inline bool compiled_math_node_cacheable(
    const CompiledMathNodeKind kind) noexcept {
  return kind == CompiledMathNodeKind::SimpleGuardCdf ||
         kind == CompiledMathNodeKind::IntegralZeroToCurrent ||
         kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw ||
         kind == CompiledMathNodeKind::UnionKernelDensity ||
         kind == CompiledMathNodeKind::UnionKernelCdf ||
         kind == CompiledMathNodeKind::UnionKernelMultiSubsetDensity ||
         kind == CompiledMathNodeKind::UnionKernelMultiSubsetCdf;
}

inline bool compiled_math_load_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value);

inline void compiled_math_store_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    double value);

inline double compiled_math_node_time(const CompiledMathNode &node,
                                      const CompiledMathWorkspace &workspace) {
  return workspace.time_values[static_cast<std::size_t>(node.time_id)];
}

inline double compiled_math_source_node_time(
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace) {
  const double time = compiled_math_node_time(node, workspace);
  if (node.aux_id == semantic::kInvalidIndex) {
    return time;
  }
  return std::min(time, workspace.time(node.aux_id));
}

struct CompiledTimedUpperBound {
  bool found{false};
  double time{std::numeric_limits<double>::infinity()};
  double normalizer{0.0};
};

inline double compiled_math_timed_upper_bound_term_time(
    const CompiledMathTimedUpperBoundTerm &term,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace) {
  const auto time_id = term.time_id;
  if (time_id != semantic::kInvalidIndex && workspace.has_time(time_id)) {
    return workspace.time(time_id);
  }
  return compiled_math_node_time(node, workspace);
}

inline double compiled_math_timed_upper_bound_term_normalizer(
    const CompiledMathTimedUpperBoundTerm &term,
    const CompiledMathWorkspace &workspace) {
  return workspace.values[static_cast<std::size_t>(term.normalizer_node_id)];
}

inline CompiledTimedUpperBound compiled_math_upper_bound_from_terms(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace) {
  CompiledTimedUpperBound best;
  if (node.aux_id == semantic::kInvalidIndex ||
      node.aux2_id == semantic::kInvalidIndex) {
    return best;
  }
  const auto offset = static_cast<std::size_t>(node.aux_id);
  const auto size = static_cast<std::size_t>(node.aux2_id);
  for (std::size_t i = 0; i < size; ++i) {
    const auto &term = program.timed_upper_bound_terms[offset + i];
    const double fact_time =
        compiled_math_timed_upper_bound_term_time(term, node, workspace);
    const double normalizer =
        compiled_math_timed_upper_bound_term_normalizer(term, workspace);
    if (std::isfinite(fact_time) && normalizer > 0.0 &&
        (!best.found || fact_time < best.time)) {
      best.found = true;
      best.time = fact_time;
      best.normalizer = normalizer;
    }
  }
  return best;
}

inline CompiledTimedUpperBound compiled_math_expr_upper_bound_for_node(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const CompiledSourceView *condition_evaluator) {
  CompiledTimedUpperBound best =
      compiled_math_upper_bound_from_terms(program, node, workspace);
  if (condition_evaluator != nullptr) {
    ExactTimedExprUpperBound sequence_upper;
    if (condition_evaluator->expr_upper_bound_for(
            node.subject_id,
            &sequence_upper) &&
        std::isfinite(sequence_upper.time) &&
        sequence_upper.normalizer > 0.0 &&
        (!best.found || sequence_upper.time < best.time)) {
      best.found = true;
      best.time = sequence_upper.time;
      best.normalizer = sequence_upper.normalizer;
    }
  }
  return best;
}

inline bool compiled_math_load_node_cache(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value) {
  if (evaluator == nullptr || !compiled_math_node_cacheable(node.kind)) {
    return false;
  }
  return compiled_math_load_node_cache_entry(
      program, node, evaluator, workspace, value);
}

inline bool compiled_math_load_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value) {
  if (evaluator == nullptr) {
    return false;
  }
  const auto slot = static_cast<std::size_t>(node.cache_slot);
  if (slot >= workspace.cache_valid.size() ||
      workspace.cache_valid[slot] == 0U ||
      workspace.cache_evaluators[slot] != evaluator ||
      workspace.cache_condition_ids[slot] !=
          compiled_node_condition_cache_id(program, node, evaluator, workspace) ||
      workspace.cache_times[slot] !=
          compiled_math_node_time(node, workspace)) {
    return false;
  }
  *value = workspace.cache_values[slot];
  return true;
}

inline CompiledSourceView *compiled_math_node_evaluator(
    const CompiledMathNode &node,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *workspace) {
  if (node.source_view_id == 0 ||
      node.source_view_id == semantic::kInvalidIndex) {
    return parent;
  }
  if (workspace == nullptr) {
    throw std::runtime_error(
        "compiled node requires a planned source view but no workspace was supplied");
  }
  return workspace->source_view_evaluator(node.source_view_id, parent);
}

inline void compiled_math_store_node_cache(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    const double value) {
  if (evaluator == nullptr || !compiled_math_node_cacheable(node.kind)) {
    return;
  }
  compiled_math_store_node_cache_entry(
      program, node, evaluator, workspace, value);
}

inline void compiled_math_store_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    const double value) {
  if (evaluator == nullptr) {
    return;
  }
  const auto slot = static_cast<std::size_t>(node.cache_slot);
  if (slot >= workspace->cache_valid.size()) {
    return;
  }
  workspace->cache_valid[slot] = 1U;
  workspace->cache_evaluators[slot] = evaluator;
  workspace->cache_condition_ids[slot] =
      compiled_node_condition_cache_id(program, node, evaluator, *workspace);
  workspace->cache_times[slot] = compiled_math_node_time(node, *workspace);
  workspace->cache_values[slot] = value;
}

inline double compiled_union_kernel_child_value(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const semantic::Index child_index) {
  if (child_index >= node.children.size) {
    return 0.0;
  }
  const auto child_id = program.child_nodes[
      static_cast<std::size_t>(node.children.offset + child_index)];
  return workspace.values[static_cast<std::size_t>(child_id)];
}

inline double compiled_union_subset_conjunction_density(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const ExactVariantPlan &plan,
    const ExactExprUnionSubset &subset,
    const semantic::Index child_count) {
  if (subset.child_positions.empty()) {
    return 0.0;
  }
  if (subset.child_positions.size == 1U) {
    const auto active_pos = plan.expr_union_subset_child_positions[
        static_cast<std::size_t>(subset.child_positions.offset)];
    return compiled_union_kernel_child_value(
        program, node, workspace, active_pos);
  }
  double total = 0.0;
  for (semantic::Index i = 0; i < subset.child_positions.size; ++i) {
    const auto active_pos = plan.expr_union_subset_child_positions[
        static_cast<std::size_t>(subset.child_positions.offset + i)];
    double term =
        compiled_union_kernel_child_value(program, node, workspace, active_pos);
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (semantic::Index j = 0; j < subset.child_positions.size; ++j) {
      if (i == j) {
        continue;
      }
      const auto other_pos = plan.expr_union_subset_child_positions[
          static_cast<std::size_t>(subset.child_positions.offset + j)];
      term *= compiled_union_kernel_child_value(
          program,
          node,
          workspace,
          child_count + other_pos);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    total += term;
  }
  return clean_signed_value(total);
}

inline double evaluate_compiled_union_kernel_density(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const ExactVariantPlan &plan,
    CompiledSourceView *evaluator,
    const bool multi_subset_only) {
  (void)evaluator;
  const auto expr_pos = static_cast<std::size_t>(node.subject_id);
  if (expr_pos >= plan.expr_kernels.size()) {
    return 0.0;
  }
  const auto &kernel = plan.expr_kernels[expr_pos];
  if (kernel.children.empty()) {
    return 0.0;
  }
  const auto child_count = kernel.children.size;
  if (!multi_subset_only && child_count == 1U) {
    return compiled_union_kernel_child_value(program, node, workspace, 0);
  }
  double total = 0.0;
  for (semantic::Index i = 0; i < kernel.union_subset_span.size; ++i) {
    const auto &subset = plan.expr_union_subsets[
        static_cast<std::size_t>(kernel.union_subset_span.offset + i)];
    if (multi_subset_only && subset.child_positions.size <= 1U) {
      continue;
    }
    const double value =
        compiled_union_subset_conjunction_density(
            program,
            node,
            workspace,
            plan,
            subset,
            child_count);
    total += static_cast<double>(subset.sign) * value;
  }
  return clean_signed_value(total);
}

inline double evaluate_compiled_union_kernel_cdf(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    CompiledSourceView *evaluator) {
  (void)evaluator;
  double total = 0.0;
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    total += compiled_union_kernel_child_value(program, node, workspace, i);
  }
  return clamp_probability(total);
}

inline double evaluate_compiled_math_root(
    const CompiledMathProgram &program,
    semantic::Index root_id,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace = nullptr);

inline double evaluate_compiled_node_span(
    const CompiledMathProgram &program,
    CompiledMathIndexSpan schedule,
    semantic::Index result_node_id,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace = nullptr,
    const std::vector<semantic::Index> *schedule_nodes = nullptr,
    bool workspace_prepared = false);

inline double evaluate_compiled_integral_kernel(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace,
    double lower,
    double upper);

inline bool compiled_math_outcome_gate_open_for_node(
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const CompiledSourceView *parent) {
  if (node.kind != CompiledMathNodeKind::OutcomeSubsetUnused) {
    throw std::runtime_error("integral outcome gate contains a non-gate node");
  }
  if (workspace.used_outcomes == nullptr) {
    return true;
  }
  if (parent == nullptr) {
    throw std::runtime_error(
        "integral outcome gate requires a compiled source view");
  }
  const auto offset = static_cast<std::size_t>(node.subject_id);
  const auto size = static_cast<std::size_t>(node.aux_id);
  const auto &indices = parent->plan().compiled_outcome_gate_indices;
  if (offset + size > indices.size()) {
    throw std::runtime_error(
        "compiled outcome-used gate points outside the plan");
  }
  for (std::size_t i = 0; i < size; ++i) {
    const auto outcome_idx = indices[offset + i];
    if (outcome_idx != semantic::kInvalidIndex &&
        static_cast<std::size_t>(outcome_idx) <
            workspace.used_outcomes->size() &&
        (*workspace.used_outcomes)[static_cast<std::size_t>(outcome_idx)] !=
            0U) {
      return false;
    }
  }
  return true;
}

inline bool compiled_math_integral_outcome_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan outcome_gate_nodes,
    const CompiledMathWorkspace &workspace,
    const CompiledSourceView *parent) {
  for (semantic::Index i = 0; i < outcome_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_outcome_gate_nodes[
            static_cast<std::size_t>(outcome_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (!compiled_math_outcome_gate_open_for_node(
            gate_node,
            workspace,
            parent)) {
      return false;
    }
  }
  return true;
}

inline bool compiled_math_integral_time_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan time_gate_nodes,
    const CompiledMathWorkspace &workspace) {
  for (semantic::Index i = 0; i < time_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_time_gate_nodes[
            static_cast<std::size_t>(time_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (gate_node.kind != CompiledMathNodeKind::TimeGate) {
      throw std::runtime_error("integral time gate contains a non-gate node");
    }
    if (workspace.time_values[static_cast<std::size_t>(gate_node.time_id)] <
        workspace.time_values[static_cast<std::size_t>(gate_node.aux_id)]) {
      return false;
    }
  }
  return true;
}

inline double compiled_math_integral_factor_value_for_node(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace) {
  if (!compiled_math_is_integral_node(node.kind)) {
    throw std::runtime_error("integral factor contains a non-integral node");
  }
  const double current_time = compiled_math_node_time(node, *workspace);
  if (!(current_time > 0.0)) {
    return 0.0;
  }
  double value =
      evaluate_compiled_integral_kernel(
          program,
          node,
          workspace,
          parent,
          scenario,
          eval_workspace,
          0.0,
          current_time);
  if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent) {
    value = clamp_probability(value);
  }
  return value;
}

inline bool compiled_math_apply_expr_upper_factor(
    const CompiledMathProgram &program,
    const CompiledMathIntegralExprUpperFactor &factor,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *eval_workspace,
    double *product) {
  const auto &node =
      program.nodes[static_cast<std::size_t>(factor.node_id)];
  if (node.kind != CompiledMathNodeKind::ExprUpperBoundDensity &&
      node.kind != CompiledMathNodeKind::ExprUpperBoundCdf) {
    throw std::runtime_error(
        "integral expression upper factor contains a non-upper node");
  }
  const auto upper =
      compiled_math_expr_upper_bound_for_node(
          program,
          node,
          *workspace,
          compiled_math_node_evaluator(node, parent, eval_workspace));
  const bool has_upper = upper.found && upper.normalizer > 0.0;
  const double current_time = compiled_math_node_time(node, *workspace);
  if (factor.mode == CompiledMathIntegralExprUpperMode::AfterOne) {
    if (!has_upper || current_time < upper.time) {
      *product = 0.0;
      return false;
    }
    return true;
  }
  if (!has_upper) {
    return true;
  }
  if (current_time >= upper.time) {
    *product = 0.0;
    return false;
  }
  *product /= upper.normalizer;
  return *product != 0.0;
}

inline double compiled_math_source_product_channel_time(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathWorkspace &workspace) {
  const double time =
      workspace.time_values[static_cast<std::size_t>(channel.time_id)];
  if (channel.time_cap_id == semantic::kInvalidIndex) {
    return time;
  }
  return std::min(
      time,
      workspace.time_values[static_cast<std::size_t>(channel.time_cap_id)]);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_finish_base_fill(
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t fill_mask) {
  const double start_prob = 1.0 - q;
  ExactSourceChannels::SourceProductScalarFill fill;
  fill.mask = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    fill.pdf = start_prob * safe_density(base_pdf);
  }
  if ((fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    const double cdf = clamp_probability(start_prob * clamp_probability(base_cdf));
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      fill.cdf = cdf;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      fill.survival = clamp_probability(1.0 - cdf);
    }
  }
  return fill;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_lognormal_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double m = loaded.params[0];
  const double s = loaded.params[1];
  return compiled_math_source_product_finish_base_fill(
      need_pdf ? lognormal_pdf_fast(x, m, s) : 0.0,
      need_cdf ? lognormal_cdf_fast(x, m, s) : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_gamma_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double shape = loaded.params[0];
  const double rate = loaded.params[1];
  const double scale = 1.0 / rate;
  return compiled_math_source_product_finish_base_fill(
      need_pdf ? gamma_pdf_fast_rate(x, shape, rate) : 0.0,
      need_cdf ? R::pgamma(x, shape, scale, 1, 0) : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_exgauss_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double mu = loaded.params[0];
  const double sigma = loaded.params[1];
  const double tau = loaded.params[2];
  const double lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
  const double lower_survival = 1.0 - lower_cdf;
  if (!(lower_survival > 0.0)) {
    ExactSourceChannels::SourceProductScalarFill fill;
    fill.mask = fill_mask;
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0;
    }
    return fill;
  }
  return compiled_math_source_product_finish_base_fill(
      need_pdf ? exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival : 0.0,
      need_cdf
          ? (exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) /
                lower_survival
          : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_lba_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return compiled_math_source_product_finish_base_fill(
      need_pdf
          ? lba_pdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      need_cdf
          ? lba_cdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_rdm_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return compiled_math_source_product_finish_base_fill(
      need_pdf
          ? rdm_pdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      need_cdf
          ? rdm_cdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_direct_leaf_fill(
    const std::uint8_t leaf_dist_kind,
    const ExactLoadedLeafInput &loaded,
    const double time_since_onset,
    const std::uint8_t fill_mask) {
  const double x = time_since_onset - loaded.t0;
  if (!(x > 0.0)) {
    ExactSourceChannels::SourceProductScalarFill fill;
    fill.mask = fill_mask;
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0;
    }
    return fill;
  }
  switch (static_cast<leaf::DistKind>(leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    return compiled_math_source_product_lognormal_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::Gamma:
    return compiled_math_source_product_gamma_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::Exgauss:
    return compiled_math_source_product_exgauss_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::LBA:
    return compiled_math_source_product_lba_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::RDM:
    return compiled_math_source_product_rdm_leaf_fill(
        loaded, x, fill_mask);
  }
  return ExactSourceChannels::SourceProductScalarFill{};
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_direct_leaf_fill(
    const CompiledMathSourceProductChannel &channel,
    const ExactLoadedLeafInput &loaded,
    const double time,
    const std::uint8_t fill_mask) {
  return compiled_math_source_product_direct_leaf_fill(
      channel.leaf_dist_kind,
      loaded,
      time - channel.leaf_onset_abs_value,
      fill_mask);
}

inline double compiled_math_source_product_finish_base_scalar(
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t channel_mask) {
  const double start_prob = 1.0 - q;
  if (channel_mask == kLeafChannelPdf) {
    return start_prob * safe_density(base_pdf);
  }
  const double cdf = clamp_probability(start_prob * clamp_probability(base_cdf));
  if (channel_mask == kLeafChannelCdf) {
    return cdf;
  }
  return channel_mask == kLeafChannelSurvival
             ? clamp_probability(1.0 - cdf)
             : 0.0;
}

inline double compiled_math_source_product_direct_leaf_scalar(
    const CompiledMathSourceProductChannel &channel,
    const ExactLoadedLeafInput &loaded,
    const double time,
    const std::uint8_t channel_mask) {
  const double x = time - channel.leaf_onset_abs_value - loaded.t0;
  if (!(x > 0.0)) {
    return channel_mask == kLeafChannelSurvival ? 1.0 : 0.0;
  }
  const bool need_pdf = channel_mask == kLeafChannelPdf;
  switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
  case leaf::DistKind::Lognormal: {
    const double m = loaded.params[0];
    const double s = loaded.params[1];
    return compiled_math_source_product_finish_base_scalar(
        need_pdf ? lognormal_pdf_fast(x, m, s) : 0.0,
        need_pdf ? 0.0 : lognormal_cdf_fast(x, m, s),
        loaded.q,
        channel_mask);
  }
  case leaf::DistKind::Gamma: {
    const double shape = loaded.params[0];
    const double rate = loaded.params[1];
    const double scale = 1.0 / rate;
    return compiled_math_source_product_finish_base_scalar(
        need_pdf ? gamma_pdf_fast_rate(x, shape, rate) : 0.0,
        need_pdf ? 0.0 : R::pgamma(x, shape, scale, 1, 0),
        loaded.q,
        channel_mask);
  }
  case leaf::DistKind::Exgauss: {
    const double mu = loaded.params[0];
    const double sigma = loaded.params[1];
    const double tau = loaded.params[2];
    const double lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
    const double lower_survival = 1.0 - lower_cdf;
    if (!(lower_survival > 0.0)) {
      return channel_mask == kLeafChannelSurvival ? 1.0 : 0.0;
    }
    return compiled_math_source_product_finish_base_scalar(
        need_pdf ? exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival : 0.0,
        need_pdf
            ? 0.0
            : (exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) /
                  lower_survival,
        loaded.q,
        channel_mask);
  }
  case leaf::DistKind::LBA:
    return compiled_math_source_product_finish_base_scalar(
        need_pdf
            ? lba_pdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3])
            : 0.0,
        need_pdf
            ? 0.0
            : lba_cdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3]),
        loaded.q,
        channel_mask);
  case leaf::DistKind::RDM:
    return compiled_math_source_product_finish_base_scalar(
        need_pdf
            ? rdm_pdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3])
            : 0.0,
        need_pdf
            ? 0.0
            : rdm_cdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3]),
        loaded.q,
        channel_mask);
  }
  return 0.0;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_impossible_fill(const std::uint8_t mask) {
  ExactSourceChannels::SourceProductScalarFill fill;
  fill.mask = mask;
  return fill;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_certain_fill(const std::uint8_t mask) {
  ExactSourceChannels::SourceProductScalarFill fill;
  fill.mask = mask;
  if ((mask & kLeafChannelCdf) != 0U) {
    fill.cdf = 1.0;
  }
  if ((mask & kLeafChannelSurvival) != 0U) {
    fill.survival = 0.0;
  }
  return fill;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_forced_fill(
    const ExactRelation relation,
    const std::uint8_t mask) {
  if (relation == ExactRelation::Before || relation == ExactRelation::At) {
    return compiled_math_source_product_certain_fill(mask);
  }
  if (relation == ExactRelation::After) {
    ExactSourceChannels::SourceProductScalarFill fill;
    fill.mask = mask;
    if ((mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0;
    }
    return fill;
  }
  return compiled_math_source_product_impossible_fill(mask);
}

inline bool compiled_math_source_product_bounds_have_overlay(
    const CompiledSourceBoundPlan &bounds) noexcept {
  return bounds.has_condition_exact ||
         bounds.has_condition_lower ||
         bounds.has_condition_upper;
}

inline bool compiled_math_source_product_relation_forces_fill(
    const ExactRelation relation,
    const std::uint8_t fill_mask) noexcept {
  return relation != ExactRelation::Unknown &&
         !(relation == ExactRelation::At &&
           (fill_mask & kLeafChannelPdf) != 0U);
}

inline ExactRelation compiled_math_source_product_program_relation(
    const CompiledMathSourceProductProgram &source_program) noexcept {
  return source_program.has_static_source_view_relation
             ? static_cast<ExactRelation>(
                   source_program.static_source_view_relation)
             : ExactRelation::Unknown;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_conditionalize(
    const ExactSourceChannels::SourceProductScalarFill uncond,
    const ExactSourceChannels::SourceProductScalarFill lower,
    const std::uint8_t fill_mask) {
  const double surv_lb = lower.survival;
  if (!std::isfinite(surv_lb) || !(surv_lb > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    out.pdf = safe_density(uncond.pdf / surv_lb);
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    out.cdf = clamp_probability((uncond.cdf + surv_lb - 1.0) / surv_lb);
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    out.survival = clamp_probability(uncond.survival / surv_lb);
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_conditionalize_between(
    const ExactSourceChannels::SourceProductScalarFill uncond,
    const ExactSourceChannels::SourceProductScalarFill lower,
    const ExactSourceChannels::SourceProductScalarFill upper,
    const std::uint8_t fill_mask) {
  const double mass = upper.cdf - lower.cdf;
  if (!std::isfinite(mass) || !(mass > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    out.pdf = safe_density(uncond.pdf / mass);
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    out.cdf = clamp_probability((uncond.cdf - lower.cdf) / mass);
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    out.survival = clamp_probability((upper.cdf - uncond.cdf) / mass);
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask);

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_leaf_fill(
    const CompiledMathSourceProductProgram &source_program,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  return compiled_math_source_product_direct_leaf_fill(
      source_program.leaf_dist_kind,
      source_channels->source_product_leaf_input(source_program.leaf_index),
      current_time - source_program.leaf_onset_abs_value,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_exact_gate_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !compiled_math_source_product_bounds_have_overlay(
          source_program.bounds)) {
    return compiled_math_source_product_forced_fill(relation, fill_mask);
  }
  const auto bounds =
      source_channels->source_product_resolved_bounds(
          source_program.source_id,
          source_program.bounds,
          workspace);
  if (const double *exact_time = bounds.exact_time) {
    if (!(current_time >= *exact_time)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    return compiled_math_source_product_certain_fill(fill_mask);
  }
  return compiled_math_source_product_program_fill(
      program,
      source_program.child_program_id,
      workspace,
      source_channels,
      current_time,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_conditioned_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  const bool has_condition_overlay =
      compiled_math_source_product_bounds_have_overlay(
          source_program.bounds);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !has_condition_overlay) {
    return compiled_math_source_product_forced_fill(relation, fill_mask);
  }

  const auto bounds =
      source_channels->source_product_resolved_bounds(
          source_program.source_id,
          source_program.bounds,
          workspace);
  const double lower_bound = bounds.lower;
  const double upper_bound = bounds.upper;
  if (std::isfinite(upper_bound) && !(upper_bound > lower_bound)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  if (const double *exact_time = bounds.exact_time) {
    if (lower_bound > 0.0 && !(*exact_time > lower_bound)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    if (std::isfinite(upper_bound) && !(*exact_time < upper_bound)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    if (!(current_time >= *exact_time)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    return compiled_math_source_product_certain_fill(fill_mask);
  }

  std::uint8_t uncond_mask = fill_mask;
  if (std::isfinite(upper_bound) &&
      (fill_mask & kLeafChannelSurvival) != 0U) {
    uncond_mask |= kLeafChannelCdf;
  }
  const auto uncond =
      compiled_math_source_product_program_fill(
          program,
          source_program.child_program_id,
          workspace,
          source_channels,
          current_time,
          uncond_mask);
  if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
    return uncond;
  }

  std::uint8_t lower_mask = kLeafChannelSurvival;
  if (std::isfinite(upper_bound)) {
    lower_mask = kLeafChannelCdf;
  }
  auto lower = compiled_math_source_product_impossible_fill(lower_mask);
  if (lower_bound > 0.0) {
    lower =
        compiled_math_source_product_program_fill(
            program,
            source_program.child_program_id,
            workspace,
            source_channels,
            lower_bound,
            lower_mask);
  }
  if (!std::isfinite(upper_bound)) {
    return compiled_math_source_product_conditionalize(
        uncond, lower, fill_mask);
  }
  const auto upper =
      compiled_math_source_product_program_fill(
          program,
          source_program.child_program_id,
          workspace,
          source_channels,
          upper_bound,
          kLeafChannelCdf);
  if (!std::isfinite(upper.cdf - lower.cdf) ||
      !(upper.cdf - lower.cdf > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  if (current_time >= upper_bound) {
    return compiled_math_source_product_certain_fill(fill_mask);
  }
  if (current_time <= lower_bound) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  return compiled_math_source_product_conditionalize_between(
      uncond, lower, upper, fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_onset_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const double upper = current_time - source_program.leaf_onset_lag;
  if (!(upper > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }

  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const std::uint8_t shifted_mask =
      (need_pdf ? kLeafChannelPdf : 0U) |
      (need_cdf ? kLeafChannelCdf : 0U);
  const auto &loaded =
      source_channels->source_product_leaf_input(source_program.leaf_index);
  auto shifted_fill = [&](const double source_time) {
    return compiled_math_source_product_direct_leaf_fill(
        source_program.leaf_dist_kind,
        loaded,
        current_time - source_time - source_program.leaf_onset_lag,
        shifted_mask);
  };

  const auto &onset_program =
      program.integral_kernel_source_product_programs[
          static_cast<std::size_t>(
              source_program.onset_source_program_id)];
  const auto onset_bounds =
      source_channels->source_product_resolved_bounds(
          onset_program.source_id,
          source_program.onset_bounds,
          workspace);
  if (const double *exact_time = onset_bounds.exact_time) {
    const auto shifted = shifted_fill(*exact_time);
    ExactSourceChannels::SourceProductScalarFill out;
    out.mask = fill_mask;
    if (need_pdf) {
      out.pdf = shifted.pdf;
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = shifted.cdf;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability(1.0 - shifted.cdf);
    }
    return out;
  }

  const auto batch = quadrature::build_finite_batch(0.0, upper);
  double pdf = 0.0;
  double cdf = 0.0;
  for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
    const double u = batch.nodes.nodes[i];
    const auto source =
        compiled_math_source_product_program_fill(
            program,
            source_program.onset_source_program_id,
            workspace,
            source_channels,
            u,
            kLeafChannelPdf);
    if (!(source.pdf > 0.0)) {
      continue;
    }
    const auto shifted = shifted_fill(u);
    const double weight = batch.nodes.weights[i] * source.pdf;
    if (need_pdf) {
      pdf += weight * shifted.pdf;
    }
    if (need_cdf) {
      cdf += weight * shifted.cdf;
    }
  }

  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if (need_pdf) {
    out.pdf = safe_density(pdf);
  }
  if (need_cdf) {
    const double clamped = clamp_probability(cdf);
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamped;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability(1.0 - clamped);
    }
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_pool_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const int n_members = static_cast<int>(source_program.member_programs.size);
  const int k = static_cast<int>(source_program.pool_k);
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const std::uint8_t member_mask =
      kLeafChannelCdf | kLeafChannelSurvival |
      (need_pdf ? kLeafChannelPdf : 0U);

  double *scratch =
      workspace->source_product_scratch.data() +
      static_cast<std::size_t>(
          source_program.source_product_scratch_offset);
  double *member_pdf = scratch;
  double *member_cdf = member_pdf + n_members;
  double *member_survival = member_cdf + n_members;
  const int width = n_members + 1;
  const std::size_t table_size =
      static_cast<std::size_t>(width) * static_cast<std::size_t>(width);
  double *prefix = member_survival + n_members;
  double *suffix = prefix + table_size;
  std::fill(prefix, prefix + table_size, 0.0);
  std::fill(suffix, suffix + table_size, 0.0);

  for (semantic::Index i = 0; i < source_program.member_programs.size; ++i) {
    const auto member_program_id =
        program.integral_kernel_source_product_program_members[
            static_cast<std::size_t>(
                source_program.member_programs.offset + i)];
    const auto member =
        compiled_math_source_product_program_fill(
            program,
            member_program_id,
            workspace,
            source_channels,
            current_time,
            member_mask);
    const auto pos = static_cast<std::size_t>(i);
    member_pdf[pos] = member.pdf;
    member_cdf[pos] = member.cdf;
    member_survival[pos] = member.survival;
  }

  const auto idx = [width](const int row, const int col) {
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
           static_cast<std::size_t>(col);
  };
  prefix[idx(0, 0)] = 1.0;
  for (int i = 0; i < n_members; ++i) {
    for (int m = 0; m <= i; ++m) {
      prefix[idx(i + 1, m)] +=
          prefix[idx(i, m)] * member_survival[static_cast<std::size_t>(i)];
      prefix[idx(i + 1, m + 1)] +=
          prefix[idx(i, m)] * member_cdf[static_cast<std::size_t>(i)];
    }
  }
  suffix[idx(n_members, 0)] = 1.0;
  for (int i = n_members - 1; i >= 0; --i) {
    const int count = n_members - i - 1;
    for (int m = 0; m <= count; ++m) {
      suffix[idx(i, m)] +=
          suffix[idx(i + 1, m)] *
          member_survival[static_cast<std::size_t>(i)];
      suffix[idx(i, m + 1)] +=
          suffix[idx(i + 1, m)] *
          member_cdf[static_cast<std::size_t>(i)];
    }
  }

  double pdf = 0.0;
  if (need_pdf) {
    for (int i = 0; i < n_members; ++i) {
      double others_exact = 0.0;
      for (int left = 0; left < k; ++left) {
        const int right = (k - 1) - left;
        if (right < 0 || right > (n_members - i - 1)) {
          continue;
        }
        others_exact +=
            prefix[idx(i, left)] * suffix[idx(i + 1, right)];
      }
      pdf += member_pdf[static_cast<std::size_t>(i)] * others_exact;
    }
  }
  double surv = 0.0;
  if (need_cdf) {
    for (int m = 0; m < k; ++m) {
      surv += prefix[idx(n_members, m)];
    }
  }

  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if (need_pdf) {
    out.pdf = safe_density(pdf);
  }
  if (need_cdf) {
    const double clamped_survival = clamp_probability(surv);
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamped_survival;
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamp_probability(1.0 - clamped_survival);
    }
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  if (source_product_program_id == semantic::kInvalidIndex) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  const auto &source_program =
      program.integral_kernel_source_product_programs[
          static_cast<std::size_t>(source_product_program_id)];
  switch (source_program.kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
    return compiled_math_source_product_impossible_fill(fill_mask);
  case CompiledMathSourceProductProgramKind::ConstantOne:
    return compiled_math_source_product_certain_fill(fill_mask);
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    return compiled_math_source_product_program_leaf_fill(
        source_program, source_channels, current_time, fill_mask);
  case CompiledMathSourceProductProgramKind::ExactGate:
    return compiled_math_source_product_program_exact_gate_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  case CompiledMathSourceProductProgramKind::Conditioned:
    return compiled_math_source_product_program_conditioned_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
    return compiled_math_source_product_program_onset_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    return compiled_math_source_product_program_pool_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  }
  return compiled_math_source_product_impossible_fill(fill_mask);
}

inline void compiled_math_store_source_product_program_fill(
    CompiledMathWorkspace *workspace,
    const std::size_t program_pos,
    const ExactSourceChannels::SourceProductScalarFill &fill) {
  if (workspace->source_product_program_epoch[program_pos] !=
      workspace->source_product_program_current_epoch) {
    workspace->source_product_program_epoch[program_pos] =
        workspace->source_product_program_current_epoch;
    workspace->source_product_program_valid_mask[program_pos] = 0U;
  }
  workspace->source_product_program_valid_mask[program_pos] |= fill.mask;
  if ((fill.mask & kLeafChannelPdf) != 0U) {
    workspace->source_product_program_pdf[program_pos] = fill.pdf;
  }
  if ((fill.mask & kLeafChannelCdf) != 0U) {
    workspace->source_product_program_cdf[program_pos] = fill.cdf;
  }
  if ((fill.mask & kLeafChannelSurvival) != 0U) {
    workspace->source_product_program_survival[program_pos] = fill.survival;
  }
}

inline double compiled_math_source_product_fill_value(
    const ExactSourceChannels::SourceProductScalarFill &fill,
    const std::uint8_t channel_mask) {
  if (channel_mask == kLeafChannelPdf) {
    return fill.pdf;
  }
  if (channel_mask == kLeafChannelCdf) {
    return fill.cdf;
  }
  if (channel_mask == kLeafChannelSurvival) {
    return fill.survival;
  }
  return 0.0;
}

inline std::uint8_t compiled_math_source_node_fill_mask(
    const CompiledMathNodeKind kind) noexcept {
  const auto value_mask = compiled_math_source_factor_channel_mask(kind);
  if ((value_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    return kLeafChannelCdf | kLeafChannelSurvival;
  }
  return value_mask;
}

inline double compiled_math_source_node_value(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace) {
  auto *source_channels = evaluator->source_channels();
  const auto program_pos =
      static_cast<std::size_t>(node.source_program_id);
  const auto value_mask = compiled_math_source_factor_channel_mask(node.kind);
  const auto fill_mask = compiled_math_source_node_fill_mask(node.kind);
  if (workspace->source_product_program_epoch[program_pos] !=
          workspace->source_product_program_current_epoch ||
      (workspace->source_product_program_valid_mask[program_pos] &
       fill_mask) != fill_mask) {
    const auto fill =
        compiled_math_source_product_program_fill(
            program,
            node.source_program_id,
            workspace,
            source_channels,
            compiled_math_source_node_time(node, *workspace),
            fill_mask);
    compiled_math_store_source_product_program_fill(
        workspace,
        program_pos,
        fill);
  }
  if (value_mask == kLeafChannelPdf) {
    return safe_density(workspace->source_product_program_pdf[program_pos]);
  }
  if (value_mask == kLeafChannelCdf) {
    return clamp_probability(
        workspace->source_product_program_cdf[program_pos]);
  }
  if (value_mask == kLeafChannelSurvival) {
    return clamp_probability(
        workspace->source_product_program_survival[program_pos]);
  }
  return 0.0;
}

inline double compiled_math_source_product_value_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels) {
  double product = 1.0;
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    double value = op.constant_value;
    const auto channel_mask = op.value_channel_mask;
    if (channel_mask == 0U) {
      if (value == 0.0) {
        return 0.0;
      }
    } else {
      const auto channel_pos =
          static_cast<std::size_t>(op.source_product_channel_id);
      const auto &channel =
          program.integral_kernel_source_product_channels[channel_pos];
      if (!op.cache_result) {
        const double time =
            compiled_math_source_product_channel_time(channel, *workspace);
        if (source_channels->source_product_direct_leaf_available(
                op.source_product_channel_id)) {
          value = batch_source_product_direct_leaf_value_for_lane(
              channel,
              source_channels->source_product_leaf_input(channel.leaf_index),
              time,
              channel_mask);
        } else {
          const auto fill =
              compiled_math_source_product_program_fill(
                  program,
                  op.source_product_program_id,
                  workspace,
                  source_channels,
                  time,
                  channel_mask);
          value = compiled_math_source_product_fill_value(fill, channel_mask);
        }
      } else {
        const auto program_pos =
            static_cast<std::size_t>(op.source_product_program_id);
        if (workspace->source_product_program_epoch[program_pos] !=
                workspace->source_product_program_current_epoch ||
            (workspace->source_product_program_valid_mask[program_pos] &
             channel_mask) == 0U) {
          const double time =
              compiled_math_source_product_channel_time(channel, *workspace);
          const auto fill_mask = op.fill_channel_mask;
          const auto fill =
              source_channels->source_product_direct_leaf_available(
                  op.source_product_channel_id)
                  ? compiled_math_source_product_direct_leaf_fill(
                        channel,
                        source_channels->source_product_leaf_input(
                            channel.leaf_index),
                        time,
                        fill_mask)
                  : compiled_math_source_product_program_fill(
                        program,
                        op.source_product_program_id,
                        workspace,
                        source_channels,
                        time,
                        fill_mask);
          compiled_math_store_source_product_program_fill(
              workspace,
              program_pos,
              fill);
        }
        if (channel_mask == kLeafChannelPdf) {
          value = workspace->source_product_program_pdf[program_pos];
        } else if (channel_mask == kLeafChannelCdf) {
          value = workspace->source_product_program_cdf[program_pos];
        } else {
          value = workspace->source_product_program_survival[program_pos];
        }
      }
    }
    product *= value;
    if (product == 0.0) {
      return 0.0;
    }
  }
  return product;
}

inline double evaluate_compiled_source_product_integral_kernel(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *eval_workspace,
    const double lower,
    const double upper) {
  auto *source_channels = parent == nullptr ? nullptr : parent->source_channels();
  if (!(upper > lower)) {
    return 0.0;
  }
  if (source_channels == nullptr) {
    return 0.0;
  }
  auto &integral_workspace = workspace->integral_workspace_for(program);
  const auto batch = quadrature::build_finite_batch(lower, upper);
  const std::vector<std::uint8_t> *term_outcome_open_ptr = nullptr;
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProductSum) {
    integral_workspace.integral_term_open.assign(
        static_cast<std::size_t>(kernel.source_product_terms.size), 1U);
    for (semantic::Index term_idx = 0;
         term_idx < kernel.source_product_terms.size;
         ++term_idx) {
      const auto &term =
          program.integral_kernel_source_product_terms[
              static_cast<std::size_t>(
                  kernel.source_product_terms.offset + term_idx)];
      integral_workspace.integral_term_open[static_cast<std::size_t>(term_idx)] =
          compiled_math_integral_outcome_gates_open(
              program,
              term.outcome_gate_nodes,
              integral_workspace,
              parent)
              ? 1U
              : 0U;
    }
    term_outcome_open_ptr = &integral_workspace.integral_term_open;
  }
  double sum = 0.0;
  CompiledMathWorkspace::RebindableTimeBinding bind(
      &integral_workspace,
      kernel.bind_time_id);
  for (std::size_t sample_idx = 0;
       sample_idx < quadrature::kDefaultFiniteOrder;
       ++sample_idx) {
    bind.set(batch.nodes.nodes[sample_idx]);
    integral_workspace.reset_source_product_program_cache();
    double value = 0.0;
    if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
      value = compiled_math_source_product_value_for_ops(
          program,
          kernel.source_product_ops,
          &integral_workspace,
          source_channels);
      if (batch_source_product_debug_enabled()) {
        batch_debug_compare_source_product_ops_to_scalar(
            program,
            kernel.source_product_ops,
            &integral_workspace,
            source_channels,
            value);
      }
    } else {
      double total = 0.0;
      for (semantic::Index term_idx = 0;
           term_idx < kernel.source_product_terms.size;
           ++term_idx) {
        const auto &term =
            program.integral_kernel_source_product_terms[
                static_cast<std::size_t>(
                    kernel.source_product_terms.offset + term_idx)];
        const auto term_pos = static_cast<std::size_t>(term_idx);
        if (term_outcome_open_ptr != nullptr &&
            term_pos < term_outcome_open_ptr->size() &&
            (*term_outcome_open_ptr)[term_pos] == 0U) {
          continue;
        }
        if (!compiled_math_integral_time_gates_open(
                program, term.time_gate_nodes, integral_workspace)) {
          continue;
        }
        double product = term.sign;
        for (semantic::Index i = 0; i < term.integral_factor_nodes.size; ++i) {
          const auto node_id =
              program.integral_kernel_integral_factor_nodes[
                  static_cast<std::size_t>(
                      term.integral_factor_nodes.offset + i)];
          const auto &integral_node =
              program.nodes[static_cast<std::size_t>(node_id)];
          product *= compiled_math_integral_factor_value_for_node(
              program,
              integral_node,
              &integral_workspace,
              parent,
              nullptr,
              eval_workspace);
          if (product == 0.0) {
            product = 0.0;
            break;
          }
        }
        if (product == 0.0) {
          continue;
        }
        for (semantic::Index i = 0; i < term.expr_upper_factors.size; ++i) {
          const auto &factor =
              program.integral_kernel_expr_upper_factors[
                  static_cast<std::size_t>(
                      term.expr_upper_factors.offset + i)];
          if (!compiled_math_apply_expr_upper_factor(
                  program,
                  factor,
                  &integral_workspace,
                  parent,
                  eval_workspace,
                  &product)) {
            break;
          }
        }
        if (product == 0.0) {
          continue;
        }
        const double source_product =
            compiled_math_source_product_value_for_ops(
                program,
                term.source_product_ops,
                &integral_workspace,
                source_channels);
        if (batch_source_product_debug_enabled()) {
          batch_debug_compare_source_product_ops_to_scalar(
              program,
              term.source_product_ops,
              &integral_workspace,
              source_channels,
              source_product);
        }
        product *= source_product;
        if (product == 0.0) {
          continue;
        }
        total += product;
      }
      value = kernel.clean_signed_source_sum ? clean_signed_value(total) : total;
    }
    if (value == 0.0) {
      continue;
    }
    sum += batch.nodes.weights[sample_idx] * value;
  }
  if (batch_finite_integral_debug_enabled()) {
    batch_debug_compare_finite_integral_to_scalar(
        program,
        kernel,
        workspace,
        source_channels,
        lower,
        upper,
        sum);
  }
  return sum;
}

inline double evaluate_compiled_integral_kernel(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace,
    const double lower,
    const double upper) {
  if (node.integral_kernel_slot == semantic::kInvalidIndex) {
    throw std::runtime_error(
        "compiled integral node has no planned integral kernel");
  }
  const auto kernel_pos =
      static_cast<std::size_t>(node.integral_kernel_slot);
  if (kernel_pos >= program.integral_kernels.size()) {
    throw std::runtime_error(
        "compiled integral node points outside integral kernels");
  }
  const auto &kernel = program.integral_kernels[kernel_pos];
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct ||
      kernel.kind == CompiledMathIntegralKernelKind::SourceProductSum) {
    return evaluate_compiled_source_product_integral_kernel(
        program,
        kernel,
        workspace,
        parent,
        eval_workspace,
        lower,
        upper);
  }
  throw std::runtime_error("compiled integral kernel has no planned evaluator");
}

inline double evaluate_compiled_node_span(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan schedule,
    const semantic::Index result_node_id,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace,
    const std::vector<semantic::Index> *schedule_nodes,
    const bool workspace_prepared) {
  (void)scenario;
  if (result_node_id == semantic::kInvalidIndex) {
    return 0.0;
  }
  if (!workspace_prepared) {
    workspace->ensure_size(program);
  }
  const auto &nodes =
      schedule_nodes == nullptr ? program.root_schedule_nodes : *schedule_nodes;
  for (semantic::Index schedule_idx = 0; schedule_idx < schedule.size;
       ++schedule_idx) {
    const auto node_id = nodes[
        static_cast<std::size_t>(schedule.offset + schedule_idx)];
    const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
    auto *condition_evaluator = parent;
    if (node.source_view_id != 0 &&
        node.source_view_id != semantic::kInvalidIndex) {
      if (eval_workspace == nullptr) {
        throw std::runtime_error(
            "compiled node requires a planned source view but no workspace was supplied");
      }
      condition_evaluator =
          eval_workspace->source_view_evaluator(node.source_view_id, parent);
    }
    const bool cacheable = compiled_math_node_cacheable(node.kind);
    double value = 0.0;
    if (condition_evaluator != nullptr &&
        cacheable &&
        compiled_math_load_node_cache_entry(
            program,
            node, condition_evaluator, *workspace, &value)) {
      workspace->values[static_cast<std::size_t>(node.cache_slot)] = value;
      continue;
    }
    switch (node.kind) {
    case CompiledMathNodeKind::Constant:
      value = node.constant;
      break;
    case CompiledMathNodeKind::SourcePdf:
    case CompiledMathNodeKind::SourceCdf:
    case CompiledMathNodeKind::SourceSurvival:
      value = compiled_math_source_node_value(
          program, node, condition_evaluator, workspace);
      break;
    case CompiledMathNodeKind::TimeGate: {
      const auto child_id =
          program.child_nodes[static_cast<std::size_t>(node.children.offset)];
      const double raw =
          workspace->values[static_cast<std::size_t>(child_id)];
      value = workspace->time(node.time_id) >= workspace->time(node.aux_id)
                  ? raw
                  : 0.0;
      break;
    }
    case CompiledMathNodeKind::ExprDensity:
    case CompiledMathNodeKind::ExprCdf:
    case CompiledMathNodeKind::ExprSurvival:
      throw std::runtime_error(
          "compiled math root contains an interpreter expression node");
    case CompiledMathNodeKind::SimpleGuardDensity: {
      if (condition_evaluator == nullptr) {
        value = 0.0;
        break;
      }
      const auto &kernel = parent->plan().expr_kernels[
          static_cast<std::size_t>(node.subject_id)];
      if (!kernel.simple_event_guard || kernel.has_unless) {
        throw std::runtime_error(
            "compiled simple guard density reached a non-simple guard");
      }
      const auto child_id =
          program.child_nodes[static_cast<std::size_t>(node.children.offset)];
      const double raw =
          workspace->values[static_cast<std::size_t>(child_id)];
      if (condition_evaluator->source_order_known_before(
              kernel.guard_blocker_source_id,
              kernel.guard_ref_source_id,
              node.condition_id,
              workspace)) {
        value = 0.0;
        break;
      }
      const auto upper =
          compiled_math_upper_bound_from_terms(program, node, *workspace);
      if (upper.found) {
        const double current_time = compiled_math_node_time(node, *workspace);
        if (!(upper.normalizer > 0.0) ||
            current_time >= upper.time) {
          value = 0.0;
          break;
        }
        value = clean_signed_value(raw / upper.normalizer);
        break;
      }
      value = clean_signed_value(raw);
      break;
    }
    case CompiledMathNodeKind::SimpleGuardCdf: {
      if (condition_evaluator == nullptr) {
        value = 0.0;
        break;
      }
      const auto &kernel = parent->plan().expr_kernels[
          static_cast<std::size_t>(node.subject_id)];
      if (!kernel.simple_event_guard || kernel.has_unless) {
        throw std::runtime_error(
            "compiled simple guard cdf reached a non-simple guard");
      }
      const auto child_id =
          program.child_nodes[static_cast<std::size_t>(node.children.offset)];
      const double raw =
          workspace->values[static_cast<std::size_t>(child_id)];
      if (condition_evaluator->source_order_known_before(
              kernel.guard_blocker_source_id,
              kernel.guard_ref_source_id,
              node.condition_id,
              workspace)) {
        value = 0.0;
        break;
      }
      const auto upper =
          compiled_math_upper_bound_from_terms(program, node, *workspace);
      if (upper.found) {
        const double current_time = compiled_math_node_time(node, *workspace);
        if (!(upper.normalizer > 0.0)) {
          value = 0.0;
          break;
        }
        if (current_time >= upper.time) {
          value = 1.0;
          break;
        }
        value = clamp_probability(raw / upper.normalizer);
        break;
      }
      value = clamp_probability(raw);
      break;
    }
    case CompiledMathNodeKind::UnionKernelDensity:
      value = parent == nullptr
                  ? 0.0
                  : evaluate_compiled_union_kernel_density(
                        program,
                        node,
                        *workspace,
                        parent->plan(),
                        condition_evaluator,
                        false);
      break;
    case CompiledMathNodeKind::UnionKernelMultiSubsetDensity:
      value = parent == nullptr
                  ? 0.0
                  : evaluate_compiled_union_kernel_density(
                        program,
                        node,
                        *workspace,
                        parent->plan(),
                        condition_evaluator,
                        true);
      break;
    case CompiledMathNodeKind::UnionKernelCdf:
      value = evaluate_compiled_union_kernel_cdf(
          program, node, *workspace, condition_evaluator);
      break;
    case CompiledMathNodeKind::UnionKernelMultiSubsetCdf: {
      if (condition_evaluator == nullptr ||
          node.integral_kernel_slot == semantic::kInvalidIndex) {
        value = 0.0;
        break;
      }
      const double current_time = compiled_math_node_time(node, *workspace);
      if (!(current_time > 0.0)) {
        value = 0.0;
        break;
      }
      value = evaluate_compiled_integral_kernel(
          program,
          node,
          workspace,
          parent,
          scenario,
          eval_workspace,
          0.0,
          current_time);
      value = std::isfinite(value) ? value : 0.0;
      break;
    }
    case CompiledMathNodeKind::OutcomeSubsetUnused: {
      value = 1.0;
      if (workspace->used_outcomes != nullptr) {
        const auto offset = static_cast<std::size_t>(node.subject_id);
        const auto size = static_cast<std::size_t>(node.aux_id);
        const auto &indices = parent->plan().compiled_outcome_gate_indices;
        if (offset + size > indices.size()) {
          throw std::runtime_error(
              "compiled outcome-used gate points outside the plan");
        }
        for (std::size_t i = 0; i < size; ++i) {
          const auto outcome_idx = indices[offset + i];
          if (outcome_idx != semantic::kInvalidIndex &&
              static_cast<std::size_t>(outcome_idx) <
                  workspace->used_outcomes->size() &&
              (*workspace->used_outcomes)[
                  static_cast<std::size_t>(outcome_idx)] != 0U) {
            value = 0.0;
            break;
          }
        }
      }
      break;
    }
    case CompiledMathNodeKind::IntegralZeroToCurrent: {
      if (condition_evaluator == nullptr ||
          node.integral_kernel_slot == semantic::kInvalidIndex) {
        value = 0.0;
        break;
      }
      const double current_time = compiled_math_node_time(node, *workspace);
      if (!(current_time > 0.0)) {
        value = 0.0;
        break;
      }
      value = clamp_probability(
          evaluate_compiled_integral_kernel(
              program,
              node,
              workspace,
              parent,
              scenario,
              eval_workspace,
              0.0,
              current_time));
      break;
    }
    case CompiledMathNodeKind::IntegralZeroToCurrentRaw: {
      if (condition_evaluator == nullptr ||
          node.integral_kernel_slot == semantic::kInvalidIndex) {
        value = 0.0;
        break;
      }
      const double current_time = compiled_math_node_time(node, *workspace);
      if (!(current_time > 0.0)) {
        value = 0.0;
        break;
      }
      value = evaluate_compiled_integral_kernel(
          program,
          node,
          workspace,
          parent,
          scenario,
          eval_workspace,
          0.0,
          current_time);
      value = std::isfinite(value) ? clean_signed_value(value) : 0.0;
      break;
    }
    case CompiledMathNodeKind::ExprUpperBoundDensity:
    case CompiledMathNodeKind::ExprUpperBoundCdf: {
      const auto child_id =
          program.child_nodes[static_cast<std::size_t>(node.children.offset)];
      const double raw =
          workspace->values[static_cast<std::size_t>(child_id)];
      const auto best_upper =
          compiled_math_expr_upper_bound_for_node(
              program,
              node,
              *workspace,
              condition_evaluator);
      if (!best_upper.found) {
        value = raw;
        break;
      }
      const double current_time = compiled_math_node_time(node, *workspace);
      if (node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
        if (!(best_upper.normalizer > 0.0)) {
          value = 0.0;
        } else if (current_time >= best_upper.time) {
          value = 1.0;
        } else {
          value = clamp_probability(raw / best_upper.normalizer);
        }
      } else {
        if (!(best_upper.normalizer > 0.0) ||
            current_time >= best_upper.time) {
          value = 0.0;
        } else {
          value = safe_density(raw / best_upper.normalizer);
        }
      }
      break;
    }
    case CompiledMathNodeKind::Product: {
      value = 1.0;
      for (semantic::Index i = 0; i < node.children.size; ++i) {
        const auto child_id = program.child_nodes[
            static_cast<std::size_t>(node.children.offset + i)];
        value *= workspace->values[static_cast<std::size_t>(child_id)];
        if (!std::isfinite(value) || value == 0.0) {
          value = 0.0;
          break;
        }
      }
      break;
    }
    case CompiledMathNodeKind::Sum:
    case CompiledMathNodeKind::CleanSignedSum: {
      double total = 0.0;
      for (semantic::Index i = 0; i < node.children.size; ++i) {
        const auto child_id = program.child_nodes[
            static_cast<std::size_t>(node.children.offset + i)];
        total += workspace->values[static_cast<std::size_t>(child_id)];
      }
      value = node.kind == CompiledMathNodeKind::CleanSignedSum
                  ? clean_signed_value(total)
                  : total;
      break;
    }
    case CompiledMathNodeKind::ClampProbability: {
      const auto child_id =
          program.child_nodes[static_cast<std::size_t>(node.children.offset)];
      value = clamp_probability(
          workspace->values[static_cast<std::size_t>(child_id)]);
      break;
    }
    case CompiledMathNodeKind::Complement: {
      const auto child_id =
          program.child_nodes[static_cast<std::size_t>(node.children.offset)];
      value = clamp_probability(
          1.0 - workspace->values[static_cast<std::size_t>(child_id)]);
      break;
    }
    case CompiledMathNodeKind::Negate: {
      const auto child_id =
          program.child_nodes[static_cast<std::size_t>(node.children.offset)];
      value = clean_signed_value(
          -workspace->values[static_cast<std::size_t>(child_id)]);
      break;
    }
    }
    if (condition_evaluator != nullptr && cacheable) {
      compiled_math_store_node_cache_entry(
          program,
          node,
          condition_evaluator,
          workspace,
          value);
    }
    workspace->values[static_cast<std::size_t>(node.cache_slot)] = value;
  }
  return workspace->values[static_cast<std::size_t>(result_node_id)];
}

inline double evaluate_compiled_math_root(
    const CompiledMathProgram &program,
    const semantic::Index root_id,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace) {
  if (root_id == semantic::kInvalidIndex) {
    return 0.0;
  }
  const auto &root = program.roots[static_cast<std::size_t>(root_id)];
  return evaluate_compiled_node_span(
      program,
      root.schedule,
      root.node_id,
      workspace,
      parent,
      scenario,
      eval_workspace,
      nullptr,
      true);
}

struct CompiledMathBatchLane {
  CompiledMathWorkspace *workspace{nullptr};
  CompiledSourceView *parent{nullptr};
  CompiledSourceView *scenario{nullptr};
  CompiledEvalWorkspace *eval_workspace{nullptr};
};

struct CompiledMathLaneBatchState {
  std::size_t lane_count{0};
  std::size_t time_slot_count{0};
  BatchTimeSlotView time_slots{};
  const std::vector<CompiledMathWorkspace *> *workspaces_by_lane{nullptr};
  const std::vector<ExactSourceChannels *> *source_channels_by_lane{nullptr};
  const std::vector<CompiledSourceView *> *parents_by_lane{nullptr};
  const std::vector<CompiledEvalWorkspace *> *eval_workspaces_by_lane{nullptr};
  const std::vector<const std::vector<std::uint8_t> *> *used_outcomes_by_lane{
      nullptr};

  CompiledMathWorkspace *workspace(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(lane);
    if (workspaces_by_lane == nullptr || pos >= workspaces_by_lane->size()) {
      return nullptr;
    }
    return (*workspaces_by_lane)[pos];
  }

  ExactSourceChannels *source_channels(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(lane);
    if (source_channels_by_lane == nullptr ||
        pos >= source_channels_by_lane->size()) {
      return nullptr;
    }
    return (*source_channels_by_lane)[pos];
  }

  CompiledSourceView *parent(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(lane);
    if (parents_by_lane == nullptr || pos >= parents_by_lane->size()) {
      return nullptr;
    }
    return (*parents_by_lane)[pos];
  }

  CompiledEvalWorkspace *eval_workspace(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(lane);
    if (eval_workspaces_by_lane == nullptr ||
        pos >= eval_workspaces_by_lane->size()) {
      return nullptr;
    }
    return (*eval_workspaces_by_lane)[pos];
  }

  const std::vector<std::uint8_t> *used_outcomes(
      const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(lane);
    if (used_outcomes_by_lane == nullptr ||
        pos >= used_outcomes_by_lane->size()) {
      return nullptr;
    }
    return (*used_outcomes_by_lane)[pos];
  }

  double time(const semantic::Index time_id,
              const semantic::Index lane) const noexcept {
    return time_slots.get(time_id, lane);
  }
};

struct CompiledMathBatchSourceProductScratch {
  std::vector<semantic::Index> lanes_a;
  std::vector<semantic::Index> lanes_b;
  std::vector<semantic::Index> lanes_c;
  std::vector<semantic::Index> lanes_d;
  std::vector<double> times;
  std::vector<double> onset_upper;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> exact;
  std::vector<std::uint8_t> has_exact;
  std::vector<std::uint8_t> mask_a;
  std::vector<double> pdf_a;
  std::vector<double> cdf_a;
  std::vector<double> survival_a;
  std::vector<std::uint8_t> mask_b;
  std::vector<double> pdf_b;
  std::vector<double> cdf_b;
  std::vector<double> survival_b;
  std::vector<std::uint8_t> mask_c;
  std::vector<double> pdf_c;
  std::vector<double> cdf_c;
  std::vector<double> survival_c;
  std::vector<double> pool_member_pdf;
  std::vector<double> pool_member_cdf;
  std::vector<double> pool_member_survival;
  std::vector<double> pool_prefix;
  std::vector<double> pool_suffix;

  void ensure_size(const std::size_t lane_capacity,
                   const std::size_t active_capacity) {
    if (lanes_a.size() < active_capacity) {
      lanes_a.resize(active_capacity);
    }
    if (lanes_b.size() < active_capacity) {
      lanes_b.resize(active_capacity);
    }
    if (lanes_c.size() < active_capacity) {
      lanes_c.resize(active_capacity);
    }
    if (lanes_d.size() < active_capacity) {
      lanes_d.resize(active_capacity);
    }
    if (times.size() < lane_capacity) {
      times.resize(lane_capacity, 0.0);
    }
    if (onset_upper.size() < lane_capacity) {
      onset_upper.resize(lane_capacity, 0.0);
    }
    if (lower.size() < lane_capacity) {
      lower.resize(lane_capacity, 0.0);
    }
    if (upper.size() < lane_capacity) {
      upper.resize(lane_capacity, 0.0);
    }
    if (exact.size() < lane_capacity) {
      exact.resize(lane_capacity, 0.0);
    }
    if (has_exact.size() < lane_capacity) {
      has_exact.resize(lane_capacity, 0U);
    }
    if (mask_a.size() < lane_capacity) {
      mask_a.resize(lane_capacity, 0U);
      pdf_a.resize(lane_capacity, 0.0);
      cdf_a.resize(lane_capacity, 0.0);
      survival_a.resize(lane_capacity, 1.0);
    }
    if (mask_b.size() < lane_capacity) {
      mask_b.resize(lane_capacity, 0U);
      pdf_b.resize(lane_capacity, 0.0);
      cdf_b.resize(lane_capacity, 0.0);
      survival_b.resize(lane_capacity, 1.0);
    }
    if (mask_c.size() < lane_capacity) {
      mask_c.resize(lane_capacity, 0U);
      pdf_c.resize(lane_capacity, 0.0);
      cdf_c.resize(lane_capacity, 0.0);
      survival_c.resize(lane_capacity, 1.0);
    }
  }

  void ensure_pool_size(const std::size_t lane_capacity,
                        const std::size_t member_count) {
    const auto member_value_count = lane_capacity * member_count;
    if (pool_member_pdf.size() < member_value_count) {
      pool_member_pdf.resize(member_value_count, 0.0);
      pool_member_cdf.resize(member_value_count, 0.0);
      pool_member_survival.resize(member_value_count, 1.0);
    }
    const auto width = member_count + 1U;
    const auto table_value_count = lane_capacity * width * width;
    if (pool_prefix.size() < table_value_count) {
      pool_prefix.resize(table_value_count, 0.0);
      pool_suffix.resize(table_value_count, 0.0);
    }
  }
};

struct CompiledMathBatchIntegralFrame {
  std::vector<double> time_values;
  std::vector<std::uint8_t> time_valid;
  std::vector<CompiledMathWorkspace *> workspaces_by_lane;
  std::vector<ExactSourceChannels *> source_channels_by_lane;
  std::vector<CompiledSourceView *> parents_by_lane;
  std::vector<CompiledEvalWorkspace *> eval_workspaces_by_lane;
  std::vector<const std::vector<std::uint8_t> *> used_outcomes_by_lane;
  std::vector<semantic::Index> active_lanes;
  std::vector<semantic::Index> parent_lanes;
  std::vector<semantic::Index> term_lanes_a;
  std::vector<semantic::Index> term_lanes_b;
  std::vector<semantic::Index> source_lanes_a;
  std::vector<semantic::Index> source_lanes_b;
  std::vector<double> weights;
  std::vector<double> values;
  std::vector<double> products;
  std::vector<double> factor_values;
  std::vector<double> source_values;
  std::vector<double> source_leaf_values;
  std::vector<double> bind_times;
  std::vector<double> lower_by_lane;
  std::vector<double> upper_by_lane;
  std::vector<std::uint8_t> source_program_valid_mask;
  std::vector<double> source_program_pdf;
  std::vector<double> source_program_cdf;
  std::vector<double> source_program_survival;
  std::vector<semantic::Index> source_program_cache_slot_by_program;
  std::vector<semantic::Index> cached_source_programs;
  std::vector<CompiledMathBatchSourceProductScratch>
      source_product_batch_scratch;

  void ensure_lane_capacity(const std::size_t time_slot_count,
                            const std::size_t lane_capacity) {
    const auto time_value_count = time_slot_count * lane_capacity;
    if (time_values.size() < time_value_count) {
      time_values.resize(time_value_count, 0.0);
    }
    if (time_valid.size() < time_value_count) {
      time_valid.resize(time_value_count, 0U);
    }
    if (workspaces_by_lane.size() < lane_capacity) {
      workspaces_by_lane.resize(lane_capacity, nullptr);
    }
    if (source_channels_by_lane.size() < lane_capacity) {
      source_channels_by_lane.resize(lane_capacity, nullptr);
    }
    if (parents_by_lane.size() < lane_capacity) {
      parents_by_lane.resize(lane_capacity, nullptr);
    }
    if (eval_workspaces_by_lane.size() < lane_capacity) {
      eval_workspaces_by_lane.resize(lane_capacity, nullptr);
    }
    if (used_outcomes_by_lane.size() < lane_capacity) {
      used_outcomes_by_lane.resize(lane_capacity, nullptr);
    }
    if (active_lanes.size() < lane_capacity) {
      active_lanes.resize(lane_capacity);
    }
    if (parent_lanes.size() < lane_capacity) {
      parent_lanes.resize(lane_capacity);
    }
    if (term_lanes_a.size() < lane_capacity) {
      term_lanes_a.resize(lane_capacity);
    }
    if (term_lanes_b.size() < lane_capacity) {
      term_lanes_b.resize(lane_capacity);
    }
    if (source_lanes_a.size() < lane_capacity) {
      source_lanes_a.resize(lane_capacity);
    }
    if (source_lanes_b.size() < lane_capacity) {
      source_lanes_b.resize(lane_capacity);
    }
    if (weights.size() < lane_capacity) {
      weights.resize(lane_capacity, 0.0);
    }
    if (values.size() < lane_capacity) {
      values.resize(lane_capacity, 0.0);
    }
    if (products.size() < lane_capacity) {
      products.resize(lane_capacity, 0.0);
    }
    if (factor_values.size() < lane_capacity) {
      factor_values.resize(lane_capacity, 0.0);
    }
    if (source_values.size() < lane_capacity) {
      source_values.resize(lane_capacity, 0.0);
    }
    if (source_leaf_values.size() < lane_capacity) {
      source_leaf_values.resize(lane_capacity, 0.0);
    }
    if (bind_times.size() < lane_capacity) {
      bind_times.resize(lane_capacity, 0.0);
    }
    if (lower_by_lane.size() < lane_capacity) {
      lower_by_lane.resize(lane_capacity, 0.0);
    }
    if (upper_by_lane.size() < lane_capacity) {
      upper_by_lane.resize(lane_capacity, 0.0);
    }
  }

  void reset_source_product_cache(
      const std::size_t program_count,
      const std::size_t lane_capacity,
      const std::vector<semantic::Index> &program_ids) {
    source_program_cache_slot_by_program.assign(
        program_count,
        semantic::kInvalidIndex);
    std::size_t cache_slot_count = 0U;
    for (const auto program_id : program_ids) {
      if (program_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(program_id) >= program_count) {
        continue;
      }
      auto &slot =
          source_program_cache_slot_by_program[
              static_cast<std::size_t>(program_id)];
      if (slot != semantic::kInvalidIndex) {
        continue;
      }
      slot = static_cast<semantic::Index>(cache_slot_count++);
    }
    const auto count = cache_slot_count * lane_capacity;
    if (source_program_valid_mask.size() < count) {
      source_program_valid_mask.resize(count, 0U);
      source_program_pdf.resize(count, 0.0);
      source_program_cdf.resize(count, 0.0);
      source_program_survival.resize(count, 1.0);
    }
    std::fill(
        source_program_valid_mask.begin(),
        source_program_valid_mask.begin() + static_cast<std::ptrdiff_t>(count),
        0U);
  }

  semantic::Index source_program_cache_slot(
      const semantic::Index program_id) const noexcept {
    if (program_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(program_id) >=
            source_program_cache_slot_by_program.size()) {
      return semantic::kInvalidIndex;
    }
    return source_program_cache_slot_by_program[
        static_cast<std::size_t>(program_id)];
  }

  CompiledMathBatchSourceProductScratch &source_product_scratch_layer(
      const std::size_t depth,
      const std::size_t lane_capacity,
      const std::size_t active_capacity) {
    constexpr std::size_t kMaxSourceProductScratchDepth = 32U;
    if (source_product_batch_scratch.size() < kMaxSourceProductScratchDepth) {
      source_product_batch_scratch.resize(kMaxSourceProductScratchDepth);
    } else if (source_product_batch_scratch.size() <= depth) {
      source_product_batch_scratch.resize(depth + 1U);
    }
    auto &scratch = source_product_batch_scratch[depth];
    scratch.ensure_size(lane_capacity, active_capacity);
    return scratch;
  }
};

struct CompiledMathBatchWorkspace {
  std::vector<CompiledSourceView *> condition_evaluators;
  std::vector<semantic::Index> active_lanes;
  std::vector<std::size_t> active_lane_indices;
  std::vector<std::size_t> compact_lane_indices;
  std::vector<CompiledMathWorkspace *> workspaces_by_lane;
  std::vector<ExactSourceChannels *> source_channels_by_lane;
  std::vector<CompiledSourceView *> parents_by_lane;
  std::vector<CompiledEvalWorkspace *> eval_workspaces_by_lane;
  std::vector<const std::vector<std::uint8_t> *> used_outcomes_by_lane;
  std::vector<double> parent_time_values;
  std::vector<std::uint8_t> parent_time_valid;
  std::vector<ExactLoadedLeafInput> parent_leaf_inputs;
  std::vector<double> lower_by_lane;
  std::vector<double> upper_by_lane;
  std::vector<double> batch_values;
  BatchFiniteIntegralWorkspace finite_integral_workspace;
  std::vector<CompiledMathBatchIntegralFrame> integral_frames;
};

inline CompiledMathBatchIntegralFrame &compiled_math_batch_integral_frame(
    CompiledMathBatchWorkspace *workspace,
    const std::size_t depth) {
  constexpr std::size_t kMaxBatchIntegralDepth = 18U;
  if (workspace->integral_frames.size() < kMaxBatchIntegralDepth) {
    workspace->integral_frames.resize(kMaxBatchIntegralDepth);
  } else if (workspace->integral_frames.size() <= depth) {
    workspace->integral_frames.resize(depth + 1U);
  }
  return workspace->integral_frames[depth];
}

inline std::size_t compiled_math_program_time_slot_count(
    const CompiledMathProgram &program) {
  semantic::Index max_time = static_cast<semantic::Index>(3);
  auto include_time = [&max_time](const semantic::Index time_id) {
    if (time_id != semantic::kInvalidIndex && time_id > max_time) {
      max_time = time_id;
    }
  };
  for (const auto &node : program.nodes) {
    include_time(node.time_id);
    if (node.kind == CompiledMathNodeKind::TimeGate ||
        compiled_math_is_source_value_node(node.kind)) {
      include_time(node.aux_id);
    }
  }
  for (const auto &kernel : program.integral_kernels) {
    include_time(kernel.bind_time_id);
  }
  for (const auto &channel : program.integral_kernel_source_product_channels) {
    include_time(channel.time_id);
    include_time(channel.time_cap_id);
  }
  for (const auto &term : program.source_condition_bound_terms) {
    include_time(term.time_id);
  }
  for (const auto &term : program.timed_upper_bound_terms) {
    include_time(term.time_id);
  }
  return static_cast<std::size_t>(max_time) + 1U;
}

inline bool compiled_math_batch_source_program_supported(
    const CompiledMathProgram &program,
    semantic::Index source_product_program_id,
    std::size_t depth);

inline bool compiled_math_batch_kernel_supported(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    std::size_t depth);

inline bool compiled_math_batch_source_product_ops_supported(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    const std::size_t depth) {
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.value_channel_mask == 0U) {
      continue;
    }
    if (op.source_product_program_id == semantic::kInvalidIndex ||
        !compiled_math_batch_source_program_supported(
            program,
            op.source_product_program_id,
            depth + 1U)) {
      return false;
    }
  }
  return true;
}

inline bool compiled_math_batch_source_program_supported(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    const std::size_t depth) {
  if (source_product_program_id == semantic::kInvalidIndex) {
    return true;
  }
  if (depth > 32U) {
    return false;
  }
  const auto pos = static_cast<std::size_t>(source_product_program_id);
  if (pos >= program.integral_kernel_source_product_programs.size()) {
    return false;
  }
  const auto &source_program =
      program.integral_kernel_source_product_programs[pos];
  switch (source_program.kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
  case CompiledMathSourceProductProgramKind::ConstantOne:
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    return true;
  case CompiledMathSourceProductProgramKind::ExactGate:
  case CompiledMathSourceProductProgramKind::Conditioned:
    return compiled_math_batch_source_program_supported(
        program,
        source_program.child_program_id,
        depth + 1U);
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
    return source_program.leaf_index != semantic::kInvalidIndex &&
           compiled_math_batch_source_program_supported(
               program,
               source_program.onset_source_program_id,
               depth + 1U);
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    if (source_program.member_programs.empty()) {
      return false;
    }
    for (semantic::Index i = 0; i < source_program.member_programs.size; ++i) {
      const auto member_program_id =
          program.integral_kernel_source_product_program_members[
              static_cast<std::size_t>(
                  source_program.member_programs.offset + i)];
      if (!compiled_math_batch_source_program_supported(
              program,
              member_program_id,
              depth + 1U)) {
        return false;
      }
    }
    return true;
  }
  return false;
}

inline bool compiled_math_batch_expr_upper_supported(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan expr_upper_factors) {
  for (semantic::Index i = 0; i < expr_upper_factors.size; ++i) {
    const auto &factor =
        program.integral_kernel_expr_upper_factors[
            static_cast<std::size_t>(expr_upper_factors.offset + i)];
    const auto &node =
        program.nodes[static_cast<std::size_t>(factor.node_id)];
    if (node.kind != CompiledMathNodeKind::ExprUpperBoundDensity &&
        node.kind != CompiledMathNodeKind::ExprUpperBoundCdf) {
      return false;
    }
    if (node.aux_id == semantic::kInvalidIndex ||
        node.aux2_id == semantic::kInvalidIndex ||
        node.aux2_id == 0) {
      continue;
    }
    const auto offset = static_cast<std::size_t>(node.aux_id);
    const auto size = static_cast<std::size_t>(node.aux2_id);
    if (offset + size > program.timed_upper_bound_terms.size()) {
      return false;
    }
    for (std::size_t term_pos = 0; term_pos < size; ++term_pos) {
      const auto &term =
          program.timed_upper_bound_terms[offset + term_pos];
      if (term.normalizer_node_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(term.normalizer_node_id) >=
              program.nodes.size()) {
        return false;
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_integral_node_supported(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const std::size_t depth) {
  if (!compiled_math_is_integral_node(node.kind) ||
      node.integral_kernel_slot == semantic::kInvalidIndex ||
      depth > 16U) {
    return false;
  }
  const auto kernel_pos = static_cast<std::size_t>(node.integral_kernel_slot);
  if (kernel_pos >= program.integral_kernels.size()) {
    return false;
  }
  if (!compiled_math_batch_kernel_supported(
          program,
          program.integral_kernels[kernel_pos],
          depth + 1U)) {
    return false;
  }
  return true;
}

inline bool compiled_math_batch_kernel_supported(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const std::size_t depth) {
  if (kernel.bind_time_id == semantic::kInvalidIndex || depth > 16U) {
    return false;
  }
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
    return compiled_math_batch_source_product_ops_supported(
        program,
        kernel.source_product_ops,
        depth + 1U);
  }
  if (kernel.kind != CompiledMathIntegralKernelKind::SourceProductSum) {
    return false;
  }
  for (semantic::Index term_idx = 0;
       term_idx < kernel.source_product_terms.size;
       ++term_idx) {
    const auto &term =
        program.integral_kernel_source_product_terms[
            static_cast<std::size_t>(
                kernel.source_product_terms.offset + term_idx)];
    if (!compiled_math_batch_expr_upper_supported(
            program,
            term.expr_upper_factors) ||
        !compiled_math_batch_source_product_ops_supported(
            program,
            term.source_product_ops,
            depth + 1U)) {
      return false;
    }
    for (semantic::Index i = 0; i < term.integral_factor_nodes.size; ++i) {
      const auto node_id =
          program.integral_kernel_integral_factor_nodes[
              static_cast<std::size_t>(
                  term.integral_factor_nodes.offset + i)];
      if (!compiled_math_batch_integral_node_supported(
              program,
              program.nodes[static_cast<std::size_t>(node_id)],
              depth + 1U)) {
        return false;
      }
    }
  }
  return true;
}

inline bool compiled_math_source_product_ops_need_batch_cache(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops) {
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.cache_result) {
      return true;
    }
  }
  return false;
}

inline void compiled_math_add_cached_source_program_id(
    const semantic::Index program_id,
    std::vector<semantic::Index> *program_ids) {
  if (program_id == semantic::kInvalidIndex || program_ids == nullptr) {
    return;
  }
  for (const auto existing : *program_ids) {
    if (existing == program_id) {
      return;
    }
  }
  program_ids->push_back(program_id);
}

inline void compiled_math_collect_cached_source_programs_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    std::vector<semantic::Index> *program_ids) {
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.cache_result) {
      compiled_math_add_cached_source_program_id(
          op.source_product_program_id,
          program_ids);
    }
  }
}

inline void compiled_math_collect_kernel_cached_source_programs(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    std::vector<semantic::Index> *program_ids) {
  if (program_ids == nullptr) {
    return;
  }
  program_ids->clear();
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
    compiled_math_collect_cached_source_programs_for_ops(
        program,
        kernel.source_product_ops,
        program_ids);
    return;
  }
  if (kernel.kind != CompiledMathIntegralKernelKind::SourceProductSum) {
    return;
  }
  for (semantic::Index term_idx = 0;
       term_idx < kernel.source_product_terms.size;
       ++term_idx) {
    const auto &term =
        program.integral_kernel_source_product_terms[
            static_cast<std::size_t>(
                kernel.source_product_terms.offset + term_idx)];
    compiled_math_collect_cached_source_programs_for_ops(
        program,
        term.source_product_ops,
        program_ids);
  }
}

inline bool compiled_math_kernel_needs_batch_source_cache(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel) {
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
    return compiled_math_source_product_ops_need_batch_cache(
        program,
        kernel.source_product_ops);
  }
  if (kernel.kind != CompiledMathIntegralKernelKind::SourceProductSum) {
    return false;
  }
  for (semantic::Index term_idx = 0;
       term_idx < kernel.source_product_terms.size;
       ++term_idx) {
    const auto &term =
        program.integral_kernel_source_product_terms[
            static_cast<std::size_t>(
                kernel.source_product_terms.offset + term_idx)];
    if (compiled_math_source_product_ops_need_batch_cache(
            program,
            term.source_product_ops)) {
      return true;
    }
  }
  return false;
}

inline void compiled_math_batch_write_impossible_fill_lane(
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  mask_by_lane[lane_pos] = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    pdf_by_lane[lane_pos] = 0.0;
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    cdf_by_lane[lane_pos] = 0.0;
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[lane_pos] = 1.0;
  }
}

inline void compiled_math_batch_write_certain_fill_lane(
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  mask_by_lane[lane_pos] = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    pdf_by_lane[lane_pos] = 0.0;
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    cdf_by_lane[lane_pos] = 1.0;
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[lane_pos] = 0.0;
  }
}

inline void compiled_math_batch_write_forced_fill_lane(
    const ExactRelation relation,
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  if (relation == ExactRelation::Before || relation == ExactRelation::At) {
    compiled_math_batch_write_certain_fill_lane(
        lane, fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
        survival_by_lane);
    return;
  }
  compiled_math_batch_write_impossible_fill_lane(
      lane, fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
      survival_by_lane);
  if (relation == ExactRelation::After &&
      (fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[static_cast<std::size_t>(lane)] = 1.0;
  }
}

inline void compiled_math_batch_write_fill_for_lanes(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const std::uint8_t fill_mask,
    const bool certain,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  for (std::size_t i = 0; i < lane_count; ++i) {
    if (certain) {
      compiled_math_batch_write_certain_fill_lane(
          lanes[i], fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
          survival_by_lane);
    } else {
      compiled_math_batch_write_impossible_fill_lane(
          lanes[i], fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
          survival_by_lane);
    }
  }
}

inline void compiled_math_batch_copy_fill_lane(
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    const std::uint8_t *src_mask_by_lane,
    const double *src_pdf_by_lane,
    const double *src_cdf_by_lane,
    const double *src_survival_by_lane,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  mask_by_lane[lane_pos] = src_mask_by_lane[lane_pos] | fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    pdf_by_lane[lane_pos] = src_pdf_by_lane[lane_pos];
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    cdf_by_lane[lane_pos] = src_cdf_by_lane[lane_pos];
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[lane_pos] = src_survival_by_lane[lane_pos];
  }
}

template <typename LeafValueFn>
inline bool compiled_math_batch_program_leaf_values_for_kind(
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    LeafValueFn leaf_value_fn) {
  if (source_program.leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr) {
      return false;
    }
    const auto lane_pos = static_cast<std::size_t>(lane);
    const auto &loaded =
        source_channels->source_product_leaf_input(source_program.leaf_index);
    const double x =
        current_time_by_lane[lane_pos] -
        source_program.leaf_onset_abs_value -
        loaded.t0;
    values_by_lane[lane_pos] = leaf_value_fn(loaded, x, channel_mask);
  }
  return true;
}

template <typename LeafValueFn>
inline bool compiled_math_batch_leaf_values_for_kind_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    LeafValueFn leaf_value_fn) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr) {
      return false;
    }
    const auto lane_pos = static_cast<std::size_t>(lane);
    const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
    const double x =
        current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded.t0;
    values_by_lane[lane_pos] = leaf_value_fn(loaded, x, channel_mask);
  }
  return true;
}

inline bool compiled_math_batch_program_lognormal_leaf_values(
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (source_program.leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const auto &loaded =
          source_channels->source_product_leaf_input(source_program.leaf_index);
      const double x =
          current_time_by_lane[lane_pos] -
          source_program.leaf_onset_abs_value -
          loaded.t0;
      const double m = loaded.params[0];
      const double s = loaded.params[1];
      if (!(x > 0.0) || !std::isfinite(m) || !std::isfinite(s) ||
          s <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded.q,
                      channel_mask);
        continue;
      }
      scratch->lanes_a[valid_count] = lane;
      scratch->lower[valid_count] = x;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(source_program.leaf_index);
        const double z =
            (scratch->upper[i] - loaded.params[0]) / loaded.params[1];
        scratch->upper[i] = -0.5 * z * z;
      }
      vvexp(scratch->upper.data(), scratch->upper.data(), &n);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(source_program.leaf_index);
        const double base_pdf =
            kInvSqrtTwoPi * scratch->upper[i] /
            (scratch->lower[i] * loaded.params[1]);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                0.0,
                loaded.q,
                channel_mask);
      }
    } else {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(source_program.leaf_index);
        scratch->upper[i] =
            (scratch->upper[i] - loaded.params[0]) / loaded.params[1];
      }
      compiled_math_batch_normal_cdf_from_z(
          scratch->upper.data(),
          scratch->lower.data(),
          valid_count);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(source_program.leaf_index);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                scratch->upper[i],
                loaded.q,
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return compiled_math_batch_program_leaf_values_for_kind(
      source_program,
      state,
      lanes,
      lane_count,
      current_time_by_lane,
      channel_mask,
      values_by_lane,
      batch_source_product_lognormal_leaf_value);
}

inline bool compiled_math_batch_lognormal_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded.t0;
      const double m = loaded.params[0];
      const double s = loaded.params[1];
      if (!(x > 0.0) || !std::isfinite(m) || !std::isfinite(s) ||
          s <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded.q,
                      channel_mask);
        continue;
      }
      scratch->lanes_a[valid_count] = lane;
      scratch->lower[valid_count] = x;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        auto *source_channels = state.source_channels(lane);
        const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
        const double z =
            (scratch->upper[i] - loaded.params[0]) / loaded.params[1];
        scratch->upper[i] = -0.5 * z * z;
      }
      vvexp(scratch->upper.data(), scratch->upper.data(), &n);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        auto *source_channels = state.source_channels(lane);
        const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
        const double base_pdf =
            kInvSqrtTwoPi * scratch->upper[i] /
            (scratch->lower[i] * loaded.params[1]);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                0.0,
                loaded.q,
                channel_mask);
      }
    } else {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        auto *source_channels = state.source_channels(lane);
        const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
        scratch->upper[i] =
            (scratch->upper[i] - loaded.params[0]) / loaded.params[1];
      }
      compiled_math_batch_normal_cdf_from_z(
          scratch->upper.data(),
          scratch->lower.data(),
          valid_count);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        auto *source_channels = state.source_channels(lane);
        const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                scratch->upper[i],
                loaded.q,
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return compiled_math_batch_leaf_values_for_kind_from_times(
      leaf_index,
      leaf_onset_abs_value,
      state,
      lanes,
      lane_count,
      current_time_by_lane,
      channel_mask,
      values_by_lane,
      batch_source_product_lognormal_leaf_value);
}

inline bool compiled_math_batch_gamma_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr && channel_mask == kLeafChannelPdf) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded.t0;
      const double shape = loaded.params[0];
      const double rate = loaded.params[1];
      if (!(x > 0.0) || !std::isfinite(shape) || shape <= 0.0 ||
          !std::isfinite(rate) || rate <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded.q,
                      channel_mask);
        continue;
      }
      scratch->lanes_a[valid_count] = lane;
      scratch->lower[valid_count] = x;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    for (std::size_t i = 0; i < valid_count; ++i) {
      const auto lane = scratch->lanes_a[i];
      auto *source_channels = state.source_channels(lane);
      const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
      const double shape = loaded.params[0];
      const double rate = loaded.params[1];
      scratch->upper[i] =
          shape * std::log(rate) + (shape - 1.0) * scratch->upper[i] -
          rate * scratch->lower[i] - std::lgamma(shape);
    }
    vvexp(scratch->upper.data(), scratch->upper.data(), &n);
    for (std::size_t i = 0; i < valid_count; ++i) {
      const auto lane = scratch->lanes_a[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      auto *source_channels = state.source_channels(lane);
      const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
      values_by_lane[lane_pos] =
          batch_source_product_finish_base_value(
              scratch->upper[i],
              0.0,
              loaded.q,
              channel_mask);
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return compiled_math_batch_leaf_values_for_kind_from_times(
      leaf_index,
      leaf_onset_abs_value,
      state,
      lanes,
      lane_count,
      current_time_by_lane,
      channel_mask,
      values_by_lane,
      batch_source_product_gamma_leaf_value);
}

inline bool compiled_math_batch_exgauss_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded.t0;
      const double mu = loaded.params[0];
      const double sigma = loaded.params[1];
      const double tau = loaded.params[2];
      if (!(x > 0.0) || !std::isfinite(mu) || !std::isfinite(sigma) ||
          sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded.q,
                      channel_mask);
        continue;
      }
      const double inv_tau = 1.0 / tau;
      const double sigma_over_tau = sigma * inv_tau;
      const double z0 = -mu / sigma;
      scratch->lanes_a[valid_count] = lane;
      scratch->pdf_a[valid_count] = x;
      scratch->upper[valid_count] = z0;
      scratch->exact[valid_count] = z0 - sigma_over_tau;
      scratch->lower[valid_count] =
          sigma * sigma * inv_tau * inv_tau * 0.5 + mu * inv_tau;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    const auto n = static_cast<int>(valid_count);
    vvexp(scratch->lower.data(), scratch->lower.data(), &n);
    compiled_math_batch_normal_cdf_from_z(
        scratch->upper.data(),
        scratch->cdf_a.data(),
        valid_count);
    compiled_math_batch_normal_cdf_from_z(
        scratch->exact.data(),
        scratch->cdf_a.data(),
        valid_count);
    std::size_t live_count = 0U;
    for (std::size_t i = 0; i < valid_count; ++i) {
      const auto lane = scratch->lanes_a[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      auto *source_channels = state.source_channels(lane);
      const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
      const double lower_cdf =
          clamp_probability(scratch->upper[i] - scratch->lower[i] * scratch->exact[i]);
      const double lower_survival = 1.0 - lower_cdf;
      if (!(lower_survival > 0.0)) {
        values_by_lane[lane_pos] =
            batch_source_product_before_onset_value(channel_mask);
        continue;
      }
      const double mu = loaded.params[0];
      const double sigma = loaded.params[1];
      const double tau = loaded.params[2];
      const double inv_tau = 1.0 / tau;
      const double sigma_over_tau = sigma * inv_tau;
      const double x = scratch->pdf_a[i];
      const double z = (x - mu) / sigma;
      scratch->lanes_b[live_count] = lane;
      scratch->pdf_b[live_count] = x;
      scratch->cdf_b[live_count] = lower_cdf;
      scratch->survival_b[live_count] = lower_survival;
      scratch->lower[live_count] =
          sigma * sigma * inv_tau * inv_tau * 0.5 -
          (x - mu) * inv_tau;
      scratch->exact[live_count] = z - sigma_over_tau;
      if (channel_mask != kLeafChannelPdf) {
        scratch->upper[live_count] = z;
      }
      ++live_count;
    }
    if (live_count == 0U) {
      return true;
    }
    const auto live_n = static_cast<int>(live_count);
    vvexp(scratch->lower.data(), scratch->lower.data(), &live_n);
    compiled_math_batch_normal_cdf_from_z(
        scratch->exact.data(),
        scratch->cdf_a.data(),
        live_count);
    if (channel_mask != kLeafChannelPdf) {
      compiled_math_batch_normal_cdf_from_z(
          scratch->upper.data(),
          scratch->cdf_a.data(),
          live_count);
    }
    for (std::size_t i = 0; i < live_count; ++i) {
      const auto lane = scratch->lanes_b[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      auto *source_channels = state.source_channels(lane);
      const auto &loaded = source_channels->source_product_leaf_input(leaf_index);
      if (channel_mask == kLeafChannelPdf) {
        const double base_pdf =
            (scratch->lower[i] * scratch->exact[i] / loaded.params[2]) /
            scratch->survival_b[i];
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                0.0,
                loaded.q,
                channel_mask);
      } else {
        const double raw_cdf =
            clamp_probability(scratch->upper[i] -
                              scratch->lower[i] * scratch->exact[i]);
        const double base_cdf =
            (raw_cdf - scratch->cdf_b[i]) / scratch->survival_b[i];
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                base_cdf,
                loaded.q,
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return compiled_math_batch_leaf_values_for_kind_from_times(
      leaf_index,
      leaf_onset_abs_value,
      state,
      lanes,
      lane_count,
      current_time_by_lane,
      channel_mask,
      values_by_lane,
      batch_source_product_exgauss_leaf_value);
}

inline bool compiled_math_batch_leaf_values_from_times(
    const std::uint8_t leaf_dist_kind,
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  switch (static_cast<leaf::DistKind>(leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    return compiled_math_batch_lognormal_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  case leaf::DistKind::Gamma:
    return compiled_math_batch_gamma_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  case leaf::DistKind::Exgauss:
    return compiled_math_batch_exgauss_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  case leaf::DistKind::LBA:
    return compiled_math_batch_leaf_values_for_kind_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        batch_source_product_lba_leaf_value);
  case leaf::DistKind::RDM:
    return compiled_math_batch_leaf_values_for_kind_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        batch_source_product_rdm_leaf_value);
  }
  return false;
}

inline bool compiled_math_batch_program_leaf_values(
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  return compiled_math_batch_leaf_values_from_times(
      source_program.leaf_dist_kind,
      source_program.leaf_index,
      source_program.leaf_onset_abs_value,
      state,
      lanes,
      lane_count,
      current_time_by_lane,
      channel_mask,
      values_by_lane,
      scratch);
}

inline bool compiled_math_batch_source_product_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane);

inline bool compiled_math_batch_source_product_program_leaf_fill(
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  auto *scratch =
      frame == nullptr
          ? nullptr
          : &frame->source_product_scratch_layer(
                scratch_depth,
                state.lane_count,
                lane_count);
  for (std::size_t i = 0; i < lane_count; ++i) {
    mask_by_lane[static_cast<std::size_t>(lanes[i])] = fill_mask;
  }
  if ((fill_mask & kLeafChannelPdf) != 0U &&
      !compiled_math_batch_program_leaf_values(
          source_program,
          state,
          lanes,
          lane_count,
          current_time_by_lane,
          kLeafChannelPdf,
          pdf_by_lane,
          scratch)) {
    return false;
  }
  if ((fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    if (!compiled_math_batch_program_leaf_values(
            source_program,
            state,
            lanes,
            lane_count,
            current_time_by_lane,
            kLeafChannelCdf,
            cdf_by_lane,
            scratch)) {
      return false;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      for (std::size_t i = 0; i < lane_count; ++i) {
        const auto lane_pos = static_cast<std::size_t>(lanes[i]);
        survival_by_lane[lane_pos] =
            clamp_probability(1.0 - cdf_by_lane[lane_pos]);
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_resolve_source_bounds_for_source(
    const semantic::Index source_id,
    const CompiledSourceBoundPlan &bounds,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *lower_by_lane,
    double *upper_by_lane,
    double *exact_by_lane,
    std::uint8_t *has_exact_by_lane) {
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr) {
      return false;
    }
    const auto lane_pos = static_cast<std::size_t>(lane);
    source_channels->source_product_resolved_bound_values_from_time_slots(
        source_id,
        bounds,
        state.time_slots.values,
        state.time_slots.lane_stride,
        lane,
        state.time_slots.valid,
        &lower_by_lane[lane_pos],
        &upper_by_lane[lane_pos],
        &exact_by_lane[lane_pos],
        &has_exact_by_lane[lane_pos]);
  }
  return true;
}

inline bool compiled_math_batch_resolve_source_bounds(
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *lower_by_lane,
    double *upper_by_lane,
    double *exact_by_lane,
    std::uint8_t *has_exact_by_lane) {
  return compiled_math_batch_resolve_source_bounds_for_source(
      source_program.source_id,
      source_program.bounds,
      state,
      lanes,
      lane_count,
      lower_by_lane,
      upper_by_lane,
      exact_by_lane,
      has_exact_by_lane);
}

inline bool compiled_math_batch_source_product_program_exact_gate_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !compiled_math_source_product_bounds_have_overlay(
          source_program.bounds)) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      compiled_math_batch_write_forced_fill_lane(
          relation,
          lanes[i],
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
    }
    return true;
  }
  auto &scratch = frame->source_product_scratch_layer(
      scratch_depth, state.lane_count, lane_count);
  if (!compiled_math_batch_resolve_source_bounds(
          source_program,
          state,
          lanes,
          lane_count,
          scratch.lower.data(),
          scratch.upper.data(),
          scratch.exact.data(),
          scratch.has_exact.data())) {
    return false;
  }
  std::size_t child_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    if (scratch.has_exact[lane_pos] != 0U) {
      if (current_time_by_lane[lane_pos] >= scratch.exact[lane_pos]) {
        compiled_math_batch_write_certain_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      } else {
        compiled_math_batch_write_impossible_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      }
      continue;
    }
    scratch.lanes_a[child_count++] = lane;
  }
  if (child_count == 0U) {
    return true;
  }
  return compiled_math_batch_source_product_program_fill(
      program,
      source_program.child_program_id,
      state,
      scratch.lanes_a.data(),
      child_count,
      current_time_by_lane,
      fill_mask,
      frame,
      scratch_depth + 1U,
      mask_by_lane,
      pdf_by_lane,
      cdf_by_lane,
      survival_by_lane);
}

inline bool compiled_math_batch_source_product_program_conditioned_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  const bool has_condition_overlay =
      compiled_math_source_product_bounds_have_overlay(
          source_program.bounds);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !has_condition_overlay) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      compiled_math_batch_write_forced_fill_lane(
          relation,
          lanes[i],
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
    }
    return true;
  }

  auto &scratch = frame->source_product_scratch_layer(
      scratch_depth, state.lane_count, lane_count);
  if (!compiled_math_batch_resolve_source_bounds(
          source_program,
          state,
          lanes,
          lane_count,
          scratch.lower.data(),
          scratch.upper.data(),
          scratch.exact.data(),
          scratch.has_exact.data())) {
    return false;
  }

  std::size_t unresolved_count = 0U;
  bool any_finite_upper = false;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower_bound = scratch.lower[lane_pos];
    const double upper_bound = scratch.upper[lane_pos];
    if (std::isfinite(upper_bound) && !(upper_bound > lower_bound)) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    if (scratch.has_exact[lane_pos] != 0U) {
      const double exact_time = scratch.exact[lane_pos];
      if ((lower_bound > 0.0 && !(exact_time > lower_bound)) ||
          (std::isfinite(upper_bound) && !(exact_time < upper_bound)) ||
          !(current_time_by_lane[lane_pos] >= exact_time)) {
        compiled_math_batch_write_impossible_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      } else {
        compiled_math_batch_write_certain_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      }
      continue;
    }
    if (std::isfinite(upper_bound)) {
      any_finite_upper = true;
    }
    scratch.lanes_a[unresolved_count++] = lane;
  }
  if (unresolved_count == 0U) {
    return true;
  }

  std::uint8_t uncond_mask = fill_mask;
  if (any_finite_upper && (fill_mask & kLeafChannelSurvival) != 0U) {
    uncond_mask |= kLeafChannelCdf;
  }
  if (!compiled_math_batch_source_product_program_fill(
          program,
          source_program.child_program_id,
          state,
          scratch.lanes_a.data(),
          unresolved_count,
          current_time_by_lane,
          uncond_mask,
          frame,
          scratch_depth + 1U,
          scratch.mask_a.data(),
          scratch.pdf_a.data(),
          scratch.cdf_a.data(),
          scratch.survival_a.data())) {
    return false;
  }

  std::size_t lower_cdf_count = 0U;
  std::size_t lower_survival_count = 0U;
  std::size_t upper_count = 0U;
  for (std::size_t i = 0; i < unresolved_count; ++i) {
    const auto lane = scratch.lanes_a[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower_bound = scratch.lower[lane_pos];
    const double upper_bound = scratch.upper[lane_pos];
    if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
      continue;
    }
    if (std::isfinite(upper_bound)) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          kLeafChannelCdf,
          scratch.mask_b.data(),
          scratch.pdf_b.data(),
          scratch.cdf_b.data(),
          scratch.survival_b.data());
      if (lower_bound > 0.0) {
        scratch.times[lane_pos] = lower_bound;
        scratch.lanes_b[lower_cdf_count++] = lane;
      }
      scratch.times[lane_pos] = upper_bound;
      scratch.lanes_d[upper_count++] = lane;
    } else {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          kLeafChannelSurvival,
          scratch.mask_b.data(),
          scratch.pdf_b.data(),
          scratch.cdf_b.data(),
          scratch.survival_b.data());
      if (lower_bound > 0.0) {
        scratch.times[lane_pos] = lower_bound;
        scratch.lanes_c[lower_survival_count++] = lane;
      }
    }
  }

  if (lower_cdf_count != 0U &&
      !compiled_math_batch_source_product_program_fill(
          program,
          source_program.child_program_id,
          state,
          scratch.lanes_b.data(),
          lower_cdf_count,
          scratch.times.data(),
          kLeafChannelCdf,
          frame,
          scratch_depth + 1U,
          scratch.mask_b.data(),
          scratch.pdf_b.data(),
          scratch.cdf_b.data(),
          scratch.survival_b.data())) {
    return false;
  }
  if (lower_survival_count != 0U &&
      !compiled_math_batch_source_product_program_fill(
          program,
          source_program.child_program_id,
          state,
          scratch.lanes_c.data(),
          lower_survival_count,
          scratch.times.data(),
          kLeafChannelSurvival,
          frame,
          scratch_depth + 1U,
          scratch.mask_b.data(),
          scratch.pdf_b.data(),
          scratch.cdf_b.data(),
          scratch.survival_b.data())) {
    return false;
  }
  if (upper_count != 0U &&
      !compiled_math_batch_source_product_program_fill(
          program,
          source_program.child_program_id,
          state,
          scratch.lanes_d.data(),
          upper_count,
          scratch.times.data(),
          kLeafChannelCdf,
          frame,
          scratch_depth + 1U,
          scratch.mask_c.data(),
          scratch.pdf_c.data(),
          scratch.cdf_c.data(),
          scratch.survival_c.data())) {
    return false;
  }

  for (std::size_t i = 0; i < unresolved_count; ++i) {
    const auto lane = scratch.lanes_a[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower_bound = scratch.lower[lane_pos];
    const double upper_bound = scratch.upper[lane_pos];
    if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
      compiled_math_batch_copy_fill_lane(
          lane,
          fill_mask,
          scratch.mask_a.data(),
          scratch.pdf_a.data(),
          scratch.cdf_a.data(),
          scratch.survival_a.data(),
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    if (!std::isfinite(upper_bound)) {
      const double surv_lb = scratch.survival_b[lane_pos];
      if (!std::isfinite(surv_lb) || !(surv_lb > 0.0)) {
        compiled_math_batch_write_impossible_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
        continue;
      }
      mask_by_lane[lane_pos] = fill_mask;
      if ((fill_mask & kLeafChannelPdf) != 0U) {
        pdf_by_lane[lane_pos] =
            safe_density(scratch.pdf_a[lane_pos] / surv_lb);
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        cdf_by_lane[lane_pos] =
            clamp_probability(
                (scratch.cdf_a[lane_pos] + surv_lb - 1.0) / surv_lb);
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        survival_by_lane[lane_pos] =
            clamp_probability(scratch.survival_a[lane_pos] / surv_lb);
      }
      continue;
    }
    const double mass =
        scratch.cdf_c[lane_pos] - scratch.cdf_b[lane_pos];
    if (!std::isfinite(mass) || !(mass > 0.0)) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    if (current_time_by_lane[lane_pos] >= upper_bound) {
      compiled_math_batch_write_certain_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    if (current_time_by_lane[lane_pos] <= lower_bound) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    mask_by_lane[lane_pos] = fill_mask;
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      pdf_by_lane[lane_pos] = safe_density(scratch.pdf_a[lane_pos] / mass);
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      cdf_by_lane[lane_pos] =
          clamp_probability(
              (scratch.cdf_a[lane_pos] - scratch.cdf_b[lane_pos]) / mass);
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      survival_by_lane[lane_pos] =
          clamp_probability(
              (scratch.cdf_c[lane_pos] - scratch.cdf_a[lane_pos]) / mass);
    }
  }
  return true;
}

inline bool compiled_math_batch_source_product_program_onset_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  if (source_program.onset_source_program_id == semantic::kInvalidIndex ||
      source_program.leaf_index == semantic::kInvalidIndex ||
      frame == nullptr || current_time_by_lane == nullptr) {
    return false;
  }
  const auto onset_pos =
      static_cast<std::size_t>(source_program.onset_source_program_id);
  if (onset_pos >= program.integral_kernel_source_product_programs.size()) {
    return false;
  }
  const auto &onset_program =
      program.integral_kernel_source_product_programs[onset_pos];
  auto &scratch = frame->source_product_scratch_layer(
      scratch_depth, state.lane_count, lane_count);

  std::size_t possible_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double upper =
        current_time_by_lane[lane_pos] - source_program.leaf_onset_lag;
    if (!(upper > 0.0)) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    scratch.onset_upper[lane_pos] = upper;
    scratch.lanes_a[possible_count++] = lane;
  }
  if (possible_count == 0U) {
    return true;
  }

  if (!compiled_math_batch_resolve_source_bounds_for_source(
          onset_program.source_id,
          source_program.onset_bounds,
          state,
          scratch.lanes_a.data(),
          possible_count,
          scratch.lower.data(),
          scratch.upper.data(),
          scratch.exact.data(),
          scratch.has_exact.data())) {
    return false;
  }

  std::size_t exact_count = 0U;
  std::size_t unresolved_count = 0U;
  for (std::size_t i = 0; i < possible_count; ++i) {
    const auto lane = scratch.lanes_a[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    if (scratch.has_exact[lane_pos] == 0U) {
      scratch.lanes_c[unresolved_count++] = lane;
      continue;
    }
    scratch.times[lane_pos] =
        current_time_by_lane[lane_pos] -
        scratch.exact[lane_pos] -
        source_program.leaf_onset_lag +
        source_program.leaf_onset_abs_value;
    scratch.lanes_b[exact_count++] = lane;
  }

  if (exact_count != 0U &&
      !compiled_math_batch_source_product_program_leaf_fill(
          source_program,
          state,
          scratch.lanes_b.data(),
          exact_count,
          scratch.times.data(),
          fill_mask,
          frame,
          scratch_depth + 1U,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane)) {
    return false;
  }
  if (unresolved_count == 0U) {
    return true;
  }

  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const std::uint8_t shifted_mask =
      static_cast<std::uint8_t>(
          (need_pdf ? kLeafChannelPdf : 0U) |
          (need_cdf ? kLeafChannelCdf : 0U));
  for (std::size_t i = 0; i < unresolved_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(scratch.lanes_c[i]);
    scratch.pdf_c[lane_pos] = 0.0;
    scratch.cdf_c[lane_pos] = 0.0;
  }

  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  for (std::size_t sample = 0U; sample < quadrature::kDefaultFiniteOrder;
       ++sample) {
    for (std::size_t i = 0U; i < unresolved_count; ++i) {
      const auto lane = scratch.lanes_c[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double scale = 0.5 * scratch.onset_upper[lane_pos];
      const double source_time =
          scale * (rule.nodes[sample] + 1.0);
      scratch.times[lane_pos] = source_time;
      scratch.lower[lane_pos] = scale * rule.weights[sample];
      scratch.upper[lane_pos] =
          current_time_by_lane[lane_pos] -
          source_time -
          source_program.leaf_onset_lag +
          source_program.leaf_onset_abs_value;
    }
    if (!compiled_math_batch_source_product_program_fill(
            program,
            source_program.onset_source_program_id,
            state,
            scratch.lanes_c.data(),
            unresolved_count,
            scratch.times.data(),
            kLeafChannelPdf,
            frame,
            scratch_depth + 1U,
            scratch.mask_a.data(),
            scratch.pdf_a.data(),
            scratch.cdf_a.data(),
            scratch.survival_a.data())) {
      return false;
    }
    if (!compiled_math_batch_source_product_program_leaf_fill(
            source_program,
            state,
            scratch.lanes_c.data(),
            unresolved_count,
            scratch.upper.data(),
            shifted_mask,
            frame,
            scratch_depth + 1U,
            scratch.mask_b.data(),
            scratch.pdf_b.data(),
            scratch.cdf_b.data(),
            scratch.survival_b.data())) {
      return false;
    }
    for (std::size_t i = 0U; i < unresolved_count; ++i) {
      const auto lane = scratch.lanes_c[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double source_pdf = scratch.pdf_a[lane_pos];
      if (!(source_pdf > 0.0)) {
        continue;
      }
      const double weight = scratch.lower[lane_pos] * source_pdf;
      if (need_pdf) {
        scratch.pdf_c[lane_pos] += weight * scratch.pdf_b[lane_pos];
      }
      if (need_cdf) {
        scratch.cdf_c[lane_pos] += weight * scratch.cdf_b[lane_pos];
      }
    }
  }

  for (std::size_t i = 0U; i < unresolved_count; ++i) {
    const auto lane = scratch.lanes_c[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    mask_by_lane[lane_pos] = fill_mask;
    if (need_pdf) {
      pdf_by_lane[lane_pos] = safe_density(scratch.pdf_c[lane_pos]);
    }
    if (need_cdf) {
      const double cdf = clamp_probability(scratch.cdf_c[lane_pos]);
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        cdf_by_lane[lane_pos] = cdf;
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        survival_by_lane[lane_pos] = clamp_probability(1.0 - cdf);
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_source_product_program_pool_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto member_count =
      static_cast<std::size_t>(source_program.member_programs.size);
  if (member_count == 0U || frame == nullptr ||
      current_time_by_lane == nullptr) {
    return false;
  }
  const auto k = static_cast<std::size_t>(source_program.pool_k);
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const std::uint8_t member_mask =
      static_cast<std::uint8_t>(
          kLeafChannelCdf | kLeafChannelSurvival |
          (need_pdf ? kLeafChannelPdf : 0U));

  auto &scratch = frame->source_product_scratch_layer(
      scratch_depth, state.lane_count, lane_count);
  scratch.ensure_pool_size(state.lane_count, member_count);
  const auto member_idx =
      [lane_stride = state.lane_count](const std::size_t member,
                                       const std::size_t lane_pos) {
        return member * lane_stride + lane_pos;
      };

  for (std::size_t member = 0U; member < member_count; ++member) {
    const auto member_program_id =
        program.integral_kernel_source_product_program_members[
            static_cast<std::size_t>(
                source_program.member_programs.offset) +
            member];
    if (member_program_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(member_program_id) >=
            program.integral_kernel_source_product_programs.size()) {
      return false;
    }
    if (!compiled_math_batch_source_product_program_fill(
            program,
            member_program_id,
            state,
            lanes,
            lane_count,
            current_time_by_lane,
            member_mask,
            frame,
            scratch_depth + 1U,
            scratch.mask_a.data(),
            scratch.pdf_a.data(),
            scratch.cdf_a.data(),
            scratch.survival_a.data())) {
      return false;
    }
    for (std::size_t i = 0U; i < lane_count; ++i) {
      const auto lane_pos = static_cast<std::size_t>(lanes[i]);
      const auto pos = member_idx(member, lane_pos);
      scratch.pool_member_pdf[pos] = scratch.pdf_a[lane_pos];
      scratch.pool_member_cdf[pos] = scratch.cdf_a[lane_pos];
      scratch.pool_member_survival[pos] = scratch.survival_a[lane_pos];
    }
  }

  const auto width = member_count + 1U;
  const auto table_stride = width * width;
  const auto table_idx =
      [table_stride, width](const std::size_t lane_pos,
                            const std::size_t row,
                            const std::size_t col) {
        return lane_pos * table_stride + row * width + col;
      };

  for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
    const auto lane = lanes[lane_i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    auto prefix = scratch.pool_prefix.data() + lane_pos * table_stride;
    auto suffix = scratch.pool_suffix.data() + lane_pos * table_stride;
    std::fill(prefix, prefix + table_stride, 0.0);
    std::fill(suffix, suffix + table_stride, 0.0);

    scratch.pool_prefix[table_idx(lane_pos, 0U, 0U)] = 1.0;
    for (std::size_t member = 0U; member < member_count; ++member) {
      const double member_cdf =
          scratch.pool_member_cdf[member_idx(member, lane_pos)];
      const double member_survival =
          scratch.pool_member_survival[member_idx(member, lane_pos)];
      for (std::size_t m = 0U; m <= member; ++m) {
        scratch.pool_prefix[table_idx(lane_pos, member + 1U, m)] +=
            scratch.pool_prefix[table_idx(lane_pos, member, m)] *
            member_survival;
        scratch.pool_prefix[table_idx(lane_pos, member + 1U, m + 1U)] +=
            scratch.pool_prefix[table_idx(lane_pos, member, m)] *
            member_cdf;
      }
    }

    scratch.pool_suffix[table_idx(lane_pos, member_count, 0U)] = 1.0;
    for (std::size_t rev = 0U; rev < member_count; ++rev) {
      const std::size_t member = member_count - 1U - rev;
      const std::size_t count = member_count - member - 1U;
      const double member_cdf =
          scratch.pool_member_cdf[member_idx(member, lane_pos)];
      const double member_survival =
          scratch.pool_member_survival[member_idx(member, lane_pos)];
      for (std::size_t m = 0U; m <= count; ++m) {
        scratch.pool_suffix[table_idx(lane_pos, member, m)] +=
            scratch.pool_suffix[table_idx(lane_pos, member + 1U, m)] *
            member_survival;
        scratch.pool_suffix[table_idx(lane_pos, member, m + 1U)] +=
            scratch.pool_suffix[table_idx(lane_pos, member + 1U, m)] *
            member_cdf;
      }
    }

    mask_by_lane[lane_pos] = fill_mask;
    if (need_pdf) {
      double pdf = 0.0;
      for (std::size_t member = 0U; member < member_count; ++member) {
        double others_exact = 0.0;
        for (std::size_t left = 0U; left < k; ++left) {
          if (left > member) {
            continue;
          }
          const auto right = k - 1U - left;
          if (right > member_count - member - 1U) {
            continue;
          }
          others_exact +=
              scratch.pool_prefix[table_idx(lane_pos, member, left)] *
              scratch.pool_suffix[table_idx(lane_pos, member + 1U, right)];
        }
        pdf += scratch.pool_member_pdf[member_idx(member, lane_pos)] *
               others_exact;
      }
      pdf_by_lane[lane_pos] = safe_density(pdf);
    }
    if (need_cdf) {
      double survival = 0.0;
      const auto survival_end = std::min(k, member_count + 1U);
      for (std::size_t m = 0U; m < survival_end; ++m) {
        survival +=
            scratch.pool_prefix[table_idx(lane_pos, member_count, m)];
      }
      const double clamped_survival = clamp_probability(survival);
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        survival_by_lane[lane_pos] = clamped_survival;
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        cdf_by_lane[lane_pos] = clamp_probability(1.0 - clamped_survival);
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_source_product_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  if (lane_count == 0U) {
    return true;
  }
  if (source_product_program_id == semantic::kInvalidIndex) {
    compiled_math_batch_write_fill_for_lanes(
        lanes,
        lane_count,
        fill_mask,
        false,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
    return true;
  }
  const auto pos = static_cast<std::size_t>(source_product_program_id);
  if (pos >= program.integral_kernel_source_product_programs.size() ||
      frame == nullptr || current_time_by_lane == nullptr ||
      mask_by_lane == nullptr || pdf_by_lane == nullptr ||
      cdf_by_lane == nullptr || survival_by_lane == nullptr) {
    return false;
  }
  const auto &source_program =
      program.integral_kernel_source_product_programs[pos];
  switch (source_program.kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
    compiled_math_batch_write_fill_for_lanes(
        lanes,
        lane_count,
        fill_mask,
        false,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
    return true;
  case CompiledMathSourceProductProgramKind::ConstantOne:
    compiled_math_batch_write_fill_for_lanes(
        lanes,
        lane_count,
        fill_mask,
        true,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
    return true;
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    return compiled_math_batch_source_product_program_leaf_fill(
        source_program,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        fill_mask,
        frame,
        scratch_depth,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
  case CompiledMathSourceProductProgramKind::ExactGate:
    return compiled_math_batch_source_product_program_exact_gate_fill(
        program,
        source_program,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        fill_mask,
        frame,
        scratch_depth,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
  case CompiledMathSourceProductProgramKind::Conditioned:
    return compiled_math_batch_source_product_program_conditioned_fill(
        program,
        source_program,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        fill_mask,
        frame,
        scratch_depth,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
    return compiled_math_batch_source_product_program_onset_fill(
        program,
        source_program,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        fill_mask,
        frame,
        scratch_depth,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    return compiled_math_batch_source_product_program_pool_fill(
        program,
        source_program,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        fill_mask,
        frame,
        scratch_depth,
        mask_by_lane,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane);
  }
  return false;
}

inline void compiled_math_batch_store_source_product_program_fill(
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t lane_stride,
    const semantic::Index cache_slot,
    const semantic::Index lane,
    const ExactSourceChannels::SourceProductScalarFill &fill) {
  const auto pos =
      static_cast<std::size_t>(cache_slot) * lane_stride +
      static_cast<std::size_t>(lane);
  frame->source_program_valid_mask[pos] |= fill.mask;
  if ((fill.mask & kLeafChannelPdf) != 0U) {
    frame->source_program_pdf[pos] = fill.pdf;
  }
  if ((fill.mask & kLeafChannelCdf) != 0U) {
    frame->source_program_cdf[pos] = fill.cdf;
  }
  if ((fill.mask & kLeafChannelSurvival) != 0U) {
    frame->source_program_survival[pos] = fill.survival;
  }
}

inline double compiled_math_batch_source_product_cached_value(
    const CompiledMathBatchIntegralFrame &frame,
    const std::size_t lane_stride,
    const semantic::Index cache_slot,
    const semantic::Index lane,
    const std::uint8_t channel_mask) {
  const auto pos =
      static_cast<std::size_t>(cache_slot) * lane_stride +
      static_cast<std::size_t>(lane);
  if (channel_mask == kLeafChannelPdf) {
    return frame.source_program_pdf[pos];
  }
  if (channel_mask == kLeafChannelCdf) {
    return frame.source_program_cdf[pos];
  }
  if (channel_mask == kLeafChannelSurvival) {
    return frame.source_program_survival[pos];
  }
  return 0.0;
}

inline double compiled_math_batch_source_product_array_value(
    const semantic::Index lane,
    const std::uint8_t channel_mask,
    const double *pdf_by_lane,
    const double *cdf_by_lane,
    const double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  if (channel_mask == kLeafChannelPdf) {
    return pdf_by_lane[lane_pos];
  }
  if (channel_mask == kLeafChannelCdf) {
    return cdf_by_lane[lane_pos];
  }
  if (channel_mask == kLeafChannelSurvival) {
    return survival_by_lane[lane_pos];
  }
  return 0.0;
}

template <typename LeafValueFn>
inline void compiled_math_batch_direct_leaf_values_for_kind(
    const CompiledMathSourceProductChannel &channel,
    const semantic::Index source_product_channel_id,
    const CompiledMathLaneBatchState &state,
    const BatchActiveLaneSpan active,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    LeafValueFn leaf_value_fn,
    bool *supported) {
  if (channel.leaf_index == semantic::kInvalidIndex) {
    *supported = false;
    return;
  }
  for (std::size_t i = 0; i < active.size; ++i) {
    const auto lane = active.lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr ||
        !source_channels->source_product_direct_leaf_available(
            source_product_channel_id)) {
      *supported = false;
      return;
    }
    const auto &loaded =
        source_channels->source_product_leaf_input(channel.leaf_index);
    const double x =
        batch_source_product_channel_time(channel, state.time_slots, lane) -
        channel.leaf_onset_abs_value - loaded.t0;
    values_by_lane[static_cast<std::size_t>(lane)] =
        leaf_value_fn(loaded, x, channel_mask);
  }
}

inline bool compiled_math_batch_direct_leaf_values(
    const CompiledMathSourceProductChannel &channel,
    const semantic::Index source_product_channel_id,
    const CompiledMathLaneBatchState &state,
    const BatchActiveLaneSpan active,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (channel.leaf_index == semantic::kInvalidIndex) {
    return false;
  }
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, active.size);
    for (std::size_t i = 0; i < active.size; ++i) {
      const auto lane = active.lanes[i];
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr ||
          !source_channels->source_product_direct_leaf_available(
              source_product_channel_id)) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      scratch->times[lane_pos] =
          batch_source_product_channel_time(channel, state.time_slots, lane);
    }
    return compiled_math_batch_leaf_values_from_times(
        channel.leaf_dist_kind,
        channel.leaf_index,
        channel.leaf_onset_abs_value,
        state,
        active.lanes,
        active.size,
        scratch->times.data(),
        channel_mask,
        values_by_lane,
        scratch);
  }
  bool supported = true;
  switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    compiled_math_batch_direct_leaf_values_for_kind(
        channel,
        source_product_channel_id,
        state,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_lognormal_leaf_value,
        &supported);
    break;
  case leaf::DistKind::Gamma:
    compiled_math_batch_direct_leaf_values_for_kind(
        channel,
        source_product_channel_id,
        state,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_gamma_leaf_value,
        &supported);
    break;
  case leaf::DistKind::Exgauss:
    compiled_math_batch_direct_leaf_values_for_kind(
        channel,
        source_product_channel_id,
        state,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_exgauss_leaf_value,
        &supported);
    break;
  case leaf::DistKind::LBA:
    compiled_math_batch_direct_leaf_values_for_kind(
        channel,
        source_product_channel_id,
        state,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_lba_leaf_value,
        &supported);
    break;
  case leaf::DistKind::RDM:
    compiled_math_batch_direct_leaf_values_for_kind(
        channel,
        source_product_channel_id,
        state,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_rdm_leaf_value,
        &supported);
    break;
  }
  return supported;
}

inline bool compiled_math_batch_source_product_value_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    const CompiledMathLaneBatchState &state,
    const BatchActiveLaneSpan active,
    CompiledMathBatchIntegralFrame *frame,
    double *out_by_lane) {
  if (frame == nullptr || out_by_lane == nullptr) {
    return false;
  }
  for (std::size_t i = 0; i < active.size; ++i) {
    const auto lane = active.lanes[i];
    out_by_lane[static_cast<std::size_t>(lane)] = 1.0;
    frame->source_lanes_a[i] = lane;
  }
  std::size_t current_count = active.size;
  auto *current_lanes = &frame->source_lanes_a;
  auto *next_lanes = &frame->source_lanes_b;

  for (semantic::Index op_idx = 0; op_idx < source_product_ops.size;
       ++op_idx) {
    if (current_count == 0U) {
      break;
    }
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + op_idx)];
    const auto channel_mask = op.value_channel_mask;
    if (channel_mask == 0U) {
      if (op.constant_value == 0.0) {
        for (std::size_t i = 0; i < current_count; ++i) {
          out_by_lane[static_cast<std::size_t>((*current_lanes)[i])] = 0.0;
        }
        current_count = 0U;
        break;
      }
      std::size_t next_count = 0U;
      for (std::size_t i = 0; i < current_count; ++i) {
        const auto lane = (*current_lanes)[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double product = out_by_lane[lane_pos] * op.constant_value;
        out_by_lane[lane_pos] =
            std::isfinite(product) && product != 0.0 ? product : 0.0;
        if (out_by_lane[lane_pos] != 0.0) {
          (*next_lanes)[next_count++] = lane;
        }
      }
      std::swap(current_lanes, next_lanes);
      current_count = next_count;
      continue;
    }

    if (op.source_product_channel_id == semantic::kInvalidIndex ||
        op.source_product_program_id == semantic::kInvalidIndex) {
      return false;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];

    bool direct_leaf_for_all_lanes =
        channel.leaf_index != semantic::kInvalidIndex;
    if (direct_leaf_for_all_lanes) {
      for (std::size_t i = 0; i < current_count; ++i) {
        auto *source_channels = state.source_channels((*current_lanes)[i]);
        if (source_channels == nullptr) {
          return false;
        }
        if (!source_channels->source_product_direct_leaf_available(
                op.source_product_channel_id)) {
          direct_leaf_for_all_lanes = false;
          break;
        }
      }
    }

    if (!op.cache_result && direct_leaf_for_all_lanes) {
      const BatchActiveLaneSpan current_active{
          current_lanes->data(),
          current_count};
      if (!compiled_math_batch_direct_leaf_values(
              channel,
              op.source_product_channel_id,
              state,
              current_active,
              channel_mask,
              frame->source_leaf_values.data(),
              &frame->source_product_scratch_layer(
                  0U,
                  state.lane_count,
                  current_count))) {
        return false;
      }
      std::size_t next_count = 0U;
      for (std::size_t i = 0; i < current_count; ++i) {
        const auto lane = (*current_lanes)[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double product =
            out_by_lane[lane_pos] * frame->source_leaf_values[lane_pos];
        out_by_lane[lane_pos] =
            std::isfinite(product) && product != 0.0 ? product : 0.0;
        if (out_by_lane[lane_pos] != 0.0) {
          (*next_lanes)[next_count++] = lane;
        }
      }
      std::swap(current_lanes, next_lanes);
      current_count = next_count;
      continue;
    }

    auto &fill_scratch =
        frame->source_product_scratch_layer(
            0U,
            state.lane_count,
            current_count);
    semantic::Index cache_slot = semantic::kInvalidIndex;
    if (!op.cache_result) {
      for (std::size_t i = 0; i < current_count; ++i) {
        const auto lane = (*current_lanes)[i];
        frame->bind_times[static_cast<std::size_t>(lane)] =
            batch_source_product_channel_time(
                channel,
                state.time_slots,
                lane);
      }
      if (!compiled_math_batch_source_product_program_fill(
              program,
              op.source_product_program_id,
              state,
              current_lanes->data(),
              current_count,
              frame->bind_times.data(),
              channel_mask,
              frame,
              1U,
              fill_scratch.mask_a.data(),
              fill_scratch.pdf_a.data(),
              fill_scratch.cdf_a.data(),
              fill_scratch.survival_a.data())) {
        return false;
      }
    } else {
      cache_slot =
          frame->source_program_cache_slot(op.source_product_program_id);
      if (cache_slot == semantic::kInvalidIndex) {
        return false;
      }
      std::size_t missing_count = 0U;
      for (std::size_t i = 0; i < current_count; ++i) {
        const auto lane = (*current_lanes)[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const auto program_pos =
            static_cast<std::size_t>(cache_slot) * state.lane_count +
            lane_pos;
        if ((frame->source_program_valid_mask[program_pos] &
             channel_mask) != 0U) {
          continue;
        }
        frame->bind_times[lane_pos] =
            batch_source_product_channel_time(
                channel,
                state.time_slots,
                lane);
        fill_scratch.lanes_a[missing_count++] = lane;
      }
      if (missing_count != 0U) {
        if (!compiled_math_batch_source_product_program_fill(
                program,
                op.source_product_program_id,
                state,
                fill_scratch.lanes_a.data(),
                missing_count,
                frame->bind_times.data(),
                op.fill_channel_mask,
                frame,
                1U,
                fill_scratch.mask_a.data(),
                fill_scratch.pdf_a.data(),
                fill_scratch.cdf_a.data(),
                fill_scratch.survival_a.data())) {
          return false;
        }
        for (std::size_t i = 0; i < missing_count; ++i) {
          const auto lane = fill_scratch.lanes_a[i];
          const auto lane_pos = static_cast<std::size_t>(lane);
          ExactSourceChannels::SourceProductScalarFill fill;
          fill.mask = fill_scratch.mask_a[lane_pos];
          fill.pdf = fill_scratch.pdf_a[lane_pos];
          fill.cdf = fill_scratch.cdf_a[lane_pos];
          fill.survival = fill_scratch.survival_a[lane_pos];
          compiled_math_batch_store_source_product_program_fill(
              frame,
              state.lane_count,
              cache_slot,
              lane,
              fill);
        }
      }
    }

    std::size_t next_count = 0U;
    for (std::size_t i = 0; i < current_count; ++i) {
      const auto lane = (*current_lanes)[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double value =
          op.cache_result
              ? compiled_math_batch_source_product_cached_value(
                    *frame,
                    state.lane_count,
                    cache_slot,
                    lane,
                    channel_mask)
              : compiled_math_batch_source_product_array_value(
                    lane,
                    channel_mask,
                    fill_scratch.pdf_a.data(),
                    fill_scratch.cdf_a.data(),
                    fill_scratch.survival_a.data());
      const double product = out_by_lane[lane_pos] * value;
      out_by_lane[lane_pos] =
          std::isfinite(product) && product != 0.0 ? product : 0.0;
      if (out_by_lane[lane_pos] != 0.0) {
        (*next_lanes)[next_count++] = lane;
      }
    }
    std::swap(current_lanes, next_lanes);
    current_count = next_count;
  }
  return true;
}

inline bool compiled_math_batch_outcome_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan outcome_gate_nodes,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane) {
  const auto *used_outcomes = state.used_outcomes(lane);
  if (used_outcomes == nullptr) {
    return true;
  }
  auto *parent = state.parent(lane);
  if (parent == nullptr) {
    return false;
  }
  for (semantic::Index i = 0; i < outcome_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_outcome_gate_nodes[
            static_cast<std::size_t>(outcome_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (gate_node.kind != CompiledMathNodeKind::OutcomeSubsetUnused) {
      return false;
    }
    const auto offset = static_cast<std::size_t>(gate_node.subject_id);
    const auto size = static_cast<std::size_t>(gate_node.aux_id);
    const auto &indices = parent->plan().compiled_outcome_gate_indices;
    if (offset + size > indices.size()) {
      return false;
    }
    for (std::size_t j = 0; j < size; ++j) {
      const auto outcome_idx = indices[offset + j];
      if (outcome_idx != semantic::kInvalidIndex &&
          static_cast<std::size_t>(outcome_idx) < used_outcomes->size() &&
          (*used_outcomes)[static_cast<std::size_t>(outcome_idx)] != 0U) {
        return false;
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_time_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan time_gate_nodes,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane) {
  for (semantic::Index i = 0; i < time_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_time_gate_nodes[
            static_cast<std::size_t>(time_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (gate_node.kind != CompiledMathNodeKind::TimeGate) {
      return false;
    }
    if (state.time(gate_node.time_id, lane) <
        state.time(gate_node.aux_id, lane)) {
      return false;
    }
  }
  return true;
}

inline double compiled_math_batch_timed_upper_bound_term_time(
    const CompiledMathTimedUpperBoundTerm &term,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane) {
  if (term.time_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(term.time_id) < state.time_slot_count &&
      state.time_slots.has(term.time_id, lane)) {
    return state.time(term.time_id, lane);
  }
  return state.time(node.time_id, lane);
}

inline bool compiled_math_batch_upper_bound_from_terms(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    CompiledTimedUpperBound *upper) {
  if (upper == nullptr) {
    return false;
  }
  if (node.aux_id == semantic::kInvalidIndex ||
      node.aux2_id == semantic::kInvalidIndex ||
      node.aux2_id == 0) {
    return true;
  }
  auto *workspace = state.workspace(lane);
  if (workspace == nullptr) {
    return false;
  }
  const auto offset = static_cast<std::size_t>(node.aux_id);
  const auto size = static_cast<std::size_t>(node.aux2_id);
  if (offset + size > program.timed_upper_bound_terms.size()) {
    return false;
  }
  for (std::size_t i = 0; i < size; ++i) {
    const auto &term = program.timed_upper_bound_terms[offset + i];
    if (term.normalizer_node_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(term.normalizer_node_id) >=
            workspace->values.size()) {
      return false;
    }
    const double fact_time =
        compiled_math_batch_timed_upper_bound_term_time(
            term,
            node,
            state,
            lane);
    const double normalizer =
        workspace->values[static_cast<std::size_t>(
            term.normalizer_node_id)];
    if (std::isfinite(fact_time) && normalizer > 0.0 &&
        (!upper->found || fact_time < upper->time)) {
      upper->found = true;
      upper->time = fact_time;
      upper->normalizer = normalizer;
    }
  }
  return true;
}

inline bool compiled_math_batch_expr_upper_bound_for_node(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const CompiledSourceView *condition_evaluator,
    CompiledTimedUpperBound *upper) {
  if (upper == nullptr) {
    return false;
  }
  *upper = CompiledTimedUpperBound{};
  if (!compiled_math_batch_upper_bound_from_terms(
          program,
          node,
          state,
          lane,
          upper)) {
    return false;
  }
  ExactTimedExprUpperBound sequence_upper;
  if (condition_evaluator != nullptr &&
      condition_evaluator->expr_upper_bound_for(
          node.subject_id,
          &sequence_upper) &&
      std::isfinite(sequence_upper.time) &&
      sequence_upper.normalizer > 0.0 &&
      (!upper->found || sequence_upper.time < upper->time)) {
    upper->found = true;
    upper->time = sequence_upper.time;
    upper->normalizer = sequence_upper.normalizer;
  }
  return true;
}

inline CompiledSourceView *compiled_math_batch_node_evaluator_for_lane(
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    bool *supported) {
  if (node.source_view_id == 0 ||
      node.source_view_id == semantic::kInvalidIndex) {
    return state.parent(lane);
  }
  auto *eval_workspace = state.eval_workspace(lane);
  if (eval_workspace == nullptr) {
    *supported = false;
    return nullptr;
  }
  return eval_workspace->source_view_evaluator(
      node.source_view_id,
      state.parent(lane));
}

inline bool compiled_math_batch_apply_expr_upper_factor_to_lanes(
    const CompiledMathProgram &program,
    const CompiledMathIntegralExprUpperFactor &factor,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *products_by_lane,
    semantic::Index *next_lanes,
    std::size_t *next_count) {
  if (products_by_lane == nullptr || next_lanes == nullptr ||
      next_count == nullptr) {
    return false;
  }
  const auto &node =
      program.nodes[static_cast<std::size_t>(factor.node_id)];
  if (node.kind != CompiledMathNodeKind::ExprUpperBoundDensity &&
      node.kind != CompiledMathNodeKind::ExprUpperBoundCdf) {
    return false;
  }
  std::size_t out_count = 0U;
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const auto lane = lanes[lane_pos];
    const auto lane_index = static_cast<std::size_t>(lane);
    bool supported = true;
    auto *condition_evaluator =
        compiled_math_batch_node_evaluator_for_lane(
            node,
            state,
            lane,
            &supported);
    if (!supported) {
      return false;
    }
    CompiledTimedUpperBound upper;
    if (!compiled_math_batch_expr_upper_bound_for_node(
            program,
            node,
            state,
            lane,
            condition_evaluator,
            &upper)) {
      return false;
    }
    const bool has_upper = upper.found && upper.normalizer > 0.0;
    const double current_time = state.time(node.time_id, lane);
    if (factor.mode == CompiledMathIntegralExprUpperMode::AfterOne) {
      if (!has_upper || current_time < upper.time) {
        products_by_lane[lane_index] = 0.0;
        continue;
      }
    } else if (has_upper) {
      if (current_time >= upper.time) {
        products_by_lane[lane_index] = 0.0;
        continue;
      }
      products_by_lane[lane_index] /= upper.normalizer;
    }
    if (products_by_lane[lane_index] != 0.0) {
      next_lanes[out_count++] = lane;
    }
  }
  *next_count = out_count;
  return true;
}

inline double compiled_math_batch_direct_leaf_channel_time(
    const CompiledMathSourceProductChannel &channel,
    const semantic::Index bind_time_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const double bind_time) noexcept {
  double time =
      channel.time_id == bind_time_id
          ? bind_time
          : state.time(channel.time_id, lane);
  if (channel.time_cap_id != semantic::kInvalidIndex) {
    const double cap =
        channel.time_cap_id == bind_time_id
            ? bind_time
            : state.time(channel.time_cap_id, lane);
    time = std::min(time, cap);
  }
  return time;
}

inline bool compiled_math_batch_lognormal_values_at_bind_times(
    const CompiledMathSourceProductChannel &channel,
    const semantic::Index source_product_channel_id,
    const semantic::Index bind_time_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *bind_times_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    double *scratch_x,
    double *scratch_y,
    semantic::Index *scratch_lanes) {
  if (channel.leaf_index == semantic::kInvalidIndex ||
      lanes == nullptr || bind_times_by_lane == nullptr ||
      values_by_lane == nullptr) {
    return false;
  }
  (void)source_product_channel_id;
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch_x != nullptr && scratch_y != nullptr &&
      scratch_lanes != nullptr) {
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const auto &loaded =
          source_channels->source_product_leaf_input(channel.leaf_index);
      const double time =
          compiled_math_batch_direct_leaf_channel_time(
              channel,
              bind_time_id,
              state,
              lane,
              bind_times_by_lane[lane_pos]);
      const double x = time - channel.leaf_onset_abs_value - loaded.t0;
      const double m = loaded.params[0];
      const double s = loaded.params[1];
      if (!(x > 0.0) || !std::isfinite(m) || !std::isfinite(s) ||
          s <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded.q,
                      channel_mask);
        continue;
      }
      scratch_lanes[valid_count] = lane;
      scratch_x[valid_count] = x;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch_y, scratch_x, &n);
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch_lanes[i];
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(channel.leaf_index);
        const double z = (scratch_y[i] - loaded.params[0]) / loaded.params[1];
        scratch_y[i] = -0.5 * z * z;
      }
      vvexp(scratch_y, scratch_y, &n);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch_lanes[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(channel.leaf_index);
        const double base_pdf =
            kInvSqrtTwoPi * scratch_y[i] /
            (scratch_x[i] * loaded.params[1]);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                0.0,
                loaded.q,
                channel_mask);
      }
    } else {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch_lanes[i];
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(channel.leaf_index);
        scratch_y[i] = (scratch_y[i] - loaded.params[0]) / loaded.params[1];
      }
      compiled_math_batch_normal_cdf_from_z(scratch_y, scratch_x, valid_count);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch_lanes[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        auto *source_channels = state.source_channels(lane);
        const auto &loaded =
            source_channels->source_product_leaf_input(channel.leaf_index);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                scratch_y[i],
                loaded.q,
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch_x;
  (void)scratch_y;
  (void)scratch_lanes;
#endif
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr) {
      return false;
    }
    const auto lane_pos = static_cast<std::size_t>(lane);
    const auto &loaded =
        source_channels->source_product_leaf_input(channel.leaf_index);
    const double time =
        compiled_math_batch_direct_leaf_channel_time(
            channel,
            bind_time_id,
            state,
            lane,
            bind_times_by_lane[lane_pos]);
    const double x = time - channel.leaf_onset_abs_value - loaded.t0;
    values_by_lane[lane_pos] =
        batch_source_product_lognormal_leaf_value(
            loaded,
            x,
            channel_mask);
  }
  return true;
}

template <typename LeafValueFn>
inline bool compiled_math_batch_direct_leaf_values_at_bind_times_for_kind(
    const CompiledMathSourceProductChannel &channel,
    const semantic::Index source_product_channel_id,
    const semantic::Index bind_time_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *bind_times_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    double *scratch_x,
    double *scratch_y,
    semantic::Index *scratch_lanes,
    LeafValueFn leaf_value_fn) {
  if (channel.leaf_index == semantic::kInvalidIndex ||
      lanes == nullptr || bind_times_by_lane == nullptr ||
      values_by_lane == nullptr) {
    return false;
  }
  (void)source_product_channel_id;
  (void)scratch_x;
  (void)scratch_y;
  (void)scratch_lanes;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr) {
      return false;
    }
    const auto lane_pos = static_cast<std::size_t>(lane);
    const auto &loaded =
        source_channels->source_product_leaf_input(channel.leaf_index);
    const double time =
        compiled_math_batch_direct_leaf_channel_time(
            channel,
            bind_time_id,
            state,
            lane,
            bind_times_by_lane[lane_pos]);
    const double x = time - channel.leaf_onset_abs_value - loaded.t0;
    values_by_lane[lane_pos] = leaf_value_fn(loaded, x, channel_mask);
  }
  return true;
}

inline bool compiled_math_batch_direct_leaf_values_at_bind_times(
    const CompiledMathSourceProductChannel &channel,
    const semantic::Index source_product_channel_id,
    const semantic::Index bind_time_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *bind_times_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    double *scratch_x,
    double *scratch_y,
    semantic::Index *scratch_lanes,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (channel.leaf_index == semantic::kInvalidIndex ||
      lanes == nullptr || bind_times_by_lane == nullptr ||
      values_by_lane == nullptr) {
    return false;
  }
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr) {
        return false;
      }
      if (!source_channels->source_product_direct_leaf_available(
              source_product_channel_id)) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      scratch->times[lane_pos] =
          compiled_math_batch_direct_leaf_channel_time(
              channel,
              bind_time_id,
              state,
              lane,
              bind_times_by_lane[lane_pos]);
    }
    return compiled_math_batch_leaf_values_from_times(
        channel.leaf_dist_kind,
        channel.leaf_index,
        channel.leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        scratch->times.data(),
        channel_mask,
        values_by_lane,
        scratch);
  }
  switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    return compiled_math_batch_lognormal_values_at_bind_times(
        channel,
        source_product_channel_id,
        bind_time_id,
        state,
        lanes,
        lane_count,
        bind_times_by_lane,
        channel_mask,
        values_by_lane,
        scratch_x,
        scratch_y,
        scratch_lanes);
  case leaf::DistKind::Gamma:
    return compiled_math_batch_direct_leaf_values_at_bind_times_for_kind(
        channel,
        source_product_channel_id,
        bind_time_id,
        state,
        lanes,
        lane_count,
        bind_times_by_lane,
        channel_mask,
        values_by_lane,
        scratch_x,
        scratch_y,
        scratch_lanes,
        batch_source_product_gamma_leaf_value);
  case leaf::DistKind::Exgauss:
    return compiled_math_batch_direct_leaf_values_at_bind_times_for_kind(
        channel,
        source_product_channel_id,
        bind_time_id,
        state,
        lanes,
        lane_count,
        bind_times_by_lane,
        channel_mask,
        values_by_lane,
        scratch_x,
        scratch_y,
        scratch_lanes,
        batch_source_product_exgauss_leaf_value);
  case leaf::DistKind::LBA:
    return compiled_math_batch_direct_leaf_values_at_bind_times_for_kind(
        channel,
        source_product_channel_id,
        bind_time_id,
        state,
        lanes,
        lane_count,
        bind_times_by_lane,
        channel_mask,
        values_by_lane,
        scratch_x,
        scratch_y,
        scratch_lanes,
        batch_source_product_lba_leaf_value);
  case leaf::DistKind::RDM:
    return compiled_math_batch_direct_leaf_values_at_bind_times_for_kind(
        channel,
        source_product_channel_id,
        bind_time_id,
        state,
        lanes,
        lane_count,
        bind_times_by_lane,
        channel_mask,
        values_by_lane,
        scratch_x,
        scratch_y,
        scratch_lanes,
        batch_source_product_rdm_leaf_value);
  }
  return false;
}

inline bool evaluate_compiled_integral_kernel_lane_batch_direct_leaf(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const BatchActiveLaneSpan parent_active,
    const double *lower_by_lane,
    const double *upper_by_lane,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane) {
  if (batch_workspace == nullptr || lower_by_lane == nullptr ||
      upper_by_lane == nullptr || out_by_parent_lane == nullptr ||
      !batch_finite_integral_kernel_direct_leaf_supported(program, kernel)) {
    return false;
  }
  const auto max_leaf =
      batch_source_product_max_leaf_index(program, kernel.source_product_ops);
  if (max_leaf == semantic::kInvalidIndex) {
    return false;
  }
  (void)max_leaf;
  for (semantic::Index op_idx = 0; op_idx < kernel.source_product_ops.size;
       ++op_idx) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(
                kernel.source_product_ops.offset + op_idx)];
    if (op.value_channel_mask == 0U) {
      continue;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];
    if (channel.time_id != kernel.bind_time_id &&
        static_cast<std::size_t>(channel.time_id) >=
            parent_state.time_slot_count) {
      return false;
    }
    if (channel.time_cap_id != semantic::kInvalidIndex &&
        channel.time_cap_id != kernel.bind_time_id &&
        static_cast<std::size_t>(channel.time_cap_id) >=
            parent_state.time_slot_count) {
      return false;
    }
  }
  for (std::size_t i = 0; i < parent_active.size; ++i) {
    auto *source_channels = parent_state.source_channels(parent_active.lanes[i]);
    if (!batch_source_product_ops_have_available_direct_leaf_inputs(
            program,
            kernel.source_product_ops,
            source_channels)) {
      return false;
    }
  }

  auto &frame = compiled_math_batch_integral_frame(batch_workspace, depth);
  frame.ensure_lane_capacity(
      parent_state.time_slot_count,
      parent_state.lane_count);
  std::size_t eligible_count = 0U;
  for (std::size_t active_pos = 0; active_pos < parent_active.size;
       ++active_pos) {
    const auto lane = parent_active.lanes[active_pos];
    const auto lane_pos = static_cast<std::size_t>(lane);
    out_by_parent_lane[lane_pos] = 0.0;
    const double lower = lower_by_lane[lane_pos];
    const double upper = upper_by_lane[lane_pos];
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      continue;
    }
    frame.lower_by_lane[lane_pos] = 0.5 * (upper - lower);
    frame.upper_by_lane[lane_pos] = 0.5 * (upper + lower);
    frame.active_lanes[eligible_count++] = lane;
  }
  if (eligible_count == 0U) {
    return true;
  }

  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  for (std::size_t q = 0; q < quadrature::kDefaultFiniteOrder; ++q) {
    std::size_t current_count = eligible_count;
    for (std::size_t i = 0; i < eligible_count; ++i) {
      const auto lane = frame.active_lanes[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double scale = frame.lower_by_lane[lane_pos];
      frame.source_lanes_a[i] = lane;
      frame.products[lane_pos] = 1.0;
      frame.weights[lane_pos] = scale * rule.weights[q];
      frame.bind_times[lane_pos] =
          frame.upper_by_lane[lane_pos] + scale * rule.nodes[q];
    }
    auto *current_lanes = &frame.source_lanes_a;
    auto *next_lanes = &frame.source_lanes_b;

    for (semantic::Index op_idx = 0; op_idx < kernel.source_product_ops.size;
         ++op_idx) {
      if (current_count == 0U) {
        break;
      }
      const auto &op =
          program.integral_kernel_source_product_ops[
              static_cast<std::size_t>(
                  kernel.source_product_ops.offset + op_idx)];
      if (op.value_channel_mask == 0U) {
        if (op.constant_value == 0.0) {
          current_count = 0U;
          break;
        }
        std::size_t next_count = 0U;
        for (std::size_t i = 0; i < current_count; ++i) {
          const auto lane = (*current_lanes)[i];
          const auto lane_pos = static_cast<std::size_t>(lane);
          const double product =
              frame.products[lane_pos] * op.constant_value;
          frame.products[lane_pos] =
              std::isfinite(product) && product != 0.0 ? product : 0.0;
          if (frame.products[lane_pos] != 0.0) {
            (*next_lanes)[next_count++] = lane;
          }
        }
        std::swap(current_lanes, next_lanes);
        current_count = next_count;
        continue;
      }

      const auto &channel =
          program.integral_kernel_source_product_channels[
              static_cast<std::size_t>(op.source_product_channel_id)];
      if (!compiled_math_batch_direct_leaf_values_at_bind_times(
              channel,
              op.source_product_channel_id,
              kernel.bind_time_id,
              parent_state,
              current_lanes->data(),
              current_count,
              frame.bind_times.data(),
              op.value_channel_mask,
              frame.source_leaf_values.data(),
              frame.source_values.data(),
              frame.factor_values.data(),
              frame.term_lanes_a.data(),
              &frame.source_product_scratch_layer(
                  0U,
                  parent_state.lane_count,
                  current_count))) {
        return false;
      }
      std::size_t next_count = 0U;
      for (std::size_t i = 0; i < current_count; ++i) {
        const auto lane = (*current_lanes)[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double product =
            frame.products[lane_pos] *
            frame.source_leaf_values[lane_pos];
        frame.products[lane_pos] =
            std::isfinite(product) && product != 0.0 ? product : 0.0;
        if (frame.products[lane_pos] != 0.0) {
          (*next_lanes)[next_count++] = lane;
        }
      }
      std::swap(current_lanes, next_lanes);
      current_count = next_count;
    }

    for (std::size_t i = 0; i < current_count; ++i) {
      const auto lane = (*current_lanes)[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      out_by_parent_lane[lane_pos] +=
          frame.weights[lane_pos] * frame.products[lane_pos];
    }
  }
  return true;
}

inline bool evaluate_compiled_integral_kernel_lane_batch(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const BatchActiveLaneSpan parent_active,
    const double *lower_by_lane,
    const double *upper_by_lane,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane) {
  if (batch_workspace == nullptr || lower_by_lane == nullptr ||
      upper_by_lane == nullptr || out_by_parent_lane == nullptr ||
      !compiled_math_batch_kernel_supported(program, kernel, 0U)) {
    return false;
  }
  for (std::size_t i = 0; i < parent_active.size; ++i) {
    out_by_parent_lane[static_cast<std::size_t>(parent_active.lanes[i])] =
        0.0;
  }
  if (parent_active.size == 0U) {
    return true;
  }
  if (evaluate_compiled_integral_kernel_lane_batch_direct_leaf(
          program,
          kernel,
          parent_state,
          parent_active,
          lower_by_lane,
          upper_by_lane,
          batch_workspace,
          depth,
          out_by_parent_lane)) {
    return true;
  }

  const auto time_slot_count =
      std::max(
          parent_state.time_slot_count,
          static_cast<std::size_t>(kernel.bind_time_id) + 1U);
  auto &frame = compiled_math_batch_integral_frame(batch_workspace, depth);
  const auto child_capacity =
      parent_active.size * quadrature::kDefaultFiniteOrder;
  frame.ensure_lane_capacity(time_slot_count, child_capacity);
  if (compiled_math_kernel_needs_batch_source_cache(program, kernel)) {
    compiled_math_collect_kernel_cached_source_programs(
        program,
        kernel,
        &frame.cached_source_programs);
    frame.reset_source_product_cache(
        program.integral_kernel_source_product_programs.size(),
        child_capacity,
        frame.cached_source_programs);
  }

  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  std::size_t child_count = 0U;
  for (std::size_t parent_pos = 0; parent_pos < parent_active.size;
       ++parent_pos) {
    const auto parent_lane = parent_active.lanes[parent_pos];
    const auto parent_lane_pos = static_cast<std::size_t>(parent_lane);
    const double lower = lower_by_lane[parent_lane_pos];
    const double upper = upper_by_lane[parent_lane_pos];
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      continue;
    }
    const double scale = 0.5 * (upper - lower);
    const double shift = 0.5 * (upper + lower);
    for (std::size_t q = 0; q < quadrature::kDefaultFiniteOrder; ++q) {
      const auto child_lane = static_cast<semantic::Index>(child_count);
      frame.active_lanes[child_count] = child_lane;
      frame.parent_lanes[child_count] = parent_lane;
      frame.weights[child_count] = scale * rule.weights[q];
      frame.workspaces_by_lane[child_count] =
          parent_state.workspace(parent_lane);
      frame.source_channels_by_lane[child_count] =
          parent_state.source_channels(parent_lane);
      frame.parents_by_lane[child_count] = parent_state.parent(parent_lane);
      frame.eval_workspaces_by_lane[child_count] =
          parent_state.eval_workspace(parent_lane);
      frame.used_outcomes_by_lane[child_count] =
          parent_state.used_outcomes(parent_lane);
      for (std::size_t time_slot = 0; time_slot < time_slot_count;
           ++time_slot) {
        double value = 0.0;
        if (time_slot < parent_state.time_slot_count) {
          value =
              parent_state.time_slots.get(
                  static_cast<semantic::Index>(time_slot),
                  parent_lane);
        }
        frame.time_values[time_slot * child_capacity + child_count] = value;
        frame.time_valid[time_slot * child_capacity + child_count] =
            time_slot < parent_state.time_slot_count &&
                    parent_state.time_slots.has(
                        static_cast<semantic::Index>(time_slot),
                        parent_lane)
                ? 1U
                : 0U;
      }
      frame.time_values[
          static_cast<std::size_t>(kernel.bind_time_id) * child_capacity +
          child_count] = shift + scale * rule.nodes[q];
      frame.time_valid[
          static_cast<std::size_t>(kernel.bind_time_id) * child_capacity +
          child_count] = 1U;
      ++child_count;
    }
  }
  if (child_count == 0U) {
    return true;
  }

  const CompiledMathLaneBatchState child_state{
      child_capacity,
      time_slot_count,
      BatchTimeSlotView{
          frame.time_values.data(),
          child_capacity,
          frame.time_valid.data()},
      &frame.workspaces_by_lane,
      &frame.source_channels_by_lane,
      &frame.parents_by_lane,
      &frame.eval_workspaces_by_lane,
      &frame.used_outcomes_by_lane};
  const BatchActiveLaneSpan child_active{frame.active_lanes.data(), child_count};

  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
    if (!compiled_math_batch_source_product_value_for_ops(
            program,
            kernel.source_product_ops,
            child_state,
            child_active,
            &frame,
            frame.values.data())) {
      return false;
    }
  } else if (kernel.kind == CompiledMathIntegralKernelKind::SourceProductSum) {
    for (std::size_t i = 0; i < child_count; ++i) {
      frame.values[i] = 0.0;
    }
    for (semantic::Index term_idx = 0;
         term_idx < kernel.source_product_terms.size;
         ++term_idx) {
      const auto &term =
          program.integral_kernel_source_product_terms[
              static_cast<std::size_t>(
                  kernel.source_product_terms.offset + term_idx)];
      std::size_t term_count = 0U;
      for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
        const auto lane = frame.active_lanes[child_pos];
        if (!compiled_math_batch_outcome_gates_open(
                program,
                term.outcome_gate_nodes,
                child_state,
                lane) ||
            !compiled_math_batch_time_gates_open(
                program,
                term.time_gate_nodes,
                child_state,
                lane)) {
          continue;
        }
        frame.term_lanes_a[term_count++] = lane;
        frame.products[static_cast<std::size_t>(lane)] = term.sign;
      }
      auto *current_lanes = &frame.term_lanes_a;
      auto *next_lanes = &frame.term_lanes_b;

      for (semantic::Index i = 0; i < term.integral_factor_nodes.size;
           ++i) {
        if (term_count == 0U) {
          break;
        }
        const auto node_id =
            program.integral_kernel_integral_factor_nodes[
                static_cast<std::size_t>(
                    term.integral_factor_nodes.offset + i)];
        const auto &factor_node =
            program.nodes[static_cast<std::size_t>(node_id)];
        if (factor_node.integral_kernel_slot == semantic::kInvalidIndex) {
          return false;
        }
        const auto &factor_kernel =
            program.integral_kernels[
                static_cast<std::size_t>(
                    factor_node.integral_kernel_slot)];
        for (std::size_t lane_pos = 0; lane_pos < term_count; ++lane_pos) {
          const auto lane = (*current_lanes)[lane_pos];
          const auto lane_index = static_cast<std::size_t>(lane);
          const double current_time = child_state.time(
              factor_node.time_id,
              lane);
          frame.lower_by_lane[lane_index] = 0.0;
          frame.upper_by_lane[lane_index] = current_time;
        }
        const BatchActiveLaneSpan factor_active{
            current_lanes->data(),
            term_count};
        if (!evaluate_compiled_integral_kernel_lane_batch(
                program,
                factor_node,
                factor_kernel,
                child_state,
                factor_active,
                frame.lower_by_lane.data(),
                frame.upper_by_lane.data(),
                batch_workspace,
                depth + 1U,
                frame.factor_values.data())) {
          return false;
        }
        std::size_t next_count = 0U;
        for (std::size_t lane_pos = 0; lane_pos < term_count; ++lane_pos) {
          const auto lane = (*current_lanes)[lane_pos];
          const auto lane_index = static_cast<std::size_t>(lane);
          double factor_value = frame.factor_values[lane_index];
          if (factor_node.kind == CompiledMathNodeKind::IntegralZeroToCurrent) {
            factor_value = clamp_probability(factor_value);
          }
          const double product =
              frame.products[lane_index] * factor_value;
          frame.products[lane_index] =
              std::isfinite(product) && product != 0.0 ? product : 0.0;
          if (frame.products[lane_index] != 0.0) {
            (*next_lanes)[next_count++] = lane;
          }
        }
        std::swap(current_lanes, next_lanes);
        term_count = next_count;
      }

      for (semantic::Index i = 0; i < term.expr_upper_factors.size; ++i) {
        if (term_count == 0U) {
          break;
        }
        const auto &factor =
            program.integral_kernel_expr_upper_factors[
                static_cast<std::size_t>(
                    term.expr_upper_factors.offset + i)];
        std::size_t next_count = 0U;
        if (!compiled_math_batch_apply_expr_upper_factor_to_lanes(
                program,
                factor,
                child_state,
                current_lanes->data(),
                term_count,
                frame.products.data(),
                next_lanes->data(),
                &next_count)) {
          return false;
        }
        std::swap(current_lanes, next_lanes);
        term_count = next_count;
      }

      if (term_count == 0U) {
        continue;
      }
      const BatchActiveLaneSpan source_active{
          current_lanes->data(),
          term_count};
      if (!compiled_math_batch_source_product_value_for_ops(
              program,
              term.source_product_ops,
              child_state,
              source_active,
              &frame,
              frame.source_values.data())) {
        return false;
      }
      for (std::size_t lane_pos = 0; lane_pos < term_count; ++lane_pos) {
        const auto lane = (*current_lanes)[lane_pos];
        const auto lane_index = static_cast<std::size_t>(lane);
        const double product =
            frame.products[lane_index] * frame.source_values[lane_index];
        if (std::isfinite(product) && product != 0.0) {
          frame.values[lane_index] += product;
        }
      }
    }
    if (kernel.clean_signed_source_sum) {
      for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
        const auto lane = frame.active_lanes[child_pos];
        const auto lane_index = static_cast<std::size_t>(lane);
        frame.values[lane_index] = clean_signed_value(frame.values[lane_index]);
      }
    }
  } else {
    return false;
  }

  for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
    const auto child_lane = frame.active_lanes[child_pos];
    const double value = frame.values[static_cast<std::size_t>(child_lane)];
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    out_by_parent_lane[
        static_cast<std::size_t>(frame.parent_lanes[child_pos])] +=
        frame.weights[child_pos] * value;
  }
  (void)node;
  return true;
}

inline bool evaluate_compiled_integral_node_batch_lane_native(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathIntegralKernel &kernel,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *node_values) {
  if (batch_workspace == nullptr || node_values == nullptr ||
      !compiled_math_batch_kernel_supported(program, kernel, 0U)) {
    return false;
  }

  const std::size_t time_slot_count =
      compiled_math_program_time_slot_count(program);
  auto &handled_lane_indices = batch_workspace->active_lane_indices;
  auto &compact_lane_indices = batch_workspace->compact_lane_indices;
  auto &active_lanes = batch_workspace->active_lanes;
  auto &lower_by_lane = batch_workspace->lower_by_lane;
  auto &upper_by_lane = batch_workspace->upper_by_lane;
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &source_channels_by_lane = batch_workspace->source_channels_by_lane;
  auto &parents_by_lane = batch_workspace->parents_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;
  auto &batch_values = batch_workspace->batch_values;

  handled_lane_indices.clear();
  compact_lane_indices.clear();
  active_lanes.clear();
  lower_by_lane.clear();
  upper_by_lane.clear();

  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    auto *workspace = lanes[lane].workspace;
    auto *parent = lanes[lane].parent;
    if (workspace == nullptr || lane >= condition_evaluators.size() ||
        condition_evaluators[lane] == nullptr || parent == nullptr ||
        parent->source_channels() == nullptr) {
      handled_lane_indices.push_back(lane);
      (*node_values)[lane] = 0.0;
      continue;
    }
    const auto compact_lane =
        static_cast<semantic::Index>(compact_lane_indices.size());
    handled_lane_indices.push_back(lane);
    compact_lane_indices.push_back(lane);
    lower_by_lane.push_back(0.0);
    const double upper = compiled_math_node_time(node, *workspace);
    upper_by_lane.push_back(upper);
    if (std::isfinite(upper) && upper > 0.0) {
      active_lanes.push_back(compact_lane);
    }
  }

  if (handled_lane_indices.empty()) {
    return false;
  }
  if (compact_lane_indices.empty()) {
    return true;
  }
  const std::size_t compact_count = compact_lane_indices.size();
  parent_time_values.assign(time_slot_count * compact_count, 0.0);
  parent_time_valid.assign(time_slot_count * compact_count, 0U);
  workspaces_by_lane.assign(compact_count, nullptr);
  source_channels_by_lane.assign(compact_count, nullptr);
  parents_by_lane.assign(compact_count, nullptr);
  eval_workspaces_by_lane.assign(compact_count, nullptr);
  used_outcomes_by_lane.assign(compact_count, nullptr);
  batch_values.assign(compact_count, 0.0);

  for (std::size_t compact_lane = 0; compact_lane < compact_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    auto *workspace = lanes[lane].workspace;
    for (std::size_t time_slot = 0; time_slot < time_slot_count;
         ++time_slot) {
      const auto time_pos = time_slot * compact_count + compact_lane;
      if (time_slot < workspace->time_values.size()) {
        parent_time_values[time_pos] = workspace->time_values[time_slot];
      }
      parent_time_valid[time_pos] =
          workspace->has_time(static_cast<semantic::Index>(time_slot))
              ? 1U
              : 0U;
    }
    workspaces_by_lane[compact_lane] = workspace;
    source_channels_by_lane[compact_lane] =
        lanes[lane].parent->source_channels();
    parents_by_lane[compact_lane] = lanes[lane].parent;
    eval_workspaces_by_lane[compact_lane] = lanes[lane].eval_workspace;
    used_outcomes_by_lane[compact_lane] = workspace->used_outcomes;
  }

  const CompiledMathLaneBatchState parent_state{
      compact_count,
      time_slot_count,
      BatchTimeSlotView{
          parent_time_values.data(),
          compact_count,
          parent_time_valid.data()},
      &workspaces_by_lane,
      &source_channels_by_lane,
      &parents_by_lane,
      &eval_workspaces_by_lane,
      &used_outcomes_by_lane};
  const BatchActiveLaneSpan active{
      active_lanes.data(),
      active_lanes.size()};
  if (!evaluate_compiled_integral_kernel_lane_batch(
          program,
          node,
          kernel,
          parent_state,
          active,
          lower_by_lane.data(),
          upper_by_lane.data(),
          batch_workspace,
          0U,
          batch_values.data())) {
    return false;
  }

  for (std::size_t compact_lane = 0; compact_lane < compact_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    double value = batch_values[compact_lane];
    if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent) {
      value = clamp_probability(value);
    } else if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw) {
      value = std::isfinite(value) ? clean_signed_value(value) : 0.0;
    } else {
      value = std::isfinite(value) ? value : 0.0;
    }
    (*node_values)[lane] = value;
  }
  return true;
}

inline bool evaluate_compiled_integral_node_batch_direct_leaf(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathIntegralKernel &kernel,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *node_values) {
  if (batch_workspace == nullptr || node_values == nullptr ||
      !batch_finite_integral_kernel_direct_leaf_supported(program, kernel)) {
    return false;
  }
  const auto max_leaf =
      batch_source_product_max_leaf_index(program, kernel.source_product_ops);
  const std::size_t leaf_count =
      max_leaf == semantic::kInvalidIndex
          ? 0U
          : static_cast<std::size_t>(max_leaf) + 1U;
  const auto max_time =
      batch_source_product_max_time_id_for_ops(program, kernel.source_product_ops);
  const std::size_t time_slot_count =
      static_cast<std::size_t>(
          std::max(max_time, kernel.bind_time_id)) +
      1U;

  auto &active_lanes = batch_workspace->active_lanes;
  auto &active_lane_indices = batch_workspace->active_lane_indices;
  auto &lower_by_lane = batch_workspace->lower_by_lane;
  auto &upper_by_lane = batch_workspace->upper_by_lane;
  active_lanes.clear();
  active_lane_indices.clear();
  lower_by_lane.clear();
  upper_by_lane.clear();

  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    auto *workspace = lanes[lane].workspace;
    auto *parent = lanes[lane].parent;
    if (workspace == nullptr || parent == nullptr ||
        lane >= condition_evaluators.size() ||
        condition_evaluators[lane] == nullptr ||
        workspace->time_values.size() < time_slot_count) {
      continue;
    }
    auto *source_channels = parent->source_channels();
    if (source_channels == nullptr ||
        !batch_source_product_ops_have_available_direct_leaf_inputs(
            program,
            kernel.source_product_ops,
            source_channels)) {
      continue;
    }
    const double upper = compiled_math_node_time(node, *workspace);
    if (!(upper > 0.0) || !std::isfinite(upper)) {
      continue;
    }
    active_lane_indices.push_back(lane);
    active_lanes.push_back(
        static_cast<semantic::Index>(active_lanes.size()));
    lower_by_lane.push_back(0.0);
    upper_by_lane.push_back(upper);
  }

  if (active_lanes.empty()) {
    return false;
  }

  const std::size_t batch_lane_count = active_lanes.size();
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_leaf_inputs = batch_workspace->parent_leaf_inputs;
  auto &batch_values = batch_workspace->batch_values;
  parent_time_values.assign(time_slot_count * batch_lane_count, 0.0);
  parent_leaf_inputs.resize(leaf_count * batch_lane_count);
  batch_values.assign(batch_lane_count, 0.0);

  for (std::size_t batch_lane = 0; batch_lane < batch_lane_count;
       ++batch_lane) {
    const auto lane = active_lane_indices[batch_lane];
    auto *workspace = lanes[lane].workspace;
    auto *source_channels = lanes[lane].parent->source_channels();
    for (std::size_t time_slot = 0; time_slot < time_slot_count;
         ++time_slot) {
      parent_time_values[time_slot * batch_lane_count + batch_lane] =
          workspace->time_values[time_slot];
    }
    for (std::size_t leaf = 0; leaf < leaf_count; ++leaf) {
      parent_leaf_inputs[leaf * batch_lane_count + batch_lane] =
          source_channels->source_product_leaf_input(
              static_cast<semantic::Index>(leaf));
    }
  }

  const BatchSourceProductContext parent_context{
      batch_lane_count,
      BatchTimeSlotView{parent_time_values.data(), batch_lane_count},
      BatchLeafInputView{parent_leaf_inputs.data(), batch_lane_count},
      0};
  const BatchActiveLaneSpan active{
      active_lanes.data(),
      batch_lane_count};
  if (!batch_finite_integral_source_product_direct_leaf(
          program,
          kernel,
          parent_context,
          active,
          lower_by_lane.data(),
          upper_by_lane.data(),
          time_slot_count,
          leaf_count,
          &batch_workspace->finite_integral_workspace,
          batch_values.data())) {
    return false;
  }

  for (std::size_t batch_lane = 0; batch_lane < batch_lane_count;
       ++batch_lane) {
    const auto lane = active_lane_indices[batch_lane];
    double value = batch_values[batch_lane];
    if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent) {
      value = clamp_probability(value);
    } else if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw) {
      value = std::isfinite(value) ? clean_signed_value(value) : 0.0;
    } else {
      value = std::isfinite(value) ? value : 0.0;
    }
    (*node_values)[lane] = value;
  }
  return true;
}

inline void evaluate_compiled_integral_node_batch(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *node_values) {
  node_values->assign(lanes.size(), 0.0);
  if (node.integral_kernel_slot == semantic::kInvalidIndex) {
    return;
  }
  const auto kernel_pos =
      static_cast<std::size_t>(node.integral_kernel_slot);
  if (kernel_pos >= program.integral_kernels.size()) {
    throw std::runtime_error(
        "compiled integral node points outside integral kernels");
  }
  const auto &kernel = program.integral_kernels[kernel_pos];
  std::vector<std::uint8_t> batched(lanes.size(), 0U);
  if (evaluate_compiled_integral_node_batch_lane_native(
          program,
          node,
          kernel,
          lanes,
          condition_evaluators,
          batch_workspace,
          node_values)) {
    for (const auto lane : batch_workspace->active_lane_indices) {
      batched[lane] = 1U;
    }
  } else if (evaluate_compiled_integral_node_batch_direct_leaf(
          program,
          node,
          kernel,
          lanes,
          condition_evaluators,
          batch_workspace,
          node_values)) {
    for (const auto lane : batch_workspace->active_lane_indices) {
      batched[lane] = 1U;
    }
  }

  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    if (batched[lane] != 0U) {
      continue;
    }
    auto *workspace = lanes[lane].workspace;
    if (workspace == nullptr || lane >= condition_evaluators.size() ||
        condition_evaluators[lane] == nullptr ||
        node.integral_kernel_slot == semantic::kInvalidIndex) {
      (*node_values)[lane] = 0.0;
      continue;
    }
    const double current_time = compiled_math_node_time(node, *workspace);
    if (!(current_time > 0.0)) {
      (*node_values)[lane] = 0.0;
      continue;
    }
    throw std::runtime_error(
        "compiled batch integral node has no batch implementation");
  }
}

inline bool prepare_compiled_source_node_batch_base(
    const CompiledMathProgram &program,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  if (batch_workspace == nullptr) {
    return false;
  }
  auto &compact_lane_indices = batch_workspace->compact_lane_indices;
  auto &active_lanes = batch_workspace->active_lanes;
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;

  compact_lane_indices.clear();
  active_lanes.clear();
  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    if (lanes[lane].workspace == nullptr) {
      return false;
    }
    compact_lane_indices.push_back(lane);
    active_lanes.push_back(static_cast<semantic::Index>(lane));
  }
  const std::size_t lane_count = compact_lane_indices.size();
  if (lane_count == 0U) {
    return true;
  }

  parent_time_values.resize(time_slot_count * lane_count);
  parent_time_valid.resize(time_slot_count * lane_count);
  workspaces_by_lane.resize(lane_count);
  eval_workspaces_by_lane.resize(lane_count);
  used_outcomes_by_lane.resize(lane_count);

  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    auto *workspace = lanes[lane].workspace;
    for (std::size_t time_slot = 0; time_slot < time_slot_count;
         ++time_slot) {
      const auto time_pos = time_slot * lane_count + compact_lane;
      if (time_slot < workspace->time_values.size()) {
        parent_time_values[time_pos] = workspace->time_values[time_slot];
      }
      parent_time_valid[time_pos] =
          workspace->has_time(static_cast<semantic::Index>(time_slot))
              ? 1U
              : 0U;
    }
    workspaces_by_lane[compact_lane] = workspace;
    eval_workspaces_by_lane[compact_lane] = lanes[lane].eval_workspace;
    used_outcomes_by_lane[compact_lane] = workspace->used_outcomes;
  }
  (void)program;
  return true;
}

inline bool evaluate_compiled_source_node_batch_prepared(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  if (batch_workspace == nullptr ||
      node.source_program_id == semantic::kInvalidIndex) {
    return false;
  }

  const auto value_mask = compiled_math_source_factor_channel_mask(node.kind);
  const auto fill_mask = compiled_math_source_node_fill_mask(node.kind);
  if (value_mask == 0U || fill_mask == 0U) {
    return false;
  }

  auto &compact_lane_indices = batch_workspace->compact_lane_indices;
  auto &active_lanes = batch_workspace->active_lanes;
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &source_channels_by_lane = batch_workspace->source_channels_by_lane;
  auto &parents_by_lane = batch_workspace->parents_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;
  const std::size_t lane_count = compact_lane_indices.size();
  if (lane_count == 0U) {
    return true;
  }

  source_channels_by_lane.resize(lane_count);
  parents_by_lane.resize(lane_count);
  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    if (lane >= condition_evaluators.size() ||
        condition_evaluators[lane] == nullptr ||
        condition_evaluators[lane]->source_channels() == nullptr) {
      return false;
    }
    source_channels_by_lane[compact_lane] =
        condition_evaluators[lane]->source_channels();
    parents_by_lane[compact_lane] = condition_evaluators[lane];
  }

  auto &frame = compiled_math_batch_integral_frame(batch_workspace, 0U);
  frame.ensure_lane_capacity(time_slot_count, lane_count);
  auto &scratch =
      frame.source_product_scratch_layer(0U, lane_count, lane_count);
  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    scratch.times[compact_lane] =
        compiled_math_source_node_time(node, *lanes[lane].workspace);
    scratch.mask_a[compact_lane] = 0U;
    scratch.pdf_a[compact_lane] = 0.0;
    scratch.cdf_a[compact_lane] = 0.0;
    scratch.survival_a[compact_lane] = 1.0;
  }

  const CompiledMathLaneBatchState state{
      lane_count,
      time_slot_count,
      BatchTimeSlotView{
          parent_time_values.data(),
          lane_count,
          parent_time_valid.data()},
      &workspaces_by_lane,
      &source_channels_by_lane,
      &parents_by_lane,
      &eval_workspaces_by_lane,
      &used_outcomes_by_lane};
  const bool supported = compiled_math_batch_source_product_program_fill(
      program,
      node.source_program_id,
      state,
      active_lanes.data(),
      lane_count,
      scratch.times.data(),
      fill_mask,
      &frame,
      1U,
      scratch.mask_a.data(),
      scratch.pdf_a.data(),
      scratch.cdf_a.data(),
      scratch.survival_a.data());
  if (!supported) {
    return false;
  }

  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    const double raw =
        value_mask == kLeafChannelPdf
            ? safe_density(scratch.pdf_a[compact_lane])
            : (value_mask == kLeafChannelCdf
                   ? clamp_probability(scratch.cdf_a[compact_lane])
                   : clamp_probability(scratch.survival_a[compact_lane]));
    lanes[lane].workspace->values[
        static_cast<std::size_t>(node.cache_slot)] = raw;
  }
  return true;
}

inline CompiledMathLaneBatchState prepare_compiled_regular_node_batch_state(
    const CompiledMathProgram &program,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &parents_by_lane = batch_workspace->parents_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;

  const std::size_t lane_count = lanes.size();
  parent_time_values.assign(time_slot_count * lane_count, 0.0);
  parent_time_valid.assign(time_slot_count * lane_count, 0U);
  workspaces_by_lane.assign(lane_count, nullptr);
  parents_by_lane.assign(lane_count, nullptr);
  eval_workspaces_by_lane.assign(lane_count, nullptr);
  used_outcomes_by_lane.assign(lane_count, nullptr);

  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    auto *workspace = lanes[lane].workspace;
    if (workspace != nullptr) {
      for (std::size_t time_slot = 0; time_slot < time_slot_count;
           ++time_slot) {
        const auto time_pos = time_slot * lane_count + lane;
        if (time_slot < workspace->time_values.size()) {
          parent_time_values[time_pos] = workspace->time_values[time_slot];
        }
        parent_time_valid[time_pos] =
            workspace->has_time(static_cast<semantic::Index>(time_slot))
                ? 1U
                : 0U;
      }
      used_outcomes_by_lane[lane] = workspace->used_outcomes;
    }
    workspaces_by_lane[lane] = workspace;
    parents_by_lane[lane] = lanes[lane].parent;
    eval_workspaces_by_lane[lane] = lanes[lane].eval_workspace;
  }

  return CompiledMathLaneBatchState{
      lane_count,
      time_slot_count,
      BatchTimeSlotView{
          parent_time_values.data(),
          lane_count,
          parent_time_valid.data()},
      &workspaces_by_lane,
      nullptr,
      &parents_by_lane,
      &eval_workspaces_by_lane,
      &used_outcomes_by_lane};
}

inline double compiled_math_batch_union_kernel_density_for_lane(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const ExactVariantPlan &plan,
    const bool multi_subset_only) {
  const auto expr_pos = static_cast<std::size_t>(node.subject_id);
  if (expr_pos >= plan.expr_kernels.size()) {
    return 0.0;
  }
  const auto &kernel = plan.expr_kernels[expr_pos];
  if (kernel.children.empty()) {
    return 0.0;
  }
  const auto child_count = kernel.children.size;
  auto child_value = [&](const semantic::Index child_index) {
    if (child_index >= node.children.size) {
      return 0.0;
    }
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + child_index)];
    return workspace.values[static_cast<std::size_t>(child_id)];
  };
  if (!multi_subset_only && child_count == 1U) {
    return child_value(0);
  }

  double total = 0.0;
  for (semantic::Index subset_idx = 0;
       subset_idx < kernel.union_subset_span.size;
       ++subset_idx) {
    const auto &subset = plan.expr_union_subsets[
        static_cast<std::size_t>(
            kernel.union_subset_span.offset + subset_idx)];
    if (multi_subset_only && subset.child_positions.size <= 1U) {
      continue;
    }
    if (subset.child_positions.empty()) {
      continue;
    }
    double subset_value = 0.0;
    if (subset.child_positions.size == 1U) {
      const auto active_pos = plan.expr_union_subset_child_positions[
          static_cast<std::size_t>(subset.child_positions.offset)];
      subset_value = child_value(active_pos);
    } else {
      for (semantic::Index i = 0; i < subset.child_positions.size; ++i) {
        const auto active_pos = plan.expr_union_subset_child_positions[
            static_cast<std::size_t>(subset.child_positions.offset + i)];
        double term = child_value(active_pos);
        if (!std::isfinite(term) || term == 0.0) {
          continue;
        }
        for (semantic::Index j = 0; j < subset.child_positions.size; ++j) {
          if (i == j) {
            continue;
          }
          const auto other_pos = plan.expr_union_subset_child_positions[
              static_cast<std::size_t>(subset.child_positions.offset + j)];
          term *= child_value(child_count + other_pos);
          if (!std::isfinite(term) || term == 0.0) {
            break;
          }
        }
        subset_value += term;
      }
      subset_value = clean_signed_value(subset_value);
    }
    total += static_cast<double>(subset.sign) * subset_value;
  }
  return clean_signed_value(total);
}

inline double compiled_math_batch_union_kernel_cdf_for_lane(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace) {
  double total = 0.0;
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + i)];
    total += workspace.values[static_cast<std::size_t>(child_id)];
  }
  return clamp_probability(total);
}

inline void evaluate_compiled_regular_node_batch(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  if (compiled_math_is_source_value_node(node.kind)) {
    throw std::runtime_error(
        "compiled batch source node reached regular-node evaluator");
  }
  const auto state =
      prepare_compiled_regular_node_batch_state(
          program,
          lanes,
          time_slot_count,
          batch_workspace);
  const auto write_value = [&](const std::size_t lane_pos,
                               const double value) {
    auto *workspace = lanes[lane_pos].workspace;
    if (workspace != nullptr) {
      workspace->values[static_cast<std::size_t>(node.cache_slot)] = value;
    }
  };
  const auto first_child_id = [&]() {
    return program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
  };
  const auto child_value = [&](const CompiledMathWorkspace &workspace,
                               const semantic::Index child_index) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + child_index)];
    return workspace.values[static_cast<std::size_t>(child_id)];
  };

  switch (node.kind) {
  case CompiledMathNodeKind::Constant:
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      write_value(lane_pos, node.constant);
    }
    return;
  case CompiledMathNodeKind::SourcePdf:
  case CompiledMathNodeKind::SourceCdf:
  case CompiledMathNodeKind::SourceSurvival:
    throw std::runtime_error(
        "compiled batch source node reached regular-node evaluator");
  case CompiledMathNodeKind::ExprDensity:
  case CompiledMathNodeKind::ExprCdf:
  case CompiledMathNodeKind::ExprSurvival:
    throw std::runtime_error(
        "compiled math root contains an interpreter expression node");
  case CompiledMathNodeKind::UnionKernelMultiSubsetCdf:
  case CompiledMathNodeKind::IntegralZeroToCurrent:
  case CompiledMathNodeKind::IntegralZeroToCurrentRaw:
    throw std::runtime_error("integral node was not handled by batch path");
  case CompiledMathNodeKind::TimeGate: {
    const auto child_id = first_child_id();
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace == nullptr) {
        continue;
      }
      const auto lane = static_cast<semantic::Index>(lane_pos);
      const double raw =
          workspace->values[static_cast<std::size_t>(child_id)];
      write_value(
          lane_pos,
          state.time(node.time_id, lane) >= state.time(node.aux_id, lane)
              ? raw
              : 0.0);
    }
    return;
  }
  case CompiledMathNodeKind::SimpleGuardDensity:
  case CompiledMathNodeKind::SimpleGuardCdf: {
    const auto child_id = first_child_id();
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      auto *parent = lanes[lane_pos].parent;
      auto *condition_evaluator = condition_evaluators[lane_pos];
      if (workspace == nullptr) {
        continue;
      }
      if (condition_evaluator == nullptr || parent == nullptr) {
        write_value(lane_pos, 0.0);
        continue;
      }
      const auto &kernel = parent->plan().expr_kernels[
          static_cast<std::size_t>(node.subject_id)];
      if (!kernel.simple_event_guard || kernel.has_unless) {
        throw std::runtime_error(
            "compiled simple guard reached a non-simple guard");
      }
      const double raw =
          workspace->values[static_cast<std::size_t>(child_id)];
      if (condition_evaluator->source_order_known_before(
              kernel.guard_blocker_source_id,
              kernel.guard_ref_source_id,
              node.condition_id,
              workspace)) {
        write_value(lane_pos, 0.0);
        continue;
      }
      const auto lane = static_cast<semantic::Index>(lane_pos);
      CompiledTimedUpperBound upper;
      if (!compiled_math_batch_upper_bound_from_terms(
              program,
              node,
              state,
              lane,
              &upper)) {
        throw std::runtime_error(
            "compiled batch simple guard upper bound is unsupported");
      }
      if (!upper.found) {
        write_value(
            lane_pos,
            node.kind == CompiledMathNodeKind::SimpleGuardCdf
                ? clamp_probability(raw)
                : clean_signed_value(raw));
        continue;
      }
      const double current_time = state.time(node.time_id, lane);
      double value = 0.0;
      if (node.kind == CompiledMathNodeKind::SimpleGuardCdf) {
        if (!(upper.normalizer > 0.0)) {
          value = 0.0;
        } else if (current_time >= upper.time) {
          value = 1.0;
        } else {
          value = clamp_probability(raw / upper.normalizer);
        }
      } else {
        if (!(upper.normalizer > 0.0) || current_time >= upper.time) {
          value = 0.0;
        } else {
          value = clean_signed_value(raw / upper.normalizer);
        }
      }
      write_value(lane_pos, value);
    }
    return;
  }
  case CompiledMathNodeKind::UnionKernelDensity:
  case CompiledMathNodeKind::UnionKernelMultiSubsetDensity: {
    const bool multi_subset_only =
        node.kind == CompiledMathNodeKind::UnionKernelMultiSubsetDensity;
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      auto *parent = lanes[lane_pos].parent;
      if (workspace == nullptr) {
        continue;
      }
      write_value(
          lane_pos,
          parent == nullptr
              ? 0.0
              : compiled_math_batch_union_kernel_density_for_lane(
                    program,
                    node,
                    *workspace,
                    parent->plan(),
                    multi_subset_only));
    }
    return;
  }
  case CompiledMathNodeKind::UnionKernelCdf:
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace != nullptr) {
        write_value(
            lane_pos,
            compiled_math_batch_union_kernel_cdf_for_lane(
                program,
                node,
                *workspace));
      }
    }
    return;
  case CompiledMathNodeKind::OutcomeSubsetUnused:
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      auto *parent = lanes[lane_pos].parent;
      if (workspace == nullptr) {
        continue;
      }
      double value = 1.0;
      if (workspace->used_outcomes != nullptr) {
        if (parent == nullptr) {
          throw std::runtime_error(
              "compiled outcome-used gate has no parent source view");
        }
        const auto offset = static_cast<std::size_t>(node.subject_id);
        const auto size = static_cast<std::size_t>(node.aux_id);
        const auto &indices = parent->plan().compiled_outcome_gate_indices;
        if (offset + size > indices.size()) {
          throw std::runtime_error(
              "compiled outcome-used gate points outside the plan");
        }
        for (std::size_t i = 0; i < size; ++i) {
          const auto outcome_idx = indices[offset + i];
          if (outcome_idx != semantic::kInvalidIndex &&
              static_cast<std::size_t>(outcome_idx) <
                  workspace->used_outcomes->size() &&
              (*workspace->used_outcomes)[
                  static_cast<std::size_t>(outcome_idx)] != 0U) {
            value = 0.0;
            break;
          }
        }
      }
      write_value(lane_pos, value);
    }
    return;
  case CompiledMathNodeKind::ExprUpperBoundDensity:
  case CompiledMathNodeKind::ExprUpperBoundCdf: {
    const auto child_id = first_child_id();
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      auto *condition_evaluator = condition_evaluators[lane_pos];
      if (workspace == nullptr) {
        continue;
      }
      const double raw =
          workspace->values[static_cast<std::size_t>(child_id)];
      const auto lane = static_cast<semantic::Index>(lane_pos);
      CompiledTimedUpperBound best_upper;
      if (!compiled_math_batch_expr_upper_bound_for_node(
              program,
              node,
              state,
              lane,
              condition_evaluator,
              &best_upper)) {
        throw std::runtime_error(
            "compiled batch expr upper bound is unsupported");
      }
      if (!best_upper.found) {
        write_value(lane_pos, raw);
        continue;
      }
      const double current_time = state.time(node.time_id, lane);
      double value = 0.0;
      if (node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
        if (!(best_upper.normalizer > 0.0)) {
          value = 0.0;
        } else if (current_time >= best_upper.time) {
          value = 1.0;
        } else {
          value = clamp_probability(raw / best_upper.normalizer);
        }
      } else {
        if (!(best_upper.normalizer > 0.0) ||
            current_time >= best_upper.time) {
          value = 0.0;
        } else {
          value = safe_density(raw / best_upper.normalizer);
        }
      }
      write_value(lane_pos, value);
    }
    return;
  }
  case CompiledMathNodeKind::Product:
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace == nullptr) {
        continue;
      }
      double value = 1.0;
      for (semantic::Index i = 0; i < node.children.size; ++i) {
        value *= child_value(*workspace, i);
        if (!std::isfinite(value) || value == 0.0) {
          value = 0.0;
          break;
        }
      }
      write_value(lane_pos, value);
    }
    return;
  case CompiledMathNodeKind::Sum:
  case CompiledMathNodeKind::CleanSignedSum:
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace == nullptr) {
        continue;
      }
      double total = 0.0;
      for (semantic::Index i = 0; i < node.children.size; ++i) {
        total += child_value(*workspace, i);
      }
      write_value(
          lane_pos,
          node.kind == CompiledMathNodeKind::CleanSignedSum
              ? clean_signed_value(total)
              : total);
    }
    return;
  case CompiledMathNodeKind::ClampProbability: {
    const auto child_id = first_child_id();
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace != nullptr) {
        write_value(
            lane_pos,
            clamp_probability(
                workspace->values[static_cast<std::size_t>(child_id)]));
      }
    }
    return;
  }
  case CompiledMathNodeKind::Complement: {
    const auto child_id = first_child_id();
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace != nullptr) {
        write_value(
            lane_pos,
            clamp_probability(
                1.0 -
                workspace->values[static_cast<std::size_t>(child_id)]));
      }
    }
    return;
  }
  case CompiledMathNodeKind::Negate: {
    const auto child_id = first_child_id();
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace != nullptr) {
        write_value(
            lane_pos,
            clean_signed_value(
                -workspace->values[static_cast<std::size_t>(child_id)]));
      }
    }
    return;
  }
  }
}

inline void evaluate_compiled_node_span_batch(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan schedule,
    const semantic::Index result_node_id,
    const std::vector<CompiledMathBatchLane> &lanes,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane,
    const std::vector<semantic::Index> *schedule_nodes = nullptr,
    const bool workspace_prepared = false) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (result_node_id == semantic::kInvalidIndex || lanes.empty()) {
    return;
  }
  if (batch_workspace == nullptr) {
    throw std::runtime_error("compiled math batch evaluation requires a workspace");
  }
  if (!workspace_prepared) {
    for (const auto &lane : lanes) {
      if (lane.workspace != nullptr) {
        lane.workspace->ensure_size(program);
      }
    }
  }
  auto &condition_evaluators = batch_workspace->condition_evaluators;
  condition_evaluators.assign(lanes.size(), nullptr);
  std::vector<double> integral_values;
  const std::size_t time_slot_count =
      compiled_math_program_time_slot_count(program);
  const auto &nodes =
      schedule_nodes == nullptr ? program.root_schedule_nodes : *schedule_nodes;
  bool source_node_batch_base_prepared = false;
  bool source_node_batch_base_supported = false;
  for (semantic::Index schedule_idx = 0; schedule_idx < schedule.size;
       ++schedule_idx) {
    const auto node_id = nodes[
        static_cast<std::size_t>(schedule.offset + schedule_idx)];
    const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      auto *condition_evaluator = lanes[lane].parent;
      if (node.source_view_id != 0 &&
          node.source_view_id != semantic::kInvalidIndex) {
        if (lanes[lane].eval_workspace == nullptr) {
          throw std::runtime_error(
              "compiled node requires a planned source view but no workspace was supplied");
        }
        condition_evaluator =
            lanes[lane].eval_workspace->source_view_evaluator(
                node.source_view_id,
                lanes[lane].parent);
      }
      condition_evaluators[lane] = condition_evaluator;
    }

    if (compiled_math_is_source_value_node(node.kind) &&
        node.source_program_id != semantic::kInvalidIndex) {
      if (!source_node_batch_base_prepared) {
        source_node_batch_base_supported =
            prepare_compiled_source_node_batch_base(
                program,
                lanes,
                time_slot_count,
                batch_workspace);
        source_node_batch_base_prepared = true;
      }
      if (source_node_batch_base_supported &&
          evaluate_compiled_source_node_batch_prepared(
              program,
              node,
              lanes,
              condition_evaluators,
              time_slot_count,
              batch_workspace)) {
        continue;
      }
      throw std::runtime_error(
          "compiled batch source node has no prepared source-program batch implementation");
    }

    if (compiled_math_is_integral_node(node.kind)) {
      evaluate_compiled_integral_node_batch(
          program,
          node,
          lanes,
          condition_evaluators,
          batch_workspace,
          &integral_values);
      for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
        auto *workspace = lanes[lane].workspace;
        if (workspace != nullptr) {
          workspace->values[static_cast<std::size_t>(node.cache_slot)] =
              integral_values[lane];
        }
      }
      continue;
    }

    evaluate_compiled_regular_node_batch(
        program,
        node,
        lanes,
        condition_evaluators,
        time_slot_count,
        batch_workspace);
  }

  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    auto *workspace = lanes[lane].workspace;
    (*out_by_lane)[lane] =
        workspace == nullptr
            ? 0.0
            : workspace->values[static_cast<std::size_t>(result_node_id)];
  }
}

inline void evaluate_compiled_math_root_batch(
    const CompiledMathProgram &program,
    const semantic::Index root_id,
    const std::vector<CompiledMathBatchLane> &lanes,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (root_id == semantic::kInvalidIndex || lanes.empty()) {
    return;
  }
  const auto &root = program.roots[static_cast<std::size_t>(root_id)];
  evaluate_compiled_node_span_batch(
      program,
      root.schedule,
      root.node_id,
      lanes,
      batch_workspace,
      out_by_lane,
      nullptr,
      true);
}

} // namespace detail
} // namespace accumulatr::eval
