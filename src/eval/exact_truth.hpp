#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "exact_source_channels.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

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

  leaf::EventChannels source_channels_for_slot(
      const semantic::Index source_channel_slot,
      const double t,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask = kLeafChannelAll) const {
    const auto slot = static_cast<std::size_t>(source_channel_slot);
    if (slot >= plan_.compiled_math.source_channel_plans.size()) {
      return leaf::EventChannels::impossible();
    }
    const auto &channel_plan =
        plan_.compiled_math.source_channel_plans[slot];
    const auto &request = channel_plan.request;
    const auto relation =
        exact_compiled_source_view_relation(
            plan_, source_view_id_, request.source_id);
    if (relation != ExactRelation::Unknown &&
        !channel_plan.has_source_condition_overlay) {
      return forced_channels(relation);
    }
    return source_channels_->evaluate_source_channel_plan_at(
        channel_plan, t, workspace, channel_mask);
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

inline semantic::Index compiled_source_channel_condition_cache_id(
    const CompiledMathProgram &program,
    const CompiledSourceChannelPlan &channel_plan,
    const CompiledMathWorkspace &workspace) {
  if (!channel_plan.source_view_condition_cache_dynamic) {
    return channel_plan.source_view_condition_static_cache_id;
  }
  const auto source_view_id =
      channel_plan.request.source_view_id == semantic::kInvalidIndex
          ? 0
          : channel_plan.request.source_view_id;
  std::size_t seed = 0x9e3779b97f4a7c15ULL;
  exact_hash_condition_part(&seed, static_cast<std::size_t>(source_view_id));
  exact_hash_condition_part(
      &seed, static_cast<std::size_t>(channel_plan.request.condition_id));
  for (semantic::Index i = 0;
       i < channel_plan.condition_time_dependencies.size;
       ++i) {
    const auto time_id =
        program.condition_cache_time_dependencies[
            static_cast<std::size_t>(
                channel_plan.condition_time_dependencies.offset + i)];
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

inline leaf::EventChannels compiled_math_store_source_channels(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace) {
  if (workspace == nullptr) {
    return leaf::EventChannels{0.0, 0.0, 0.0};
  }
  const double time = compiled_math_source_node_time(node, *workspace);
  const auto slot = static_cast<std::size_t>(node.source_channel_slot);
  if (node.source_channel_slot == semantic::kInvalidIndex ||
      slot >= workspace->source_channel_valid.size()) {
    throw std::runtime_error(
        "compiled source node reached evaluation without a source-channel slot");
  }
  const auto &channel_plan = program.source_channel_plans[slot];
  const auto condition_cache_id =
      evaluator == nullptr
          ? 0
          : compiled_source_channel_condition_cache_id(
                program, channel_plan, *workspace);
  if (workspace->source_channel_valid[slot] != 0U &&
      workspace->source_channel_evaluators[slot] == evaluator &&
      workspace->source_channel_condition_ids[slot] == condition_cache_id &&
      workspace->source_channel_times[slot] == time) {
    return leaf::EventChannels{
        workspace->source_channel_pdf[slot],
        workspace->source_channel_cdf[slot],
        workspace->source_channel_survival[slot]};
  }

  const auto channels =
      evaluator == nullptr
          ? leaf::EventChannels{0.0, 0.0, 0.0}
          : evaluator->source_channels_for_slot(
                node.source_channel_slot,
                time,
                workspace,
                channel_plan.required_channels);
  workspace->source_channel_valid[slot] = 1U;
  workspace->source_channel_evaluators[slot] = evaluator;
  workspace->source_channel_condition_ids[slot] = condition_cache_id;
  workspace->source_channel_times[slot] = time;
  workspace->source_channel_pdf[slot] = channels.pdf;
  workspace->source_channel_cdf[slot] = channels.cdf;
  workspace->source_channel_survival[slot] = channels.survival;
  return channels;
}

inline leaf::EventChannels compiled_math_source_channels(
    const CompiledMathNode &node,
    const CompiledMathWorkspace *workspace) {
  if (workspace == nullptr ||
      node.source_channel_slot == semantic::kInvalidIndex) {
    throw std::runtime_error(
        "compiled source value has no planned source-channel slot");
  }
  const auto slot = static_cast<std::size_t>(node.source_channel_slot);
  if (slot >= workspace->source_channel_valid.size() ||
      workspace->source_channel_valid[slot] == 0U) {
    throw std::runtime_error(
        "compiled source value reached before source-channel load");
  }
  return leaf::EventChannels{
      workspace->source_channel_pdf[slot],
      workspace->source_channel_cdf[slot],
      workspace->source_channel_survival[slot]};
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
      need_pdf ? R::dlnorm(x, m, s, 0) : 0.0,
      need_cdf ? R::plnorm(x, m, s, 1, 0) : 0.0,
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
  const double scale = 1.0 / loaded.params[1];
  return compiled_math_source_product_finish_base_fill(
      need_pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
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
    const CompiledMathSourceProductChannel &channel,
    const ExactLoadedLeafInput &loaded,
    const double time,
    const std::uint8_t fill_mask) {
  const double x = time - channel.leaf_onset_abs_value - loaded.t0;
  if (!(x > 0.0)) {
    ExactSourceChannels::SourceProductScalarFill fill;
    fill.mask = fill_mask;
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0;
    }
    return fill;
  }
  switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
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
        need_pdf ? R::dlnorm(x, m, s, 0) : 0.0,
        need_pdf ? 0.0 : R::plnorm(x, m, s, 1, 0),
        loaded.q,
        channel_mask);
  }
  case leaf::DistKind::Gamma: {
    const double shape = loaded.params[0];
    const double scale = 1.0 / loaded.params[1];
    return compiled_math_source_product_finish_base_scalar(
        need_pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
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

inline void compiled_math_store_source_product_fill(
    CompiledMathWorkspace *workspace,
    const std::size_t channel_pos,
    const ExactSourceChannels::SourceProductScalarFill &fill) {
  if (workspace->source_product_channel_epoch[channel_pos] !=
      workspace->source_product_channel_current_epoch) {
    workspace->source_product_channel_epoch[channel_pos] =
        workspace->source_product_channel_current_epoch;
    workspace->source_product_channel_valid_mask[channel_pos] = 0U;
  }
  workspace->source_product_channel_valid_mask[channel_pos] |= fill.mask;
  if ((fill.mask & kLeafChannelPdf) != 0U) {
    workspace->source_product_channel_pdf[channel_pos] = fill.pdf;
  }
  if ((fill.mask & kLeafChannelCdf) != 0U) {
    workspace->source_product_channel_cdf[channel_pos] = fill.cdf;
  }
  if ((fill.mask & kLeafChannelSurvival) != 0U) {
    workspace->source_product_channel_survival[channel_pos] = fill.survival;
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

inline double compiled_math_source_product_value_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const semantic::Index parent_source_view_id) {
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
          value = compiled_math_source_product_direct_leaf_scalar(
              channel,
              source_channels->source_product_leaf_input(channel.leaf_index),
              time,
              channel_mask);
        } else {
          const auto fill =
              source_channels->evaluate_source_product_channel_fill(
                  channel,
                  channel_mask,
                  time,
                  parent_source_view_id,
                  workspace);
          value = compiled_math_source_product_fill_value(fill, channel_mask);
        }
      } else {
        if (workspace->source_product_channel_epoch[channel_pos] !=
                workspace->source_product_channel_current_epoch ||
            (workspace->source_product_channel_valid_mask[channel_pos] &
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
                  : source_channels->evaluate_source_product_channel_fill(
                        channel,
                        fill_mask,
                        time,
                        parent_source_view_id,
                        workspace);
          compiled_math_store_source_product_fill(
              workspace,
              channel_pos,
              fill);
        }
        if (channel_mask == kLeafChannelPdf) {
          value = workspace->source_product_channel_pdf[channel_pos];
        } else if (channel_mask == kLeafChannelCdf) {
          value = workspace->source_product_channel_cdf[channel_pos];
        } else {
          value = workspace->source_product_channel_survival[channel_pos];
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
  if (!(upper > lower)) {
    return 0.0;
  }
  auto *source_channels = parent == nullptr ? nullptr : parent->source_channels();
  if (source_channels == nullptr) {
    return 0.0;
  }
  const auto parent_source_view_id =
      parent == nullptr ? semantic::Index{0} : parent->source_view_id();
  const auto source_product_source_view_id =
      kernel.source_view_id == 0 || kernel.source_view_id == semantic::kInvalidIndex
          ? parent_source_view_id
          : kernel.source_view_id;

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
    integral_workspace.reset_source_product_channel_cache();
    double value = 0.0;
    if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
      value = compiled_math_source_product_value_for_ops(
            program,
            kernel.source_product_ops,
            &integral_workspace,
            source_channels,
            source_product_source_view_id);
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
        product *= compiled_math_source_product_value_for_ops(
            program,
            term.source_product_ops,
            &integral_workspace,
            source_channels,
            source_product_source_view_id);
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
    case CompiledMathNodeKind::SourceChannelLoad:
      (void)compiled_math_store_source_channels(
          program, node, condition_evaluator, workspace);
      value = 0.0;
      break;
    case CompiledMathNodeKind::SourcePdf:
      value = safe_density(
          workspace->source_channel_pdf[
              static_cast<std::size_t>(node.source_channel_slot)]);
      break;
    case CompiledMathNodeKind::SourceCdf:
      value = clamp_probability(
          workspace->source_channel_cdf[
              static_cast<std::size_t>(node.source_channel_slot)]);
      break;
    case CompiledMathNodeKind::SourceSurvival:
      value = clamp_probability(
          workspace->source_channel_survival[
              static_cast<std::size_t>(node.source_channel_slot)]);
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

} // namespace detail
} // namespace accumulatr::eval
