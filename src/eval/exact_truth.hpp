#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "batch_source_product.hpp"
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
          value = compiled_math_source_product_direct_leaf_scalar(
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
  const std::vector<ExactSourceChannels *> *source_channels_by_lane{nullptr};
  const std::vector<CompiledSourceView *> *parents_by_lane{nullptr};
  const std::vector<CompiledEvalWorkspace *> *eval_workspaces_by_lane{nullptr};
  const std::vector<const std::vector<std::uint8_t> *> *used_outcomes_by_lane{
      nullptr};

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

struct CompiledMathBatchIntegralFrame {
  std::vector<double> time_values;
  std::vector<ExactSourceChannels *> source_channels_by_lane;
  std::vector<CompiledSourceView *> parents_by_lane;
  std::vector<CompiledEvalWorkspace *> eval_workspaces_by_lane;
  std::vector<const std::vector<std::uint8_t> *> used_outcomes_by_lane;
  std::vector<semantic::Index> active_lanes;
  std::vector<semantic::Index> parent_lanes;
  std::vector<semantic::Index> term_lanes_a;
  std::vector<semantic::Index> term_lanes_b;
  std::vector<double> weights;
  std::vector<double> values;
  std::vector<double> products;
  std::vector<double> factor_values;
  std::vector<double> source_values;
  std::vector<double> lower_by_lane;
  std::vector<double> upper_by_lane;
  std::vector<std::uint8_t> source_program_valid_mask;
  std::vector<double> source_program_pdf;
  std::vector<double> source_program_cdf;
  std::vector<double> source_program_survival;

  void ensure_lane_capacity(const std::size_t time_slot_count,
                            const std::size_t lane_capacity) {
    const auto time_value_count = time_slot_count * lane_capacity;
    if (time_values.size() < time_value_count) {
      time_values.resize(time_value_count, 0.0);
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
    if (lower_by_lane.size() < lane_capacity) {
      lower_by_lane.resize(lane_capacity, 0.0);
    }
    if (upper_by_lane.size() < lane_capacity) {
      upper_by_lane.resize(lane_capacity, 0.0);
    }
  }

  void reset_source_product_cache(const std::size_t program_count,
                                  const std::size_t lane_capacity) {
    const auto count = program_count * lane_capacity;
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
};

struct CompiledMathBatchWorkspace {
  std::vector<CompiledSourceView *> condition_evaluators;
  std::vector<semantic::Index> active_lanes;
  std::vector<std::size_t> active_lane_indices;
  std::vector<std::size_t> compact_lane_indices;
  std::vector<ExactSourceChannels *> source_channels_by_lane;
  std::vector<CompiledSourceView *> parents_by_lane;
  std::vector<CompiledEvalWorkspace *> eval_workspaces_by_lane;
  std::vector<const std::vector<std::uint8_t> *> used_outcomes_by_lane;
  std::vector<double> parent_time_values;
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
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    return false;
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
    if (node.aux_id != semantic::kInvalidIndex ||
        node.aux2_id != semantic::kInvalidIndex) {
      return false;
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
  return kernel_pos < program.integral_kernels.size() &&
         compiled_math_batch_kernel_supported(
             program,
             program.integral_kernels[kernel_pos],
             depth + 1U);
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

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_batch_source_product_program_fill_lane(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const double current_time,
    const std::uint8_t fill_mask,
    bool *supported);

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_batch_source_product_program_leaf_fill_lane(
    const CompiledMathSourceProductProgram &source_program,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask,
    bool *supported) {
  if (source_channels == nullptr ||
      source_program.leaf_index == semantic::kInvalidIndex) {
    *supported = false;
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  return compiled_math_source_product_direct_leaf_fill(
      source_program.leaf_dist_kind,
      source_channels->source_product_leaf_input(source_program.leaf_index),
      current_time - source_program.leaf_onset_abs_value,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_batch_source_product_program_exact_gate_fill_lane(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const double current_time,
    const std::uint8_t fill_mask,
    bool *supported) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !compiled_math_source_product_bounds_have_overlay(
          source_program.bounds)) {
    return compiled_math_source_product_forced_fill(relation, fill_mask);
  }
  auto *source_channels = state.source_channels(lane);
  if (source_channels == nullptr) {
    *supported = false;
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  const auto bounds =
      source_channels->source_product_resolved_bounds_from_time_slots(
          source_program.source_id,
          source_program.bounds,
          state.time_slots.values,
          state.time_slots.lane_stride,
          lane);
  if (const double *exact_time = bounds.exact_time) {
    if (!(current_time >= *exact_time)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    return compiled_math_source_product_certain_fill(fill_mask);
  }
  return compiled_math_batch_source_product_program_fill_lane(
      program,
      source_program.child_program_id,
      state,
      lane,
      current_time,
      fill_mask,
      supported);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_batch_source_product_program_conditioned_fill_lane(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const double current_time,
    const std::uint8_t fill_mask,
    bool *supported) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  const bool has_condition_overlay =
      compiled_math_source_product_bounds_have_overlay(
          source_program.bounds);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !has_condition_overlay) {
    return compiled_math_source_product_forced_fill(relation, fill_mask);
  }

  auto *source_channels = state.source_channels(lane);
  if (source_channels == nullptr) {
    *supported = false;
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  const auto bounds =
      source_channels->source_product_resolved_bounds_from_time_slots(
          source_program.source_id,
          source_program.bounds,
          state.time_slots.values,
          state.time_slots.lane_stride,
          lane);
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
      compiled_math_batch_source_product_program_fill_lane(
          program,
          source_program.child_program_id,
          state,
          lane,
          current_time,
          uncond_mask,
          supported);
  if (!*supported) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
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
        compiled_math_batch_source_product_program_fill_lane(
            program,
            source_program.child_program_id,
            state,
            lane,
            lower_bound,
            lower_mask,
            supported);
    if (!*supported) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
  }
  if (!std::isfinite(upper_bound)) {
    return compiled_math_source_product_conditionalize(
        uncond,
        lower,
        fill_mask);
  }
  const auto upper =
      compiled_math_batch_source_product_program_fill_lane(
          program,
          source_program.child_program_id,
          state,
          lane,
          upper_bound,
          kLeafChannelCdf,
          supported);
  if (!*supported) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
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
      uncond,
      lower,
      upper,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_batch_source_product_program_fill_lane(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const double current_time,
    const std::uint8_t fill_mask,
    bool *supported) {
  if (source_product_program_id == semantic::kInvalidIndex) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  const auto pos = static_cast<std::size_t>(source_product_program_id);
  if (pos >= program.integral_kernel_source_product_programs.size()) {
    *supported = false;
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  const auto &source_program =
      program.integral_kernel_source_product_programs[pos];
  switch (source_program.kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
    return compiled_math_source_product_impossible_fill(fill_mask);
  case CompiledMathSourceProductProgramKind::ConstantOne:
    return compiled_math_source_product_certain_fill(fill_mask);
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    return compiled_math_batch_source_product_program_leaf_fill_lane(
        source_program,
        state.source_channels(lane),
        current_time,
        fill_mask,
        supported);
  case CompiledMathSourceProductProgramKind::ExactGate:
    return compiled_math_batch_source_product_program_exact_gate_fill_lane(
        program,
        source_program,
        state,
        lane,
        current_time,
        fill_mask,
        supported);
  case CompiledMathSourceProductProgramKind::Conditioned:
    return compiled_math_batch_source_product_program_conditioned_fill_lane(
        program,
        source_program,
        state,
        lane,
        current_time,
        fill_mask,
        supported);
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    *supported = false;
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  *supported = false;
  return compiled_math_source_product_impossible_fill(fill_mask);
}

inline void compiled_math_batch_store_source_product_program_fill(
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t lane_stride,
    const semantic::Index program_id,
    const semantic::Index lane,
    const ExactSourceChannels::SourceProductScalarFill &fill) {
  const auto pos =
      static_cast<std::size_t>(program_id) * lane_stride +
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
    const semantic::Index program_id,
    const semantic::Index lane,
    const std::uint8_t channel_mask) {
  const auto pos =
      static_cast<std::size_t>(program_id) * lane_stride +
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
    out_by_lane[static_cast<std::size_t>(active.lanes[i])] = 1.0;
  }

  for (semantic::Index op_idx = 0; op_idx < source_product_ops.size;
       ++op_idx) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + op_idx)];
    const auto channel_mask = op.value_channel_mask;
    if (channel_mask == 0U) {
      if (op.constant_value == 0.0) {
        for (std::size_t i = 0; i < active.size; ++i) {
          out_by_lane[static_cast<std::size_t>(active.lanes[i])] = 0.0;
        }
        break;
      }
      for (std::size_t i = 0; i < active.size; ++i) {
        const auto lane_pos = static_cast<std::size_t>(active.lanes[i]);
        if (out_by_lane[lane_pos] == 0.0) {
          continue;
        }
        out_by_lane[lane_pos] *=
            op.constant_value;
      }
      continue;
    }

    if (op.source_product_channel_id == semantic::kInvalidIndex ||
        op.source_product_program_id == semantic::kInvalidIndex) {
      return false;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];

    for (std::size_t i = 0; i < active.size; ++i) {
      const auto lane = active.lanes[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      if (out_by_lane[lane_pos] == 0.0) {
        continue;
      }
      double value = 0.0;
      bool supported = true;
      auto *source_channels = state.source_channels(lane);
      if (source_channels == nullptr) {
        return false;
      }
      if (!op.cache_result) {
        const double time =
            batch_source_product_channel_time(
                channel,
                state.time_slots,
                lane);
        if (source_channels->source_product_direct_leaf_available(
                op.source_product_channel_id)) {
          value = compiled_math_source_product_direct_leaf_scalar(
              channel,
              source_channels->source_product_leaf_input(channel.leaf_index),
              time,
              channel_mask);
        } else {
          const auto fill =
              compiled_math_batch_source_product_program_fill_lane(
                  program,
                  op.source_product_program_id,
                  state,
                  lane,
                  time,
                  channel_mask,
                  &supported);
          if (!supported) {
            return false;
          }
          value = compiled_math_source_product_fill_value(fill, channel_mask);
        }
      } else {
        const auto program_pos =
            static_cast<std::size_t>(op.source_product_program_id) *
            state.lane_count +
            lane_pos;
        if ((frame->source_program_valid_mask[program_pos] &
             channel_mask) == 0U) {
          const double time =
              batch_source_product_channel_time(
                  channel,
                  state.time_slots,
                  lane);
          const auto fill =
              source_channels->source_product_direct_leaf_available(
                  op.source_product_channel_id)
                  ? compiled_math_source_product_direct_leaf_fill(
                        channel,
                        source_channels->source_product_leaf_input(
                            channel.leaf_index),
                        time,
                        op.fill_channel_mask)
                  : compiled_math_batch_source_product_program_fill_lane(
                        program,
                        op.source_product_program_id,
                        state,
                        lane,
                        time,
                        op.fill_channel_mask,
                        &supported);
          if (!supported) {
            return false;
          }
          compiled_math_batch_store_source_product_program_fill(
              frame,
              state.lane_count,
              op.source_product_program_id,
              lane,
              fill);
        }
        value = compiled_math_batch_source_product_cached_value(
            *frame,
            state.lane_count,
            op.source_product_program_id,
            lane,
            channel_mask);
      }

      const double product = out_by_lane[lane_pos] * value;
      out_by_lane[lane_pos] =
          std::isfinite(product) && product != 0.0 ? product : 0.0;
    }
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

inline bool compiled_math_batch_apply_expr_upper_factor(
    const CompiledMathProgram &program,
    const CompiledMathIntegralExprUpperFactor &factor,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    double *product) {
  const auto &node =
      program.nodes[static_cast<std::size_t>(factor.node_id)];
  if (node.kind != CompiledMathNodeKind::ExprUpperBoundDensity &&
      node.kind != CompiledMathNodeKind::ExprUpperBoundCdf) {
    return false;
  }
  if (node.aux_id != semantic::kInvalidIndex ||
      node.aux2_id != semantic::kInvalidIndex) {
    return false;
  }
  CompiledTimedUpperBound upper;
  auto *condition_evaluator =
      compiled_math_node_evaluator(
          node,
          state.parent(lane),
          state.eval_workspace(lane));
  ExactTimedExprUpperBound sequence_upper;
  if (condition_evaluator != nullptr &&
      condition_evaluator->expr_upper_bound_for(
          node.subject_id,
          &sequence_upper) &&
      std::isfinite(sequence_upper.time) &&
      sequence_upper.normalizer > 0.0) {
    upper.found = true;
    upper.time = sequence_upper.time;
    upper.normalizer = sequence_upper.normalizer;
  }
  const bool has_upper = upper.found && upper.normalizer > 0.0;
  const double current_time = state.time(node.time_id, lane);
  if (factor.mode == CompiledMathIntegralExprUpperMode::AfterOne) {
    if (!has_upper || current_time < upper.time) {
      *product = 0.0;
      return true;
    }
    return true;
  }
  if (!has_upper) {
    return true;
  }
  if (current_time >= upper.time) {
    *product = 0.0;
    return true;
  }
  *product /= upper.normalizer;
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

  const auto time_slot_count =
      std::max(
          parent_state.time_slot_count,
          static_cast<std::size_t>(kernel.bind_time_id) + 1U);
  auto &frame = compiled_math_batch_integral_frame(batch_workspace, depth);
  const auto child_capacity =
      parent_active.size * quadrature::kDefaultFiniteOrder;
  frame.ensure_lane_capacity(time_slot_count, child_capacity);
  frame.reset_source_product_cache(
      program.integral_kernel_source_product_programs.size(),
      child_capacity);

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
      }
      frame.time_values[
          static_cast<std::size_t>(kernel.bind_time_id) * child_capacity +
          child_count] = shift + scale * rule.nodes[q];
      ++child_count;
    }
  }
  if (child_count == 0U) {
    return true;
  }

  const CompiledMathLaneBatchState child_state{
      child_capacity,
      time_slot_count,
      BatchTimeSlotView{frame.time_values.data(), child_capacity},
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
        for (std::size_t lane_pos = 0; lane_pos < term_count; ++lane_pos) {
          const auto lane = (*current_lanes)[lane_pos];
          const auto lane_index = static_cast<std::size_t>(lane);
          if (!compiled_math_batch_apply_expr_upper_factor(
                  program,
                  factor,
                  child_state,
                  lane,
                  &frame.products[lane_index])) {
            return false;
          }
          if (frame.products[lane_index] != 0.0) {
            (*next_lanes)[next_count++] = lane;
          }
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
      kernel.kind != CompiledMathIntegralKernelKind::SourceProductSum ||
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
    if (workspace->time_values.size() < time_slot_count) {
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
      parent_time_values[time_slot * compact_count + compact_lane] =
          workspace->time_values[time_slot];
    }
    source_channels_by_lane[compact_lane] =
        lanes[lane].parent->source_channels();
    parents_by_lane[compact_lane] = lanes[lane].parent;
    eval_workspaces_by_lane[compact_lane] = lanes[lane].eval_workspace;
    used_outcomes_by_lane[compact_lane] = workspace->used_outcomes;
  }

  const CompiledMathLaneBatchState parent_state{
      compact_count,
      time_slot_count,
      BatchTimeSlotView{parent_time_values.data(), compact_count},
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
    double value = evaluate_compiled_integral_kernel(
        program,
        node,
        workspace,
        lanes[lane].parent,
        lanes[lane].scenario,
        lanes[lane].eval_workspace,
        0.0,
        current_time);
    if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent) {
      value = clamp_probability(value);
    } else if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw) {
      value = std::isfinite(value) ? clean_signed_value(value) : 0.0;
    } else {
      value = std::isfinite(value) ? value : 0.0;
    }
    (*node_values)[lane] = value;
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
  const auto &nodes =
      schedule_nodes == nullptr ? program.root_schedule_nodes : *schedule_nodes;
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

    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      auto *workspace = lanes[lane].workspace;
      auto *parent = lanes[lane].parent;
      auto *condition_evaluator = condition_evaluators[lane];
      if (workspace == nullptr) {
        continue;
      }
      double value = 0.0;
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
      case CompiledMathNodeKind::UnionKernelMultiSubsetCdf:
      case CompiledMathNodeKind::IntegralZeroToCurrent:
      case CompiledMathNodeKind::IntegralZeroToCurrentRaw:
        throw std::runtime_error("integral node was not handled by batch path");
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
      workspace->values[static_cast<std::size_t>(node.cache_slot)] = value;
    }
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
