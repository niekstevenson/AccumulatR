#pragma once

#include <algorithm>
#include <array>
#include <limits>

#include "exact_source_channels.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool relation_view_with_overlay(const RelationView &base,
                                       const ExactRelationTemplate &overlay,
                                       RelationView *out) {
  for (std::size_t i = 0; i < overlay.source_ids.size(); ++i) {
    const auto source_id = overlay.source_ids[i];
    const auto relation = overlay.relations[i];
    const auto inherited = base.relation_for(source_id);
    if (inherited != ExactRelation::Unknown && inherited != relation) {
      return false;
    }
  }
  if (out != nullptr) {
    *out = overlay.empty() ? base : base.with_overlay(&overlay);
  }
  return true;
}

class CompiledSourceView {
public:
  explicit CompiledSourceView(const ExactVariantPlan &plan)
      : plan_(plan) {}

  void reset(ExactSourceChannels *source_channels,
             const RelationView relation_view = {}) {
    source_channels_ = source_channels;
    relation_view_ = relation_view;
  }

  void invalidate_cache() noexcept {}

  ExactSourceChannels *source_channels() const {
    return source_channels_;
  }

  const ExactVariantPlan &plan() const {
    return plan_;
  }

  const RelationView &relation_view() const {
    return relation_view_;
  }

  leaf::EventChannels source_channels_for_slot(
      const semantic::Index source_channel_slot,
      const double t,
      const CompiledMathWorkspace *workspace) const {
    const auto slot = static_cast<std::size_t>(source_channel_slot);
    if (slot >= plan_.compiled_math.source_channel_plans.size()) {
      return leaf::EventChannels::impossible();
    }
    const auto &channel_plan =
        plan_.compiled_math.source_channel_plans[slot];
    const auto &request = channel_plan.request;
    const auto relation = relation_view_.relation_for(request.source_id);
    if (relation != ExactRelation::Unknown &&
        !channel_plan.has_source_condition_overlay) {
      return forced_channels(relation);
    }
    return source_channels_->evaluate_source_channel_plan_at(
        channel_plan, t, workspace);
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

  bool has_guard_upper_bound_overlay(
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    (void)workspace;
    return source_channels_->has_guard_upper_bound_overlay(
        condition_id);
  }

private:
  bool relation_view_knows_before(const semantic::Index before_source_id,
                                  const semantic::Index after_source_id) const {
    if (before_source_id == semantic::kInvalidIndex ||
        after_source_id == semantic::kInvalidIndex ||
        before_source_id == after_source_id) {
      return false;
    }
    const auto before_relation = relation_view_.relation_for(before_source_id);
    const auto after_relation = relation_view_.relation_for(after_source_id);
    return before_relation != ExactRelation::Unknown &&
           after_relation != ExactRelation::Unknown &&
           before_relation < after_relation;
  }

public:
  semantic::Index condition_cache_id() const noexcept {
    return 0;
  }

private:
  const ExactVariantPlan &plan_;
  ExactSourceChannels *source_channels_{nullptr};
  RelationView relation_view_;
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
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace) {
  const auto evaluator_id = exact_condition_cache_id(evaluator);
  if (node.condition_id == 0 ||
      node.condition_id == semantic::kInvalidIndex) {
    return evaluator_id;
  }
  std::size_t seed = 0x9e3779b97f4a7c15ULL;
  exact_hash_condition_part(&seed, static_cast<std::size_t>(evaluator_id));
  exact_hash_condition_part(&seed, static_cast<std::size_t>(node.condition_id));
  const auto condition_pos = static_cast<std::size_t>(node.condition_id);
  if (condition_pos >= workspace.condition_time_dependency_spans.size()) {
    return static_cast<semantic::Index>(
        (seed & static_cast<std::size_t>(0x3fffffffU)) + 1U);
  }
  const auto span = workspace.condition_time_dependency_spans[condition_pos];
  for (semantic::Index i = 0; i < span.size; ++i) {
    const auto time_id =
        workspace.condition_time_dependency_ids[
            static_cast<std::size_t>(span.offset + i)];
    const auto time_pos = static_cast<std::size_t>(time_id);
    if (time_pos >= workspace.time_valid.size() ||
        workspace.time_valid[time_pos] == 0U) {
      continue;
    }
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
    source_view_parents.assign(source_view_count, nullptr);
    source_view_condition_ids.assign(source_view_count, 0);
    source_view_valid.assign(source_view_count, 0U);
  }

  void reset(ExactSourceChannels *source_channels,
             const RelationView relation_view = {}) {
    evaluator.reset(source_channels, relation_view);
    compiled_math.reset_cache();
    reset_source_views();
  }

  void reset_planned_caches() {}

  void reset_source_views() {
    std::fill(source_view_valid.begin(), source_view_valid.end(), 0U);
  }

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
    const auto condition_id = exact_condition_cache_id(parent);
    if (source_view_valid[pos] != 0U &&
        source_view_parents[pos] == parent &&
        source_view_condition_ids[pos] == condition_id) {
      return &source_view_evaluators[pos];
    }

    RelationView view;
    if (!relation_view_with_overlay(
            parent->relation_view(),
            parent->plan().compiled_source_views[pos],
            &view)) {
      source_view_valid[pos] = 0U;
      return nullptr;
    }
    auto &evaluator_for_view = source_view_evaluators[pos];
    evaluator_for_view.reset(parent->source_channels(), view);
    source_view_parents[pos] = parent;
    source_view_condition_ids[pos] = condition_id;
    source_view_valid[pos] = 1U;
    return &evaluator_for_view;
  }

  CompiledSourceView evaluator;
  CompiledMathWorkspace compiled_math;
  std::vector<CompiledSourceView> source_view_evaluators;
  std::vector<const CompiledSourceView *> source_view_parents;
  std::vector<semantic::Index> source_view_condition_ids;
  std::vector<std::uint8_t> source_view_valid;
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
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value);

inline void compiled_math_store_node_cache_entry(
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    double value);

inline double compiled_math_node_time(const CompiledMathNode &node,
                                      const CompiledMathWorkspace &workspace) {
  return workspace.time(node.time_id);
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

inline const CompiledMathCondition *compiled_math_condition_by_id(
    const CompiledMathProgram &program,
    const semantic::Index condition_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return nullptr;
  }
  const auto pos = static_cast<std::size_t>(condition_id - 1U);
  if (pos >= program.conditions.size()) {
    return nullptr;
  }
  return &program.conditions[pos];
}

inline double compiled_math_fact_time(
    const CompiledMathCondition &condition,
    const std::size_t fact_pos,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace) {
  if (fact_pos >= condition.fact_time_ids.size()) {
    return compiled_math_node_time(node, workspace);
  }
  const auto time_id = condition.fact_time_ids[fact_pos];
  if (time_id != semantic::kInvalidIndex && workspace.has_time(time_id)) {
    return workspace.time(time_id);
  }
  return compiled_math_node_time(node, workspace);
}

inline double compiled_math_fact_normalizer(
    const CompiledMathCondition &condition,
    const std::size_t fact_pos,
    const CompiledMathWorkspace &workspace) {
  if (fact_pos >= condition.fact_normalizer_node_ids.size()) {
    return 0.0;
  }
  const auto node_id = condition.fact_normalizer_node_ids[fact_pos];
  if (node_id == semantic::kInvalidIndex) {
    return 0.0;
  }
  const auto pos = static_cast<std::size_t>(node_id);
  if (pos >= workspace.values.size()) {
    return 0.0;
  }
  return workspace.values[pos];
}

inline CompiledTimedUpperBound compiled_math_upper_bound_from_facts(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const std::vector<semantic::Index> &fact_indices) {
  CompiledTimedUpperBound best;
  const auto *condition =
      compiled_math_condition_by_id(program, node.condition_id);
  if (condition == nullptr ||
      node.aux_id == semantic::kInvalidIndex ||
      node.aux2_id == semantic::kInvalidIndex) {
    return best;
  }
  const auto offset = static_cast<std::size_t>(node.aux_id);
  const auto size = static_cast<std::size_t>(node.aux2_id);
  if (offset + size > fact_indices.size()) {
    throw std::runtime_error(
        "compiled upper-bound node points outside condition fact slots");
  }
  for (std::size_t i = 0; i < size; ++i) {
    const auto fact_pos = static_cast<std::size_t>(fact_indices[offset + i]);
    const double fact_time =
        compiled_math_fact_time(*condition, fact_pos, node, workspace);
    const double normalizer =
        compiled_math_fact_normalizer(*condition, fact_pos, workspace);
    if (std::isfinite(fact_time) && normalizer > 0.0 &&
        (!best.found || fact_time < best.time)) {
      best.found = true;
      best.time = fact_time;
      best.normalizer = normalizer;
    }
  }
  return best;
}

inline bool compiled_math_load_node_cache(
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value) {
  if (evaluator == nullptr || !compiled_math_node_cacheable(node.kind)) {
    return false;
  }
  return compiled_math_load_node_cache_entry(
      node, evaluator, workspace, value);
}

inline bool compiled_math_load_node_cache_entry(
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
          compiled_node_condition_cache_id(node, evaluator, workspace) ||
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
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    const double value) {
  if (evaluator == nullptr || !compiled_math_node_cacheable(node.kind)) {
    return;
  }
  compiled_math_store_node_cache_entry(node, evaluator, workspace, value);
}

inline void compiled_math_store_node_cache_entry(
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
      compiled_node_condition_cache_id(node, evaluator, *workspace);
  workspace->cache_times[slot] = compiled_math_node_time(node, *workspace);
  workspace->cache_values[slot] = value;
}

inline leaf::EventChannels compiled_math_store_source_channels(
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

  const auto condition_cache_id =
      evaluator == nullptr
          ? 0
          : compiled_node_condition_cache_id(node, evaluator, *workspace);
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
                node.source_channel_slot, time, workspace);
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

inline double evaluate_compiled_lazy_child(
    const CompiledMathProgram &program,
    const semantic::Index node_id,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace) {
  if (node_id == semantic::kInvalidIndex) {
    return 0.0;
  }
  const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
  switch (node.kind) {
  case CompiledMathNodeKind::Constant:
    return node.constant;
  case CompiledMathNodeKind::RootValue: {
    auto *condition_evaluator =
        compiled_math_node_evaluator(node, parent, eval_workspace);
    if (node.subject_id == semantic::kInvalidIndex ||
        condition_evaluator == nullptr) {
      return 0.0;
    }
    const double value =
        evaluate_compiled_math_root(
            program,
            node.subject_id,
            workspace,
            condition_evaluator,
            scenario,
            eval_workspace);
    return std::isfinite(value) ? value : 0.0;
  }
  case CompiledMathNodeKind::LazyProduct: {
    double product = 1.0;
    for (semantic::Index i = 0; i < node.children.size; ++i) {
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset + i)];
      const double child =
          evaluate_compiled_lazy_child(
              program,
              child_id,
              workspace,
              parent,
              scenario,
              eval_workspace);
      product *= child;
      if (!std::isfinite(product) || product == 0.0) {
        return 0.0;
      }
    }
    return product;
  }
  default:
    throw std::runtime_error(
        "lazy compiled product contains a non-root child");
  }
}

inline double compiled_math_source_value_for_node(
    const CompiledMathNode &node,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *eval_workspace) {
  auto *condition_evaluator =
      compiled_math_node_evaluator(node, parent, eval_workspace);
  if (condition_evaluator == nullptr) {
    return 0.0;
  }
  const auto channels =
      compiled_math_store_source_channels(node, condition_evaluator, workspace);
  switch (node.kind) {
  case CompiledMathNodeKind::SourcePdf:
    return safe_density(channels.pdf);
  case CompiledMathNodeKind::SourceCdf:
    return clamp_probability(channels.cdf);
  case CompiledMathNodeKind::SourceSurvival:
    return clamp_probability(channels.survival);
  default:
    break;
  }
  throw std::runtime_error("source-product integral contains a non-source node");
}

inline double evaluate_compiled_source_product_integral_kernel(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *eval_workspace,
    const double lower,
    const double upper) {
  auto &integral_workspace = workspace->integral_workspace_for(program);
  return quadrature::integrate_finite_default(
      [&](const double u) {
        CompiledMathWorkspace::TimeBinding bind(
            &integral_workspace,
            kernel.bind_time_id,
            u);
        if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
          double product = 1.0;
          for (semantic::Index i = 0; i < kernel.source_value_nodes.size; ++i) {
            const auto node_id =
                program.integral_kernel_source_value_nodes[
                    static_cast<std::size_t>(
                        kernel.source_value_nodes.offset + i)];
            const auto &source_node =
                program.nodes[static_cast<std::size_t>(node_id)];
            product *= compiled_math_source_value_for_node(
                source_node,
                &integral_workspace,
                parent,
                eval_workspace);
            if (!std::isfinite(product) || product == 0.0) {
              return 0.0;
            }
          }
          return std::isfinite(product) ? product : 0.0;
        }
        double total = 0.0;
        for (semantic::Index term_idx = 0;
             term_idx < kernel.source_product_terms.size;
             ++term_idx) {
          const auto &term =
              program.integral_kernel_source_product_terms[
                  static_cast<std::size_t>(
                      kernel.source_product_terms.offset + term_idx)];
          double product = term.sign;
          for (semantic::Index i = 0; i < term.source_value_nodes.size; ++i) {
            const auto node_id =
                program.integral_kernel_source_value_nodes[
                    static_cast<std::size_t>(
                        term.source_value_nodes.offset + i)];
            const auto &source_node =
                program.nodes[static_cast<std::size_t>(node_id)];
            product *= compiled_math_source_value_for_node(
                source_node,
                &integral_workspace,
                parent,
                eval_workspace);
            if (!std::isfinite(product) || product == 0.0) {
              product = 0.0;
              break;
            }
          }
          total += product;
        }
        if (!std::isfinite(total)) {
          return 0.0;
        }
        return kernel.clean_signed_source_sum ? clean_signed_value(total)
                                              : total;
      },
      lower,
      upper);
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
  if (kernel.kind == CompiledMathIntegralKernelKind::Subgraph) {
    auto &integral_workspace = workspace->integral_workspace_for(program);
    return quadrature::integrate_finite_default(
        [&](const double u) {
          CompiledMathWorkspace::TimeBinding bind(
              &integral_workspace,
              kernel.bind_time_id,
              u);
          const double term = evaluate_compiled_node_span(
              program,
              kernel.subgraph,
              kernel.result_node_id,
              &integral_workspace,
              parent,
              scenario,
              eval_workspace,
              nullptr,
              true);
          return std::isfinite(term) ? term : 0.0;
        },
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
    auto *condition_evaluator =
        compiled_math_node_evaluator(node, parent, eval_workspace);
    double value = 0.0;
    if (compiled_math_load_node_cache(
            node,
            condition_evaluator,
            *workspace,
            &value)) {
      workspace->values[static_cast<std::size_t>(node.cache_slot)] = value;
      continue;
    }
    switch (node.kind) {
    case CompiledMathNodeKind::Constant:
      value = node.constant;
      break;
    case CompiledMathNodeKind::SourceChannelLoad:
      (void)compiled_math_store_source_channels(
          node, condition_evaluator, workspace);
      value = 0.0;
      break;
    case CompiledMathNodeKind::SourcePdf:
      value = safe_density(
          compiled_math_source_channels(node, workspace).pdf);
      break;
    case CompiledMathNodeKind::SourceCdf:
      value = clamp_probability(
          compiled_math_source_channels(node, workspace).cdf);
      break;
    case CompiledMathNodeKind::SourceSurvival:
      value = clamp_probability(
          compiled_math_source_channels(node, workspace).survival);
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
      const auto *condition =
          compiled_math_condition_by_id(program, node.condition_id);
      const auto upper =
          condition == nullptr
              ? CompiledTimedUpperBound{}
              : compiled_math_upper_bound_from_facts(
                    program,
                    node,
                    *workspace,
                    condition->guard_upper_fact_lookup.fact_indices);
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
      const auto *condition =
          compiled_math_condition_by_id(program, node.condition_id);
      const auto upper =
          condition == nullptr
              ? CompiledTimedUpperBound{}
              : compiled_math_upper_bound_from_facts(
                    program,
                    node,
                    *workspace,
                    condition->guard_upper_fact_lookup.fact_indices);
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
    case CompiledMathNodeKind::RootValue: {
      if (node.subject_id == semantic::kInvalidIndex ||
          condition_evaluator == nullptr) {
        value = 0.0;
        break;
      }
      value = evaluate_compiled_math_root(
          program,
          node.subject_id,
          workspace,
          condition_evaluator,
          scenario,
          eval_workspace);
      value = std::isfinite(value) ? value : 0.0;
      break;
    }
    case CompiledMathNodeKind::LazyProduct: {
      value = 1.0;
      for (semantic::Index i = 0; i < node.children.size; ++i) {
        const auto child_id = program.child_nodes[
            static_cast<std::size_t>(node.children.offset + i)];
        value *= evaluate_compiled_lazy_child(
            program,
            child_id,
            workspace,
            parent,
            scenario,
            eval_workspace);
        if (!std::isfinite(value) || value == 0.0) {
          value = 0.0;
          break;
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
      const auto *condition =
          compiled_math_condition_by_id(program, node.condition_id);
      const auto upper =
          condition == nullptr
              ? CompiledTimedUpperBound{}
              : compiled_math_upper_bound_from_facts(
                    program,
                    node,
                    *workspace,
                    condition->expr_upper_fact_indices);
      if (!upper.found) {
        value = raw;
        break;
      }
      const double current_time = compiled_math_node_time(node, *workspace);
      if (node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
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
          value = safe_density(raw / upper.normalizer);
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
    compiled_math_store_node_cache(
        node,
        condition_evaluator,
        workspace,
        value);
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
      eval_workspace);
}

} // namespace detail
} // namespace accumulatr::eval
