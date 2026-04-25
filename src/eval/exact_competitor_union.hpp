#pragma once

#include <algorithm>
#include <deque>
#include <memory>
#include <optional>

#include "exact_truth.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactRuntimePreparedCompetitorSubset {
  semantic::Index block_index{0};
  semantic::Index subset_index{0};
  double fixed_mass{0.0};
};

struct ExactRuntimePreparedCompetitorCache {
  std::vector<semantic::Index> block_offsets;
  std::vector<ExactRuntimePreparedCompetitorSubset> subsets;

  void clear() {
    block_offsets.clear();
    subsets.clear();
  }
};

struct ExactStepWorkspace {
  explicit ExactStepWorkspace(const ExactVariantPlan &plan)
      : oracle(plan),
        target_evaluator(plan),
        competitor_evaluator(plan),
        target_workspace(plan),
        competitor_workspace(plan) {}

  void reset(const ParamView &params,
             const int first_param_row,
             const ExactTriggerState &trigger_state,
             const ExactSequenceState &sequence_state,
             const double observed_time) {
    oracle.reset(
        params, first_param_row, trigger_state, sequence_state, observed_time);
    target_evaluator.reset(&oracle, RelationView{});
    competitor_evaluator.reset(&oracle, RelationView{});
    target_workspace.reset(&oracle, RelationView{});
    competitor_workspace.reset(&oracle, RelationView{});
    condition_frames.clear();
    next_condition_frame_id = 1;
  }

  ExactSourceOracle oracle;
  ForcedExprEvaluator target_evaluator;
  ForcedExprEvaluator competitor_evaluator;
  ForcedExprWorkspace target_workspace;
  ForcedExprWorkspace competitor_workspace;
  ExactRuntimePreparedCompetitorCache base_competitor_cache;
  ExactRuntimePreparedCompetitorCache conditioned_competitor_cache;
  std::vector<double> term_values;
  std::vector<double> group_tails;
  std::vector<double> group_non_wins;
  std::deque<ExactConditionFrame> condition_frames;
  semantic::Index next_condition_frame_id{1};
};

struct ExactEvaluatorConditionGuard {
  ExactEvaluatorConditionGuard(ForcedExprEvaluator *evaluator,
                               const ExactConditionFrame *frame)
      : evaluator_(evaluator) {
    if (evaluator_ != nullptr) {
      previous_ = evaluator_->set_condition_frame(frame);
    }
  }

  ExactEvaluatorConditionGuard(const ExactEvaluatorConditionGuard &) = delete;
  ExactEvaluatorConditionGuard &operator=(
      const ExactEvaluatorConditionGuard &) = delete;

  ~ExactEvaluatorConditionGuard() {
    if (evaluator_ != nullptr) {
      evaluator_->set_condition_frame(previous_);
    }
  }

private:
  ForcedExprEvaluator *evaluator_{nullptr};
  const ExactConditionFrame *previous_{nullptr};
};

inline const ExactConditionFrame *store_condition_frame(
    ExactStepWorkspace *workspace,
    ExactConditionFrame frame) {
  frame.id = workspace->next_condition_frame_id++;
  workspace->condition_frames.push_back(std::move(frame));
  return &workspace->condition_frames.back();
}

inline ExactConditionFrame copy_condition_frame(
    const ExactConditionFrame *parent) {
  ExactConditionFrame frame;
  frame.parent = parent;
  frame.impossible = parent != nullptr && parent->impossible;
  frame.has_source_order_facts =
      parent != nullptr && parent->has_source_order_facts;
  frame.has_guard_upper_bounds =
      parent != nullptr && parent->has_guard_upper_bounds;
  frame.has_expr_upper_bounds =
      parent != nullptr && parent->has_expr_upper_bounds;
  if (parent != nullptr) {
    frame.source_exact_times = parent->source_exact_times;
    frame.source_lower_bounds = parent->source_lower_bounds;
    frame.source_upper_bounds = parent->source_upper_bounds;
  }
  return frame;
}

inline void ensure_source_condition_size(std::vector<double> *values,
                                         const semantic::Index source_id) {
  if (source_id == semantic::kInvalidIndex) {
    return;
  }
  const auto size = static_cast<std::size_t>(source_id + 1);
  if (values->size() < size) {
    values->resize(size, std::numeric_limits<double>::quiet_NaN());
  }
}

inline const ExactConditionFrame *append_exact_source_frame(
    ExactStepWorkspace *workspace,
    const ExactConditionFrame *parent,
    const semantic::Index source_id,
    const double time) {
  if (source_id == semantic::kInvalidIndex || !std::isfinite(time)) {
    return parent;
  }
  auto frame = copy_condition_frame(parent);
  ensure_source_condition_size(&frame.source_exact_times, source_id);
  frame.source_exact_times[static_cast<std::size_t>(source_id)] = time;
  frame.exact_times.push_back(ExactTimedSourceFact{source_id, time});
  return store_condition_frame(workspace, std::move(frame));
}

inline const ExactConditionFrame *append_runtime_condition_frame(
    ExactStepWorkspace *workspace,
    const ExactConditionFrame *parent,
    const ExactRuntimeTermCondition &condition,
    const double time,
    ForcedExprEvaluator *normalizer_evaluator) {
  if (runtime_condition_empty(condition)) {
    return parent;
  }
  auto frame = copy_condition_frame(parent);
  if (runtime_condition_order_contradiction(condition)) {
    frame.impossible = true;
    return store_condition_frame(workspace, std::move(frame));
  }
  if (condition.exact_source_id != semantic::kInvalidIndex) {
    ensure_source_condition_size(
        &frame.source_exact_times, condition.exact_source_id);
    frame.source_exact_times[
        static_cast<std::size_t>(condition.exact_source_id)] = time;
    frame.exact_times.push_back(
        ExactTimedSourceFact{condition.exact_source_id, time});
  }
  for (const auto source_id : condition.upper_bound_source_ids) {
    ensure_source_condition_size(&frame.source_upper_bounds, source_id);
    auto &upper =
        frame.source_upper_bounds[static_cast<std::size_t>(source_id)];
    upper = std::isfinite(upper) ? std::min(upper, time) : time;
    frame.upper_bounds.push_back(ExactTimedSourceFact{source_id, time});
  }
  for (const auto source_id : condition.lower_bound_source_ids) {
    ensure_source_condition_size(&frame.source_lower_bounds, source_id);
    auto &lower =
        frame.source_lower_bounds[static_cast<std::size_t>(source_id)];
    lower = std::isfinite(lower) ? std::max(lower, time) : time;
    frame.lower_bounds.push_back(ExactTimedSourceFact{source_id, time});
  }
  for (const auto &fact : condition.source_order_facts) {
    append_runtime_source_order_fact(
        &frame.source_order_facts,
        fact.before_source_id,
        fact.after_source_id);
  }
  frame.has_source_order_facts =
      frame.has_source_order_facts || !frame.source_order_facts.empty();

  const ExactConditionFrame *current =
      store_condition_frame(workspace, std::move(frame));
  for (const auto expr_id : condition.upper_bound_expr_ids) {
    if (normalizer_evaluator == nullptr) {
      auto impossible = copy_condition_frame(current);
      impossible.impossible = true;
      return store_condition_frame(workspace, std::move(impossible));
    }
    const ExactEvaluatorConditionGuard evaluator_frame(
        normalizer_evaluator, current);
    const double normalizer = normalizer_evaluator->expr_cdf(expr_id);
    if (!(normalizer > 0.0)) {
      auto impossible = copy_condition_frame(current);
      impossible.impossible = true;
      return store_condition_frame(workspace, std::move(impossible));
    }
    auto next = copy_condition_frame(current);
    next.has_expr_upper_bounds = true;
    next.expr_upper_bounds.push_back(
        ExactTimedExprUpperBound{expr_id, time, normalizer});
    current = store_condition_frame(workspace, std::move(next));
  }

  for (const auto &fact : condition.guard_upper_bound_facts) {
    if (normalizer_evaluator == nullptr) {
      auto impossible = copy_condition_frame(current);
      impossible.impossible = true;
      return store_condition_frame(workspace, std::move(impossible));
    }
    const ExactEvaluatorConditionGuard evaluator_frame(
        normalizer_evaluator, current);
    const double normalizer = normalizer_evaluator->expr_cdf(fact.expr_id);
    if (!(normalizer > 0.0)) {
      auto impossible = copy_condition_frame(current);
      impossible.impossible = true;
      return store_condition_frame(workspace, std::move(impossible));
    }
    auto next = copy_condition_frame(current);
    next.has_guard_upper_bounds = true;
    next.guard_upper_bounds.push_back(
        ExactTimedGuardUpperBound{
            fact.ref_source_id,
            fact.blocker_source_id,
            time,
            normalizer});
    current = store_condition_frame(workspace, std::move(next));
  }
  return current;
}

struct ExactRuntimeConditionGuards {
  bool impossible{false};
  std::unique_ptr<ExactSourceOracle::ExactTimeOverlayGuard> exact_time;
  std::vector<std::unique_ptr<ExactSourceOracle::SourceUpperBoundOverlayGuard>>
      upper_bounds;
  std::vector<std::unique_ptr<ExactSourceOracle::SourceLowerBoundOverlayGuard>>
      lower_bounds;
  std::vector<std::unique_ptr<ExactSourceOracle::ExprUpperBoundOverlayGuard>>
      expr_upper_bounds;
  std::vector<std::unique_ptr<ExactSourceOracle::GuardUpperBoundOverlayGuard>>
      guard_upper_bounds;
  std::unique_ptr<ExactSourceOracle::SourceOrderOverlayGuard> source_order;
};

inline ExactRuntimeConditionGuards apply_runtime_condition(
    ExactSourceOracle *oracle,
    const ExactRuntimeTermCondition &condition,
    const double time,
    ForcedExprEvaluator *expr_evaluator = nullptr) {
  ExactRuntimeConditionGuards guards;
  if (oracle == nullptr || runtime_condition_empty(condition)) {
    return guards;
  }
  if (runtime_condition_order_contradiction(condition)) {
    guards.impossible = true;
    return guards;
  }
  if (condition.exact_source_id != semantic::kInvalidIndex) {
    guards.exact_time =
        std::make_unique<ExactSourceOracle::ExactTimeOverlayGuard>(
            oracle, condition.exact_source_id, time);
  }
  guards.upper_bounds.reserve(condition.upper_bound_source_ids.size());
  for (const auto source_id : condition.upper_bound_source_ids) {
    guards.upper_bounds.push_back(
        std::make_unique<ExactSourceOracle::SourceUpperBoundOverlayGuard>(
            oracle, source_id, time));
  }
  guards.lower_bounds.reserve(condition.lower_bound_source_ids.size());
  for (const auto source_id : condition.lower_bound_source_ids) {
    guards.lower_bounds.push_back(
        std::make_unique<ExactSourceOracle::SourceLowerBoundOverlayGuard>(
            oracle, source_id, time));
  }
  if (!condition.source_order_facts.empty()) {
    guards.source_order =
        std::make_unique<ExactSourceOracle::SourceOrderOverlayGuard>(
            oracle, condition.source_order_facts);
  }
  guards.expr_upper_bounds.reserve(condition.upper_bound_expr_ids.size());
  for (const auto expr_id : condition.upper_bound_expr_ids) {
    const double normalizer =
        expr_evaluator != nullptr ? expr_evaluator->expr_cdf(expr_id) : 0.0;
    if (!(normalizer > 0.0)) {
      guards.impossible = true;
      return guards;
    }
    guards.expr_upper_bounds.push_back(
        std::make_unique<ExactSourceOracle::ExprUpperBoundOverlayGuard>(
            oracle, expr_id, time, normalizer));
  }
  guards.guard_upper_bounds.reserve(condition.guard_upper_bound_facts.size());
  for (const auto &fact : condition.guard_upper_bound_facts) {
    const double normalizer =
        expr_evaluator != nullptr ? expr_evaluator->expr_cdf(fact.expr_id) : 0.0;
    if (!(normalizer > 0.0)) {
      guards.impossible = true;
      return guards;
    }
    guards.guard_upper_bounds.push_back(
        std::make_unique<ExactSourceOracle::GuardUpperBoundOverlayGuard>(
            oracle,
            fact.ref_source_id,
            fact.blocker_source_id,
            time,
            normalizer));
  }
  return guards;
}

inline double competitor_subset_win_mass(
    const ExactCompetitorSubsetPlan &subset_plan,
    const semantic::Index target_active_source_id,
    const double readiness_upper,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t) {
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  double win_prob = 0.0;
  if (subset_plan.outcome_indices.size() == 1U &&
      subset_plan.singleton_expr_root != semantic::kInvalidIndex) {
    win_prob = evaluator->expr_cdf(subset_plan.singleton_expr_root);
  } else {
    if (evaluator->oracle()->has_guard_upper_bound_overlay(
            evaluator->condition_frame())) {
      for (const auto outcome_idx : subset_plan.outcome_indices) {
        const auto expr_root =
            evaluator->plan().outcomes[static_cast<std::size_t>(outcome_idx)]
                .expr_root;
        const double marginal = evaluator->expr_cdf(expr_root);
        if (!(marginal > 0.0)) {
          return 0.0;
        }
      }
    }
    for (const auto &scenario : subset_plan.scenarios) {
      const auto scenario_view =
          make_exact_scenario_runtime_view(evaluator->plan(), scenario);
      win_prob += scenario_truth_cdf(scenario_view, evaluator, t, workspace);
      if (!std::isfinite(win_prob)) {
        return 1.0;
      }
    }
    win_prob = clamp_probability(win_prob);
  }
  double same_active_t = 0.0;
  double same_active_ready = 0.0;
  for (const auto &scenario : subset_plan.scenarios) {
    if (scenario.active_source_id != target_active_source_id) {
      continue;
    }
    const auto scenario_view =
        make_exact_scenario_runtime_view(evaluator->plan(), scenario);
    same_active_t +=
        same_active_win_mass(scenario_view, evaluator, t, t, workspace);
    same_active_ready +=
        same_active_win_mass(
            scenario_view, evaluator, t, readiness_upper, workspace);
  }
  return clamp_probability(win_prob - same_active_t + same_active_ready);
}

inline double competitor_non_win_probability(
    const ExactTargetCompetitorPlan &competitor_plan,
    const semantic::Index target_active_source_id,
    const double readiness_upper,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    const std::vector<std::uint8_t> *used_outcomes = nullptr) {
  double non_win = 1.0;
  for (const auto &block : competitor_plan.blocks) {
    double union_prob = 0.0;
    for (const auto &subset : block.subsets) {
      if (used_outcomes != nullptr) {
        bool skip_subset = false;
        for (const auto outcome_idx : subset.outcome_indices) {
          if ((*used_outcomes)[static_cast<std::size_t>(outcome_idx)] != 0U) {
            skip_subset = true;
            break;
          }
        }
        if (skip_subset) {
          continue;
        }
      }
      union_prob += static_cast<double>(subset.inclusion_sign) *
                    competitor_subset_win_mass(
                        subset,
                        target_active_source_id,
                        readiness_upper,
                        evaluator,
                        workspace,
                        t);
    }
    union_prob = clamp_probability(union_prob);
    non_win *= clamp_probability(1.0 - union_prob);
    if (!(non_win > 0.0)) {
      return 0.0;
    }
  }
  return clamp_probability(non_win);
}

inline double runtime_competitor_subset_win_at_t(
    const ExactRuntimeCompetitorSubsetPlan &subset,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t) {
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  if (subset.outcome_indices.size() == 1U &&
      subset.singleton_expr_root != semantic::kInvalidIndex) {
    return evaluator->expr_cdf(subset.singleton_expr_root);
  }

  if (evaluator->oracle()->has_guard_upper_bound_overlay(
          evaluator->condition_frame())) {
    for (const auto outcome_idx : subset.outcome_indices) {
      const auto expr_root =
          evaluator->plan().outcomes[static_cast<std::size_t>(outcome_idx)]
              .expr_root;
      const double marginal = evaluator->expr_cdf(expr_root);
      if (!(marginal > 0.0)) {
        return 0.0;
      }
    }
  }
  double win_prob = 0.0;
  for (const auto &scenario : subset.scenarios) {
    win_prob += runtime_scenario_truth_cdf(
        scenario,
        evaluator,
        t,
        workspace);
    if (!std::isfinite(win_prob)) {
      return 1.0;
    }
  }
  return clamp_probability(win_prob);
}

inline double runtime_same_active_mass(
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const ExactRuntimeScenarioSubsetView &subset_view,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    const double ready_upper) {
  double total = 0.0;
  for (const auto scenario_idx : subset_view.same_active_scenario_indices) {
    total += runtime_same_active_win_mass(
        subset.scenarios[static_cast<std::size_t>(scenario_idx)],
        evaluator,
        t,
        ready_upper,
        workspace);
  }
  return total;
}

inline void prepare_runtime_competitor_non_win(
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactRuntimeScenarioCompetitorView &scenario_competitor,
    const std::vector<std::uint8_t> *used_outcomes,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    ExactRuntimePreparedCompetitorCache *cache,
    const ExactRuntimePreparedCompetitorCache *base_cache = nullptr,
    const ExactRuntimeCompetitorSubsetMask *mask_a = nullptr,
    const ExactRuntimeCompetitorSubsetMask *mask_b = nullptr) {
  auto &offsets = cache->block_offsets;
  auto &prepared = cache->subsets;
  cache->clear();
  offsets.reserve(runtime_outcome.competitor_blocks.size() + 1U);
  offsets.push_back(0);

  std::size_t base_pos = 0;
  for (std::size_t block_idx = 0;
       block_idx < runtime_outcome.competitor_blocks.size();
       ++block_idx) {
    const auto &block = runtime_outcome.competitor_blocks[block_idx];
    const auto &block_view = scenario_competitor.blocks[block_idx];
    for (std::size_t subset_idx = 0; subset_idx < block.subsets.size();
         ++subset_idx) {
      const auto &subset = block.subsets[subset_idx];
      if (runtime_subset_is_used(subset, used_outcomes)) {
        continue;
      }
      const bool has_subset_mask = mask_a != nullptr || mask_b != nullptr;
      const bool affected =
          !has_subset_mask ||
          runtime_competitor_subset_mask_affected(mask_a, block_idx, subset_idx) ||
          runtime_competitor_subset_mask_affected(mask_b, block_idx, subset_idx);
      if (base_cache != nullptr && !affected &&
          base_pos < base_cache->subsets.size() &&
          base_cache->subsets[base_pos].block_index ==
              static_cast<semantic::Index>(block_idx) &&
          base_cache->subsets[base_pos].subset_index ==
              static_cast<semantic::Index>(subset_idx)) {
        prepared.push_back(base_cache->subsets[base_pos]);
        ++base_pos;
        continue;
      }
      if (base_cache != nullptr) {
        ++base_pos;
      }
      const auto &subset_view = block_view.subsets[subset_idx];
      const double win_prob =
          runtime_competitor_subset_win_at_t(subset, evaluator, workspace, t);
      const double same_active_t =
          runtime_same_active_mass(subset, subset_view, evaluator, workspace, t, t);
      prepared.push_back(ExactRuntimePreparedCompetitorSubset{
          static_cast<semantic::Index>(block_idx),
          static_cast<semantic::Index>(subset_idx),
          win_prob - same_active_t});
    }
    offsets.push_back(static_cast<semantic::Index>(prepared.size()));
  }
}

inline double runtime_competitor_non_win_probability(
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactRuntimeScenarioCompetitorView &scenario_competitor,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    const double readiness_upper,
    const ExactRuntimePreparedCompetitorCache &cache) {
  double non_win = 1.0;
  const auto &offsets = cache.block_offsets;
  const auto &prepared = cache.subsets;
  for (std::size_t block_idx = 0;
       block_idx < runtime_outcome.competitor_blocks.size();
       ++block_idx) {
    double union_prob = 0.0;
    const auto begin = static_cast<std::size_t>(offsets[block_idx]);
    const auto end = static_cast<std::size_t>(offsets[block_idx + 1U]);
    for (std::size_t i = begin; i < end; ++i) {
      const auto &prepared_subset = prepared[i];
      const auto subset_idx = static_cast<std::size_t>(
          prepared_subset.subset_index);
      const auto &subset =
          runtime_outcome.competitor_blocks[block_idx].subsets[subset_idx];
      const auto &subset_view =
          scenario_competitor.blocks[block_idx].subsets[subset_idx];
      const double same_active_ready = runtime_same_active_mass(
          subset,
          subset_view,
          evaluator,
          workspace,
          t,
          readiness_upper);
      const double subset_mass =
          clamp_probability(prepared_subset.fixed_mass + same_active_ready);
      union_prob += static_cast<double>(subset.inclusion_sign) * subset_mass;
    }
    union_prob = clamp_probability(union_prob);
    non_win *= clamp_probability(1.0 - union_prob);
    if (!(non_win > 0.0)) {
      return 0.0;
    }
  }
  return clamp_probability(non_win);
}

inline ExactSequenceState advance_exact_sequence_state(
    const ExactVariantPlan &plan,
    const ExactSequenceState &sequence_state,
    const ExactTransitionScenario &scenario,
    const double observed_time) {
  ExactSequenceState next_state = sequence_state;
  next_state.lower_bound = observed_time;
  if (next_state.exact_times.empty()) {
    next_state.exact_times.assign(
        static_cast<std::size_t>(plan.source_count),
        std::numeric_limits<double>::quiet_NaN());
  }
  next_state.exact_times[static_cast<std::size_t>(scenario.active_source_id)] =
      observed_time;
  return next_state;
}

inline double evaluate_scenario_probability(
    const ExactTargetCompetitorPlan &competitor_plan,
    const ExactScenarioRuntimeView &scenario_view,
    const double observed_time,
    const std::vector<std::uint8_t> *used_outcomes,
    ForcedExprEvaluator *target_evaluator,
    ForcedExprWorkspace *target_workspace,
    ForcedExprEvaluator *competitor_evaluator,
    ForcedExprWorkspace *competitor_workspace) {
  const double active_pdf = source_pdf_at(
      target_evaluator->oracle(),
      scenario_view.scenario.active_source_id,
      observed_time);
  const double tail = after_survival(
      scenario_view, target_evaluator, observed_time, target_workspace);
  if (!(active_pdf > 0.0) || !(tail > 0.0)) {
    return 0.0;
  }

  const auto competitor_non_win = [&](const double readiness_upper) {
    return competitor_non_win_probability(
        competitor_plan,
        scenario_view.scenario.active_source_id,
        readiness_upper,
        competitor_evaluator,
        competitor_workspace,
        observed_time,
        used_outcomes);
  };

  if (!scenario_view.has_readiness()) {
    return active_pdf * tail * competitor_non_win(0.0);
  }

  double total = 0.0;
  const double initial_ready =
      readiness_cdf(scenario_view, target_evaluator, 0.0, target_workspace);
  if (initial_ready > 0.0) {
    total += initial_ready * active_pdf * tail * competitor_non_win(0.0);
  }

  const auto batch = quadrature::build_finite_batch(0.0, observed_time);
  for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
    const double readiness_time = batch.nodes.nodes[i];
    double value = readiness_density(
        scenario_view, target_evaluator, readiness_time, target_workspace);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    value *= active_pdf * tail;
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    value *= competitor_non_win(readiness_time);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    total += batch.nodes.weights[i] * value;
  }
  return total;
}

inline double evaluate_runtime_scenario_probability(
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactRuntimeScenarioFormula &scenario_formula,
    const ExactRuntimeScenarioCompetitorView &scenario_competitor,
    const double observed_time,
    const std::vector<std::uint8_t> *used_outcomes,
    ForcedExprEvaluator *target_evaluator,
    ForcedExprWorkspace *target_workspace,
    ForcedExprEvaluator *competitor_evaluator,
    ForcedExprWorkspace *competitor_workspace,
    ExactStepWorkspace *step_workspace) {
  const double active_pdf = source_pdf_at(
      target_evaluator->oracle(),
      scenario_formula.active_source_id,
      observed_time);
  if (!(active_pdf > 0.0)) {
    return 0.0;
  }

  const ExactConditionFrame *scenario_frame = append_exact_source_frame(
      step_workspace,
      nullptr,
      scenario_formula.active_source_id,
      observed_time);
  const ExactEvaluatorConditionGuard target_scenario_frame(
      target_evaluator, scenario_frame);
  const ExactEvaluatorConditionGuard competitor_scenario_frame(
      competitor_evaluator, scenario_frame);
  std::optional<ExactRuntimeTermCondition> ranked_tail_competitor_condition;
  const ExactRuntimeTermCondition *tail_competitor_condition =
      &scenario_formula.tail_competitor_condition;
  if (used_outcomes != nullptr) {
    ranked_tail_competitor_condition =
        filter_runtime_condition_for_competitors(
            scenario_formula.tail_condition,
            runtime_outcome,
            used_outcomes,
            target_evaluator->plan());
    tail_competitor_condition = &*ranked_tail_competitor_condition;
  }
  const bool has_tail_competitor_condition =
      !runtime_condition_empty(*tail_competitor_condition);

  auto prepare_competitor = [&]() {
    prepare_runtime_competitor_non_win(
        runtime_outcome,
        scenario_competitor,
        used_outcomes,
        competitor_evaluator,
        competitor_workspace,
        observed_time,
        &step_workspace->base_competitor_cache);
  };

  auto competitor_non_win_prepared =
      [&](const double readiness_upper,
          const ExactRuntimePreparedCompetitorCache &cache) {
    return runtime_competitor_non_win_probability(
        runtime_outcome,
        scenario_competitor,
        competitor_evaluator,
        competitor_workspace,
        observed_time,
        readiness_upper,
        cache);
  };

  bool base_competitor_prepared = false;
  const auto competitor_non_win = [&](const double readiness_upper) {
    if (!base_competitor_prepared) {
      prepare_competitor();
      base_competitor_prepared = true;
    }
    return competitor_non_win_prepared(
        readiness_upper,
        step_workspace->base_competitor_cache);
  };

  const auto conditioned_competitor_non_win =
      [&](const double readiness_upper,
          const ExactRuntimeCompetitorSubsetMask *mask_a,
          const ExactRuntimeCompetitorSubsetMask *mask_b,
          const ExactConditionFrame *condition_frame) {
    if (!base_competitor_prepared) {
      prepare_competitor();
      base_competitor_prepared = true;
    }
    const ExactEvaluatorConditionGuard competitor_frame(
        competitor_evaluator, condition_frame);
    prepare_runtime_competitor_non_win(
        runtime_outcome,
        scenario_competitor,
        used_outcomes,
        competitor_evaluator,
        competitor_workspace,
        observed_time,
        &step_workspace->conditioned_competitor_cache,
        &step_workspace->base_competitor_cache,
        mask_a,
        mask_b);
    return competitor_non_win_prepared(
        readiness_upper,
        step_workspace->conditioned_competitor_cache);
  };

  const auto competitor_non_win_with_conditions =
      [&](const double readiness_upper,
          const ExactRuntimeTermCondition &term_condition,
          const double term_time,
          ForcedExprEvaluator *term_evaluator,
          const ExactRuntimeCompetitorSubsetMask *term_mask) {
        const bool has_term_condition =
            !runtime_condition_empty(term_condition);
        if (!has_term_condition && !has_tail_competitor_condition) {
          return competitor_non_win(readiness_upper);
        }
        const auto *tail_mask =
            has_tail_competitor_condition && used_outcomes == nullptr
                ? &scenario_formula.tail_competitor_subset_mask
                : nullptr;
        const auto *prepared_term_mask =
            has_term_condition && used_outcomes == nullptr ? term_mask : nullptr;
        const ExactConditionFrame *condition_frame = scenario_frame;
        if (has_term_condition) {
          condition_frame = append_runtime_condition_frame(
              step_workspace,
              condition_frame,
              term_condition,
              term_time,
              term_evaluator);
          if (condition_frame->impossible) {
            return 0.0;
          }
        }

        if (has_tail_competitor_condition) {
          condition_frame = append_runtime_condition_frame(
              step_workspace,
              condition_frame,
              *tail_competitor_condition,
              observed_time,
              target_evaluator);
          if (condition_frame->impossible) {
            return 0.0;
          }
        }
        return conditioned_competitor_non_win(
            readiness_upper,
            prepared_term_mask,
            tail_mask,
            condition_frame);
      };

  const auto tail_with_condition =
      [&](const ExactRuntimeTermCondition &tail_condition,
          const double condition_time,
          ForcedExprEvaluator *condition_evaluator) {
        if (runtime_condition_empty(tail_condition)) {
          return runtime_after_survival(
              scenario_formula,
              target_evaluator,
              observed_time,
              target_workspace);
        }
        const ExactConditionFrame *condition_frame =
            append_runtime_condition_frame(
                step_workspace,
                scenario_frame,
                tail_condition,
                condition_time,
                condition_evaluator);
        if (condition_frame->impossible) {
          return 0.0;
        }
        const ExactEvaluatorConditionGuard target_frame(
            target_evaluator, condition_frame);
        return runtime_after_survival(
            scenario_formula,
            target_evaluator,
            observed_time,
            target_workspace);
      };

  if (!scenario_formula.has_readiness) {
    const double tail = runtime_after_survival(
        scenario_formula, target_evaluator, observed_time, target_workspace);
    if (!(tail > 0.0)) {
      return 0.0;
    }
    const double non_win = competitor_non_win_with_conditions(
        0.0,
        ExactRuntimeTermCondition{},
        observed_time,
        target_evaluator,
        nullptr);
    return active_pdf * tail * non_win;
  }

  double total = 0.0;
  const double base_tail = runtime_after_survival(
      scenario_formula, target_evaluator, observed_time, target_workspace);
  const double initial_ready = runtime_readiness_cdf(
      scenario_formula, target_evaluator, 0.0, target_workspace);
  if (initial_ready > 0.0 && base_tail > 0.0) {
    const double non_win = competitor_non_win_with_conditions(
        0.0,
        ExactRuntimeTermCondition{},
        observed_time,
        target_evaluator,
        nullptr);
    total += initial_ready * active_pdf * base_tail * non_win;
  }

  const auto &readiness_terms = scenario_formula.readiness_density.sum_terms;
  std::vector<ExactRuntimeTermCondition> ranked_term_competitor_conditions;
  if (used_outcomes != nullptr) {
    ranked_term_competitor_conditions.reserve(readiness_terms.size());
    for (const auto &term : readiness_terms) {
      ranked_term_competitor_conditions.push_back(
          filter_runtime_condition_for_competitors(
              term.condition,
              runtime_outcome,
              used_outcomes,
              target_evaluator->plan()));
    }
  }
  const auto term_competitor_condition =
      [&](const std::size_t term_idx,
          const ExactRuntimeProductTerm &term) -> const ExactRuntimeTermCondition & {
        return used_outcomes != nullptr
                   ? ranked_term_competitor_conditions[term_idx]
                   : term.competitor_condition;
      };

  if (!scenario_formula.has_conditioned_readiness_terms &&
      !has_tail_competitor_condition) {
    if (!(base_tail > 0.0)) {
      return total;
    }
    const auto batch = quadrature::build_finite_batch(0.0, observed_time);
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double readiness_time = batch.nodes.nodes[i];
      double value = runtime_readiness_density(
          scenario_formula, target_evaluator, readiness_time, target_workspace);
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      value *= active_pdf * base_tail;
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      value *= competitor_non_win(readiness_time);
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      total += batch.nodes.weights[i] * value;
    }
    return total;
  }

  const auto batch = quadrature::build_finite_batch(0.0, observed_time);
  for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
    const double readiness_time = batch.nodes.nodes[i];
    const auto time_guard =
        target_evaluator->oracle()->conditional_time_guard(readiness_time);
    ForcedExprEvaluator *scenario_evaluator = nullptr;
    if (scenario_formula.readiness_density.requires_scenario) {
      scenario_evaluator = prepare_runtime_scenario_evaluator(
          scenario_formula,
          target_evaluator,
          target_workspace);
      if (scenario_evaluator == nullptr) {
        continue;
      }
    }
    const auto condition_evaluator =
        scenario_evaluator != nullptr ? scenario_evaluator : target_evaluator;
    auto &term_values = step_workspace->term_values;
    term_values.assign(readiness_terms.size(), 0.0);
    bool has_nonzero_term = false;
    for (std::size_t term_idx = 0; term_idx < readiness_terms.size();
         ++term_idx) {
      const auto &term = readiness_terms[term_idx];
      if (term.condition_impossible) {
        continue;
      }
      const double value = evaluate_runtime_factors(
          term.factors,
          target_evaluator,
          scenario_evaluator);
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      term_values[term_idx] = value;
      has_nonzero_term = true;
    }
    if (!has_nonzero_term) {
      continue;
    }

    if (used_outcomes == nullptr) {
      auto &group_tails = step_workspace->group_tails;
      auto &group_non_wins = step_workspace->group_non_wins;
      group_tails.assign(
          scenario_formula.condition_groups.size(),
          std::numeric_limits<double>::quiet_NaN());
      group_non_wins.assign(
          scenario_formula.condition_groups.size(),
          std::numeric_limits<double>::quiet_NaN());
      const auto ensure_group = [&](const std::size_t group_idx) {
        if (std::isfinite(group_tails[group_idx])) {
          return;
        }
        const auto &group = scenario_formula.condition_groups[group_idx];
        const double group_tail =
            runtime_condition_empty(group.tail_condition)
                ? base_tail
                : tail_with_condition(
                      group.tail_condition,
                      readiness_time,
                      condition_evaluator);
        group_tails[group_idx] = group_tail;
        group_non_wins[group_idx] =
            group_tail > 0.0
                ? competitor_non_win_with_conditions(
                      readiness_time,
                      group.competitor_condition,
                      readiness_time,
                      condition_evaluator,
                      &group.competitor_subset_mask)
                : 0.0;
      };
      for (std::size_t term_idx = 0; term_idx < readiness_terms.size();
           ++term_idx) {
        const auto &term = readiness_terms[term_idx];
        double value = term_values[term_idx];
        if (value == 0.0 ||
            term.context_group_index == semantic::kInvalidIndex) {
          continue;
        }
        const auto group_idx =
            static_cast<std::size_t>(term.context_group_index);
        ensure_group(group_idx);
        const double term_tail = group_tails[group_idx];
        const double term_non_win = group_non_wins[group_idx];
        if (!(term_tail > 0.0)) {
          continue;
        }
        value *= active_pdf * term_tail;
        if (!std::isfinite(value) || value == 0.0) {
          continue;
        }
        value *= term_non_win;
        if (!std::isfinite(value) || value == 0.0) {
          continue;
        }
        total += batch.nodes.weights[i] * value;
      }
      continue;
    }

    for (std::size_t term_idx = 0; term_idx < readiness_terms.size();
         ++term_idx) {
      const auto &term = readiness_terms[term_idx];
      double value = term_values[term_idx];
      if (value == 0.0) {
        continue;
      }

      double term_tail = 0.0;
      double term_non_win = 0.0;
      term_tail = runtime_condition_empty(term.tail_condition)
                      ? base_tail
                      : tail_with_condition(
                            term.tail_condition,
                            readiness_time,
                            condition_evaluator);
      term_non_win = competitor_non_win_with_conditions(
          readiness_time,
          term_competitor_condition(term_idx, term),
          readiness_time,
          condition_evaluator,
          nullptr);
      if (!(term_tail > 0.0)) {
        continue;
      }
      value *= active_pdf * term_tail;
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      value *= term_non_win;
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      total += batch.nodes.weights[i] * value;
    }
  }
  return total;
}

inline ExactStepResult evaluate_exact_step(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const ExactSequenceState &sequence_state,
    const semantic::Index target_idx,
    const double observed_time,
    const std::vector<std::uint8_t> *used_outcomes = nullptr,
    const bool collect_successors = false,
    ExactStepWorkspace *workspace = nullptr) {
  ExactStepResult result;
  const auto target_pos = static_cast<std::size_t>(target_idx);
  std::optional<ExactStepWorkspace> local_workspace;
  if (workspace == nullptr) {
    local_workspace.emplace(plan);
    workspace = &*local_workspace;
  }
  auto &step_workspace = *workspace;
  step_workspace.reset(
      params, first_param_row, trigger_state, sequence_state, observed_time);
  auto &oracle = step_workspace.oracle;
  auto &target_evaluator = step_workspace.target_evaluator;
  auto &competitor_evaluator = step_workspace.competitor_evaluator;
  auto &target_workspace = step_workspace.target_workspace;
  auto &competitor_workspace = step_workspace.competitor_workspace;
  const auto &runtime_outcome = plan.runtime.outcomes[target_pos];

  for (std::size_t scenario_idx = 0;
       scenario_idx < plan.outcomes[target_pos].scenarios.size();
       ++scenario_idx) {
    const auto &scenario =
        plan.outcomes[target_pos].scenarios[scenario_idx];
    const auto &runtime_scenario =
        runtime_outcome.scenarios[scenario_idx];
    if (collect_successors && !scenario_supports_ranked_sequence(scenario)) {
      throw std::logic_error(
          "ranked exact step requested for unsupported latent-readiness scenario");
    }

    RelationView competitor_view;
    if (!relation_view_with_overlay(
            RelationView{}, runtime_scenario.relation_template, &competitor_view)) {
      continue;
    }
    competitor_evaluator.reset(&oracle, competitor_view);

    const double scenario_prob = evaluate_runtime_scenario_probability(
        runtime_outcome,
        runtime_scenario,
        runtime_outcome.competitor_by_scenario[scenario_idx],
        observed_time,
        used_outcomes,
        &target_evaluator,
        &target_workspace,
        &competitor_evaluator,
        &competitor_workspace,
        &step_workspace);
    if (!(scenario_prob > 0.0)) {
      continue;
    }
    result.total_probability += scenario_prob;
    if (collect_successors) {
      result.branches.push_back(ExactStepBranch{
          scenario_prob,
          advance_exact_sequence_state(
              plan, sequence_state, scenario, observed_time)});
    }
  }

  return result;
}

} // namespace detail
} // namespace accumulatr::eval
