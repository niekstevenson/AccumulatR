#pragma once

#include <limits>
#include <optional>
#include <utility>

#include "eval_query.hpp"
#include "exact_step_distribution.hpp"
#include "quadrature.hpp"
#include "trial_data.hpp"
#include "../compile/exact_evaluation_program_lowering.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactTrialColumns {
  std::vector<const int *> labels;
  std::vector<const double *> times;
};

struct ExactTrialView {
  semantic::Index variant_index{semantic::kInvalidIndex};
  R_xlen_t row{0};
  int rank_count{0};
  const ExactTrialColumns *columns{nullptr};
};

inline ExactTrialColumns make_exact_trial_columns(
    SEXP dataSEXP,
    const PreparedTrialLayout &layout) {
  ExactTrialColumns columns;
  const int max_rank = layout.max_rank;
  columns.labels.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  columns.times.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  for (int rank = 1; rank <= max_rank; ++rank) {
    columns.labels[static_cast<std::size_t>(rank)] =
        INTEGER(trusted_data_column(
            dataSEXP,
            layout.label_cols[static_cast<std::size_t>(rank)]));
    columns.times[static_cast<std::size_t>(rank)] =
        REAL(trusted_data_column(
            dataSEXP,
            layout.time_cols[static_cast<std::size_t>(rank)]));
  }
  return columns;
}

inline semantic::Index exact_trial_view_outcome_code(const ExactTrialView &view,
                                                     const std::size_t rank_idx) {
  return view.columns->labels[rank_idx + 1U][view.row];
}

inline double exact_trial_view_rt(const ExactTrialView &view,
                                  const std::size_t rank_idx) {
  return view.columns->times[rank_idx + 1U][view.row];
}

inline ExactTrialView read_exact_observation_view(
    const PreparedDataView &table,
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const PreparedTrialLayout &layout,
    const std::size_t trial_index,
    const ExactTrialColumns &columns) {
  const auto &trial_row = layout.trials[trial_index];
  const auto row = static_cast<R_xlen_t>(trial_row.start_row);
  const int component_code = table.component[row];
  ExactTrialView obs;
  obs.variant_index =
      variant_index_by_component_code[static_cast<std::size_t>(component_code)];
  obs.row = row;
  obs.columns = &columns;
  for (int rank = 1; rank <= layout.max_rank; ++rank) {
    const auto *rank_labels = columns.labels[static_cast<std::size_t>(rank)];
    const auto *rank_times = columns.times[static_cast<std::size_t>(rank)];
    if (integer_cell_is_na(rank_labels, row) &&
        Rcpp::NumericVector::is_na(rank_times[row])) {
      break;
    }
    ++obs.rank_count;
  }
  return obs;
}

inline void build_exact_plan_cache(
    const compile::CompiledModel &compiled,
    const std::unordered_map<std::string, semantic::Index> &component_code_by_id,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_component_codes,
    const std::size_t n_outcome_codes,
    std::vector<semantic::Index> *variant_index_by_component_code,
    std::vector<ExactVariantPlan> *plans,
    std::vector<ExactComplexityMetrics> *complexity_metrics_by_variant = nullptr) {
  variant_index_by_component_code->assign(
      n_component_codes + 1U,
      semantic::kInvalidIndex);
  plans->clear();
  plans->reserve(compiled.variants.size());
  if (complexity_metrics_by_variant != nullptr) {
    complexity_metrics_by_variant->clear();
    complexity_metrics_by_variant->reserve(compiled.variants.size());
  }
  for (const auto &variant : compiled.variants) {
    const auto plan_index = static_cast<semantic::Index>(plans->size());
    const auto component_it = component_code_by_id.find(variant.component_id);
    if (component_it == component_code_by_id.end()) {
      throw std::runtime_error(
          "exact evaluator found no prepared component code for '" +
          variant.component_id + "'");
    }
    (*variant_index_by_component_code)[static_cast<std::size_t>(component_it->second)] =
        plan_index;
    auto evaluation_program =
        accumulatr::compile::lower_exact_evaluation_program(
            variant,
            outcome_code_by_label);
    ExactComplexityMetrics metrics;
    plans->push_back(
        make_exact_variant_plan(
            std::move(evaluation_program),
            n_outcome_codes,
            complexity_metrics_by_variant == nullptr ? nullptr : &metrics));
    if (complexity_metrics_by_variant != nullptr) {
      complexity_metrics_by_variant->push_back(metrics);
    }
  }
}

inline double exact_unranked_target_density(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const semantic::Index target_idx,
    const double rt,
    ExactStepWorkspace *workspace = nullptr) {
  if (target_idx == semantic::kInvalidIndex ||
      !std::isfinite(rt) ||
      !(rt > 0.0)) {
    return 0.0;
  }
  std::optional<ExactStepWorkspace> local_workspace;
  if (workspace == nullptr) {
    local_workspace.emplace(plan);
    workspace = &*local_workspace;
  }
  double total = 0.0;
  for (const auto &compiled_state : plan.trigger_state_table.states) {
    const auto trigger_state =
        exact_compiled_trigger_state_view(
            plan, params, first_param_row, compiled_state);
    if (!(trigger_state.weight > 0.0)) {
      continue;
    }
    const auto step =
        evaluate_exact_step_distribution(
            plan,
            params,
            first_param_row,
            trigger_state,
            workspace->initial_state,
            target_idx,
            rt,
            nullptr,
            false,
            workspace);
    total += trigger_state.weight * step.total_probability;
  }
  return std::isfinite(total) && total > 0.0 ? total : 0.0;
}

inline double exact_outcome_probability_between(
    const ExactVariantPlan &plan,
    const ParamView &params,
    int first_param_row,
    semantic::Index target_idx,
    double lower,
    double upper,
    ExactStepWorkspace *workspace);

inline double exact_finite_outcome_probability(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const semantic::Index target_idx,
    ExactStepWorkspace *workspace) {
  return exact_outcome_probability_between(
      plan,
      params,
      first_param_row,
      target_idx,
      0.0,
      R_PosInf,
      workspace);
}

inline double exact_outcome_probability_between(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const semantic::Index target_idx,
    double lower,
    const double upper,
    ExactStepWorkspace *workspace) {
  if (target_idx == semantic::kInvalidIndex) {
    return 0.0;
  }
  if (!std::isfinite(lower)) {
    return 0.0;
  }
  lower = std::max(0.0, lower);
  if (std::isfinite(upper) && !(upper > lower)) {
    return 0.0;
  }
  const double probability = std::isfinite(upper)
      ? quadrature::integrate_finite_default(
            lower,
            upper,
            [&](const double rt) {
              return exact_unranked_target_density(
                  plan,
                  params,
                  first_param_row,
                  target_idx,
                  rt,
                  workspace);
            })
      : integrate_to_infinity(
            [&](const double offset) {
              return exact_unranked_target_density(
                  plan,
                  params,
                  first_param_row,
                  target_idx,
                  lower + offset,
                  workspace);
            });
  return std::isfinite(probability) ? clamp_probability(probability) : 0.0;
}

inline double exact_all_outcome_probability_between(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const std::vector<semantic::Index> &outcome_codes,
    const double lower,
    const double upper,
    ExactStepWorkspace *workspace) {
  double total = 0.0;
  for (const auto outcome_code : outcome_codes) {
    if (outcome_code == semantic::kInvalidIndex ||
        static_cast<std::size_t>(outcome_code) >=
            plan.outcome_index_by_code.size()) {
      continue;
    }
    const auto target_idx =
        plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
    total += exact_outcome_probability_between(
        plan,
        params,
        first_param_row,
        target_idx,
        lower,
        upper,
        workspace);
  }
  return std::isfinite(total) ? clamp_probability(total) : 0.0;
}

inline double exact_terminal_no_response_probability(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row) {
  if (!plan.no_response.direct_leaf_failure_product ||
      plan.no_response.leaf_indices.empty()) {
    return 0.0;
  }
  double total = 0.0;
  for (const auto &compiled_state : plan.trigger_state_table.states) {
    const auto trigger_state =
        exact_compiled_trigger_state_view(
            plan, params, first_param_row, compiled_state);
    if (!(trigger_state.weight > 0.0)) {
      continue;
    }
    double product = trigger_state.weight;
    for (const auto leaf_index : plan.no_response.leaf_indices) {
      product *= clamp_probability(
          exact_leaf_q_for_trigger_state(
              plan.program,
              params,
              first_param_row + static_cast<int>(leaf_index),
              trigger_state,
              leaf_index));
      if (!(product > 0.0)) {
        break;
      }
    }
    total += product;
  }
  return std::isfinite(total) ? clamp_probability(total) : 0.0;
}

inline double exact_loglik_for_trial(const ExactVariantPlan &plan,
                                     const ParamView &params,
                                     const int first_param_row,
                                     const semantic::Index outcome_code,
                                     const double rt,
                                     const double min_ll,
                                     ExactStepWorkspace *workspace = nullptr) {
  const auto target_idx =
      plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
  const double total =
      exact_unranked_target_density(
          plan, params, first_param_row, target_idx, rt, workspace);
  return total > 0.0 ? std::log(total) : min_ll;
}

inline void advance_exact_sequence_state(
    ExactSequenceState *state,
    const ExactCompiledTransitionPlan &transition,
    const ExactVariantPlan &plan,
    const double observed_time,
    const std::vector<double> *ready_expr_normalizers = nullptr) {
  if (state == nullptr) {
    return;
  }
  state->lower_bound = observed_time;
  if (transition.release_source_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(transition.release_source_id) <
          state->exact_times.size()) {
    state->exact_times[static_cast<std::size_t>(transition.release_source_id)] =
        observed_time;
  }
  for (const auto source_id : transition.readiness_source_ids) {
    if (source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(source_id) >= state->upper_bounds.size()) {
      continue;
    }
    auto &upper = state->upper_bounds[static_cast<std::size_t>(source_id)];
    upper = std::isfinite(upper) ? std::min(upper, observed_time)
                                 : observed_time;
  }
  for (std::size_t i = 0; i < transition.readiness_expr_ids.size(); ++i) {
    const auto expr_id = transition.readiness_expr_ids[i];
    if (expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(expr_id) >=
            state->expr_upper_bounds.size() ||
        static_cast<std::size_t>(expr_id) >=
            state->expr_upper_normalizers.size()) {
      continue;
    }
    if (ready_expr_normalizers == nullptr ||
        i >= ready_expr_normalizers->size()) {
      continue;
    }
    const double normalizer = (*ready_expr_normalizers)[i];
    if (!(normalizer > 0.0) || !std::isfinite(normalizer)) {
      continue;
    }
    auto &upper =
        state->expr_upper_bounds[static_cast<std::size_t>(expr_id)];
    if (std::isfinite(upper) && upper <= observed_time) {
      continue;
    }
    upper = observed_time;
    state->expr_upper_normalizers[static_cast<std::size_t>(expr_id)] =
        normalizer;
  }
}

inline bool exact_sequence_states_equal(const ExactSequenceState &lhs,
                                        const ExactSequenceState &rhs) {
  if (lhs.lower_bound != rhs.lower_bound ||
      lhs.exact_times.size() != rhs.exact_times.size() ||
      lhs.upper_bounds.size() != rhs.upper_bounds.size()) {
    return false;
  }
  if (lhs.expr_upper_bounds.size() != rhs.expr_upper_bounds.size() ||
      lhs.expr_upper_normalizers.size() !=
          rhs.expr_upper_normalizers.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.exact_times.size(); ++i) {
    const bool lhs_na = std::isnan(lhs.exact_times[i]);
    const bool rhs_na = std::isnan(rhs.exact_times[i]);
    if (lhs_na || rhs_na) {
      if (lhs_na != rhs_na) {
        return false;
      }
      continue;
    }
    if (lhs.exact_times[i] != rhs.exact_times[i]) {
      return false;
    }
  }
  for (std::size_t i = 0; i < lhs.upper_bounds.size(); ++i) {
    if (lhs.upper_bounds[i] != rhs.upper_bounds[i]) {
      return false;
    }
  }
  for (std::size_t i = 0; i < lhs.expr_upper_bounds.size(); ++i) {
    if (lhs.expr_upper_bounds[i] != rhs.expr_upper_bounds[i]) {
      return false;
    }
  }
  for (std::size_t i = 0; i < lhs.expr_upper_normalizers.size(); ++i) {
    if (lhs.expr_upper_normalizers[i] != rhs.expr_upper_normalizers[i]) {
      return false;
    }
  }
  return true;
}

inline void append_ranked_frontier_entry(
    std::vector<ExactRankedFrontierEntry> *frontier,
    const std::vector<ExactSequenceState> &states,
    const semantic::Index state_index,
    const double probability) {
  if (!(probability > 0.0)) {
    return;
  }
  const auto state_pos = static_cast<std::size_t>(state_index);
  for (auto &existing : *frontier) {
    if (exact_sequence_states_equal(
            states[static_cast<std::size_t>(existing.state_index)],
            states[state_pos])) {
      existing.probability += probability;
      return;
    }
  }
  frontier->push_back(ExactRankedFrontierEntry{probability, state_index});
}

inline ExactSequenceState &ranked_sequence_state_slot(
    std::vector<ExactSequenceState> *states,
    const ExactVariantPlan &plan,
    const std::size_t index) {
  while (states->size() <= index) {
    states->push_back(make_exact_sequence_state(plan));
  }
  return (*states)[index];
}

inline void exact_sequence_ready_expr_normalizers(
    const ExactVariantPlan &plan,
    const ExactCompiledTransitionPlan &transition,
    ExactStepWorkspace *workspace,
    std::vector<double> *normalizers) {
  normalizers->clear();
  normalizers->reserve(transition.readiness_expr_ids.size());
  for (const auto expr_id : transition.readiness_expr_ids) {
    double normalizer = 0.0;
    if (expr_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(expr_id) <
            plan.sequence.expr_cdf_roots.size()) {
      const auto root_id =
          plan.sequence.expr_cdf_roots[static_cast<std::size_t>(expr_id)];
      if (root_id != semantic::kInvalidIndex) {
        normalizer =
            evaluate_compiled_math_root(
                plan.compiled_math,
                root_id,
                &workspace->target_workspace.compiled_math,
                &workspace->target_evaluator,
                nullptr,
                &workspace->target_workspace);
      }
    }
    normalizers->push_back(
        std::isfinite(normalizer) ? clamp_probability(normalizer) : 0.0);
  }
}

inline double exact_ranked_trigger_probability(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const ExactTrialView &obs,
    ExactStepWorkspace *step_workspace) {
  auto &used_outcomes = step_workspace->ranked_used_outcomes;
  auto &frontier = step_workspace->ranked_frontier;
  auto &next_frontier = step_workspace->ranked_next_frontier;
  auto &states = step_workspace->ranked_states;
  auto &next_states = step_workspace->ranked_next_states;
  used_outcomes.assign(plan.compiled_outcomes.size(), 0U);
  frontier.clear();
  next_frontier.clear();
  ranked_sequence_state_slot(&states, plan, 0) = step_workspace->initial_state;
  frontier.push_back(ExactRankedFrontierEntry{1.0, 0});

  for (std::size_t rank_idx = 0;
       rank_idx < static_cast<std::size_t>(obs.rank_count);
       ++rank_idx) {
    const auto outcome_code = exact_trial_view_outcome_code(obs, rank_idx);
    const auto target_outcome_index =
        plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
    const auto target_idx = static_cast<std::size_t>(target_outcome_index);
    if (used_outcomes[target_idx] != 0U) {
      return 0.0;
    }

    const auto &compiled_outcome = plan.compiled_outcomes[target_idx];
    next_frontier.clear();
    std::size_t next_state_count = 0;
    for (const auto &entry : frontier) {
      if (!(entry.probability > 0.0)) {
        continue;
      }
      const auto entry_state_index =
          static_cast<std::size_t>(entry.state_index);
      const ExactStepDistributionView step =
          evaluate_exact_step_distribution(
              plan,
              params,
              first_param_row,
              trigger_state,
              states[entry_state_index],
              target_outcome_index,
              exact_trial_view_rt(obs, rank_idx),
              &used_outcomes,
              true,
              step_workspace);
      if (step.transition_probabilities == nullptr) {
        continue;
      }
      for (std::size_t transition_idx = 0;
           transition_idx < step.transition_probabilities->size();
           ++transition_idx) {
        const double transition_probability =
            (*step.transition_probabilities)[transition_idx];
        if (!(transition_probability > 0.0)) {
          continue;
        }
        const auto candidate_index = next_state_count++;
        auto &candidate_state =
            ranked_sequence_state_slot(&next_states, plan, candidate_index);
        candidate_state = states[entry_state_index];
        const auto &transition = compiled_outcome.transitions[transition_idx];
        exact_sequence_ready_expr_normalizers(
            plan,
            transition,
            step_workspace,
            &step_workspace->ready_expr_normalizers);
        advance_exact_sequence_state(
            &candidate_state,
            transition,
            plan,
            exact_trial_view_rt(obs, rank_idx),
            &step_workspace->ready_expr_normalizers);
        append_ranked_frontier_entry(
            &next_frontier,
            next_states,
            static_cast<semantic::Index>(candidate_index),
            entry.probability * transition_probability);
      }
    }
    if (next_frontier.empty()) {
      return 0.0;
    }
    used_outcomes[target_idx] = 1U;
    frontier.swap(next_frontier);
    states.swap(next_states);
  }

  double total = 0.0;
  for (const auto &entry : frontier) {
    total += entry.probability;
  }
  return total;
}

inline double exact_ranked_loglik_for_trial(const ExactVariantPlan &plan,
                                            const ParamView &params,
                                            const int first_param_row,
                                            const ExactTrialView &obs,
                                            const double min_ll,
                                            ExactStepWorkspace *workspace = nullptr) {
  std::optional<ExactStepWorkspace> local_workspace;
  if (workspace == nullptr) {
    local_workspace.emplace(plan);
    workspace = &*local_workspace;
  }
  double total = 0.0;
  for (const auto &compiled_state : plan.trigger_state_table.states) {
    const auto trigger_state =
        exact_compiled_trigger_state_view(
            plan, params, first_param_row, compiled_state);
    if (!(trigger_state.weight > 0.0)) {
      continue;
    }
    total += trigger_state.weight *
             exact_ranked_trigger_probability(
                 plan,
                 params,
                 first_param_row,
                 trigger_state,
                 obs,
                 workspace);
  }
  if (!std::isfinite(total) || !(total > 0.0)) {
    return min_ll;
  }
  return std::log(total);
}

template <typename TrialSink>
inline void evaluate_exact_trials_cached(
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const int *ok,
    TrialSink &&sink) {
  const double *onset =
      layout.onset_col >= 0
          ? REAL(trusted_data_column(dataSEXP, layout.onset_col))
          : nullptr;
  ParamView params(paramsSEXP, onset);
  const auto table = read_prepared_data_view(dataSEXP, layout);
  const auto columns = make_exact_trial_columns(dataSEXP, layout);
  ExactStepWorkspacePool workspace_pool(plans.size());
  std::size_t param_row = 0;
  for (std::size_t trial_index = 0; trial_index < layout.trials.size(); ++trial_index) {
    const auto row = static_cast<R_xlen_t>(layout.trials[trial_index].start_row);
    const auto variant_index =
        variant_index_by_component_code[
            static_cast<std::size_t>(table.component[row])];

    const auto &plan = plans[static_cast<std::size_t>(variant_index)];
    const auto leaf_count =
        static_cast<std::size_t>(plan.program.layout.n_leaves);
    double value = min_ll;
    if (trial_is_selected(ok, trial_index)) {
      const auto obs = read_exact_observation_view(
          table,
          variant_index_by_component_code,
          layout,
          trial_index,
          columns);
      auto &workspace = workspace_pool.get(plans, variant_index);
      value =
          obs.rank_count == 1
              ? exact_loglik_for_trial(
                    plan,
                    params,
                    static_cast<int>(param_row),
                    exact_trial_view_outcome_code(obs, 0U),
                    exact_trial_view_rt(obs, 0U),
                    min_ll,
                    &workspace)
              : exact_ranked_loglik_for_trial(
                    plan,
                    params,
                    static_cast<int>(param_row),
                    obs,
                    min_ll,
                    &workspace);
    }
    sink(trial_index, value);
    param_row += leaf_count;
  }
}

} // namespace detail
} // namespace accumulatr::eval
