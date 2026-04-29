#pragma once

#include <array>
#include <limits>
#include <optional>

#include "eval_query.hpp"
#include "exact_competitor_union.hpp"
#include "trial_data.hpp"

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
  const auto &span = layout.spans[trial_index];
  const auto row = static_cast<R_xlen_t>(span.start_row);
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
    std::vector<ExactVariantPlan> *plans) {
  variant_index_by_component_code->assign(
      n_component_codes + 1U,
      semantic::kInvalidIndex);
  plans->clear();
  plans->reserve(compiled.variants.size());
  for (const auto &variant : compiled.variants) {
    if (variant.backend != compile::BackendKind::Exact) {
      continue;
    }
    const auto plan_index = static_cast<semantic::Index>(plans->size());
    const auto component_it = component_code_by_id.find(variant.component_id);
    if (component_it == component_code_by_id.end()) {
      throw std::runtime_error(
          "exact evaluator found no prepared component code for '" +
          variant.component_id + "'");
    }
    (*variant_index_by_component_code)[static_cast<std::size_t>(component_it->second)] =
        plan_index;
    plans->push_back(
      make_exact_variant_plan(
          runtime::lower_exact_variant(variant),
          outcome_code_by_label,
          n_outcome_codes));
  }
}

inline double exact_density_program_value(
    const ExactVariantPlan &plan,
    const ParamView &params,
    int first_param_row,
    semantic::Index target_idx,
    double rt,
    ExactStepWorkspace *workspace = nullptr);

inline double exact_no_response_program_value(
    const ExactVariantPlan &plan,
    const ParamView &params,
    int first_param_row,
    ExactStepWorkspace *workspace = nullptr);

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
  return exact_density_program_value(
      plan,
      params,
      first_param_row,
      target_idx,
      rt,
      workspace);
}

inline double exact_simple_race_leaf_q(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const semantic::Index leaf_index) {
  const auto leaf_pos = static_cast<std::size_t>(leaf_index);
  const auto &program = plan.lowered.program;
  const auto trigger_index = program.leaf_trigger_index[leaf_pos];
  if (trigger_index != semantic::kInvalidIndex) {
    const auto trigger_pos = static_cast<std::size_t>(trigger_index);
    if (static_cast<semantic::TriggerKind>(program.trigger_kind[trigger_pos]) ==
            semantic::TriggerKind::Shared &&
        trigger_state.shared_started[trigger_pos] <= 1U) {
      return trigger_state.shared_started[trigger_pos] == 1U ? 0.0 : 1.0;
    }
  }
  return params.q(first_param_row + leaf_index);
}

inline leaf::EventChannels exact_simple_race_leaf_channels_at(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const semantic::Index leaf_index,
    const double rt) {
  const auto leaf_pos = static_cast<std::size_t>(leaf_index);
  const auto &program = plan.lowered.program;
  const auto &desc = program.leaf_descriptors[leaf_pos];
  std::array<double, 8> local_params{};
  const int row = first_param_row + leaf_index;
  const int n_local = std::min<int>(desc.param_count, 8);
  for (int j = 0; j < n_local; ++j) {
    local_params[static_cast<std::size_t>(j)] = params.p(row, j);
  }
  const double q =
      exact_simple_race_leaf_q(plan, params, first_param_row, trigger_state, leaf_index);
  return standard_leaf_channels(
      desc.dist_kind,
      local_params.data(),
      n_local,
      q,
      params.t0(row),
      rt - desc.onset_abs_value);
}

inline double exact_terminal_leaf_survival_probability(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const semantic::Index source_id) {
  const auto source_pos = static_cast<std::size_t>(source_id);
  const auto &kernel = plan.source_kernels[source_pos];
  const auto leaf_pos = static_cast<std::size_t>(kernel.leaf_index);
  const auto &program = plan.lowered.program;
  const auto trigger_index = program.leaf_trigger_index[leaf_pos];
  if (trigger_index != semantic::kInvalidIndex) {
    const auto trigger_pos = static_cast<std::size_t>(trigger_index);
    if (static_cast<semantic::TriggerKind>(program.trigger_kind[trigger_pos]) ==
            semantic::TriggerKind::Shared &&
        trigger_state.shared_started[trigger_pos] <= 1U) {
      return trigger_state.shared_started[trigger_pos] == 1U ? 0.0 : 1.0;
    }
  }
  return clamp_probability(params.q(first_param_row + kernel.leaf_index));
}

inline double evaluate_exact_probability_trigger_op(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ExactProbabilityProgramSet &programs,
    const ParamView &params,
    const int first_param_row,
    const double observed_time,
    const ExactProbabilityOp &op,
    const ExactTriggerState &trigger_state,
    ExactStepWorkspace *workspace) {
  switch (op.kind) {
  case ExactProbabilityOpKind::Constant:
    return op.constant;

  case ExactProbabilityOpKind::Top1LeafRaceDensity: {
    double term = 1.0;
    for (semantic::Index i = 0; i < op.source_span.size; ++i) {
      const auto source_id =
          programs.source_ids[
              static_cast<std::size_t>(op.source_span.offset + i)];
      const auto &kernel =
          plan.source_kernels[static_cast<std::size_t>(source_id)];
      const auto channels = exact_simple_race_leaf_channels_at(
          plan,
          params,
          first_param_row,
          trigger_state,
          kernel.leaf_index,
          observed_time);
      const double factor =
          source_id == op.target_source_id
              ? safe_density(channels.pdf)
              : clamp_probability(channels.survival);
      term *= factor;
      if (!(term > 0.0)) {
        break;
      }
    }
    return term > 0.0 ? term : 0.0;
  }

  case ExactProbabilityOpKind::TerminalNoResponseProbability: {
    double terminal_survival = 1.0;
    for (semantic::Index i = 0; i < op.source_span.size; ++i) {
      const auto source_id =
          programs.source_ids[
              static_cast<std::size_t>(op.source_span.offset + i)];
      const double source_survival =
          exact_terminal_leaf_survival_probability(
              plan,
              params,
              first_param_row,
              trigger_state,
              source_id);
      terminal_survival *= source_survival;
      if (!(terminal_survival > 0.0)) {
        break;
      }
    }
    return clamp_probability(terminal_survival);
  }

  case ExactProbabilityOpKind::GenericTransitionDensity: {
    const ExactStepDistributionView step = evaluate_exact_step_distribution(
        plan,
        params,
        first_param_row,
        trigger_state,
        workspace->initial_state,
        op.outcome_index,
        observed_time,
        nullptr,
        false,
        workspace);
    return step.total_probability;
  }

  case ExactProbabilityOpKind::GenericTransitionProbability: {
    const double probability = integrate_to_infinity(
        [&](const double rt) {
          const ExactStepDistributionView step =
              evaluate_exact_step_distribution(
                  plan,
                  params,
                  first_param_row,
                  trigger_state,
                  workspace->initial_state,
                  op.outcome_index,
                  rt,
                  nullptr,
                  false,
                  workspace);
          return step.total_probability;
        });
    return clamp_probability(probability);
  }

  case ExactProbabilityOpKind::Integral: {
    const auto child =
        program.child_ops[static_cast<std::size_t>(op.children.offset)];
    const auto &child_op = program.ops[static_cast<std::size_t>(child)];
    const double probability = integrate_to_infinity(
        [&](const double rt) {
          return evaluate_exact_probability_trigger_op(
              plan,
              program,
              programs,
              params,
              first_param_row,
              rt,
              child_op,
              trigger_state,
              workspace);
        });
    return op.value_kind == ExactProbabilityValueKind::Probability
               ? clamp_probability(probability)
               : probability;
  }

  case ExactProbabilityOpKind::Log: {
    const auto child =
        program.child_ops[static_cast<std::size_t>(op.children.offset)];
    const double probability = evaluate_exact_probability_trigger_op(
        plan,
        program,
        programs,
        params,
        first_param_row,
        observed_time,
        program.ops[static_cast<std::size_t>(child)],
        trigger_state,
        workspace);
    return std::isfinite(probability) && probability > 0.0
               ? std::log(probability)
               : R_NegInf;
  }

  case ExactProbabilityOpKind::WeightedTriggerSum:
  case ExactProbabilityOpKind::RankedTransitionSequence:
    return 0.0;
  }
  return 0.0;
}

inline double evaluate_exact_probability_weighted_trigger_sum(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ExactProbabilityProgramSet &programs,
    const ParamView &params,
    const int first_param_row,
    const double observed_time,
    const ExactProbabilityOp &root,
    ExactStepWorkspace *workspace) {
  if (program.root_child != semantic::kInvalidIndex) {
    const auto &child_op =
        program.ops[static_cast<std::size_t>(program.root_child)];
    if (!program.requires_trigger_enumeration) {
      const double total = evaluate_exact_probability_trigger_op(
          plan,
          program,
          programs,
          params,
          first_param_row,
          observed_time,
          child_op,
          workspace->default_trigger_state,
          workspace);
      return root.value_kind == ExactProbabilityValueKind::Probability
                 ? clamp_probability(total)
                 : (std::isfinite(total) && total > 0.0 ? total : 0.0);
    }

    double total = 0.0;
    for (const auto &compiled_state : plan.trigger_state_table.states) {
      const auto state =
          exact_compiled_trigger_state_view(
              plan, params, first_param_row, compiled_state);
      if (!(state.weight > 0.0)) {
        continue;
      }
      total += state.weight * evaluate_exact_probability_trigger_op(
          plan,
          program,
          programs,
          params,
          first_param_row,
          observed_time,
          child_op,
          state,
          workspace);
    }
    return root.value_kind == ExactProbabilityValueKind::Probability
               ? clamp_probability(total)
               : (std::isfinite(total) && total > 0.0 ? total : 0.0);
  }

  if (!program.requires_trigger_enumeration) {
    double total = 0.0;
    for (semantic::Index i = 0; i < root.children.size; ++i) {
      const auto child =
          program.child_ops[
              static_cast<std::size_t>(root.children.offset + i)];
      total += evaluate_exact_probability_trigger_op(
          plan,
          program,
          programs,
          params,
          first_param_row,
          observed_time,
          program.ops[static_cast<std::size_t>(child)],
          workspace->default_trigger_state,
          workspace);
    }
    return root.value_kind == ExactProbabilityValueKind::Probability
               ? clamp_probability(total)
               : (std::isfinite(total) && total > 0.0 ? total : 0.0);
  }

  double total = 0.0;
  for (const auto &compiled_state : plan.trigger_state_table.states) {
    const auto state =
        exact_compiled_trigger_state_view(
            plan, params, first_param_row, compiled_state);
    if (!(state.weight > 0.0)) {
      continue;
    }
    double state_value = 0.0;
    for (semantic::Index i = 0; i < root.children.size; ++i) {
      const auto child =
          program.child_ops[
              static_cast<std::size_t>(root.children.offset + i)];
      state_value += evaluate_exact_probability_trigger_op(
          plan,
          program,
          programs,
          params,
          first_param_row,
          observed_time,
          program.ops[static_cast<std::size_t>(child)],
          state,
          workspace);
    }
    total += state.weight * state_value;
  }
  return root.value_kind == ExactProbabilityValueKind::Probability
             ? clamp_probability(total)
             : (std::isfinite(total) && total > 0.0 ? total : 0.0);
}

inline double evaluate_exact_probability_program(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ParamView &params,
    const int first_param_row,
    const double observed_time,
    ExactStepWorkspace *workspace = nullptr) {
  std::optional<ExactStepWorkspace> local_workspace;
  if (workspace == nullptr) {
    local_workspace.emplace(plan);
    workspace = &*local_workspace;
  }
  const auto evaluate_program_op =
      [&](const auto &self,
          const ExactProbabilityOp &op,
          const double op_time) -> double {
    switch (op.kind) {
    case ExactProbabilityOpKind::WeightedTriggerSum:
      return evaluate_exact_probability_weighted_trigger_sum(
          plan,
          program,
          plan.probability_programs,
          params,
          first_param_row,
          op_time,
          op,
          workspace);
    case ExactProbabilityOpKind::Integral: {
      const auto child =
          program.child_ops[static_cast<std::size_t>(op.children.offset)];
      const auto &child_op = program.ops[static_cast<std::size_t>(child)];
      const double probability = integrate_to_infinity(
          [&](const double rt) {
            return self(self, child_op, rt);
          });
      return op.value_kind == ExactProbabilityValueKind::Probability
                 ? clamp_probability(probability)
                 : probability;
    }
    case ExactProbabilityOpKind::Log: {
      const auto child =
          program.child_ops[static_cast<std::size_t>(op.children.offset)];
      const double probability =
          self(self, program.ops[static_cast<std::size_t>(child)], op_time);
      return std::isfinite(probability) && probability > 0.0
                 ? std::log(probability)
                 : R_NegInf;
    }
    case ExactProbabilityOpKind::Top1LeafRaceDensity:
    case ExactProbabilityOpKind::TerminalNoResponseProbability:
    case ExactProbabilityOpKind::GenericTransitionDensity:
    case ExactProbabilityOpKind::GenericTransitionProbability:
      return evaluate_exact_probability_trigger_op(
          plan,
          program,
          plan.probability_programs,
          params,
          first_param_row,
          op_time,
          op,
          workspace->default_trigger_state,
          workspace);
    case ExactProbabilityOpKind::Constant:
      return op.constant;
    case ExactProbabilityOpKind::RankedTransitionSequence:
      return 0.0;
    }
    return 0.0;
  };
  const auto &root = program.ops[static_cast<std::size_t>(program.root)];
  switch (root.kind) {
  case ExactProbabilityOpKind::WeightedTriggerSum:
    return evaluate_exact_probability_weighted_trigger_sum(
        plan,
        program,
        plan.probability_programs,
        params,
        first_param_row,
        observed_time,
        root,
        workspace);
  case ExactProbabilityOpKind::Integral:
  case ExactProbabilityOpKind::Log:
    return evaluate_program_op(evaluate_program_op, root, observed_time);
  case ExactProbabilityOpKind::Top1LeafRaceDensity:
  case ExactProbabilityOpKind::TerminalNoResponseProbability:
  case ExactProbabilityOpKind::Constant:
  case ExactProbabilityOpKind::GenericTransitionDensity:
  case ExactProbabilityOpKind::GenericTransitionProbability:
  case ExactProbabilityOpKind::RankedTransitionSequence:
    return evaluate_program_op(evaluate_program_op, root, observed_time);
  }
  return 0.0;
}

inline double exact_density_program_value(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const semantic::Index target_idx,
    const double rt,
    ExactStepWorkspace *workspace) {
  const auto program_index =
      plan.probability_programs.density_by_outcome[
          static_cast<std::size_t>(target_idx)];
  const auto &program =
      plan.probability_programs.programs[
          static_cast<std::size_t>(program_index)];
  return evaluate_exact_probability_program(
      plan,
      program,
      params,
      first_param_row,
      rt,
      workspace);
}

inline double exact_no_response_program_value(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    ExactStepWorkspace *workspace) {
  const auto &program =
      plan.probability_programs.programs[
          static_cast<std::size_t>(
              plan.probability_programs.no_response_probability)];
  return evaluate_exact_probability_program(
      plan,
      program,
      params,
      first_param_row,
      NA_REAL,
      workspace);
}

inline double exact_finite_probability_program_value(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const semantic::Index target_idx,
    ExactStepWorkspace *workspace) {
  const auto program_index =
      plan.probability_programs.finite_probability_by_outcome[
          static_cast<std::size_t>(target_idx)];
  const auto &program =
      plan.probability_programs.programs[
          static_cast<std::size_t>(program_index)];
  return evaluate_exact_probability_program(
      plan,
      program,
      params,
      first_param_row,
      NA_REAL,
      workspace);
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
    const ExactRuntimeScenarioTransitionPlan &transition,
    const ExactVariantPlan &plan,
    const double observed_time,
    const std::vector<double> *ready_expr_normalizers = nullptr) {
  if (state == nullptr) {
    return;
  }
  state->lower_bound = observed_time;
  if (transition.active_source_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(transition.active_source_id) <
          state->exact_times.size()) {
    state->exact_times[static_cast<std::size_t>(transition.active_source_id)] =
        observed_time;
  }
  for (semantic::Index i = 0; i < transition.before_source_span.size; ++i) {
    const auto source_id =
        plan.scenario_source_ids[
            static_cast<std::size_t>(
                transition.before_source_span.offset + i)];
    if (source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(source_id) >= state->upper_bounds.size()) {
      continue;
    }
    auto &upper = state->upper_bounds[static_cast<std::size_t>(source_id)];
    upper = std::isfinite(upper) ? std::min(upper, observed_time)
                                 : observed_time;
  }
  for (semantic::Index i = 0; i < transition.ready_expr_span.size; ++i) {
    const auto expr_id =
        plan.scenario_expr_ids[
            static_cast<std::size_t>(
                transition.ready_expr_span.offset + i)];
    if (expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(expr_id) >=
            state->expr_upper_bounds.size() ||
        static_cast<std::size_t>(expr_id) >=
            state->expr_upper_normalizers.size()) {
      continue;
    }
    const std::size_t ready_pos = static_cast<std::size_t>(i);
    if (ready_expr_normalizers == nullptr ||
        ready_pos >= ready_expr_normalizers->size()) {
      continue;
    }
    const double normalizer = (*ready_expr_normalizers)[ready_pos];
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
    const ExactRuntimeScenarioTransitionPlan &transition,
    ExactStepWorkspace *workspace,
    std::vector<double> *normalizers) {
  normalizers->clear();
  normalizers->reserve(
      static_cast<std::size_t>(transition.ready_expr_span.size));
  for (semantic::Index i = 0; i < transition.ready_expr_span.size; ++i) {
    const auto expr_id =
        plan.scenario_expr_ids[
            static_cast<std::size_t>(
                transition.ready_expr_span.offset + i)];
    double normalizer = 0.0;
    if (expr_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(expr_id) <
            plan.sequence_expr_cdf_roots.size()) {
      const auto root_id =
          plan.sequence_expr_cdf_roots[static_cast<std::size_t>(expr_id)];
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
  used_outcomes.assign(plan.outcomes.size(), 0U);
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

    const auto &successors =
        plan.runtime.outcomes[target_idx].successor_distribution;
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
        const auto &transition = successors.transitions[transition_idx];
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

inline SEXP evaluate_exact_trials_cached(
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const int *ok = nullptr) {
  ParamView params(paramsSEXP);
  const auto table = read_prepared_data_view(dataSEXP, layout);
  const int max_rank = layout.max_rank;
  ExactTrialColumns columns;
  columns.labels.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  columns.times.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  columns.labels[1] =
      INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
  columns.times[1] =
      REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
  for (int rank = 2; rank <= max_rank; ++rank) {
    columns.labels[static_cast<std::size_t>(rank)] =
        INTEGER(trusted_data_column(
            dataSEXP,
            layout.label_cols[static_cast<std::size_t>(rank)]));
    columns.times[static_cast<std::size_t>(rank)] =
        REAL(trusted_data_column(
            dataSEXP,
            layout.time_cols[static_cast<std::size_t>(rank)]));
  }
  Rcpp::NumericVector loglik(layout.spans.size(), min_ll);
  std::vector<runtime::TrialBlock> blocks;
  ExactStepWorkspacePool workspace_pool(plans.size());
  std::size_t param_row = 0;
  runtime::TrialBlock current_block;
  bool have_block = false;
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
    const auto variant_index =
        variant_index_by_component_code[
            static_cast<std::size_t>(table.component[row])];
    if (!have_block || current_block.variant_index != variant_index) {
      if (have_block) {
        blocks.push_back(current_block);
      }
      current_block.variant_index = variant_index;
      current_block.start_row = static_cast<int>(trial_index);
      current_block.row_count = 1;
      have_block = true;
    } else {
      ++current_block.row_count;
    }

    const auto &plan = plans[static_cast<std::size_t>(variant_index)];
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    if (!trial_is_selected(ok, trial_index)) {
      param_row += leaf_count;
      continue;
    }
    const auto obs = read_exact_observation_view(
        table,
        variant_index_by_component_code,
        layout,
        trial_index,
        columns);
    auto &workspace = workspace_pool.get(plans, variant_index);
    loglik[static_cast<R_xlen_t>(trial_index)] =
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
    param_row += leaf_count;
  }
  if (have_block) {
    blocks.push_back(current_block);
  }

  Rcpp::List block_list(blocks.size());
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    const auto &block = blocks[i];
    const auto &plan = plans[static_cast<std::size_t>(block.variant_index)];
    block_list[i] = Rcpp::List::create(
        Rcpp::Named("component_id") =
            plan.lowered.component_id,
        Rcpp::Named("start_row") = block.start_row + 1,
        Rcpp::Named("row_count") = block.row_count);
  }

  const double total_loglik = aggregate_trial_loglik(loglik, layout);

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("blocks") = block_list);
}

} // namespace detail
} // namespace accumulatr::eval
