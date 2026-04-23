#pragma once

#include "exact_step.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactTrialView {
  semantic::Index variant_index{semantic::kInvalidIndex};
  R_xlen_t row{0};
  int rank_count{0};
  const std::vector<Rcpp::IntegerVector> *labels{nullptr};
  const std::vector<Rcpp::NumericVector> *times{nullptr};
};

inline semantic::Index exact_trial_view_outcome_code(const ExactTrialView &view,
                                                     const std::size_t rank_idx) {
  return (*(view.labels))[rank_idx + 1U][view.row];
}

inline double exact_trial_view_rt(const ExactTrialView &view,
                                  const std::size_t rank_idx) {
  return (*(view.times))[rank_idx + 1U][view.row];
}

inline ExactTrialView read_exact_observation_view(
    const Rcpp::DataFrame &data,
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const PreparedTrialLayout &layout,
    const std::size_t trial_index,
    const std::vector<Rcpp::IntegerVector> &labels,
    const std::vector<Rcpp::NumericVector> &times) {
  const auto table = read_prepared_data_view(data);
  const auto &span = layout.spans.at(trial_index);
  const auto row = static_cast<R_xlen_t>(span.start_row);
  const int component_code = table.component[row];
  ExactTrialView obs;
  obs.variant_index =
      variant_index_by_component_code[static_cast<std::size_t>(component_code)];
  obs.row = row;
  obs.labels = &labels;
  obs.times = &times;
  for (int rank = 1; rank <= layout.max_rank; ++rank) {
    const auto &rank_labels = labels[static_cast<std::size_t>(rank)];
    const auto &rank_times = times[static_cast<std::size_t>(rank)];
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

inline double exact_loglik_for_trial(const ExactVariantPlan &plan,
                                     const ParamView &params,
                                     const int first_param_row,
                                     const semantic::Index outcome_code,
                                     const double rt,
                                     const PreparedTrialLayout *layout,
                                     const double min_ll) {
  struct ExactUnrankedTrialKernel {
    const ExactVariantPlan &plan;
    const ParamView &params;
    int first_param_row;
    semantic::Index target_idx;
    std::vector<ExactTriggerState> trigger_states;
    ExactSequenceState initial_state;
    ExactStepWorkspace step_workspace;

    ExactUnrankedTrialKernel(const ExactVariantPlan &plan_,
                             const ParamView &params_,
                             const int first_param_row_,
                             const semantic::Index target_idx_)
        : plan(plan_),
          params(params_),
          first_param_row(first_param_row_),
          target_idx(target_idx_),
          trigger_states(
              enumerate_trigger_states(plan_, params_, first_param_row_)),
          initial_state(make_exact_sequence_state(plan_)),
          step_workspace(plan_) {}

    double density(const double rt,
                   const quadrature::FiniteBatch *step_batch = nullptr) {
      if (!std::isfinite(rt) || !(rt > 0.0)) {
        return 0.0;
      }
      double total = 0.0;
      for (const auto &trigger_state : trigger_states) {
        if (!(trigger_state.weight > 0.0)) {
          continue;
        }
        const ExactStepResult step = evaluate_exact_step(
            plan,
            params,
            first_param_row,
            trigger_state,
            initial_state,
            target_idx,
            rt,
            nullptr,
            false,
            step_batch,
            &step_workspace);
        total += trigger_state.weight * step.total_probability;
      }
      return std::isfinite(total) && total > 0.0 ? total : 0.0;
    }

    double loglik(const double rt,
                  const quadrature::FiniteBatch *step_batch,
                  const double min_ll) {
      const double total = density(rt, step_batch);
      return total > 0.0 ? std::log(total) : min_ll;
    }
  };

  const auto target_idx =
      plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
  const auto *step_batch = layout == nullptr ? nullptr : find_finite_batch(*layout, rt);
  ExactUnrankedTrialKernel kernel(plan, params, first_param_row, target_idx);
  return kernel.loglik(rt, step_batch, min_ll);
}

inline double exact_ranked_sequence_probability(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const ExactTrialView &obs,
    const std::size_t rank_idx,
    const ExactSequenceState &sequence_state,
    std::vector<std::uint8_t> *used_outcomes,
    const PreparedTrialLayout *layout,
    ExactStepWorkspace *step_workspace) {
  if (rank_idx >= static_cast<std::size_t>(obs.rank_count)) {
    return 1.0;
  }

  const auto outcome_code = exact_trial_view_outcome_code(obs, rank_idx);
  const auto target_outcome_index =
      plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
  const auto target_idx = static_cast<std::size_t>(target_outcome_index);
  if ((*used_outcomes)[target_idx] != 0U) {
    return 0.0;
  }
  const ExactStepResult step = evaluate_exact_step(
      plan,
      params,
      first_param_row,
      trigger_state,
      sequence_state,
      target_outcome_index,
      exact_trial_view_rt(obs, rank_idx),
      used_outcomes,
      true,
      layout == nullptr ? nullptr
                        : find_finite_batch(*layout, exact_trial_view_rt(obs, rank_idx)),
      step_workspace);
  double total = 0.0;
  for (const auto &branch : step.branches) {
    (*used_outcomes)[target_idx] = 1U;
    const double tail_prob = exact_ranked_sequence_probability(
        plan,
        params,
        first_param_row,
        trigger_state,
        obs,
        rank_idx + 1U,
        branch.next_state,
        used_outcomes,
        layout,
        step_workspace);
    (*used_outcomes)[target_idx] = 0U;
    total += branch.probability * tail_prob;
  }
  return total;
}

inline double exact_ranked_loglik_for_trial(const ExactVariantPlan &plan,
                                            const ParamView &params,
                                            const int first_param_row,
                                            const ExactTrialView &obs,
                                            const PreparedTrialLayout *layout,
                                            const double min_ll) {
  if (obs.rank_count == 1) {
    return exact_loglik_for_trial(
        plan,
        params,
        first_param_row,
        exact_trial_view_outcome_code(obs, 0U),
        exact_trial_view_rt(obs, 0U),
        layout,
        min_ll);
  }

  const auto trigger_states = enumerate_trigger_states(plan, params, first_param_row);
  double total = 0.0;
  ExactStepWorkspace step_workspace(plan);
  for (const auto &trigger_state : trigger_states) {
    if (!(trigger_state.weight > 0.0)) {
      continue;
    }
    std::vector<std::uint8_t> used_outcomes(plan.outcomes.size(), 0U);
    ExactSequenceState state = make_exact_sequence_state(plan);
    total += trigger_state.weight *
             exact_ranked_sequence_probability(
                 plan,
                 params,
                 first_param_row,
                 trigger_state,
                 obs,
                 0U,
                 state,
                 &used_outcomes,
                 layout,
                 &step_workspace);
  }
  if (!std::isfinite(total) || !(total > 0.0)) {
    return min_ll;
  }
  return std::log(total);
}

inline bool exact_plan_supports_observation(const ExactVariantPlan &plan,
                                            const ExactTrialView &obs) {
  for (std::size_t rank_idx = 0; rank_idx < static_cast<std::size_t>(obs.rank_count);
       ++rank_idx) {
    const auto outcome_code = exact_trial_view_outcome_code(obs, rank_idx);
    if (outcome_code <= 0 ||
        static_cast<std::size_t>(outcome_code) >= plan.outcome_index_by_code.size() ||
        plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)] ==
            semantic::kInvalidIndex) {
      return false;
    }
  }
  return true;
}

inline Rcpp::NumericVector evaluate_exact_loglik_queries_cached(
    const std::vector<ExactVariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectLoglikQuery> &queries,
    const double min_ll,
    const std::vector<std::vector<int>> *row_maps = nullptr) {
  Rcpp::NumericVector out(queries.size(), min_ll);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    const int *row_map =
        row_maps != nullptr && query.row_map_index != semantic::kInvalidIndex
            ? row_maps->at(static_cast<std::size_t>(query.row_map_index)).data()
            : nullptr;
    ParamView params(paramsSEXP, row_map);
    const int first_param_row =
        row_map == nullptr
            ? static_cast<int>(
                  layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row)
            : 0;
    const auto target_idx =
        plan.outcome_index_by_code[static_cast<std::size_t>(query.outcome_code)];
    if (target_idx == semantic::kInvalidIndex) {
      out[static_cast<R_xlen_t>(i)] = min_ll;
      continue;
    }
    const auto *step_batch = find_finite_batch(layout, query.rt);
    struct ExactUnrankedQueryKernel {
      const ExactVariantPlan &plan;
      const ParamView &params;
      int first_param_row;
      semantic::Index target_idx;
      std::vector<ExactTriggerState> trigger_states;
      ExactSequenceState initial_state;
      ExactStepWorkspace step_workspace;

      ExactUnrankedQueryKernel(const ExactVariantPlan &plan_,
                               const ParamView &params_,
                               const int first_param_row_,
                               const semantic::Index target_idx_)
          : plan(plan_),
            params(params_),
            first_param_row(first_param_row_),
            target_idx(target_idx_),
            trigger_states(
                enumerate_trigger_states(plan_, params_, first_param_row_)),
            initial_state(make_exact_sequence_state(plan_)),
            step_workspace(plan_) {}

      double loglik(const double rt,
                    const quadrature::FiniteBatch *step_batch,
                    const double min_ll) {
        if (!std::isfinite(rt) || !(rt > 0.0)) {
          return min_ll;
        }
        double total = 0.0;
        for (const auto &trigger_state : trigger_states) {
          if (!(trigger_state.weight > 0.0)) {
            continue;
          }
          const ExactStepResult step = evaluate_exact_step(
              plan,
              params,
              first_param_row,
              trigger_state,
              initial_state,
              target_idx,
              rt,
              nullptr,
              false,
              step_batch,
              &step_workspace);
          total += trigger_state.weight * step.total_probability;
        }
        return std::isfinite(total) && total > 0.0 ? std::log(total) : min_ll;
      }
    } kernel(plan, params, first_param_row, target_idx);
    out[static_cast<R_xlen_t>(i)] = kernel.loglik(query.rt, step_batch, min_ll);
  }
  return out;
}

inline Rcpp::NumericVector evaluate_exact_probability_queries_cached(
    const std::vector<ExactVariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectProbabilityQuery> &queries,
    const std::vector<std::vector<int>> *row_maps = nullptr) {
  Rcpp::NumericVector out(queries.size(), 0.0);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    const int *row_map =
        row_maps != nullptr && query.row_map_index != semantic::kInvalidIndex
            ? row_maps->at(static_cast<std::size_t>(query.row_map_index)).data()
            : nullptr;
    ParamView params(paramsSEXP, row_map);
    const int first_param_row =
        row_map == nullptr
            ? static_cast<int>(
                  layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row)
            : 0;
    const auto target_idx =
        plan.outcome_index_by_code[static_cast<std::size_t>(query.outcome_code)];
    if (target_idx == semantic::kInvalidIndex) {
      out[static_cast<R_xlen_t>(i)] = 0.0;
      continue;
    }
    struct ExactUnrankedProbabilityKernel {
      const ExactVariantPlan &plan;
      const ParamView &params;
      int first_param_row;
      semantic::Index target_idx;
      std::vector<ExactTriggerState> trigger_states;
      ExactSequenceState initial_state;
      ExactStepWorkspace step_workspace;

      ExactUnrankedProbabilityKernel(const ExactVariantPlan &plan_,
                                     const ParamView &params_,
                                     const int first_param_row_,
                                     const semantic::Index target_idx_)
          : plan(plan_),
            params(params_),
            first_param_row(first_param_row_),
            target_idx(target_idx_),
            trigger_states(
                enumerate_trigger_states(plan_, params_, first_param_row_)),
            initial_state(make_exact_sequence_state(plan_)),
            step_workspace(plan_) {}

      double density(const double rt) {
        if (!std::isfinite(rt) || !(rt > 0.0)) {
          return 0.0;
        }
        double total = 0.0;
        for (const auto &trigger_state : trigger_states) {
          if (!(trigger_state.weight > 0.0)) {
            continue;
          }
          const ExactStepResult step = evaluate_exact_step(
              plan,
              params,
              first_param_row,
              trigger_state,
              initial_state,
              target_idx,
              rt,
              nullptr,
              false,
              nullptr,
              &step_workspace);
          total += trigger_state.weight * step.total_probability;
        }
        return std::isfinite(total) && total > 0.0 ? total : 0.0;
      }
    } kernel(plan, params, first_param_row, target_idx);
    out[static_cast<R_xlen_t>(i)] = integrate_to_infinity(
        [&](const double rt) {
          return kernel.density(rt);
        });
  }
  return out;
}

inline Rcpp::List evaluate_exact_outcome_probabilities_cached(
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP) {
  ParamView params(paramsSEXP);
  Rcpp::DataFrame data(dataSEXP);
  const auto queries =
      collapse_probability_queries(data, variant_index_by_component_code, layout);
  const auto blocks = build_variant_blocks(queries);

  Rcpp::NumericVector prob(queries.size());
  std::size_t param_row = 0;
  for (const auto &block : blocks) {
    const auto &plan = plans[static_cast<std::size_t>(block.variant_index)];
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    for (int local = 0; local < block.row_count; ++local) {
      const auto query_idx = static_cast<std::size_t>(block.start_row + local);
      const auto target_idx = plan.outcome_index_by_code[static_cast<std::size_t>(
          queries[query_idx].outcome_code)];
      struct ExactUnrankedProbabilityKernel {
        const ExactVariantPlan &plan;
        const ParamView &params;
        int first_param_row;
        semantic::Index target_idx;
        std::vector<ExactTriggerState> trigger_states;
        ExactSequenceState initial_state;
        ExactStepWorkspace step_workspace;

        ExactUnrankedProbabilityKernel(const ExactVariantPlan &plan_,
                                       const ParamView &params_,
                                       const int first_param_row_,
                                       const semantic::Index target_idx_)
            : plan(plan_),
              params(params_),
              first_param_row(first_param_row_),
              target_idx(target_idx_),
              trigger_states(
                  enumerate_trigger_states(plan_, params_, first_param_row_)),
              initial_state(make_exact_sequence_state(plan_)),
              step_workspace(plan_) {}

        double density(const double rt) {
          if (!std::isfinite(rt) || !(rt > 0.0)) {
            return 0.0;
          }
          double total = 0.0;
          for (const auto &trigger_state : trigger_states) {
            if (!(trigger_state.weight > 0.0)) {
              continue;
            }
            const ExactStepResult step = evaluate_exact_step(
                plan,
                params,
                first_param_row,
                trigger_state,
                initial_state,
                target_idx,
                rt,
                nullptr,
                false,
                nullptr,
                &step_workspace);
            total += trigger_state.weight * step.total_probability;
          }
          return std::isfinite(total) && total > 0.0 ? total : 0.0;
        }
      } kernel(plan, params, static_cast<int>(param_row), target_idx);
      prob[static_cast<R_xlen_t>(query_idx)] = integrate_to_infinity(
          [&](const double rt) {
            return kernel.density(rt);
          });
      param_row += leaf_count;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = prob,
      Rcpp::Named("n_trials") = static_cast<int>(queries.size()),
      Rcpp::Named("n_variants") = static_cast<int>(plans.size()));
}

inline SEXP evaluate_exact_trials_cached(
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll) {
  ParamView params(paramsSEXP);
  Rcpp::DataFrame data(dataSEXP);
  const int max_rank = layout.max_rank;
  std::vector<Rcpp::IntegerVector> labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> times(static_cast<std::size_t>(max_rank + 1));
  labels[1] = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
  times[1] = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  for (int rank = 2; rank <= max_rank; ++rank) {
    labels[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::IntegerVector>(data[("R" + std::to_string(rank)).c_str()]);
    times[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::NumericVector>(data[("rt" + std::to_string(rank)).c_str()]);
  }
  Rcpp::NumericVector loglik(layout.spans.size(), min_ll);
  std::vector<runtime::TrialBlock> blocks;
  std::size_t param_row = 0;
  runtime::TrialBlock current_block;
  bool have_block = false;
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    const auto obs = read_exact_observation_view(
        data,
        variant_index_by_component_code,
        layout,
        trial_index,
        labels,
        times);
    const auto variant_index = obs.variant_index;
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
    if (!exact_plan_supports_observation(plan, obs)) {
      loglik[static_cast<R_xlen_t>(trial_index)] = min_ll;
      param_row += leaf_count;
      continue;
    }
    loglik[static_cast<R_xlen_t>(trial_index)] = exact_ranked_loglik_for_trial(
        plan,
        params,
        static_cast<int>(param_row),
        obs,
        &layout,
        min_ll);
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

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    total_loglik += static_cast<double>(loglik[i]);
  }

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("blocks") = block_list);
}

} // namespace detail
} // namespace accumulatr::eval
