#pragma once

#include "exact_step.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

inline std::vector<ExactTrialObservation> collapse_exact_observations(
    const Rcpp::DataFrame &data,
    const std::unordered_map<std::string, semantic::Index> &variant_index_by_component) {
  const auto n_rows = data.nrows();
  if (n_rows == 0) {
    return {};
  }
  if (!data.containsElementNamed("R") || !data.containsElementNamed("rt")) {
    throw std::runtime_error("exact evaluator data must include columns R and rt");
  }

  const int max_rank = detect_rank_count(data, "exact evaluator");

  std::vector<Rcpp::CharacterVector> labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> times(static_cast<std::size_t>(max_rank + 1));
  labels[1] = Rcpp::as<Rcpp::CharacterVector>(data["R"]);
  times[1] = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  for (int rank = 2; rank <= max_rank; ++rank) {
    labels[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::CharacterVector>(data[("R" + std::to_string(rank)).c_str()]);
    times[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::NumericVector>(data[("rt" + std::to_string(rank)).c_str()]);
  }

  const auto table = read_trial_table_view(data);

  std::vector<ExactTrialObservation> out;
  out.reserve(static_cast<std::size_t>(n_rows));
  int last_trial = NA_INTEGER;
  for (R_xlen_t i = 0; i < n_rows; ++i) {
    const int trial_id = table.trial[i];
    if (i > 0 && trial_id == last_trial) {
      continue;
    }
    last_trial = trial_id;

    const std::string component_id = component_id_at_row(table, i);
    const auto it = variant_index_by_component.find(component_id);
    if (it == variant_index_by_component.end()) {
      throw std::runtime_error(
          "exact evaluator found no exact variant for component '" + component_id + "'");
    }

    ExactTrialObservation obs;
    obs.variant_index = it->second;
    obs.component_id = component_id;
    double prev_rt = -std::numeric_limits<double>::infinity();
    std::unordered_map<std::string, std::uint8_t> seen_labels;
    for (int rank = 1; rank <= max_rank; ++rank) {
      const auto &rank_labels = labels[static_cast<std::size_t>(rank)];
      const auto &rank_times = times[static_cast<std::size_t>(rank)];
      const SEXP label_sexp = rank_labels[i];
      const bool label_missing = label_sexp == NA_STRING;
      const double rt = rank_times[i];
      const bool time_missing = Rcpp::NumericVector::is_na(rt);
      if (label_missing && time_missing) {
        break;
      }
      if (label_missing || time_missing || !std::isfinite(rt) || !(rt > prev_rt)) {
        obs.valid = false;
        break;
      }
      const auto label = Rcpp::as<std::string>(rank_labels[i]);
      if (seen_labels.find(label) != seen_labels.end()) {
        obs.valid = false;
        break;
      }
      seen_labels.emplace(label, 1U);
      obs.ranks.push_back(ExactObservedRank{label, rt});
      prev_rt = rt;
    }
    if (obs.ranks.empty()) {
      obs.valid = false;
    }
    out.push_back(std::move(obs));
  }
  return out;
}

inline double exact_loglik_for_trial(const ExactVariantPlan &plan,
                                     const ParamView &params,
                                     const int first_param_row,
                                     const std::string &label,
                                     const double rt,
                                     const double min_ll) {
  if (!std::isfinite(rt) || !(rt > 0.0)) {
    return min_ll;
  }
  const auto target_it = plan.outcome_index.find(label);
  if (target_it == plan.outcome_index.end()) {
    return min_ll;
  }
  const auto target_idx = target_it->second;
  const auto trigger_states = enumerate_trigger_states(plan, params, first_param_row);
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
        ExactSequenceState{},
        target_idx,
        rt);
    total += trigger_state.weight * step.total_probability;
  }

  if (!std::isfinite(total) || !(total > 0.0)) {
    return min_ll;
  }
  return std::log(total);
}

inline double exact_ranked_sequence_probability(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const std::vector<ExactObservedRank> &ranks,
    const std::size_t rank_idx,
    const ExactSequenceState &sequence_state,
    std::vector<std::uint8_t> *used_outcomes) {
  if (rank_idx >= ranks.size()) {
    return 1.0;
  }

  const auto &rank = ranks[rank_idx];
  const auto target_it = plan.outcome_index.find(rank.outcome_label);
  if (target_it == plan.outcome_index.end()) {
    return 0.0;
  }
  const auto target_idx = static_cast<std::size_t>(target_it->second);
  if ((*used_outcomes)[target_idx] != 0U) {
    return 0.0;
  }
  const ExactStepResult step = evaluate_exact_step(
      plan,
      params,
      first_param_row,
      trigger_state,
      sequence_state,
      target_it->second,
      rank.rt,
      used_outcomes,
      true);
  double total = 0.0;
  for (const auto &branch : step.branches) {
    (*used_outcomes)[target_idx] = 1U;
    const double tail_prob = exact_ranked_sequence_probability(
        plan,
        params,
        first_param_row,
        trigger_state,
        ranks,
        rank_idx + 1U,
        branch.next_state,
        used_outcomes);
    (*used_outcomes)[target_idx] = 0U;
    total += branch.probability * tail_prob;
  }
  return total;
}

inline double exact_ranked_loglik_for_trial(const ExactVariantPlan &plan,
                                            const ParamView &params,
                                            const int first_param_row,
                                            const ExactTrialObservation &obs,
                                            const double min_ll) {
  if (!obs.valid || obs.ranks.empty()) {
    return min_ll;
  }
  if (obs.ranks.size() == 1U) {
    return exact_loglik_for_trial(
        plan,
        params,
        first_param_row,
        obs.ranks.front().outcome_label,
        obs.ranks.front().rt,
        min_ll);
  }

  const auto trigger_states = enumerate_trigger_states(plan, params, first_param_row);
  double total = 0.0;
  for (const auto &trigger_state : trigger_states) {
    if (!(trigger_state.weight > 0.0)) {
      continue;
    }
    std::vector<std::uint8_t> used_outcomes(plan.outcomes.size(), 0U);
    ExactSequenceState state;
    total += trigger_state.weight *
             exact_ranked_sequence_probability(
                 plan, params, first_param_row, trigger_state, obs.ranks, 0U, state, &used_outcomes);
  }
  if (!std::isfinite(total) || !(total > 0.0)) {
    return min_ll;
  }
  return std::log(total);
}

inline Rcpp::List evaluate_exact_outcome_probabilities(
    const compile::CompiledModel &compiled,
    SEXP paramsSEXP,
    SEXP dataSEXP) {
  std::unordered_map<std::string, semantic::Index> variant_index_by_component;
  std::vector<runtime::LoweredExactVariant> lowered_variants;
  lowered_variants.reserve(compiled.variants.size());

  for (const auto &variant : compiled.variants) {
    if (variant.backend != compile::BackendKind::Exact) {
      continue;
    }
    variant_index_by_component.emplace(
        variant.component_id,
        static_cast<semantic::Index>(lowered_variants.size()));
    lowered_variants.push_back(runtime::lower_exact_variant(variant));
  }
  if (lowered_variants.empty()) {
    throw std::runtime_error("exact probability evaluator found no exact variants");
  }

  ParamView params(paramsSEXP);
  Rcpp::DataFrame data(dataSEXP);
  const auto queries = collapse_probability_queries(data, variant_index_by_component);
  const auto blocks = build_variant_blocks(queries);

  std::vector<ExactVariantPlan> plans;
  plans.reserve(lowered_variants.size());
  for (const auto &variant : lowered_variants) {
    plans.push_back(make_exact_variant_plan(variant));
  }

  Rcpp::NumericVector prob(queries.size());
  std::size_t param_row = 0;
  for (const auto &block : blocks) {
    const auto &plan = plans[static_cast<std::size_t>(block.variant_index)];
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    for (int local = 0; local < block.row_count; ++local) {
      const auto query_idx = static_cast<std::size_t>(block.start_row + local);
      if (param_row + leaf_count > static_cast<std::size_t>(params.matrix.nrow())) {
        throw std::runtime_error(
            "parameter rows do not match exact probability query layout");
      }
      prob[static_cast<R_xlen_t>(query_idx)] = integrate_to_infinity(
          [&](const double rt) {
            const double log_density = exact_loglik_for_trial(
                plan,
                params,
                static_cast<int>(param_row),
                queries[query_idx].outcome_label,
                rt,
                -std::numeric_limits<double>::infinity());
            if (!std::isfinite(log_density)) {
              return 0.0;
            }
            const double density = std::exp(log_density);
            return std::isfinite(density) && density > 0.0 ? density : 0.0;
          });
      param_row += leaf_count;
    }
  }

  if (param_row != static_cast<std::size_t>(params.matrix.nrow())) {
    throw std::runtime_error(
        "parameter rows contain unused entries for the exact probability layout");
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = prob,
      Rcpp::Named("n_trials") = static_cast<int>(queries.size()),
      Rcpp::Named("n_variants") = static_cast<int>(plans.size()));
}

inline SEXP evaluate_exact_trials(const compile::CompiledModel &compiled,
                                  SEXP paramsSEXP,
                                  SEXP dataSEXP,
                                  const double min_ll) {
  std::unordered_map<std::string, semantic::Index> variant_index_by_component;
  std::vector<runtime::LoweredExactVariant> lowered_variants;
  lowered_variants.reserve(compiled.variants.size());

  for (const auto &variant : compiled.variants) {
    if (variant.backend != compile::BackendKind::Exact) {
      continue;
    }
    variant_index_by_component.emplace(
        variant.component_id,
        static_cast<semantic::Index>(lowered_variants.size()));
    lowered_variants.push_back(runtime::lower_exact_variant(variant));
  }
  if (lowered_variants.empty()) {
    throw std::runtime_error("exact evaluator found no exact variants");
  }

  ParamView params(paramsSEXP);
  Rcpp::DataFrame data(dataSEXP);
  const auto observations = collapse_exact_observations(data, variant_index_by_component);
  const auto blocks = build_variant_blocks(observations);

  std::vector<ExactVariantPlan> plans;
  plans.reserve(lowered_variants.size());
  for (const auto &variant : lowered_variants) {
    plans.push_back(make_exact_variant_plan(variant));
  }

  Rcpp::NumericVector loglik(observations.size(), min_ll);
  std::size_t param_row = 0;
  for (const auto &block : blocks) {
    const auto &plan = plans[static_cast<std::size_t>(block.variant_index)];
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    if (!plan.ranked_supported) {
      for (int i = 0; i < block.row_count; ++i) {
        const auto row = static_cast<std::size_t>(block.start_row + i);
        if (observations[row].ranks.size() > 1U) {
          throw std::runtime_error(
              "exact ranked kernel does not support latent prerequisite timings yet");
        }
      }
    }
    for (int i = 0; i < block.row_count; ++i) {
      const auto row = block.start_row + i;
      if (param_row + leaf_count > static_cast<std::size_t>(params.matrix.nrow())) {
        throw std::runtime_error(
            "parameter rows do not match exact trial layout");
      }
      loglik[row] = exact_ranked_loglik_for_trial(
          plan,
          params,
          static_cast<int>(param_row),
          observations[static_cast<std::size_t>(row)],
          min_ll);
      param_row += leaf_count;
    }
  }

  if (param_row != static_cast<std::size_t>(params.matrix.nrow())) {
    throw std::runtime_error(
        "parameter rows contain unused entries for the exact trial layout");
  }

  Rcpp::List block_list(blocks.size());
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    const auto &block = blocks[i];
    block_list[i] = Rcpp::List::create(
        Rcpp::Named("component_id") =
            lowered_variants[static_cast<std::size_t>(block.variant_index)].component_id,
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
