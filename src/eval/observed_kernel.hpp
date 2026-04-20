#pragma once

#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "direct_kernel.hpp"
#include "exact_kernel.hpp"

namespace accumulatr::eval {
namespace detail {

enum class ObservedRtKind : std::uint8_t {
  Keep = 0,
  Missing = 1
};

enum class ObservedTrialKind : std::uint8_t {
  Finite = 0,
  Ranked = 1,
  MissingRt = 2,
  MissingAll = 3,
  UnsupportedMissingLabel = 4
};

struct ObservedBranch {
  std::string semantic_label;
  double weight{1.0};
  ObservedRtKind rt_kind{ObservedRtKind::Keep};
};

struct ComponentObservationPlan {
  compile::BackendKind observed_backend{compile::BackendKind::Exact};
  compile::BackendKind semantic_backend{compile::BackendKind::Exact};
  std::unordered_map<std::string, std::vector<ObservedBranch>> keep;
  std::unordered_map<std::string, std::vector<ObservedBranch>> missing_rt;
  std::vector<ObservedBranch> finite_observed_branches;
};

struct ObservedTrialPlan {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index start_row{semantic::kInvalidIndex};
  semantic::Index end_row{semantic::kInvalidIndex};
  std::string component_id;
  std::string observed_label;
  double observed_rt{NA_REAL};
  ObservedTrialKind kind{ObservedTrialKind::Finite};
};

struct ObservedRecord {
  semantic::Index trial_index{semantic::kInvalidIndex};
  compile::BackendKind backend{compile::BackendKind::Exact};
  std::string semantic_label;
  double weight{1.0};
};

inline bool backend_is_semantic_direct(const compile::CompiledVariant &variant) {
  if (variant.backend == compile::BackendKind::Direct) {
    return true;
  }
  for (const auto &reason : variant.backend_reasons) {
    if (reason != "outcome remapping" && reason != "guess outcome") {
      return false;
    }
  }
  return true;
}

inline bool string_cell_is_na(const Rcpp::CharacterVector &column,
                              const R_xlen_t row) {
  return STRING_ELT(column, row) == NA_STRING;
}

inline bool string_vector_contains(const std::vector<std::string> &values,
                                   const std::string &needle) {
  return std::find(values.begin(), values.end(), needle) != values.end();
}

inline std::vector<std::string> as_string_vector_checked(SEXP value) {
  if (Rf_isNull(value)) {
    return {};
  }
  if (TYPEOF(value) != STRSXP) {
    throw std::runtime_error("observation options must use character vectors");
  }
  Rcpp::CharacterVector chr(value);
  std::vector<std::string> out;
  out.reserve(chr.size());
  for (R_xlen_t i = 0; i < chr.size(); ++i) {
    if (STRING_ELT(chr, i) == NA_STRING) {
      out.emplace_back();
    } else {
      out.push_back(Rcpp::as<std::string>(chr[i]));
    }
  }
  return out;
}

inline std::string as_optional_string(SEXP value) {
  if (Rf_isNull(value) || TYPEOF(value) != STRSXP) {
    return {};
  }
  Rcpp::CharacterVector chr(value);
  if (chr.size() == 0 || STRING_ELT(chr, 0) == NA_STRING) {
    return {};
  }
  return Rcpp::as<std::string>(chr[0]);
}

inline bool component_matches_outcome(SEXP options_sexp,
                                      const std::string &component_id) {
  if (Rf_isNull(options_sexp)) {
    return false;
  }
  Rcpp::List options(options_sexp);
  if (!options.containsElementNamed("component") || Rf_isNull(options["component"])) {
    return false;
  }
  const auto component_ids = as_string_vector_checked(options["component"]);
  return string_vector_contains(component_ids, component_id);
}

inline bool outcome_is_unrestricted(SEXP options_sexp) {
  if (Rf_isNull(options_sexp)) {
    return true;
  }
  Rcpp::List options(options_sexp);
  return !options.containsElementNamed("component") || Rf_isNull(options["component"]);
}

inline semantic::Index resolve_outcome_def_index(const Rcpp::List &outcomes,
                                                 const Rcpp::CharacterVector &outcome_names,
                                                 const std::string &semantic_label,
                                                 const std::string &component_id) {
  std::vector<semantic::Index> matches;
  for (R_xlen_t i = 0; i < outcomes.size(); ++i) {
    if (i < outcome_names.size() &&
        !string_cell_is_na(outcome_names, i) &&
        Rcpp::as<std::string>(outcome_names[i]) == semantic_label) {
      matches.push_back(static_cast<semantic::Index>(i));
    }
  }
  if (matches.empty()) {
    return semantic::kInvalidIndex;
  }
  for (const auto idx : matches) {
    const Rcpp::List outcome(outcomes[static_cast<R_xlen_t>(idx)]);
    const SEXP options_sexp =
        outcome.containsElementNamed("options") ? outcome["options"] : R_NilValue;
    if (component_matches_outcome(options_sexp, component_id)) {
      return idx;
    }
  }
  for (const auto idx : matches) {
    const Rcpp::List outcome(outcomes[static_cast<R_xlen_t>(idx)]);
    const SEXP options_sexp =
        outcome.containsElementNamed("options") ? outcome["options"] : R_NilValue;
    if (outcome_is_unrestricted(options_sexp)) {
      return idx;
    }
  }
  return semantic::kInvalidIndex;
}

inline void append_branch(
    std::unordered_map<std::string, std::vector<ObservedBranch>> *map,
    const std::string &observed_label,
    const ObservedBranch &branch,
    std::vector<ObservedBranch> *finite_observed_branches) {
  auto &bucket = (*map)[observed_label];
  bucket.push_back(branch);
  finite_observed_branches->push_back(branch);
}

inline std::unordered_map<std::string, ComponentObservationPlan>
build_component_observation_plans(const Rcpp::List &prep,
                                  const compile::CompiledModel &compiled) {
  if (!prep.containsElementNamed("outcomes") || Rf_isNull(prep["outcomes"])) {
    throw std::runtime_error("prep must contain outcomes for observed likelihood evaluation");
  }

  std::unordered_map<std::string, ComponentObservationPlan> plans;
  plans.reserve(compiled.variants.size());
  for (const auto &variant : compiled.variants) {
    ComponentObservationPlan plan;
    plan.observed_backend = variant.backend;
    plan.semantic_backend =
        backend_is_semantic_direct(variant) ? compile::BackendKind::Direct
                                            : compile::BackendKind::Exact;
    plans.emplace(variant.component_id, std::move(plan));
  }

  Rcpp::List outcomes(prep["outcomes"]);
  const SEXP names_sexp = outcomes.names();
  const Rcpp::CharacterVector outcome_names =
      Rf_isNull(names_sexp) ? Rcpp::CharacterVector()
                            : Rcpp::CharacterVector(names_sexp);

  std::vector<std::string> semantic_labels;
  semantic_labels.reserve(outcomes.size());
  for (R_xlen_t i = 0; i < outcome_names.size(); ++i) {
    if (string_cell_is_na(outcome_names, i)) {
      continue;
    }
    const auto label = Rcpp::as<std::string>(outcome_names[i]);
    if (!string_vector_contains(semantic_labels, label)) {
      semantic_labels.push_back(label);
    }
  }

  for (auto &entry : plans) {
    const auto &component_id = entry.first;
    auto &plan = entry.second;
    const auto variant_it = std::find_if(
        compiled.variants.begin(),
        compiled.variants.end(),
        [&component_id](const compile::CompiledVariant &variant) {
          return variant.component_id == component_id;
        });
    if (variant_it == compiled.variants.end()) {
      throw std::runtime_error(
          "observation plan found no compiled variant for component '" +
          component_id + "'");
    }
    std::vector<std::string> surviving_labels;
    surviving_labels.reserve(variant_it->model.outcomes.size());
    for (const auto &outcome : variant_it->model.outcomes) {
      surviving_labels.push_back(outcome.label);
    }

    for (const auto &semantic_label : semantic_labels) {
      if (!string_vector_contains(surviving_labels, semantic_label)) {
        continue;
      }
      const auto outcome_idx =
          resolve_outcome_def_index(outcomes, outcome_names, semantic_label, component_id);
      if (outcome_idx == semantic::kInvalidIndex) {
        continue;
      }

      const Rcpp::List outcome(outcomes[static_cast<R_xlen_t>(outcome_idx)]);
      const Rcpp::List options =
          outcome.containsElementNamed("options") && !Rf_isNull(outcome["options"])
              ? Rcpp::List(outcome["options"])
              : Rcpp::List();

      std::vector<std::pair<std::string, ObservedBranch>> branches;
      branches.push_back({semantic_label,
                          ObservedBranch{semantic_label, 1.0, ObservedRtKind::Keep}});

      if (options.containsElementNamed("guess") && !Rf_isNull(options["guess"])) {
        const Rcpp::List guess(options["guess"]);
        if (!guess.containsElementNamed("labels") ||
            !guess.containsElementNamed("weights") ||
            Rf_isNull(guess["labels"]) ||
            Rf_isNull(guess["weights"])) {
          throw std::runtime_error(
              "guess option requires labels and weights");
        }
        const auto guess_labels = as_string_vector_checked(guess["labels"]);
        const Rcpp::NumericVector guess_weights(guess["weights"]);
        const std::string rt_policy =
            guess.containsElementNamed("rt_policy") && !Rf_isNull(guess["rt_policy"])
                ? as_optional_string(guess["rt_policy"])
                : "keep";
        if (guess_labels.size() != static_cast<std::size_t>(guess_weights.size())) {
          throw std::runtime_error(
              "guess option requires labels and weights of equal length");
        }
        if (rt_policy != "keep" && rt_policy != "na") {
          throw std::runtime_error("guess rt_policy must be 'keep' or 'na'");
        }
        branches.clear();
        branches.reserve(guess_labels.size());
        for (R_xlen_t i = 0; i < guess_weights.size(); ++i) {
          branches.push_back(
              {guess_labels[static_cast<std::size_t>(i)],
               ObservedBranch{
                   semantic_label,
                   guess_weights[i],
                   rt_policy == "na" ? ObservedRtKind::Missing
                                      : ObservedRtKind::Keep}});
        }
      }

      bool map_missing = false;
      std::string mapped_label;
      if (options.containsElementNamed("map_outcome_to") &&
          !Rf_isNull(options["map_outcome_to"])) {
        if (TYPEOF(options["map_outcome_to"]) != STRSXP) {
          throw std::runtime_error("map_outcome_to must be character or NA");
        }
        Rcpp::CharacterVector target(options["map_outcome_to"]);
        if (target.size() > 0 && STRING_ELT(target, 0) == NA_STRING) {
          map_missing = true;
        } else if (target.size() > 0) {
          mapped_label = Rcpp::as<std::string>(target[0]);
        }
      }

      for (auto &branch : branches) {
        const double weight = branch.second.weight;
        if (!std::isfinite(weight) || !(weight > 0.0)) {
          continue;
        }

        std::string observed_label = branch.first;
        if (map_missing) {
          observed_label.clear();
          branch.second.rt_kind = ObservedRtKind::Missing;
        } else if (!mapped_label.empty()) {
          observed_label = mapped_label;
        }

        if (observed_label.empty()) {
          continue;
        }

        if (branch.second.rt_kind == ObservedRtKind::Keep) {
          append_branch(
              &plan.keep,
              observed_label,
              branch.second,
              &plan.finite_observed_branches);
        } else {
          append_branch(
              &plan.missing_rt,
              observed_label,
              branch.second,
              &plan.finite_observed_branches);
        }
      }
    }
  }

  return plans;
}

inline std::string read_string_cell(const Rcpp::CharacterVector &column,
                                    const R_xlen_t row) {
  if (string_cell_is_na(column, row)) {
    return {};
  }
  return Rcpp::as<std::string>(column[row]);
}

inline std::vector<ObservedTrialPlan> build_observed_trial_plan(
    const Rcpp::DataFrame &data,
    const std::unordered_map<std::string, ComponentObservationPlan> &component_plans) {
  const auto n_rows = data.nrows();
  if (n_rows == 0) {
    return {};
  }
  if (!data.containsElementNamed("R") || !data.containsElementNamed("rt")) {
    throw std::runtime_error("observed likelihood evaluator requires columns R and rt");
  }

  int max_rank = 1;
  for (int rank = 2;; ++rank) {
    const std::string r_col = "R" + std::to_string(rank);
    const std::string rt_col = "rt" + std::to_string(rank);
    const bool has_r = data.containsElementNamed(r_col.c_str());
    const bool has_rt = data.containsElementNamed(rt_col.c_str());
    if (has_r != has_rt) {
      throw std::runtime_error(
          "ranked columns must appear as matched Rk/rtk pairs");
    }
    if (!has_r) {
      break;
    }
    max_rank = rank;
  }

  const auto trial =
      data.containsElementNamed("trial")
          ? Rcpp::as<Rcpp::IntegerVector>(data["trial"])
          : Rcpp::seq_len(n_rows);
  const auto label = Rcpp::as<Rcpp::CharacterVector>(data["R"]);
  const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  const bool has_component = data.containsElementNamed("component");
  const auto component =
      has_component ? Rcpp::as<Rcpp::CharacterVector>(data["component"])
                    : Rcpp::CharacterVector(n_rows, NA_STRING);

  std::vector<Rcpp::CharacterVector> rank_labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> rank_times(static_cast<std::size_t>(max_rank + 1));
  for (int rank = 2; rank <= max_rank; ++rank) {
    rank_labels[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::CharacterVector>(data[("R" + std::to_string(rank)).c_str()]);
    rank_times[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::NumericVector>(data[("rt" + std::to_string(rank)).c_str()]);
  }

  std::vector<ObservedTrialPlan> out;
  int last_trial = NA_INTEGER;
  for (R_xlen_t i = 0; i < n_rows; ++i) {
    if (i > 0 && trial[i] == last_trial) {
      continue;
    }
    const auto start = static_cast<semantic::Index>(i);
    semantic::Index end = start;
    while (static_cast<R_xlen_t>(end + 1) < n_rows && trial[end + 1] == trial[i]) {
      ++end;
    }
    last_trial = trial[i];

    std::string component_id = "__default__";
    if (has_component && !string_cell_is_na(component, i)) {
      component_id = Rcpp::as<std::string>(component[i]);
    } else if (component_plans.size() == 1U) {
      component_id = component_plans.begin()->first;
    }
    const auto component_it = component_plans.find(component_id);
    if (component_it == component_plans.end()) {
      throw std::runtime_error(
          "observed likelihood evaluator found no variant for component '" +
          component_id + "'");
    }

    const std::string observed_label = read_string_cell(label, i);
    const double observed_rt = rt[i];
    bool has_ranked = false;
    for (int rank = 2; rank <= max_rank; ++rank) {
      const bool label_missing =
          STRING_ELT(rank_labels[static_cast<std::size_t>(rank)], i) == NA_STRING;
      const bool rt_missing =
          Rcpp::NumericVector::is_na(rank_times[static_cast<std::size_t>(rank)][i]);
      if (!label_missing || !rt_missing) {
        has_ranked = true;
        break;
      }
    }

    ObservedTrialPlan plan;
    plan.trial_index = static_cast<semantic::Index>(out.size());
    plan.start_row = start;
    plan.end_row = end;
    plan.component_id = component_id;
    plan.observed_label = observed_label;
    plan.observed_rt = observed_rt;

    if (has_ranked) {
      if (component_it->second.observed_backend != compile::BackendKind::Exact) {
        throw std::runtime_error("ranked observations require the exact backend");
      }
      plan.kind = ObservedTrialKind::Ranked;
    } else if (observed_label.empty() && Rcpp::NumericVector::is_na(observed_rt)) {
      plan.kind = ObservedTrialKind::MissingAll;
    } else if (!observed_label.empty() && Rcpp::NumericVector::is_na(observed_rt)) {
      plan.kind = ObservedTrialKind::MissingRt;
    } else if (observed_label.empty()) {
      plan.kind = ObservedTrialKind::UnsupportedMissingLabel;
    } else {
      plan.kind = ObservedTrialKind::Finite;
    }

    out.push_back(std::move(plan));
  }
  return out;
}

inline Rcpp::NumericMatrix subset_param_rows(const Rcpp::NumericMatrix &params,
                                             const std::vector<int> &rows) {
  Rcpp::NumericMatrix out(static_cast<R_xlen_t>(rows.size()), params.ncol());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    for (int j = 0; j < params.ncol(); ++j) {
      out(static_cast<R_xlen_t>(i), j) = params(rows[i], j);
    }
  }
  out.attr("dimnames") = params.attr("dimnames");
  return out;
}

inline Rcpp::DataFrame copy_trial_blocks(
    const Rcpp::DataFrame &data,
    const std::vector<ObservedTrialPlan> &trials,
    const std::vector<semantic::Index> &trial_indices,
    const std::vector<std::string> *replacement_labels,
    Rcpp::NumericMatrix *out_params,
    const Rcpp::NumericMatrix &params) {
  int max_rank = 1;
  for (int rank = 2;; ++rank) {
    const std::string r_col = "R" + std::to_string(rank);
    const std::string rt_col = "rt" + std::to_string(rank);
    const bool has_r = data.containsElementNamed(r_col.c_str());
    const bool has_rt = data.containsElementNamed(rt_col.c_str());
    if (!has_r || !has_rt) {
      break;
    }
    max_rank = rank;
  }

  const auto source_trial =
      data.containsElementNamed("trial")
          ? Rcpp::as<Rcpp::IntegerVector>(data["trial"])
          : Rcpp::seq_len(data.nrows());
  const auto source_r = Rcpp::as<Rcpp::CharacterVector>(data["R"]);
  const auto source_rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  const bool has_component = data.containsElementNamed("component");
  const auto source_component =
      has_component ? Rcpp::as<Rcpp::CharacterVector>(data["component"])
                    : Rcpp::CharacterVector(data.nrows(), NA_STRING);
  std::vector<Rcpp::CharacterVector> rank_labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> rank_times(static_cast<std::size_t>(max_rank + 1));
  for (int rank = 2; rank <= max_rank; ++rank) {
    rank_labels[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::CharacterVector>(data[("R" + std::to_string(rank)).c_str()]);
    rank_times[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::NumericVector>(data[("rt" + std::to_string(rank)).c_str()]);
  }

  std::vector<int> row_indices;
  for (const auto trial_index : trial_indices) {
    const auto &trial = trials[static_cast<std::size_t>(trial_index)];
    for (semantic::Index row = trial.start_row; row <= trial.end_row; ++row) {
      row_indices.push_back(static_cast<int>(row));
    }
  }
  *out_params = subset_param_rows(params, row_indices);

  Rcpp::IntegerVector new_trial(row_indices.size());
  Rcpp::CharacterVector new_r(row_indices.size(), NA_STRING);
  Rcpp::NumericVector new_rt(row_indices.size(), NA_REAL);
  Rcpp::CharacterVector new_component(row_indices.size(), NA_STRING);
  std::vector<Rcpp::CharacterVector> new_rank_labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> new_rank_times(static_cast<std::size_t>(max_rank + 1));
  for (int rank = 2; rank <= max_rank; ++rank) {
    new_rank_labels[static_cast<std::size_t>(rank)] =
        Rcpp::CharacterVector(row_indices.size(), NA_STRING);
    new_rank_times[static_cast<std::size_t>(rank)] =
        Rcpp::NumericVector(row_indices.size(), NA_REAL);
  }

  std::size_t out_row = 0;
  for (std::size_t i = 0; i < trial_indices.size(); ++i) {
    const auto &trial = trials[static_cast<std::size_t>(trial_indices[i])];
    for (semantic::Index row = trial.start_row; row <= trial.end_row; ++row, ++out_row) {
      new_trial[static_cast<R_xlen_t>(out_row)] = static_cast<int>(i + 1U);
      new_r[static_cast<R_xlen_t>(out_row)] = source_r[static_cast<R_xlen_t>(row)];
      new_rt[static_cast<R_xlen_t>(out_row)] = source_rt[static_cast<R_xlen_t>(row)];
      if (has_component) {
        new_component[static_cast<R_xlen_t>(out_row)] =
            source_component[static_cast<R_xlen_t>(row)];
      }
      for (int rank = 2; rank <= max_rank; ++rank) {
        new_rank_labels[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(out_row)] =
            rank_labels[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(row)];
        new_rank_times[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(out_row)] =
            rank_times[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(row)];
      }
    }
    if (replacement_labels != nullptr) {
      const auto replacement = (*replacement_labels)[i];
      out_row -= static_cast<std::size_t>(trial.end_row - trial.start_row + 1);
      new_r[static_cast<R_xlen_t>(out_row)] = replacement;
      out_row += static_cast<std::size_t>(trial.end_row - trial.start_row + 1);
    }
  }

  Rcpp::List out;
  Rcpp::CharacterVector names;
  out.push_back(new_trial, "trial");
  out.push_back(new_r, "R");
  out.push_back(new_rt, "rt");
  if (has_component) {
    out.push_back(new_component, "component");
  }
  for (int rank = 2; rank <= max_rank; ++rank) {
    out.push_back(new_rank_labels[static_cast<std::size_t>(rank)],
                  ("R" + std::to_string(rank)).c_str());
    out.push_back(new_rank_times[static_cast<std::size_t>(rank)],
                  ("rt" + std::to_string(rank)).c_str());
  }
  out.attr("class") = "data.frame";
  out.attr("row.names") =
      Rcpp::IntegerVector::create(NA_INTEGER, -static_cast<int>(row_indices.size()));
  return Rcpp::DataFrame(out);
}

inline Rcpp::NumericVector evaluate_loglik_records(
    const compile::CompiledModel &compiled,
    const semantic::SemanticModel &model,
    const std::vector<ObservedTrialPlan> &trials,
    const Rcpp::DataFrame &data,
    const Rcpp::NumericMatrix &params,
    const std::vector<ObservedRecord> &records,
    const double min_ll) {
  Rcpp::NumericVector out(records.size(), R_NegInf);
  if (records.empty()) {
    return out;
  }

  for (const auto backend : {compile::BackendKind::Direct, compile::BackendKind::Exact}) {
    std::vector<semantic::Index> idxs;
    idxs.reserve(records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
      if (records[i].backend == backend) {
        idxs.push_back(static_cast<semantic::Index>(i));
      }
    }
    if (idxs.empty()) {
      continue;
    }

    std::vector<semantic::Index> trial_indices;
    std::vector<std::string> labels;
    trial_indices.reserve(idxs.size());
    labels.reserve(idxs.size());
    for (const auto idx : idxs) {
      trial_indices.push_back(records[static_cast<std::size_t>(idx)].trial_index);
      labels.push_back(records[static_cast<std::size_t>(idx)].semantic_label);
    }

    Rcpp::NumericMatrix subset_params;
    const auto subset_data =
        copy_trial_blocks(data, trials, trial_indices, &labels, &subset_params, params);
    SEXP eval_sexp = R_NilValue;
    if (backend == compile::BackendKind::Direct) {
      eval_sexp = evaluate_direct_trials(
          compiled, model, subset_params, subset_data, min_ll);
    } else {
      eval_sexp = evaluate_exact_trials(
          compiled, subset_params, subset_data, min_ll);
    }
    Rcpp::List eval(eval_sexp);
    const auto loglik = Rcpp::as<Rcpp::NumericVector>(eval["loglik"]);
    for (std::size_t i = 0; i < idxs.size(); ++i) {
      out[static_cast<R_xlen_t>(idxs[i])] = loglik[static_cast<R_xlen_t>(i)];
    }
  }

  return out;
}

inline Rcpp::NumericVector evaluate_probability_records(
    const compile::CompiledModel &compiled,
    const semantic::SemanticModel &model,
    const std::vector<ObservedTrialPlan> &trials,
    const Rcpp::DataFrame &data,
    const Rcpp::NumericMatrix &params,
    const std::vector<ObservedRecord> &records) {
  Rcpp::NumericVector out(records.size(), 0.0);
  if (records.empty()) {
    return out;
  }

  for (const auto backend : {compile::BackendKind::Direct, compile::BackendKind::Exact}) {
    std::vector<semantic::Index> idxs;
    idxs.reserve(records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
      if (records[i].backend == backend) {
        idxs.push_back(static_cast<semantic::Index>(i));
      }
    }
    if (idxs.empty()) {
      continue;
    }

    std::vector<semantic::Index> trial_indices;
    std::vector<std::string> labels;
    trial_indices.reserve(idxs.size());
    labels.reserve(idxs.size());
    for (const auto idx : idxs) {
      trial_indices.push_back(records[static_cast<std::size_t>(idx)].trial_index);
      labels.push_back(records[static_cast<std::size_t>(idx)].semantic_label);
    }

    Rcpp::NumericMatrix subset_params;
    const auto subset_data =
        copy_trial_blocks(data, trials, trial_indices, &labels, &subset_params, params);
    SEXP eval_sexp = R_NilValue;
    if (backend == compile::BackendKind::Direct) {
      eval_sexp = evaluate_direct_outcome_probabilities(
          compiled, model, subset_params, subset_data);
    } else {
      eval_sexp = evaluate_exact_outcome_probabilities(
          compiled, subset_params, subset_data);
    }
    Rcpp::List eval(eval_sexp);
    const auto prob = Rcpp::as<Rcpp::NumericVector>(eval["probability"]);
    for (std::size_t i = 0; i < idxs.size(); ++i) {
      out[static_cast<R_xlen_t>(idxs[i])] = prob[static_cast<R_xlen_t>(i)];
    }
  }

  return out;
}

inline double logsumexp_records(const std::vector<double> &values) {
  if (values.empty()) {
    return R_NegInf;
  }
  double anchor = R_NegInf;
  for (const auto value : values) {
    if (std::isfinite(value) && value > anchor) {
      anchor = value;
    }
  }
  if (!std::isfinite(anchor)) {
    return R_NegInf;
  }
  double sum = 0.0;
  for (const auto value : values) {
    if (std::isfinite(value)) {
      sum += std::exp(value - anchor);
    }
  }
  return sum > 0.0 ? anchor + std::log(sum) : R_NegInf;
}

inline SEXP evaluate_observed_trials(const Rcpp::List &prep,
                                     const semantic::SemanticModel &model,
                                     const compile::CompiledModel &compiled,
                                     SEXP paramsSEXP,
                                     SEXP dataSEXP,
                                     const double min_ll) {
  Rcpp::DataFrame data(dataSEXP);
  Rcpp::NumericMatrix params(Rcpp::as<Rcpp::NumericMatrix>(paramsSEXP));
  const auto component_plans = build_component_observation_plans(prep, compiled);
  const auto trials = build_observed_trial_plan(data, component_plans);

  Rcpp::NumericVector loglik(trials.size(), min_ll);

  std::vector<semantic::Index> ranked_trial_indices;
  ranked_trial_indices.reserve(trials.size());
  std::vector<ObservedRecord> finite_records;
  std::vector<ObservedRecord> missing_rt_records;
  std::vector<ObservedRecord> missing_all_records;

  for (const auto &trial : trials) {
    const auto component_it = component_plans.find(trial.component_id);
    if (component_it == component_plans.end()) {
      throw std::runtime_error(
          "observed likelihood evaluator found no observation plan for component '" +
          trial.component_id + "'");
    }
    const auto &component_plan = component_it->second;

    switch (trial.kind) {
    case ObservedTrialKind::Ranked:
      ranked_trial_indices.push_back(trial.trial_index);
      break;
    case ObservedTrialKind::Finite: {
      const auto keep_it = component_plan.keep.find(trial.observed_label);
      if (keep_it != component_plan.keep.end()) {
        for (const auto &branch : keep_it->second) {
          finite_records.push_back(ObservedRecord{
              trial.trial_index,
              component_plan.semantic_backend,
              branch.semantic_label,
              branch.weight});
        }
      }
      break;
    }
    case ObservedTrialKind::MissingRt: {
      const auto missing_it = component_plan.missing_rt.find(trial.observed_label);
      if (missing_it != component_plan.missing_rt.end()) {
        for (const auto &branch : missing_it->second) {
          missing_rt_records.push_back(ObservedRecord{
              trial.trial_index,
              component_plan.semantic_backend,
              branch.semantic_label,
              branch.weight});
        }
      }
      break;
    }
    case ObservedTrialKind::MissingAll:
      for (const auto &branch : component_plan.finite_observed_branches) {
        missing_all_records.push_back(ObservedRecord{
            trial.trial_index,
            component_plan.semantic_backend,
            branch.semantic_label,
            branch.weight});
      }
      break;
    case ObservedTrialKind::UnsupportedMissingLabel:
      break;
    }
  }

  if (!ranked_trial_indices.empty()) {
    Rcpp::NumericMatrix subset_params;
    const auto subset_data =
        copy_trial_blocks(data, trials, ranked_trial_indices, nullptr, &subset_params, params);
    Rcpp::List eval(
        evaluate_exact_trials(compiled, subset_params, subset_data, min_ll));
    const auto ranked_loglik = Rcpp::as<Rcpp::NumericVector>(eval["loglik"]);
    for (std::size_t i = 0; i < ranked_trial_indices.size(); ++i) {
      loglik[static_cast<R_xlen_t>(ranked_trial_indices[i])] =
          ranked_loglik[static_cast<R_xlen_t>(i)];
    }
  }

  if (!finite_records.empty()) {
    const auto donor_loglik = evaluate_loglik_records(
        compiled, model, trials, data, params, finite_records, min_ll);
    for (const auto &trial : trials) {
      if (trial.kind != ObservedTrialKind::Finite) {
        continue;
      }
      std::vector<double> values;
      for (std::size_t i = 0; i < finite_records.size(); ++i) {
        if (finite_records[i].trial_index != trial.trial_index) {
          continue;
        }
        if (!(finite_records[i].weight > 0.0)) {
          continue;
        }
        values.push_back(std::log(finite_records[i].weight) + donor_loglik[static_cast<R_xlen_t>(i)]);
      }
      const auto value = logsumexp_records(values);
      loglik[static_cast<R_xlen_t>(trial.trial_index)] =
          std::isfinite(value) ? value : min_ll;
    }
  }

  if (!missing_rt_records.empty()) {
    const auto donor_prob = evaluate_probability_records(
        compiled, model, trials, data, params, missing_rt_records);
    for (const auto &trial : trials) {
      if (trial.kind != ObservedTrialKind::MissingRt) {
        continue;
      }
      double total = 0.0;
      for (std::size_t i = 0; i < missing_rt_records.size(); ++i) {
        if (missing_rt_records[i].trial_index != trial.trial_index) {
          continue;
        }
        total += missing_rt_records[i].weight * donor_prob[static_cast<R_xlen_t>(i)];
      }
      loglik[static_cast<R_xlen_t>(trial.trial_index)] =
          std::isfinite(total) && total > 0.0 ? std::log(total) : min_ll;
    }
  }

  if (!missing_all_records.empty()) {
    const auto finite_prob = evaluate_probability_records(
        compiled, model, trials, data, params, missing_all_records);
    for (const auto &trial : trials) {
      if (trial.kind != ObservedTrialKind::MissingAll) {
        continue;
      }
      double total_finite = 0.0;
      for (std::size_t i = 0; i < missing_all_records.size(); ++i) {
        if (missing_all_records[i].trial_index != trial.trial_index) {
          continue;
        }
        total_finite += missing_all_records[i].weight * finite_prob[static_cast<R_xlen_t>(i)];
      }
      const double missing_prob =
          std::max(0.0, 1.0 - (std::isfinite(total_finite) ? total_finite : 1.0));
      loglik[static_cast<R_xlen_t>(trial.trial_index)] =
          missing_prob > 0.0 ? std::log(missing_prob) : min_ll;
    }
  } else {
    for (const auto &trial : trials) {
      if (trial.kind == ObservedTrialKind::MissingAll) {
        loglik[static_cast<R_xlen_t>(trial.trial_index)] = 0.0;
      }
    }
  }

  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    if (!std::isfinite(loglik[i])) {
      loglik[i] = min_ll;
    }
  }

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    total_loglik += static_cast<double>(loglik[i]);
  }

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("n_trials") = static_cast<int>(trials.size()));
}

} // namespace detail
} // namespace accumulatr::eval
