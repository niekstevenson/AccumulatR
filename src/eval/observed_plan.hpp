#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../compile/project_semantic.hpp"
#include "../semantic/model.hpp"

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
    plan.semantic_backend = variant.semantic_backend;
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

inline Rcpp::List to_r_list(const ObservedBranch &branch) {
  return Rcpp::List::create(
      Rcpp::Named("semantic_label") = branch.semantic_label,
      Rcpp::Named("weight") = branch.weight,
      Rcpp::Named("rt_kind") = static_cast<int>(branch.rt_kind));
}

inline ObservedBranch observed_branch_from_r(const Rcpp::List &branch) {
  ObservedBranch out;
  out.semantic_label = Rcpp::as<std::string>(branch["semantic_label"]);
  out.weight = Rcpp::as<double>(branch["weight"]);
  out.rt_kind = static_cast<ObservedRtKind>(Rcpp::as<int>(branch["rt_kind"]));
  return out;
}

inline Rcpp::List to_r_branch_map(
    const std::unordered_map<std::string, std::vector<ObservedBranch>> &branches) {
  Rcpp::List out(branches.size());
  Rcpp::CharacterVector names(branches.size());
  std::size_t i = 0;
  for (const auto &[label, values] : branches) {
    Rcpp::List bucket(values.size());
    for (std::size_t j = 0; j < values.size(); ++j) {
      bucket[static_cast<R_xlen_t>(j)] = to_r_list(values[j]);
    }
    out[static_cast<R_xlen_t>(i)] = bucket;
    names[static_cast<R_xlen_t>(i)] = label;
    ++i;
  }
  out.attr("names") = names;
  return out;
}

inline std::unordered_map<std::string, std::vector<ObservedBranch>>
branch_map_from_r(SEXP value) {
  std::unordered_map<std::string, std::vector<ObservedBranch>> out;
  if (Rf_isNull(value)) {
    return out;
  }
  Rcpp::List branches(value);
  const auto names = Rcpp::CharacterVector(branches.names());
  for (R_xlen_t i = 0; i < branches.size(); ++i) {
    const auto label = Rcpp::as<std::string>(names[i]);
    const Rcpp::List bucket(branches[i]);
    std::vector<ObservedBranch> values;
    values.reserve(bucket.size());
    for (R_xlen_t j = 0; j < bucket.size(); ++j) {
      values.push_back(observed_branch_from_r(Rcpp::List(bucket[j])));
    }
    out.emplace(label, std::move(values));
  }
  return out;
}

inline Rcpp::List to_r_list(const ComponentObservationPlan &plan) {
  Rcpp::List finite(plan.finite_observed_branches.size());
  for (std::size_t i = 0; i < plan.finite_observed_branches.size(); ++i) {
    finite[static_cast<R_xlen_t>(i)] = to_r_list(plan.finite_observed_branches[i]);
  }
  return Rcpp::List::create(
      Rcpp::Named("observed_backend") = static_cast<int>(plan.observed_backend),
      Rcpp::Named("semantic_backend") = static_cast<int>(plan.semantic_backend),
      Rcpp::Named("keep") = to_r_branch_map(plan.keep),
      Rcpp::Named("missing_rt") = to_r_branch_map(plan.missing_rt),
      Rcpp::Named("finite_observed_branches") = finite);
}

inline ComponentObservationPlan component_plan_from_r(const Rcpp::List &plan) {
  ComponentObservationPlan out;
  out.observed_backend =
      static_cast<compile::BackendKind>(Rcpp::as<int>(plan["observed_backend"]));
  out.semantic_backend =
      static_cast<compile::BackendKind>(Rcpp::as<int>(plan["semantic_backend"]));
  out.keep = branch_map_from_r(plan["keep"]);
  out.missing_rt = branch_map_from_r(plan["missing_rt"]);
  const Rcpp::List finite(plan["finite_observed_branches"]);
  out.finite_observed_branches.reserve(finite.size());
  for (R_xlen_t i = 0; i < finite.size(); ++i) {
    out.finite_observed_branches.push_back(
        observed_branch_from_r(Rcpp::List(finite[i])));
  }
  return out;
}

inline Rcpp::List to_r_observed_plan(
    const std::unordered_map<std::string, ComponentObservationPlan> &plans) {
  Rcpp::List out(plans.size());
  Rcpp::CharacterVector names(plans.size());
  std::size_t i = 0;
  for (const auto &[component_id, plan] : plans) {
    out[static_cast<R_xlen_t>(i)] = to_r_list(plan);
    names[static_cast<R_xlen_t>(i)] = component_id;
    ++i;
  }
  out.attr("names") = names;
  return out;
}

inline std::unordered_map<std::string, ComponentObservationPlan>
observed_plan_from_r(SEXP value) {
  std::unordered_map<std::string, ComponentObservationPlan> out;
  if (Rf_isNull(value)) {
    return out;
  }
  Rcpp::List plans(value);
  const auto names = Rcpp::CharacterVector(plans.names());
  for (R_xlen_t i = 0; i < plans.size(); ++i) {
    const auto component_id = Rcpp::as<std::string>(names[i]);
    out.emplace(component_id, component_plan_from_r(Rcpp::List(plans[i])));
  }
  return out;
}

} // namespace detail
} // namespace accumulatr::eval
