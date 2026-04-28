#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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
  MissingRt = 1,
  MissingAll = 2
};

enum class MissingAllProbabilityMode : std::uint8_t {
  FiniteComplement = 0,
  DirectNoResponse = 1
};

enum class MissingRtProbabilityMode : std::uint8_t {
  BranchIntegral = 0,
  FiniteComplementOfNoResponse = 1
};

enum class ObservationPlanValueKind : std::uint8_t {
  Probability = 0,
  Log = 1
};

enum class ObservationPlanOpKind : std::uint8_t {
  Constant = 0,
  LogDensity = 1,
  FiniteOutcomeProbability = 2,
  NoResponseProbability = 3,
  Complement = 4,
  WeightedSum = 5,
  Log = 6
};

struct ObservationIndexSpan {
  semantic::Index offset{0};
  semantic::Index size{0};
};

struct ObservationPlanOp {
  ObservationPlanOpKind kind{ObservationPlanOpKind::Constant};
  ObservationPlanValueKind value_kind{ObservationPlanValueKind::Probability};
  semantic::Index semantic_code{semantic::kInvalidIndex};
  compile::BackendKind backend{compile::BackendKind::Exact};
  double weight{1.0};
  double constant{0.0};
  ObservationIndexSpan children{};
};

struct ObservationProbabilityPlan {
  std::vector<ObservationPlanOp> ops;
  std::vector<semantic::Index> child_ops;
  semantic::Index root{semantic::kInvalidIndex};
  ObservationPlanValueKind value_kind{ObservationPlanValueKind::Probability};

  bool empty() const {
    return root == semantic::kInvalidIndex || ops.empty();
  }
};

struct ObservedBranch {
  semantic::Index semantic_code{semantic::kInvalidIndex};
  double weight{1.0};
  ObservedRtKind rt_kind{ObservedRtKind::Keep};
};

struct ComponentObservationPlan {
  bool present{false};
  compile::BackendKind observed_backend{compile::BackendKind::Exact};
  compile::BackendKind semantic_backend{compile::BackendKind::Exact};
  std::vector<std::vector<ObservedBranch>> keep_by_code;
  std::vector<std::vector<ObservedBranch>> missing_rt_by_code;
  std::vector<MissingRtProbabilityMode> missing_rt_mode_by_code;
  std::vector<ObservedBranch> finite_observed_branches;
  std::vector<ObservedBranch> missing_all_branches;
  MissingAllProbabilityMode missing_all_mode{
      MissingAllProbabilityMode::FiniteComplement};
  std::vector<ObservationProbabilityPlan> finite_log_plans_by_code;
  std::vector<ObservationProbabilityPlan> missing_rt_probability_plans_by_code;
  std::vector<ObservationProbabilityPlan> missing_rt_log_plans_by_code;
  ObservationProbabilityPlan missing_all_probability_plan;
  ObservationProbabilityPlan missing_all_log_plan;
};

struct ObservedRecord {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  compile::BackendKind backend{compile::BackendKind::Exact};
  semantic::Index semantic_code{semantic::kInvalidIndex};
  double observed_rt{NA_REAL};
  double weight{1.0};
  semantic::Index row_map_index{semantic::kInvalidIndex};
};

inline semantic::Index append_observation_plan_op(
    ObservationProbabilityPlan *plan,
    ObservationPlanOp op,
    const std::vector<semantic::Index> &children = {}) {
  if (plan == nullptr) {
    return semantic::kInvalidIndex;
  }
  op.children.offset =
      static_cast<semantic::Index>(plan->child_ops.size());
  op.children.size = static_cast<semantic::Index>(children.size());
  plan->child_ops.insert(plan->child_ops.end(), children.begin(), children.end());
  const auto op_index = static_cast<semantic::Index>(plan->ops.size());
  plan->ops.push_back(op);
  return op_index;
}

inline ObservationProbabilityPlan make_empty_observation_plan(
    const ObservationPlanValueKind value_kind) {
  ObservationProbabilityPlan plan;
  plan.value_kind = value_kind;
  return plan;
}

inline ObservationProbabilityPlan make_weighted_log_density_plan(
    const std::vector<ObservedBranch> &branches,
    const compile::BackendKind backend) {
  ObservationProbabilityPlan plan =
      make_empty_observation_plan(ObservationPlanValueKind::Log);
  std::vector<semantic::Index> children;
  children.reserve(branches.size());
  for (const auto &branch : branches) {
    if (!(std::isfinite(branch.weight) && branch.weight > 0.0)) {
      continue;
    }
    ObservationPlanOp op;
    op.kind = ObservationPlanOpKind::LogDensity;
    op.value_kind = ObservationPlanValueKind::Log;
    op.semantic_code = branch.semantic_code;
    op.backend = backend;
    op.weight = branch.weight;
    children.push_back(append_observation_plan_op(&plan, op));
  }
  if (children.empty()) {
    return plan;
  }
  ObservationPlanOp root;
  root.kind = ObservationPlanOpKind::WeightedSum;
  root.value_kind = ObservationPlanValueKind::Log;
  plan.root = append_observation_plan_op(&plan, root, children);
  return plan;
}

inline ObservationProbabilityPlan make_weighted_probability_plan(
    const std::vector<ObservedBranch> &branches,
    const compile::BackendKind backend) {
  ObservationProbabilityPlan plan =
      make_empty_observation_plan(ObservationPlanValueKind::Probability);
  std::vector<semantic::Index> children;
  children.reserve(branches.size());
  for (const auto &branch : branches) {
    if (!(std::isfinite(branch.weight) && branch.weight > 0.0)) {
      continue;
    }
    ObservationPlanOp op;
    op.kind = ObservationPlanOpKind::FiniteOutcomeProbability;
    op.value_kind = ObservationPlanValueKind::Probability;
    op.semantic_code = branch.semantic_code;
    op.backend = backend;
    op.weight = branch.weight;
    children.push_back(append_observation_plan_op(&plan, op));
  }
  if (children.empty()) {
    return plan;
  }
  ObservationPlanOp root;
  root.kind = ObservationPlanOpKind::WeightedSum;
  root.value_kind = ObservationPlanValueKind::Probability;
  plan.root = append_observation_plan_op(&plan, root, children);
  return plan;
}

inline ObservationProbabilityPlan make_no_response_probability_plan() {
  ObservationProbabilityPlan plan =
      make_empty_observation_plan(ObservationPlanValueKind::Probability);
  ObservationPlanOp root;
  root.kind = ObservationPlanOpKind::NoResponseProbability;
  root.value_kind = ObservationPlanValueKind::Probability;
  plan.root = append_observation_plan_op(&plan, root);
  return plan;
}

inline ObservationProbabilityPlan make_complement_probability_plan(
    const ObservationProbabilityPlan &inner) {
  if (inner.empty()) {
    ObservationProbabilityPlan one =
        make_empty_observation_plan(ObservationPlanValueKind::Probability);
    ObservationPlanOp root;
    root.kind = ObservationPlanOpKind::Constant;
    root.value_kind = ObservationPlanValueKind::Probability;
    root.constant = 1.0;
    one.root = append_observation_plan_op(&one, root);
    return one;
  }
  ObservationProbabilityPlan plan = inner;
  plan.value_kind = ObservationPlanValueKind::Probability;
  ObservationPlanOp root;
  root.kind = ObservationPlanOpKind::Complement;
  root.value_kind = ObservationPlanValueKind::Probability;
  plan.root = append_observation_plan_op(&plan, root, {inner.root});
  return plan;
}

inline ObservationProbabilityPlan wrap_observation_plan_log(
    const ObservationProbabilityPlan &inner) {
  if (inner.empty()) {
    return make_empty_observation_plan(ObservationPlanValueKind::Log);
  }
  ObservationProbabilityPlan plan = inner;
  plan.value_kind = ObservationPlanValueKind::Log;
  ObservationPlanOp root;
  root.kind = ObservationPlanOpKind::Log;
  root.value_kind = ObservationPlanValueKind::Log;
  plan.root = append_observation_plan_op(&plan, root, {inner.root});
  return plan;
}

inline ObservationProbabilityPlan make_missing_rt_probability_plan(
    const ComponentObservationPlan &component_plan,
    const std::size_t observed_code) {
  if (observed_code < component_plan.missing_rt_mode_by_code.size() &&
      component_plan.missing_rt_mode_by_code[observed_code] ==
          MissingRtProbabilityMode::FiniteComplementOfNoResponse) {
    return make_complement_probability_plan(make_no_response_probability_plan());
  }
  std::vector<ObservedBranch> branches;
  if (observed_code < component_plan.keep_by_code.size()) {
    branches.insert(
        branches.end(),
        component_plan.keep_by_code[observed_code].begin(),
        component_plan.keep_by_code[observed_code].end());
  }
  if (observed_code < component_plan.missing_rt_by_code.size()) {
    branches.insert(
        branches.end(),
        component_plan.missing_rt_by_code[observed_code].begin(),
        component_plan.missing_rt_by_code[observed_code].end());
  }
  return make_weighted_probability_plan(
      branches,
      component_plan.semantic_backend);
}

inline ObservationProbabilityPlan make_missing_all_probability_plan(
    const ComponentObservationPlan &component_plan) {
  if (component_plan.missing_all_mode ==
      MissingAllProbabilityMode::DirectNoResponse) {
    return make_no_response_probability_plan();
  }
  return make_complement_probability_plan(
      make_weighted_probability_plan(
          component_plan.finite_observed_branches,
          component_plan.semantic_backend));
}

inline void compile_component_observation_probability_plans(
    ComponentObservationPlan *component_plan) {
  if (component_plan == nullptr) {
    return;
  }
  const auto n_codes = component_plan->keep_by_code.size();
  component_plan->finite_log_plans_by_code.assign(
      n_codes,
      make_empty_observation_plan(ObservationPlanValueKind::Log));
  component_plan->missing_rt_probability_plans_by_code.assign(
      n_codes,
      make_empty_observation_plan(ObservationPlanValueKind::Probability));
  component_plan->missing_rt_log_plans_by_code.assign(
      n_codes,
      make_empty_observation_plan(ObservationPlanValueKind::Log));
  for (std::size_t observed_code = 1; observed_code < n_codes; ++observed_code) {
    component_plan->finite_log_plans_by_code[observed_code] =
        make_weighted_log_density_plan(
            component_plan->keep_by_code[observed_code],
            component_plan->semantic_backend);
    component_plan->missing_rt_probability_plans_by_code[observed_code] =
        make_missing_rt_probability_plan(*component_plan, observed_code);
    component_plan->missing_rt_log_plans_by_code[observed_code] =
        wrap_observation_plan_log(
            component_plan->missing_rt_probability_plans_by_code[observed_code]);
  }
  component_plan->missing_all_probability_plan =
      make_missing_all_probability_plan(*component_plan);
  component_plan->missing_all_log_plan =
      wrap_observation_plan_log(component_plan->missing_all_probability_plan);
}

inline void compile_observation_probability_plans(
    std::vector<ComponentObservationPlan> *plans) {
  if (plans == nullptr) {
    return;
  }
  for (auto &plan : *plans) {
    if (plan.present) {
      compile_component_observation_probability_plans(&plan);
    }
  }
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
  if (!options.containsElementNamed("component") ||
      Rf_isNull(options["component"])) {
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
  return !options.containsElementNamed("component") ||
         Rf_isNull(options["component"]);
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

inline void append_branch(std::vector<std::vector<ObservedBranch>> *buckets,
                          const semantic::Index observed_code,
                          const ObservedBranch &branch,
                          std::vector<ObservedBranch> *finite_observed_branches) {
  if (observed_code <= 0 ||
      observed_code >= static_cast<semantic::Index>(buckets->size())) {
    throw std::runtime_error("observation plan produced an invalid observed outcome code");
  }
  (*buckets)[static_cast<std::size_t>(observed_code)].push_back(branch);
  finite_observed_branches->push_back(branch);
}

inline std::vector<ComponentObservationPlan> build_component_observation_plans(
    const Rcpp::List &prep,
    const compile::CompiledModel &compiled,
    const std::unordered_map<std::string, semantic::Index> &component_code_by_id,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_component_codes,
    const std::size_t n_outcome_codes) {
  if (!prep.containsElementNamed("outcomes") || Rf_isNull(prep["outcomes"])) {
    throw std::runtime_error(
        "prep must contain outcomes for observed likelihood evaluation");
  }

  std::vector<ComponentObservationPlan> plans(n_component_codes + 1U);
  for (const auto &variant : compiled.variants) {
    const auto component_it = component_code_by_id.find(variant.component_id);
    if (component_it == component_code_by_id.end()) {
      throw std::runtime_error(
          "observation plan found no prepared component code for '" +
          variant.component_id + "'");
    }
    auto &plan = plans[static_cast<std::size_t>(component_it->second)];
    plan.present = true;
    plan.observed_backend = variant.backend;
    plan.semantic_backend = variant.semantic_backend;
    plan.keep_by_code.assign(n_outcome_codes + 1U, {});
    plan.missing_rt_by_code.assign(n_outcome_codes + 1U, {});
    plan.missing_rt_mode_by_code.assign(
        n_outcome_codes + 1U,
        MissingRtProbabilityMode::BranchIntegral);
  }

  Rcpp::List outcomes(prep["outcomes"]);
  const SEXP names_sexp = outcomes.names();
  const Rcpp::CharacterVector outcome_names =
      Rf_isNull(names_sexp) ? Rcpp::CharacterVector()
                            : Rcpp::CharacterVector(names_sexp);

  for (const auto &variant : compiled.variants) {
    const auto component_it = component_code_by_id.find(variant.component_id);
    auto &plan = plans[static_cast<std::size_t>(component_it->second)];

    for (const auto &semantic_outcome : variant.model.outcomes) {
      const auto &semantic_label = semantic_outcome.label;
      const auto outcome_idx = resolve_outcome_def_index(
          outcomes,
          outcome_names,
          semantic_label,
          variant.component_id);
      if (outcome_idx == semantic::kInvalidIndex) {
        throw std::runtime_error(
            "observation plan found no outcome definition for '" +
            semantic_label + "'");
      }

      const auto semantic_code_it = outcome_code_by_label.find(semantic_label);
      if (semantic_code_it == outcome_code_by_label.end()) {
        throw std::runtime_error(
            "observation plan found no prepared outcome code for '" +
            semantic_label + "'");
      }

      const Rcpp::List outcome(outcomes[static_cast<R_xlen_t>(outcome_idx)]);
      const Rcpp::List options =
          outcome.containsElementNamed("options") && !Rf_isNull(outcome["options"])
              ? Rcpp::List(outcome["options"])
              : Rcpp::List();

      std::vector<std::pair<semantic::Index, ObservedBranch>> branches;
      branches.push_back({semantic_code_it->second,
                          ObservedBranch{
                              semantic_code_it->second,
                              1.0,
                              ObservedRtKind::Keep}});

      if (options.containsElementNamed("guess") && !Rf_isNull(options["guess"])) {
        const Rcpp::List guess(options["guess"]);
        if (!guess.containsElementNamed("labels") ||
            !guess.containsElementNamed("weights") ||
            Rf_isNull(guess["labels"]) ||
            Rf_isNull(guess["weights"])) {
          throw std::runtime_error("guess option requires labels and weights");
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
          const auto observed_code_it =
              outcome_code_by_label.find(guess_labels[static_cast<std::size_t>(i)]);
          if (observed_code_it == outcome_code_by_label.end()) {
            throw std::runtime_error(
                "guess option refers to unknown observed label '" +
                guess_labels[static_cast<std::size_t>(i)] + "'");
          }
          branches.push_back(
              {observed_code_it->second,
               ObservedBranch{
                   semantic_code_it->second,
                   guess_weights[i],
                   rt_policy == "na" ? ObservedRtKind::Missing
                                      : ObservedRtKind::Keep}});
        }
      }

      bool map_missing = false;
      semantic::Index mapped_code = semantic::kInvalidIndex;
      if (options.containsElementNamed("map_outcome_to") &&
          !Rf_isNull(options["map_outcome_to"])) {
        if (TYPEOF(options["map_outcome_to"]) != STRSXP) {
          throw std::runtime_error("map_outcome_to must be character or NA");
        }
        Rcpp::CharacterVector target(options["map_outcome_to"]);
        if (target.size() > 0 && STRING_ELT(target, 0) == NA_STRING) {
          map_missing = true;
        } else if (target.size() > 0) {
          const auto mapped_label = Rcpp::as<std::string>(target[0]);
          const auto mapped_it = outcome_code_by_label.find(mapped_label);
          if (mapped_it == outcome_code_by_label.end()) {
            throw std::runtime_error(
                "map_outcome_to refers to unknown observed label '" +
                mapped_label + "'");
          }
          mapped_code = mapped_it->second;
        }
      }

      for (auto &branch : branches) {
        if (!std::isfinite(branch.second.weight) || !(branch.second.weight > 0.0)) {
          continue;
        }

        semantic::Index observed_code = branch.first;
        if (map_missing) {
          plan.missing_all_branches.push_back(branch.second);
          observed_code = semantic::kInvalidIndex;
          branch.second.rt_kind = ObservedRtKind::Missing;
        } else if (mapped_code != semantic::kInvalidIndex) {
          observed_code = mapped_code;
        }

        if (observed_code == semantic::kInvalidIndex) {
          continue;
        }

        if (branch.second.rt_kind == ObservedRtKind::Keep) {
          append_branch(
              &plan.keep_by_code,
              observed_code,
              branch.second,
              &plan.finite_observed_branches);
        } else {
          append_branch(
              &plan.missing_rt_by_code,
              observed_code,
              branch.second,
              &plan.finite_observed_branches);
        }
      }
    }
  }

  return plans;
}

} // namespace detail
} // namespace accumulatr::eval
