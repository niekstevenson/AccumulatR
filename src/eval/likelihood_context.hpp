#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../compile/prep_to_semantic.hpp"
#include "../compile/project_semantic.hpp"
#include "../semantic/model.hpp"
#include "exact_sequence.hpp"
#include "observed_plan.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct NativeLikelihoodContext {
  semantic::SemanticModel model;
  std::vector<ComponentObservationPlan> observed_plans_by_component_code;
  bool observed_identity{false};
  compile::BackendKind identity_backend{compile::BackendKind::Exact};
  bool ranked_supported{true};
  std::vector<semantic::Index> exact_variant_index_by_component_code;
  std::vector<ExactVariantPlan> exact_plans;
};

inline void compile_component_weight_parameter_layout(
    semantic::SemanticModel *model) {
  if (model == nullptr) {
    return;
  }
  std::unordered_map<std::string, int> index_by_name;
  for (auto &component : model->components) {
    component.weight_param_index = -1;
    if (component.weight_name.empty()) {
      continue;
    }
    auto it = index_by_name.find(component.weight_name);
    if (it == index_by_name.end()) {
      const int next = static_cast<int>(index_by_name.size());
      it = index_by_name.emplace(component.weight_name, next).first;
    }
    component.weight_param_index = it->second;
  }
  model->component_weight_param_count =
      static_cast<int>(index_by_name.size());
}

inline void compile_exact_plan_leaf_row_offsets(
    const semantic::SemanticModel &model,
    std::vector<ExactVariantPlan> *plans) {
  if (plans == nullptr) {
    return;
  }
  std::unordered_map<std::string, int> global_leaf_index;
  global_leaf_index.reserve(model.leaves.size());
  for (std::size_t i = 0; i < model.leaves.size(); ++i) {
    global_leaf_index.emplace(model.leaves[i].id, static_cast<int>(i));
  }
  for (auto &plan : *plans) {
    plan.leaf_row_offsets.clear();
    plan.leaf_row_offsets.reserve(plan.lowered.leaf_ids.size());
    for (const auto &leaf_id : plan.lowered.leaf_ids) {
      plan.leaf_row_offsets.push_back(global_leaf_index.at(leaf_id));
    }
  }
}

inline bool observation_plans_are_identity(
    const std::vector<ComponentObservationPlan> &plans,
    compile::BackendKind *identity_backend = nullptr) {
  compile::BackendKind backend = compile::BackendKind::Exact;
  bool have_backend = false;
  for (std::size_t component_code = 1; component_code < plans.size(); ++component_code) {
    const auto &plan = plans[component_code];
    if (!plan.present) {
      continue;
    }
    if (plan.observed_backend != plan.semantic_backend) {
      return false;
    }
    for (std::size_t outcome_code = 1; outcome_code < plan.missing_rt_by_code.size();
         ++outcome_code) {
      if (!plan.missing_rt_by_code[outcome_code].empty()) {
        return false;
      }
    }
    for (std::size_t outcome_code = 1; outcome_code < plan.keep_by_code.size();
         ++outcome_code) {
      const auto &branches = plan.keep_by_code[outcome_code];
      if (branches.empty()) {
        continue;
      }
      if (branches.size() != 1U) {
        return false;
      }
      const auto &branch = branches.front();
      if (branch.semantic_code != static_cast<semantic::Index>(outcome_code) ||
          branch.rt_kind != ObservedRtKind::Keep ||
          std::fabs(branch.weight - 1.0) > 1e-12) {
        return false;
      }
    }
    if (!have_backend) {
      backend = plan.semantic_backend;
      have_backend = true;
    } else if (backend != plan.semantic_backend) {
      return false;
    }
  }
  if (identity_backend != nullptr) {
    *identity_backend = backend;
  }
  return true;
}

inline bool observed_branches_cover_exact_outcomes(
    const ExactVariantPlan &exact_plan,
    const std::vector<ObservedBranch> &branches) {
  if (exact_plan.outcomes.empty()) {
    return false;
  }

  std::vector<double> covered_weight(exact_plan.outcomes.size(), 0.0);
  for (const auto &branch : branches) {
    if (!(std::isfinite(branch.weight) && branch.weight > 0.0) ||
        branch.semantic_code == semantic::kInvalidIndex ||
        static_cast<std::size_t>(branch.semantic_code) >=
            exact_plan.outcome_index_by_code.size()) {
      return false;
    }
    const auto outcome_index =
        exact_plan.outcome_index_by_code[
            static_cast<std::size_t>(branch.semantic_code)];
    if (outcome_index == semantic::kInvalidIndex ||
        static_cast<std::size_t>(outcome_index) >= covered_weight.size()) {
      return false;
    }
    covered_weight[static_cast<std::size_t>(outcome_index)] += branch.weight;
  }

  for (const auto weight : covered_weight) {
    if (std::fabs(weight - 1.0) > 1e-12) {
      return false;
    }
  }
  return true;
}

inline bool observation_plan_can_use_direct_missing_all(
    const ComponentObservationPlan &observed_plan,
    const ExactVariantPlan &exact_plan) {
  if (!exact_plan.no_response.terminal_leaf_survival_product ||
      !observed_plan.missing_all_branches.empty()) {
    return false;
  }
  if (observed_plan.finite_observed_branches.empty()) {
    return false;
  }
  std::vector<ObservedBranch> branches;
  branches.reserve(observed_plan.finite_observed_branches.size());
  branches.insert(
      branches.end(),
      observed_plan.finite_observed_branches.begin(),
      observed_plan.finite_observed_branches.end());
  return observed_branches_cover_exact_outcomes(exact_plan, branches);
}

inline bool observation_code_can_use_no_response_complement(
    const ComponentObservationPlan &observed_plan,
    const ExactVariantPlan &exact_plan,
    const std::size_t observed_code) {
  if (!exact_plan.no_response.terminal_leaf_survival_product ||
      observed_code >= observed_plan.keep_by_code.size() ||
      observed_code >= observed_plan.missing_rt_by_code.size()) {
    return false;
  }
  std::vector<ObservedBranch> branches;
  branches.reserve(
      observed_plan.keep_by_code[observed_code].size() +
      observed_plan.missing_rt_by_code[observed_code].size());
  branches.insert(
      branches.end(),
      observed_plan.keep_by_code[observed_code].begin(),
      observed_plan.keep_by_code[observed_code].end());
  branches.insert(
      branches.end(),
      observed_plan.missing_rt_by_code[observed_code].begin(),
      observed_plan.missing_rt_by_code[observed_code].end());
  return observed_branches_cover_exact_outcomes(exact_plan, branches);
}

inline void compile_observation_probability_modes(
    std::vector<ComponentObservationPlan> *observed_plans,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans) {
  if (observed_plans == nullptr) {
    return;
  }
  for (std::size_t component_code = 1;
       component_code < observed_plans->size() &&
       component_code < exact_variant_index_by_component_code.size();
       ++component_code) {
    auto &observed_plan = (*observed_plans)[component_code];
    observed_plan.missing_all_mode =
        MissingAllProbabilityMode::FiniteComplement;
    if (observed_plan.missing_rt_mode_by_code.size() <
        observed_plan.keep_by_code.size()) {
      observed_plan.missing_rt_mode_by_code.assign(
          observed_plan.keep_by_code.size(),
          MissingRtProbabilityMode::BranchIntegral);
    } else {
      std::fill(
          observed_plan.missing_rt_mode_by_code.begin(),
          observed_plan.missing_rt_mode_by_code.end(),
          MissingRtProbabilityMode::BranchIntegral);
    }
    if (!observed_plan.present) {
      continue;
    }
    const auto variant_index =
        exact_variant_index_by_component_code[component_code];
    if (variant_index == semantic::kInvalidIndex ||
        static_cast<std::size_t>(variant_index) >= exact_plans.size()) {
      continue;
    }
    if (observation_plan_can_use_direct_missing_all(
            observed_plan,
            exact_plans[static_cast<std::size_t>(variant_index)])) {
      observed_plan.missing_all_mode =
          MissingAllProbabilityMode::DirectNoResponse;
    }
    for (std::size_t observed_code = 1;
         observed_code < observed_plan.missing_rt_mode_by_code.size();
         ++observed_code) {
      if (observation_code_can_use_no_response_complement(
              observed_plan,
              exact_plans[static_cast<std::size_t>(variant_index)],
              observed_code)) {
        observed_plan.missing_rt_mode_by_code[observed_code] =
            MissingRtProbabilityMode::FiniteComplementOfNoResponse;
      }
    }
  }
}

inline std::vector<std::string> prepared_outcome_labels(const Rcpp::List &prep) {
  if (!prep.containsElementNamed("outcomes") || Rf_isNull(prep["outcomes"])) {
    throw std::runtime_error("prep must contain outcomes");
  }
  const Rcpp::List outcomes(prep["outcomes"]);
  const SEXP names_sexp = outcomes.names();
  if (Rf_isNull(names_sexp)) {
    throw std::runtime_error("prep outcomes must be named");
  }
  const Rcpp::CharacterVector outcome_names(names_sexp);
  std::vector<std::string> labels;
  labels.reserve(outcome_names.size());
  for (R_xlen_t i = 0; i < outcome_names.size(); ++i) {
    if (STRING_ELT(outcome_names, i) == NA_STRING) {
      throw std::runtime_error("prep outcomes must not have missing names");
    }
    labels.push_back(Rcpp::as<std::string>(outcome_names[i]));
  }
  return labels;
}

inline std::vector<std::string> prepared_component_ids(const Rcpp::List &prep) {
  if (!prep.containsElementNamed("components") || Rf_isNull(prep["components"])) {
    return {"__default__"};
  }
  const Rcpp::List components(prep["components"]);
  if (!components.containsElementNamed("ids") || Rf_isNull(components["ids"])) {
    return {"__default__"};
  }
  const Rcpp::CharacterVector ids(components["ids"]);
  if (ids.size() == 0) {
    return {"__default__"};
  }
  std::vector<std::string> out;
  out.reserve(ids.size());
  for (R_xlen_t i = 0; i < ids.size(); ++i) {
    if (STRING_ELT(ids, i) == NA_STRING) {
      throw std::runtime_error("prep component ids must not be missing");
    }
    out.push_back(Rcpp::as<std::string>(ids[i]));
  }
  return out;
}

inline std::unordered_map<std::string, semantic::Index> make_code_map(
    const std::vector<std::string> &labels) {
  std::unordered_map<std::string, semantic::Index> out;
  out.reserve(labels.size());
  for (std::size_t i = 0; i < labels.size(); ++i) {
    out.emplace(labels[i], static_cast<semantic::Index>(i + 1U));
  }
  return out;
}

inline NativeLikelihoodContext build_native_likelihood_context(
    const Rcpp::List &prep) {
  NativeLikelihoodContext ctx;
  ctx.model = compile::compile_prep(prep);
  compile_component_weight_parameter_layout(&ctx.model);
  const auto compiled = compile::project_semantic_model(ctx.model);

  const auto outcome_labels = prepared_outcome_labels(prep);
  const auto component_ids = prepared_component_ids(prep);
  const auto outcome_code_by_label = make_code_map(outcome_labels);
  const auto component_code_by_id = make_code_map(component_ids);

  ctx.observed_plans_by_component_code = build_component_observation_plans(
      prep,
      compiled,
      component_code_by_id,
      outcome_code_by_label,
      component_ids.size(),
      outcome_labels.size());
  build_exact_plan_cache(
      compiled,
      component_code_by_id,
      outcome_code_by_label,
      component_ids.size(),
      outcome_labels.size(),
      &ctx.exact_variant_index_by_component_code,
      &ctx.exact_plans);
  compile_exact_plan_leaf_row_offsets(ctx.model, &ctx.exact_plans);
  compile_observation_probability_modes(
      &ctx.observed_plans_by_component_code,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans);
  compile_observation_probability_plans(&ctx.observed_plans_by_component_code);
  ctx.observed_identity = observation_plans_are_identity(
      ctx.observed_plans_by_component_code,
      &ctx.identity_backend);
  for (const auto &plan : ctx.exact_plans) {
    if (!plan.ranked_supported) {
      ctx.ranked_supported = false;
      break;
    }
  }
  return ctx;
}

inline PreparedTrialLayout build_native_trial_layout(SEXP dataSEXP) {
  return build_prepared_trial_layout(
      Rcpp::DataFrame(dataSEXP),
      "prepared data");
}

template <typename T>
inline Rcpp::XPtr<T> checked_xptr(SEXP value, const char *what) {
  if (TYPEOF(value) != EXTPTRSXP) {
    throw std::runtime_error(std::string(what) + " must be an external pointer");
  }
  Rcpp::XPtr<T> ptr(value);
  if (ptr.get() == nullptr) {
    throw std::runtime_error(std::string(what) + " is null");
  }
  return ptr;
}

inline const NativeLikelihoodContext &likelihood_context_from_xptr(SEXP value) {
  return *checked_xptr<NativeLikelihoodContext>(value, "likelihood context");
}

inline const PreparedTrialLayout &trial_layout_from_xptr(SEXP value) {
  return *checked_xptr<PreparedTrialLayout>(value, "prepared trial layout");
}

} // namespace detail
} // namespace accumulatr::eval
