#pragma once

#include <Rcpp.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../compile/prep_to_semantic.hpp"
#include "../compile/project_semantic.hpp"
#include "../semantic/model.hpp"
#include "direct_kernel.hpp"
#include "exact_sequence.hpp"
#include "observed_plan.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct NativeLikelihoodContext {
  semantic::SemanticModel model;
  std::vector<ComponentObservationPlan> observed_plans_by_component_code;
  bool observed_identity{false};
  compile::BackendKind identity_backend{compile::BackendKind::Direct};
  bool ranked_supported{true};
  std::vector<semantic::Index> direct_variant_index_by_component_code;
  std::vector<VariantPlan> direct_plans;
  std::vector<semantic::Index> exact_variant_index_by_component_code;
  std::vector<ExactVariantPlan> exact_plans;
};

inline bool observation_plans_are_identity(
    const std::vector<ComponentObservationPlan> &plans,
    compile::BackendKind *identity_backend = nullptr) {
  compile::BackendKind backend = compile::BackendKind::Direct;
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
  ctx.observed_identity = observation_plans_are_identity(
      ctx.observed_plans_by_component_code,
      &ctx.identity_backend);
  build_direct_plan_cache(
      compiled,
      component_code_by_id,
      outcome_code_by_label,
      component_ids.size(),
      outcome_labels.size(),
      &ctx.direct_variant_index_by_component_code,
      &ctx.direct_plans);
  build_exact_plan_cache(
      compiled,
      component_code_by_id,
      outcome_code_by_label,
      component_ids.size(),
      outcome_labels.size(),
      &ctx.exact_variant_index_by_component_code,
      &ctx.exact_plans);
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
