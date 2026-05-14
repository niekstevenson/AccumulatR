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
#include "observation_component_mixture.hpp"
#include "observation_model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct NativeLikelihoodContext {
  ComponentMixturePlan component_mixture;
  std::vector<ComponentObservationPlan> observation_plans_by_component_code;
  bool observation_is_identity{false};
  std::size_t outcome_count{0U};
  std::vector<semantic::Index> exact_variant_index_by_component_code;
  std::vector<ExactVariantPlan> exact_plans;
  std::vector<std::vector<int>> exact_leaf_row_offsets_by_variant;
};

inline ComponentMixturePlan build_component_mixture_plan(
    const semantic::SemanticModel &model,
    const std::unordered_map<std::string, semantic::Index> &component_code_by_id,
    const std::vector<ComponentObservationPlan> &observation_plans_by_component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code) {
  ComponentMixturePlan plan;
  plan.mode = model.component_mode == "sample"
                  ? ComponentMixtureMode::Sample
                  : ComponentMixtureMode::Fixed;
  plan.component_by_code.resize(observation_plans_by_component_code.size());

  std::unordered_map<std::string, int> index_by_name;
  auto weight_param_index = [&index_by_name](const std::string &name) {
    if (name.empty()) {
      return -1;
    }
    auto it = index_by_name.find(name);
    if (it == index_by_name.end()) {
      it = index_by_name.emplace(
          name,
          static_cast<int>(index_by_name.size())).first;
    }
    return it->second;
  };

  auto add_component = [&](const std::string &id,
                           const double fixed_weight,
                           const std::string &weight_name) {
    const auto code =
        component_code_by_id.at(id);
    auto &entry = plan.component_by_code[static_cast<std::size_t>(code)];
    entry.fixed_weight = fixed_weight;
    entry.weight_param_index = weight_param_index(weight_name);
    if (id == model.component_reference) {
      plan.reference_component_code = code;
    }
  };

  if (model.components.empty()) {
    add_component("__default__", 1.0, "");
  } else {
    for (const auto &component : model.components) {
      add_component(component.id, component.weight, component.weight_name);
    }
  }

  plan.weight_param_count = static_cast<int>(index_by_name.size());
  for (std::size_t component_code = 1;
       component_code < observation_plans_by_component_code.size() &&
       component_code < exact_variant_index_by_component_code.size();
       ++component_code) {
    if (!observation_plans_by_component_code[component_code].present ||
        exact_variant_index_by_component_code[component_code] ==
            semantic::kInvalidIndex) {
      continue;
    }
    plan.present_component_codes.push_back(
        static_cast<semantic::Index>(component_code));
  }
  if (plan.reference_component_code == semantic::kInvalidIndex &&
      !plan.present_component_codes.empty()) {
    plan.reference_component_code = plan.present_component_codes.front();
  }
  return plan;
}

inline std::vector<std::vector<int>> make_exact_leaf_row_offsets_by_variant(
    const semantic::SemanticModel &model,
    const compile::CompiledModel &compiled) {
  std::vector<std::vector<int>> offsets_by_variant;
  offsets_by_variant.reserve(compiled.variants.size());
  std::unordered_map<std::string, int> global_leaf_index;
  global_leaf_index.reserve(model.leaves.size());
  for (std::size_t i = 0; i < model.leaves.size(); ++i) {
    global_leaf_index.emplace(model.leaves[i].id, static_cast<int>(i));
  }
  for (const auto &variant : compiled.variants) {
    auto &offsets = offsets_by_variant.emplace_back();
    offsets.reserve(variant.model.leaves.size());
    for (const auto &leaf : variant.model.leaves) {
      offsets.push_back(global_leaf_index.at(leaf.id));
    }
  }
  return offsets_by_variant;
}

inline bool observation_plans_are_identity(
    const std::vector<ComponentObservationPlan> &plans) {
  for (std::size_t component_code = 1; component_code < plans.size(); ++component_code) {
    const auto &plan = plans[component_code];
    if (!plan.present) {
      continue;
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
  const auto model = compile::compile_prep(prep);
  const auto compiled = compile::project_semantic_model(model);

  const auto outcome_labels = prepared_outcome_labels(prep);
  const auto component_ids = prepared_component_ids(prep);
  ctx.outcome_count = outcome_labels.size();
  const auto outcome_code_by_label = make_code_map(outcome_labels);
  const auto component_code_by_id = make_code_map(component_ids);

  ctx.observation_plans_by_component_code = build_component_observation_plans(
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
  ctx.component_mixture = build_component_mixture_plan(
      model,
      component_code_by_id,
      ctx.observation_plans_by_component_code,
      ctx.exact_variant_index_by_component_code);
  for (std::size_t component_code = 0;
       component_code < ctx.observation_plans_by_component_code.size() &&
       component_code < ctx.exact_variant_index_by_component_code.size();
       ++component_code) {
    const auto variant_index =
        ctx.exact_variant_index_by_component_code[component_code];
    if (variant_index == semantic::kInvalidIndex ||
        static_cast<std::size_t>(variant_index) >= ctx.exact_plans.size()) {
      continue;
    }
    ctx.observation_plans_by_component_code[component_code]
        .direct_no_response =
        ctx.exact_plans[static_cast<std::size_t>(variant_index)]
            .no_response.direct_leaf_failure_product &&
        !ctx.exact_plans[static_cast<std::size_t>(variant_index)]
             .no_response.leaf_indices.empty();
  }
  ctx.exact_leaf_row_offsets_by_variant =
      make_exact_leaf_row_offsets_by_variant(model, compiled);
  compile_observation_probability_plans(&ctx.observation_plans_by_component_code);
  ctx.observation_is_identity = observation_plans_are_identity(
      ctx.observation_plans_by_component_code);
  prune_observation_planning_state(&ctx.observation_plans_by_component_code);
  return ctx;
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

} // namespace detail
} // namespace accumulatr::eval
