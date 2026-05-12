#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <vector>

#include "observation_model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct TrustedParamMatrix {
  const double *base{nullptr};
  int nrow{0};
  int component_weight_start{0};

  TrustedParamMatrix(SEXP paramsSEXP, const int component_weight_param_count)
      : base(REAL(paramsSEXP)),
        nrow(Rf_nrows(paramsSEXP)),
        component_weight_start(
            Rf_ncols(paramsSEXP) - component_weight_param_count) {}

  inline double component_weight(const int row,
                                 const int weight_param_index) const {
    return base[
        static_cast<R_xlen_t>(component_weight_start + weight_param_index) *
            nrow +
        row];
  }

};

inline std::vector<double> resolve_component_weights(
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &component_codes,
    const TrustedParamMatrix &params,
    const int row) {
  std::vector<double> weights(component_codes.size(), 0.0);
  if (component_codes.empty()) {
    return weights;
  }

  if (model.component_mode != "sample") {
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      weights[i] = model.components[component_index].weight;
    }
  } else {
    const auto reference_id =
        !model.component_reference.empty() && std::any_of(
                                                 component_codes.begin(),
                                                 component_codes.end(),
                                                 [&](const auto code) {
                                                   return model.components
                                                              [static_cast<std::size_t>(code - 1)]
                                                                  .id ==
                                                          model.component_reference;
                                                 })
            ? model.component_reference
            : model.components[static_cast<std::size_t>(component_codes.front() - 1)]
                  .id;
    double sum_nonref = 0.0;
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      const auto &component = model.components[component_index];
      if (component.id == reference_id) {
        continue;
      }
      double weight = component.weight;
      if (component.weight_param_index >= 0) {
        weight = params.component_weight(row, component.weight_param_index);
      }
      weights[i] = weight;
      sum_nonref += weight;
    }
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      const auto &component = model.components[component_index];
      if (component.id != reference_id) {
        continue;
      }
      double ref_weight = 1.0 - sum_nonref;
      weights[i] = ref_weight;
      break;
    }
  }

  double total = 0.0;
  for (const auto weight : weights) {
    total += weight;
  }
  const double inv_total = 1.0 / total;
  for (auto &weight : weights) {
    weight *= inv_total;
  }
  return weights;
}

struct TrialComponentChoice {
  semantic::Index component_code{semantic::kInvalidIndex};
  double weight{0.0};
};

inline std::vector<TrialComponentChoice> resolve_trial_components(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const semantic::SemanticModel &model,
    const PreparedTrialLayout &layout,
    const int *component,
    const TrustedParamMatrix &params,
    const std::size_t trial_index) {
  const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
  std::vector<TrialComponentChoice> resolved;

  if (!integer_cell_is_na(component, row)) {
    resolved.push_back(TrialComponentChoice{
        static_cast<semantic::Index>(component[row]),
        1.0});
    return resolved;
  }

  std::vector<semantic::Index> component_codes;
  const auto max_component_code = std::min(
      component_plans_by_code.size(),
      model.components.size() + 1U);
  for (std::size_t component_code = 1; component_code < max_component_code;
       ++component_code) {
    if (component_plans_by_code[component_code].present) {
      component_codes.push_back(static_cast<semantic::Index>(component_code));
    }
  }
  if (component_codes.empty()) {
    return resolved;
  }

  const auto weights = resolve_component_weights(
      model,
      component_codes,
      params,
      static_cast<int>(row));
  resolved.reserve(component_codes.size());
  for (std::size_t i = 0; i < component_codes.size(); ++i) {
    if (!(weights[i] > 0.0)) {
      continue;
    }
    resolved.push_back(TrialComponentChoice{component_codes[i], weights[i]});
  }
  return resolved;
}

inline semantic::Index resolve_variant_index_by_component_code(
    const semantic::Index component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code) {
  const auto idx = static_cast<std::size_t>(component_code);
  return exact_variant_index_by_component_code[idx];
}

} // namespace detail
} // namespace accumulatr::eval
