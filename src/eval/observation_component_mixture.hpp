#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "observation_model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct TrustedParamMatrix {
  const double *base{nullptr};
  int nrow{0};
  int component_weight_start{0};

  TrustedParamMatrix(SEXP paramsSEXP, const int weight_param_count)
      : base(REAL(paramsSEXP)),
        nrow(Rf_nrows(paramsSEXP)),
        component_weight_start(
            Rf_ncols(paramsSEXP) - weight_param_count) {}

  inline double component_weight(const int row,
                                 const int weight_param_index) const {
    return base[
        static_cast<R_xlen_t>(component_weight_start + weight_param_index) *
            nrow +
        row];
  }
};

enum class ComponentMixtureMode : std::uint8_t {
  Fixed = 0,
  Sample = 1
};

struct ComponentMixtureEntry {
  double fixed_weight{0.0};
  int weight_param_index{-1};
};

struct ComponentMixturePlan {
  ComponentMixtureMode mode{ComponentMixtureMode::Fixed};
  semantic::Index reference_component_code{semantic::kInvalidIndex};
  int weight_param_count{0};
  std::vector<ComponentMixtureEntry> component_by_code;
  std::vector<semantic::Index> present_component_codes;
};

inline std::vector<double> resolve_component_weights(
    const ComponentMixturePlan &mixture,
    const std::vector<semantic::Index> &component_codes,
    const TrustedParamMatrix &params,
    const int row) {
  std::vector<double> weights(component_codes.size(), 0.0);
  if (component_codes.empty()) {
    return weights;
  }

  if (mixture.mode != ComponentMixtureMode::Sample) {
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      weights[i] =
          mixture.component_by_code[static_cast<std::size_t>(component_codes[i])]
              .fixed_weight;
    }
  } else {
    const auto reference_it =
        std::find(
            component_codes.begin(),
            component_codes.end(),
            mixture.reference_component_code);
    const bool reference_available = reference_it != component_codes.end();
    const auto reference_code =
        reference_available ? mixture.reference_component_code
                            : semantic::kInvalidIndex;
    double sum_nonref = 0.0;
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto code = component_codes[i];
      if (reference_available && code == reference_code) {
        continue;
      }
      const auto &component =
          mixture.component_by_code[static_cast<std::size_t>(code)];
      double weight = component.fixed_weight;
      if (component.weight_param_index >= 0) {
        weight = params.component_weight(row, component.weight_param_index);
      }
      weights[i] = weight;
      sum_nonref += weight;
    }
    if (reference_available) {
      for (std::size_t i = 0; i < component_codes.size(); ++i) {
        if (component_codes[i] != reference_code) {
          continue;
        }
        weights[i] = 1.0 - sum_nonref;
        break;
      }
    }
  }

  double total = 0.0;
  for (const auto weight : weights) {
    total += weight;
  }
  if (!(std::isfinite(total) && total > 0.0)) {
    const double uniform = 1.0 / static_cast<double>(weights.size());
    for (auto &weight : weights) {
      weight = uniform;
    }
    return weights;
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
    const ComponentMixturePlan &mixture,
    const PreparedTrialLayout &layout,
    const int *component,
    const TrustedParamMatrix &params,
    const std::size_t trial_index) {
  const auto row = static_cast<R_xlen_t>(layout.trials[trial_index].start_row);
  std::vector<TrialComponentChoice> resolved;

  if (!integer_cell_is_na(component, row)) {
    resolved.push_back(TrialComponentChoice{
        static_cast<semantic::Index>(component[row]),
        1.0});
    return resolved;
  }

  if (mixture.present_component_codes.empty()) {
    return resolved;
  }

  const auto weights = resolve_component_weights(
      mixture,
      mixture.present_component_codes,
      params,
      static_cast<int>(row));
  resolved.reserve(mixture.present_component_codes.size());
  for (std::size_t i = 0; i < mixture.present_component_codes.size(); ++i) {
    if (!(weights[i] > 0.0)) {
      continue;
    }
    resolved.push_back(
        TrialComponentChoice{mixture.present_component_codes[i], weights[i]});
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
