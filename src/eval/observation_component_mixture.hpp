#pragma once

#include <Rcpp.h>

#include <cmath>
#include <cstdint>
#include <vector>

#include "observation_model.hpp"

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

inline void resolve_component_weights(
    const ComponentMixturePlan &mixture,
    const TrustedParamMatrix &params,
    const int row,
    std::vector<double> *weights) {
  const auto &component_codes = mixture.present_component_codes;
  weights->resize(component_codes.size());

  if (mixture.mode != ComponentMixtureMode::Sample) {
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      (*weights)[i] =
          mixture.component_by_code[static_cast<std::size_t>(component_codes[i])]
              .fixed_weight;
    }
    return;
  }

  double sum_nonref = 0.0;
  for (std::size_t i = 0; i < component_codes.size(); ++i) {
    const auto code = component_codes[i];
    if (code == mixture.reference_component_code) {
      continue;
    }
    const auto &component =
        mixture.component_by_code[static_cast<std::size_t>(code)];
    const double weight = component.weight_param_index >= 0
                              ? params.component_weight(
                                    row, component.weight_param_index)
                              : component.fixed_weight;
    if (!std::isfinite(weight) || weight < 0.0) {
      weights->assign(component_codes.size(), 0.0);
      return;
    }
    (*weights)[i] = weight;
    sum_nonref += weight;
  }
  if (!std::isfinite(sum_nonref) || sum_nonref > 1.0) {
    weights->assign(component_codes.size(), 0.0);
    return;
  }
  for (std::size_t i = 0; i < component_codes.size(); ++i) {
    if (component_codes[i] == mixture.reference_component_code) {
      (*weights)[i] = 1.0 - sum_nonref;
      break;
    }
  }
}

} // namespace detail
} // namespace accumulatr::eval
