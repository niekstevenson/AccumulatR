#pragma once

#include <Rcpp.h>

#include <vector>

#include "../semantic/model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ParamView {
  const double *base;
  int nrow;
  const int *row_map{nullptr};

  explicit ParamView(SEXP paramsSEXP)
      : base(REAL(paramsSEXP)), nrow(Rf_nrows(paramsSEXP)) {}

  ParamView(SEXP paramsSEXP, const int *row_map_)
      : base(REAL(paramsSEXP)),
        nrow(Rf_nrows(paramsSEXP)),
        row_map(row_map_) {}

  inline int physical_row(const int row) const {
    return row_map == nullptr ? row : row_map[row];
  }

  inline double q(const int row) const {
    return base[physical_row(row)];
  }

  inline double t0(const int row) const {
    return base[nrow + physical_row(row)];
  }

  inline double p(const int row, const int slot) const {
    return base[(slot + 2) * nrow + physical_row(row)];
  }
};

struct ProbabilityQuery {
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index outcome_code{semantic::kInvalidIndex};
};

struct EvalLoglikQuery {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index outcome_code{semantic::kInvalidIndex};
  double rt{NA_REAL};
  semantic::Index row_map_index{semantic::kInvalidIndex};
};

struct EvalProbabilityQuery {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index outcome_code{semantic::kInvalidIndex};
  semantic::Index row_map_index{semantic::kInvalidIndex};
};

inline std::vector<ProbabilityQuery> collapse_probability_queries(
    const Rcpp::DataFrame &data,
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const PreparedTrialLayout &layout) {
  if (layout.spans.empty()) {
    return {};
  }

  const auto table = read_prepared_data_view(data);
  Rcpp::IntegerVector outcome = Rcpp::as<Rcpp::IntegerVector>(data["R"]);

  std::vector<ProbabilityQuery> out;
  out.reserve(layout.spans.size());
  for (const auto &span : layout.spans) {
    const auto row = static_cast<R_xlen_t>(span.start_row);
    const int component_code = table.component[row];
    out.push_back(ProbabilityQuery{
        variant_index_by_component_code[static_cast<std::size_t>(component_code)],
        outcome[row]});
  }
  return out;
}

} // namespace detail
} // namespace accumulatr::eval
