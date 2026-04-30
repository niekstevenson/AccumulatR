#pragma once

#include <Rcpp.h>

#include <cstring>

#include "../runtime/layout.hpp"

namespace accumulatr::eval {
namespace detail {

struct PreparedTrialSpan {
  semantic::Index start_row{0};
  semantic::Index end_row{-1};
};

struct PreparedTrialSpansView {
  const int *start_rows{nullptr};
  const int *end_rows{nullptr};
  R_xlen_t n{0};

  std::size_t size() const {
    return static_cast<std::size_t>(n);
  }

  PreparedTrialSpan operator[](const std::size_t index) const {
    return PreparedTrialSpan{
        static_cast<semantic::Index>(start_rows[index] - 1),
        static_cast<semantic::Index>(end_rows[index] - 1)};
  }
};

struct PreparedRankColumnView {
  const int *cols{nullptr};

  int operator[](const std::size_t rank) const {
    return cols[rank - 1U] - 1;
  }
};

struct PreparedTrialLayout {
  PreparedTrialSpansView spans;
  int max_rank{1};
  int trial_col{-1};
  int component_col{-1};
  int accumulator_col{-1};
  PreparedRankColumnView label_cols;
  PreparedRankColumnView time_cols;
};

struct PreparedDataView {
  const int *trial{nullptr};
  const int *component{nullptr};
  R_xlen_t n_rows{0};
};

inline bool trial_is_selected(const int *ok,
                              const std::size_t trial_index) {
  return ok == nullptr || ok[static_cast<R_xlen_t>(trial_index)] == TRUE;
}

inline double aggregate_trial_loglik(
    const Rcpp::NumericVector &loglik,
    SEXP expandSEXP) {
  if (expandSEXP == R_NilValue || XLENGTH(expandSEXP) == 0) {
    double total_loglik = 0.0;
    for (R_xlen_t i = 0; i < loglik.size(); ++i) {
      total_loglik += static_cast<double>(loglik[i]);
    }
    return total_loglik;
  }

  double total_loglik = 0.0;
  const int *expand = INTEGER(expandSEXP);
  for (R_xlen_t i = 0; i < XLENGTH(expandSEXP); ++i) {
    total_loglik += static_cast<double>(
        loglik[static_cast<R_xlen_t>(expand[i] - 1)]);
  }
  return total_loglik;
}

inline SEXP trusted_data_column(SEXP dataSEXP, const int column_index) {
  return VECTOR_ELT(dataSEXP, column_index);
}

inline SEXP trusted_data_attr(SEXP dataSEXP, const char *name) {
  return Rf_getAttrib(dataSEXP, Rf_install(name));
}

inline int trusted_named_integer(SEXP valuesSEXP, const char *name) {
  const SEXP names = Rf_getAttrib(valuesSEXP, R_NamesSymbol);
  const R_xlen_t n = XLENGTH(valuesSEXP);
  const int *values = INTEGER(valuesSEXP);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return values[i];
    }
  }
  return NA_INTEGER;
}

inline PreparedTrialLayout read_prepared_trial_layout(SEXP dataSEXP) {
  PreparedTrialLayout layout;

  const SEXP spansSEXP = trusted_data_attr(dataSEXP, "trial_spans");
  const SEXP dimsSEXP = Rf_getAttrib(spansSEXP, R_DimSymbol);
  const int *dims = INTEGER(dimsSEXP);
  const R_xlen_t n_trials = static_cast<R_xlen_t>(dims[0]);
  const int *spans = INTEGER(spansSEXP);
  layout.spans.start_rows = spans;
  layout.spans.end_rows = spans + n_trials;
  layout.spans.n = n_trials;

  const SEXP layoutColsSEXP = trusted_data_attr(dataSEXP, "layout_cols");
  layout.trial_col = trusted_named_integer(layoutColsSEXP, "trial") - 1;
  layout.component_col = trusted_named_integer(layoutColsSEXP, "component") - 1;
  layout.accumulator_col = trusted_named_integer(layoutColsSEXP, "accumulator") - 1;

  layout.label_cols.cols = INTEGER(trusted_data_attr(dataSEXP, "label_cols"));
  layout.time_cols.cols = INTEGER(trusted_data_attr(dataSEXP, "time_cols"));
  layout.max_rank = INTEGER(trusted_data_attr(dataSEXP, "max_rank"))[0];

  return layout;
}

inline PreparedDataView read_prepared_data_view(
    SEXP dataSEXP,
    const PreparedTrialLayout &layout) {
  const SEXP component = trusted_data_column(dataSEXP, layout.component_col);
  return PreparedDataView{
      INTEGER(trusted_data_column(dataSEXP, layout.trial_col)),
      INTEGER(component),
      XLENGTH(component)};
}

inline bool integer_cell_is_na(const int *column,
                               const R_xlen_t row) {
  return column[row] == NA_INTEGER;
}

} // namespace detail
} // namespace accumulatr::eval
