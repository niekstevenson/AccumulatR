#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstring>

#include "../runtime/layout.hpp"

namespace accumulatr::eval {
namespace detail {

struct PreparedTrialRow {
  semantic::Index start_row{0};
};

struct PreparedTrialRowsView {
  const int *start_rows{nullptr};
  R_xlen_t n{0};

  std::size_t size() const {
    return static_cast<std::size_t>(n);
  }

  PreparedTrialRow operator[](const std::size_t index) const {
    return PreparedTrialRow{
        static_cast<semantic::Index>(start_rows[index] - 1)};
  }
};

struct PreparedRankColumnView {
  const int *cols{nullptr};

  int operator[](const std::size_t rank) const {
    return cols[rank - 1U] - 1;
  }
};

struct PreparedBoundColumnView {
  int lt_col{-1};
  int ut_col{-1};
  int lc_col{-1};
  int uc_col{-1};
  bool present{false};
};

struct PreparedTrialLayout {
  PreparedTrialRowsView trials;
  int max_rank{1};
  int component_col{-1};
  int onset_col{-1};
  PreparedRankColumnView label_cols;
  PreparedRankColumnView time_cols;
  PreparedBoundColumnView bounds;
};

struct PreparedDataView {
  const int *component{nullptr};
  R_xlen_t n_rows{0};
};

struct PreparedBoundDataView {
  const double *lt{nullptr};
  const double *ut{nullptr};
  const double *lc{nullptr};
  const double *uc{nullptr};
};

struct ObservationBounds {
  bool active{false};
  bool interval{false};
  bool censored{false};
  bool truncation_active{false};
  bool right_censored{false};
  bool left_censored{false};
  double trunc_lower{0.0};
  double trunc_upper{R_PosInf};
  double event_lower{0.0};
  double event_upper{R_PosInf};
};

inline bool trial_is_selected(const int *ok,
                              const std::size_t trial_index) {
  return ok == nullptr || ok[static_cast<R_xlen_t>(trial_index)] == TRUE;
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

inline int trusted_named_column_index(SEXP valuesSEXP, const char *name) {
  const int value = trusted_named_integer(valuesSEXP, name);
  return value == NA_INTEGER ? -1 : value - 1;
}

inline PreparedTrialLayout read_prepared_trial_layout(
    SEXP dataSEXP) {
  PreparedTrialLayout layout;

  const SEXP startsSEXP = trusted_data_attr(dataSEXP, "trials_start_rows");
  layout.trials.start_rows = INTEGER(startsSEXP);
  layout.trials.n = XLENGTH(startsSEXP);

  const SEXP layoutColsSEXP = trusted_data_attr(dataSEXP, "layout_cols");
  layout.component_col = trusted_named_integer(layoutColsSEXP, "component") - 1;
  layout.onset_col = trusted_named_integer(layoutColsSEXP, "onset") - 1;
  layout.bounds.lt_col = trusted_named_column_index(layoutColsSEXP, "LT");
  layout.bounds.ut_col = trusted_named_column_index(layoutColsSEXP, "UT");
  layout.bounds.lc_col = trusted_named_column_index(layoutColsSEXP, "LC");
  layout.bounds.uc_col = trusted_named_column_index(layoutColsSEXP, "UC");
  layout.bounds.present =
      layout.bounds.lt_col >= 0 ||
      layout.bounds.ut_col >= 0 ||
      layout.bounds.lc_col >= 0 ||
      layout.bounds.uc_col >= 0;

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
      INTEGER(component),
      XLENGTH(component)};
}

inline PreparedBoundDataView read_prepared_bound_data_view(
    SEXP dataSEXP,
    const PreparedTrialLayout &layout) {
  PreparedBoundDataView view;
  if (!layout.bounds.present) {
    return view;
  }
  if (layout.bounds.lt_col >= 0) {
    view.lt = REAL(trusted_data_column(dataSEXP, layout.bounds.lt_col));
  }
  if (layout.bounds.ut_col >= 0) {
    view.ut = REAL(trusted_data_column(dataSEXP, layout.bounds.ut_col));
  }
  if (layout.bounds.lc_col >= 0) {
    view.lc = REAL(trusted_data_column(dataSEXP, layout.bounds.lc_col));
  }
  if (layout.bounds.uc_col >= 0) {
    view.uc = REAL(trusted_data_column(dataSEXP, layout.bounds.uc_col));
  }
  return view;
}

inline double bound_cell_or_default(const double *column,
                                    const R_xlen_t row,
                                    const double default_value) {
  if (column == nullptr || Rcpp::NumericVector::is_na(column[row])) {
    return default_value;
  }
  return column[row];
}

inline double bound_cell_or_na(const double *column,
                               const R_xlen_t row) {
  if (column == nullptr || Rcpp::NumericVector::is_na(column[row])) {
    return NA_REAL;
  }
  return column[row];
}

inline ObservationBounds observation_bounds_for_row(
    const PreparedBoundDataView &view,
    const R_xlen_t row,
    const double observed_rt) {
  ObservationBounds bounds;
  bounds.trunc_lower = bound_cell_or_default(view.lt, row, 0.0);
  bounds.trunc_upper = bound_cell_or_default(view.ut, row, R_PosInf);
  const double lc = bound_cell_or_na(view.lc, row);
  const double uc = bound_cell_or_na(view.uc, row);
  bounds.left_censored = !Rcpp::NumericVector::is_na(lc);
  bounds.right_censored = !Rcpp::NumericVector::is_na(uc);
  bounds.censored = bounds.left_censored || bounds.right_censored;
  bounds.truncation_active =
      bounds.trunc_lower > 0.0 || std::isfinite(bounds.trunc_upper);
  bounds.active =
      bounds.truncation_active ||
      (Rcpp::NumericVector::is_na(observed_rt) && bounds.censored);
  bounds.interval =
      Rcpp::NumericVector::is_na(observed_rt) && bounds.active;
  bounds.event_lower =
      bounds.right_censored
          ? std::max(bounds.trunc_lower, uc)
          : bounds.trunc_lower;
  bounds.event_upper =
      bounds.left_censored
          ? std::min(bounds.trunc_upper, lc)
          : bounds.trunc_upper;
  return bounds;
}

inline bool integer_cell_is_na(const int *column,
                               const R_xlen_t row) {
  return column[row] == NA_INTEGER;
}

} // namespace detail
} // namespace accumulatr::eval
