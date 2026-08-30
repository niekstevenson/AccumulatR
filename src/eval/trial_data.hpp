#pragma once

#include <Rcpp.h>

#include <cmath>
#include <cstring>

namespace accumulatr::eval {
namespace detail {

struct PreparedRankColumnView {
  const int *cols{nullptr};

  int operator[](const std::size_t rank) const {
    return cols[rank - 1U] - 1;
  }
};

struct PreparedObservationColumnView {
  int lt{-1};
  int ut{-1};
  int lc{-1};
  int uc{-1};
  int missingness{-1};

  bool present() const noexcept {
    return lt >= 0;
  }
};

struct PreparedTrialLayout {
  int max_rank{1};
  int component_col{-1};
  int onset_col{-1};
  PreparedRankColumnView label_cols;
  PreparedRankColumnView time_cols;
  PreparedObservationColumnView observation;
};

struct PreparedObservationDataView {
  const double *lt{nullptr};
  const double *ut{nullptr};
  const double *lc{nullptr};
  const double *uc{nullptr};
  const int *missingness{nullptr};
};

struct ObservationBounds {
  double trunc_lower{0.0};
  double trunc_upper{R_PosInf};
  double censor_lower{0.0};
  double censor_upper{R_PosInf};
  int missingness{NA_INTEGER};

  bool truncates() const noexcept {
    return trunc_lower > 0.0 || std::isfinite(trunc_upper);
  }

  bool censored() const noexcept {
    return missingness != NA_INTEGER;
  }
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

inline int trusted_named_column(SEXP valuesSEXP, const char *name) {
  const int value = trusted_named_integer(valuesSEXP, name);
  return value == NA_INTEGER ? -1 : value - 1;
}

inline PreparedTrialLayout read_prepared_trial_layout(
    SEXP dataSEXP) {
  PreparedTrialLayout layout;

  const SEXP layoutColsSEXP = trusted_data_attr(dataSEXP, "layout_cols");
  layout.component_col = trusted_named_column(layoutColsSEXP, "component");
  layout.onset_col = trusted_named_column(layoutColsSEXP, "onset");
  layout.observation.lt = trusted_named_column(layoutColsSEXP, "LT");
  layout.observation.ut = trusted_named_column(layoutColsSEXP, "UT");
  layout.observation.lc = trusted_named_column(layoutColsSEXP, "LC");
  layout.observation.uc = trusted_named_column(layoutColsSEXP, "UC");
  layout.observation.missingness =
      trusted_named_column(layoutColsSEXP, "missingness");

  layout.label_cols.cols = INTEGER(trusted_data_attr(dataSEXP, "label_cols"));
  layout.time_cols.cols = INTEGER(trusted_data_attr(dataSEXP, "time_cols"));
  layout.max_rank = INTEGER(trusted_data_attr(dataSEXP, "max_rank"))[0];

  return layout;
}

inline PreparedObservationDataView read_prepared_observation_data_view(
    SEXP dataSEXP,
    const PreparedTrialLayout &layout) {
  return PreparedObservationDataView{
      REAL(trusted_data_column(dataSEXP, layout.observation.lt)),
      REAL(trusted_data_column(dataSEXP, layout.observation.ut)),
      REAL(trusted_data_column(dataSEXP, layout.observation.lc)),
      REAL(trusted_data_column(dataSEXP, layout.observation.uc)),
      INTEGER(trusted_data_column(
          dataSEXP, layout.observation.missingness))};
}

inline ObservationBounds observation_bounds_for_row(
    const PreparedObservationDataView &view,
    const R_xlen_t row) noexcept {
  return ObservationBounds{
      view.lt[row],
      view.ut[row],
      view.lc[row],
      view.uc[row],
      view.missingness[row]};
}

inline bool integer_cell_is_na(const int *column,
                               const R_xlen_t row) {
  return column[row] == NA_INTEGER;
}

} // namespace detail
} // namespace accumulatr::eval
