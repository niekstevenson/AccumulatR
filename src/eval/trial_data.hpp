#pragma once

#include <Rcpp.h>

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

struct PreparedTrialLayout {
  PreparedTrialRowsView trials;
  const int *trial_weights{nullptr};
  R_xlen_t n_trial_weights{0};
  int max_rank{1};
  int component_col{-1};
  PreparedRankColumnView label_cols;
  PreparedRankColumnView time_cols;
};

struct PreparedDataView {
  const int *component{nullptr};
  R_xlen_t n_rows{0};
};

inline bool trial_is_selected(const int *ok,
                              const std::size_t trial_index) {
  return ok == nullptr || ok[static_cast<R_xlen_t>(trial_index)] == TRUE;
}

inline double prepared_trial_weight(
    const PreparedTrialLayout &layout,
    const std::size_t trial_index) {
  if (layout.trial_weights == nullptr ||
      static_cast<R_xlen_t>(trial_index) >= layout.n_trial_weights) {
    return 1.0;
  }
  return static_cast<double>(
      layout.trial_weights[static_cast<R_xlen_t>(trial_index)]);
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

inline PreparedTrialLayout read_prepared_trial_layout(
    SEXP dataSEXP,
    SEXP trialWeightsSEXP = R_NilValue) {
  PreparedTrialLayout layout;

  const SEXP startsSEXP = trusted_data_attr(dataSEXP, "trial_start_rows");
  layout.trials.start_rows = INTEGER(startsSEXP);
  layout.trials.n = XLENGTH(startsSEXP);

  SEXP weightsSEXP = trialWeightsSEXP;
  if (weightsSEXP == R_NilValue || XLENGTH(weightsSEXP) == 0) {
    weightsSEXP = trusted_data_attr(dataSEXP, "trial_weights");
  }
  if (weightsSEXP != R_NilValue && XLENGTH(weightsSEXP) > 0) {
    layout.trial_weights = INTEGER(weightsSEXP);
    layout.n_trial_weights = XLENGTH(weightsSEXP);
  }

  const SEXP layoutColsSEXP = trusted_data_attr(dataSEXP, "layout_cols");
  layout.component_col = trusted_named_integer(layoutColsSEXP, "component") - 1;

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

inline bool integer_cell_is_na(const int *column,
                               const R_xlen_t row) {
  return column[row] == NA_INTEGER;
}

} // namespace detail
} // namespace accumulatr::eval
