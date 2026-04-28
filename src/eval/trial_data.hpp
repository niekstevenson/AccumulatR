#pragma once

#include <Rcpp.h>

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "../runtime/layout.hpp"

namespace accumulatr::eval {
namespace detail {

struct PreparedTrialSpan {
  semantic::Index start_row{0};
  semantic::Index end_row{-1};
};

struct PreparedTrialLayout {
  std::vector<PreparedTrialSpan> spans;
  int max_rank{1};
  int trial_col{-1};
  int component_col{-1};
  int accumulator_col{-1};
  std::vector<int> label_cols;
  std::vector<int> time_cols;
};

struct PreparedDataView {
  const int *trial{nullptr};
  const int *component{nullptr};
  R_xlen_t n_rows{0};
};

struct TrialEvalControl {
  std::vector<int> weights;
  std::vector<unsigned char> selected;
};

inline TrialEvalControl build_trial_eval_control(
    const std::size_t n_trials,
    const Rcpp::LogicalVector &ok,
    const Rcpp::IntegerVector &expand) {
  TrialEvalControl out;

  if (expand.size() > 0) {
    out.weights.assign(n_trials, 0);
    for (R_xlen_t i = 0; i < expand.size(); ++i) {
      ++out.weights[static_cast<std::size_t>(expand[i] - 1)];
    }
  }

  if (ok.size() > 0) {
    bool any_false = false;
    for (std::size_t i = 0; i < n_trials; ++i) {
      if (ok[static_cast<R_xlen_t>(i)] != TRUE) {
        any_false = true;
        break;
      }
    }
    if (any_false) {
      out.selected.assign(n_trials, 0U);
      for (std::size_t i = 0; i < n_trials; ++i) {
        out.selected[i] = static_cast<unsigned char>(
            ok[static_cast<R_xlen_t>(i)] == TRUE);
      }
    }
  }
  return out;
}

inline bool trial_is_selected(const std::vector<unsigned char> *selected,
                              const std::size_t trial_index) {
  return selected == nullptr || (*selected)[trial_index] != 0U;
}

inline double aggregate_trial_loglik(
    const Rcpp::NumericVector &loglik,
    const TrialEvalControl &control) {
  if (control.weights.empty()) {
    double total_loglik = 0.0;
    for (R_xlen_t i = 0; i < loglik.size(); ++i) {
      total_loglik += static_cast<double>(loglik[i]);
    }
    return total_loglik;
  }

  double total_loglik = 0.0;
  for (std::size_t i = 0; i < control.weights.size(); ++i) {
    if (control.weights[i] == 0) {
      continue;
    }
    total_loglik += static_cast<double>(control.weights[i]) *
                    static_cast<double>(loglik[static_cast<R_xlen_t>(i)]);
  }
  return total_loglik;
}

inline SEXP trusted_data_column(SEXP dataSEXP, const char *name) {
  const SEXP names = Rf_getAttrib(dataSEXP, R_NamesSymbol);
  const R_xlen_t n_cols = XLENGTH(dataSEXP);
  for (R_xlen_t i = 0; i < n_cols; ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(dataSEXP, i);
    }
  }
  return R_NilValue;
}

inline SEXP trusted_data_column(SEXP dataSEXP, const int column_index) {
  return VECTOR_ELT(dataSEXP, column_index);
}

inline int trusted_data_column_index(SEXP dataSEXP, const char *name) {
  const SEXP names = Rf_getAttrib(dataSEXP, R_NamesSymbol);
  const R_xlen_t n_cols = XLENGTH(dataSEXP);
  for (R_xlen_t i = 0; i < n_cols; ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

inline PreparedDataView read_prepared_data_view(SEXP dataSEXP) {
  const auto n_rows = XLENGTH(trusted_data_column(dataSEXP, "component"));
  return PreparedDataView{
      INTEGER(trusted_data_column(dataSEXP, "trial")),
      INTEGER(trusted_data_column(dataSEXP, "component")),
      n_rows};
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

inline PreparedDataView read_prepared_data_view(const Rcpp::DataFrame &data) {
  return read_prepared_data_view(static_cast<SEXP>(data));
}

inline bool integer_cell_is_na(const Rcpp::IntegerVector &column,
                               const R_xlen_t row) {
  return column[row] == NA_INTEGER;
}

inline bool integer_cell_is_na(const int *column,
                               const R_xlen_t row) {
  return column[row] == NA_INTEGER;
}

inline int detect_rank_count(const Rcpp::DataFrame &data,
                             const std::string &context) {
  const std::string prefix = context.empty() ? "" : (context + " ");
  int max_rank = 1;
  for (int rank = 2;; ++rank) {
    const auto r_col = "R" + std::to_string(rank);
    const auto rt_col = "rt" + std::to_string(rank);
    const bool has_r = data.containsElementNamed(r_col.c_str());
    const bool has_rt = data.containsElementNamed(rt_col.c_str());
    if (has_r != has_rt) {
      throw std::runtime_error(
          prefix + "ranked columns must appear as matched Rk/rtk pairs");
    }
    if (!has_r) {
      break;
    }
    max_rank = rank;
  }
  return max_rank;
}

inline PreparedTrialLayout build_prepared_trial_layout(
    const Rcpp::DataFrame &data,
    const std::string &context = std::string()) {
  PreparedTrialLayout layout;
  layout.max_rank = detect_rank_count(data, context);
  const SEXP dataSEXP = static_cast<SEXP>(data);
  layout.trial_col = trusted_data_column_index(dataSEXP, "trial");
  layout.component_col = trusted_data_column_index(dataSEXP, "component");
  layout.accumulator_col = trusted_data_column_index(dataSEXP, "accumulator");
  layout.label_cols.assign(static_cast<std::size_t>(layout.max_rank + 1), -1);
  layout.time_cols.assign(static_cast<std::size_t>(layout.max_rank + 1), -1);
  layout.label_cols[1] = trusted_data_column_index(dataSEXP, "R");
  layout.time_cols[1] = trusted_data_column_index(dataSEXP, "rt");
  for (int rank = 2; rank <= layout.max_rank; ++rank) {
    layout.label_cols[static_cast<std::size_t>(rank)] =
        trusted_data_column_index(dataSEXP, ("R" + std::to_string(rank)).c_str());
    layout.time_cols[static_cast<std::size_t>(rank)] =
        trusted_data_column_index(dataSEXP, ("rt" + std::to_string(rank)).c_str());
  }
  const auto n_rows = data.nrows();
  if (n_rows == 0) {
    return layout;
  }
  const auto table = read_prepared_data_view(dataSEXP, layout);
  semantic::Index start = 0;
  int last_trial = table.trial[0];
  for (R_xlen_t row = 1; row < n_rows; ++row) {
    if (table.trial[row] == last_trial) {
      continue;
    }
    layout.spans.push_back(
        PreparedTrialSpan{start, static_cast<semantic::Index>(row - 1)});
    start = static_cast<semantic::Index>(row);
    last_trial = table.trial[row];
  }
  layout.spans.push_back(
      PreparedTrialSpan{start, static_cast<semantic::Index>(n_rows - 1)});

  return layout;
}

template <typename T>
inline std::vector<runtime::TrialBlock> build_variant_blocks(
    const std::vector<T> &records) {
  std::vector<runtime::TrialBlock> blocks;
  if (records.empty()) {
    return blocks;
  }
  runtime::TrialBlock current;
  current.variant_index = records.front().variant_index;
  current.start_row = 0;
  current.row_count = 1;
  for (std::size_t i = 1; i < records.size(); ++i) {
    if (records[i].variant_index == current.variant_index) {
      ++current.row_count;
      continue;
    }
    blocks.push_back(current);
    current.variant_index = records[i].variant_index;
    current.start_row = static_cast<int>(i);
    current.row_count = 1;
  }
  blocks.push_back(current);
  return blocks;
}

} // namespace detail
} // namespace accumulatr::eval
