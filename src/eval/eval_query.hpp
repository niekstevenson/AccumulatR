#pragma once

#include <Rcpp.h>

#include <vector>

#include "../semantic/model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ParamView {
  const double *base;
  const double *onset{nullptr};
  int nrow;
  const int *row_map{nullptr};
  int row_offset{0};

  explicit ParamView(SEXP paramsSEXP)
      : base(REAL(paramsSEXP)), nrow(Rf_nrows(paramsSEXP)) {}

  ParamView(SEXP paramsSEXP, const double *onset_)
      : base(REAL(paramsSEXP)), onset(onset_), nrow(Rf_nrows(paramsSEXP)) {}

  ParamView(SEXP paramsSEXP, const double *onset_, const int *row_map_)
      : base(REAL(paramsSEXP)),
        onset(onset_),
        nrow(Rf_nrows(paramsSEXP)),
        row_map(row_map_) {}

  ParamView(SEXP paramsSEXP,
            const double *onset_,
            const int *row_map_,
            const int row_offset_)
      : base(REAL(paramsSEXP)),
        onset(onset_),
        nrow(Rf_nrows(paramsSEXP)),
        row_map(row_map_),
        row_offset(row_offset_) {}

  inline int physical_row(const int row) const {
    return row_offset + (row_map == nullptr ? row : row_map[row]);
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

  inline double onset_abs(const int row, const double fallback) const {
    return onset == nullptr ? fallback : onset[physical_row(row)];
  }
};

} // namespace detail
} // namespace accumulatr::eval
