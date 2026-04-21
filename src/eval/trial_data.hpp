#pragma once

#include <Rcpp.h>

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
};

struct PreparedDataView {
  Rcpp::IntegerVector trial;
  Rcpp::IntegerVector component;
};

inline PreparedDataView read_prepared_data_view(const Rcpp::DataFrame &data) {
  const auto n_rows = data.nrows();
  return PreparedDataView{
      data.containsElementNamed("trial")
          ? Rcpp::as<Rcpp::IntegerVector>(data["trial"])
          : Rcpp::seq_len(n_rows),
      Rcpp::as<Rcpp::IntegerVector>(data["component"])};
}

inline bool integer_cell_is_na(const Rcpp::IntegerVector &column,
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
  const auto n_rows = data.nrows();
  if (n_rows == 0) {
    return layout;
  }
  const auto table = read_prepared_data_view(data);
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
