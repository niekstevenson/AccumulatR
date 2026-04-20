#pragma once

#include <Rcpp.h>

#include <string>
#include <vector>

#include "../runtime/layout.hpp"

namespace accumulatr::eval {
namespace detail {

struct TrialTableView {
  Rcpp::IntegerVector trial;
  bool has_component{false};
  Rcpp::CharacterVector component;
};

inline TrialTableView read_trial_table_view(const Rcpp::DataFrame &data) {
  const auto n_rows = data.nrows();
  const bool has_component = data.containsElementNamed("component");
  return TrialTableView{
      data.containsElementNamed("trial")
          ? Rcpp::as<Rcpp::IntegerVector>(data["trial"])
          : Rcpp::seq_len(n_rows),
      has_component,
      has_component ? Rcpp::as<Rcpp::CharacterVector>(data["component"])
                    : Rcpp::CharacterVector(n_rows, NA_STRING)};
}

inline std::string component_id_at_row(const TrialTableView &table,
                                       const R_xlen_t row,
                                       const std::string &default_component = "__default__") {
  if (table.has_component && STRING_ELT(table.component, row) != NA_STRING) {
    return Rcpp::as<std::string>(table.component[row]);
  }
  return default_component;
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
