#pragma once

#include <Rcpp.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "../compile/project_semantic.hpp"

namespace accumulatr::runtime {

struct RuntimeLayout {
  int n_leaves{0};
  int n_pools{0};
  int n_outcomes{0};
  int n_params{0};
  int n_triggers{0};
};

struct TrialBlock {
  int variant_index{-1};
  int start_row{0};
  int row_count{0};
};

struct ParameterLayout {
  std::vector<semantic::Index> leaf_param_offsets;
  std::vector<semantic::Index> leaf_param_slots;
  std::vector<semantic::Index> leaf_q_slots;
  std::vector<semantic::Index> leaf_t0_slots;
};

struct LeafRuntimeDescriptor {
  std::uint8_t dist_kind{0};
  std::uint8_t onset_kind{0};
  std::uint8_t onset_source_kind{0};
  semantic::Index onset_source_index{semantic::kInvalidIndex};
  semantic::Index onset_source_id{semantic::kInvalidIndex};
  double onset_lag{0.0};
  double onset_abs_value{0.0};
  semantic::Index trigger_index{semantic::kInvalidIndex};
  semantic::Index param_offset{0};
  int param_count{0};
};

namespace detail {

class SlotAllocator {
public:
  semantic::Index slot_for(const std::string &key) {
    if (key.empty()) {
      return semantic::kInvalidIndex;
    }
    const auto it = slot_by_key_.find(key);
    if (it != slot_by_key_.end()) {
      return it->second;
    }
    const auto slot = static_cast<semantic::Index>(keys_.size());
    keys_.push_back(key);
    slot_by_key_.emplace(key, slot);
    return slot;
  }

  const std::vector<std::string> &keys() const noexcept {
    return keys_;
  }

private:
  std::unordered_map<std::string, semantic::Index> slot_by_key_;
  std::vector<std::string> keys_;
};

template <typename T>
Rcpp::IntegerVector as_integer_vector(const std::vector<T> &values) {
  Rcpp::IntegerVector out(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = static_cast<int>(values[i]);
  }
  return out;
}

inline Rcpp::NumericVector as_numeric_vector(
    const std::vector<double> &values) {
  Rcpp::NumericVector out(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = values[i];
  }
  return out;
}

inline Rcpp::List runtime_layout_to_r_list(const RuntimeLayout &layout) {
  return Rcpp::List::create(
      Rcpp::Named("n_leaves") = layout.n_leaves,
      Rcpp::Named("n_pools") = layout.n_pools,
      Rcpp::Named("n_outcomes") = layout.n_outcomes,
      Rcpp::Named("n_params") = layout.n_params,
      Rcpp::Named("n_triggers") = layout.n_triggers);
}

inline Rcpp::List parameter_layout_to_r_list(
    const ParameterLayout &layout) {
  return Rcpp::List::create(
      Rcpp::Named("leaf_param_offsets") =
          as_integer_vector(layout.leaf_param_offsets),
      Rcpp::Named("leaf_param_slots") =
          as_integer_vector(layout.leaf_param_slots),
      Rcpp::Named("leaf_q_slots") = as_integer_vector(layout.leaf_q_slots),
      Rcpp::Named("leaf_t0_slots") =
          as_integer_vector(layout.leaf_t0_slots));
}

} // namespace detail

} // namespace accumulatr::runtime
