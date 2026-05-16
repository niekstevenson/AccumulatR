#pragma once

#include <cstdint>
#include <vector>

#include "layout.hpp"

namespace accumulatr::runtime {

struct ExactEvaluationProgram {
  RuntimeLayout layout{};

  std::vector<LeafRuntimeDescriptor> leaf_descriptors;
  std::vector<std::uint8_t> leaf_dist_kind;

  std::vector<std::uint8_t> onset_kind;
  std::vector<std::uint8_t> onset_source_kind;
  std::vector<semantic::Index> onset_source_index;
  std::vector<semantic::Index> onset_source_ids;
  std::vector<double> onset_lag;
  std::vector<double> onset_abs_value;

  std::vector<semantic::Index> leaf_trigger_index;

  std::vector<semantic::Index> trigger_member_offsets;
  std::vector<semantic::Index> trigger_member_indices;

  std::vector<semantic::Index> pool_k;
  std::vector<semantic::Index> pool_member_offsets;
  std::vector<semantic::Index> pool_member_indices;
  std::vector<std::uint8_t> pool_member_kind;
  std::vector<semantic::Index> pool_member_source_ids;

  std::vector<std::uint8_t> expr_kind;
  std::vector<semantic::Index> expr_arg_offsets;
  std::vector<semantic::Index> expr_args;
  std::vector<semantic::Index> expr_ref_child;
  std::vector<semantic::Index> expr_blocker_child;
  std::vector<semantic::Index> expr_source_index;
  std::vector<std::uint8_t> expr_source_kind;
  std::vector<semantic::Index> expr_source_ids;
  std::vector<int> expr_event_k;

  std::vector<semantic::Index> outcome_expr_root;
  std::vector<semantic::Index> outcome_codes;
  std::vector<semantic::Index> outcome_competitor_offsets;
  std::vector<semantic::Index> outcome_competitor_expr_roots;
  std::vector<semantic::Index> outcome_competitor_indices;
};

} // namespace accumulatr::runtime
