#pragma once

#include <cstdint>
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

} // namespace accumulatr::runtime
