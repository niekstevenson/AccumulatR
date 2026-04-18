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

} // namespace accumulatr::runtime
