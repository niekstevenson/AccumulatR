#pragma once

#include <cstdint>
#include <vector>

#include "../semantic/model.hpp"

namespace accumulatr::runtime {

struct RuntimeLayout {
  int n_leaves{0};
  int n_pools{0};
  int n_outcomes{0};
  int n_triggers{0};
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
  int param_count{0};
};

} // namespace accumulatr::runtime
