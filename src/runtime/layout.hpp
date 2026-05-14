#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "../semantic/model.hpp"

namespace accumulatr::runtime {

struct RuntimeLayout {
  int n_leaves{0};
  int n_pools{0};
  int n_outcomes{0};
  int n_params{0};
  int n_triggers{0};
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

} // namespace detail

} // namespace accumulatr::runtime
