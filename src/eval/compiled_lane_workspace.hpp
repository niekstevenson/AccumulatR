#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "compiled_math_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline constexpr std::size_t kExactExpandedLaneTileSize = 256U;

inline semantic::Index compiled_frame_lane(
    const semantic::Index *lanes,
    const std::size_t position) noexcept {
  return lanes == nullptr ? static_cast<semantic::Index>(position)
                          : lanes[position];
}

template <std::size_t ExpansionWidth>
constexpr std::size_t exact_expanded_parent_lane_tile_size() noexcept {
  static_assert(ExpansionWidth > 0U);
  static_assert(ExpansionWidth <= kExactExpandedLaneTileSize);
  return kExactExpandedLaneTileSize / ExpansionWidth;
}

class ExactLaneSourceState;

struct CompiledLaneFrame {
  explicit CompiledLaneFrame(const CompiledMathProgram &program)
      : node_count(program.nodes.size()),
        integral_cache_count(program.integral_kernels.size()),
        time_count(static_cast<std::size_t>(program.time_slot_count)),
        source_program_count(static_cast<std::size_t>(
            program.source_program_cache_count)),
        time_valid(time_count, 0U) {}

  void ensure(const std::size_t required_lanes) {
    if (required_lanes <= stride) {
      lane_count = required_lanes;
      return;
    }

    stride = required_lanes;
    lane_count = required_lanes;

    node_values.resize(node_count * stride, 0.0);
    cache_epoch.resize(integral_cache_count * stride, 0U);
    cache_values.resize(integral_cache_count * stride, 0.0);

    source_program_epoch.resize(source_program_count * stride, 0U);
    source_program_valid_mask.resize(source_program_count * stride, 0U);
    source_program_pdf.resize(source_program_count * stride, 0.0);
    source_program_cdf.resize(source_program_count * stride, 0.0);
    source_program_survival.resize(source_program_count * stride, 1.0);

    time_values.resize(time_count * stride, 0.0);
    const auto zero =
        static_cast<std::size_t>(CompiledMathTimeSlot::Zero);
    std::fill_n(time_values.data() + zero * stride, stride, 0.0);
    source_lanes.resize(stride, 0);
    used_outcomes.resize(stride, nullptr);
  }

  void begin(const std::size_t required_lanes) {
    ensure(required_lanes);
    has_sequence_history = false;
    source_lanes_identity = true;
    advance_epoch(&cache_current_epoch, &cache_epoch);
    advance_epoch(&source_program_current_epoch, &source_program_epoch);
    std::fill(time_valid.begin(), time_valid.end(), std::uint8_t{0U});
    const auto zero =
        static_cast<std::size_t>(CompiledMathTimeSlot::Zero);
    time_valid[zero] = 1U;
  }

  std::size_t node_pos(const semantic::Index node_id,
                       const std::size_t lane) const noexcept {
    return static_cast<std::size_t>(node_id) * stride + lane;
  }

  std::size_t source_program_pos(const semantic::Index program_id,
                                 const std::size_t lane) const noexcept {
    return static_cast<std::size_t>(program_id) * stride + lane;
  }

  std::size_t integral_cache_pos(const semantic::Index cache_id,
                                 const std::size_t lane) const noexcept {
    return static_cast<std::size_t>(cache_id) * stride + lane;
  }

  std::size_t time_pos(const semantic::Index time_id,
                       const std::size_t lane) const noexcept {
    return static_cast<std::size_t>(time_id) * stride + lane;
  }

  double value(const semantic::Index node_id,
               const std::size_t lane) const noexcept {
    return node_values[node_pos(node_id, lane)];
  }

  double &value(const semantic::Index node_id,
                const std::size_t lane) noexcept {
    return node_values[node_pos(node_id, lane)];
  }

  const double *values_for(const semantic::Index node_id) const noexcept {
    return node_values.data() + static_cast<std::size_t>(node_id) * stride;
  }

  double *values_for(const semantic::Index node_id) noexcept {
    return node_values.data() + static_cast<std::size_t>(node_id) * stride;
  }

  bool has_time(const semantic::Index time_id,
                const std::size_t lane) const noexcept {
    (void)lane;
    return time_valid[static_cast<std::size_t>(time_id)] != 0U;
  }

  double time(const semantic::Index time_id,
              const std::size_t lane) const noexcept {
    return time_values[time_pos(time_id, lane)];
  }

  void copy_lanes_from(const CompiledLaneFrame &parent,
                       const std::size_t *parent_lanes,
                       const std::size_t count) noexcept {
    source_lanes_identity = false;
    const auto zero =
        static_cast<std::size_t>(CompiledMathTimeSlot::Zero);
    for (std::size_t time_id = 0; time_id < time_count; ++time_id) {
      time_valid[time_id] = parent.time_valid[time_id];
      if (time_id == zero) {
        continue;
      }
      if (time_valid[time_id] == 0U) {
        continue;
      }
      const auto parent_offset = time_id * parent.stride;
      const auto offset = time_id * stride;
      for (std::size_t lane = 0; lane < count; ++lane) {
        const auto parent_pos = parent_offset + parent_lanes[lane];
        time_values[offset + lane] = parent.time_values[parent_pos];
      }
    }
    for (std::size_t lane = 0; lane < count; ++lane) {
      source_lanes[lane] = static_cast<semantic::Index>(
          parent.source_lanes_identity
              ? parent_lanes[lane]
              : parent.source_lanes[parent_lanes[lane]]);
    }
    if (parent.has_sequence_history) {
      for (std::size_t lane = 0; lane < count; ++lane) {
        used_outcomes[lane] = parent.used_outcomes[parent_lanes[lane]];
      }
    }
    has_sequence_history = parent.has_sequence_history;
  }

  void set_time_plane(const semantic::Index time_id,
                      const double *values,
                      const std::size_t count) noexcept {
    const auto offset = static_cast<std::size_t>(time_id) * stride;
    std::copy_n(values, count, time_values.data() + offset);
    time_valid[static_cast<std::size_t>(time_id)] = 1U;
  }

  std::size_t lane_count{0U};
  std::size_t stride{0U};
  std::size_t node_count{0U};
  std::size_t integral_cache_count{0U};
  std::size_t time_count{0U};
  std::size_t source_program_count{0U};
  bool has_sequence_history{false};
  bool source_lanes_identity{true};

  std::vector<double> node_values;
  std::vector<std::uint32_t> cache_epoch;
  std::vector<double> cache_values;
  std::uint32_t cache_current_epoch{1U};

  std::vector<std::uint32_t> source_program_epoch;
  std::vector<std::uint8_t> source_program_valid_mask;
  std::vector<double> source_program_pdf;
  std::vector<double> source_program_cdf;
  std::vector<double> source_program_survival;
  std::uint32_t source_program_current_epoch{1U};

  std::vector<double> time_values;
  std::vector<std::uint8_t> time_valid;
  std::vector<semantic::Index> source_lanes;
  std::vector<const std::vector<std::uint8_t> *> used_outcomes;

private:
  template <typename T>
  static void advance_epoch(std::uint32_t *epoch,
                            std::vector<T> *epochs) {
    ++*epoch;
    if (*epoch == 0U) {
      *epoch = 1U;
      std::fill(epochs->begin(), epochs->end(), T{0});
    }
  }
};

struct CompiledLaneWorkspace {
  explicit CompiledLaneWorkspace(const CompiledMathProgram &program)
      : program_(&program), top_(program) {}

  CompiledLaneFrame &top(const std::size_t lane_count) {
    top_.begin(lane_count);
    return top_;
  }

  CompiledLaneFrame &integral_frame(const std::size_t depth,
                                    const std::size_t lane_count) {
    while (integral_frames_.size() <= depth) {
      integral_frames_.push_back(
          std::make_unique<CompiledLaneFrame>(*program_));
    }
    auto &frame = *integral_frames_[depth];
    frame.begin(lane_count);
    return frame;
  }

  void set_source_state(ExactLaneSourceState *source_state) noexcept {
    source_state_ = source_state;
  }

  std::size_t source_lane(const CompiledLaneFrame &frame,
                          const std::size_t lane) const noexcept {
    return frame.source_lanes_identity
               ? lane
               : static_cast<std::size_t>(frame.source_lanes[lane]);
  }

  ExactLaneSourceState &source_state() const noexcept {
    return *source_state_;
  }

private:
  const CompiledMathProgram *program_{nullptr};
  CompiledLaneFrame top_;
  std::vector<std::unique_ptr<CompiledLaneFrame>> integral_frames_;
  ExactLaneSourceState *source_state_{nullptr};
};

} // namespace detail
} // namespace accumulatr::eval
