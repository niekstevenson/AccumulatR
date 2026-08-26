#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "eval_query.hpp"
#include "compiled_lane_workspace.hpp"
#include "exact_planner.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"
#include "../leaf/dist_kind.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactLaneLeafBatchView {
  inline int physical_row(const std::size_t lane) const noexcept {
    return physical_rows[lane];
  }

  inline double param(const std::size_t slot,
                      const int row) const noexcept {
    return base[(slot + 2U) * static_cast<std::size_t>(nrow) +
                static_cast<std::size_t>(row)];
  }

  inline double q(const std::size_t lane,
                  const int row) const noexcept {
    if (trigger_index == semantic::kInvalidIndex) {
      return 0.0;
    }
    const auto *started_values =
        shared_started == nullptr ? uniform_shared_started
                                  : shared_started[lane];
    if (started_values != nullptr) {
      const auto started =
          started_values[static_cast<std::size_t>(trigger_index)];
      if (started <= 1U) {
        return started == 1U ? 0.0 : 1.0;
      }
    }
    return base[row];
  }

  inline double t0(const int row) const noexcept {
    return base[nrow + row];
  }

  inline double onset(const int row) const noexcept {
    return onset_values == nullptr
               ? onset_abs_value
               : onset_values[row];
  }

  const double *base{nullptr};
  const double *onset_values{nullptr};
  const int *physical_rows{nullptr};
  const std::uint8_t *const *shared_started{nullptr};
  const std::uint8_t *uniform_shared_started{nullptr};
  semantic::Index trigger_index{semantic::kInvalidIndex};
  int nrow{0};
  double onset_abs_value{0.0};
};

inline double exact_leaf_q_for_trigger_state(
    const std::vector<semantic::Index> &leaf_trigger_index,
    const std::uint8_t *shared_started,
    const semantic::Index leaf_index,
    const double fallback) {
  const auto trigger_index =
      leaf_trigger_index[static_cast<std::size_t>(leaf_index)];
  if (trigger_index != semantic::kInvalidIndex &&
      shared_started != nullptr) {
    const auto started =
        shared_started[static_cast<std::size_t>(trigger_index)];
    if (started <= 1U) {
      return started == 1U ? 0.0 : 1.0;
    }
  }
  return fallback;
}

inline double exact_trigger_lane_q(
    const ObservationLaneBatchView lanes,
    const std::size_t lane,
    const semantic::Index leaf) noexcept {
  return lanes.q(leaf, lane);
}

template <typename Lane>
inline double exact_trigger_lane_q(
    const Lane *lanes,
    const std::size_t lane,
    const semantic::Index leaf) noexcept {
  return lanes[lane].params.q(leaf);
}

class ExactLaneSourceState {
public:
  ExactLaneSourceState(const ExactVariantPlan &plan,
                       const std::size_t capacity)
      : plan_(plan),
        capacity_(capacity),
        physical_rows_(plan.leaf_descriptors.size() * capacity, 0),
        physical_row_epoch_(plan.leaf_descriptors.size(), 0U),
        row_maps_(capacity, nullptr),
        row_offsets_(capacity, 0),
        shared_started_(capacity, nullptr),
        sequence_states_(capacity, nullptr),
        initial_sequence_states_(capacity, nullptr),
        active_row_maps_(row_maps_.data()),
        active_row_offsets_(row_offsets_.data()),
        active_shared_started_(shared_started_.data()),
        active_sequence_states_(sequence_states_.data()) {}

  void bind_matrix(const ParamMatrixView &matrix,
                   const std::size_t lane_count) noexcept {
    base_ = matrix.base;
    onset_ = matrix.onset;
    nrow_ = matrix.nrow;
    active_physical_rows_ = nullptr;
    active_row_maps_ = row_maps_.data();
    active_row_offsets_ = row_offsets_.data();
    active_lane_count_ = lane_count;
    advance_physical_row_epoch();
    active_shared_started_ = shared_started_.data();
    uniform_shared_started_ = nullptr;
    active_sequence_states_ = sequence_states_.data();
  }

  void bind_initial_batch(
      const ObservationLaneBatchView &lanes,
      const std::uint8_t *shared_started) noexcept {
    base_ = lanes.matrix->base;
    onset_ = lanes.matrix->onset;
    nrow_ = lanes.matrix->nrow;
    if (lanes.physical_rows != nullptr) {
      active_physical_rows_ = lanes.physical_rows;
      active_physical_row_stride_ = lanes.physical_row_stride;
    } else {
      active_physical_rows_ = nullptr;
      active_row_maps_ = lanes.row_maps;
      active_row_offsets_ = lanes.row_offsets;
      active_lane_count_ = lanes.size;
      advance_physical_row_epoch();
    }
    active_shared_started_ = nullptr;
    uniform_shared_started_ = shared_started;
    active_sequence_states_ = initial_sequence_states_.data();
  }

  void bind_lane(const std::size_t lane,
                 const ExactSequenceState &sequence_state,
                 const ParamView &params,
                 const std::uint8_t *shared_started) noexcept {
    sequence_states_[lane] = &sequence_state;
    row_maps_[lane] = params.row_map;
    row_offsets_[lane] = params.row_offset;
    shared_started_[lane] = shared_started;
  }

  ExactLaneLeafBatchView leaf_batch(
      const semantic::Index leaf_index) const noexcept {
    const auto leaf = static_cast<std::size_t>(leaf_index);
    const int *rows = nullptr;
    if (active_physical_rows_ != nullptr) {
      rows = active_physical_rows_ + leaf * active_physical_row_stride_;
    } else {
      auto *resolved = physical_rows_.data() + leaf * capacity_;
      if (physical_row_epoch_[leaf] != physical_row_current_epoch_) {
        for (std::size_t lane = 0U; lane < active_lane_count_; ++lane) {
          const auto *row_map = active_row_maps_[lane];
          resolved[lane] = active_row_offsets_[lane] +
                           (row_map == nullptr
                                ? static_cast<int>(leaf)
                                : row_map[leaf]);
        }
        physical_row_epoch_[leaf] = physical_row_current_epoch_;
      }
      rows = resolved;
    }
    return ExactLaneLeafBatchView{
        base_,
        onset_,
        rows,
        active_shared_started_,
        uniform_shared_started_,
        plan_.leaf_trigger_index[leaf],
        nrow_,
        plan_.leaf_descriptors[leaf].onset_abs_value};
  }

  bool sequence_bounds_inactive(
      const std::size_t lane,
      const semantic::Index source_id) const noexcept {
    const auto *sequence = active_sequence_states_[lane];
    if (sequence == nullptr || !sequence->has_history) {
      return true;
    }
    if (sequence->lower_bound > 0.0) {
      return false;
    }
    if (source_id == semantic::kInvalidIndex) {
      return true;
    }
    const auto source = static_cast<std::size_t>(source_id);
    return !(source < sequence->exact_times.size() &&
             std::isfinite(sequence->exact_times[source])) &&
           !(source < sequence->upper_bounds.size() &&
             std::isfinite(sequence->upper_bounds[source]));
  }

  double sequence_lower_bound(const std::size_t lane) const noexcept {
    const auto *sequence = active_sequence_states_[lane];
    return sequence != nullptr && sequence->has_history
               ? sequence->lower_bound
               : 0.0;
  }

  double sequence_exact_time(
      const std::size_t lane,
      const semantic::Index source_id) const noexcept {
    const auto *sequence = active_sequence_states_[lane];
    if (sequence == nullptr || !sequence->has_history ||
        source_id == semantic::kInvalidIndex) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const auto source = static_cast<std::size_t>(source_id);
    return source < sequence->exact_times.size()
               ? sequence->exact_times[source]
               : std::numeric_limits<double>::quiet_NaN();
  }

  double sequence_upper_bound(
      const std::size_t lane,
      const semantic::Index source_id) const noexcept {
    const auto *sequence = active_sequence_states_[lane];
    if (sequence == nullptr || !sequence->has_history ||
        source_id == semantic::kInvalidIndex) {
      return std::numeric_limits<double>::infinity();
    }
    const auto source = static_cast<std::size_t>(source_id);
    return source < sequence->upper_bounds.size()
               ? sequence->upper_bounds[source]
               : std::numeric_limits<double>::infinity();
  }

  bool expr_upper_bound(const std::size_t lane,
                        const semantic::Index expr_id,
                        ExactTimedExprUpperBound *out) const noexcept {
    const auto *sequence = active_sequence_states_[lane];
    if (sequence == nullptr || expr_id == semantic::kInvalidIndex) {
      return false;
    }
    const auto expr = static_cast<std::size_t>(expr_id);
    if (expr >= sequence->expr_upper_bounds.size() ||
        expr >= sequence->expr_upper_normalizers.size()) {
      return false;
    }
    const double time = sequence->expr_upper_bounds[expr];
    const double normalizer = sequence->expr_upper_normalizers[expr];
    if (!std::isfinite(time) || !(normalizer > 0.0)) {
      return false;
    }
    *out = ExactTimedExprUpperBound{expr_id, time, normalizer};
    return true;
  }

private:
  void advance_physical_row_epoch() noexcept {
    ++physical_row_current_epoch_;
    if (physical_row_current_epoch_ == 0U) {
      physical_row_current_epoch_ = 1U;
      std::fill(
          physical_row_epoch_.begin(), physical_row_epoch_.end(), 0U);
    }
  }

  const ExactVariantPlan &plan_;
  std::size_t capacity_{0U};
  const double *base_{nullptr};
  const double *onset_{nullptr};
  int nrow_{0};
  const std::uint8_t *uniform_shared_started_{nullptr};
  mutable std::vector<int> physical_rows_;
  mutable std::vector<std::uint32_t> physical_row_epoch_;
  std::uint32_t physical_row_current_epoch_{1U};
  std::vector<const int *> row_maps_;
  std::vector<int> row_offsets_;
  std::vector<const std::uint8_t *> shared_started_;
  std::vector<const ExactSequenceState *> sequence_states_;
  std::vector<const ExactSequenceState *> initial_sequence_states_;
  const int *active_physical_rows_{nullptr};
  std::size_t active_physical_row_stride_{0U};
  const int *const *active_row_maps_{nullptr};
  const int *active_row_offsets_{nullptr};
  std::size_t active_lane_count_{0U};
  const std::uint8_t *const *active_shared_started_{nullptr};
  const ExactSequenceState *const *active_sequence_states_{nullptr};
};

template <typename Lanes>
inline void exact_compiled_trigger_state_weights_lanes(
    const ExactVariantPlan &plan,
    const Lanes lanes,
    const std::size_t lane_count,
    const ExactCompiledTriggerState &compiled_state,
    std::vector<double> *weights) {
  weights->assign(lane_count, compiled_state.fixed_weight);
  const auto &table = plan.trigger_state_table;
  for (semantic::Index i = 0; i < compiled_state.weight_terms.size; ++i) {
    const auto &term =
        table.weight_terms[
            static_cast<std::size_t>(
                compiled_state.weight_terms.offset + i)];
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      double &weight = (*weights)[lane];
      if (!(weight > 0.0)) {
        weight = 0.0;
        continue;
      }
      const double q =
          clamp_probability(
              exact_trigger_lane_q(lanes, lane, term.leaf_index));
      weight *= term.shared_started == 0U ? q : (1.0 - q);
      if (!(weight > 0.0)) {
        weight = 0.0;
      }
    }
  }
}

inline const std::uint8_t *exact_compiled_trigger_shared_started(
    const ExactVariantPlan &plan,
    const ExactCompiledTriggerState &compiled_state) {
  const auto shared_offset =
      static_cast<std::size_t>(compiled_state.shared_started_offset);
  return shared_offset < plan.trigger_state_table.shared_started_values.size()
             ? plan.trigger_state_table.shared_started_values.data() +
                   shared_offset
             : nullptr;
}

} // namespace detail
} // namespace accumulatr::eval
