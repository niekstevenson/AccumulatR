#pragma once

#include <Rcpp.h>

#include <cstddef>
#include <vector>

#include "../semantic/model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ParamMatrixView {
  const double *base{nullptr};
  const double *onset{nullptr};
  int nrow{0};

  explicit ParamMatrixView(SEXP paramsSEXP,
                           const double *onset_ = nullptr)
      : base(REAL(paramsSEXP)), onset(onset_), nrow(Rf_nrows(paramsSEXP)) {}
};

struct ParamView {
  const ParamMatrixView *matrix{nullptr};
  const int *row_map{nullptr};
  int row_offset{0};

  ParamView() = default;

  ParamView(const ParamMatrixView &matrix_,
            const int *row_map_,
            const int row_offset_)
      : matrix(&matrix_),
        row_map(row_map_),
        row_offset(row_offset_) {}

  inline int physical_row(const int row) const {
    return row_offset + row_map[row];
  }

  inline double q(const int row) const {
    return matrix->base[physical_row(row)];
  }
};

struct ObservationLaneView {
  ParamView params;
  double observed_time{NA_REAL};
};

struct ObservationLaneBatchView {
  int physical_row(const semantic::Index leaf,
                   const std::size_t lane) const noexcept {
    if (physical_rows != nullptr) {
      return physical_rows[
          static_cast<std::size_t>(leaf) * physical_row_stride + lane];
    }
    return row_offsets[lane] + row_maps[lane][leaf];
  }

  double q(const semantic::Index leaf,
           const std::size_t lane) const noexcept {
    return matrix->base[physical_row(leaf, lane)];
  }

  ObservationLaneView operator[](
      const std::size_t lane) const noexcept {
    return ObservationLaneView{
        ParamView(*matrix, row_maps[lane], row_offsets[lane]),
        observed_times[lane]};
  }

  ObservationLaneBatchView operator+(
      const std::size_t offset) const noexcept {
    return ObservationLaneBatchView{
        matrix,
        row_maps + offset,
        row_offsets + offset,
        observed_times + offset,
        physical_rows == nullptr ? nullptr : physical_rows + offset,
        physical_row_stride,
        size - offset};
  }

  const ParamMatrixView *matrix{nullptr};
  const int *const *row_maps{nullptr};
  const int *row_offsets{nullptr};
  const double *observed_times{nullptr};
  const int *physical_rows{nullptr};
  std::size_t physical_row_stride{0U};
  std::size_t size{0U};
};

struct ObservationLaneBatch {
  void clear() noexcept {
    row_maps.clear();
    row_offsets.clear();
    observed_times.clear();
    physical_rows.clear();
    physical_row_stride = 0U;
  }

  void reserve(const std::size_t count) {
    row_maps.reserve(count);
    row_offsets.reserve(count);
    observed_times.reserve(count);
  }

  void emplace_back(const int *row_map,
                    const int row_offset,
                    const double observed_time) {
    physical_rows.clear();
    physical_row_stride = 0U;
    row_maps.push_back(row_map);
    row_offsets.push_back(row_offset);
    observed_times.push_back(observed_time);
  }

  void materialize_physical_rows(const std::size_t leaf_count) {
    const auto lane_count = size();
    physical_row_stride = lane_count;
    physical_rows.resize(leaf_count * lane_count);
    for (std::size_t leaf = 0U; leaf < leaf_count; ++leaf) {
      auto *rows = physical_rows.data() + leaf * lane_count;
      for (std::size_t lane = 0U; lane < lane_count; ++lane) {
        rows[lane] = row_offsets[lane] + row_maps[lane][leaf];
      }
    }
  }

  std::size_t size() const noexcept {
    return row_offsets.size();
  }

  bool empty() const noexcept {
    return row_offsets.empty();
  }

  ObservationLaneBatchView view(
      const ParamMatrixView &matrix) const noexcept {
    return view(matrix, 0U, size());
  }

  ObservationLaneBatchView view(
      const ParamMatrixView &matrix,
      const std::size_t begin,
      const std::size_t count) const noexcept {
    return ObservationLaneBatchView{
        &matrix,
        row_maps.data() + begin,
        row_offsets.data() + begin,
        observed_times.data() + begin,
        physical_rows.empty() ? nullptr : physical_rows.data() + begin,
        physical_row_stride,
        count};
  }

  std::vector<const int *> row_maps;
  std::vector<int> row_offsets;
  std::vector<double> observed_times;
  std::vector<int> physical_rows;
  std::size_t physical_row_stride{0U};
};

} // namespace detail
} // namespace accumulatr::eval
