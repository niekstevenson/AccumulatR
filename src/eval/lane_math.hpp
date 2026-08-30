#pragma once

#include <cstddef>

namespace accumulatr::eval::detail {

void exp_lanes(const double *input,
               double *output,
               std::size_t size) noexcept;

void log_lanes(const double *input,
               double *output,
               std::size_t size) noexcept;

} // namespace accumulatr::eval::detail
