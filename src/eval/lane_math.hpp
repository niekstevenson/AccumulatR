#pragma once

#include <cstddef>

namespace accumulatr::eval::detail {

void log_lanes(const double *input,
               double *output,
               std::size_t size) noexcept;

} // namespace accumulatr::eval::detail
