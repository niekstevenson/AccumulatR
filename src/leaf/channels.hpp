#pragma once

#include <cmath>

namespace accumulatr::leaf {

struct EventChannels {
  double pdf{0.0};
  double cdf{0.0};
  double survival{1.0};

  static constexpr EventChannels impossible() noexcept {
    return EventChannels{0.0, 0.0, 1.0};
  }

  static constexpr EventChannels certain() noexcept {
    return EventChannels{0.0, 1.0, 0.0};
  }

  bool finite() const noexcept {
    return std::isfinite(pdf) && std::isfinite(cdf) && std::isfinite(survival);
  }
};

struct LeafChannels {
  EventChannels event{};
  EventChannels given_start{};
  bool has_given_start{false};
};

} // namespace accumulatr::leaf
