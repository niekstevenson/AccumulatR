#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "../semantic/model.hpp"

namespace accumulatr::eval {
namespace detail {

enum class ExactRelation : std::uint8_t {
  Unknown = 0,
  Before = 1,
  At = 2,
  After = 3
};

struct ExactIndexSpan {
  semantic::Index offset{0};
  semantic::Index size{0};

  [[nodiscard]] bool empty() const noexcept {
    return size == 0;
  }
};

struct ExactSourceOrderFact {
  semantic::Index before_source_id{semantic::kInvalidIndex};
  semantic::Index after_source_id{semantic::kInvalidIndex};
};

struct ExactTimedExprUpperBound {
  semantic::Index expr_id{semantic::kInvalidIndex};
  double time{std::numeric_limits<double>::quiet_NaN()};
  double normalizer{0.0};
};

struct ExactRelationTemplate {
  std::vector<semantic::Index> source_ids;
  std::vector<ExactRelation> relations;

  bool empty() const noexcept {
    return source_ids.empty();
  }
};

inline double clean_signed_value(const double value,
                                 const double eps = 1e-15) {
  if (!std::isfinite(value)) {
    return 0.0;
  }
  return std::fabs(value) <= eps ? 0.0 : value;
}

} // namespace detail
} // namespace accumulatr::eval
