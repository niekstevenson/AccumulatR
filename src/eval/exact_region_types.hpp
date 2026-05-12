#pragma once

#include <cstdint>
#include <vector>

#include "compiled_math_types.hpp"
#include "exact_common_types.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactOrderRegionTimeOrder {
  semantic::Index before_time_id{semantic::kInvalidIndex};
  semantic::Index after_time_id{semantic::kInvalidIndex};
  bool strict{true};
};

struct ExactOrderRegionSourceTime {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index time_id{semantic::kInvalidIndex};
  bool inclusive{false};
};

struct ExactOrderRegionExprValueFactor {
  semantic::Index expr_id{semantic::kInvalidIndex};
  semantic::Index time_id{semantic::kInvalidIndex};
  bool before{true};
  bool inclusive{false};
  bool density{false};
};

enum class ExactRegionVarKind : std::uint8_t {
  SourceTime = 0,
  ExprTime = 1,
  Time = 2
};

struct ExactRegionVar {
  ExactRegionVarKind kind{ExactRegionVarKind::Time};
  semantic::Index id{semantic::kInvalidIndex};
};

enum class ExactRegionAtomKind : std::uint8_t {
  SourceExact = 0,
  SourceLower = 1,
  SourceUpper = 2,
  ExprBefore = 3,
  ExprNotBefore = 4,
  ExprDensity = 5,
  TimeOrder = 6,
  OutcomeUnused = 7,
  OutcomeUsed = 8
};

enum class ExactRegionEqualityMass : std::uint8_t {
  MeasureZero = 0,
  PositiveMass = 1
};

enum class ExactRegionEqualityOrigin : std::uint8_t {
  ContinuousBoundary = 0,
  ModelTie = 1,
  SharedLatentIdentity = 2
};

struct ExactRegionAtom {
  ExactRegionAtomKind kind{ExactRegionAtomKind::SourceLower};
  ExactRegionVar lhs{};
  ExactRegionVar rhs{};
  std::vector<semantic::Index> outcome_indices;
  bool inclusive{false};
  bool strict{true};
};

struct ExactRegionEquality {
  semantic::Index lhs_time_id{semantic::kInvalidIndex};
  semantic::Index rhs_time_id{semantic::kInvalidIndex};
  ExactRegionEqualityMass mass{ExactRegionEqualityMass::MeasureZero};
  ExactRegionEqualityOrigin origin{
      ExactRegionEqualityOrigin::ContinuousBoundary};
};

struct ExactRegionCell {
  double sign{1.0};
  std::vector<ExactRegionAtom> atoms;
  std::vector<ExactRegionEquality> equalities;
  bool impossible{false};
};

struct ExactOrderRegionExpr {
  std::vector<ExactRegionCell> terms;
};

struct ExactOrderRegionBuilder {
  semantic::Index next_time_id{
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) + 1U};
};

struct ExactOrderRegionTimeClosure {
  std::vector<semantic::Index> time_ids;
  std::vector<std::uint8_t> relation;
  std::vector<semantic::Index> representative;
  bool impossible{false};
};

} // namespace detail
} // namespace accumulatr::eval
