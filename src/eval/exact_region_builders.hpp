#pragma once

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include "exact_region_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline ExactOrderRegionExpr exact_order_region_zero() {
  return {};
}

inline ExactOrderRegionExpr exact_order_region_one() {
  ExactOrderRegionExpr out;
  out.terms.push_back(ExactRegionCell{});
  return out;
}

inline ExactRegionVar exact_region_source_time_var(
    const semantic::Index source_id) noexcept {
  return ExactRegionVar{ExactRegionVarKind::SourceTime, source_id};
}

inline ExactRegionVar exact_region_expr_time_var(
    const semantic::Index expr_id) noexcept {
  return ExactRegionVar{ExactRegionVarKind::ExprTime, expr_id};
}

inline ExactRegionVar exact_region_time_var(
    const semantic::Index time_id) noexcept {
  return ExactRegionVar{ExactRegionVarKind::Time, time_id};
}

inline bool exact_region_var_less(
    const ExactRegionVar &lhs,
    const ExactRegionVar &rhs) noexcept {
  if (lhs.kind != rhs.kind) {
    return lhs.kind < rhs.kind;
  }
  return lhs.id < rhs.id;
}

inline bool exact_region_var_equal(
    const ExactRegionVar &lhs,
    const ExactRegionVar &rhs) noexcept {
  return lhs.kind == rhs.kind && lhs.id == rhs.id;
}

inline ExactRegionAtom exact_region_source_atom(
    const ExactRegionAtomKind kind,
    const semantic::Index source_id,
    const semantic::Index time_id,
    const bool inclusive = false) {
  return ExactRegionAtom{
      kind,
      exact_region_source_time_var(source_id),
      exact_region_time_var(time_id),
      {},
      inclusive,
      true};
}

inline ExactRegionAtom exact_region_expr_atom(
    const ExactRegionAtomKind kind,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool inclusive = false) {
  return ExactRegionAtom{
      kind,
      exact_region_expr_time_var(expr_id),
      exact_region_time_var(time_id),
      {},
      inclusive,
      true};
}

inline ExactRegionAtom exact_region_time_order_atom(
    const semantic::Index before_time_id,
    const semantic::Index after_time_id,
    const bool strict) {
  return ExactRegionAtom{
      ExactRegionAtomKind::TimeOrder,
      exact_region_time_var(before_time_id),
      exact_region_time_var(after_time_id),
      {},
      false,
      strict};
}

inline ExactRegionAtom exact_region_outcome_atom(
    const ExactRegionAtomKind kind,
    std::vector<semantic::Index> outcome_indices) {
  outcome_indices.erase(
      std::remove(outcome_indices.begin(),
                  outcome_indices.end(),
                  semantic::kInvalidIndex),
      outcome_indices.end());
  std::sort(outcome_indices.begin(), outcome_indices.end());
  outcome_indices.erase(
      std::unique(outcome_indices.begin(), outcome_indices.end()),
      outcome_indices.end());
  return ExactRegionAtom{
      kind,
      ExactRegionVar{},
      ExactRegionVar{},
      std::move(outcome_indices),
      false,
      true};
}

inline bool exact_region_atom_less(
    const ExactRegionAtom &lhs,
    const ExactRegionAtom &rhs) noexcept {
  if (lhs.kind != rhs.kind) {
    return lhs.kind < rhs.kind;
  }
  if (!exact_region_var_equal(lhs.lhs, rhs.lhs)) {
    return exact_region_var_less(lhs.lhs, rhs.lhs);
  }
  if (!exact_region_var_equal(lhs.rhs, rhs.rhs)) {
    return exact_region_var_less(lhs.rhs, rhs.rhs);
  }
  if (lhs.outcome_indices != rhs.outcome_indices) {
    return lhs.outcome_indices < rhs.outcome_indices;
  }
  if (lhs.inclusive != rhs.inclusive) {
    return lhs.inclusive < rhs.inclusive;
  }
  return lhs.strict < rhs.strict;
}

inline bool exact_region_atom_equal(
    const ExactRegionAtom &lhs,
    const ExactRegionAtom &rhs) noexcept {
  return lhs.kind == rhs.kind &&
         exact_region_var_equal(lhs.lhs, rhs.lhs) &&
         exact_region_var_equal(lhs.rhs, rhs.rhs) &&
         lhs.outcome_indices == rhs.outcome_indices &&
         lhs.inclusive == rhs.inclusive &&
         lhs.strict == rhs.strict;
}

inline int exact_region_equality_origin_rank(
    const ExactRegionEqualityOrigin origin) noexcept {
  return static_cast<int>(origin);
}

inline bool exact_region_equality_positive(
    const ExactRegionEquality &equality) noexcept {
  return equality.mass == ExactRegionEqualityMass::PositiveMass ||
         equality.origin ==
             ExactRegionEqualityOrigin::SharedLatentIdentity;
}

inline ExactRegionEquality exact_region_normalized_equality(
    ExactRegionEquality equality) {
  if (equality.lhs_time_id > equality.rhs_time_id) {
    std::swap(equality.lhs_time_id, equality.rhs_time_id);
  }
  if (equality.origin ==
      ExactRegionEqualityOrigin::SharedLatentIdentity) {
    equality.mass = ExactRegionEqualityMass::PositiveMass;
  }
  return equality;
}

inline bool exact_region_equality_less(
    const ExactRegionEquality &lhs,
    const ExactRegionEquality &rhs) noexcept {
  if (lhs.lhs_time_id != rhs.lhs_time_id) {
    return lhs.lhs_time_id < rhs.lhs_time_id;
  }
  if (lhs.rhs_time_id != rhs.rhs_time_id) {
    return lhs.rhs_time_id < rhs.rhs_time_id;
  }
  if (lhs.mass != rhs.mass) {
    return lhs.mass < rhs.mass;
  }
  return lhs.origin < rhs.origin;
}

inline void exact_region_normalize_equalities(ExactRegionCell *cell) {
  std::vector<ExactRegionEquality> normalized;
  normalized.reserve(cell->equalities.size());
  for (auto equality : cell->equalities) {
    if (equality.lhs_time_id == semantic::kInvalidIndex ||
        equality.rhs_time_id == semantic::kInvalidIndex) {
      continue;
    }
    equality = exact_region_normalized_equality(equality);
    if (equality.lhs_time_id == equality.rhs_time_id &&
        equality.mass == ExactRegionEqualityMass::MeasureZero) {
      continue;
    }
    normalized.push_back(equality);
  }
  std::sort(
      normalized.begin(), normalized.end(), exact_region_equality_less);
  std::vector<ExactRegionEquality> unique;
  unique.reserve(normalized.size());
  for (const auto &equality : normalized) {
    if (!unique.empty() &&
        unique.back().lhs_time_id == equality.lhs_time_id &&
        unique.back().rhs_time_id == equality.rhs_time_id) {
      auto &back = unique.back();
      if (exact_region_equality_positive(equality)) {
        back.mass = ExactRegionEqualityMass::PositiveMass;
      }
      if (exact_region_equality_origin_rank(equality.origin) >
          exact_region_equality_origin_rank(back.origin)) {
        back.origin = equality.origin;
      }
      continue;
    }
    unique.push_back(equality);
  }
  cell->equalities = std::move(unique);
}

inline void exact_region_append_equality(
    ExactRegionCell *cell,
    const semantic::Index lhs_time_id,
    const semantic::Index rhs_time_id,
    const ExactRegionEqualityMass mass,
    const ExactRegionEqualityOrigin origin) {
  if (cell->impossible ||
      lhs_time_id == semantic::kInvalidIndex ||
      rhs_time_id == semantic::kInvalidIndex) {
    return;
  }
  ExactRegionEquality equality{
      lhs_time_id, rhs_time_id, mass, origin};
  equality = exact_region_normalized_equality(equality);
  if (equality.lhs_time_id == equality.rhs_time_id &&
      equality.mass == ExactRegionEqualityMass::MeasureZero) {
    return;
  }
  cell->equalities.push_back(equality);
}

inline std::vector<ExactOrderRegionSourceTime>
exact_region_source_atoms(
    const ExactRegionCell &cell,
    const ExactRegionAtomKind kind) {
  std::vector<ExactOrderRegionSourceTime> out;
  for (const auto &atom : cell.atoms) {
    if (atom.kind == kind &&
        atom.lhs.kind == ExactRegionVarKind::SourceTime &&
        atom.rhs.kind == ExactRegionVarKind::Time) {
      out.push_back(
          ExactOrderRegionSourceTime{
              atom.lhs.id,
              atom.rhs.id,
              atom.inclusive});
    }
  }
  return out;
}

inline std::vector<ExactOrderRegionSourceTime>
exact_region_exact_source_atoms(const ExactRegionCell &cell) {
  return exact_region_source_atoms(cell, ExactRegionAtomKind::SourceExact);
}

inline std::vector<ExactOrderRegionSourceTime>
exact_region_lower_source_atoms(const ExactRegionCell &cell) {
  return exact_region_source_atoms(cell, ExactRegionAtomKind::SourceLower);
}

inline std::vector<ExactOrderRegionSourceTime>
exact_region_upper_source_atoms(const ExactRegionCell &cell) {
  return exact_region_source_atoms(cell, ExactRegionAtomKind::SourceUpper);
}

inline std::vector<ExactOrderRegionExprValueFactor>
exact_region_expr_atoms(const ExactRegionCell &cell) {
  std::vector<ExactOrderRegionExprValueFactor> out;
  for (const auto &atom : cell.atoms) {
    if (atom.lhs.kind != ExactRegionVarKind::ExprTime ||
        atom.rhs.kind != ExactRegionVarKind::Time) {
      continue;
    }
    if (atom.kind == ExactRegionAtomKind::ExprBefore) {
      out.push_back(
	          ExactOrderRegionExprValueFactor{
	              atom.lhs.id,
	              atom.rhs.id,
	              true,
	              atom.inclusive,
	              false});
	    } else if (atom.kind == ExactRegionAtomKind::ExprNotBefore) {
	      out.push_back(
	          ExactOrderRegionExprValueFactor{
	              atom.lhs.id,
	              atom.rhs.id,
	              false,
	              atom.inclusive,
	              false});
	    } else if (atom.kind == ExactRegionAtomKind::ExprDensity) {
	      out.push_back(
	          ExactOrderRegionExprValueFactor{
	              atom.lhs.id,
	              atom.rhs.id,
	              true,
	              false,
	              true});
    }
  }
  return out;
}

inline semantic::Index exact_region_expr_density_time(
    const ExactRegionCell &cell,
    const semantic::Index expr_id) {
  for (const auto &factor : exact_region_expr_atoms(cell)) {
    if (factor.expr_id == expr_id && factor.density) {
      return factor.time_id;
    }
  }
  return semantic::kInvalidIndex;
}

inline std::vector<ExactOrderRegionTimeOrder>
exact_region_time_order_atoms(const ExactRegionCell &cell) {
  std::vector<ExactOrderRegionTimeOrder> out;
  for (const auto &atom : cell.atoms) {
    if (atom.kind == ExactRegionAtomKind::TimeOrder &&
        atom.lhs.kind == ExactRegionVarKind::Time &&
        atom.rhs.kind == ExactRegionVarKind::Time) {
      out.push_back(
          ExactOrderRegionTimeOrder{atom.lhs.id, atom.rhs.id, atom.strict});
    }
  }
  return out;
}

inline std::vector<std::vector<semantic::Index>>
exact_region_outcome_atoms(
    const ExactRegionCell &cell,
    const ExactRegionAtomKind kind) {
  std::vector<std::vector<semantic::Index>> out;
  for (const auto &atom : cell.atoms) {
    if (atom.kind == kind && !atom.outcome_indices.empty()) {
      out.push_back(atom.outcome_indices);
    }
  }
  return out;
}

inline void exact_region_append_atom(
    ExactRegionCell *cell,
    ExactRegionAtom atom) {
  if (atom.kind == ExactRegionAtomKind::OutcomeUnused ||
      atom.kind == ExactRegionAtomKind::OutcomeUsed) {
    if (atom.outcome_indices.empty()) {
      return;
    }
  } else if (atom.lhs.id == semantic::kInvalidIndex ||
             atom.rhs.id == semantic::kInvalidIndex) {
    return;
  }
  cell->atoms.push_back(std::move(atom));
}

inline bool exact_region_time_is_latent_variable(
    const semantic::Index time_id) noexcept {
  return time_id >
         static_cast<semantic::Index>(CompiledMathTimeSlot::Zero);
}

inline semantic::Index exact_order_region_new_time(
    ExactOrderRegionBuilder *builder) {
  return builder->next_time_id++;
}

inline void exact_order_region_reserve_time(
    ExactOrderRegionBuilder *builder,
    const semantic::Index time_id) {
  if (builder == nullptr ||
      time_id == semantic::kInvalidIndex ||
      !exact_region_time_is_latent_variable(time_id)) {
    return;
  }
  if (builder->next_time_id <= time_id) {
    builder->next_time_id = time_id + 1U;
  }
}

inline bool exact_order_region_time_is_special(
    const semantic::Index time_id) noexcept {
  return time_id ==
             static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) ||
         time_id ==
             static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
}

} // namespace detail
} // namespace accumulatr::eval
