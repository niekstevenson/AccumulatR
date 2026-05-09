#pragma once

#include "exact_types.hpp"

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

inline bool exact_order_region_time_is_special(
    const semantic::Index time_id) noexcept {
  return time_id ==
             static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) ||
         time_id ==
             static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
}

inline void exact_order_region_append_time_id(
    std::vector<semantic::Index> *time_ids,
    const semantic::Index time_id) {
  if (time_id == semantic::kInvalidIndex) {
    return;
  }
  if (std::find(time_ids->begin(), time_ids->end(), time_id) ==
      time_ids->end()) {
    time_ids->push_back(time_id);
  }
}

inline std::size_t exact_order_region_time_pos(
    const ExactOrderRegionTimeClosure &closure,
    const semantic::Index time_id) {
  const auto it =
      std::find(closure.time_ids.begin(), closure.time_ids.end(), time_id);
  return it == closure.time_ids.end()
             ? closure.time_ids.size()
             : static_cast<std::size_t>(it - closure.time_ids.begin());
}

inline std::uint8_t &exact_order_region_time_relation_slot(
    ExactOrderRegionTimeClosure *closure,
    const std::size_t before_pos,
    const std::size_t after_pos) {
  return closure->relation[
      before_pos * closure->time_ids.size() + after_pos];
}

inline std::uint8_t exact_order_region_time_relation_at(
    const ExactOrderRegionTimeClosure &closure,
    const std::size_t before_pos,
    const std::size_t after_pos) {
  return closure.relation[
      before_pos * closure.time_ids.size() + after_pos];
}

inline void exact_order_region_set_time_relation(
    ExactOrderRegionTimeClosure *closure,
    const semantic::Index before_time_id,
    const semantic::Index after_time_id,
    const bool strict) {
  const auto before_pos =
      exact_order_region_time_pos(*closure, before_time_id);
  const auto after_pos =
      exact_order_region_time_pos(*closure, after_time_id);
  if (before_pos >= closure->time_ids.size() ||
      after_pos >= closure->time_ids.size()) {
    return;
  }
  auto &slot =
      exact_order_region_time_relation_slot(closure, before_pos, after_pos);
  const std::uint8_t value = strict ? 2U : 1U;
  if (slot < value) {
    slot = value;
  }
}

inline ExactOrderRegionTimeClosure exact_order_region_build_time_closure(
    const ExactRegionCell &term) {
  ExactOrderRegionTimeClosure closure;
  exact_order_region_append_time_id(
      &closure.time_ids,
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero));
  for (const auto &item : exact_region_exact_source_atoms(term)) {
    exact_order_region_append_time_id(&closure.time_ids, item.time_id);
  }
  for (const auto &item : exact_region_lower_source_atoms(term)) {
    exact_order_region_append_time_id(&closure.time_ids, item.time_id);
  }
  for (const auto &item : exact_region_upper_source_atoms(term)) {
    exact_order_region_append_time_id(&closure.time_ids, item.time_id);
  }
  for (const auto &factor : exact_region_expr_atoms(term)) {
    exact_order_region_append_time_id(&closure.time_ids, factor.time_id);
  }
  for (const auto &order : exact_region_time_order_atoms(term)) {
    exact_order_region_append_time_id(&closure.time_ids, order.before_time_id);
    exact_order_region_append_time_id(&closure.time_ids, order.after_time_id);
  }
  for (const auto &equality : term.equalities) {
    exact_order_region_append_time_id(&closure.time_ids, equality.lhs_time_id);
    exact_order_region_append_time_id(&closure.time_ids, equality.rhs_time_id);
  }
  std::sort(closure.time_ids.begin(), closure.time_ids.end());
  const auto n = closure.time_ids.size();
  closure.relation.assign(n * n, 0U);
  closure.representative = closure.time_ids;
  for (std::size_t i = 0; i < n; ++i) {
    exact_order_region_time_relation_slot(&closure, i, i) = 1U;
  }
  const auto zero_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero);
  for (const auto time_id : closure.time_ids) {
    exact_order_region_set_time_relation(
        &closure, zero_time_id, time_id, false);
  }
  for (const auto &equality : term.equalities) {
    if (equality.origin !=
            ExactRegionEqualityOrigin::SharedLatentIdentity ||
        equality.lhs_time_id == equality.rhs_time_id) {
      continue;
    }
    exact_order_region_set_time_relation(
        &closure, equality.lhs_time_id, equality.rhs_time_id, false);
    exact_order_region_set_time_relation(
        &closure, equality.rhs_time_id, equality.lhs_time_id, false);
  }
  for (const auto &order : exact_region_time_order_atoms(term)) {
    if (order.before_time_id == order.after_time_id) {
      if (order.strict) {
        closure.impossible = true;
      }
      continue;
    }
    exact_order_region_set_time_relation(
        &closure, order.before_time_id, order.after_time_id, order.strict);
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      const auto ik = exact_order_region_time_relation_at(closure, i, k);
      if (ik == 0U) {
        continue;
      }
      for (std::size_t j = 0; j < n; ++j) {
        const auto kj = exact_order_region_time_relation_at(closure, k, j);
        if (kj == 0U) {
          continue;
        }
        const std::uint8_t composed = (ik == 2U || kj == 2U) ? 2U : 1U;
        auto &ij = exact_order_region_time_relation_slot(&closure, i, j);
        if (ij < composed) {
          ij = composed;
        }
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (exact_order_region_time_relation_at(closure, i, i) == 2U) {
      closure.impossible = true;
      return closure;
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    semantic::Index rep = closure.time_ids[i];
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j ||
          exact_order_region_time_is_special(closure.time_ids[i]) ||
          exact_order_region_time_is_special(closure.time_ids[j])) {
        continue;
      }
      const auto ij = exact_order_region_time_relation_at(closure, i, j);
      const auto ji = exact_order_region_time_relation_at(closure, j, i);
      if (ij == 1U && ji == 1U) {
        rep = std::min(rep, closure.time_ids[j]);
      }
    }
    closure.representative[i] = rep;
  }
  bool changed = true;
  while (changed) {
    changed = false;
    for (std::size_t i = 0; i < n; ++i) {
      auto rep = closure.representative[i];
      const auto rep_pos = exact_order_region_time_pos(closure, rep);
      if (rep_pos < n && closure.representative[rep_pos] != rep) {
        rep = closure.representative[rep_pos];
      }
      if (closure.representative[i] != rep) {
        closure.representative[i] = rep;
        changed = true;
      }
    }
  }
  for (const auto &equality : term.equalities) {
    const auto lhs_pos =
        exact_order_region_time_pos(closure, equality.lhs_time_id);
    const auto rhs_pos =
        exact_order_region_time_pos(closure, equality.rhs_time_id);
    const auto lhs = lhs_pos >= n
                         ? equality.lhs_time_id
                         : closure.representative[lhs_pos];
    const auto rhs = rhs_pos >= n
                         ? equality.rhs_time_id
                         : closure.representative[rhs_pos];
    if (lhs == rhs) {
      continue;
    }
    if (equality.origin ==
        ExactRegionEqualityOrigin::ModelTie) {
      const auto lhs_original_pos =
          exact_order_region_time_pos(closure, equality.lhs_time_id);
      const auto rhs_original_pos =
          exact_order_region_time_pos(closure, equality.rhs_time_id);
      const bool strict_before =
          lhs_original_pos < n &&
          rhs_original_pos < n &&
          exact_order_region_time_relation_at(
              closure, lhs_original_pos, rhs_original_pos) == 2U;
      const bool strict_after =
          lhs_original_pos < n &&
          rhs_original_pos < n &&
          exact_order_region_time_relation_at(
              closure, rhs_original_pos, lhs_original_pos) == 2U;
      if (!strict_before && !strict_after) {
        continue;
      }
    }
    closure.impossible = true;
    return closure;
  }
  return closure;
}

inline semantic::Index exact_order_region_canonical_time(
    const ExactOrderRegionTimeClosure &closure,
    const semantic::Index time_id) {
  const auto pos = exact_order_region_time_pos(closure, time_id);
  return pos >= closure.time_ids.size() ? time_id : closure.representative[pos];
}

inline int exact_order_region_time_relation(
    const ExactRegionCell &term,
    const semantic::Index before_time_id,
    const semantic::Index after_time_id) {
  if (before_time_id == semantic::kInvalidIndex ||
      after_time_id == semantic::kInvalidIndex) {
    return 0;
  }
  const auto closure = exact_order_region_build_time_closure(term);
  const auto before_pos =
      exact_order_region_time_pos(closure, before_time_id);
  const auto after_pos =
      exact_order_region_time_pos(closure, after_time_id);
  if (before_pos >= closure.time_ids.size() ||
      after_pos >= closure.time_ids.size()) {
    return 0;
  }
  return static_cast<int>(
      exact_order_region_time_relation_at(closure, before_pos, after_pos));
}

inline bool exact_order_region_time_known_before_or_equal(
    const ExactRegionCell &term,
    const semantic::Index before_time_id,
    const semantic::Index after_time_id) {
  return exact_order_region_time_relation(
             term, before_time_id, after_time_id) != 0;
}

inline void exact_order_region_append_time_order(
    ExactRegionCell *term,
    const semantic::Index before_time_id,
    const semantic::Index after_time_id,
    const bool strict = true) {
  if (term->impossible ||
      before_time_id == semantic::kInvalidIndex ||
      after_time_id == semantic::kInvalidIndex) {
    return;
  }
  if (before_time_id == after_time_id) {
    if (strict) {
      term->impossible = true;
    }
    return;
  }
  const auto forward =
      exact_order_region_time_relation(*term, before_time_id, after_time_id);
  if (forward == 2 || (forward == 1 && !strict)) {
    return;
  }
  const auto reverse =
      exact_order_region_time_relation(*term, after_time_id, before_time_id);
  if (reverse != 0) {
    if (!strict && reverse == 1) {
      exact_region_append_equality(
          term,
          before_time_id,
          after_time_id,
          ExactRegionEqualityMass::MeasureZero,
          ExactRegionEqualityOrigin::ContinuousBoundary);
    }
    term->impossible = true;
    return;
  }
  exact_region_append_atom(
      term,
      exact_region_time_order_atom(before_time_id, after_time_id, strict));
}

inline semantic::Index exact_order_region_exact_time_for_source(
    const ExactRegionCell &term,
    const semantic::Index source_id) {
  for (const auto &exact : exact_region_exact_source_atoms(term)) {
    if (exact.source_id == source_id) {
      return exact.time_id;
    }
  }
  return semantic::kInvalidIndex;
}

inline bool exact_order_region_contains_source_time(
    const std::vector<ExactOrderRegionSourceTime> &items,
    const semantic::Index source_id,
    const semantic::Index time_id,
    const bool inclusive) {
  return std::any_of(
      items.begin(),
      items.end(),
      [&](const ExactOrderRegionSourceTime &item) {
        return item.source_id == source_id &&
               item.time_id == time_id &&
               item.inclusive == inclusive;
      });
}

inline void exact_order_region_remove_source_bounds(
    ExactRegionCell *term,
    const semantic::Index source_id) {
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return (atom.kind == ExactRegionAtomKind::SourceLower ||
                    atom.kind == ExactRegionAtomKind::SourceUpper) &&
                   atom.lhs.kind == ExactRegionVarKind::SourceTime &&
                   atom.lhs.id == source_id;
          }),
      term->atoms.end());
}

inline void exact_order_region_remove_expr_value_atoms(
    ExactRegionCell *term,
    const semantic::Index expr_id) {
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return (atom.kind == ExactRegionAtomKind::ExprBefore ||
                    atom.kind == ExactRegionAtomKind::ExprNotBefore) &&
                   atom.lhs.kind == ExactRegionVarKind::ExprTime &&
                   atom.lhs.id == expr_id;
          }),
	      term->atoms.end());
}

inline bool exact_order_region_remove_expr_value_atom(
    ExactRegionCell *term,
    const ExactOrderRegionExprValueFactor &factor) {
  if (factor.density) {
    return false;
  }
  const auto kind =
      factor.before ? ExactRegionAtomKind::ExprBefore
                    : ExactRegionAtomKind::ExprNotBefore;
  const auto before = term->atoms.size();
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return atom.kind == kind &&
                   atom.lhs.kind == ExactRegionVarKind::ExprTime &&
                   atom.lhs.id == factor.expr_id &&
                   atom.rhs.kind == ExactRegionVarKind::Time &&
                   atom.rhs.id == factor.time_id &&
                   atom.inclusive == factor.inclusive;
          }),
      term->atoms.end());
  return term->atoms.size() != before;
}

inline void exact_order_region_replace_time(
    ExactRegionCell *term,
    const semantic::Index from_time_id,
    const semantic::Index to_time_id);

inline void exact_order_region_append_exact(
    ExactRegionCell *term,
    const semantic::Index source_id,
    const semantic::Index time_id) {
  if (term->impossible ||
      source_id == semantic::kInvalidIndex ||
      time_id == semantic::kInvalidIndex) {
    return;
  }
  const auto existing_time =
      exact_order_region_exact_time_for_source(*term, source_id);
  if (existing_time != semantic::kInvalidIndex) {
    if (existing_time != time_id) {
      exact_region_append_equality(
          term,
          existing_time,
          time_id,
          ExactRegionEqualityMass::PositiveMass,
          ExactRegionEqualityOrigin::SharedLatentIdentity);
      exact_order_region_replace_time(term, time_id, existing_time);
    }
    return;
  }
  exact_region_append_atom(
      term,
      exact_region_source_atom(
          ExactRegionAtomKind::SourceExact, source_id, time_id));
  for (const auto &bound : exact_region_lower_source_atoms(*term)) {
    if (bound.source_id == source_id) {
      if (bound.time_id == time_id && !bound.inclusive) {
        term->impossible = true;
        return;
      }
      exact_order_region_append_time_order(
          term, bound.time_id, time_id, !bound.inclusive);
    }
  }
  for (const auto &bound : exact_region_upper_source_atoms(*term)) {
    if (bound.source_id == source_id) {
      if (bound.time_id == time_id) {
        term->impossible = true;
        return;
      }
      exact_order_region_append_time_order(
          term, time_id, bound.time_id, !bound.inclusive);
    }
  }
  exact_order_region_remove_source_bounds(term, source_id);
}

inline void exact_order_region_append_lower(
    ExactRegionCell *term,
    const semantic::Index source_id,
    const semantic::Index time_id,
    const bool inclusive = false) {
  if (term->impossible ||
      source_id == semantic::kInvalidIndex ||
      time_id == semantic::kInvalidIndex) {
    return;
  }
  const auto exact_time =
      exact_order_region_exact_time_for_source(*term, source_id);
  if (exact_time != semantic::kInvalidIndex) {
    if (exact_time == time_id && !inclusive) {
      term->impossible = true;
      return;
    }
    exact_order_region_append_time_order(
        term, time_id, exact_time, !inclusive);
    return;
  }
  const auto upper = exact_region_upper_source_atoms(*term);
  for (const auto &bound : upper) {
    if (bound.source_id == source_id && bound.time_id == time_id) {
      term->impossible = true;
      return;
    }
  }
  const auto lower = exact_region_lower_source_atoms(*term);
  if (!exact_order_region_contains_source_time(
          lower, source_id, time_id, inclusive)) {
    exact_region_append_atom(
        term,
        exact_region_source_atom(
            ExactRegionAtomKind::SourceLower, source_id, time_id, inclusive));
  }
}

inline void exact_order_region_append_upper(
    ExactRegionCell *term,
    const semantic::Index source_id,
    const semantic::Index time_id,
    const bool inclusive = false) {
  if (term->impossible ||
      source_id == semantic::kInvalidIndex ||
      time_id == semantic::kInvalidIndex) {
    return;
  }
  const auto exact_time =
      exact_order_region_exact_time_for_source(*term, source_id);
  if (exact_time != semantic::kInvalidIndex) {
    if (exact_time == time_id && !inclusive) {
      term->impossible = true;
      return;
    }
    exact_order_region_append_time_order(
        term, exact_time, time_id, !inclusive);
    return;
  }
  const auto lower = exact_region_lower_source_atoms(*term);
  for (const auto &bound : lower) {
    if (bound.source_id == source_id && bound.time_id == time_id) {
      term->impossible = true;
      return;
    }
  }
  const auto upper = exact_region_upper_source_atoms(*term);
  if (!exact_order_region_contains_source_time(
          upper, source_id, time_id, inclusive)) {
    exact_region_append_atom(
        term,
        exact_region_source_atom(
            ExactRegionAtomKind::SourceUpper, source_id, time_id, inclusive));
  }
}

inline void exact_order_region_append_expr_factor(
    ExactRegionCell *term,
    const semantic::Index expr_id,
    const semantic::Index time_id,
    const bool before,
    const bool density = false,
    const bool inclusive = false) {
  if (term->impossible ||
      expr_id == semantic::kInvalidIndex ||
      time_id == semantic::kInvalidIndex) {
    return;
  }
  const auto existing_density_time =
      exact_region_expr_density_time(*term, expr_id);
  if (density) {
    if (existing_density_time != semantic::kInvalidIndex) {
      if (existing_density_time != time_id) {
        exact_region_append_equality(
            term,
            existing_density_time,
            time_id,
            ExactRegionEqualityMass::PositiveMass,
            ExactRegionEqualityOrigin::SharedLatentIdentity);
        exact_order_region_replace_time(term, time_id, existing_density_time);
      }
      return;
    }
    const auto existing_factors = exact_region_expr_atoms(*term);
    for (const auto &factor : existing_factors) {
      if (factor.expr_id != expr_id || factor.density) {
        continue;
      }
      if (factor.before) {
        exact_order_region_append_time_order(
            term, time_id, factor.time_id, !factor.inclusive);
      } else {
        exact_order_region_append_time_order(
            term, factor.time_id, time_id, !factor.inclusive);
      }
      if (term->impossible) {
        return;
      }
    }
    exact_order_region_remove_expr_value_atoms(term, expr_id);
    exact_region_append_atom(
        term,
        exact_region_expr_atom(
            ExactRegionAtomKind::ExprDensity, expr_id, time_id));
    return;
  }
  if (existing_density_time != semantic::kInvalidIndex) {
    if (before) {
      exact_order_region_append_time_order(
          term, existing_density_time, time_id, !inclusive);
    } else {
      exact_order_region_append_time_order(
          term, time_id, existing_density_time, !inclusive);
    }
    return;
  }
  for (const auto &factor : exact_region_expr_atoms(*term)) {
    if (factor.expr_id == expr_id &&
        factor.time_id == time_id &&
        !factor.density) {
      if (factor.before != before) {
        term->impossible = true;
        return;
      }
      if (factor.inclusive == inclusive) {
        return;
      }
	    }
  }
  exact_region_append_atom(
      term,
      exact_region_expr_atom(
          before ? ExactRegionAtomKind::ExprBefore
                 : ExactRegionAtomKind::ExprNotBefore,
          expr_id,
          time_id,
          inclusive));
}

inline void exact_order_region_append_outcome_gate(
    ExactRegionCell *term,
    const std::vector<semantic::Index> &outcome_indices) {
  exact_region_append_atom(
      term,
      exact_region_outcome_atom(
          ExactRegionAtomKind::OutcomeUnused, outcome_indices));
}

inline void exact_order_region_append_outcome_used_gate(
    ExactRegionCell *term,
    const std::vector<semantic::Index> &outcome_indices) {
  exact_region_append_atom(
      term,
      exact_region_outcome_atom(
          ExactRegionAtomKind::OutcomeUsed, outcome_indices));
}

inline void exact_order_region_replace_time(
    ExactRegionCell *term,
    const semantic::Index from_time_id,
    const semantic::Index to_time_id) {
  if (from_time_id == to_time_id ||
      from_time_id == semantic::kInvalidIndex ||
      to_time_id == semantic::kInvalidIndex) {
    return;
  }
  for (auto &atom : term->atoms) {
    if (atom.lhs.kind == ExactRegionVarKind::Time &&
        atom.lhs.id == from_time_id) {
      atom.lhs.id = to_time_id;
    }
    if (atom.rhs.kind == ExactRegionVarKind::Time &&
        atom.rhs.id == from_time_id) {
      atom.rhs.id = to_time_id;
    }
  }
  for (auto &equality : term->equalities) {
    if (equality.lhs_time_id == from_time_id) {
      equality.lhs_time_id = to_time_id;
    }
    if (equality.rhs_time_id == from_time_id) {
      equality.rhs_time_id = to_time_id;
    }
  }
  exact_region_normalize_equalities(term);
}

inline ExactRegionCell exact_order_region_intersect_terms(
    const ExactRegionCell &lhs,
    const ExactRegionCell &rhs) {
  ExactRegionCell out = lhs;
  out.sign *= rhs.sign;
  if (out.impossible || rhs.impossible) {
    out.impossible = true;
    return out;
  }
  ExactRegionCell rhs_term = rhs;
  for (const auto &exact : exact_region_exact_source_atoms(rhs_term)) {
    const auto existing_time =
        exact_order_region_exact_time_for_source(out, exact.source_id);
    if (existing_time != semantic::kInvalidIndex &&
        existing_time != exact.time_id) {
      exact_region_append_equality(
          &out,
          existing_time,
          exact.time_id,
          ExactRegionEqualityMass::PositiveMass,
          ExactRegionEqualityOrigin::SharedLatentIdentity);
      exact_order_region_replace_time(
          &rhs_term, exact.time_id, existing_time);
    }
  }
  for (const auto &factor : exact_region_expr_atoms(rhs_term)) {
    if (!factor.density) {
      continue;
    }
    const auto existing_time =
        exact_region_expr_density_time(out, factor.expr_id);
    if (existing_time != semantic::kInvalidIndex &&
        existing_time != factor.time_id) {
      exact_region_append_equality(
          &out,
          existing_time,
          factor.time_id,
          ExactRegionEqualityMass::PositiveMass,
          ExactRegionEqualityOrigin::SharedLatentIdentity);
      exact_order_region_replace_time(
          &rhs_term, factor.time_id, existing_time);
    }
  }
  for (const auto &equality : rhs_term.equalities) {
    exact_region_append_equality(
        &out,
        equality.lhs_time_id,
        equality.rhs_time_id,
        equality.mass,
        equality.origin);
  }
  for (const auto &order : exact_region_time_order_atoms(rhs_term)) {
    exact_order_region_append_time_order(
        &out, order.before_time_id, order.after_time_id, order.strict);
  }
  for (const auto &exact : exact_region_exact_source_atoms(rhs_term)) {
    exact_order_region_append_exact(&out, exact.source_id, exact.time_id);
  }
  for (const auto &bound : exact_region_lower_source_atoms(rhs_term)) {
    exact_order_region_append_lower(
        &out, bound.source_id, bound.time_id, bound.inclusive);
  }
  for (const auto &bound : exact_region_upper_source_atoms(rhs_term)) {
    exact_order_region_append_upper(
        &out, bound.source_id, bound.time_id, bound.inclusive);
  }
	  for (const auto &factor : exact_region_expr_atoms(rhs_term)) {
	    exact_order_region_append_expr_factor(
	        &out,
	        factor.expr_id,
	        factor.time_id,
	        factor.before,
	        factor.density,
	        factor.inclusive);
	  }
  for (const auto &outcome_indices :
       exact_region_outcome_atoms(rhs_term, ExactRegionAtomKind::OutcomeUnused)) {
    exact_order_region_append_outcome_gate(&out, outcome_indices);
  }
  for (const auto &outcome_indices :
       exact_region_outcome_atoms(rhs_term, ExactRegionAtomKind::OutcomeUsed)) {
    exact_order_region_append_outcome_used_gate(&out, outcome_indices);
  }
  return out;
}

inline ExactOrderRegionExpr exact_order_region_conjoin(
    const ExactOrderRegionExpr &lhs,
    const ExactOrderRegionExpr &rhs) {
  ExactOrderRegionExpr out;
  out.terms.reserve(lhs.terms.size() * rhs.terms.size());
  for (const auto &lhs_term : lhs.terms) {
    for (const auto &rhs_term : rhs.terms) {
      auto term = exact_order_region_intersect_terms(lhs_term, rhs_term);
      if (!term.impossible && term.sign != 0.0) {
        out.terms.push_back(std::move(term));
      }
    }
  }
  return out;
}

inline void exact_order_region_append_expr(ExactOrderRegionExpr *dst,
                                           ExactOrderRegionExpr src) {
  dst->terms.insert(
      dst->terms.end(),
      std::make_move_iterator(src.terms.begin()),
      std::make_move_iterator(src.terms.end()));
}

inline ExactOrderRegionExpr exact_order_region_with_outcome_gate(
    ExactOrderRegionExpr expr,
    const std::vector<semantic::Index> &outcome_indices) {
  if (outcome_indices.empty()) {
    return expr;
  }
  for (auto &term : expr.terms) {
    exact_order_region_append_outcome_gate(&term, outcome_indices);
  }
  return expr;
}

inline ExactOrderRegionExpr exact_order_region_with_outcome_used_gate(
    ExactOrderRegionExpr expr,
    const std::vector<semantic::Index> &outcome_indices) {
  if (outcome_indices.empty()) {
    return expr;
  }
  for (auto &term : expr.terms) {
    exact_order_region_append_outcome_used_gate(&term, outcome_indices);
  }
  return expr;
}

inline ExactOrderRegionExpr exact_order_region_union(
    const ExactOrderRegionExpr &lhs,
    const ExactOrderRegionExpr &rhs);

inline bool exact_order_region_source_time_less(
    const ExactOrderRegionSourceTime &lhs,
    const ExactOrderRegionSourceTime &rhs) noexcept {
  if (lhs.source_id != rhs.source_id) {
    return lhs.source_id < rhs.source_id;
  }
  if (lhs.time_id != rhs.time_id) {
    return lhs.time_id < rhs.time_id;
  }
  return lhs.inclusive < rhs.inclusive;
}

inline bool exact_order_region_expr_factor_less(
    const ExactOrderRegionExprValueFactor &lhs,
    const ExactOrderRegionExprValueFactor &rhs) noexcept {
  if (lhs.expr_id != rhs.expr_id) {
    return lhs.expr_id < rhs.expr_id;
  }
  if (lhs.time_id != rhs.time_id) {
    return lhs.time_id < rhs.time_id;
  }
  if (lhs.before != rhs.before) {
    return lhs.before < rhs.before;
  }
  if (lhs.inclusive != rhs.inclusive) {
    return lhs.inclusive < rhs.inclusive;
  }
  return lhs.density < rhs.density;
}

inline bool exact_order_region_time_order_less(
    const ExactOrderRegionTimeOrder &lhs,
    const ExactOrderRegionTimeOrder &rhs) noexcept {
  if (lhs.before_time_id != rhs.before_time_id) {
    return lhs.before_time_id < rhs.before_time_id;
  }
  if (lhs.after_time_id != rhs.after_time_id) {
    return lhs.after_time_id < rhs.after_time_id;
  }
  return lhs.strict < rhs.strict;
}

inline void exact_order_region_unique_source_times(
    std::vector<ExactOrderRegionSourceTime> *items) {
  std::sort(items->begin(), items->end(), exact_order_region_source_time_less);
  std::vector<ExactOrderRegionSourceTime> out;
  out.reserve(items->size());
  for (const auto &item : *items) {
    if (!out.empty() &&
        out.back().source_id == item.source_id &&
        out.back().time_id == item.time_id) {
      out.back().inclusive = out.back().inclusive && item.inclusive;
      continue;
    }
    out.push_back(item);
  }
  *items = std::move(out);
}

inline void exact_order_region_rewrite_term_times(
    ExactRegionCell *term,
    const ExactOrderRegionTimeClosure &closure) {
  for (auto &atom : term->atoms) {
    if (atom.lhs.kind == ExactRegionVarKind::Time) {
      atom.lhs.id = exact_order_region_canonical_time(closure, atom.lhs.id);
    }
    if (atom.rhs.kind == ExactRegionVarKind::Time) {
      atom.rhs.id = exact_order_region_canonical_time(closure, atom.rhs.id);
    }
  }
  for (auto &equality : term->equalities) {
    equality.lhs_time_id =
        exact_order_region_canonical_time(closure, equality.lhs_time_id);
    equality.rhs_time_id =
        exact_order_region_canonical_time(closure, equality.rhs_time_id);
  }
  exact_region_normalize_equalities(term);
}

inline std::uint8_t exact_order_region_canonical_relation(
    const ExactOrderRegionTimeClosure &closure,
    const semantic::Index before_time_id,
    const semantic::Index after_time_id) {
  std::uint8_t relation = 0U;
  for (std::size_t i = 0; i < closure.time_ids.size(); ++i) {
    if (closure.representative[i] != before_time_id) {
      continue;
    }
    for (std::size_t j = 0; j < closure.time_ids.size(); ++j) {
      if (closure.representative[j] != after_time_id) {
        continue;
      }
      relation = std::max(
          relation,
          exact_order_region_time_relation_at(closure, i, j));
    }
  }
  return relation;
}

inline std::vector<semantic::Index> exact_order_region_canonical_time_ids(
    const ExactOrderRegionTimeClosure &closure) {
  std::vector<semantic::Index> out;
  out.reserve(closure.representative.size());
  for (const auto time_id : closure.representative) {
    if (std::find(out.begin(), out.end(), time_id) == out.end()) {
      out.push_back(time_id);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

inline bool exact_order_region_canonical_order_redundant(
    const ExactOrderRegionTimeClosure &closure,
    const std::vector<semantic::Index> &time_ids,
    const semantic::Index before_time_id,
    const semantic::Index after_time_id,
    const std::uint8_t relation) {
  for (const auto middle_time_id : time_ids) {
    if (middle_time_id == before_time_id ||
        middle_time_id == after_time_id) {
      continue;
    }
    const auto left =
        exact_order_region_canonical_relation(
            closure, before_time_id, middle_time_id);
    const auto right =
        exact_order_region_canonical_relation(
            closure, middle_time_id, after_time_id);
    if (left == 0U || right == 0U) {
      continue;
    }
    const std::uint8_t composed =
        (left == 2U || right == 2U) ? 2U : 1U;
    if (composed >= relation) {
      return true;
    }
  }
  return false;
}

inline std::vector<ExactOrderRegionTimeOrder>
exact_order_region_canonical_time_orders(
    const ExactOrderRegionTimeClosure &closure) {
  std::vector<ExactOrderRegionTimeOrder> out;
  const auto time_ids = exact_order_region_canonical_time_ids(closure);
  const auto zero_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero);
  for (const auto before_time_id : time_ids) {
    if (before_time_id == zero_time_id) {
      continue;
    }
    for (const auto after_time_id : time_ids) {
      if (before_time_id == after_time_id || after_time_id == zero_time_id) {
        continue;
      }
      const auto relation =
          exact_order_region_canonical_relation(
              closure, before_time_id, after_time_id);
      if (relation == 0U) {
        continue;
      }
      if (exact_order_region_canonical_order_redundant(
              closure,
              time_ids,
              before_time_id,
              after_time_id,
              relation)) {
        continue;
      }
      out.push_back(
          ExactOrderRegionTimeOrder{
              before_time_id, after_time_id, relation == 2U});
    }
  }
  std::sort(out.begin(), out.end(), exact_order_region_time_order_less);
  return out;
}

inline bool exact_order_region_exact_sources_consistent(
    const ExactRegionCell &term) {
  std::vector<ExactOrderRegionSourceTime> exact =
      exact_region_exact_source_atoms(term);
  std::sort(exact.begin(), exact.end(), exact_order_region_source_time_less);
  for (std::size_t i = 1; i < exact.size(); ++i) {
    if (exact[i - 1U].source_id == exact[i].source_id &&
        exact[i - 1U].time_id != exact[i].time_id) {
      return false;
    }
  }
  return true;
}

inline bool exact_order_region_has_independent_source_density_tie(
    const ExactRegionCell &term) {
  std::vector<ExactOrderRegionSourceTime> exact =
      exact_region_exact_source_atoms(term);
  std::sort(
      exact.begin(),
      exact.end(),
      [](const ExactOrderRegionSourceTime &lhs,
         const ExactOrderRegionSourceTime &rhs) {
        if (lhs.time_id != rhs.time_id) {
          return lhs.time_id < rhs.time_id;
        }
        return lhs.source_id < rhs.source_id;
      });
  for (std::size_t i = 1; i < exact.size(); ++i) {
    if (exact[i - 1U].time_id == exact[i].time_id &&
        exact[i - 1U].source_id != exact[i].source_id) {
      return true;
    }
  }
  return false;
}

inline void exact_order_region_canonicalize_term(ExactRegionCell *term) {
  exact_region_normalize_equalities(term);
  auto closure = exact_order_region_build_time_closure(*term);
  if (closure.impossible) {
    term->impossible = true;
    return;
  }
  exact_order_region_rewrite_term_times(term, closure);
  closure = exact_order_region_build_time_closure(*term);
  if (closure.impossible ||
      !exact_order_region_exact_sources_consistent(*term) ||
      exact_order_region_has_independent_source_density_tie(*term)) {
    term->impossible = true;
    return;
  }
  const auto remove_source_atom =
      [&](const ExactRegionAtomKind kind,
          const semantic::Index source_id,
          const semantic::Index time_id) {
        term->atoms.erase(
            std::remove_if(
                term->atoms.begin(),
                term->atoms.end(),
                [&](const ExactRegionAtom &atom) {
                  return atom.kind == kind &&
                         atom.lhs.kind == ExactRegionVarKind::SourceTime &&
                         atom.lhs.id == source_id &&
                         atom.rhs.kind == ExactRegionVarKind::Time &&
                         atom.rhs.id == time_id;
                }),
            term->atoms.end());
      };

  for (const auto &source : exact_region_lower_source_atoms(*term)) {
    bool dominated = false;
    for (const auto &other : exact_region_lower_source_atoms(*term)) {
      if (source.source_id != other.source_id ||
          source.time_id == other.time_id) {
        continue;
      }
      if (exact_order_region_time_known_before_or_equal(
              *term, source.time_id, other.time_id)) {
        dominated = true;
        break;
      }
    }
    if (dominated) {
      remove_source_atom(
          ExactRegionAtomKind::SourceLower, source.source_id, source.time_id);
    }
  }
  for (const auto &source : exact_region_upper_source_atoms(*term)) {
    bool dominated = false;
    for (const auto &other : exact_region_upper_source_atoms(*term)) {
      if (source.source_id != other.source_id ||
          source.time_id == other.time_id) {
        continue;
      }
      if (exact_order_region_time_known_before_or_equal(
              *term, other.time_id, source.time_id)) {
        dominated = true;
        break;
      }
    }
    if (dominated) {
      remove_source_atom(
          ExactRegionAtomKind::SourceUpper, source.source_id, source.time_id);
    }
  }

  std::sort(term->atoms.begin(), term->atoms.end(), exact_region_atom_less);
  std::vector<ExactRegionAtom> unique_atoms;
  unique_atoms.reserve(term->atoms.size());
  for (const auto &atom : term->atoms) {
    if (!unique_atoms.empty() &&
        atom.kind == unique_atoms.back().kind &&
        exact_region_var_equal(atom.lhs, unique_atoms.back().lhs) &&
        exact_region_var_equal(atom.rhs, unique_atoms.back().rhs) &&
        atom.outcome_indices == unique_atoms.back().outcome_indices) {
      auto &back = unique_atoms.back();
      if (atom.kind == ExactRegionAtomKind::SourceLower ||
          atom.kind == ExactRegionAtomKind::SourceUpper) {
        back.inclusive = back.inclusive && atom.inclusive;
      } else if (atom.kind == ExactRegionAtomKind::TimeOrder) {
        back.strict = back.strict || atom.strict;
      }
      continue;
    }
    unique_atoms.push_back(atom);
  }

  std::vector<semantic::Index> unused_outcomes;
  std::vector<ExactRegionAtom> normalized_atoms;
  normalized_atoms.reserve(unique_atoms.size());
  for (const auto &atom : unique_atoms) {
    if (atom.kind == ExactRegionAtomKind::OutcomeUnused) {
      unused_outcomes.insert(
          unused_outcomes.end(),
          atom.outcome_indices.begin(),
          atom.outcome_indices.end());
      continue;
    }
    normalized_atoms.push_back(atom);
  }
  std::sort(unused_outcomes.begin(), unused_outcomes.end());
  unused_outcomes.erase(
      std::unique(unused_outcomes.begin(), unused_outcomes.end()),
      unused_outcomes.end());
  if (!unused_outcomes.empty()) {
    normalized_atoms.push_back(
        exact_region_outcome_atom(
            ExactRegionAtomKind::OutcomeUnused, unused_outcomes));
  }
  for (const auto &atom : normalized_atoms) {
    if (atom.kind != ExactRegionAtomKind::OutcomeUsed) {
      continue;
    }
    const bool fully_blocked =
        std::all_of(
            atom.outcome_indices.begin(),
            atom.outcome_indices.end(),
            [&](const semantic::Index outcome_id) {
              return std::binary_search(
                  unused_outcomes.begin(),
                  unused_outcomes.end(),
                  outcome_id);
            });
    if (fully_blocked) {
      term->impossible = true;
      return;
    }
  }
  std::sort(
      normalized_atoms.begin(),
      normalized_atoms.end(),
      exact_region_atom_less);
  term->atoms = std::move(normalized_atoms);
}

inline bool exact_order_region_remove_exact_source_time(
    ExactRegionCell *term,
    const semantic::Index source_id,
    const semantic::Index time_id) {
  const auto before = term->atoms.size();
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return atom.kind == ExactRegionAtomKind::SourceExact &&
                   atom.lhs.kind == ExactRegionVarKind::SourceTime &&
                   atom.lhs.id == source_id &&
                   atom.rhs.kind == ExactRegionVarKind::Time &&
                   atom.rhs.id == time_id;
          }),
      term->atoms.end());
  return term->atoms.size() != before;
}

inline bool exact_order_region_remove_lower_source_time(
    ExactRegionCell *term,
    const semantic::Index source_id,
    const semantic::Index time_id) {
  const auto before = term->atoms.size();
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return atom.kind == ExactRegionAtomKind::SourceLower &&
                   atom.lhs.kind == ExactRegionVarKind::SourceTime &&
                   atom.lhs.id == source_id &&
                   atom.rhs.kind == ExactRegionVarKind::Time &&
                   atom.rhs.id == time_id;
          }),
      term->atoms.end());
  return term->atoms.size() != before;
}

inline bool exact_order_region_remove_upper_source_time(
    ExactRegionCell *term,
    const semantic::Index source_id,
    const semantic::Index time_id) {
  const auto before = term->atoms.size();
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return atom.kind == ExactRegionAtomKind::SourceUpper &&
                   atom.lhs.kind == ExactRegionVarKind::SourceTime &&
                   atom.lhs.id == source_id &&
                   atom.rhs.kind == ExactRegionVarKind::Time &&
                   atom.rhs.id == time_id;
          }),
      term->atoms.end());
  return term->atoms.size() != before;
}

inline bool exact_order_region_remove_expr_density_time(
    ExactRegionCell *term,
    const semantic::Index expr_id,
    const semantic::Index time_id) {
  const auto before = term->atoms.size();
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return atom.kind == ExactRegionAtomKind::ExprDensity &&
                   atom.lhs.kind == ExactRegionVarKind::ExprTime &&
                   atom.lhs.id == expr_id &&
                   atom.rhs.kind == ExactRegionVarKind::Time &&
                   atom.rhs.id == time_id;
          }),
      term->atoms.end());
  return term->atoms.size() != before;
}

inline void exact_order_region_remove_orders_touching_time(
    ExactRegionCell *term,
    const semantic::Index time_id) {
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return atom.kind == ExactRegionAtomKind::TimeOrder &&
                   ((atom.lhs.kind == ExactRegionVarKind::Time &&
                     atom.lhs.id == time_id) ||
                    (atom.rhs.kind == ExactRegionVarKind::Time &&
                     atom.rhs.id == time_id));
          }),
      term->atoms.end());
  term->equalities.erase(
      std::remove_if(
          term->equalities.begin(),
          term->equalities.end(),
          [&](const ExactRegionEquality &equality) {
            return equality.lhs_time_id == time_id ||
                   equality.rhs_time_id == time_id;
          }),
      term->equalities.end());
}

inline bool exact_order_region_same_constraints(
    const ExactRegionCell &lhs,
    const ExactRegionCell &rhs);

enum class ExactOrderRegionAxisKind : std::uint8_t {
  Source = 0,
  Expr = 1
};

struct ExactOrderRegionAxis {
  ExactOrderRegionAxisKind kind{ExactOrderRegionAxisKind::Source};
  semantic::Index id{semantic::kInvalidIndex};
};

inline bool exact_order_region_axis_equal(
    const ExactOrderRegionAxis &lhs,
    const ExactOrderRegionAxis &rhs) noexcept {
  return lhs.kind == rhs.kind && lhs.id == rhs.id;
}

inline bool exact_order_region_axis_less(
    const ExactOrderRegionAxis &lhs,
    const ExactOrderRegionAxis &rhs) noexcept {
  if (lhs.kind != rhs.kind) {
    return lhs.kind < rhs.kind;
  }
  return lhs.id < rhs.id;
}

enum class ExactOrderRegionAxisSegmentKind : std::uint8_t {
  None = 0,
  Head = 1,
  Tail = 2,
  Interval = 3,
  All = 4
};

struct ExactOrderRegionAxisSegment {
  ExactOrderRegionAxisSegmentKind kind{
      ExactOrderRegionAxisSegmentKind::None};
  ExactOrderRegionAxis axis{};
  semantic::Index exact_time_id{semantic::kInvalidIndex};
  semantic::Index lower_time_id{semantic::kInvalidIndex};
  semantic::Index upper_time_id{semantic::kInvalidIndex};
  bool lower_inclusive{false};
  bool upper_inclusive{false};
};

inline bool exact_order_region_single_adjacent_interval_for_time(
    const ExactRegionCell &term,
    const semantic::Index exact_time_id,
    semantic::Index *lower_time_id,
    semantic::Index *upper_time_id) {
  if (exact_time_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto closure = exact_order_region_build_time_closure(term);
  if (closure.impossible) {
    return false;
  }
  *lower_time_id = semantic::kInvalidIndex;
  *upper_time_id = semantic::kInvalidIndex;
  const auto canonical_time_ids =
      exact_order_region_canonical_time_ids(closure);
  for (const auto candidate_time_id : canonical_time_ids) {
    if (candidate_time_id == exact_time_id) {
      continue;
    }
    const auto before_exact =
        exact_order_region_canonical_relation(
            closure, candidate_time_id, exact_time_id);
    if (before_exact != 0U) {
      bool dominated = false;
      for (const auto other_time_id : canonical_time_ids) {
        if (other_time_id == candidate_time_id ||
            other_time_id == exact_time_id) {
          continue;
        }
        const auto other_before_exact =
            exact_order_region_canonical_relation(
                closure, other_time_id, exact_time_id);
        const auto candidate_before_other =
            exact_order_region_canonical_relation(
                closure, candidate_time_id, other_time_id);
        if (other_before_exact != 0U && candidate_before_other != 0U) {
          dominated = true;
          break;
        }
      }
      if (!dominated) {
        if (*lower_time_id != semantic::kInvalidIndex) {
          return false;
        }
        *lower_time_id = candidate_time_id;
      }
    }
    const auto after_exact =
        exact_order_region_canonical_relation(
            closure, exact_time_id, candidate_time_id);
    if (after_exact != 0U) {
      bool dominated = false;
      for (const auto other_time_id : canonical_time_ids) {
        if (other_time_id == candidate_time_id ||
            other_time_id == exact_time_id) {
          continue;
        }
        const auto exact_before_other =
            exact_order_region_canonical_relation(
                closure, exact_time_id, other_time_id);
        const auto other_before_candidate =
            exact_order_region_canonical_relation(
                closure, other_time_id, candidate_time_id);
        if (exact_before_other != 0U && other_before_candidate != 0U) {
          dominated = true;
          break;
        }
      }
      if (!dominated) {
        if (*upper_time_id != semantic::kInvalidIndex) {
          return false;
        }
        *upper_time_id = candidate_time_id;
      }
    }
  }
  return *lower_time_id != semantic::kInvalidIndex &&
         *upper_time_id != semantic::kInvalidIndex;
}

inline bool exact_order_region_simple_source_axis_segment(
    const ExactRegionCell &term,
    const semantic::Index source_id,
    ExactOrderRegionAxisSegment *out) {
  *out = ExactOrderRegionAxisSegment{};
  out->axis = ExactOrderRegionAxis{ExactOrderRegionAxisKind::Source, source_id};

  semantic::Index exact_time_id{semantic::kInvalidIndex};
  for (const auto &exact : exact_region_exact_source_atoms(term)) {
    if (exact.source_id != source_id) {
      continue;
    }
    if (exact_time_id != semantic::kInvalidIndex) {
      return false;
    }
    exact_time_id = exact.time_id;
  }

  semantic::Index lower_bound{semantic::kInvalidIndex};
  semantic::Index upper_bound{semantic::kInvalidIndex};
  for (const auto &item : exact_region_lower_source_atoms(term)) {
    if (item.source_id != source_id) {
      continue;
    }
    if (lower_bound != semantic::kInvalidIndex) {
      return false;
    }
    lower_bound = item.time_id;
    out->lower_inclusive = item.inclusive;
  }
  for (const auto &item : exact_region_upper_source_atoms(term)) {
    if (item.source_id != source_id) {
      continue;
    }
    if (upper_bound != semantic::kInvalidIndex) {
      return false;
    }
    upper_bound = item.time_id;
    out->upper_inclusive = item.inclusive;
  }

  if (exact_time_id != semantic::kInvalidIndex &&
      lower_bound == semantic::kInvalidIndex &&
      upper_bound == semantic::kInvalidIndex &&
      exact_order_region_single_adjacent_interval_for_time(
          term, exact_time_id, &out->lower_time_id, &out->upper_time_id)) {
    out->kind = ExactOrderRegionAxisSegmentKind::Interval;
    out->exact_time_id = exact_time_id;
    return true;
  }
  if (exact_time_id == semantic::kInvalidIndex &&
      lower_bound != semantic::kInvalidIndex &&
      upper_bound == semantic::kInvalidIndex) {
    out->kind = ExactOrderRegionAxisSegmentKind::Tail;
    out->lower_time_id = lower_bound;
    return true;
  }
  if (exact_time_id == semantic::kInvalidIndex &&
      upper_bound != semantic::kInvalidIndex &&
      lower_bound == semantic::kInvalidIndex) {
    out->kind = ExactOrderRegionAxisSegmentKind::Head;
    out->upper_time_id = upper_bound;
    return true;
  }
  return false;
}

inline bool exact_order_region_simple_expr_axis_segment(
    const ExactRegionCell &term,
    const semantic::Index expr_id,
    ExactOrderRegionAxisSegment *out) {
  *out = ExactOrderRegionAxisSegment{};
  out->axis = ExactOrderRegionAxis{ExactOrderRegionAxisKind::Expr, expr_id};

  semantic::Index density_time_id{semantic::kInvalidIndex};
  ExactOrderRegionExprValueFactor value_factor;
  bool has_value_factor = false;
  for (const auto &factor : exact_region_expr_atoms(term)) {
    if (factor.expr_id != expr_id) {
      continue;
    }
    if (factor.density) {
      if (density_time_id != semantic::kInvalidIndex) {
        return false;
      }
      density_time_id = factor.time_id;
      continue;
    }
    if (has_value_factor) {
      return false;
    }
    value_factor = factor;
    has_value_factor = true;
  }

  if (density_time_id != semantic::kInvalidIndex && !has_value_factor &&
      exact_order_region_single_adjacent_interval_for_time(
          term, density_time_id, &out->lower_time_id, &out->upper_time_id)) {
    out->kind = ExactOrderRegionAxisSegmentKind::Interval;
    out->exact_time_id = density_time_id;
    return true;
  }
  if (density_time_id == semantic::kInvalidIndex && has_value_factor) {
    if (value_factor.before) {
      out->kind = ExactOrderRegionAxisSegmentKind::Head;
      out->upper_time_id = value_factor.time_id;
      out->upper_inclusive = value_factor.inclusive;
    } else {
      out->kind = ExactOrderRegionAxisSegmentKind::Tail;
      out->lower_time_id = value_factor.time_id;
      out->lower_inclusive = value_factor.inclusive;
    }
    return true;
  }
  return false;
}

inline bool exact_order_region_simple_axis_segment(
    const ExactRegionCell &term,
    const ExactOrderRegionAxis &axis,
    ExactOrderRegionAxisSegment *out) {
  if (axis.kind == ExactOrderRegionAxisKind::Source) {
    return exact_order_region_simple_source_axis_segment(term, axis.id, out);
  }
  return exact_order_region_simple_expr_axis_segment(term, axis.id, out);
}

inline ExactRegionCell exact_order_region_axis_segment_rest(
    ExactRegionCell term,
    const ExactOrderRegionAxisSegment &segment) {
  switch (segment.kind) {
  case ExactOrderRegionAxisSegmentKind::Head:
    if (segment.axis.kind == ExactOrderRegionAxisKind::Source) {
      exact_order_region_remove_upper_source_time(
          &term, segment.axis.id, segment.upper_time_id);
    } else {
      ExactOrderRegionExprValueFactor factor{
          segment.axis.id,
          segment.upper_time_id,
          true,
          segment.upper_inclusive,
          false};
      exact_order_region_remove_expr_value_atom(&term, factor);
    }
    break;
  case ExactOrderRegionAxisSegmentKind::Tail:
    if (segment.axis.kind == ExactOrderRegionAxisKind::Source) {
      exact_order_region_remove_lower_source_time(
          &term, segment.axis.id, segment.lower_time_id);
    } else {
      ExactOrderRegionExprValueFactor factor{
          segment.axis.id,
          segment.lower_time_id,
          false,
          segment.lower_inclusive,
          false};
      exact_order_region_remove_expr_value_atom(&term, factor);
    }
    break;
  case ExactOrderRegionAxisSegmentKind::Interval:
    if (segment.axis.kind == ExactOrderRegionAxisKind::Source) {
      exact_order_region_remove_exact_source_time(
          &term, segment.axis.id, segment.exact_time_id);
    } else {
      exact_order_region_remove_expr_density_time(
          &term, segment.axis.id, segment.exact_time_id);
    }
    exact_order_region_remove_orders_touching_time(
        &term, segment.exact_time_id);
    break;
  case ExactOrderRegionAxisSegmentKind::All:
  case ExactOrderRegionAxisSegmentKind::None:
    break;
  }
  exact_order_region_canonicalize_term(&term);
  return term;
}

inline bool exact_order_region_merge_adjacent_axis_segments(
    const ExactOrderRegionAxisSegment &lhs,
    const ExactOrderRegionAxisSegment &rhs,
    ExactOrderRegionAxisSegment *out) {
  if (!exact_order_region_axis_equal(lhs.axis, rhs.axis)) {
    return false;
  }
  *out = ExactOrderRegionAxisSegment{};
  out->axis = lhs.axis;
  if (lhs.kind == ExactOrderRegionAxisSegmentKind::All ||
      rhs.kind == ExactOrderRegionAxisSegmentKind::All) {
    out->kind = ExactOrderRegionAxisSegmentKind::All;
    return true;
  }
  if (lhs.kind == ExactOrderRegionAxisSegmentKind::Head &&
      rhs.kind == ExactOrderRegionAxisSegmentKind::Tail &&
      lhs.upper_time_id == rhs.lower_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::All;
    return true;
  }
  if (rhs.kind == ExactOrderRegionAxisSegmentKind::Head &&
      lhs.kind == ExactOrderRegionAxisSegmentKind::Tail &&
      rhs.upper_time_id == lhs.lower_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::All;
    return true;
  }
  if (lhs.kind == ExactOrderRegionAxisSegmentKind::Tail &&
      rhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      lhs.lower_time_id == rhs.upper_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::Tail;
    out->lower_time_id = rhs.lower_time_id;
    out->lower_inclusive = true;
    return true;
  }
  if (rhs.kind == ExactOrderRegionAxisSegmentKind::Tail &&
      lhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      rhs.lower_time_id == lhs.upper_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::Tail;
    out->lower_time_id = lhs.lower_time_id;
    out->lower_inclusive = true;
    return true;
  }
  if (lhs.kind == ExactOrderRegionAxisSegmentKind::Head &&
      rhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      lhs.upper_time_id == rhs.lower_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::Head;
    out->upper_time_id = rhs.upper_time_id;
    out->upper_inclusive = true;
    return true;
  }
  if (rhs.kind == ExactOrderRegionAxisSegmentKind::Head &&
      lhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      rhs.upper_time_id == lhs.lower_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::Head;
    out->upper_time_id = lhs.upper_time_id;
    out->upper_inclusive = true;
    return true;
  }
  if (lhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      rhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      lhs.upper_time_id == rhs.lower_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::Interval;
    out->exact_time_id = lhs.exact_time_id;
    out->lower_time_id = lhs.lower_time_id;
    out->upper_time_id = rhs.upper_time_id;
    return true;
  }
  if (lhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      rhs.kind == ExactOrderRegionAxisSegmentKind::Interval &&
      rhs.upper_time_id == lhs.lower_time_id) {
    out->kind = ExactOrderRegionAxisSegmentKind::Interval;
    out->exact_time_id = lhs.exact_time_id;
    out->lower_time_id = rhs.lower_time_id;
    out->upper_time_id = lhs.upper_time_id;
    return true;
  }
  return false;
}

inline ExactRegionCell exact_order_region_with_axis_segment(
    ExactRegionCell rest,
    const ExactOrderRegionAxisSegment &segment) {
  switch (segment.kind) {
  case ExactOrderRegionAxisSegmentKind::Head:
    if (segment.axis.kind == ExactOrderRegionAxisKind::Source) {
      exact_order_region_append_upper(
          &rest,
          segment.axis.id,
          segment.upper_time_id,
          segment.upper_inclusive);
    } else {
      exact_order_region_append_expr_factor(
          &rest,
          segment.axis.id,
          segment.upper_time_id,
          true,
          false,
          segment.upper_inclusive);
    }
    break;
  case ExactOrderRegionAxisSegmentKind::Tail:
    if (segment.axis.kind == ExactOrderRegionAxisKind::Source) {
      exact_order_region_append_lower(
          &rest,
          segment.axis.id,
          segment.lower_time_id,
          segment.lower_inclusive);
    } else {
      exact_order_region_append_expr_factor(
          &rest,
          segment.axis.id,
          segment.lower_time_id,
          false,
          false,
          segment.lower_inclusive);
    }
    break;
  case ExactOrderRegionAxisSegmentKind::Interval:
    if (segment.axis.kind == ExactOrderRegionAxisKind::Source) {
      exact_order_region_append_exact(
          &rest, segment.axis.id, segment.exact_time_id);
    } else {
      exact_order_region_append_expr_factor(
          &rest,
          segment.axis.id,
          segment.exact_time_id,
          true,
          true);
    }
    exact_order_region_append_time_order(
        &rest,
        segment.lower_time_id,
        segment.exact_time_id,
        !segment.lower_inclusive);
    exact_order_region_append_time_order(
        &rest,
        segment.exact_time_id,
        segment.upper_time_id,
        !segment.upper_inclusive);
    break;
  case ExactOrderRegionAxisSegmentKind::All:
  case ExactOrderRegionAxisSegmentKind::None:
    break;
  }
  exact_order_region_canonicalize_term(&rest);
  return rest;
}

inline ExactRegionCell exact_order_region_merged_axis_segment_rest(
    const ExactOrderRegionAxisSegment &lhs_segment,
    const ExactRegionCell &lhs_rest,
    const ExactOrderRegionAxisSegment &rhs_segment,
    const ExactRegionCell &rhs_rest,
    const ExactOrderRegionAxisSegment &merged) {
  if (merged.kind == ExactOrderRegionAxisSegmentKind::Tail) {
    if (lhs_segment.kind == ExactOrderRegionAxisSegmentKind::Tail) {
      return lhs_rest;
    }
    if (rhs_segment.kind == ExactOrderRegionAxisSegmentKind::Tail) {
      return rhs_rest;
    }
  }
  if (merged.kind == ExactOrderRegionAxisSegmentKind::Head) {
    if (lhs_segment.kind == ExactOrderRegionAxisSegmentKind::Head) {
      return lhs_rest;
    }
    if (rhs_segment.kind == ExactOrderRegionAxisSegmentKind::Head) {
      return rhs_rest;
    }
  }
  return lhs_rest;
}

inline std::vector<ExactOrderRegionAxis> exact_order_region_pair_axes(
    const ExactRegionCell &lhs,
    const ExactRegionCell &rhs) {
  std::vector<ExactOrderRegionAxis> out;
  const auto append = [&](const ExactOrderRegionAxis axis) {
    if (axis.id == semantic::kInvalidIndex) {
      return;
    }
    if (std::find_if(
            out.begin(),
            out.end(),
            [&](const ExactOrderRegionAxis &existing) {
              return exact_order_region_axis_equal(existing, axis);
            }) == out.end()) {
      out.push_back(axis);
    }
  };
  const auto append_source = [&](const semantic::Index source_id) {
    append(ExactOrderRegionAxis{ExactOrderRegionAxisKind::Source, source_id});
  };
  const auto append_expr = [&](const semantic::Index expr_id) {
    append(ExactOrderRegionAxis{ExactOrderRegionAxisKind::Expr, expr_id});
  };
  for (const auto &item : exact_region_exact_source_atoms(lhs)) {
    append_source(item.source_id);
  }
  for (const auto &item : exact_region_lower_source_atoms(lhs)) {
    append_source(item.source_id);
  }
  for (const auto &item : exact_region_upper_source_atoms(lhs)) {
    append_source(item.source_id);
  }
  for (const auto &factor : exact_region_expr_atoms(lhs)) {
    append_expr(factor.expr_id);
  }
  for (const auto &item : exact_region_exact_source_atoms(rhs)) {
    append_source(item.source_id);
  }
  for (const auto &item : exact_region_lower_source_atoms(rhs)) {
    append_source(item.source_id);
  }
  for (const auto &item : exact_region_upper_source_atoms(rhs)) {
    append_source(item.source_id);
  }
  for (const auto &factor : exact_region_expr_atoms(rhs)) {
    append_expr(factor.expr_id);
  }
  std::sort(out.begin(), out.end(), exact_order_region_axis_less);
  return out;
}

inline bool exact_order_region_try_merge_axis_partition_step(
    ExactRegionCell *lhs,
    ExactRegionCell *rhs,
    const bool allow_density_segments) {
  if (lhs->sign != rhs->sign || lhs->impossible || rhs->impossible) {
    return false;
  }
  for (const auto axis : exact_order_region_pair_axes(*lhs, *rhs)) {
    ExactOrderRegionAxisSegment lhs_segment;
    ExactOrderRegionAxisSegment rhs_segment;
    if (!exact_order_region_simple_axis_segment(*lhs, axis, &lhs_segment) ||
        !exact_order_region_simple_axis_segment(*rhs, axis, &rhs_segment) ||
        lhs_segment.kind == ExactOrderRegionAxisSegmentKind::None ||
        rhs_segment.kind == ExactOrderRegionAxisSegmentKind::None) {
      continue;
    }
    // Interval segments carry density binders; merging them is projection
    // algebra, not pure region canonicalization.
    if (!allow_density_segments &&
        (lhs_segment.kind == ExactOrderRegionAxisSegmentKind::Interval ||
         rhs_segment.kind == ExactOrderRegionAxisSegmentKind::Interval)) {
      continue;
    }
    auto lhs_rest =
        exact_order_region_axis_segment_rest(*lhs, lhs_segment);
    auto rhs_rest =
        exact_order_region_axis_segment_rest(*rhs, rhs_segment);
    if (!exact_order_region_same_constraints(lhs_rest, rhs_rest)) {
      continue;
    }
    ExactOrderRegionAxisSegment merged;
    if (!exact_order_region_merge_adjacent_axis_segments(
            lhs_segment, rhs_segment, &merged)) {
      continue;
    }
    auto rest = exact_order_region_merged_axis_segment_rest(
        lhs_segment, lhs_rest, rhs_segment, rhs_rest, merged);
    auto merged_term =
        exact_order_region_with_axis_segment(std::move(rest), merged);
    merged_term.sign = lhs->sign;
    *lhs = std::move(merged_term);
    rhs->sign = 0.0;
    return true;
  }
  return false;
}

inline std::vector<ExactRegionEquality> exact_region_constraint_equalities(
    ExactRegionCell term) {
  exact_region_normalize_equalities(&term);
  std::vector<ExactRegionEquality> out;
  out.reserve(term.equalities.size());
  for (const auto &equality : term.equalities) {
    if (equality.lhs_time_id == equality.rhs_time_id) {
      continue;
    }
    if (equality.origin ==
        ExactRegionEqualityOrigin::SharedLatentIdentity) {
      continue;
    }
    out.push_back(equality);
  }
  std::sort(out.begin(), out.end(), exact_region_equality_less);
  return out;
}

inline bool exact_order_region_same_constraints(
    const ExactRegionCell &lhs,
    const ExactRegionCell &rhs) {
  auto lhs_term = lhs;
  auto rhs_term = rhs;
  exact_order_region_canonicalize_term(&lhs_term);
  exact_order_region_canonicalize_term(&rhs_term);
  const auto lhs_closure = exact_order_region_build_time_closure(lhs_term);
  const auto rhs_closure = exact_order_region_build_time_closure(rhs_term);
  const auto lhs_orders =
      lhs_closure.impossible
          ? std::vector<ExactOrderRegionTimeOrder>{}
          : exact_order_region_canonical_time_orders(lhs_closure);
  const auto rhs_orders =
      rhs_closure.impossible
          ? std::vector<ExactOrderRegionTimeOrder>{}
          : exact_order_region_canonical_time_orders(rhs_closure);
  auto lhs_atoms = lhs_term.atoms;
  auto rhs_atoms = rhs_term.atoms;
  const auto lhs_equalities =
      exact_region_constraint_equalities(lhs_term);
  const auto rhs_equalities =
      exact_region_constraint_equalities(rhs_term);
  lhs_atoms.erase(
      std::remove_if(
          lhs_atoms.begin(),
          lhs_atoms.end(),
          [](const ExactRegionAtom &atom) {
            return atom.kind == ExactRegionAtomKind::TimeOrder;
          }),
      lhs_atoms.end());
  rhs_atoms.erase(
      std::remove_if(
          rhs_atoms.begin(),
          rhs_atoms.end(),
          [](const ExactRegionAtom &atom) {
            return atom.kind == ExactRegionAtomKind::TimeOrder;
          }),
      rhs_atoms.end());
  return lhs_term.impossible == rhs_term.impossible &&
         lhs_atoms.size() == rhs_atoms.size() &&
         lhs_orders.size() == rhs_orders.size() &&
         lhs_equalities.size() == rhs_equalities.size() &&
         std::equal(
             lhs_atoms.begin(),
             lhs_atoms.end(),
             rhs_atoms.begin(),
             exact_region_atom_equal) &&
         std::equal(
             lhs_equalities.begin(),
             lhs_equalities.end(),
             rhs_equalities.begin(),
             [](const auto &a, const auto &b) {
               return a.lhs_time_id == b.lhs_time_id &&
                      a.rhs_time_id == b.rhs_time_id &&
                      a.mass == b.mass &&
                      a.origin == b.origin;
             }) &&
         std::equal(
             lhs_orders.begin(), lhs_orders.end(),
             rhs_orders.begin(),
             [](const auto &a, const auto &b) {
               return a.before_time_id == b.before_time_id &&
                      a.after_time_id == b.after_time_id &&
                      a.strict == b.strict;
             });
}

inline ExactOrderRegionExpr exact_order_region_simplify(
    ExactOrderRegionExpr expr);

inline bool exact_region_atom_boolean(
    const ExactRegionAtom &atom) noexcept {
  return atom.kind == ExactRegionAtomKind::SourceLower ||
         atom.kind == ExactRegionAtomKind::SourceUpper ||
         atom.kind == ExactRegionAtomKind::ExprBefore ||
         atom.kind == ExactRegionAtomKind::ExprNotBefore ||
         atom.kind == ExactRegionAtomKind::TimeOrder ||
         atom.kind == ExactRegionAtomKind::OutcomeUnused ||
         atom.kind == ExactRegionAtomKind::OutcomeUsed;
}

inline std::vector<ExactRegionAtom> exact_order_region_boolean_atoms(
    const ExactRegionCell &term) {
  std::vector<ExactRegionAtom> atoms;
  atoms.reserve(term.atoms.size());
  for (const auto &atom : term.atoms) {
    if (exact_region_atom_boolean(atom)) {
      atoms.push_back(atom);
    }
  }
  std::sort(atoms.begin(), atoms.end(), exact_region_atom_less);
  atoms.erase(
      std::unique(atoms.begin(), atoms.end(), exact_region_atom_equal),
      atoms.end());
  return atoms;
}

inline std::vector<ExactRegionAtom> exact_order_region_density_atoms(
    const ExactRegionCell &term) {
  std::vector<ExactRegionAtom> atoms;
  atoms.reserve(term.atoms.size());
  for (const auto &atom : term.atoms) {
    if (atom.kind == ExactRegionAtomKind::SourceExact ||
        atom.kind == ExactRegionAtomKind::ExprDensity) {
      atoms.push_back(atom);
    }
  }
  std::sort(atoms.begin(), atoms.end(), exact_region_atom_less);
  atoms.erase(
      std::unique(atoms.begin(), atoms.end(), exact_region_atom_equal),
      atoms.end());
  return atoms;
}

inline void exact_order_region_append_atom(
    ExactRegionCell *term,
    const ExactRegionAtom &atom) {
  switch (atom.kind) {
  case ExactRegionAtomKind::SourceExact:
    exact_order_region_append_exact(term, atom.lhs.id, atom.rhs.id);
    break;
  case ExactRegionAtomKind::SourceLower:
    exact_order_region_append_lower(
        term, atom.lhs.id, atom.rhs.id, atom.inclusive);
    break;
  case ExactRegionAtomKind::SourceUpper:
    exact_order_region_append_upper(
        term, atom.lhs.id, atom.rhs.id, atom.inclusive);
    break;
	  case ExactRegionAtomKind::ExprBefore:
	    exact_order_region_append_expr_factor(
	        term, atom.lhs.id, atom.rhs.id, true, false, atom.inclusive);
	    break;
	  case ExactRegionAtomKind::ExprNotBefore:
	    exact_order_region_append_expr_factor(
	        term, atom.lhs.id, atom.rhs.id, false, false, atom.inclusive);
	    break;
  case ExactRegionAtomKind::ExprDensity:
    exact_order_region_append_expr_factor(
        term, atom.lhs.id, atom.rhs.id, true, true);
    break;
  case ExactRegionAtomKind::TimeOrder:
    exact_order_region_append_time_order(
        term, atom.lhs.id, atom.rhs.id, atom.strict);
    break;
  case ExactRegionAtomKind::OutcomeUnused:
    exact_order_region_append_outcome_gate(term, atom.outcome_indices);
    break;
  case ExactRegionAtomKind::OutcomeUsed:
    exact_order_region_append_outcome_used_gate(term, atom.outcome_indices);
    break;
  }
}

inline void exact_order_region_append_atom_complement(
    ExactRegionCell *term,
    const ExactRegionAtom &atom) {
  switch (atom.kind) {
  case ExactRegionAtomKind::SourceExact:
    break;
  case ExactRegionAtomKind::SourceLower:
    exact_order_region_append_upper(
        term, atom.lhs.id, atom.rhs.id, !atom.inclusive);
    break;
  case ExactRegionAtomKind::SourceUpper:
    exact_order_region_append_lower(
        term, atom.lhs.id, atom.rhs.id, !atom.inclusive);
    break;
	  case ExactRegionAtomKind::ExprBefore:
	    exact_order_region_append_expr_factor(
	        term, atom.lhs.id, atom.rhs.id, false, false, !atom.inclusive);
	    break;
	  case ExactRegionAtomKind::ExprNotBefore:
	    exact_order_region_append_expr_factor(
	        term, atom.lhs.id, atom.rhs.id, true, false, !atom.inclusive);
	    break;
  case ExactRegionAtomKind::ExprDensity:
    break;
  case ExactRegionAtomKind::TimeOrder:
    exact_order_region_append_time_order(
        term, atom.rhs.id, atom.lhs.id, !atom.strict);
    break;
  case ExactRegionAtomKind::OutcomeUnused:
    exact_order_region_append_outcome_used_gate(term, atom.outcome_indices);
    break;
  case ExactRegionAtomKind::OutcomeUsed:
    exact_order_region_append_outcome_gate(term, atom.outcome_indices);
    break;
  }
}

inline bool exact_order_region_remove_atom(
    ExactRegionCell *term,
    const ExactRegionAtom &needle) {
  const auto before = term->atoms.size();
  term->atoms.erase(
      std::remove_if(
          term->atoms.begin(),
          term->atoms.end(),
          [&](const ExactRegionAtom &atom) {
            return exact_region_atom_equal(atom, needle);
          }),
      term->atoms.end());
  return term->atoms.size() != before;
}

inline bool exact_order_region_try_merge_boolean_partition_side(
    ExactRegionCell *lhs,
    ExactRegionCell *rhs) {
  for (const auto &atom : exact_order_region_boolean_atoms(*lhs)) {
    auto rest = *lhs;
    if (!exact_order_region_remove_atom(&rest, atom)) {
      continue;
    }
    exact_order_region_canonicalize_term(&rest);
    if (rest.impossible) {
      continue;
    }

    auto complement = rest;
    exact_order_region_append_atom_complement(&complement, atom);
    exact_order_region_canonicalize_term(&complement);
    if (complement.impossible ||
        !exact_order_region_same_constraints(complement, *rhs)) {
      continue;
    }

    rest.sign = lhs->sign;
    *lhs = std::move(rest);
    rhs->sign = 0.0;
    return true;
  }
  return false;
}

inline bool exact_order_region_try_merge_boolean_partition_step(
    ExactRegionCell *lhs,
    ExactRegionCell *rhs) {
  if (lhs->sign != rhs->sign || lhs->impossible || rhs->impossible) {
    return false;
  }
  return exact_order_region_try_merge_boolean_partition_side(lhs, rhs) ||
         exact_order_region_try_merge_boolean_partition_side(rhs, lhs);
}

inline bool exact_order_region_try_minimize_union_step(
    std::vector<ExactRegionCell> *terms,
    const bool allow_density_segments) {
  for (std::size_t i = 0; i < terms->size(); ++i) {
    if ((*terms)[i].sign == 0.0) {
      continue;
    }
    for (std::size_t j = 0; j < terms->size(); ++j) {
      if (i == j || (*terms)[j].sign == 0.0) {
        continue;
      }
      if (exact_order_region_try_merge_boolean_partition_step(
              &(*terms)[i], &(*terms)[j]) ||
          exact_order_region_try_merge_axis_partition_step(
              &(*terms)[i], &(*terms)[j], allow_density_segments)) {
        return true;
      }
    }
  }
  return false;
}

inline void exact_order_region_minimize_union_cells(
    std::vector<ExactRegionCell> *terms,
    const bool allow_density_segments) {
  while (exact_order_region_try_minimize_union_step(
      terms, allow_density_segments)) {
  }
}

inline bool exact_order_region_equality_implied(
    const ExactRegionCell &term,
    const ExactRegionEquality &equality) {
  auto with_equality = term;
  exact_region_append_equality(
      &with_equality,
      equality.lhs_time_id,
      equality.rhs_time_id,
      equality.mass,
      equality.origin);
  exact_order_region_canonicalize_term(&with_equality);
  if (with_equality.impossible) {
    return false;
  }
  return exact_order_region_same_constraints(term, with_equality);
}

inline ExactOrderRegionExpr exact_order_region_equality_complement(
    const ExactRegionCell &term,
    const ExactRegionEquality &equality) {
  ExactOrderRegionExpr out;
  if (equality.lhs_time_id == equality.rhs_time_id ||
      equality.origin ==
          ExactRegionEqualityOrigin::SharedLatentIdentity) {
    return out;
  }
  if (equality.mass == ExactRegionEqualityMass::MeasureZero) {
    out.terms.push_back(term);
    return out;
  }
  const auto append_live = [&out](ExactRegionCell cell) {
    exact_order_region_canonicalize_term(&cell);
    if (!cell.impossible && cell.sign != 0.0) {
      out.terms.push_back(std::move(cell));
    }
  };

  auto lhs_before_rhs = term;
  exact_order_region_append_time_order(
      &lhs_before_rhs, equality.lhs_time_id, equality.rhs_time_id);
  append_live(std::move(lhs_before_rhs));

  auto rhs_before_lhs = term;
  exact_order_region_append_time_order(
      &rhs_before_lhs, equality.rhs_time_id, equality.lhs_time_id);
  append_live(std::move(rhs_before_lhs));
  return out;
}

inline ExactRegionCell exact_order_region_canonicalized(
    ExactRegionCell term) {
  exact_order_region_canonicalize_term(&term);
  return term;
}

inline bool exact_order_region_term_subset_of(
    const ExactRegionCell &lhs,
    const ExactRegionCell &rhs) {
  auto lhs_term = exact_order_region_canonicalized(lhs);
  auto intersection =
      exact_order_region_canonicalized(
          exact_order_region_intersect_terms(lhs_term, rhs));
  return !lhs_term.impossible &&
         !intersection.impossible &&
         exact_order_region_same_constraints(lhs_term, intersection);
}

inline void exact_order_region_append_live_term(
    ExactOrderRegionExpr *expr,
    ExactRegionCell term) {
  exact_order_region_canonicalize_term(&term);
  if (!term.impossible && term.sign != 0.0) {
    expr->terms.push_back(std::move(term));
  }
}

inline bool exact_region_density_subject_present(
    const ExactRegionCell &term,
    const ExactRegionAtom &atom) {
  if (atom.kind == ExactRegionAtomKind::SourceExact) {
    return exact_order_region_exact_time_for_source(term, atom.lhs.id) !=
           semantic::kInvalidIndex;
  }
  if (atom.kind == ExactRegionAtomKind::ExprDensity) {
    return exact_region_expr_density_time(term, atom.lhs.id) !=
           semantic::kInvalidIndex;
  }
  return true;
}

inline bool exact_region_missing_density_atom_has_positive_overlap(
    const ExactRegionCell &term,
    const ExactRegionAtom &atom) {
  if (!exact_region_time_is_latent_variable(atom.rhs.id)) {
    return false;
  }
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  if (atom.kind == ExactRegionAtomKind::SourceExact) {
    for (const auto &bound : exact_region_lower_source_atoms(term)) {
      if (bound.source_id == atom.lhs.id &&
          exact_order_region_time_known_before_or_equal(
              term, observed_time_id, bound.time_id)) {
        return false;
      }
    }
    return true;
  }
  if (atom.kind == ExactRegionAtomKind::ExprDensity) {
    for (const auto &factor : exact_region_expr_atoms(term)) {
      if (factor.expr_id == atom.lhs.id &&
          !factor.density &&
          !factor.before &&
          exact_order_region_time_known_before_or_equal(
              term, observed_time_id, factor.time_id)) {
        return false;
      }
    }
    return true;
  }
  return false;
}

inline std::vector<ExactRegionAtom> exact_region_missing_density_atoms(
    const ExactRegionCell &minuend,
    const ExactRegionCell &subtrahend,
    bool *zero_measure_overlap) {
  *zero_measure_overlap = false;
  std::vector<ExactRegionAtom> out;
  for (const auto &atom : subtrahend.atoms) {
    if (atom.kind != ExactRegionAtomKind::SourceExact &&
        atom.kind != ExactRegionAtomKind::ExprDensity) {
      continue;
    }
    if (exact_region_density_subject_present(minuend, atom)) {
      continue;
    }
    if (!exact_region_missing_density_atom_has_positive_overlap(
            minuend, atom)) {
      *zero_measure_overlap = true;
      return {};
    }
    out.push_back(atom);
  }
  std::sort(out.begin(), out.end(), exact_region_atom_less);
  out.erase(
      std::unique(out.begin(), out.end(), exact_region_atom_equal),
      out.end());
  return out;
}

inline ExactOrderRegionExpr exact_order_region_lift_cell_for_density_atom(
    const ExactRegionCell &term,
    const ExactRegionAtom &atom) {
  ExactOrderRegionExpr out;
  if (atom.kind == ExactRegionAtomKind::SourceExact) {
    auto head_branch = term;
    exact_order_region_append_upper(
        &head_branch, atom.lhs.id, atom.rhs.id);
    exact_order_region_append_live_term(&out, std::move(head_branch));

    auto exact_branch = term;
    exact_order_region_append_exact(
        &exact_branch, atom.lhs.id, atom.rhs.id);
    exact_order_region_append_live_term(&out, std::move(exact_branch));

    auto tail_branch = term;
    exact_order_region_append_lower(
        &tail_branch, atom.lhs.id, atom.rhs.id);
    exact_order_region_append_live_term(&out, std::move(tail_branch));
    return out;
  }
  if (atom.kind == ExactRegionAtomKind::ExprDensity) {
    auto before_branch = term;
    exact_order_region_append_expr_factor(
        &before_branch, atom.lhs.id, atom.rhs.id, true);
    exact_order_region_append_live_term(&out, std::move(before_branch));

    auto exact_branch = term;
    exact_order_region_append_expr_factor(
        &exact_branch, atom.lhs.id, atom.rhs.id, true, true);
    exact_order_region_append_live_term(&out, std::move(exact_branch));

    auto after_branch = term;
    exact_order_region_append_expr_factor(
        &after_branch, atom.lhs.id, atom.rhs.id, false);
    exact_order_region_append_live_term(&out, std::move(after_branch));
    return out;
  }
  out.terms.push_back(term);
  return out;
}

inline ExactOrderRegionExpr exact_order_region_subtract_term(
    const ExactRegionCell &minuend,
    const ExactRegionCell &subtrahend) {
  ExactOrderRegionExpr out;
  auto minuend_term = exact_order_region_canonicalized(minuend);
  auto subtrahend_term = exact_order_region_canonicalized(subtrahend);
  if (minuend_term.impossible || minuend_term.sign == 0.0) {
    return out;
  }
  bool zero_measure_overlap = false;
  const auto missing_density_atoms =
      exact_region_missing_density_atoms(
          minuend_term, subtrahend_term, &zero_measure_overlap);
  if (zero_measure_overlap) {
    out.terms.push_back(std::move(minuend_term));
    return out;
  }
  if (!missing_density_atoms.empty()) {
    ExactOrderRegionExpr exact_overlap;
    exact_overlap.terms.push_back(std::move(minuend_term));
    for (const auto &atom : missing_density_atoms) {
      ExactOrderRegionExpr next_exact_overlap;
      for (const auto &term : exact_overlap.terms) {
        const auto lifted =
            exact_order_region_lift_cell_for_density_atom(term, atom);
        for (const auto &part : lifted.terms) {
          if (exact_region_density_subject_present(part, atom)) {
            next_exact_overlap.terms.push_back(part);
          } else {
            out.terms.push_back(part);
          }
        }
      }
      exact_overlap = std::move(next_exact_overlap);
      if (exact_overlap.terms.empty()) {
        break;
      }
    }
    for (const auto &term : exact_overlap.terms) {
      exact_order_region_append_expr(
          &out,
          exact_order_region_subtract_term(term, subtrahend_term));
    }
    return exact_order_region_simplify(std::move(out));
  }
  auto intersection =
      exact_order_region_canonicalized(
          exact_order_region_intersect_terms(minuend_term, subtrahend_term));
  if (intersection.impossible || intersection.sign == 0.0) {
    out.terms.push_back(std::move(minuend_term));
    return out;
  }
  if (exact_order_region_term_subset_of(minuend_term, subtrahend_term)) {
    return out;
  }

  auto prefix = minuend_term;
  for (const auto &atom : exact_order_region_boolean_atoms(subtrahend_term)) {
    auto residual = prefix;
    exact_order_region_append_atom_complement(&residual, atom);
    exact_order_region_append_live_term(&out, std::move(residual));

    exact_order_region_append_atom(&prefix, atom);
    exact_order_region_canonicalize_term(&prefix);
    if (prefix.impossible || prefix.sign == 0.0) {
      return out;
    }
  }

  for (const auto &atom : exact_order_region_density_atoms(subtrahend_term)) {
    if (!exact_region_density_subject_present(prefix, atom)) {
      const auto lifted =
          exact_order_region_lift_cell_for_density_atom(prefix, atom);
      for (const auto &term : lifted.terms) {
        exact_order_region_append_expr(
            &out,
            exact_order_region_subtract_term(term, subtrahend_term));
      }
      return exact_order_region_simplify(std::move(out));
    }
    exact_order_region_append_atom(&prefix, atom);
    exact_order_region_canonicalize_term(&prefix);
    if (prefix.impossible || prefix.sign == 0.0) {
      return out;
    }
  }

  for (const auto &equality :
       exact_region_constraint_equalities(subtrahend_term)) {
    if (exact_order_region_equality_implied(prefix, equality)) {
      continue;
    }
    exact_order_region_append_expr(
        &out,
        exact_order_region_equality_complement(prefix, equality));

    exact_region_append_equality(
        &prefix,
        equality.lhs_time_id,
        equality.rhs_time_id,
        equality.mass,
        equality.origin);
    exact_order_region_canonicalize_term(&prefix);
    if (prefix.impossible || prefix.sign == 0.0) {
      return out;
    }
  }

  return exact_order_region_simplify(std::move(out));
}

inline ExactOrderRegionExpr exact_order_region_subtract_expr(
    const ExactOrderRegionExpr &minuend,
    const ExactRegionCell &subtrahend) {
  ExactOrderRegionExpr out;
  for (const auto &term : minuend.terms) {
    exact_order_region_append_expr(
        &out,
        exact_order_region_subtract_term(term, subtrahend));
  }
  return out;
}

inline ExactOrderRegionExpr exact_order_region_union(
    const ExactOrderRegionExpr &lhs,
    const ExactOrderRegionExpr &rhs) {
  ExactOrderRegionExpr out = exact_order_region_simplify(lhs);
  for (const auto &rhs_term : exact_order_region_simplify(rhs).terms) {
    ExactOrderRegionExpr residual;
    residual.terms.push_back(rhs_term);
    for (const auto &existing : out.terms) {
      residual = exact_order_region_subtract_expr(residual, existing);
      if (residual.terms.empty()) {
        break;
      }
    }
    exact_order_region_append_expr(&out, std::move(residual));
  }
  return exact_order_region_simplify(std::move(out));
}

inline ExactOrderRegionExpr exact_order_region_minimize_positive_union(
    ExactOrderRegionExpr expr) {
  expr = exact_order_region_simplify(std::move(expr));
  exact_order_region_minimize_union_cells(&expr.terms, true);
  for (std::size_t i = 0; i < expr.terms.size(); ++i) {
    if (expr.terms[i].sign <= 0.0 || expr.terms[i].impossible) {
      continue;
    }
    for (std::size_t j = 0; j < expr.terms.size(); ++j) {
      if (i == j ||
          expr.terms[j].sign <= 0.0 ||
          expr.terms[j].impossible) {
        continue;
      }
      if (exact_order_region_term_subset_of(expr.terms[j], expr.terms[i])) {
        expr.terms[j].sign = 0.0;
      }
    }
  }
  return exact_order_region_simplify(std::move(expr));
}

inline ExactOrderRegionExpr exact_order_region_simplify(
    ExactOrderRegionExpr expr) {
  std::vector<ExactRegionCell> simplified;
  simplified.reserve(expr.terms.size());
  for (auto &term : expr.terms) {
    if (term.impossible || term.sign == 0.0) {
      continue;
    }
    exact_order_region_canonicalize_term(&term);
    bool merged = false;
    for (auto &existing : simplified) {
      if (exact_order_region_same_constraints(existing, term)) {
        existing.sign += term.sign;
        merged = true;
        break;
      }
    }
    if (!merged) {
      simplified.push_back(std::move(term));
    }
  }
  exact_order_region_minimize_union_cells(&simplified, false);
  for (std::size_t i = 0; i < simplified.size(); ++i) {
    if (simplified[i].sign == 0.0) {
      continue;
    }
    for (std::size_t j = i + 1; j < simplified.size(); ++j) {
      if (simplified[j].sign == 0.0) {
        continue;
      }
      if (exact_order_region_same_constraints(simplified[i], simplified[j])) {
        simplified[i].sign += simplified[j].sign;
        simplified[j].sign = 0.0;
      }
    }
  }
  simplified.erase(
      std::remove_if(
          simplified.begin(),
          simplified.end(),
          [](const ExactRegionCell &term) {
            return term.sign == 0.0;
          }),
      simplified.end());
  expr.terms = std::move(simplified);
  return expr;
}


} // namespace detail
} // namespace accumulatr::eval
