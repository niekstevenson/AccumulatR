#pragma once

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "exact_region_builders.hpp"

namespace accumulatr::eval {
namespace detail {

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

inline bool exact_order_region_source_bounds_consistent(
    const ExactRegionCell &term) {
  const auto lower_bounds = exact_region_lower_source_atoms(term);
  const auto upper_bounds = exact_region_upper_source_atoms(term);
  for (const auto &lower : lower_bounds) {
    for (const auto &upper : upper_bounds) {
      if (lower.source_id != upper.source_id) {
        continue;
      }
      if (exact_order_region_time_known_before_or_equal(
              term, upper.time_id, lower.time_id)) {
        return false;
      }
    }
  }
  return true;
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
      !exact_order_region_source_bounds_consistent(*term) ||
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

} // namespace detail
} // namespace accumulatr::eval
