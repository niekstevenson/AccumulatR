#pragma once

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <vector>

#include "exact_region_closure.hpp"

namespace accumulatr::eval {
namespace detail {

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
