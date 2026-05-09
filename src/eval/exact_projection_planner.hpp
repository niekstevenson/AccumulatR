#pragma once

#include "exact_types.hpp"
#include "exact_region.hpp"
#include "exact_compiled_math_lowering.hpp"

namespace accumulatr::eval {
namespace detail {

inline void exact_order_region_reduce_bounds(
    const ExactRegionCell &term,
    std::vector<semantic::Index> *times,
    const bool keep_latest) {
  std::vector<semantic::Index> reduced;
  for (const auto time_id : *times) {
    bool dominated = false;
    for (const auto other_time_id : *times) {
      if (time_id == other_time_id) {
        continue;
      }
      if (keep_latest) {
        if (exact_order_region_time_known_before_or_equal(
                term, time_id, other_time_id)) {
          dominated = true;
          break;
        }
      } else {
        if (exact_order_region_time_known_before_or_equal(
                term, other_time_id, time_id)) {
          dominated = true;
          break;
        }
      }
    }
    if (!dominated &&
        std::find(reduced.begin(), reduced.end(), time_id) == reduced.end()) {
      reduced.push_back(time_id);
    }
  }
  *times = std::move(reduced);
}

inline semantic::Index exact_order_region_source_interval_node(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const semantic::Index lower_time_id,
    const semantic::Index upper_time_id,
    const semantic::Index source_view_id) {
  const auto upper_cdf =
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourceCdf,
          source_id,
          0,
          upper_time_id,
          source_view_id);
  const auto lower_cdf =
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourceCdf,
          source_id,
          0,
          lower_time_id,
          source_view_id);
  const auto interval =
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::CleanSignedSum,
          std::vector<semantic::Index>{
              upper_cdf,
              compiled_math_unary_node(
                  &plan->compiled_math,
                  CompiledMathNodeKind::Negate,
                  lower_cdf)},
          CompiledMathValueKind::Scalar);
  return compiled_math_time_gate_node(
      &plan->compiled_math,
      interval,
      upper_time_id,
      lower_time_id,
      CompiledMathValueKind::Scalar);
}

inline semantic::Index exact_order_region_source_cdf_min_node(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id) {
  if (upper_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (upper_time_ids.size() == 1U) {
    return compiled_math_source_node(
        &plan->compiled_math,
        CompiledMathNodeKind::SourceCdf,
        source_id,
        0,
        upper_time_ids.front(),
        source_view_id);
  }
  std::vector<semantic::Index> candidates;
  candidates.reserve(upper_time_ids.size());
  for (std::size_t i = 0; i < upper_time_ids.size(); ++i) {
    auto node =
        compiled_math_source_node(
            &plan->compiled_math,
            CompiledMathNodeKind::SourceCdf,
            source_id,
            0,
            upper_time_ids[i],
            source_view_id);
    for (std::size_t j = 0; j < upper_time_ids.size(); ++j) {
      if (i == j) {
        continue;
      }
      node =
          (j < i ? compiled_math_strict_time_gate_node
                 : compiled_math_time_gate_node)(
              &plan->compiled_math,
              node,
              upper_time_ids[j],
              upper_time_ids[i],
              CompiledMathValueKind::Scalar);
    }
    candidates.push_back(node);
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Sum,
      std::move(candidates),
      CompiledMathValueKind::Scalar);
}

inline semantic::Index exact_order_region_source_survival_max_node(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const semantic::Index source_view_id) {
  if (lower_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (lower_time_ids.size() == 1U) {
    return compiled_math_source_node(
        &plan->compiled_math,
        CompiledMathNodeKind::SourceSurvival,
        source_id,
        0,
        lower_time_ids.front(),
        source_view_id);
  }
  std::vector<semantic::Index> candidates;
  candidates.reserve(lower_time_ids.size());
  for (std::size_t i = 0; i < lower_time_ids.size(); ++i) {
    auto node =
        compiled_math_source_node(
            &plan->compiled_math,
            CompiledMathNodeKind::SourceSurvival,
            source_id,
            0,
            lower_time_ids[i],
            source_view_id);
    for (std::size_t j = 0; j < lower_time_ids.size(); ++j) {
      if (i == j) {
        continue;
      }
      node =
          (j < i ? compiled_math_strict_time_gate_node
                 : compiled_math_time_gate_node)(
              &plan->compiled_math,
              node,
              lower_time_ids[i],
              lower_time_ids[j],
              CompiledMathValueKind::Scalar);
    }
    candidates.push_back(node);
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Sum,
      std::move(candidates),
      CompiledMathValueKind::Scalar);
}

inline semantic::Index exact_order_region_source_interval_partition_node(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id) {
  if (lower_time_ids.empty()) {
    return exact_order_region_source_cdf_min_node(
        plan, source_id, upper_time_ids, source_view_id);
  }
  if (upper_time_ids.empty()) {
    return exact_order_region_source_survival_max_node(
        plan, source_id, lower_time_ids, source_view_id);
  }

  std::vector<semantic::Index> candidates;
  candidates.reserve(lower_time_ids.size() * upper_time_ids.size());
  for (std::size_t lower_idx = 0; lower_idx < lower_time_ids.size();
       ++lower_idx) {
    const auto lower_time_id = lower_time_ids[lower_idx];
    for (std::size_t upper_idx = 0; upper_idx < upper_time_ids.size();
         ++upper_idx) {
      const auto upper_time_id = upper_time_ids[upper_idx];
      auto node =
          exact_order_region_source_interval_node(
              plan,
              source_id,
              lower_time_id,
              upper_time_id,
              source_view_id);
      for (std::size_t other = 0; other < lower_time_ids.size(); ++other) {
        if (other == lower_idx) {
          continue;
        }
        node =
            (other < lower_idx ? compiled_math_strict_time_gate_node
                               : compiled_math_time_gate_node)(
                &plan->compiled_math,
                node,
                lower_time_id,
                lower_time_ids[other],
                CompiledMathValueKind::Scalar);
      }
      for (std::size_t other = 0; other < upper_time_ids.size(); ++other) {
        if (other == upper_idx) {
          continue;
        }
        node =
            (other < upper_idx ? compiled_math_strict_time_gate_node
                               : compiled_math_time_gate_node)(
                &plan->compiled_math,
                node,
                upper_time_ids[other],
                upper_time_id,
                CompiledMathValueKind::Scalar);
      }
      candidates.push_back(node);
    }
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Sum,
      std::move(candidates),
      CompiledMathValueKind::Scalar);
}

inline semantic::Index exact_order_region_expr_value_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind kind,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  return compile_expr_value_node(
      plan,
      expr_id,
      kind,
      0,
      time_id,
      source_view_id);
}

inline semantic::Index exact_order_region_expr_interval_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const semantic::Index lower_time_id,
    const semantic::Index upper_time_id,
    const semantic::Index source_view_id) {
  const auto upper_cdf =
      exact_order_region_expr_value_node(
          plan,
          expr_id,
          CompiledMathNodeKind::ExprCdf,
          upper_time_id,
          source_view_id);
  const auto lower_cdf =
      exact_order_region_expr_value_node(
          plan,
          expr_id,
          CompiledMathNodeKind::ExprCdf,
          lower_time_id,
          source_view_id);
  const auto interval =
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::CleanSignedSum,
          std::vector<semantic::Index>{
              upper_cdf,
              compiled_math_unary_node(
                  &plan->compiled_math,
                  CompiledMathNodeKind::Negate,
                  lower_cdf)},
          CompiledMathValueKind::Scalar);
  return compiled_math_time_gate_node(
      &plan->compiled_math,
      interval,
      upper_time_id,
      lower_time_id,
      CompiledMathValueKind::Scalar);
}

inline semantic::Index exact_order_region_expr_cdf_min_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id) {
  if (upper_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (upper_time_ids.size() == 1U) {
    return exact_order_region_expr_value_node(
        plan,
        expr_id,
        CompiledMathNodeKind::ExprCdf,
        upper_time_ids.front(),
        source_view_id);
  }
  std::vector<semantic::Index> candidates;
  candidates.reserve(upper_time_ids.size());
  for (std::size_t i = 0; i < upper_time_ids.size(); ++i) {
    auto node =
        exact_order_region_expr_value_node(
            plan,
            expr_id,
            CompiledMathNodeKind::ExprCdf,
            upper_time_ids[i],
            source_view_id);
    for (std::size_t j = 0; j < upper_time_ids.size(); ++j) {
      if (i == j) {
        continue;
      }
      node =
          (j < i ? compiled_math_strict_time_gate_node
                 : compiled_math_time_gate_node)(
              &plan->compiled_math,
              node,
              upper_time_ids[j],
              upper_time_ids[i],
              CompiledMathValueKind::Scalar);
    }
    candidates.push_back(node);
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Sum,
      std::move(candidates),
      CompiledMathValueKind::Scalar);
}

inline semantic::Index exact_order_region_expr_survival_max_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const semantic::Index source_view_id) {
  if (lower_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (lower_time_ids.size() == 1U) {
    return exact_order_region_expr_value_node(
        plan,
        expr_id,
        CompiledMathNodeKind::ExprSurvival,
        lower_time_ids.front(),
        source_view_id);
  }
  std::vector<semantic::Index> candidates;
  candidates.reserve(lower_time_ids.size());
  for (std::size_t i = 0; i < lower_time_ids.size(); ++i) {
    auto node =
        exact_order_region_expr_value_node(
            plan,
            expr_id,
            CompiledMathNodeKind::ExprSurvival,
            lower_time_ids[i],
            source_view_id);
    for (std::size_t j = 0; j < lower_time_ids.size(); ++j) {
      if (i == j) {
        continue;
      }
      node =
          (j < i ? compiled_math_strict_time_gate_node
                 : compiled_math_time_gate_node)(
              &plan->compiled_math,
              node,
              lower_time_ids[i],
              lower_time_ids[j],
              CompiledMathValueKind::Scalar);
    }
    candidates.push_back(node);
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Sum,
      std::move(candidates),
      CompiledMathValueKind::Scalar);
}

inline semantic::Index exact_order_region_expr_interval_partition_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id) {
  if (lower_time_ids.empty()) {
    return exact_order_region_expr_cdf_min_node(
        plan, expr_id, upper_time_ids, source_view_id);
  }
  if (upper_time_ids.empty()) {
    return exact_order_region_expr_survival_max_node(
        plan, expr_id, lower_time_ids, source_view_id);
  }

  std::vector<semantic::Index> candidates;
  candidates.reserve(lower_time_ids.size() * upper_time_ids.size());
  for (std::size_t lower_idx = 0; lower_idx < lower_time_ids.size();
       ++lower_idx) {
    const auto lower_time_id = lower_time_ids[lower_idx];
    for (std::size_t upper_idx = 0; upper_idx < upper_time_ids.size();
         ++upper_idx) {
      const auto upper_time_id = upper_time_ids[upper_idx];
      auto node =
          exact_order_region_expr_interval_node(
              plan,
              expr_id,
              lower_time_id,
              upper_time_id,
              source_view_id);
      for (std::size_t other = 0; other < lower_time_ids.size(); ++other) {
        if (other == lower_idx) {
          continue;
        }
        node =
            (other < lower_idx ? compiled_math_strict_time_gate_node
                               : compiled_math_time_gate_node)(
                &plan->compiled_math,
                node,
                lower_time_id,
                lower_time_ids[other],
                CompiledMathValueKind::Scalar);
      }
      for (std::size_t other = 0; other < upper_time_ids.size(); ++other) {
        if (other == upper_idx) {
          continue;
        }
        node =
            (other < upper_idx ? compiled_math_strict_time_gate_node
                               : compiled_math_time_gate_node)(
                &plan->compiled_math,
                node,
                upper_time_ids[other],
                upper_time_id,
                CompiledMathValueKind::Scalar);
      }
      candidates.push_back(node);
    }
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Sum,
      std::move(candidates),
      CompiledMathValueKind::Scalar);
}

enum class ExactOrderRegionDensityBinderKind : std::uint8_t {
  Source = 0,
  Expr = 1
};

struct ExactOrderRegionDensityBinder {
  ExactOrderRegionDensityBinderKind kind{
      ExactOrderRegionDensityBinderKind::Source};
  semantic::Index subject_id{semantic::kInvalidIndex};
  semantic::Index time_id{semantic::kInvalidIndex};
};

struct ExactOrderRegionProjectionBounds {
  std::vector<semantic::Index> lower_time_ids;
  std::vector<semantic::Index> upper_time_ids;
};

struct ExactOrderRegionProjectionCandidate {
  ExactOrderRegionDensityBinder binder;
  ExactOrderRegionProjectionBounds bounds;
};

inline bool exact_order_region_contains_time_id(
    const std::vector<semantic::Index> &time_ids,
    const semantic::Index time_id) {
  return std::find(time_ids.begin(), time_ids.end(), time_id) !=
         time_ids.end();
}

inline void exact_order_region_append_latent_time_id(
    std::vector<semantic::Index> *time_ids,
    const semantic::Index time_id) {
  if (!exact_region_time_is_latent_variable(time_id)) {
    return;
  }
  exact_order_region_append_time_id(time_ids, time_id);
}

inline bool exact_order_region_atom_references_time(
    const ExactRegionAtom &atom,
    const semantic::Index time_id) {
  return (atom.lhs.kind == ExactRegionVarKind::Time &&
          atom.lhs.id == time_id) ||
         (atom.rhs.kind == ExactRegionVarKind::Time &&
          atom.rhs.id == time_id);
}

inline bool exact_order_region_atom_is_density_binder(
    const ExactRegionAtom &atom,
    const ExactOrderRegionDensityBinder &binder) {
  if (binder.kind == ExactOrderRegionDensityBinderKind::Source) {
    return atom.kind == ExactRegionAtomKind::SourceExact &&
           atom.lhs.kind == ExactRegionVarKind::SourceTime &&
           atom.lhs.id == binder.subject_id &&
           atom.rhs.kind == ExactRegionVarKind::Time &&
           atom.rhs.id == binder.time_id;
  }
  return atom.kind == ExactRegionAtomKind::ExprDensity &&
         atom.lhs.kind == ExactRegionVarKind::ExprTime &&
         atom.lhs.id == binder.subject_id &&
         atom.rhs.kind == ExactRegionVarKind::Time &&
         atom.rhs.id == binder.time_id;
}

inline bool exact_order_region_density_time_used_outside_orders(
    const ExactRegionCell &term,
    const ExactOrderRegionDensityBinder &binder) {
  for (const auto &atom : term.atoms) {
    if (!exact_order_region_atom_references_time(atom, binder.time_id)) {
      continue;
    }
    if (atom.kind == ExactRegionAtomKind::TimeOrder ||
        exact_order_region_atom_is_density_binder(atom, binder)) {
      continue;
    }
    return true;
  }
  for (const auto &equality : term.equalities) {
    const bool touches =
        equality.lhs_time_id == binder.time_id ||
        equality.rhs_time_id == binder.time_id;
    if (touches && equality.lhs_time_id != equality.rhs_time_id) {
      return true;
    }
  }
  return false;
}

inline bool exact_order_region_projection_bounds(
    const ExactRegionCell &term,
    const semantic::Index time_id,
    ExactOrderRegionProjectionBounds *out) {
  *out = ExactOrderRegionProjectionBounds{};
  const auto closure = exact_order_region_build_time_closure(term);
  if (closure.impossible) {
    return false;
  }
  const auto projected_time_id =
      exact_order_region_canonical_time(closure, time_id);
  const auto zero_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero);
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto canonical_time_ids =
      exact_order_region_canonical_time_ids(closure);
  for (const auto candidate_time_id : canonical_time_ids) {
    if (candidate_time_id == projected_time_id) {
      continue;
    }
    const auto before_projected =
        exact_order_region_canonical_relation(
            closure, candidate_time_id, projected_time_id);
    if (before_projected != 0U && candidate_time_id != zero_time_id) {
      exact_order_region_append_time_id(
          &out->lower_time_ids, candidate_time_id);
    }
    const auto after_projected =
        exact_order_region_canonical_relation(
            closure, projected_time_id, candidate_time_id);
    if (after_projected != 0U) {
      exact_order_region_append_time_id(
          &out->upper_time_ids, candidate_time_id);
    }
  }
  exact_order_region_append_time_id(&out->upper_time_ids, observed_time_id);
  exact_order_region_reduce_bounds(term, &out->lower_time_ids, true);
  exact_order_region_reduce_bounds(term, &out->upper_time_ids, false);
  return true;
}

inline std::size_t exact_order_region_projection_latent_dependency_count(
    const ExactOrderRegionProjectionBounds &bounds,
    std::vector<semantic::Index> *latent_time_ids = nullptr) {
  std::vector<semantic::Index> local;
  auto &time_ids = latent_time_ids == nullptr ? local : *latent_time_ids;
  for (const auto time_id : bounds.lower_time_ids) {
    exact_order_region_append_latent_time_id(&time_ids, time_id);
  }
  for (const auto time_id : bounds.upper_time_ids) {
    exact_order_region_append_latent_time_id(&time_ids, time_id);
  }
  return time_ids.size();
}

inline bool exact_order_region_density_binder_projectable(
    const ExactRegionCell &term,
    const ExactOrderRegionDensityBinder &binder,
    const std::vector<semantic::Index> &blocked_time_ids,
    ExactOrderRegionProjectionBounds *bounds,
    std::size_t *latent_dependency_count) {
  if (!exact_region_time_is_latent_variable(binder.time_id) ||
      exact_order_region_contains_time_id(blocked_time_ids, binder.time_id) ||
      exact_order_region_density_time_used_outside_orders(term, binder)) {
    return false;
  }
  if (!exact_order_region_projection_bounds(term, binder.time_id, bounds)) {
    return false;
  }
  *latent_dependency_count =
      exact_order_region_projection_latent_dependency_count(*bounds);
  return true;
}

inline bool exact_order_region_time_has_density(
    const ExactRegionCell &term,
    const semantic::Index time_id) {
  for (const auto &exact : exact_region_exact_source_atoms(term)) {
    if (exact.time_id == time_id) {
      return true;
    }
  }
  for (const auto &factor : exact_region_expr_atoms(term)) {
    if (factor.density && factor.time_id == time_id) {
      return true;
    }
  }
  return false;
}

inline bool exact_order_region_explicit_positive_equality(
    const ExactRegionCell &term,
    const semantic::Index lhs_time_id,
    const semantic::Index rhs_time_id) {
  for (const auto &equality : term.equalities) {
    const bool same_pair =
        (equality.lhs_time_id == lhs_time_id &&
         equality.rhs_time_id == rhs_time_id) ||
        (equality.lhs_time_id == rhs_time_id &&
         equality.rhs_time_id == lhs_time_id);
    if (!same_pair) {
      continue;
    }
    if (equality.mass == ExactRegionEqualityMass::PositiveMass ||
        equality.origin == ExactRegionEqualityOrigin::SharedLatentIdentity ||
        equality.origin == ExactRegionEqualityOrigin::ModelTie) {
      return true;
    }
  }
  return false;
}

inline bool exact_order_region_cell_has_positive_measure(
    const ExactRegionCell &term) {
  auto canonical = term;
  exact_order_region_canonicalize_term(&canonical);
  if (canonical.impossible || canonical.sign == 0.0) {
    return false;
  }
  const auto closure = exact_order_region_build_time_closure(canonical);
  if (closure.impossible) {
    return false;
  }
  for (std::size_t i = 0; i < closure.time_ids.size(); ++i) {
    for (std::size_t j = i + 1U; j < closure.time_ids.size(); ++j) {
      const auto ij = exact_order_region_time_relation_at(closure, i, j);
      const auto ji = exact_order_region_time_relation_at(closure, j, i);
      if (ij != 1U || ji != 1U) {
        continue;
      }
      const auto lhs_time_id = closure.time_ids[i];
      const auto rhs_time_id = closure.time_ids[j];
      if (exact_order_region_explicit_positive_equality(
              canonical, lhs_time_id, rhs_time_id)) {
        continue;
      }
      const bool touches_special =
          exact_order_region_time_is_special(lhs_time_id) ||
          exact_order_region_time_is_special(rhs_time_id);
      if (touches_special &&
          (exact_order_region_time_has_density(canonical, lhs_time_id) ||
           exact_order_region_time_has_density(canonical, rhs_time_id))) {
        continue;
      }
      return false;
    }
  }
  return true;
}

inline bool exact_order_region_find_projection_candidate(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    const std::vector<semantic::Index> &blocked_time_ids,
    ExactOrderRegionProjectionCandidate *out) {
  bool found = false;
  std::size_t best_score = static_cast<std::size_t>(-1);
  const auto consider =
      [&](const ExactOrderRegionDensityBinder &binder) {
        ExactOrderRegionProjectionBounds bounds;
        std::size_t score = 0U;
        if (!exact_order_region_density_binder_projectable(
                term, binder, blocked_time_ids, &bounds, &score)) {
          return;
        }
        if (!found || score < best_score ||
            (score == best_score &&
             (binder.kind < out->binder.kind ||
              (binder.kind == out->binder.kind &&
               binder.subject_id < out->binder.subject_id)))) {
          out->binder = binder;
          out->bounds = std::move(bounds);
          best_score = score;
          found = true;
        }
      };
  for (const auto &exact : exact_region_exact_source_atoms(term)) {
    if (exact.source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(exact.source_id) >=
            static_cast<std::size_t>(plan.source_count)) {
      continue;
    }
    consider(
        ExactOrderRegionDensityBinder{
            ExactOrderRegionDensityBinderKind::Source,
            exact.source_id,
            exact.time_id});
  }
  for (const auto &factor : exact_region_expr_atoms(term)) {
    if (!factor.density ||
        factor.expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(factor.expr_id) >=
            plan.expr_kernels.size()) {
      continue;
    }
    consider(
        ExactOrderRegionDensityBinder{
            ExactOrderRegionDensityBinderKind::Expr,
            factor.expr_id,
            factor.time_id});
  }
  return found;
}

inline semantic::Index exact_order_region_projection_node(
    ExactVariantPlan *plan,
    const ExactOrderRegionProjectionCandidate &candidate,
    const semantic::Index source_view_id) {
  if (candidate.binder.kind == ExactOrderRegionDensityBinderKind::Source) {
    return exact_order_region_source_interval_partition_node(
        plan,
        candidate.binder.subject_id,
        candidate.bounds.lower_time_ids,
        candidate.bounds.upper_time_ids,
        source_view_id);
  }
  return exact_order_region_expr_interval_partition_node(
      plan,
      candidate.binder.subject_id,
      candidate.bounds.lower_time_ids,
      candidate.bounds.upper_time_ids,
      source_view_id);
}

inline bool exact_order_region_project_density_binder(
    ExactVariantPlan *plan,
    ExactRegionCell *term,
    const ExactOrderRegionProjectionCandidate &candidate,
    const semantic::Index source_view_id,
    std::vector<semantic::Index> *factors,
    std::vector<semantic::Index> *factor_latent_time_ids,
    std::vector<semantic::Index> *blocked_time_ids) {
  const auto node =
      exact_order_region_projection_node(plan, candidate, source_view_id);
  factors->push_back(node);
  for (const auto time_id : candidate.bounds.lower_time_ids) {
    exact_order_region_append_latent_time_id(
        factor_latent_time_ids, time_id);
    exact_order_region_append_latent_time_id(blocked_time_ids, time_id);
  }
  for (const auto time_id : candidate.bounds.upper_time_ids) {
    exact_order_region_append_latent_time_id(
        factor_latent_time_ids, time_id);
    exact_order_region_append_latent_time_id(blocked_time_ids, time_id);
  }
  bool removed = false;
  if (candidate.binder.kind == ExactOrderRegionDensityBinderKind::Source) {
    removed =
        exact_order_region_remove_exact_source_time(
            term,
            candidate.binder.subject_id,
            candidate.binder.time_id);
  } else {
    removed =
        exact_order_region_remove_expr_density_time(
            term,
            candidate.binder.subject_id,
            candidate.binder.time_id);
  }
  if (!removed) {
    return false;
  }
  exact_order_region_remove_orders_touching_time(
      term, candidate.binder.time_id);
  exact_order_region_canonicalize_term(term);
  return !term->impossible;
}

inline bool exact_order_region_project_density_binders(
    ExactVariantPlan *plan,
    ExactRegionCell *term,
    const semantic::Index source_view_id,
    std::vector<semantic::Index> *factors,
    std::vector<semantic::Index> *factor_latent_time_ids) {
  std::vector<semantic::Index> blocked_time_ids;
  while (true) {
    ExactOrderRegionProjectionCandidate candidate;
    if (!exact_order_region_find_projection_candidate(
            *plan, *term, blocked_time_ids, &candidate)) {
      break;
    }
    if (!exact_order_region_project_density_binder(
            plan,
            term,
            candidate,
            source_view_id,
            factors,
            factor_latent_time_ids,
            &blocked_time_ids)) {
      return false;
    }
  }
  return true;
}

inline bool exact_order_region_lower_term_node(
    ExactVariantPlan *plan,
    const ExactRegionCell &term,
    const semantic::Index source_view_id,
    semantic::Index *out_node_id,
    ExactRegionCell *out_residual,
    std::vector<semantic::Index> *out_factor_latent_time_ids) {
  if (term.impossible || term.sign == 0.0) {
    return false;
  }
  struct SourceBounds {
    semantic::Index exact_time_id{semantic::kInvalidIndex};
    std::vector<semantic::Index> lower_time_ids;
    std::vector<semantic::Index> upper_time_ids;
  };
  struct ExprBounds {
    semantic::Index density_time_id{semantic::kInvalidIndex};
    std::vector<semantic::Index> lower_time_ids;
    std::vector<semantic::Index> upper_time_ids;
  };
  ExactRegionCell residual = term;
  exact_order_region_canonicalize_term(&residual);
  if (residual.impossible || residual.sign == 0.0) {
    *out_node_id =
        compiled_math_constant(&plan->compiled_math, 0.0);
    if (out_residual != nullptr) {
      *out_residual = std::move(residual);
    }
    if (out_factor_latent_time_ids != nullptr) {
      out_factor_latent_time_ids->clear();
    }
    return true;
  }
  std::vector<semantic::Index> factors;
  std::vector<semantic::Index> factor_latent_time_ids;
  if (!exact_order_region_project_density_binders(
          plan,
          &residual,
          source_view_id,
          &factors,
          &factor_latent_time_ids)) {
    return false;
  }
  std::vector<SourceBounds> bounds(static_cast<std::size_t>(plan->source_count));
  std::vector<ExprBounds> expr_bounds(plan->expr_kernels.size());
  for (const auto &exact : exact_region_exact_source_atoms(residual)) {
    if (exact.source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(exact.source_id) >= bounds.size()) {
      return false;
    }
    auto &source = bounds[static_cast<std::size_t>(exact.source_id)];
    if (source.exact_time_id != semantic::kInvalidIndex &&
        source.exact_time_id != exact.time_id) {
      return false;
    }
    source.exact_time_id = exact.time_id;
  }
  for (const auto &bound : exact_region_lower_source_atoms(residual)) {
    if (bound.source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(bound.source_id) >= bounds.size()) {
      return false;
    }
    auto &source = bounds[static_cast<std::size_t>(bound.source_id)];
    if (source.exact_time_id != semantic::kInvalidIndex) {
      continue;
    }
    source.lower_time_ids.push_back(bound.time_id);
  }
  for (const auto &bound : exact_region_upper_source_atoms(residual)) {
    if (bound.source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(bound.source_id) >= bounds.size()) {
      return false;
    }
    auto &source = bounds[static_cast<std::size_t>(bound.source_id)];
    if (source.exact_time_id != semantic::kInvalidIndex) {
      continue;
    }
    source.upper_time_ids.push_back(bound.time_id);
  }
  for (const auto &factor : exact_region_expr_atoms(residual)) {
    if (factor.expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(factor.expr_id) >= expr_bounds.size()) {
      return false;
    }
    auto &expr = expr_bounds[static_cast<std::size_t>(factor.expr_id)];
    if (factor.density) {
      if (expr.density_time_id != semantic::kInvalidIndex &&
          expr.density_time_id != factor.time_id) {
        return false;
      }
      expr.density_time_id = factor.time_id;
    } else if (factor.before) {
      if (expr.density_time_id == semantic::kInvalidIndex) {
        expr.upper_time_ids.push_back(factor.time_id);
      }
    } else if (expr.density_time_id == semantic::kInvalidIndex) {
      expr.lower_time_ids.push_back(factor.time_id);
    }
  }
  for (const auto &equality : residual.equalities) {
    if (equality.lhs_time_id != equality.rhs_time_id) {
      return false;
    }
  }
  for (semantic::Index source_id = 0;
       source_id < static_cast<semantic::Index>(bounds.size());
       ++source_id) {
    auto &source = bounds[static_cast<std::size_t>(source_id)];
    exact_order_region_reduce_bounds(residual, &source.lower_time_ids, true);
    exact_order_region_reduce_bounds(residual, &source.upper_time_ids, false);
    if (source.exact_time_id != semantic::kInvalidIndex) {
      factors.push_back(
          compiled_math_source_node(
              &plan->compiled_math,
              CompiledMathNodeKind::SourcePdf,
              source_id,
              0,
              source.exact_time_id,
              source_view_id));
      continue;
    }
    const auto has_lower = !source.lower_time_ids.empty();
    const auto has_upper = !source.upper_time_ids.empty();
    if (has_lower && has_upper) {
      factors.push_back(
          exact_order_region_source_interval_partition_node(
              plan,
              source_id,
              source.lower_time_ids,
              source.upper_time_ids,
              source_view_id));
      continue;
    }
    if (!has_lower && source.upper_time_ids.size() > 1U) {
      factors.push_back(
          exact_order_region_source_cdf_min_node(
              plan, source_id, source.upper_time_ids, source_view_id));
      continue;
    }
    if (!has_upper && source.lower_time_ids.size() > 1U) {
      factors.push_back(
          exact_order_region_source_survival_max_node(
              plan, source_id, source.lower_time_ids, source_view_id));
      continue;
    }
    if (has_lower) {
      factors.push_back(
          compiled_math_source_node(
              &plan->compiled_math,
              CompiledMathNodeKind::SourceSurvival,
              source_id,
              0,
              source.lower_time_ids.front(),
              source_view_id));
    } else if (has_upper) {
      factors.push_back(
          compiled_math_source_node(
              &plan->compiled_math,
              CompiledMathNodeKind::SourceCdf,
              source_id,
              0,
              source.upper_time_ids.front(),
              source_view_id));
    }
  }
  for (semantic::Index expr_id = 0;
       expr_id < static_cast<semantic::Index>(expr_bounds.size());
       ++expr_id) {
    auto &expr = expr_bounds[static_cast<std::size_t>(expr_id)];
    exact_order_region_reduce_bounds(residual, &expr.lower_time_ids, true);
    exact_order_region_reduce_bounds(residual, &expr.upper_time_ids, false);
    if (expr.density_time_id != semantic::kInvalidIndex) {
      factors.push_back(
          exact_order_region_expr_value_node(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprDensity,
              expr.density_time_id,
              source_view_id));
      continue;
    }
    if (!expr.lower_time_ids.empty() || !expr.upper_time_ids.empty()) {
      factors.push_back(
          exact_order_region_expr_interval_partition_node(
              plan,
              expr_id,
              expr.lower_time_ids,
              expr.upper_time_ids,
              source_view_id));
    }
  }
  for (const auto &outcome_indices :
       exact_region_outcome_atoms(residual, ExactRegionAtomKind::OutcomeUnused)) {
    factors.push_back(
        compile_outcome_subset_unused_node(
            plan, outcome_indices, false));
  }
  for (const auto &outcome_indices :
       exact_region_outcome_atoms(residual, ExactRegionAtomKind::OutcomeUsed)) {
    factors.push_back(
        compile_outcome_subset_unused_node(
            plan, outcome_indices, true));
  }
  semantic::Index node =
      factors.empty()
          ? compiled_math_constant(&plan->compiled_math, 1.0)
          : compiled_math_algebra_node(
                &plan->compiled_math,
                CompiledMathNodeKind::Product,
                std::move(factors),
                CompiledMathValueKind::Scalar);
  for (const auto &order : exact_region_time_order_atoms(residual)) {
    node =
        (order.strict ? compiled_math_strict_time_gate_node
                      : compiled_math_time_gate_node)(
            &plan->compiled_math,
            node,
            order.after_time_id,
            order.before_time_id,
            CompiledMathValueKind::Scalar);
  }
  if (residual.sign < 0.0) {
    node =
        compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Negate,
            node,
            CompiledMathValueKind::Scalar);
  }
  *out_node_id = node;
  if (out_residual != nullptr) {
    *out_residual = std::move(residual);
  }
  if (out_factor_latent_time_ids != nullptr) {
    *out_factor_latent_time_ids = std::move(factor_latent_time_ids);
  }
  return true;
}

inline bool exact_order_region_lower_term_root(
    ExactVariantPlan *plan,
    const ExactRegionCell &term,
    const semantic::Index source_view_id,
    ExactOrderRegionBuilder *builder,
    semantic::Index *out_root_id) {
  if (term.impossible || term.sign == 0.0) {
    *out_root_id = compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 0.0));
    return true;
  }
  semantic::Index node{semantic::kInvalidIndex};
  ExactRegionCell residual;
  std::vector<semantic::Index> factor_latent_time_ids;
  if (!exact_order_region_lower_term_node(
          plan,
          term,
          source_view_id,
          &node,
          &residual,
          &factor_latent_time_ids)) {
    return false;
  }
  std::vector<semantic::Index> latent_time_ids;
  const auto append_latent_time = [&](const semantic::Index time_id) {
    if (exact_region_time_is_latent_variable(time_id) &&
        std::find(
            latent_time_ids.begin(),
            latent_time_ids.end(),
            time_id) == latent_time_ids.end()) {
      latent_time_ids.push_back(time_id);
    }
  };
  for (const auto time_id : factor_latent_time_ids) {
    append_latent_time(time_id);
  }
  for (const auto &exact : exact_region_exact_source_atoms(residual)) {
    append_latent_time(exact.time_id);
  }
  for (const auto &bound : exact_region_lower_source_atoms(residual)) {
    append_latent_time(bound.time_id);
  }
  for (const auto &bound : exact_region_upper_source_atoms(residual)) {
    append_latent_time(bound.time_id);
  }
  for (const auto &factor : exact_region_expr_atoms(residual)) {
    append_latent_time(factor.time_id);
  }
  for (const auto &order : exact_region_time_order_atoms(residual)) {
    append_latent_time(order.before_time_id);
    append_latent_time(order.after_time_id);
  }
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  for (auto it = latent_time_ids.rbegin();
       it != latent_time_ids.rend();
       ++it) {
    const auto integrand_root =
        compiled_math_make_root(&plan->compiled_math, node);
    node =
        compiled_math_raw_integral_zero_to_current_node(
            &plan->compiled_math,
            integrand_root,
            0,
            observed_time_id,
            0,
            *it);
  }
  *out_root_id = compiled_math_make_root(&plan->compiled_math, node);
  return true;
}



} // namespace detail
} // namespace accumulatr::eval
