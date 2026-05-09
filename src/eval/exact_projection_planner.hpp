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

struct ExactProjectionRelationOps {
  bool (*context_overlaps_expr)(const ExactVariantPlan &,
                                const ExactRegionCell &,
                                semantic::Index){nullptr};
  bool (*relation_can_collapse)(const ExactVariantPlan &,
                                semantic::Index){nullptr};
  bool (*expand_relation)(const ExactVariantPlan &,
                          const ExactOrderRegionExprValueFactor &,
                          ExactOrderRegionBuilder *,
                          ExactOrderRegionExpr *){nullptr};
};

struct ExactProjectionCost {
  semantic::Index generic_integral_nodes{0};
  semantic::Index max_integral_depth{0};
  semantic::Index integral_nodes{0};
  semantic::Index symbolic_cells{0};
  semantic::Index compiled_nodes{0};
  semantic::Index latent_times{0};
};

inline bool exact_projection_cost_less(const ExactProjectionCost &lhs,
                                       const ExactProjectionCost &rhs) {
  if (lhs.generic_integral_nodes != rhs.generic_integral_nodes) {
    return lhs.generic_integral_nodes < rhs.generic_integral_nodes;
  }
  if (lhs.max_integral_depth != rhs.max_integral_depth) {
    return lhs.max_integral_depth < rhs.max_integral_depth;
  }
  if (lhs.integral_nodes != rhs.integral_nodes) {
    return lhs.integral_nodes < rhs.integral_nodes;
  }
  if (lhs.symbolic_cells != rhs.symbolic_cells) {
    return lhs.symbolic_cells < rhs.symbolic_cells;
  }
  if (lhs.compiled_nodes != rhs.compiled_nodes) {
    return lhs.compiled_nodes < rhs.compiled_nodes;
  }
  return lhs.latent_times < rhs.latent_times;
}

inline ExactProjectionCost exact_projection_cost_sum(
    const ExactProjectionCost &lhs,
    const ExactProjectionCost &rhs) {
  return ExactProjectionCost{
      lhs.generic_integral_nodes + rhs.generic_integral_nodes,
      std::max(lhs.max_integral_depth, rhs.max_integral_depth),
      lhs.integral_nodes + rhs.integral_nodes,
      lhs.symbolic_cells + rhs.symbolic_cells,
      lhs.compiled_nodes + rhs.compiled_nodes,
      lhs.latent_times + rhs.latent_times};
}

enum class ExactProjectionPlanKind : std::uint8_t {
  Terminal = 0,
  Product = 1,
  Sum = 2
};

struct ExactProjectionFactor {
  ExactOrderRegionProjectionCandidate candidate;
};

struct ExactProjectionPlan {
  ExactProjectionPlanKind kind{ExactProjectionPlanKind::Terminal};
  ExactProjectionCost cost;
  ExactOrderRegionBuilder builder_after;
  std::vector<ExactProjectionFactor> factors;
  std::vector<ExactProjectionPlan> children;
  ExactRegionCell residual;
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

inline void exact_projection_append_latent_times(
    std::vector<semantic::Index> *dst,
    const std::vector<semantic::Index> &src) {
  for (const auto time_id : src) {
    exact_order_region_append_latent_time_id(dst, time_id);
  }
}

inline void exact_projection_append_factor_latent_times(
    std::vector<semantic::Index> *dst,
    const ExactProjectionFactor &factor) {
  exact_projection_append_latent_times(
      dst, factor.candidate.bounds.lower_time_ids);
  exact_projection_append_latent_times(
      dst, factor.candidate.bounds.upper_time_ids);
}

inline void exact_projection_append_residual_latent_times(
    std::vector<semantic::Index> *latent_time_ids,
    const ExactRegionCell &residual) {
  for (const auto &exact : exact_region_exact_source_atoms(residual)) {
    exact_order_region_append_latent_time_id(latent_time_ids, exact.time_id);
  }
  for (const auto &bound : exact_region_lower_source_atoms(residual)) {
    exact_order_region_append_latent_time_id(latent_time_ids, bound.time_id);
  }
  for (const auto &bound : exact_region_upper_source_atoms(residual)) {
    exact_order_region_append_latent_time_id(latent_time_ids, bound.time_id);
  }
  for (const auto &factor : exact_region_expr_atoms(residual)) {
    exact_order_region_append_latent_time_id(latent_time_ids, factor.time_id);
  }
  for (const auto &order : exact_region_time_order_atoms(residual)) {
    exact_order_region_append_latent_time_id(
        latent_time_ids, order.before_time_id);
    exact_order_region_append_latent_time_id(
        latent_time_ids, order.after_time_id);
  }
}

inline ExactProjectionCost exact_projection_terminal_cost(
    const ExactRegionCell &residual,
    const std::vector<semantic::Index> &projected_latent_time_ids) {
  std::vector<semantic::Index> latent_time_ids = projected_latent_time_ids;
  exact_projection_append_residual_latent_times(&latent_time_ids, residual);
  const auto latent_count =
      static_cast<semantic::Index>(latent_time_ids.size());
  const auto atom_count =
      static_cast<semantic::Index>(residual.atoms.size());
  return ExactProjectionCost{
      latent_count,
      latent_count,
      latent_count,
      1,
      static_cast<semantic::Index>(1 + atom_count + residual.equalities.size()),
      latent_count};
}

inline semantic::Index exact_projection_factor_node(
    ExactVariantPlan *plan,
    const ExactProjectionFactor &factor,
    const semantic::Index source_view_id);

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

inline std::vector<ExactOrderRegionProjectionCandidate>
exact_order_region_projection_candidates(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    const std::vector<semantic::Index> &blocked_time_ids) {
  std::vector<ExactOrderRegionProjectionCandidate> out;
  const auto consider =
      [&](const ExactOrderRegionDensityBinder &binder) {
        ExactOrderRegionProjectionBounds bounds;
        std::size_t score = 0U;
        if (!exact_order_region_density_binder_projectable(
                term, binder, blocked_time_ids, &bounds, &score)) {
          return;
        }
        out.push_back(
            ExactOrderRegionProjectionCandidate{
                binder, std::move(bounds)});
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
  std::sort(
      out.begin(),
      out.end(),
      [](const auto &lhs, const auto &rhs) {
        const auto lhs_deps =
            exact_order_region_projection_latent_dependency_count(lhs.bounds);
        const auto rhs_deps =
            exact_order_region_projection_latent_dependency_count(rhs.bounds);
        if (lhs_deps != rhs_deps) {
          return lhs_deps < rhs_deps;
        }
        if (lhs.binder.kind != rhs.binder.kind) {
          return lhs.binder.kind < rhs.binder.kind;
        }
        if (lhs.binder.subject_id != rhs.binder.subject_id) {
          return lhs.binder.subject_id < rhs.binder.subject_id;
        }
        return lhs.binder.time_id < rhs.binder.time_id;
      });
  return out;
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

inline semantic::Index exact_projection_factor_node(
    ExactVariantPlan *plan,
    const ExactProjectionFactor &factor,
    const semantic::Index source_view_id) {
  return exact_order_region_projection_node(
      plan, factor.candidate, source_view_id);
}

inline bool exact_projection_apply_density_projection(
    ExactRegionCell *term,
    const ExactOrderRegionProjectionCandidate &candidate,
    std::vector<semantic::Index> *blocked_time_ids) {
  for (const auto time_id : candidate.bounds.lower_time_ids) {
    exact_order_region_append_latent_time_id(blocked_time_ids, time_id);
  }
  for (const auto time_id : candidate.bounds.upper_time_ids) {
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

inline bool exact_projection_relation_factor_coupled(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    const ExactOrderRegionExprValueFactor &factor,
    const ExactProjectionRelationOps *ops) {
  if (ops == nullptr || factor.density) {
    return false;
  }
  if (ops->relation_can_collapse != nullptr &&
      !ops->relation_can_collapse(plan, factor.expr_id)) {
    return true;
  }
  return ops->context_overlaps_expr != nullptr &&
         ops->context_overlaps_expr(plan, term, factor.expr_id);
}

inline std::vector<ExactOrderRegionExprValueFactor>
exact_projection_coupled_relation_factors(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    const ExactProjectionRelationOps *ops) {
  std::vector<ExactOrderRegionExprValueFactor> out;
  for (const auto &factor : exact_region_expr_atoms(term)) {
    if (exact_projection_relation_factor_coupled(plan, term, factor, ops)) {
      out.push_back(factor);
    }
  }
  std::sort(out.begin(), out.end(), exact_order_region_expr_factor_less);
  out.erase(
      std::unique(
          out.begin(),
          out.end(),
          [](const auto &lhs, const auto &rhs) {
            return lhs.expr_id == rhs.expr_id &&
                   lhs.time_id == rhs.time_id &&
                   lhs.before == rhs.before &&
                   lhs.inclusive == rhs.inclusive &&
                   lhs.density == rhs.density;
          }),
      out.end());
  return out;
}

inline bool exact_projection_materialize_relation_factor(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    const ExactOrderRegionExprValueFactor &factor,
    const ExactProjectionRelationOps &ops,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (ops.expand_relation == nullptr || factor.density) {
    return false;
  }
  ExactRegionCell residual = term;
  if (!exact_order_region_remove_expr_value_atom(&residual, factor)) {
    return false;
  }
  exact_order_region_canonicalize_term(&residual);
  if (residual.impossible || residual.sign == 0.0) {
    *out = exact_order_region_zero();
    return true;
  }
  ExactOrderRegionExpr residual_expr;
  residual_expr.terms.push_back(std::move(residual));
  ExactOrderRegionExpr relation;
  if (!ops.expand_relation(plan, factor, builder, &relation)) {
    return false;
  }
  *out =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(
              exact_order_region_conjoin(
                  std::move(residual_expr), std::move(relation))));
  return true;
}

inline ExactProjectionCost exact_projection_factor_cost(
    const ExactProjectionFactor &factor) {
  const auto latent_count =
      static_cast<semantic::Index>(
          exact_order_region_projection_latent_dependency_count(
              factor.candidate.bounds));
  const auto lower_count =
      static_cast<semantic::Index>(
          factor.candidate.bounds.lower_time_ids.size());
  const auto upper_count =
      static_cast<semantic::Index>(
          factor.candidate.bounds.upper_time_ids.size());
  return ExactProjectionCost{
      0,
      0,
      latent_count,
      0,
      static_cast<semantic::Index>(1 + lower_count + upper_count),
      latent_count};
}

inline bool exact_projection_plan_cell(
    const ExactVariantPlan &plan,
    ExactRegionCell term,
    std::vector<semantic::Index> blocked_time_ids,
    std::vector<semantic::Index> projected_latent_time_ids,
    const ExactProjectionRelationOps *ops,
    ExactOrderRegionBuilder builder,
    ExactProjectionPlan *out);

inline bool exact_projection_candidate_better(
    const bool have_best,
    const ExactProjectionPlan &candidate,
    const ExactProjectionPlan &best) {
  return !have_best || exact_projection_cost_less(candidate.cost, best.cost);
}

inline bool exact_projection_make_terminal_plan(
    ExactRegionCell term,
    const std::vector<semantic::Index> &projected_latent_time_ids,
    ExactOrderRegionBuilder builder,
    ExactProjectionPlan *out) {
  exact_order_region_canonicalize_term(&term);
  if (term.impossible || term.sign == 0.0) {
    return false;
  }
  out->kind = ExactProjectionPlanKind::Terminal;
  out->residual = std::move(term);
  out->children.clear();
  out->factors.clear();
  out->builder_after = builder;
  out->cost =
      exact_projection_terminal_cost(
          out->residual, projected_latent_time_ids);
  return true;
}

inline void exact_projection_make_zero_plan(
    ExactOrderRegionBuilder builder,
    ExactProjectionPlan *out) {
  out->kind = ExactProjectionPlanKind::Sum;
  out->cost = ExactProjectionCost{};
  out->builder_after = builder;
  out->factors.clear();
  out->children.clear();
  out->residual = ExactRegionCell{};
}

inline bool exact_projection_make_project_plan(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    const ExactOrderRegionProjectionCandidate &candidate,
    const std::vector<semantic::Index> &blocked_time_ids,
    const std::vector<semantic::Index> &projected_latent_time_ids,
    const ExactProjectionRelationOps *ops,
    ExactOrderRegionBuilder builder,
    ExactProjectionPlan *out) {
  auto next_term = term;
  auto next_blocked_time_ids = blocked_time_ids;
  if (!exact_projection_apply_density_projection(
          &next_term, candidate, &next_blocked_time_ids)) {
    return false;
  }
  auto next_projected_latent_time_ids = projected_latent_time_ids;
  exact_projection_append_latent_times(
      &next_projected_latent_time_ids, candidate.bounds.lower_time_ids);
  exact_projection_append_latent_times(
      &next_projected_latent_time_ids, candidate.bounds.upper_time_ids);

  ExactProjectionPlan child;
  if (!exact_projection_plan_cell(
          plan,
          std::move(next_term),
          std::move(next_blocked_time_ids),
          std::move(next_projected_latent_time_ids),
          ops,
          builder,
          &child)) {
    return false;
  }
  ExactProjectionFactor factor{candidate};
  out->kind = ExactProjectionPlanKind::Product;
  out->factors.clear();
  out->factors.push_back(factor);
  out->children.clear();
  out->children.push_back(std::move(child));
  out->residual = ExactRegionCell{};
  out->builder_after = out->children.front().builder_after;
  out->cost =
      exact_projection_cost_sum(
          exact_projection_factor_cost(factor),
          out->children.front().cost);
  return true;
}

inline bool exact_projection_make_materialized_plan(
    const ExactVariantPlan &plan,
    const ExactRegionCell &term,
    const ExactOrderRegionExprValueFactor &factor,
    const std::vector<semantic::Index> &blocked_time_ids,
    const std::vector<semantic::Index> &projected_latent_time_ids,
    const ExactProjectionRelationOps &ops,
    ExactOrderRegionBuilder builder,
    ExactProjectionPlan *out) {
  ExactOrderRegionExpr materialized;
  if (!exact_projection_materialize_relation_factor(
          plan, term, factor, ops, &builder, &materialized)) {
    return false;
  }
  materialized = exact_order_region_simplify(std::move(materialized));
  if (materialized.terms.empty()) {
    exact_projection_make_zero_plan(builder, out);
    return true;
  }

  out->kind = ExactProjectionPlanKind::Sum;
  out->factors.clear();
  out->children.clear();
  out->residual = ExactRegionCell{};
  out->cost = ExactProjectionCost{};
  out->cost.symbolic_cells =
      static_cast<semantic::Index>(materialized.terms.size());

  bool have_child = false;
  for (auto child_term : materialized.terms) {
    ExactProjectionPlan child;
    if (!exact_projection_plan_cell(
            plan,
            std::move(child_term),
            blocked_time_ids,
            projected_latent_time_ids,
            &ops,
            builder,
            &child)) {
      continue;
    }
    builder = child.builder_after;
    out->cost = exact_projection_cost_sum(out->cost, child.cost);
    out->children.push_back(std::move(child));
    have_child = true;
  }
  if (!have_child) {
    return false;
  }
  out->builder_after = builder;
  return true;
}

inline bool exact_projection_plan_cell(
    const ExactVariantPlan &plan,
    ExactRegionCell term,
    std::vector<semantic::Index> blocked_time_ids,
    std::vector<semantic::Index> projected_latent_time_ids,
    const ExactProjectionRelationOps *ops,
    ExactOrderRegionBuilder builder,
    ExactProjectionPlan *out) {
  exact_order_region_canonicalize_term(&term);
  if (term.impossible || term.sign == 0.0) {
    exact_projection_make_zero_plan(builder, out);
    return true;
  }

  ExactProjectionPlan best;
  bool have_best = false;
  const auto coupled_relations =
      exact_projection_coupled_relation_factors(plan, term, ops);
  const auto projection_candidates =
      exact_order_region_projection_candidates(
          plan, term, blocked_time_ids);
  if (coupled_relations.empty()) {
    ExactProjectionPlan terminal;
    if (exact_projection_make_terminal_plan(
            term, projected_latent_time_ids, builder, &terminal)) {
      best = std::move(terminal);
      have_best = true;
    }
  }

  for (const auto &candidate : projection_candidates) {
    ExactProjectionPlan projected;
    if (!exact_projection_make_project_plan(
            plan,
            term,
            candidate,
            blocked_time_ids,
            projected_latent_time_ids,
            ops,
            builder,
            &projected)) {
      continue;
    }
    if (exact_projection_candidate_better(have_best, projected, best)) {
      best = std::move(projected);
      have_best = true;
    }
  }

  if (ops != nullptr && ops->expand_relation != nullptr) {
    for (const auto &factor : coupled_relations) {
      ExactProjectionPlan materialized;
      if (!exact_projection_make_materialized_plan(
              plan,
              term,
              factor,
              blocked_time_ids,
              projected_latent_time_ids,
              *ops,
              builder,
              &materialized)) {
        continue;
      }
      if (exact_projection_candidate_better(
              have_best, materialized, best)) {
        best = std::move(materialized);
        have_best = true;
      }
    }
  }

  if (!have_best) {
    return false;
  }
  *out = std::move(best);
  return true;
}

inline bool exact_projection_emit_terminal_node(
    ExactVariantPlan *plan,
    const ExactRegionCell &term,
    const semantic::Index source_view_id,
    const std::vector<ExactProjectionFactor> &projected_factors,
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
  factors.reserve(projected_factors.size());
  for (const auto &factor : projected_factors) {
    factors.push_back(
        exact_projection_factor_node(plan, factor, source_view_id));
    exact_projection_append_factor_latent_times(
        &factor_latent_time_ids, factor);
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

inline bool exact_projection_emit_plan_root(
    ExactVariantPlan *plan,
    const ExactProjectionPlan &projection_plan,
    const semantic::Index source_view_id,
    std::vector<ExactProjectionFactor> inherited_factors,
    semantic::Index *out_root_id) {
  if (projection_plan.kind == ExactProjectionPlanKind::Product) {
    inherited_factors.insert(
        inherited_factors.end(),
        projection_plan.factors.begin(),
        projection_plan.factors.end());
    if (projection_plan.children.size() != 1U) {
      return false;
    }
    return exact_projection_emit_plan_root(
        plan,
        projection_plan.children.front(),
        source_view_id,
        std::move(inherited_factors),
        out_root_id);
  }

  if (projection_plan.kind == ExactProjectionPlanKind::Sum) {
    std::vector<semantic::Index> child_nodes;
    child_nodes.reserve(projection_plan.children.size());
    for (const auto &child : projection_plan.children) {
      semantic::Index child_root{semantic::kInvalidIndex};
      if (!exact_projection_emit_plan_root(
              plan,
              child,
              source_view_id,
              inherited_factors,
              &child_root)) {
        return false;
      }
      child_nodes.push_back(
          compiled_math_root_node_id(plan->compiled_math, child_root));
    }
    const auto node =
        child_nodes.empty()
            ? compiled_math_constant(&plan->compiled_math, 0.0)
            : compiled_math_algebra_node(
                  &plan->compiled_math,
                  CompiledMathNodeKind::CleanSignedSum,
                  std::move(child_nodes),
                  CompiledMathValueKind::Scalar);
    *out_root_id = compiled_math_make_root(&plan->compiled_math, node);
    return true;
  }

  const auto &term = projection_plan.residual;
  if (term.impossible || term.sign == 0.0) {
    *out_root_id = compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 0.0));
    return true;
  }
  semantic::Index node{semantic::kInvalidIndex};
  ExactRegionCell residual;
  std::vector<semantic::Index> factor_latent_time_ids;
  if (!exact_projection_emit_terminal_node(
          plan,
          term,
          source_view_id,
          inherited_factors,
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

inline void exact_projection_apply_factor_to_metric_cell(
    ExactRegionCell *term,
    const ExactProjectionFactor &factor) {
  const auto &candidate = factor.candidate;
  if (candidate.binder.kind == ExactOrderRegionDensityBinderKind::Source) {
    exact_order_region_append_exact(
        term,
        candidate.binder.subject_id,
        candidate.binder.time_id);
  } else {
    exact_order_region_append_expr_factor(
        term,
        candidate.binder.subject_id,
        candidate.binder.time_id,
        true,
        true);
  }
  for (const auto lower_time_id : candidate.bounds.lower_time_ids) {
    exact_order_region_append_time_order(
        term, lower_time_id, candidate.binder.time_id);
  }
  for (const auto upper_time_id : candidate.bounds.upper_time_ids) {
    exact_order_region_append_time_order(
        term, candidate.binder.time_id, upper_time_id);
  }
}

inline void exact_projection_collect_metric_cells_impl(
    const ExactProjectionPlan &projection_plan,
    std::vector<ExactProjectionFactor> inherited_factors,
    ExactOrderRegionExpr *out) {
  if (projection_plan.kind == ExactProjectionPlanKind::Terminal) {
    auto term = projection_plan.residual;
    for (const auto &factor : inherited_factors) {
      exact_projection_apply_factor_to_metric_cell(&term, factor);
    }
    exact_order_region_canonicalize_term(&term);
    if (!term.impossible && term.sign != 0.0) {
      out->terms.push_back(std::move(term));
    }
    return;
  }
  if (projection_plan.kind == ExactProjectionPlanKind::Product) {
    inherited_factors.insert(
        inherited_factors.end(),
        projection_plan.factors.begin(),
        projection_plan.factors.end());
  }
  for (const auto &child : projection_plan.children) {
    exact_projection_collect_metric_cells_impl(
        child, inherited_factors, out);
  }
}

inline void exact_projection_collect_metric_cells(
    const ExactProjectionPlan &projection_plan,
    ExactOrderRegionExpr *out) {
  exact_projection_collect_metric_cells_impl(projection_plan, {}, out);
}

inline bool exact_order_region_lower_term_root(
    ExactVariantPlan *plan,
    const ExactRegionCell &term,
    const semantic::Index source_view_id,
    ExactOrderRegionBuilder *builder,
    const ExactProjectionRelationOps *ops,
    semantic::Index *out_root_id) {
  ExactProjectionPlan projection_plan;
  if (!exact_projection_plan_cell(
          *plan,
          term,
          {},
          {},
          ops,
          *builder,
          &projection_plan)) {
    return false;
  }
  *builder = projection_plan.builder_after;
  return exact_projection_emit_plan_root(
      plan,
      projection_plan,
      source_view_id,
      {},
      out_root_id);
}



} // namespace detail
} // namespace accumulatr::eval
