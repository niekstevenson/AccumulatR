#pragma once

#include <memory>
#include <unordered_map>

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
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const semantic::Index lower_time_id,
    const semantic::Index upper_time_id,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  const auto upper_cdf =
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourceCdf,
          source_id,
          condition_id,
          upper_time_id,
          source_view_id);
  const auto lower_cdf =
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourceCdf,
          source_id,
          condition_id,
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
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  if (upper_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (upper_time_ids.size() == 1U) {
    return compiled_math_source_node(
        &plan->compiled_math,
        CompiledMathNodeKind::SourceCdf,
        source_id,
        condition_id,
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
            condition_id,
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
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  if (lower_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (lower_time_ids.size() == 1U) {
    return compiled_math_source_node(
        &plan->compiled_math,
        CompiledMathNodeKind::SourceSurvival,
        source_id,
        condition_id,
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
            condition_id,
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
    ExactVariantBuildState *plan,
    const semantic::Index source_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  if (lower_time_ids.empty()) {
    return exact_order_region_source_cdf_min_node(
        plan, source_id, upper_time_ids, source_view_id, condition_id);
  }
  if (upper_time_ids.empty()) {
    return exact_order_region_source_survival_max_node(
        plan, source_id, lower_time_ids, source_view_id, condition_id);
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
              source_view_id,
              condition_id);
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
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind kind,
    const semantic::Index time_id,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  return compile_expr_value_node(
      plan,
      expr_id,
      kind,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index exact_order_region_expr_interval_node(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const semantic::Index lower_time_id,
    const semantic::Index upper_time_id,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  const auto upper_cdf =
      exact_order_region_expr_value_node(
          plan,
          expr_id,
          CompiledMathNodeKind::ExprCdf,
          upper_time_id,
          source_view_id,
          condition_id);
  const auto lower_cdf =
      exact_order_region_expr_value_node(
          plan,
          expr_id,
          CompiledMathNodeKind::ExprCdf,
          lower_time_id,
          source_view_id,
          condition_id);
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
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  if (upper_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (upper_time_ids.size() == 1U) {
    return exact_order_region_expr_value_node(
        plan,
        expr_id,
        CompiledMathNodeKind::ExprCdf,
        upper_time_ids.front(),
        source_view_id,
        condition_id);
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
            source_view_id,
            condition_id);
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
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  if (lower_time_ids.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  if (lower_time_ids.size() == 1U) {
    return exact_order_region_expr_value_node(
        plan,
        expr_id,
        CompiledMathNodeKind::ExprSurvival,
        lower_time_ids.front(),
        source_view_id,
        condition_id);
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
            source_view_id,
            condition_id);
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
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const std::vector<semantic::Index> &lower_time_ids,
    const std::vector<semantic::Index> &upper_time_ids,
    const semantic::Index source_view_id,
    const semantic::Index condition_id = 0) {
  if (lower_time_ids.empty()) {
    return exact_order_region_expr_cdf_min_node(
        plan, expr_id, upper_time_ids, source_view_id, condition_id);
  }
  if (upper_time_ids.empty()) {
    return exact_order_region_expr_survival_max_node(
        plan, expr_id, lower_time_ids, source_view_id, condition_id);
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
              source_view_id,
              condition_id);
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
  bool (*context_overlaps_expr)(const ExactVariantBuildState &,
                                const ExactRegionCell &,
                                semantic::Index){nullptr};
  bool (*relation_can_collapse)(const ExactVariantBuildState &,
                                semantic::Index){nullptr};
  bool (*expand_relation)(const ExactVariantBuildState &,
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

struct ExactProjectionPlan;
using ExactProjectionPlanPtr = std::shared_ptr<const ExactProjectionPlan>;

struct ExactProjectionPlan {
  ExactProjectionPlanKind kind{ExactProjectionPlanKind::Terminal};
  ExactProjectionCost cost;
  ExactOrderRegionBuilder builder_after;
  std::vector<ExactProjectionFactor> factors;
  std::vector<ExactProjectionPlanPtr> children;
  ExactRegionCell residual;
};

struct ExactProjectionMemoKey {
  ExactRegionCell term;
  std::vector<semantic::Index> blocked_time_ids;
  std::vector<semantic::Index> projected_latent_time_ids;
  semantic::Index builder_time_id{semantic::kInvalidIndex};
  semantic::Index relation_ops_id{0};

  bool operator==(const ExactProjectionMemoKey &other) const noexcept {
    if (term.sign != other.term.sign ||
        term.impossible != other.term.impossible ||
        term.atoms.size() != other.term.atoms.size() ||
        term.equalities.size() != other.term.equalities.size() ||
        blocked_time_ids != other.blocked_time_ids ||
        projected_latent_time_ids != other.projected_latent_time_ids ||
        builder_time_id != other.builder_time_id ||
        relation_ops_id != other.relation_ops_id) {
      return false;
    }
    for (std::size_t i = 0; i < term.atoms.size(); ++i) {
      if (!exact_region_atom_equal(term.atoms[i], other.term.atoms[i])) {
        return false;
      }
    }
    for (std::size_t i = 0; i < term.equalities.size(); ++i) {
      const auto &lhs = term.equalities[i];
      const auto &rhs = other.term.equalities[i];
      if (lhs.lhs_time_id != rhs.lhs_time_id ||
          lhs.rhs_time_id != rhs.rhs_time_id ||
          lhs.mass != rhs.mass || lhs.origin != rhs.origin) {
        return false;
      }
    }
    return true;
  }
};

struct ExactProjectionMemoKeyHash {
  std::size_t operator()(const ExactProjectionMemoKey &key) const noexcept {
    std::size_t seed = std::hash<double>{}(key.term.sign);
    combine(&seed, static_cast<std::size_t>(key.term.impossible));
    for (const auto &atom : key.term.atoms) {
      combine(&seed, static_cast<std::size_t>(atom.kind));
      combine(&seed, static_cast<std::size_t>(atom.lhs.kind));
      combine(&seed, static_cast<std::size_t>(atom.lhs.id));
      combine(&seed, static_cast<std::size_t>(atom.rhs.kind));
      combine(&seed, static_cast<std::size_t>(atom.rhs.id));
      for (const auto outcome_id : atom.outcome_indices) {
        combine(&seed, static_cast<std::size_t>(outcome_id));
      }
      combine(&seed, atom.outcome_indices.size());
      combine(&seed, static_cast<std::size_t>(atom.inclusive));
      combine(&seed, static_cast<std::size_t>(atom.strict));
    }
    combine(&seed, key.term.atoms.size());
    for (const auto &equality : key.term.equalities) {
      combine(&seed, static_cast<std::size_t>(equality.lhs_time_id));
      combine(&seed, static_cast<std::size_t>(equality.rhs_time_id));
      combine(&seed, static_cast<std::size_t>(equality.mass));
      combine(&seed, static_cast<std::size_t>(equality.origin));
    }
    combine(&seed, key.term.equalities.size());
    for (const auto time_id : key.blocked_time_ids) {
      combine(&seed, static_cast<std::size_t>(time_id));
    }
    combine(&seed, key.blocked_time_ids.size());
    for (const auto time_id : key.projected_latent_time_ids) {
      combine(&seed, static_cast<std::size_t>(time_id));
    }
    combine(&seed, key.projected_latent_time_ids.size());
    combine(&seed, static_cast<std::size_t>(key.builder_time_id));
    combine(&seed, static_cast<std::size_t>(key.relation_ops_id));
    return seed;
  }

private:
  static void combine(std::size_t *seed, const std::size_t value) noexcept {
    *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
  }
};

struct ExactProjectionMemoEntry {
  ExactProjectionPlanPtr plan;
  bool complete{false};
};

struct ExactProjectionPlannerState {
  std::vector<ExactProjectionRelationOps> relation_ops;
  std::unordered_map<
      ExactProjectionMemoKey,
      ExactProjectionMemoEntry,
      ExactProjectionMemoKeyHash>
      plans;
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
  semantic::Index nested_expr_relation_count{0};
  for (const auto &factor : exact_region_expr_atoms(residual)) {
    if (!factor.density) {
      ++nested_expr_relation_count;
    }
  }
  return ExactProjectionCost{
      latent_count,
      latent_count,
      latent_count,
      static_cast<semantic::Index>(1 + nested_expr_relation_count),
      static_cast<semantic::Index>(
          1 + atom_count + residual.equalities.size() +
          nested_expr_relation_count),
      latent_count};
}

inline semantic::Index exact_projection_factor_node(
    ExactVariantBuildState *plan,
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
    const ExactVariantBuildState &plan,
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
    ExactVariantBuildState *plan,
    const ExactOrderRegionProjectionCandidate &candidate,
    const semantic::Index source_view_id,
    const semantic::Index condition_id) {
  if (candidate.binder.kind == ExactOrderRegionDensityBinderKind::Source) {
    return exact_order_region_source_interval_partition_node(
        plan,
        candidate.binder.subject_id,
        candidate.bounds.lower_time_ids,
        candidate.bounds.upper_time_ids,
        source_view_id,
        condition_id);
  }
  return exact_order_region_expr_interval_partition_node(
      plan,
      candidate.binder.subject_id,
      candidate.bounds.lower_time_ids,
      candidate.bounds.upper_time_ids,
      source_view_id,
      condition_id);
}

inline semantic::Index exact_projection_factor_node(
    ExactVariantBuildState *plan,
    const ExactProjectionFactor &factor,
    const semantic::Index source_view_id,
    const semantic::Index condition_id) {
  return exact_order_region_projection_node(
      plan, factor.candidate, source_view_id, condition_id);
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
    const ExactVariantBuildState &plan,
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
  if (exact_region_time_is_latent_variable(factor.time_id)) {
    return true;
  }
  return ops->context_overlaps_expr != nullptr &&
         ops->context_overlaps_expr(plan, term, factor.expr_id);
}

inline std::vector<ExactOrderRegionExprValueFactor>
exact_projection_coupled_relation_factors(
    const ExactVariantBuildState &plan,
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

inline std::vector<ExactOrderRegionExprValueFactor>
exact_projection_materializable_relation_factors(
    const ExactRegionCell &term) {
  std::vector<ExactOrderRegionExprValueFactor> out;
  for (const auto &factor : exact_region_expr_atoms(term)) {
    if (!factor.density) {
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
    const ExactVariantBuildState &plan,
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

inline bool exact_projection_plan_cell_memoized(
    const ExactVariantBuildState &plan,
    ExactRegionCell term,
    std::vector<semantic::Index> blocked_time_ids,
    std::vector<semantic::Index> projected_latent_time_ids,
    const ExactProjectionRelationOps *ops,
    ExactOrderRegionBuilder builder,
    ExactProjectionPlanPtr *out);

inline bool exact_projection_candidate_better(
    const ExactProjectionPlanPtr &candidate,
    const ExactProjectionPlanPtr &best) {
  return candidate != nullptr &&
         (best == nullptr ||
          exact_projection_cost_less(candidate->cost, best->cost));
}

inline ExactProjectionPlanPtr exact_projection_make_terminal_plan(
    ExactRegionCell term,
    const std::vector<semantic::Index> &projected_latent_time_ids,
    const ExactOrderRegionBuilder builder) {
  exact_order_region_canonicalize_term(&term);
  if (term.impossible || term.sign == 0.0) {
    return nullptr;
  }
  auto out = std::make_shared<ExactProjectionPlan>();
  out->kind = ExactProjectionPlanKind::Terminal;
  out->residual = std::move(term);
  out->builder_after = builder;
  out->cost =
      exact_projection_terminal_cost(
          out->residual, projected_latent_time_ids);
  return out;
}

inline ExactProjectionPlanPtr exact_projection_make_zero_plan(
    const ExactOrderRegionBuilder builder) {
  auto out = std::make_shared<ExactProjectionPlan>();
  out->kind = ExactProjectionPlanKind::Sum;
  out->builder_after = builder;
  return out;
}

inline ExactProjectionPlanPtr exact_projection_make_project_plan(
    const ExactVariantBuildState &plan,
    const ExactRegionCell &term,
    const ExactOrderRegionProjectionCandidate &candidate,
    const std::vector<semantic::Index> &blocked_time_ids,
    const std::vector<semantic::Index> &projected_latent_time_ids,
    const ExactProjectionRelationOps *ops,
    const ExactOrderRegionBuilder builder) {
  auto next_term = term;
  auto next_blocked_time_ids = blocked_time_ids;
  if (!exact_projection_apply_density_projection(
          &next_term, candidate, &next_blocked_time_ids)) {
    return nullptr;
  }
  auto next_projected_latent_time_ids = projected_latent_time_ids;
  exact_projection_append_latent_times(
      &next_projected_latent_time_ids, candidate.bounds.lower_time_ids);
  exact_projection_append_latent_times(
      &next_projected_latent_time_ids, candidate.bounds.upper_time_ids);

  ExactProjectionPlanPtr child;
  if (!exact_projection_plan_cell_memoized(
          plan,
          std::move(next_term),
          std::move(next_blocked_time_ids),
          std::move(next_projected_latent_time_ids),
          ops,
          builder,
          &child)) {
    return nullptr;
  }
  const ExactProjectionFactor factor{candidate};
  auto out = std::make_shared<ExactProjectionPlan>();
  out->kind = ExactProjectionPlanKind::Product;
  out->factors.push_back(factor);
  out->children.push_back(std::move(child));
  out->builder_after = out->children.front()->builder_after;
  out->cost =
      exact_projection_cost_sum(
          exact_projection_factor_cost(factor),
          out->children.front()->cost);
  return out;
}

inline ExactProjectionPlanPtr exact_projection_make_materialized_plan(
    const ExactVariantBuildState &plan,
    const ExactRegionCell &term,
    const ExactOrderRegionExprValueFactor &factor,
    const std::vector<semantic::Index> &blocked_time_ids,
    const std::vector<semantic::Index> &projected_latent_time_ids,
    const ExactProjectionRelationOps &ops,
    ExactOrderRegionBuilder builder,
    const ExactProjectionCost *upper_bound) {
  ExactOrderRegionExpr materialized;
  if (!exact_projection_materialize_relation_factor(
          plan, term, factor, ops, &builder, &materialized)) {
    return nullptr;
  }
  materialized = exact_order_region_simplify(std::move(materialized));
  if (materialized.terms.empty()) {
    return exact_projection_make_zero_plan(builder);
  }

  auto out = std::make_shared<ExactProjectionPlan>();
  out->kind = ExactProjectionPlanKind::Sum;
  out->cost.symbolic_cells =
      static_cast<semantic::Index>(materialized.terms.size());
  if (upper_bound != nullptr &&
      !exact_projection_cost_less(out->cost, *upper_bound)) {
    return nullptr;
  }

  for (auto child_term : materialized.terms) {
    ExactProjectionPlanPtr child;
    if (!exact_projection_plan_cell_memoized(
            plan,
            std::move(child_term),
            blocked_time_ids,
            projected_latent_time_ids,
            &ops,
            builder,
            &child)) {
      continue;
    }
    builder = child->builder_after;
    out->cost = exact_projection_cost_sum(out->cost, child->cost);
    out->children.push_back(std::move(child));
    if (upper_bound != nullptr &&
        !exact_projection_cost_less(out->cost, *upper_bound)) {
      return nullptr;
    }
  }
  if (out->children.empty()) {
    return nullptr;
  }
  out->builder_after = builder;
  return out;
}

inline bool exact_projection_relation_ops_equal(
    const ExactProjectionRelationOps &lhs,
    const ExactProjectionRelationOps &rhs) noexcept {
  return lhs.context_overlaps_expr == rhs.context_overlaps_expr &&
         lhs.relation_can_collapse == rhs.relation_can_collapse &&
         lhs.expand_relation == rhs.expand_relation;
}

inline semantic::Index exact_projection_relation_ops_id(
    ExactProjectionPlannerState *state,
    const ExactProjectionRelationOps *ops) {
  if (ops == nullptr) {
    return 0;
  }
  for (std::size_t i = 0; i < state->relation_ops.size(); ++i) {
    if (exact_projection_relation_ops_equal(state->relation_ops[i], *ops)) {
      return static_cast<semantic::Index>(i + 1U);
    }
  }
  state->relation_ops.push_back(*ops);
  return static_cast<semantic::Index>(state->relation_ops.size());
}

inline bool exact_projection_plan_cell_memoized(
    const ExactVariantBuildState &plan,
    ExactRegionCell term,
    std::vector<semantic::Index> blocked_time_ids,
    std::vector<semantic::Index> projected_latent_time_ids,
    const ExactProjectionRelationOps *ops,
    const ExactOrderRegionBuilder builder,
    ExactProjectionPlanPtr *out) {
  exact_order_region_canonicalize_term(&term);
  std::sort(blocked_time_ids.begin(), blocked_time_ids.end());
  blocked_time_ids.erase(
      std::unique(blocked_time_ids.begin(), blocked_time_ids.end()),
      blocked_time_ids.end());
  std::sort(
      projected_latent_time_ids.begin(), projected_latent_time_ids.end());
  projected_latent_time_ids.erase(
      std::unique(
          projected_latent_time_ids.begin(),
          projected_latent_time_ids.end()),
      projected_latent_time_ids.end());

  if (plan.projection_planner == nullptr) {
    plan.projection_planner = std::make_shared<ExactProjectionPlannerState>();
  }
  auto &state = *plan.projection_planner;
  const auto materializable_relations =
      ops != nullptr && ops->expand_relation != nullptr
          ? exact_projection_materializable_relation_factors(term)
          : std::vector<ExactOrderRegionExprValueFactor>{};
  const bool builder_neutral = materializable_relations.empty();
  ExactProjectionMemoKey key{
      term,
      blocked_time_ids,
      projected_latent_time_ids,
      builder_neutral ? semantic::kInvalidIndex : builder.next_time_id,
      exact_projection_relation_ops_id(&state, ops)};

  const auto found = state.plans.find(key);
  if (found != state.plans.end()) {
    if (!found->second.complete || found->second.plan == nullptr) {
      return false;
    }
    if (builder_neutral &&
        found->second.plan->builder_after.next_time_id !=
            builder.next_time_id) {
      auto adapted = std::make_shared<ExactProjectionPlan>(*found->second.plan);
      adapted->builder_after = builder;
      *out = std::move(adapted);
    } else {
      *out = found->second.plan;
    }
    return true;
  }
  state.plans.emplace(key, ExactProjectionMemoEntry{});

  ExactProjectionPlanPtr best;
  if (term.impossible || term.sign == 0.0) {
    best = exact_projection_make_zero_plan(builder);
  } else {
    const auto coupled_relations =
        exact_projection_coupled_relation_factors(plan, term, ops);
    const auto projection_candidates =
        exact_order_region_projection_candidates(
            plan, term, blocked_time_ids);
    if (coupled_relations.empty()) {
      best = exact_projection_make_terminal_plan(
          term, projected_latent_time_ids, builder);
    }

    for (const auto &candidate : projection_candidates) {
      const ExactProjectionFactor factor{candidate};
      const auto lower_bound = exact_projection_factor_cost(factor);
      if (best != nullptr &&
          !exact_projection_cost_less(lower_bound, best->cost)) {
        continue;
      }
      auto projected =
          exact_projection_make_project_plan(
              plan,
              term,
              candidate,
              blocked_time_ids,
              projected_latent_time_ids,
              ops,
              builder);
      if (exact_projection_candidate_better(projected, best)) {
        best = std::move(projected);
      }
    }

    if (ops != nullptr && ops->expand_relation != nullptr) {
      for (const auto &factor : materializable_relations) {
        auto materialized =
            exact_projection_make_materialized_plan(
                plan,
                term,
                factor,
                blocked_time_ids,
                projected_latent_time_ids,
                *ops,
                builder,
                best == nullptr ? nullptr : &best->cost);
        if (exact_projection_candidate_better(materialized, best)) {
          best = std::move(materialized);
        }
      }
    }
  }

  auto stored = state.plans.find(key);
  stored->second.plan = best;
  stored->second.complete = true;
  if (best == nullptr) {
    return false;
  }
  *out = std::move(best);
  return true;
}

inline bool exact_projection_plan_cell(
    const ExactVariantBuildState &plan,
    ExactRegionCell term,
    std::vector<semantic::Index> blocked_time_ids,
    std::vector<semantic::Index> projected_latent_time_ids,
    const ExactProjectionRelationOps *ops,
    const ExactOrderRegionBuilder builder,
    ExactProjectionPlan *out) {
  ExactProjectionPlanPtr shared;
  if (!exact_projection_plan_cell_memoized(
          plan,
          std::move(term),
          std::move(blocked_time_ids),
          std::move(projected_latent_time_ids),
          ops,
          builder,
          &shared)) {
    return false;
  }
  *out = *shared;
  return true;
}

inline bool exact_projection_emit_terminal_node(
    ExactVariantBuildState *plan,
    const ExactRegionCell &term,
    const semantic::Index source_view_id,
    const semantic::Index condition_id,
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
        exact_projection_factor_node(plan, factor, source_view_id, condition_id));
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
              condition_id,
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
              source_view_id,
              condition_id));
      continue;
    }
    if (!has_lower && source.upper_time_ids.size() > 1U) {
      factors.push_back(
          exact_order_region_source_cdf_min_node(
              plan,
              source_id,
              source.upper_time_ids,
              source_view_id,
              condition_id));
      continue;
    }
    if (!has_upper && source.lower_time_ids.size() > 1U) {
      factors.push_back(
          exact_order_region_source_survival_max_node(
              plan,
              source_id,
              source.lower_time_ids,
              source_view_id,
              condition_id));
      continue;
    }
    if (has_lower) {
      factors.push_back(
          compiled_math_source_node(
              &plan->compiled_math,
              CompiledMathNodeKind::SourceSurvival,
              source_id,
              condition_id,
              source.lower_time_ids.front(),
              source_view_id));
    } else if (has_upper) {
      factors.push_back(
          compiled_math_source_node(
              &plan->compiled_math,
              CompiledMathNodeKind::SourceCdf,
              source_id,
              condition_id,
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
              source_view_id,
              condition_id));
      continue;
    }
    if (!expr.lower_time_ids.empty() || !expr.upper_time_ids.empty()) {
      factors.push_back(
          exact_order_region_expr_interval_partition_node(
              plan,
              expr_id,
              expr.lower_time_ids,
              expr.upper_time_ids,
              source_view_id,
              condition_id));
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

inline semantic::Index exact_projection_raw_integral_upper_time(
    const ExactRegionCell &residual,
    const semantic::Index latent_time_id) {
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  std::vector<semantic::Index> upper_time_ids;
  for (const auto &order : exact_region_time_order_atoms(residual)) {
    if (order.before_time_id == latent_time_id) {
      exact_order_region_append_time_id(&upper_time_ids, order.after_time_id);
    }
  }
  if (upper_time_ids.empty()) {
    return observed_time_id;
  }
  exact_order_region_append_time_id(&upper_time_ids, observed_time_id);
  exact_order_region_reduce_bounds(residual, &upper_time_ids, false);
  return upper_time_ids.size() == 1U ? upper_time_ids.front() : observed_time_id;
}

inline bool exact_projection_emit_plan_root(
    ExactVariantBuildState *plan,
    const ExactProjectionPlan &projection_plan,
    const semantic::Index source_view_id,
    const semantic::Index condition_id,
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
        *projection_plan.children.front(),
        source_view_id,
        condition_id,
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
              *child,
              source_view_id,
              condition_id,
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
          condition_id,
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
  for (auto it = latent_time_ids.rbegin();
       it != latent_time_ids.rend();
       ++it) {
    const auto upper_time_id =
        exact_projection_raw_integral_upper_time(residual, *it);
    const auto integrand_root =
        compiled_math_make_root(&plan->compiled_math, node);
    node =
        compiled_math_raw_integral_zero_to_current_node(
            &plan->compiled_math,
            integrand_root,
            0,
            upper_time_id,
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
        *child, inherited_factors, out);
  }
}

inline void exact_projection_collect_metric_cells(
    const ExactProjectionPlan &projection_plan,
    ExactOrderRegionExpr *out) {
  exact_projection_collect_metric_cells_impl(projection_plan, {}, out);
}

inline bool exact_order_region_lower_term_root(
    ExactVariantBuildState *plan,
    const ExactRegionCell &term,
    const semantic::Index source_view_id,
    const semantic::Index condition_id,
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
      condition_id,
      {},
      out_root_id);
}



} // namespace detail
} // namespace accumulatr::eval
