#pragma once

#include "exact_region_compiler.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool exact_expr_distribution_key_equal(
    const ExactExprDistributionKey &lhs,
    const ExactExprDistributionKey &rhs) noexcept {
  return lhs.expr_id == rhs.expr_id &&
         lhs.value_kind == rhs.value_kind &&
         lhs.condition_id == rhs.condition_id &&
         lhs.time_id == rhs.time_id &&
         lhs.source_view_id == rhs.source_view_id;
}

inline bool exact_expr_distribution_region(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index time_id,
    ExactOrderRegionBuilder *builder,
    ExactOrderRegionExpr *out) {
  if (value_kind == CompiledMathNodeKind::ExprDensity) {
    return exact_order_region_expr_at_time(plan, expr_id, time_id, builder, out);
  }
  if (value_kind == CompiledMathNodeKind::ExprCdf) {
    return exact_order_region_expr_satisfied_at_time(
        plan, expr_id, time_id, builder, out);
  }
  if (value_kind == CompiledMathNodeKind::ExprSurvival) {
    if (!exact_expr_completion_monotone(plan, expr_id)) {
      return false;
    }
    return exact_order_region_expr_not_satisfied_at_time(
        plan, expr_id, time_id, builder, out);
  }
  return false;
}

enum class ExactExprDistributionLoweringKind : std::uint8_t {
  Region = 0,
  Independent = 1
};

enum class ExactVirtualExprKind : std::uint8_t {
  True = 0,
  Expr = 1,
  And = 2,
  Or = 3
};

struct ExactVirtualExpr {
  ExactVirtualExprKind kind{ExactVirtualExprKind::True};
  semantic::Index expr_id{semantic::kInvalidIndex};
  std::vector<ExactVirtualExpr> children;
};

struct ExactExprDistributionLowering {
  ExactExprDistributionLoweringKind kind{
      ExactExprDistributionLoweringKind::Region};
  ExactExprDistributionKey key;
  ExactOrderRegionExpr region;
  ExactProjectionCost cost;
  ExactOrderRegionBuilder builder;
  ExactVirtualExpr virtual_expr;
  semantic::Index integral_upper_time_id{semantic::kInvalidIndex};
  bool integrate_density{false};
  bool complement_result{false};
};

inline ExactOrderRegionBuilder exact_expr_distribution_seed_builder(
    const ExactVariantBuildState &plan,
    const semantic::Index condition_id,
    const semantic::Index time_id) {
  ExactOrderRegionBuilder builder;
  exact_order_region_reserve_time(&builder, time_id);
  if (condition_id != 0 &&
      condition_id != semantic::kInvalidIndex) {
    const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
    if (condition_pos < plan.compiled_math.conditions.size()) {
      const auto &condition =
          plan.compiled_math.conditions[condition_pos];
      for (const auto fact_time_id : condition.fact_time_ids) {
        exact_order_region_reserve_time(&builder, fact_time_id);
      }
    }
  }
  return builder;
}

inline bool exact_expr_distribution_prepare_region(
    const ExactVariantBuildState &plan,
    const ExactExprDistributionKey &key,
    ExactExprDistributionLowering *out) {
  auto builder =
      exact_expr_distribution_seed_builder(plan, key.condition_id, key.time_id);
  ExactOrderRegionExpr region;
  if (!exact_expr_distribution_region(
          plan,
          key.expr_id,
          key.value_kind,
          key.time_id,
          &builder,
          &region)) {
    return false;
  }
  region =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(std::move(region)));

  const ExactProjectionRelationOps projection_ops{
      exact_order_region_factor_context_overlaps_expr,
      exact_order_region_expr_relation_can_collapse,
      exact_order_region_expand_relation_factor};

  ExactOrderRegionExpr planned_metric_region;
  ExactProjectionCost cost;
  for (const auto &term : region.terms) {
    ExactProjectionPlan projection_plan;
    if (!exact_projection_plan_cell(
            plan,
            term,
            {},
            {},
            &projection_ops,
            builder,
            &projection_plan)) {
      return false;
    }
    builder = projection_plan.builder_after;
    cost = exact_projection_cost_sum(cost, projection_plan.cost);
    exact_projection_collect_metric_cells(
        projection_plan, &planned_metric_region);
  }
  region =
      exact_order_region_minimize_positive_union(
          exact_order_region_simplify(std::move(planned_metric_region)));
  cost.symbolic_cells = static_cast<semantic::Index>(region.terms.size());
  out->key = key;
  out->region = std::move(region);
  out->cost = cost;
  out->builder = builder;
  out->integrate_density = false;
  out->complement_result = false;
  return true;
}

inline bool exact_expr_distribution_prepare_complement_survival(
    const ExactVariantBuildState &plan,
    const ExactExprDistributionKey &key,
    ExactExprDistributionLowering *out) {
  if (key.value_kind != CompiledMathNodeKind::ExprCdf) {
    return false;
  }
  auto survival_key = key;
  survival_key.value_kind = CompiledMathNodeKind::ExprSurvival;
  ExactExprDistributionLowering candidate;
  if (!exact_expr_distribution_prepare_region(plan, survival_key, &candidate)) {
    return false;
  }
  candidate.key = survival_key;
  candidate.complement_result = true;
  candidate.cost.compiled_nodes += 1;
  *out = std::move(candidate);
  return true;
}

inline bool exact_expr_distribution_prepare_integrated_density(
    const ExactVariantBuildState &plan,
    const ExactExprDistributionKey &key,
    ExactExprDistributionLowering *out) {
  if (key.value_kind != CompiledMathNodeKind::ExprCdf) {
    return false;
  }
  auto density_key = key;
  density_key.value_kind = CompiledMathNodeKind::ExprDensity;
  density_key.time_id =
      key.time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
          ? key.time_id
          : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
  ExactExprDistributionLowering candidate;
  if (!exact_expr_distribution_prepare_region(plan, density_key, &candidate)) {
    return false;
  }
  candidate.key = density_key;
  candidate.integrate_density = true;
  candidate.integral_upper_time_id = key.time_id;
  candidate.cost.integral_nodes += 1;
  candidate.cost.max_integral_depth =
      static_cast<semantic::Index>(candidate.cost.max_integral_depth + 1);
  *out = std::move(candidate);
  return true;
}

inline bool exact_expr_distribution_lowering_less(
    const ExactExprDistributionLowering &lhs,
    const ExactExprDistributionLowering &rhs) {
  return exact_projection_cost_less(lhs.cost, rhs.cost);
}

inline void exact_expr_distribution_collect_same_kind_children(
    const ExactVariantBuildState &plan,
    const semantic::Index expr_id,
    const semantic::ExprKind kind,
    std::vector<semantic::Index> *children) {
  const auto &program = plan.program;
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_kernels.size()) {
    return;
  }
  const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
  if (kernel.kind != kind) {
    children->push_back(expr_id);
    return;
  }
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    exact_expr_distribution_collect_same_kind_children(
        plan,
        program.expr_args[
            static_cast<std::size_t>(kernel.children.offset + i)],
        kind,
        children);
  }
}

inline bool exact_expr_distribution_children_independent(
    const ExactVariantBuildState &plan,
    const std::vector<semantic::Index> &children) {
  for (std::size_t i = 0; i < children.size(); ++i) {
    const auto lhs = children[i];
    if (lhs == semantic::kInvalidIndex ||
        static_cast<std::size_t>(lhs) >= plan.expr_supports.size()) {
      return false;
    }
    for (std::size_t j = i + 1; j < children.size(); ++j) {
      const auto rhs = children[j];
      if (rhs == semantic::kInvalidIndex ||
          static_cast<std::size_t>(rhs) >= plan.expr_supports.size() ||
          supports_overlap(
              plan.expr_supports[static_cast<std::size_t>(lhs)],
              plan.expr_supports[static_cast<std::size_t>(rhs)])) {
        return false;
      }
    }
  }
  return true;
}

inline ExactVirtualExpr exact_virtual_expr_true() {
  return ExactVirtualExpr{};
}

inline ExactVirtualExpr exact_virtual_expr_leaf(const semantic::Index expr_id) {
  ExactVirtualExpr node;
  node.kind = ExactVirtualExprKind::Expr;
  node.expr_id = expr_id;
  return node;
}

inline ExactVirtualExpr exact_virtual_expr_logical(
    const ExactVirtualExprKind kind,
    std::vector<ExactVirtualExpr> children) {
  std::vector<ExactVirtualExpr> normalized;
  normalized.reserve(children.size());
  for (auto &child : children) {
    if (child.kind == ExactVirtualExprKind::True) {
      if (kind == ExactVirtualExprKind::Or) {
        return exact_virtual_expr_true();
      }
      continue;
    }
    if (child.kind == kind) {
      normalized.insert(
          normalized.end(),
          std::make_move_iterator(child.children.begin()),
          std::make_move_iterator(child.children.end()));
      continue;
    }
    normalized.push_back(std::move(child));
  }
  if (normalized.empty()) {
    return exact_virtual_expr_true();
  }
  if (normalized.size() == 1U) {
    return std::move(normalized.front());
  }
  ExactVirtualExpr node;
  node.kind = kind;
  node.children = std::move(normalized);
  return node;
}

inline void exact_expr_distribution_append_support(
    std::vector<semantic::Index> *dst,
    const std::vector<semantic::Index> &src) {
  dst->insert(dst->end(), src.begin(), src.end());
  std::sort(dst->begin(), dst->end());
  dst->erase(std::unique(dst->begin(), dst->end()), dst->end());
}

inline std::vector<semantic::Index> exact_virtual_expr_support(
    const ExactVariantBuildState &plan,
    const ExactVirtualExpr &node) {
  if (node.kind == ExactVirtualExprKind::True) {
    return {};
  }
  if (node.kind == ExactVirtualExprKind::Expr) {
    if (node.expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(node.expr_id) >= plan.expr_supports.size()) {
      return {semantic::kInvalidIndex};
    }
    return plan.expr_supports[static_cast<std::size_t>(node.expr_id)];
  }
  std::vector<semantic::Index> support;
  for (const auto &child : node.children) {
    exact_expr_distribution_append_support(
        &support,
        exact_virtual_expr_support(plan, child));
  }
  return support;
}

inline bool exact_virtual_expr_children_independent(
    const ExactVariantBuildState &plan,
    const std::vector<ExactVirtualExpr> &children) {
  std::vector<std::vector<semantic::Index>> supports;
  supports.reserve(children.size());
  for (const auto &child : children) {
    supports.push_back(exact_virtual_expr_support(plan, child));
    if (std::find(
            supports.back().begin(),
            supports.back().end(),
            semantic::kInvalidIndex) != supports.back().end()) {
      return false;
    }
  }
  for (std::size_t i = 0; i < supports.size(); ++i) {
    for (std::size_t j = i + 1U; j < supports.size(); ++j) {
      if (supports_overlap(supports[i], supports[j])) {
        return false;
      }
    }
  }
  return true;
}

inline bool exact_virtual_expr_distribution_node(
    ExactVariantBuildState *plan,
    const ExactVirtualExpr &node,
    CompiledMathNodeKind value_kind,
    semantic::Index condition_id,
    semantic::Index time_id,
    semantic::Index source_view_id,
    semantic::Index *out_node_id);

inline bool exact_virtual_expr_product_node(
    ExactVariantBuildState *plan,
    const std::vector<ExactVirtualExpr> &children,
    const CompiledMathNodeKind child_value_kind,
    const CompiledMathValueKind product_value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id,
    semantic::Index *out_node_id) {
  std::vector<semantic::Index> factors;
  factors.reserve(children.size());
  for (const auto &child : children) {
    semantic::Index child_node{semantic::kInvalidIndex};
    if (!exact_virtual_expr_distribution_node(
            plan,
            child,
            child_value_kind,
            condition_id,
            time_id,
            source_view_id,
            &child_node)) {
      return false;
    }
    factors.push_back(child_node);
  }
  *out_node_id =
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Product,
          std::move(factors),
          product_value_kind);
  return true;
}

inline bool exact_virtual_expr_density_node(
    ExactVariantBuildState *plan,
    const ExactVirtualExprKind kind,
    const std::vector<ExactVirtualExpr> &children,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id,
    semantic::Index *out_node_id) {
  const auto other_value_kind =
      kind == ExactVirtualExprKind::And
          ? CompiledMathNodeKind::ExprCdf
          : CompiledMathNodeKind::ExprSurvival;
  std::vector<semantic::Index> terms;
  terms.reserve(children.size());
  for (std::size_t i = 0; i < children.size(); ++i) {
    std::vector<semantic::Index> factors;
    factors.reserve(children.size());
    semantic::Index active_density{semantic::kInvalidIndex};
    if (!exact_virtual_expr_distribution_node(
            plan,
            children[i],
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id,
            &active_density)) {
      return false;
    }
    factors.push_back(active_density);
    for (std::size_t j = 0; j < children.size(); ++j) {
      if (i == j) {
        continue;
      }
      semantic::Index sibling_node{semantic::kInvalidIndex};
      if (!exact_virtual_expr_distribution_node(
              plan,
              children[j],
              other_value_kind,
              condition_id,
              time_id,
              source_view_id,
              &sibling_node)) {
        return false;
      }
      factors.push_back(sibling_node);
    }
    terms.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(factors),
            CompiledMathValueKind::Density));
  }
  *out_node_id =
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::CleanSignedSum,
          std::move(terms),
          CompiledMathValueKind::Density);
  return true;
}

inline bool exact_virtual_expr_distribution_node(
    ExactVariantBuildState *plan,
    const ExactVirtualExpr &node,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id,
    semantic::Index *out_node_id) {
  if (node.kind == ExactVirtualExprKind::True) {
    *out_node_id =
        compiled_math_constant(
            &plan->compiled_math,
            value_kind == CompiledMathNodeKind::ExprSurvival ||
                    value_kind == CompiledMathNodeKind::ExprDensity
                ? 0.0
                : 1.0);
    return true;
  }
  if (node.kind == ExactVirtualExprKind::Expr) {
    *out_node_id =
        compile_expr_value_node(
            plan,
            node.expr_id,
            value_kind,
            condition_id,
            time_id,
            source_view_id);
    return true;
  }
  if (!exact_virtual_expr_children_independent(*plan, node.children)) {
    return false;
  }

  if (node.kind == ExactVirtualExprKind::And) {
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return exact_virtual_expr_product_node(
          plan,
          node.children,
          CompiledMathNodeKind::ExprCdf,
          CompiledMathValueKind::Cdf,
          condition_id,
          time_id,
          source_view_id,
          out_node_id);
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      semantic::Index cdf_product{semantic::kInvalidIndex};
      if (!exact_virtual_expr_product_node(
              plan,
              node.children,
              CompiledMathNodeKind::ExprCdf,
              CompiledMathValueKind::Cdf,
              condition_id,
              time_id,
              source_view_id,
              &cdf_product)) {
        return false;
      }
      *out_node_id =
          compiled_math_unary_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Complement,
              cdf_product,
              CompiledMathValueKind::Survival);
      return true;
    }
    if (value_kind == CompiledMathNodeKind::ExprDensity) {
      return exact_virtual_expr_density_node(
          plan,
          node.kind,
          node.children,
          condition_id,
          time_id,
          source_view_id,
          out_node_id);
    }
  }

  if (node.kind == ExactVirtualExprKind::Or) {
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return exact_virtual_expr_product_node(
          plan,
          node.children,
          CompiledMathNodeKind::ExprSurvival,
          CompiledMathValueKind::Survival,
          condition_id,
          time_id,
          source_view_id,
          out_node_id);
    }
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      semantic::Index survival_product{semantic::kInvalidIndex};
      if (!exact_virtual_expr_product_node(
              plan,
              node.children,
              CompiledMathNodeKind::ExprSurvival,
              CompiledMathValueKind::Survival,
              condition_id,
              time_id,
              source_view_id,
              &survival_product)) {
        return false;
      }
      *out_node_id =
          compiled_math_unary_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Complement,
              survival_product,
              CompiledMathValueKind::Cdf);
      return true;
    }
    if (value_kind == CompiledMathNodeKind::ExprDensity) {
      return exact_virtual_expr_density_node(
          plan,
          node.kind,
          node.children,
          condition_id,
          time_id,
          source_view_id,
          out_node_id);
    }
  }
  return false;
}

inline ExactProjectionCost exact_virtual_expr_distribution_cost(
    const ExactVirtualExpr &node,
    const CompiledMathNodeKind value_kind) {
  ExactProjectionCost cost;
  if (node.kind == ExactVirtualExprKind::True ||
      node.kind == ExactVirtualExprKind::Expr) {
    cost.compiled_nodes = 1;
    return cost;
  }

  const auto product_cost =
      [&](const CompiledMathNodeKind child_value_kind,
          const semantic::Index extra_nodes) {
        ExactProjectionCost out;
        out.compiled_nodes = extra_nodes;
        for (const auto &child : node.children) {
          out = exact_projection_cost_sum(
              out,
              exact_virtual_expr_distribution_cost(child, child_value_kind));
        }
        return out;
      };

  if (node.kind == ExactVirtualExprKind::And) {
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return product_cost(CompiledMathNodeKind::ExprCdf, 1);
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return product_cost(CompiledMathNodeKind::ExprCdf, 2);
    }
  }

  if (node.kind == ExactVirtualExprKind::Or) {
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return product_cost(CompiledMathNodeKind::ExprSurvival, 1);
    }
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return product_cost(CompiledMathNodeKind::ExprSurvival, 2);
    }
  }

  if (value_kind == CompiledMathNodeKind::ExprDensity) {
    const auto other_value_kind =
        node.kind == ExactVirtualExprKind::And
            ? CompiledMathNodeKind::ExprCdf
            : CompiledMathNodeKind::ExprSurvival;
    cost.compiled_nodes =
        static_cast<semantic::Index>(1 + node.children.size());
    for (std::size_t i = 0; i < node.children.size(); ++i) {
      cost = exact_projection_cost_sum(
          cost,
          exact_virtual_expr_distribution_cost(
              node.children[i],
              CompiledMathNodeKind::ExprDensity));
      for (std::size_t j = 0; j < node.children.size(); ++j) {
        if (i == j) {
          continue;
        }
        cost = exact_projection_cost_sum(
            cost,
            exact_virtual_expr_distribution_cost(
                node.children[j],
                other_value_kind));
      }
    }
  }
  return cost;
}

inline bool exact_expr_distribution_prepare_independent(
    const ExactVariantBuildState &plan,
    const ExactExprDistributionKey &key,
    ExactExprDistributionLowering *out) {
  if (key.condition_id != 0 || key.source_view_id != 0 ||
      key.expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(key.expr_id) >= plan.expr_kernels.size()) {
    return false;
  }
  const auto &kernel =
      plan.expr_kernels[static_cast<std::size_t>(key.expr_id)];
  if ((kernel.kind != semantic::ExprKind::And &&
       kernel.kind != semantic::ExprKind::Or) ||
      kernel.children.empty()) {
    return false;
  }
  std::vector<semantic::Index> children;
  exact_expr_distribution_collect_same_kind_children(
      plan, key.expr_id, kernel.kind, &children);
  if (children.empty() ||
      !exact_expr_distribution_children_independent(plan, children)) {
    return false;
  }

  std::vector<ExactVirtualExpr> virtual_children;
  virtual_children.reserve(children.size());
  for (const auto child : children) {
    virtual_children.push_back(exact_virtual_expr_leaf(child));
  }
  const auto virtual_kind =
      kernel.kind == semantic::ExprKind::And ? ExactVirtualExprKind::And
                                             : ExactVirtualExprKind::Or;
  auto virtual_expr =
      exact_virtual_expr_logical(virtual_kind, std::move(virtual_children));
  if (virtual_expr.kind == ExactVirtualExprKind::Expr) {
    return false;
  }
  out->kind = ExactExprDistributionLoweringKind::Independent;
  out->key = key;
  out->region = ExactOrderRegionExpr{};
  out->cost =
      exact_virtual_expr_distribution_cost(virtual_expr, key.value_kind);
  out->builder =
      exact_expr_distribution_seed_builder(plan, key.condition_id, key.time_id);
  out->virtual_expr = std::move(virtual_expr);
  out->integrate_density = false;
  out->complement_result = false;
  out->integral_upper_time_id = semantic::kInvalidIndex;
  return true;
}

inline semantic::Index compile_expr_distribution_lowering_root(
    ExactVariantBuildState *plan,
    ExactExprDistributionLowering lowering) {
  semantic::Index node{semantic::kInvalidIndex};
  if (lowering.kind == ExactExprDistributionLoweringKind::Independent) {
    if (!exact_virtual_expr_distribution_node(
            plan,
            lowering.virtual_expr,
            lowering.key.value_kind,
            lowering.key.condition_id,
            lowering.key.time_id,
            lowering.key.source_view_id,
            &node)) {
      throw std::runtime_error(
          "exact independent expression distribution lowering failed");
    }
  } else if (lowering.integrate_density) {
    const auto density_node =
        compile_expr_distribution_node(
            plan,
            lowering.key.expr_id,
            CompiledMathNodeKind::ExprDensity,
            lowering.key.condition_id,
            lowering.key.time_id,
            lowering.key.source_view_id);
    const auto integral_node =
        compile_integral_zero_to_current_node(
            plan,
            density_node,
            lowering.key.condition_id,
            lowering.integral_upper_time_id,
            lowering.key.source_view_id,
            lowering.key.time_id);
    node = integral_node;
  } else {
    exact_complexity_observe_region(plan, lowering.region);

    std::vector<semantic::Index> terms;
    terms.reserve(lowering.region.terms.size());
    const ExactProjectionRelationOps projection_ops{
        exact_order_region_factor_context_overlaps_expr,
        exact_order_region_expr_relation_can_collapse,
        exact_order_region_expand_relation_factor};
    auto builder = lowering.builder;
    for (const auto &term : lowering.region.terms) {
      semantic::Index term_root{semantic::kInvalidIndex};
      if (!exact_order_region_lower_term_root(
              plan,
              term,
              lowering.key.source_view_id,
              lowering.key.condition_id,
              &builder,
              &projection_ops,
              &term_root)) {
        throw std::runtime_error("exact expression term lowering failed");
      }
      terms.push_back(compiled_math_root_node_id(plan->compiled_math, term_root));
    }

    const auto value_kind =
        lowering.key.value_kind == CompiledMathNodeKind::ExprDensity
            ? CompiledMathValueKind::Density
            : lowering.key.value_kind == CompiledMathNodeKind::ExprSurvival
                  ? CompiledMathValueKind::Survival
                  : CompiledMathValueKind::Cdf;
    node =
        terms.empty()
            ? compiled_math_constant(&plan->compiled_math, 0.0)
            : compiled_math_algebra_node(
                  &plan->compiled_math,
                  CompiledMathNodeKind::CleanSignedSum,
                  std::move(terms),
                  value_kind);
  }

  if (lowering.complement_result) {
    node =
        compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Complement,
            node,
            CompiledMathValueKind::Cdf);
  }
  return compiled_math_make_root(&plan->compiled_math, node);
}

inline semantic::Index compile_expr_distribution_root_uncached(
    ExactVariantBuildState *plan,
    const ExactExprDistributionKey &key) {
  ExactExprDistributionLowering best;
  if (!exact_expr_distribution_prepare_region(*plan, key, &best)) {
    throw std::runtime_error("exact expression distribution lowering failed");
  }

  ExactExprDistributionLowering independent;
  if (exact_expr_distribution_prepare_independent(
          *plan, key, &independent) &&
      exact_expr_distribution_lowering_less(independent, best)) {
    best = std::move(independent);
  }

  ExactExprDistributionLowering complement_survival;
  if (exact_expr_distribution_prepare_complement_survival(
          *plan, key, &complement_survival) &&
      exact_expr_distribution_lowering_less(complement_survival, best)) {
    best = std::move(complement_survival);
  }

  ExactExprDistributionLowering integrated_density;
  if (exact_expr_distribution_prepare_integrated_density(
          *plan, key, &integrated_density) &&
      exact_expr_distribution_lowering_less(integrated_density, best)) {
    best = std::move(integrated_density);
  }

  return compile_expr_distribution_lowering_root(plan, std::move(best));
}

inline semantic::Index compile_expr_distribution_node(
    ExactVariantBuildState *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  if (value_kind == CompiledMathNodeKind::ExprSurvival &&
      !exact_expr_completion_monotone(*plan, expr_id)) {
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Complement,
        compile_expr_distribution_node(
            plan,
            expr_id,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id),
        CompiledMathValueKind::Survival);
  }

  ExactExprDistributionKey key;
  key.expr_id = expr_id;
  key.value_kind = value_kind;
  key.condition_id =
      condition_id == semantic::kInvalidIndex ? 0 : condition_id;
  key.time_id = time_id;
  key.source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;

  for (const auto &entry : plan->expr_distributions) {
    if (exact_expr_distribution_key_equal(entry.key, key) &&
        entry.root_id != semantic::kInvalidIndex) {
      return compiled_math_root_node_id(plan->compiled_math, entry.root_id);
    }
  }
  for (const auto &entry : plan->expr_distributions) {
    if (exact_expr_distribution_key_equal(entry.key, key) &&
        entry.compiling) {
      throw std::runtime_error(
          "recursive exact expression distribution compilation");
    }
  }

  const auto entry_idx = plan->expr_distributions.size();
  plan->expr_distributions.push_back(
      ExactExprDistributionPlan{key, semantic::kInvalidIndex, true});
  const auto root_id = compile_expr_distribution_root_uncached(plan, key);
  plan->expr_distributions[entry_idx].root_id = root_id;
  plan->expr_distributions[entry_idx].compiling = false;
  return compiled_math_root_node_id(plan->compiled_math, root_id);
}

} // namespace detail
} // namespace accumulatr::eval
