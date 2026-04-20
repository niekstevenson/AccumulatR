#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "prep_to_semantic.hpp"

namespace accumulatr::compile {

enum class BackendKind : std::uint8_t {
  Direct = 0,
  Exact = 1
};

struct DirectOutcomeSpec {
  std::string label;
  semantic::SourceRef source{};
};

struct VariantCapabilities {
  bool no_surviving_outcomes{false};
  bool ranked_observation{false};
  bool chained_onset{false};
  bool shared_trigger{false};
  bool outcome_remapping{false};
  bool guess_outcome{false};
  bool non_direct_outcome{false};
};

struct CompiledVariant {
  std::string component_id;
  double weight{1.0};
  std::string weight_name;
  semantic::SemanticModel model{};
  BackendKind backend{BackendKind::Exact};
  BackendKind semantic_backend{BackendKind::Exact};
  VariantCapabilities capabilities{};
  std::vector<DirectOutcomeSpec> direct_outcomes;
};

struct CompiledModel {
  std::string component_mode{"fixed"};
  std::string component_reference;
  std::vector<CompiledVariant> variants;
};

namespace detail {

inline std::string to_string(BackendKind backend) {
  return backend == BackendKind::Direct ? "direct" : "exact";
}

inline bool supports_semantic_direct(const VariantCapabilities &capabilities) {
  return !capabilities.no_surviving_outcomes &&
         !capabilities.ranked_observation &&
         !capabilities.chained_onset &&
         !capabilities.shared_trigger &&
         !capabilities.non_direct_outcome;
}

inline bool source_is_dead(const semantic::SourceRef &ref) {
  return !ref.valid();
}

inline bool source_in_component(const semantic::OutcomeSpec &outcome,
                                const std::string &component_id) {
  if (outcome.component_ids.empty()) {
    return true;
  }
  return std::find(
             outcome.component_ids.begin(),
             outcome.component_ids.end(),
             component_id) != outcome.component_ids.end();
}

enum class PoolProjectionKind : std::uint8_t {
  Dead = 0,
  Alias = 1,
  Retained = 2
};

struct PoolProjection {
  PoolProjectionKind kind{PoolProjectionKind::Dead};
  semantic::SourceRef alias_source{};
  int k{1};
  std::vector<semantic::SourceRef> members;
};

enum class ExprResultKind : std::uint8_t {
  Impossible = 0,
  TrueExpr = 1,
  Node = 2
};

struct ExprResult {
  ExprResultKind kind{ExprResultKind::Impossible};
  semantic::Index node_index{semantic::kInvalidIndex};
};

class VariantProjector {
public:
  VariantProjector(const semantic::SemanticModel &model,
                   const semantic::ComponentSpec &component)
      : model_(model),
        component_(component),
        allowed_leaves_(model.leaves.size(), false),
        live_leaves_(model.leaves.size(), false),
        pool_cache_(model.pools.size()),
        pool_ready_(model.pools.size(), false),
        pool_visiting_(model.pools.size(), false) {
    for (const auto leaf_index : component_.active_leaf_indices) {
      if (leaf_index >= 0 &&
          leaf_index < static_cast<semantic::Index>(allowed_leaves_.size())) {
        allowed_leaves_[static_cast<std::size_t>(leaf_index)] = true;
      }
    }
    live_leaves_ = allowed_leaves_;
  }

  CompiledVariant project() {
    resolve_live_leaves();

    std::vector<semantic::ExprNode> temp_expr_nodes;
    std::vector<semantic::OutcomeSpec> projected_outcomes;
    projected_outcomes.reserve(model_.outcomes.size());

    for (const auto &outcome : model_.outcomes) {
      if (!source_in_component(outcome, component_.id)) {
        continue;
      }
      const auto expr_result =
          simplify_expr(outcome.expr_root, &temp_expr_nodes);
      if (expr_result.kind == ExprResultKind::Impossible) {
        continue;
      }

      semantic::OutcomeSpec projected = outcome;
      projected.component_ids.clear();
      projected.expr_root = materialize_expr(expr_result, &temp_expr_nodes);
      projected_outcomes.push_back(std::move(projected));
    }

    std::vector<bool> reachable_leaves(model_.leaves.size(), false);
    std::vector<bool> reachable_pools(model_.pools.size(), false);
    for (const auto &outcome : projected_outcomes) {
      mark_expr_sources(
          outcome.expr_root,
          temp_expr_nodes,
          &reachable_leaves,
          &reachable_pools);
    }

    CompiledVariant variant;
    variant.component_id = component_.id;
    variant.weight = component_.weight;
    variant.weight_name = component_.weight_name;
    variant.model.component_mode = "projected";
    variant.model.component_reference = component_.id;
    variant.model.observation.mode = model_.observation.mode;
    variant.model.observation.n_outcomes =
        component_.n_outcomes_override > 0
            ? component_.n_outcomes_override
            : model_.observation.global_n_outcomes;
    variant.model.observation.global_n_outcomes =
        variant.model.observation.n_outcomes;
    variant.model.components.clear();

    std::vector<semantic::Index> leaf_map(
        model_.leaves.size(), semantic::kInvalidIndex);
    for (semantic::Index i = 0;
         i < static_cast<semantic::Index>(model_.leaves.size()); ++i) {
      if (!live_leaves_[static_cast<std::size_t>(i)] ||
          !reachable_leaves[static_cast<std::size_t>(i)]) {
        continue;
      }
      leaf_map[static_cast<std::size_t>(i)] =
          static_cast<semantic::Index>(variant.model.leaves.size());
      variant.model.leaves.push_back(model_.leaves[static_cast<std::size_t>(i)]);
      variant.model.leaves.back().trigger_index = semantic::kInvalidIndex;
    }

    std::vector<semantic::Index> pool_map(
        model_.pools.size(), semantic::kInvalidIndex);
    for (semantic::Index i = 0;
         i < static_cast<semantic::Index>(model_.pools.size()); ++i) {
      const auto &projection = pool_projection(i);
      if (projection.kind != PoolProjectionKind::Retained ||
          !reachable_pools[static_cast<std::size_t>(i)]) {
        continue;
      }
      pool_map[static_cast<std::size_t>(i)] =
          static_cast<semantic::Index>(variant.model.pools.size());
      variant.model.pools.push_back(model_.pools[static_cast<std::size_t>(i)]);
      variant.model.pools.back().members.clear();
    }

    for (semantic::Index i = 0;
         i < static_cast<semantic::Index>(model_.leaves.size()); ++i) {
      const auto new_index = leaf_map[static_cast<std::size_t>(i)];
      if (new_index == semantic::kInvalidIndex) {
        continue;
      }
      auto &leaf = variant.model.leaves[static_cast<std::size_t>(new_index)];
      leaf.onset = rewrite_onset(model_.leaves[static_cast<std::size_t>(i)].onset,
                                 leaf_map,
                                 pool_map);
    }

    for (semantic::Index i = 0;
         i < static_cast<semantic::Index>(model_.pools.size()); ++i) {
      const auto new_index = pool_map[static_cast<std::size_t>(i)];
      if (new_index == semantic::kInvalidIndex) {
        continue;
      }
      const auto &projection = pool_projection(i);
      auto &pool = variant.model.pools[static_cast<std::size_t>(new_index)];
      pool.k = projection.k;
      for (const auto &member : projection.members) {
        pool.members.push_back(remap_source(member, leaf_map, pool_map));
      }
    }

    variant.model.expr_nodes = temp_expr_nodes;
    for (auto &node : variant.model.expr_nodes) {
      if (node.kind == semantic::ExprKind::Event) {
        node.source = remap_source(node.source, leaf_map, pool_map);
      }
    }
    variant.model.outcomes = projected_outcomes;

    for (const auto &trigger : model_.triggers) {
      std::vector<semantic::Index> members;
      members.reserve(trigger.leaf_indices.size());
      for (const auto leaf_index : trigger.leaf_indices) {
        const auto new_index = leaf_map[static_cast<std::size_t>(leaf_index)];
        if (new_index != semantic::kInvalidIndex) {
          members.push_back(new_index);
        }
      }
      if (members.empty()) {
        continue;
      }

      semantic::TriggerSpec projected = trigger;
      projected.leaf_indices = std::move(members);
      if (projected.kind == semantic::TriggerKind::Shared &&
          projected.leaf_indices.size() < 2U) {
        projected.kind = semantic::TriggerKind::Independent;
      }
      const auto trigger_index =
          static_cast<semantic::Index>(variant.model.triggers.size());
      variant.model.triggers.push_back(projected);
      for (const auto leaf_index : variant.model.triggers.back().leaf_indices) {
        variant.model.leaves[static_cast<std::size_t>(leaf_index)].trigger_index =
            trigger_index;
      }
    }

    classify_backend(&variant);
    return variant;
  }

private:
  const semantic::SemanticModel &model_;
  const semantic::ComponentSpec &component_;
  std::vector<bool> allowed_leaves_;
  std::vector<bool> live_leaves_;
  mutable std::vector<PoolProjection> pool_cache_;
  mutable std::vector<bool> pool_ready_;
  mutable std::vector<bool> pool_visiting_;

  void invalidate_pool_cache() {
    std::fill(pool_ready_.begin(), pool_ready_.end(), false);
    std::fill(pool_visiting_.begin(), pool_visiting_.end(), false);
  }

  void resolve_live_leaves() {
    bool changed = true;
    while (changed) {
      changed = false;
      invalidate_pool_cache();
      for (semantic::Index i = 0;
           i < static_cast<semantic::Index>(model_.leaves.size()); ++i) {
        if (!live_leaves_[static_cast<std::size_t>(i)]) {
          continue;
        }
        if (!leaf_can_survive(i)) {
          live_leaves_[static_cast<std::size_t>(i)] = false;
          changed = true;
        }
      }
    }
    invalidate_pool_cache();
  }

  bool leaf_can_survive(const semantic::Index leaf_index) {
    if (!allowed_leaves_[static_cast<std::size_t>(leaf_index)]) {
      return false;
    }
    const auto &leaf = model_.leaves[static_cast<std::size_t>(leaf_index)];
    if (leaf.onset.kind == semantic::OnsetKind::Absolute) {
      return true;
    }
    return !source_is_dead(project_source_ref(leaf.onset.source));
  }

  const PoolProjection &pool_projection(const semantic::Index pool_index) {
    const auto idx = static_cast<std::size_t>(pool_index);
    if (pool_ready_[idx]) {
      return pool_cache_[idx];
    }
    if (pool_visiting_[idx]) {
      throw std::runtime_error(
          "cyclic pool dependency in projected semantic model");
    }

    pool_visiting_[idx] = true;
    PoolProjection projected;
    const auto &pool = model_.pools[idx];
    projected.k = pool.k;
    projected.members.reserve(pool.members.size());

    for (const auto &member : pool.members) {
      const auto projected_member = project_source_ref(member);
      if (source_is_dead(projected_member)) {
        continue;
      }
      projected.members.push_back(projected_member);
    }

    if (projected.members.size() < static_cast<std::size_t>(pool.k)) {
      projected.kind = PoolProjectionKind::Dead;
    } else if (pool.k == 1 && projected.members.size() == 1U) {
      projected.kind = PoolProjectionKind::Alias;
      projected.alias_source = projected.members.front();
    } else {
      projected.kind = PoolProjectionKind::Retained;
    }

    pool_cache_[idx] = std::move(projected);
    pool_ready_[idx] = true;
    pool_visiting_[idx] = false;
    return pool_cache_[idx];
  }

  semantic::SourceRef project_source_ref(const semantic::SourceRef &ref) {
    switch (ref.kind) {
    case semantic::SourceKind::Leaf:
      if (ref.index == semantic::kInvalidIndex) {
        return semantic::SourceRef{};
      }
      return live_leaves_[static_cast<std::size_t>(ref.index)]
                 ? ref
                 : semantic::SourceRef{};
    case semantic::SourceKind::Pool: {
      if (ref.index == semantic::kInvalidIndex) {
        return semantic::SourceRef{};
      }
      const auto &projection = pool_projection(ref.index);
      switch (projection.kind) {
      case PoolProjectionKind::Dead:
        return semantic::SourceRef{};
      case PoolProjectionKind::Alias:
        return projection.alias_source;
      case PoolProjectionKind::Retained:
        return ref;
      }
      break;
    }
    case semantic::SourceKind::Special:
      return ref;
    }
    return semantic::SourceRef{};
  }

  semantic::Index add_expr_node(std::vector<semantic::ExprNode> *expr_nodes,
                                semantic::ExprNode node) {
    expr_nodes->push_back(std::move(node));
    return static_cast<semantic::Index>(expr_nodes->size() - 1);
  }

  semantic::Index materialize_expr(const ExprResult result,
                                   std::vector<semantic::ExprNode> *expr_nodes) {
    if (result.kind == ExprResultKind::Node) {
      return result.node_index;
    }
    semantic::ExprNode node;
    node.kind = result.kind == ExprResultKind::TrueExpr
                    ? semantic::ExprKind::TrueExpr
                    : semantic::ExprKind::Impossible;
    return add_expr_node(expr_nodes, std::move(node));
  }

  ExprResult simplify_expr(const semantic::Index expr_index,
                           std::vector<semantic::ExprNode> *expr_nodes) {
    const auto &node = model_.expr_nodes[static_cast<std::size_t>(expr_index)];
    switch (node.kind) {
    case semantic::ExprKind::Event: {
      const auto projected_source = project_source_ref(node.source);
      if (source_is_dead(projected_source)) {
        return ExprResult{ExprResultKind::Impossible, semantic::kInvalidIndex};
      }
      semantic::ExprNode out = node;
      out.source = projected_source;
      return ExprResult{
          ExprResultKind::Node,
          add_expr_node(expr_nodes, std::move(out))};
    }
    case semantic::ExprKind::Impossible:
      return ExprResult{ExprResultKind::Impossible, semantic::kInvalidIndex};
    case semantic::ExprKind::TrueExpr:
      return ExprResult{ExprResultKind::TrueExpr, semantic::kInvalidIndex};
    case semantic::ExprKind::Not: {
      const auto child = simplify_expr(node.children.front(), expr_nodes);
      if (child.kind == ExprResultKind::Impossible) {
        return ExprResult{ExprResultKind::TrueExpr, semantic::kInvalidIndex};
      }
      if (child.kind == ExprResultKind::TrueExpr) {
        return ExprResult{ExprResultKind::Impossible, semantic::kInvalidIndex};
      }
      semantic::ExprNode out = node;
      out.children = {child.node_index};
      return ExprResult{
          ExprResultKind::Node,
          add_expr_node(expr_nodes, std::move(out))};
    }
    case semantic::ExprKind::And: {
      std::vector<semantic::Index> children;
      children.reserve(node.children.size());
      for (const auto child_index : node.children) {
        const auto child = simplify_expr(child_index, expr_nodes);
        if (child.kind == ExprResultKind::Impossible) {
          return ExprResult{
              ExprResultKind::Impossible, semantic::kInvalidIndex};
        }
        if (child.kind == ExprResultKind::TrueExpr) {
          continue;
        }
        children.push_back(child.node_index);
      }
      if (children.empty()) {
        return ExprResult{ExprResultKind::TrueExpr, semantic::kInvalidIndex};
      }
      if (children.size() == 1U) {
        return ExprResult{ExprResultKind::Node, children.front()};
      }
      semantic::ExprNode out = node;
      out.children = std::move(children);
      return ExprResult{
          ExprResultKind::Node,
          add_expr_node(expr_nodes, std::move(out))};
    }
    case semantic::ExprKind::Or: {
      std::vector<semantic::Index> children;
      children.reserve(node.children.size());
      for (const auto child_index : node.children) {
        const auto child = simplify_expr(child_index, expr_nodes);
        if (child.kind == ExprResultKind::TrueExpr) {
          return ExprResult{ExprResultKind::TrueExpr, semantic::kInvalidIndex};
        }
        if (child.kind == ExprResultKind::Impossible) {
          continue;
        }
        children.push_back(child.node_index);
      }
      if (children.empty()) {
        return ExprResult{
            ExprResultKind::Impossible, semantic::kInvalidIndex};
      }
      if (children.size() == 1U) {
        return ExprResult{ExprResultKind::Node, children.front()};
      }
      semantic::ExprNode out = node;
      out.children = std::move(children);
      return ExprResult{
          ExprResultKind::Node,
          add_expr_node(expr_nodes, std::move(out))};
    }
    case semantic::ExprKind::Guard: {
      const auto reference =
          simplify_expr(node.reference_child, expr_nodes);
      if (reference.kind == ExprResultKind::Impossible) {
        return ExprResult{
            ExprResultKind::Impossible, semantic::kInvalidIndex};
      }

      const auto blocker = simplify_expr(node.blocker_child, expr_nodes);
      if (blocker.kind == ExprResultKind::Impossible) {
        return reference;
      }

      std::vector<semantic::Index> unless_children;
      unless_children.reserve(node.unless_children.size());
      for (const auto child_index : node.unless_children) {
        const auto child = simplify_expr(child_index, expr_nodes);
        if (child.kind == ExprResultKind::TrueExpr) {
          return reference;
        }
        if (child.kind == ExprResultKind::Impossible) {
          continue;
        }
        unless_children.push_back(child.node_index);
      }

      if (blocker.kind == ExprResultKind::TrueExpr && unless_children.empty()) {
        return ExprResult{
            ExprResultKind::Impossible, semantic::kInvalidIndex};
      }

      semantic::ExprNode out = node;
      out.reference_child = materialize_expr(reference, expr_nodes);
      out.blocker_child = materialize_expr(blocker, expr_nodes);
      out.unless_children = std::move(unless_children);
      return ExprResult{
          ExprResultKind::Node,
          add_expr_node(expr_nodes, std::move(out))};
    }
    }

    throw std::runtime_error("unknown expression kind in projector");
  }

  void mark_source(const semantic::SourceRef &source,
                   std::vector<bool> *reachable_leaves,
                   std::vector<bool> *reachable_pools) {
    switch (source.kind) {
    case semantic::SourceKind::Leaf: {
      const auto idx = static_cast<std::size_t>(source.index);
      if ((*reachable_leaves)[idx]) {
        return;
      }
      (*reachable_leaves)[idx] = true;
      const auto &leaf = model_.leaves[idx];
      if (leaf.onset.kind != semantic::OnsetKind::Absolute) {
        const auto projected_source = project_source_ref(leaf.onset.source);
        if (!source_is_dead(projected_source)) {
          mark_source(projected_source, reachable_leaves, reachable_pools);
        }
      }
      break;
    }
    case semantic::SourceKind::Pool: {
      const auto idx = static_cast<std::size_t>(source.index);
      if ((*reachable_pools)[idx]) {
        return;
      }
      (*reachable_pools)[idx] = true;
      const auto &projection = pool_projection(source.index);
      for (const auto &member : projection.members) {
        mark_source(member, reachable_leaves, reachable_pools);
      }
      break;
    }
    case semantic::SourceKind::Special:
      break;
    }
  }

  void mark_expr_sources(const semantic::Index expr_index,
                         const std::vector<semantic::ExprNode> &expr_nodes,
                         std::vector<bool> *reachable_leaves,
                         std::vector<bool> *reachable_pools) {
    const auto &node = expr_nodes[static_cast<std::size_t>(expr_index)];
    switch (node.kind) {
    case semantic::ExprKind::Event:
      mark_source(node.source, reachable_leaves, reachable_pools);
      break;
    case semantic::ExprKind::And:
    case semantic::ExprKind::Or:
    case semantic::ExprKind::Not:
      for (const auto child_index : node.children) {
        mark_expr_sources(
            child_index,
            expr_nodes,
            reachable_leaves,
            reachable_pools);
      }
      break;
    case semantic::ExprKind::Guard:
      mark_expr_sources(
          node.reference_child,
          expr_nodes,
          reachable_leaves,
          reachable_pools);
      mark_expr_sources(
          node.blocker_child,
          expr_nodes,
          reachable_leaves,
          reachable_pools);
      for (const auto child_index : node.unless_children) {
        mark_expr_sources(
            child_index,
            expr_nodes,
            reachable_leaves,
            reachable_pools);
      }
      break;
    case semantic::ExprKind::Impossible:
    case semantic::ExprKind::TrueExpr:
      break;
    }
  }

  semantic::SourceRef remap_source(
      const semantic::SourceRef &source,
      const std::vector<semantic::Index> &leaf_map,
      const std::vector<semantic::Index> &pool_map) {
    switch (source.kind) {
    case semantic::SourceKind::Leaf:
      return semantic::SourceRef{
          semantic::SourceKind::Leaf,
          leaf_map[static_cast<std::size_t>(source.index)],
          std::string()};
    case semantic::SourceKind::Pool:
      return semantic::SourceRef{
          semantic::SourceKind::Pool,
          pool_map[static_cast<std::size_t>(source.index)],
          std::string()};
    case semantic::SourceKind::Special:
      return source;
    }
    return semantic::SourceRef{};
  }

  semantic::OnsetSpec rewrite_onset(
      const semantic::OnsetSpec &onset,
      const std::vector<semantic::Index> &leaf_map,
      const std::vector<semantic::Index> &pool_map) {
    if (onset.kind == semantic::OnsetKind::Absolute) {
      return onset;
    }

    semantic::OnsetSpec rewritten = onset;
    const auto projected_source = project_source_ref(onset.source);
    if (source_is_dead(projected_source)) {
      throw std::runtime_error("live projected leaf has dead onset source");
    }
    rewritten.source = remap_source(projected_source, leaf_map, pool_map);
    rewritten.kind = projected_source.kind == semantic::SourceKind::Leaf
                         ? semantic::OnsetKind::AfterLeaf
                         : semantic::OnsetKind::AfterPool;
    return rewritten;
  }

  void classify_backend(CompiledVariant *variant) {
    if (variant->model.outcomes.empty()) {
      variant->capabilities.no_surviving_outcomes = true;
    }
    if (variant->model.observation.n_outcomes > 1) {
      variant->capabilities.ranked_observation = true;
    }
    for (const auto &leaf : variant->model.leaves) {
      if (leaf.onset.kind != semantic::OnsetKind::Absolute) {
        variant->capabilities.chained_onset = true;
      }
    }
    for (const auto &trigger : variant->model.triggers) {
      if (trigger.kind == semantic::TriggerKind::Shared &&
          trigger.leaf_indices.size() > 1U) {
        variant->capabilities.shared_trigger = true;
      }
    }

    for (const auto &outcome : variant->model.outcomes) {
      if (outcome.mapping.maps_to_missing ||
          !outcome.mapping.observed_label.empty()) {
        variant->capabilities.outcome_remapping = true;
      }
      if (outcome.has_guess) {
        variant->capabilities.guess_outcome = true;
      }

      const auto &root =
          variant->model.expr_nodes[static_cast<std::size_t>(outcome.expr_root)];
      if (root.kind != semantic::ExprKind::Event) {
        variant->capabilities.non_direct_outcome = true;
        continue;
      }
      variant->direct_outcomes.push_back(
          DirectOutcomeSpec{outcome.label, root.source});
    }

    variant->semantic_backend =
        supports_semantic_direct(variant->capabilities)
            ? BackendKind::Direct
            : BackendKind::Exact;
    variant->backend =
        (variant->semantic_backend == BackendKind::Direct &&
         !variant->capabilities.outcome_remapping &&
         !variant->capabilities.guess_outcome)
            ? BackendKind::Direct
            : BackendKind::Exact;
    if (variant->backend == BackendKind::Exact) {
      const bool keep_direct_outcomes =
          variant->semantic_backend == BackendKind::Direct;
      if (!keep_direct_outcomes) {
        variant->direct_outcomes.clear();
      }
    }
  }
};

} // namespace detail

inline CompiledModel project_semantic_model(const semantic::SemanticModel &model) {
  CompiledModel compiled;
  compiled.component_mode = model.component_mode;
  compiled.component_reference = model.component_reference;

  if (model.components.empty()) {
    semantic::ComponentSpec component;
    component.id = "__default__";
    component.weight = 1.0;
    component.active_leaf_indices.resize(model.leaves.size());
    for (std::size_t i = 0; i < model.leaves.size(); ++i) {
      component.active_leaf_indices[i] = static_cast<semantic::Index>(i);
    }
    compiled.variants.push_back(detail::VariantProjector(model, component).project());
    return compiled;
  }

  compiled.variants.reserve(model.components.size());
  for (const auto &component : model.components) {
    compiled.variants.push_back(detail::VariantProjector(model, component).project());
  }
  return compiled;
}

inline Rcpp::List to_r_list(const CompiledModel &compiled) {
  Rcpp::List variants(compiled.variants.size());
  for (std::size_t i = 0; i < compiled.variants.size(); ++i) {
    const auto &variant = compiled.variants[i];
    Rcpp::List variant_list = detail::to_r_list(variant.model);

    Rcpp::List direct_outcomes(variant.direct_outcomes.size());
    for (std::size_t j = 0; j < variant.direct_outcomes.size(); ++j) {
      const auto &outcome = variant.direct_outcomes[j];
      std::string source_kind;
      std::string source_id;
      if (outcome.source.kind == semantic::SourceKind::Leaf) {
        source_kind = "leaf";
        source_id = variant.model.leaves[
            static_cast<std::size_t>(outcome.source.index)].id;
      } else if (outcome.source.kind == semantic::SourceKind::Pool) {
        source_kind = "pool";
        source_id = variant.model.pools[
            static_cast<std::size_t>(outcome.source.index)].id;
      } else {
        source_kind = "special";
        source_id = outcome.source.special_id;
      }
      direct_outcomes[j] = Rcpp::List::create(
          Rcpp::Named("label") = outcome.label,
          Rcpp::Named("source_kind") = source_kind,
          Rcpp::Named("source_id") = source_id);
    }

    variant_list["component_id"] = variant.component_id;
    variant_list["weight"] = variant.weight;
    variant_list["weight_name"] = variant.weight_name;
    variant_list["backend"] = detail::to_string(variant.backend);
    variant_list["semantic_backend"] = detail::to_string(variant.semantic_backend);
    variant_list["capabilities"] = Rcpp::List::create(
        Rcpp::Named("no_surviving_outcomes") =
            variant.capabilities.no_surviving_outcomes,
        Rcpp::Named("ranked_observation") =
            variant.capabilities.ranked_observation,
        Rcpp::Named("chained_onset") = variant.capabilities.chained_onset,
        Rcpp::Named("shared_trigger") = variant.capabilities.shared_trigger,
        Rcpp::Named("outcome_remapping") =
            variant.capabilities.outcome_remapping,
        Rcpp::Named("guess_outcome") = variant.capabilities.guess_outcome,
        Rcpp::Named("non_direct_outcome") =
            variant.capabilities.non_direct_outcome);
    variant_list["direct_outcomes"] = direct_outcomes;
    variants[i] = variant_list;
  }

  return Rcpp::List::create(
      Rcpp::Named("component_mode") = compiled.component_mode,
      Rcpp::Named("component_reference") = compiled.component_reference,
      Rcpp::Named("variants") = variants);
}

} // namespace accumulatr::compile
