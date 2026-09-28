#pragma once

#include <functional>
#include <utility>
#include <vector>

#include "../semantic/model.hpp"

namespace accumulatr::compile {

inline void normalize_simulation_expressions(semantic::SemanticModel &model) {
  using semantic::ExprKind;
  using semantic::Index;
  const auto &input = model.expr_nodes;
  std::vector<semantic::ExprNode> nodes;
  nodes.reserve(input.size());
  std::vector<Index> normalized(input.size(), semantic::kInvalidIndex);

  const auto append = [&](semantic::ExprNode node) {
    nodes.push_back(std::move(node));
    return static_cast<Index>(nodes.size() - 1);
  };
  const auto logical = [&](ExprKind kind, std::vector<Index> children) {
    if (children.size() == 1) return children.front();
    semantic::ExprNode node;
    node.kind = children.empty()
                    ? (kind == ExprKind::And ? ExprKind::TrueExpr
                                            : ExprKind::Impossible)
                    : kind;
    node.children = std::move(children);
    return append(std::move(node));
  };
  std::function<void(Index, ExprKind, std::vector<Index> &)> collect =
      [&](Index id, ExprKind kind, std::vector<Index> &children) {
        const auto &node = input[id];
        if (node.kind == kind) {
          for (const auto child : node.children) collect(child, kind, children);
        } else {
          children.push_back(id);
        }
      };
  std::function<Index(Index)> normalize = [&](Index id) -> Index {
    auto &cached = normalized[id];
    if (cached != semantic::kInvalidIndex) return cached;
    auto node = input[id];
    if (node.kind == ExprKind::And || node.kind == ExprKind::Or) {
      std::vector<Index> children;
      collect(id, node.kind, children);
      std::vector<Index> positive, negative;
      for (const auto child : children) {
        if (node.kind == ExprKind::And && input[child].kind == ExprKind::Not) {
          negative.push_back(child);
        } else {
          positive.push_back(normalize(child));
        }
      }
      if (!positive.empty() && !negative.empty()) {
        // Absence is checked when the conjunction completes, not forever.
        // Flatten first so nested conjunctions do not latch absence early.
        for (auto &child : negative) child = normalize(input[child].children[0]);
        node.kind = ExprKind::Guard;
        node.children.clear();
        node.reference_child = logical(ExprKind::And, std::move(positive));
        node.blocker_child = logical(ExprKind::Or, std::move(negative));
      } else {
        // Pure negation has no response-generating transition in likelihood;
        // retain its existing simulator-only absence semantics.
        for (const auto child : negative) positive.push_back(normalize(child));
        return cached = logical(node.kind, std::move(positive));
      }
    } else if (node.kind == ExprKind::Not) {
      node.children[0] = normalize(node.children[0]);
    } else if (node.kind == ExprKind::Guard) {
      node.reference_child = normalize(node.reference_child);
      node.blocker_child = normalize(node.blocker_child);
    }
    return cached = append(std::move(node));
  };
  for (auto &outcome : model.outcomes) {
    outcome.expr_root = normalize(outcome.expr_root);
    for (auto &root : outcome.competitor_expr_roots) root = normalize(root);
  }
  model.expr_nodes = std::move(nodes);
}

} // namespace accumulatr::compile
