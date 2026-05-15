#pragma once

#include <algorithm>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

class ExactExprCanonicalizer {
public:
  explicit ExactExprCanonicalizer(const runtime::ExactEvaluationProgram &program)
      : input_(program),
        canonical_(program.expr_kind.size(), semantic::kInvalidIndex) {
    output_ = input_;
    output_.expr_kind.clear();
    output_.expr_arg_offsets.clear();
    output_.expr_args.clear();
    output_.expr_ref_child.clear();
    output_.expr_blocker_child.clear();
    output_.expr_source_index.clear();
    output_.expr_source_kind.clear();
    output_.expr_source_ids.clear();
    output_.expr_event_k.clear();
    output_.expr_arg_offsets.push_back(0);
  }

  runtime::ExactEvaluationProgram run() {
    for (auto &root : output_.outcome_expr_root) {
      root = canonicalize(root);
    }
    for (auto &root : output_.outcome_competitor_expr_roots) {
      root = canonicalize(root);
    }
    return std::move(output_);
  }

private:
  const runtime::ExactEvaluationProgram &input_;
  runtime::ExactEvaluationProgram output_;
  std::vector<semantic::Index> canonical_;
  std::map<std::vector<semantic::Index>, semantic::Index> interned_;

  semantic::Index canonicalize(const semantic::Index expr_id) {
    if (expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(expr_id) >= input_.expr_kind.size()) {
      return expr_id;
    }
    auto &cached = canonical_[static_cast<std::size_t>(expr_id)];
    if (cached != semantic::kInvalidIndex) {
      return cached;
    }

    const auto kind = input_kind(expr_id);
    switch (kind) {
    case semantic::ExprKind::Event:
      cached = canonical_event(expr_id);
      break;
    case semantic::ExprKind::Impossible:
    case semantic::ExprKind::TrueExpr:
      cached = intern_constant(kind);
      break;
    case semantic::ExprKind::And:
    case semantic::ExprKind::Or:
      cached = canonical_logical(expr_id, kind);
      break;
    case semantic::ExprKind::Not:
      cached = canonical_not(expr_id);
      break;
    case semantic::ExprKind::Guard:
      cached = canonical_guard(expr_id);
      break;
    }
    return cached;
  }

  semantic::ExprKind input_kind(const semantic::Index expr_id) const {
    return static_cast<semantic::ExprKind>(
        input_.expr_kind[static_cast<std::size_t>(expr_id)]);
  }

  semantic::ExprKind output_kind(const semantic::Index expr_id) const {
    return static_cast<semantic::ExprKind>(
        output_.expr_kind[static_cast<std::size_t>(expr_id)]);
  }

  static void normalize_child_set(std::vector<semantic::Index> *children) {
    std::sort(children->begin(), children->end());
    children->erase(std::unique(children->begin(), children->end()),
                    children->end());
  }

  semantic::Index canonical_event(const semantic::Index expr_id) {
    const auto pos = static_cast<std::size_t>(expr_id);
    return intern_event(
        static_cast<semantic::SourceKind>(input_.expr_source_kind[pos]),
        input_.expr_source_index[pos],
        input_.expr_event_k[pos]);
  }

  semantic::Index canonical_not(const semantic::Index expr_id) {
    const auto pos = static_cast<std::size_t>(expr_id);
    if (input_.expr_arg_offsets[pos] + 1 != input_.expr_arg_offsets[pos + 1U]) {
      return semantic::kInvalidIndex;
    }
    const auto child =
        canonicalize(input_.expr_args[
            static_cast<std::size_t>(input_.expr_arg_offsets[pos])]);
    return intern_not(child);
  }

  semantic::Index canonical_guard(const semantic::Index expr_id) {
    const auto pos = static_cast<std::size_t>(expr_id);
    std::vector<semantic::Index> unless_children;
    unless_children.reserve(static_cast<std::size_t>(
        input_.expr_arg_offsets[pos + 1U] - input_.expr_arg_offsets[pos]));
    for (semantic::Index i = input_.expr_arg_offsets[pos];
         i < input_.expr_arg_offsets[pos + 1U];
         ++i) {
      const auto child =
          canonicalize(input_.expr_args[static_cast<std::size_t>(i)]);
      if (child != semantic::kInvalidIndex) {
        unless_children.push_back(child);
      }
    }
    normalize_child_set(&unless_children);
    return intern_guard(
        canonicalize(input_.expr_ref_child[pos]),
        canonicalize(input_.expr_blocker_child[pos]),
        unless_children);
  }

  std::vector<semantic::Index> canonical_children(
      const semantic::Index expr_id,
      const semantic::ExprKind flatten_kind) {
    std::vector<semantic::Index> out;
    const auto pos = static_cast<std::size_t>(expr_id);
    for (semantic::Index i = input_.expr_arg_offsets[pos];
         i < input_.expr_arg_offsets[pos + 1U];
         ++i) {
      const auto child =
          canonicalize(input_.expr_args[static_cast<std::size_t>(i)]);
      if (child == semantic::kInvalidIndex) {
        continue;
      }
      if (output_kind(child) == flatten_kind) {
        const auto child_pos = static_cast<std::size_t>(child);
        for (semantic::Index j = output_.expr_arg_offsets[child_pos];
             j < output_.expr_arg_offsets[child_pos + 1U];
             ++j) {
          out.push_back(output_.expr_args[static_cast<std::size_t>(j)]);
        }
      } else {
        out.push_back(child);
      }
    }
    return out;
  }

  semantic::Index canonical_logical(const semantic::Index expr_id,
                                    const semantic::ExprKind kind) {
    auto children = canonical_children(expr_id, kind);
    std::vector<semantic::Index> normalized;
    normalized.reserve(children.size());
    for (const auto child : children) {
      const auto child_kind = output_kind(child);
      if (kind == semantic::ExprKind::And) {
        if (child_kind == semantic::ExprKind::Impossible) {
          return intern_constant(semantic::ExprKind::Impossible);
        }
        if (child_kind != semantic::ExprKind::TrueExpr) {
          normalized.push_back(child);
        }
      } else {
        if (child_kind == semantic::ExprKind::TrueExpr) {
          return intern_constant(semantic::ExprKind::TrueExpr);
        }
        if (child_kind != semantic::ExprKind::Impossible) {
          normalized.push_back(child);
        }
      }
    }
    normalize_child_set(&normalized);
    if (normalized.empty()) {
      return intern_constant(kind == semantic::ExprKind::And
                                 ? semantic::ExprKind::TrueExpr
                                 : semantic::ExprKind::Impossible);
    }
    if (normalized.size() == 1U) {
      return normalized.front();
    }
    if (kind == semantic::ExprKind::Or) {
      const auto factored = factored_or(normalized);
      if (factored != semantic::kInvalidIndex) {
        return factored;
      }
    }
    return intern_logical(kind, normalized);
  }

  void collect_conjuncts(const semantic::Index expr_id,
                         std::vector<semantic::Index> *out) {
    if (output_kind(expr_id) != semantic::ExprKind::And) {
      out->push_back(expr_id);
      return;
    }
    const auto pos = static_cast<std::size_t>(expr_id);
    for (semantic::Index i = output_.expr_arg_offsets[pos];
         i < output_.expr_arg_offsets[pos + 1U];
         ++i) {
      collect_conjuncts(output_.expr_args[static_cast<std::size_t>(i)], out);
    }
  }

  semantic::Index make_conjunction(std::vector<semantic::Index> children) {
    normalize_child_set(&children);
    if (children.empty()) {
      return intern_constant(semantic::ExprKind::TrueExpr);
    }
    if (children.size() == 1U) {
      return children.front();
    }
    return intern_logical(semantic::ExprKind::And, children);
  }

  semantic::Index make_disjunction(std::vector<semantic::Index> children) {
    normalize_child_set(&children);
    if (children.empty()) {
      return intern_constant(semantic::ExprKind::Impossible);
    }
    if (children.size() == 1U) {
      return children.front();
    }
    const auto factored = factored_or(children);
    if (factored != semantic::kInvalidIndex) {
      return factored;
    }
    return intern_logical(semantic::ExprKind::Or, children);
  }

  semantic::Index factored_or(
      const std::vector<semantic::Index> &branches) {
    if (branches.size() < 2U) {
      return semantic::kInvalidIndex;
    }

    std::vector<std::vector<semantic::Index>> branch_conjuncts;
    branch_conjuncts.reserve(branches.size());
    for (const auto branch : branches) {
      std::vector<semantic::Index> conjuncts;
      collect_conjuncts(branch, &conjuncts);
      normalize_child_set(&conjuncts);
      if (conjuncts.empty()) {
        return semantic::kInvalidIndex;
      }
      branch_conjuncts.push_back(std::move(conjuncts));
    }

    auto common = branch_conjuncts.front();
    for (std::size_t i = 1; i < branch_conjuncts.size(); ++i) {
      std::vector<semantic::Index> next_common;
      std::set_intersection(
          common.begin(),
          common.end(),
          branch_conjuncts[i].begin(),
          branch_conjuncts[i].end(),
          std::back_inserter(next_common));
      common = std::move(next_common);
      if (common.empty()) {
        return semantic::kInvalidIndex;
      }
    }

    std::vector<semantic::Index> residual_branches;
    residual_branches.reserve(branch_conjuncts.size());
    for (const auto &branch : branch_conjuncts) {
      std::vector<semantic::Index> residual;
      std::set_difference(
          branch.begin(),
          branch.end(),
          common.begin(),
          common.end(),
          std::back_inserter(residual));
      if (residual.empty()) {
        return make_conjunction(common);
      }
      residual_branches.push_back(make_conjunction(std::move(residual)));
    }

    auto factored_children = common;
    factored_children.push_back(make_disjunction(std::move(residual_branches)));
    return make_conjunction(std::move(factored_children));
  }

  static void append_key_header(std::vector<semantic::Index> *key,
                                const semantic::ExprKind kind) {
    key->push_back(static_cast<semantic::Index>(kind));
  }

  semantic::Index intern_constant(const semantic::ExprKind kind) {
    std::vector<semantic::Index> key;
    append_key_header(&key, kind);
    return intern_node(
        std::move(key),
        kind,
        {},
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        semantic::SourceKind::Special,
        semantic::kInvalidIndex,
        0);
  }

  semantic::Index intern_event(const semantic::SourceKind source_kind,
                               const semantic::Index source_index,
                               const int event_k) {
    std::vector<semantic::Index> key;
    append_key_header(&key, semantic::ExprKind::Event);
    key.push_back(static_cast<semantic::Index>(source_kind));
    key.push_back(source_index);
    key.push_back(static_cast<semantic::Index>(event_k));
    return intern_node(
        std::move(key),
        semantic::ExprKind::Event,
        {},
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        source_kind,
        source_index,
        event_k);
  }

  semantic::Index intern_logical(const semantic::ExprKind kind,
                                 const std::vector<semantic::Index> &children) {
    std::vector<semantic::Index> key;
    append_key_header(&key, kind);
    key.push_back(static_cast<semantic::Index>(children.size()));
    key.insert(key.end(), children.begin(), children.end());
    return intern_node(
        std::move(key),
        kind,
        children,
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        semantic::SourceKind::Special,
        semantic::kInvalidIndex,
        0);
  }

  semantic::Index intern_not(const semantic::Index child) {
    std::vector<semantic::Index> key;
    append_key_header(&key, semantic::ExprKind::Not);
    key.push_back(child);
    return intern_node(
        std::move(key),
        semantic::ExprKind::Not,
        std::vector<semantic::Index>{child},
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        semantic::SourceKind::Special,
        semantic::kInvalidIndex,
        0);
  }

  semantic::Index intern_guard(const semantic::Index ref,
                               const semantic::Index blocker,
                               const std::vector<semantic::Index> &unless_children) {
    std::vector<semantic::Index> key;
    append_key_header(&key, semantic::ExprKind::Guard);
    key.push_back(ref);
    key.push_back(blocker);
    key.push_back(static_cast<semantic::Index>(unless_children.size()));
    key.insert(key.end(), unless_children.begin(), unless_children.end());
    return intern_node(
        std::move(key),
        semantic::ExprKind::Guard,
        unless_children,
        ref,
        blocker,
        semantic::SourceKind::Special,
        semantic::kInvalidIndex,
        0);
  }

  semantic::Index intern_node(std::vector<semantic::Index> key,
                              const semantic::ExprKind kind,
                              const std::vector<semantic::Index> &children,
                              const semantic::Index ref,
                              const semantic::Index blocker,
                              const semantic::SourceKind source_kind,
                              const semantic::Index source_index,
                              const int event_k) {
    const auto found = interned_.find(key);
    if (found != interned_.end()) {
      return found->second;
    }

    const auto expr_id =
        static_cast<semantic::Index>(output_.expr_kind.size());
    output_.expr_kind.push_back(static_cast<std::uint8_t>(kind));
    output_.expr_args.insert(
        output_.expr_args.end(), children.begin(), children.end());
    output_.expr_arg_offsets.push_back(
        static_cast<semantic::Index>(output_.expr_args.size()));
    output_.expr_ref_child.push_back(ref);
    output_.expr_blocker_child.push_back(blocker);
    output_.expr_source_index.push_back(source_index);
    output_.expr_source_kind.push_back(
        static_cast<std::uint8_t>(source_kind));
    output_.expr_source_ids.push_back(semantic::kInvalidIndex);
    output_.expr_event_k.push_back(event_k);
    interned_.emplace(std::move(key), expr_id);
    return expr_id;
  }
};

inline void canonicalize_exact_evaluation_program_expressions(
    runtime::ExactEvaluationProgram *program) {
  if (program == nullptr) {
    return;
  }
  *program = ExactExprCanonicalizer(*program).run();
}

} // namespace detail
} // namespace accumulatr::eval
