#pragma once

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "../leaf/dist_kind.hpp"

namespace accumulatr::semantic {

using Index = std::int32_t;
inline constexpr Index kInvalidIndex = -1;

enum class SourceKind : std::uint8_t {
  Leaf = 0,
  Pool = 1,
  Special = 2
};

enum class OnsetKind : std::uint8_t {
  Absolute = 0,
  AfterLeaf = 1,
  AfterPool = 2
};

enum class ExprKind : std::uint8_t {
  Event = 0,
  And = 1,
  Or = 2,
  Not = 3,
  Guard = 4,
  Impossible = 5,
  TrueExpr = 6
};

enum class TriggerKind : std::uint8_t {
  Independent = 0,
  Shared = 1
};

enum class ObservationMode : std::uint8_t {
  TopK = 0
};

struct SourceRef {
  SourceKind kind{SourceKind::Leaf};
  Index index{kInvalidIndex};
  std::string special_id;

  bool valid() const noexcept {
    return kind == SourceKind::Special ? !special_id.empty() : index != kInvalidIndex;
  }
};

struct ParamBinding {
  std::vector<std::string> dist_param_names;
  std::string q_name;
  std::string t0_name;

  bool empty() const noexcept {
    return dist_param_names.empty() && q_name.empty() && t0_name.empty();
  }
};

struct OnsetSpec {
  OnsetKind kind{OnsetKind::Absolute};
  SourceRef source{};
  double absolute_value{0.0};
  double lag{0.0};
};

struct LeafSpec {
  std::string id;
  leaf::DistKind dist{leaf::DistKind::Lognormal};
  OnsetSpec onset{};
  ParamBinding params{};
  Index trigger_index{kInvalidIndex};
};

struct PoolSpec {
  std::string id;
  int k{1};
  std::vector<SourceRef> members;
};

struct TriggerSpec {
  std::string id;
  TriggerKind kind{TriggerKind::Independent};
  std::vector<Index> leaf_indices;
  std::string q_name;
  double fixed_q{0.0};
  bool has_fixed_q{false};
};

struct ExprNode {
  ExprKind kind{ExprKind::Impossible};
  SourceRef source{};
  int event_k{0};
  std::vector<Index> children;
  Index reference_child{kInvalidIndex};
  Index blocker_child{kInvalidIndex};
  std::vector<Index> unless_children;
};

struct OutcomeMapping {
  bool maps_to_missing{false};
  std::string observed_label;
};

struct OutcomeSpec {
  std::string label;
  Index expr_root{kInvalidIndex};
  OutcomeMapping mapping{};
  std::vector<std::string> component_ids;
  bool has_guess{false};
  std::string outcome_class;
};

struct ComponentSpec {
  std::string id;
  std::vector<Index> active_leaf_indices;
  double weight{1.0};
  std::string weight_name;
  int n_outcomes_override{0};
};

struct ObservationSpec {
  ObservationMode mode{ObservationMode::TopK};
  int n_outcomes{1};
  int global_n_outcomes{1};
};

struct SemanticModel {
  std::vector<LeafSpec> leaves;
  std::vector<PoolSpec> pools;
  std::vector<TriggerSpec> triggers;
  std::vector<ExprNode> expr_nodes;
  std::vector<OutcomeSpec> outcomes;
  std::vector<ComponentSpec> components;
  ObservationSpec observation{};
  std::string component_mode{"fixed"};
  std::string component_reference;
};

struct ValidationIssue {
  std::string message;
};

inline std::vector<ValidationIssue> validate_basic(const SemanticModel &model) {
  auto add_issue = [](std::vector<ValidationIssue> *issues, std::string message) {
    issues->push_back(ValidationIssue{std::move(message)});
  };

  auto has_duplicates = [](const auto &values) {
    std::unordered_set<typename std::decay_t<decltype(values)>::value_type> seen;
    for (const auto &value : values) {
      if (!seen.insert(value).second) {
        return true;
      }
    }
    return false;
  };

  auto valid_source_ref = [&model](const SourceRef &ref) {
    if (!ref.valid()) {
      return false;
    }
    switch (ref.kind) {
    case SourceKind::Leaf:
      return ref.index >= 0 &&
             ref.index < static_cast<Index>(model.leaves.size());
    case SourceKind::Pool:
      return ref.index >= 0 &&
             ref.index < static_cast<Index>(model.pools.size());
    case SourceKind::Special:
      return !ref.special_id.empty();
    }
    return false;
  };

  auto valid_expr_index = [&model](Index index) {
    return index >= 0 && index < static_cast<Index>(model.expr_nodes.size());
  };

  std::vector<ValidationIssue> issues;

  if (model.observation.n_outcomes < 1) {
    add_issue(&issues, "observation n_outcomes must be >= 1");
  }
  if (model.observation.global_n_outcomes < 1) {
    add_issue(&issues, "observation global_n_outcomes must be >= 1");
  }

  std::vector<std::string> source_ids;
  source_ids.reserve(model.leaves.size() + model.pools.size());
  for (const auto &leaf : model.leaves) {
    if (leaf.id.empty()) {
      add_issue(&issues, "leaf id must be non-empty");
    } else {
      source_ids.push_back(leaf.id);
    }
  }
  for (const auto &pool : model.pools) {
    if (pool.id.empty()) {
      add_issue(&issues, "pool id must be non-empty");
    } else {
      source_ids.push_back(pool.id);
    }
  }
  if (has_duplicates(source_ids)) {
    add_issue(&issues, "leaf and pool ids must be unique across source nodes");
  }

  std::vector<std::string> trigger_ids;
  trigger_ids.reserve(model.triggers.size());
  for (const auto &trigger : model.triggers) {
    if (trigger.id.empty()) {
      add_issue(&issues, "trigger id must be non-empty");
    } else {
      trigger_ids.push_back(trigger.id);
    }
  }
  if (has_duplicates(trigger_ids)) {
    add_issue(&issues, "trigger ids must be unique");
  }

  std::vector<std::string> outcome_labels;
  outcome_labels.reserve(model.outcomes.size());
  for (const auto &outcome : model.outcomes) {
    if (outcome.label.empty()) {
      add_issue(&issues, "outcome label must be non-empty");
    } else {
      outcome_labels.push_back(outcome.label);
    }
  }
  if (has_duplicates(outcome_labels)) {
    add_issue(&issues, "outcome labels must be unique");
  }

  std::vector<std::string> component_ids;
  component_ids.reserve(model.components.size());
  for (const auto &component : model.components) {
    if (component.id.empty()) {
      add_issue(&issues, "component id must be non-empty");
    } else {
      component_ids.push_back(component.id);
    }
  }
  if (has_duplicates(component_ids)) {
    add_issue(&issues, "component ids must be unique");
  }

  for (Index i = 0; i < static_cast<Index>(model.leaves.size()); ++i) {
    const auto &leaf = model.leaves[static_cast<std::size_t>(i)];
    if (leaf.onset.lag < 0.0) {
      std::ostringstream ss;
      ss << "leaf[" << i << "] onset lag must be >= 0";
      add_issue(&issues, ss.str());
    }
    switch (leaf.onset.kind) {
    case OnsetKind::Absolute:
      break;
    case OnsetKind::AfterLeaf:
      if (!(leaf.onset.source.kind == SourceKind::Leaf &&
            valid_source_ref(leaf.onset.source))) {
        std::ostringstream ss;
        ss << "leaf[" << i
           << "] after-leaf onset must reference a valid leaf source";
        add_issue(&issues, ss.str());
      }
      break;
    case OnsetKind::AfterPool:
      if (!(leaf.onset.source.kind == SourceKind::Pool &&
            valid_source_ref(leaf.onset.source))) {
        std::ostringstream ss;
        ss << "leaf[" << i
           << "] after-pool onset must reference a valid pool source";
        add_issue(&issues, ss.str());
      }
      break;
    }

    if (leaf.trigger_index != kInvalidIndex &&
        (leaf.trigger_index < 0 ||
         leaf.trigger_index >= static_cast<Index>(model.triggers.size()))) {
      std::ostringstream ss;
      ss << "leaf[" << i << "] references an invalid trigger";
      add_issue(&issues, ss.str());
    }
  }

  for (Index i = 0; i < static_cast<Index>(model.pools.size()); ++i) {
    const auto &pool = model.pools[static_cast<std::size_t>(i)];
    if (pool.k < 1) {
      std::ostringstream ss;
      ss << "pool[" << i << "] k must be >= 1";
      add_issue(&issues, ss.str());
    }
    if (pool.members.empty()) {
      std::ostringstream ss;
      ss << "pool[" << i << "] must have at least one member";
      add_issue(&issues, ss.str());
    }
    if (!pool.members.empty() && pool.k > static_cast<int>(pool.members.size())) {
      std::ostringstream ss;
      ss << "pool[" << i << "] k cannot exceed member count";
      add_issue(&issues, ss.str());
    }
    for (const auto &member : pool.members) {
      if (!valid_source_ref(member)) {
        std::ostringstream ss;
        ss << "pool[" << i << "] has invalid member reference";
        add_issue(&issues, ss.str());
      }
    }
  }

  for (Index i = 0; i < static_cast<Index>(model.triggers.size()); ++i) {
    const auto &trigger = model.triggers[static_cast<std::size_t>(i)];
    if (trigger.has_fixed_q && (trigger.fixed_q < 0.0 || trigger.fixed_q > 1.0)) {
      std::ostringstream ss;
      ss << "trigger[" << i << "] fixed_q must be in [0, 1]";
      add_issue(&issues, ss.str());
    }
    for (const auto &leaf_index : trigger.leaf_indices) {
      if (leaf_index < 0 ||
          leaf_index >= static_cast<Index>(model.leaves.size())) {
        std::ostringstream ss;
        ss << "trigger[" << i << "] references an invalid leaf";
        add_issue(&issues, ss.str());
      }
    }
  }

  for (Index i = 0; i < static_cast<Index>(model.expr_nodes.size()); ++i) {
    const auto &node = model.expr_nodes[static_cast<std::size_t>(i)];
    switch (node.kind) {
    case ExprKind::Event:
      if (!valid_source_ref(node.source)) {
        std::ostringstream ss;
        ss << "expr[" << i << "] event node has invalid source";
        add_issue(&issues, ss.str());
      }
      break;
    case ExprKind::And:
    case ExprKind::Or:
      if (node.children.empty()) {
        std::ostringstream ss;
        ss << "expr[" << i << "] logical node must have children";
        add_issue(&issues, ss.str());
      }
      for (const auto child : node.children) {
        if (!valid_expr_index(child)) {
          std::ostringstream ss;
          ss << "expr[" << i << "] has invalid child reference";
          add_issue(&issues, ss.str());
        }
      }
      break;
    case ExprKind::Not:
      if (node.children.size() != 1U) {
        std::ostringstream ss;
        ss << "expr[" << i << "] not node must have exactly one child";
        add_issue(&issues, ss.str());
      } else if (!valid_expr_index(node.children.front())) {
        std::ostringstream ss;
        ss << "expr[" << i << "] not node has invalid child reference";
        add_issue(&issues, ss.str());
      }
      break;
    case ExprKind::Guard:
      if (!valid_expr_index(node.reference_child) ||
          !valid_expr_index(node.blocker_child)) {
        std::ostringstream ss;
        ss << "expr[" << i
           << "] guard node must reference valid reference and blocker children";
        add_issue(&issues, ss.str());
      }
      for (const auto child : node.unless_children) {
        if (!valid_expr_index(child)) {
          std::ostringstream ss;
          ss << "expr[" << i << "] guard node has invalid unless child";
          add_issue(&issues, ss.str());
        }
      }
      break;
    case ExprKind::Impossible:
    case ExprKind::TrueExpr:
      break;
    }
  }

  for (Index i = 0; i < static_cast<Index>(model.outcomes.size()); ++i) {
    const auto &outcome = model.outcomes[static_cast<std::size_t>(i)];
    if (!valid_expr_index(outcome.expr_root)) {
      std::ostringstream ss;
      ss << "outcome[" << i << "] has invalid expr_root";
      add_issue(&issues, ss.str());
    }
    if (has_duplicates(outcome.component_ids)) {
      std::ostringstream ss;
      ss << "outcome[" << i << "] component ids must be unique";
      add_issue(&issues, ss.str());
    }
  }

  for (Index i = 0; i < static_cast<Index>(model.components.size()); ++i) {
    const auto &component = model.components[static_cast<std::size_t>(i)];
    if (component.weight < 0.0) {
      std::ostringstream ss;
      ss << "component[" << i << "] weight must be >= 0";
      add_issue(&issues, ss.str());
    }
    if (component.n_outcomes_override < 0) {
      std::ostringstream ss;
      ss << "component[" << i
         << "] n_outcomes_override must be 0 or >= 1";
      add_issue(&issues, ss.str());
    }
    for (const auto leaf_index : component.active_leaf_indices) {
      if (leaf_index < 0 ||
          leaf_index >= static_cast<Index>(model.leaves.size())) {
        std::ostringstream ss;
        ss << "component[" << i << "] references an invalid leaf";
        add_issue(&issues, ss.str());
      }
    }
  }

  return issues;
}

inline bool is_basic_valid(const SemanticModel &model) {
  return validate_basic(model).empty();
}

} // namespace accumulatr::semantic
