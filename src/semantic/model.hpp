#pragma once

#include <cstdint>
#include <string>
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
  std::string t0_name;

  bool empty() const noexcept {
    return dist_param_names.empty() && t0_name.empty();
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
  std::vector<Index> leaf_indices;
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
  std::vector<Index> competitor_expr_roots;
  std::vector<Index> competitor_outcome_indices;
  std::vector<std::string> component_ids;
  bool has_guess{false};
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

} // namespace accumulatr::semantic
