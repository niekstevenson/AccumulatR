#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#if defined(__APPLE__)
#include <vecLib/vForce.h>
#define ACCUMULATR_USE_ACCELERATE_VFORCE 1
#endif

#include "batch_source_product.hpp"
#include "exact_source_channels.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

constexpr std::size_t kFiniteDirectLeafQuadratureTileOrder =
    quadrature::kDefaultFiniteOrder;
constexpr std::size_t kSourceVectorOnsetQuadratureTargetExpandedLanes = 1024U;

#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
inline void compiled_math_batch_normal_cdf_from_z(
    double *z_by_lane,
    double *exp_scratch,
    const std::size_t count) {
  if (z_by_lane == nullptr || exp_scratch == nullptr || count == 0U) {
    return;
  }
  if (count < 16U) {
    for (std::size_t i = 0; i < count; ++i) {
      z_by_lane[i] = normal_cdf_rational_fast(z_by_lane[i]);
    }
    return;
  }
  for (std::size_t i = 0; i < count; ++i) {
    const double z = z_by_lane[i];
    if (!std::isfinite(z)) {
      exp_scratch[i] = 0.0;
      continue;
    }
    const double x = z * kInvSqrtTwo;
    z_by_lane[i] = x;
    exp_scratch[i] = -x * x;
  }
  const auto n = static_cast<int>(count);
  vvexp(exp_scratch, exp_scratch, &n);
  for (std::size_t i = 0; i < count; ++i) {
    const double x = z_by_lane[i];
    if (!std::isfinite(x)) {
      z_by_lane[i] = x > 0.0 ? 1.0 : 0.0;
      continue;
    }
    const double abs_x = std::abs(x);
    if (abs_x < kInvSqrtTwo) {
      z_by_lane[i] =
          clamp_probability(0.5 + 0.5 * normal_erf_small(x));
      continue;
    }
    double value =
        0.5 * normal_erfc_positive_from_exp(abs_x, exp_scratch[i]);
    if (x > 0.0) {
      value = 1.0 - value;
    }
    z_by_lane[i] = clamp_probability(value);
  }
}

inline void compiled_math_batch_normal_cdf_from_finite_z(
    double *z_by_lane,
    double *exp_scratch,
    const std::size_t count) {
  if (z_by_lane == nullptr || exp_scratch == nullptr || count == 0U) {
    return;
  }
  if (count < 16U) {
    for (std::size_t i = 0; i < count; ++i) {
      const double x = z_by_lane[i] * kInvSqrtTwo;
      const double abs_x = std::abs(x);
      if (abs_x < kInvSqrtTwo) {
        z_by_lane[i] =
            clamp_probability(0.5 + 0.5 * normal_erf_small(x));
        continue;
      }
      double value =
          0.5 * normal_erfc_positive_from_exp(
                    abs_x,
                    std::exp(-abs_x * abs_x));
      if (x > 0.0) {
        value = 1.0 - value;
      }
      z_by_lane[i] = clamp_probability(value);
    }
    return;
  }
  for (std::size_t i = 0; i < count; ++i) {
    const double x = z_by_lane[i] * kInvSqrtTwo;
    z_by_lane[i] = x;
    exp_scratch[i] = -x * x;
  }
  const auto n = static_cast<int>(count);
  vvexp(exp_scratch, exp_scratch, &n);
  for (std::size_t i = 0; i < count; ++i) {
    const double x = z_by_lane[i];
    const double abs_x = std::abs(x);
    if (abs_x < kInvSqrtTwo) {
      z_by_lane[i] =
          clamp_probability(0.5 + 0.5 * normal_erf_small(x));
      continue;
    }
    double value =
        0.5 * normal_erfc_positive_from_exp(abs_x, exp_scratch[i]);
    if (x > 0.0) {
      value = 1.0 - value;
    }
    z_by_lane[i] = clamp_probability(value);
  }
}

inline void compiled_math_batch_normal_log_cdf_from_z(
    double *z_by_lane,
    double *scratch,
    const std::size_t count) {
  if (z_by_lane == nullptr || scratch == nullptr || count == 0U) {
    return;
  }
  compiled_math_batch_normal_cdf_from_z(z_by_lane, scratch, count);
  for (std::size_t i = 0; i < count; ++i) {
    z_by_lane[i] = std::log(std::max(1e-300, z_by_lane[i]));
  }
}
#endif

class CompiledSourceView {
public:
  explicit CompiledSourceView(const ExactVariantPlan &plan)
      : plan_(plan) {}

  void reset(ExactSourceChannels *source_channels,
             const semantic::Index source_view_id = 0) {
    source_channels_ = source_channels;
    source_view_id_ =
        source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  }

  void invalidate_cache() noexcept {}

  ExactSourceChannels *source_channels() const {
    return source_channels_;
  }

  const ExactVariantPlan &plan() const {
    return plan_;
  }

  bool source_order_known_before(
      const semantic::Index before_source_id,
      const semantic::Index after_source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    if (relation_view_knows_before(before_source_id, after_source_id)) {
      return true;
    }
    return source_channels_->source_order_known_before(
        before_source_id, after_source_id, condition_id, workspace);
  }

  const double *exact_time_for_source(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    return source_channels_->exact_time_for_source(
        source_id, condition_id, workspace);
  }

  bool expr_upper_bound_for(
      const semantic::Index expr_id,
      ExactTimedExprUpperBound *out) const {
    return source_channels_ != nullptr &&
           source_channels_->expr_upper_bound_for(expr_id, out);
  }

private:
  bool relation_view_knows_before(const semantic::Index before_source_id,
                                  const semantic::Index after_source_id) const {
    if (before_source_id == semantic::kInvalidIndex ||
        after_source_id == semantic::kInvalidIndex ||
        before_source_id == after_source_id) {
      return false;
    }
    const auto before_relation =
        exact_compiled_source_view_relation(
            plan_, source_view_id_, before_source_id);
    const auto after_relation =
        exact_compiled_source_view_relation(
            plan_, source_view_id_, after_source_id);
    return before_relation != ExactRelation::Unknown &&
           after_relation != ExactRelation::Unknown &&
           before_relation < after_relation;
  }

public:
  semantic::Index condition_cache_id() const noexcept {
    return source_view_id_;
  }

  semantic::Index source_view_id() const noexcept {
    return source_view_id_;
  }

private:
  const ExactVariantPlan &plan_;
  ExactSourceChannels *source_channels_{nullptr};
  semantic::Index source_view_id_{0};
};

inline semantic::Index exact_condition_cache_id(
    const CompiledSourceView *evaluator) {
  return evaluator == nullptr ? 0 : evaluator->condition_cache_id();
}

inline void exact_hash_condition_part(std::size_t *seed,
                                      const std::size_t value) noexcept {
  *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
}

inline semantic::Index compiled_node_condition_cache_id(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace) {
  const auto evaluator_id = exact_condition_cache_id(evaluator);
  if (node.condition_id == 0 ||
      node.condition_id == semantic::kInvalidIndex) {
    return evaluator_id;
  }
  const auto condition_pos = static_cast<std::size_t>(node.condition_id);
  const auto &cache_plan = program.condition_cache_plans[condition_pos];
  if (!cache_plan.dynamic) {
    return compiled_math_static_source_view_condition_cache_id(
        evaluator_id,
        node.condition_id);
  }
  std::size_t seed = 0x9e3779b97f4a7c15ULL;
  exact_hash_condition_part(&seed, static_cast<std::size_t>(evaluator_id));
  exact_hash_condition_part(&seed, static_cast<std::size_t>(node.condition_id));
  for (semantic::Index i = 0; i < cache_plan.time_dependencies.size; ++i) {
    const auto time_id =
        program.condition_cache_time_dependencies[
            static_cast<std::size_t>(
                cache_plan.time_dependencies.offset + i)];
    const auto time_pos = static_cast<std::size_t>(time_id);
    exact_hash_condition_part(&seed, static_cast<std::size_t>(time_id) + 1U);
    exact_hash_condition_part(
        &seed, std::hash<double>{}(workspace.time_values[time_pos]));
  }
  return static_cast<semantic::Index>(
      (seed & static_cast<std::size_t>(0x3fffffffU)) + 1U);
}

struct CompiledEvalWorkspace {
  explicit CompiledEvalWorkspace(const ExactVariantPlan &plan)
      : evaluator(plan),
        compiled_math(plan.compiled_math) {
    const auto source_view_count = plan.compiled_source_views.size();
    source_view_evaluators.reserve(source_view_count);
    for (std::size_t i = 0; i < source_view_count; ++i) {
      source_view_evaluators.emplace_back(plan);
    }
  }

  void reset(ExactSourceChannels *source_channels,
             const semantic::Index source_view_id = 0) {
    evaluator.reset(source_channels, source_view_id);
    if (source_views_channels_ != source_channels) {
      for (std::size_t i = 0; i < source_view_evaluators.size(); ++i) {
        source_view_evaluators[i].reset(
            source_channels, static_cast<semantic::Index>(i + 1U));
      }
      source_views_channels_ = source_channels;
    }
    compiled_math.reset_cache();
  }

  void reset_planned_caches() {}

  CompiledSourceView *source_view_evaluator(
      const semantic::Index source_view_id,
      CompiledSourceView *parent) {
    if (source_view_id == 0 || source_view_id == semantic::kInvalidIndex) {
      return parent;
    }
    if (parent == nullptr) {
      return nullptr;
    }
    const auto pos = static_cast<std::size_t>(source_view_id - 1U);
    if (pos >= parent->plan().compiled_source_views.size() ||
        pos >= source_view_evaluators.size()) {
      throw std::runtime_error("compiled source view id is outside the plan");
    }
    return &source_view_evaluators[pos];
  }

  CompiledSourceView evaluator;
  CompiledMathWorkspace compiled_math;
  std::vector<CompiledSourceView> source_view_evaluators;
  ExactSourceChannels *source_views_channels_{nullptr};
};

inline bool compiled_math_node_cacheable(
    const CompiledMathNodeKind kind) noexcept {
  return kind == CompiledMathNodeKind::SimpleGuardCdf ||
         kind == CompiledMathNodeKind::IntegralZeroToCurrent ||
         kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw ||
         kind == CompiledMathNodeKind::UnionKernelDensity ||
         kind == CompiledMathNodeKind::UnionKernelCdf ||
         kind == CompiledMathNodeKind::UnionKernelMultiSubsetDensity ||
         kind == CompiledMathNodeKind::UnionKernelMultiSubsetCdf;
}

inline bool compiled_math_load_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value);

inline void compiled_math_store_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    double value);

inline double compiled_math_node_time(const CompiledMathNode &node,
                                      const CompiledMathWorkspace &workspace) {
  return workspace.time_values[static_cast<std::size_t>(node.time_id)];
}

inline double compiled_math_source_node_time(
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace) {
  const double time = compiled_math_node_time(node, workspace);
  if (node.aux_id == semantic::kInvalidIndex) {
    return time;
  }
  return std::min(time, workspace.time(node.aux_id));
}

inline bool compiled_math_load_node_cache(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value) {
  if (evaluator == nullptr || !compiled_math_node_cacheable(node.kind)) {
    return false;
  }
  return compiled_math_load_node_cache_entry(
      program, node, evaluator, workspace, value);
}

inline bool compiled_math_load_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value) {
  if (evaluator == nullptr) {
    return false;
  }
  const auto slot = static_cast<std::size_t>(node.cache_slot);
  if (slot >= workspace.cache_valid.size() ||
      workspace.cache_valid[slot] == 0U ||
      workspace.cache_evaluators[slot] != evaluator ||
      workspace.cache_condition_ids[slot] !=
          compiled_node_condition_cache_id(program, node, evaluator, workspace) ||
      workspace.cache_times[slot] !=
          compiled_math_node_time(node, workspace)) {
    return false;
  }
  *value = workspace.cache_values[slot];
  return true;
}

inline CompiledSourceView *compiled_math_node_evaluator(
    const CompiledMathNode &node,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *workspace) {
  if (node.source_view_id == 0 ||
      node.source_view_id == semantic::kInvalidIndex) {
    return parent;
  }
  if (workspace == nullptr) {
    throw std::runtime_error(
        "compiled node requires a planned source view but no workspace was supplied");
  }
  return workspace->source_view_evaluator(node.source_view_id, parent);
}

inline void compiled_math_store_node_cache(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    const double value) {
  if (evaluator == nullptr || !compiled_math_node_cacheable(node.kind)) {
    return;
  }
  compiled_math_store_node_cache_entry(
      program, node, evaluator, workspace, value);
}

inline void compiled_math_store_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    const double value) {
  if (evaluator == nullptr) {
    return;
  }
  const auto slot = static_cast<std::size_t>(node.cache_slot);
  if (slot >= workspace->cache_valid.size()) {
    return;
  }
  workspace->cache_valid[slot] = 1U;
  workspace->cache_evaluators[slot] = evaluator;
  workspace->cache_condition_ids[slot] =
      compiled_node_condition_cache_id(program, node, evaluator, *workspace);
  workspace->cache_times[slot] = compiled_math_node_time(node, *workspace);
  workspace->cache_values[slot] = value;
}

inline double compiled_math_source_product_channel_time(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathWorkspace &workspace) {
  const double time =
      workspace.time_values[static_cast<std::size_t>(channel.time_id)];
  if (channel.time_cap_id == semantic::kInvalidIndex) {
    return time;
  }
  return std::min(
      time,
      workspace.time_values[static_cast<std::size_t>(channel.time_cap_id)]);
}

inline bool compiled_math_source_product_bounds_have_overlay(
    const CompiledSourceBoundPlan &bounds) noexcept {
  return bounds.has_condition_exact ||
         bounds.has_condition_lower ||
         bounds.has_condition_upper;
}

inline bool compiled_math_source_product_relation_forces_fill(
    const ExactRelation relation,
    const std::uint8_t fill_mask) noexcept {
  return relation != ExactRelation::Unknown &&
         !(relation == ExactRelation::At &&
           (fill_mask & kLeafChannelPdf) != 0U);
}

inline ExactRelation compiled_math_source_vector_op_relation(
    const CompiledMathSourceVectorOp &op) noexcept {
  return op.has_static_source_view_relation
             ? static_cast<ExactRelation>(op.static_source_view_relation)
             : ExactRelation::Unknown;
}

inline ExactRelation compiled_math_source_product_channel_relation(
    const CompiledMathSourceProductChannel &channel) noexcept {
  return channel.has_static_source_view_relation
             ? static_cast<ExactRelation>(channel.static_source_view_relation)
             : ExactRelation::Unknown;
}

inline std::uint8_t compiled_math_source_node_fill_mask(
    const CompiledMathNodeKind kind) noexcept {
  const auto value_mask = compiled_math_source_factor_channel_mask(kind);
  if ((value_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    return kLeafChannelCdf | kLeafChannelSurvival;
  }
  return value_mask;
}

struct CompiledMathBatchLane {
  CompiledMathWorkspace *workspace{nullptr};
  CompiledSourceView *parent{nullptr};
  CompiledSourceView *scenario{nullptr};
  CompiledEvalWorkspace *eval_workspace{nullptr};
};

struct CompiledMathLaneBatchState {
  std::size_t lane_count{0};
  std::size_t time_slot_count{0};
  BatchTimeSlotView time_slots{};
  const std::vector<CompiledMathWorkspace *> *workspaces_by_lane{nullptr};
  const std::vector<ExactSourceChannels *> *source_channels_by_lane{nullptr};
  const std::vector<CompiledSourceView *> *parents_by_lane{nullptr};
  const std::vector<CompiledEvalWorkspace *> *eval_workspaces_by_lane{nullptr};
  const std::vector<const std::vector<std::uint8_t> *> *used_outcomes_by_lane{
      nullptr};
  const double *node_values{nullptr};
  std::size_t node_value_lane_stride{0};
  const semantic::Index *node_value_lane_map{nullptr};
  const std::vector<const ExactLoadedLeafInput *> *leaf_inputs_by_lane{nullptr};
  const semantic::Index *lane_parent_map{nullptr};

  semantic::Index storage_lane(const semantic::Index lane) const noexcept {
    return lane_parent_map == nullptr
               ? lane
               : lane_parent_map[static_cast<std::size_t>(lane)];
  }

  CompiledMathWorkspace *workspace(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(storage_lane(lane));
    if (workspaces_by_lane == nullptr || pos >= workspaces_by_lane->size()) {
      return nullptr;
    }
    return (*workspaces_by_lane)[pos];
  }

  ExactSourceChannels *source_channels(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(storage_lane(lane));
    if (source_channels_by_lane == nullptr ||
        pos >= source_channels_by_lane->size()) {
      return nullptr;
    }
    return (*source_channels_by_lane)[pos];
  }

  CompiledSourceView *parent(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(storage_lane(lane));
    if (parents_by_lane == nullptr || pos >= parents_by_lane->size()) {
      return nullptr;
    }
    return (*parents_by_lane)[pos];
  }

  CompiledEvalWorkspace *eval_workspace(const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(storage_lane(lane));
    if (eval_workspaces_by_lane == nullptr ||
        pos >= eval_workspaces_by_lane->size()) {
      return nullptr;
    }
    return (*eval_workspaces_by_lane)[pos];
  }

  const std::vector<std::uint8_t> *used_outcomes(
      const semantic::Index lane) const {
    const auto pos = static_cast<std::size_t>(storage_lane(lane));
    if (used_outcomes_by_lane == nullptr ||
        pos >= used_outcomes_by_lane->size()) {
      return nullptr;
    }
    return (*used_outcomes_by_lane)[pos];
  }

  double time(const semantic::Index time_id,
              const semantic::Index lane) const noexcept {
    return time_slots.get(time_id, storage_lane(lane));
  }

  bool has_time(const semantic::Index time_id,
                const semantic::Index lane) const noexcept {
    return time_slots.has(time_id, storage_lane(lane));
  }

  double node_value(const semantic::Index node_id,
                    const semantic::Index lane) const {
    const auto node_pos = static_cast<std::size_t>(node_id);
    if (node_values != nullptr && node_value_lane_stride != 0U) {
      const auto mapped_lane =
          node_value_lane_map == nullptr
              ? storage_lane(lane)
              : node_value_lane_map[
                    static_cast<std::size_t>(storage_lane(lane))];
      return node_values[
          node_pos * node_value_lane_stride +
          static_cast<std::size_t>(mapped_lane)];
    }
    (void)node_pos;
    return 0.0;
  }
};

inline const ExactLoadedLeafInput *compiled_math_batch_leaf_input_for_lane(
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const semantic::Index leaf_index) {
  if (state.leaf_inputs_by_lane != nullptr &&
      static_cast<std::size_t>(state.storage_lane(lane)) <
          state.leaf_inputs_by_lane->size()) {
    return (*state.leaf_inputs_by_lane)[
        static_cast<std::size_t>(state.storage_lane(lane))];
  }
  auto *source_channels = state.source_channels(lane);
  if (source_channels == nullptr) {
    return nullptr;
  }
  return &source_channels->source_product_leaf_input(leaf_index);
}

struct CompiledMathBatchSourceProductScratch {
  std::vector<semantic::Index> lanes_a;
  std::vector<semantic::Index> lanes_b;
  std::vector<semantic::Index> lanes_c;
  std::vector<semantic::Index> lanes_d;
  std::vector<double> times;
  std::vector<double> onset_upper;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> exact;
  std::vector<std::uint8_t> has_exact;
  std::vector<std::uint8_t> mask_a;
  std::vector<double> pdf_a;
  std::vector<double> cdf_a;
  std::vector<double> survival_a;
  std::vector<std::uint8_t> mask_b;
  std::vector<double> pdf_b;
  std::vector<double> cdf_b;
  std::vector<double> survival_b;
  std::vector<std::uint8_t> mask_c;
  std::vector<double> pdf_c;
  std::vector<double> cdf_c;
  std::vector<double> survival_c;
  std::vector<std::uint8_t> source_vector_mask_slots;
  std::vector<double> source_vector_pdf_slots;
  std::vector<double> source_vector_cdf_slots;
  std::vector<double> source_vector_survival_slots;
  std::vector<double> pool_member_pdf;
  std::vector<double> pool_member_cdf;
  std::vector<double> pool_member_survival;
  std::vector<double> pool_prefix;
  std::vector<double> pool_suffix;
  std::vector<const double *> pool_member_pdf_ptrs;
  std::vector<const double *> pool_member_cdf_ptrs;
  std::vector<const double *> pool_member_survival_ptrs;
  std::vector<semantic::Index> expanded_lanes;
  std::vector<semantic::Index> expanded_parent_lanes;
  std::vector<semantic::Index> expanded_node_value_lane_map;
  std::vector<double> expanded_weights;
  std::vector<double> expanded_source_times;
  std::vector<double> expanded_shifted_times;
  std::vector<double> expanded_time_values;
  std::vector<std::uint8_t> expanded_time_valid;
  std::vector<CompiledMathWorkspace *> expanded_workspaces_by_lane;
  std::vector<ExactSourceChannels *> expanded_source_channels_by_lane;
  std::vector<CompiledSourceView *> expanded_parents_by_lane;
  std::vector<CompiledEvalWorkspace *> expanded_eval_workspaces_by_lane;
  std::vector<const std::vector<std::uint8_t> *>
      expanded_used_outcomes_by_lane;

  void ensure_size(const std::size_t lane_capacity,
                   const std::size_t active_capacity) {
    if (lanes_a.size() < active_capacity) {
      lanes_a.resize(active_capacity);
    }
    if (lanes_b.size() < active_capacity) {
      lanes_b.resize(active_capacity);
    }
    if (lanes_c.size() < active_capacity) {
      lanes_c.resize(active_capacity);
    }
    if (lanes_d.size() < active_capacity) {
      lanes_d.resize(active_capacity);
    }
    if (times.size() < lane_capacity) {
      times.resize(lane_capacity, 0.0);
    }
    if (onset_upper.size() < lane_capacity) {
      onset_upper.resize(lane_capacity, 0.0);
    }
    if (lower.size() < lane_capacity) {
      lower.resize(lane_capacity, 0.0);
    }
    if (upper.size() < lane_capacity) {
      upper.resize(lane_capacity, 0.0);
    }
    if (exact.size() < lane_capacity) {
      exact.resize(lane_capacity, 0.0);
    }
    if (has_exact.size() < lane_capacity) {
      has_exact.resize(lane_capacity, 0U);
    }
    if (mask_a.size() < lane_capacity) {
      mask_a.resize(lane_capacity, 0U);
      pdf_a.resize(lane_capacity, 0.0);
      cdf_a.resize(lane_capacity, 0.0);
      survival_a.resize(lane_capacity, 1.0);
    }
    if (mask_b.size() < lane_capacity) {
      mask_b.resize(lane_capacity, 0U);
      pdf_b.resize(lane_capacity, 0.0);
      cdf_b.resize(lane_capacity, 0.0);
      survival_b.resize(lane_capacity, 1.0);
    }
    if (mask_c.size() < lane_capacity) {
      mask_c.resize(lane_capacity, 0U);
      pdf_c.resize(lane_capacity, 0.0);
      cdf_c.resize(lane_capacity, 0.0);
      survival_c.resize(lane_capacity, 1.0);
    }
  }

  void ensure_pool_size(const std::size_t lane_capacity,
                        const std::size_t member_count) {
    const auto member_value_count = lane_capacity * member_count;
    if (pool_member_pdf.size() < member_value_count) {
      pool_member_pdf.resize(member_value_count, 0.0);
      pool_member_cdf.resize(member_value_count, 0.0);
      pool_member_survival.resize(member_value_count, 1.0);
    }
    const auto width = member_count + 1U;
    const auto table_value_count = lane_capacity * width * width;
    if (pool_prefix.size() < table_value_count) {
      pool_prefix.resize(table_value_count, 0.0);
      pool_suffix.resize(table_value_count, 0.0);
    }
    if (pool_member_pdf_ptrs.size() < member_count) {
      pool_member_pdf_ptrs.resize(member_count, nullptr);
      pool_member_cdf_ptrs.resize(member_count, nullptr);
      pool_member_survival_ptrs.resize(member_count, nullptr);
    }
  }

  void ensure_source_vector_slots(const std::size_t lane_capacity,
                                  const std::size_t slot_count) {
    const auto value_count = lane_capacity * slot_count;
    if (source_vector_mask_slots.size() < value_count) {
      source_vector_mask_slots.resize(value_count, 0U);
      source_vector_pdf_slots.resize(value_count, 0.0);
      source_vector_cdf_slots.resize(value_count, 0.0);
      source_vector_survival_slots.resize(value_count, 1.0);
    }
  }

  void ensure_expanded_lane_state(const std::size_t time_slot_count,
                                  const std::size_t lane_capacity) {
    if (expanded_lanes.size() < lane_capacity) {
      expanded_lanes.resize(lane_capacity);
      expanded_parent_lanes.resize(lane_capacity);
      expanded_node_value_lane_map.resize(lane_capacity);
      expanded_weights.resize(lane_capacity, 0.0);
      expanded_source_times.resize(lane_capacity, 0.0);
      expanded_shifted_times.resize(lane_capacity, 0.0);
      expanded_workspaces_by_lane.resize(lane_capacity, nullptr);
      expanded_source_channels_by_lane.resize(lane_capacity, nullptr);
      expanded_parents_by_lane.resize(lane_capacity, nullptr);
      expanded_eval_workspaces_by_lane.resize(lane_capacity, nullptr);
      expanded_used_outcomes_by_lane.resize(lane_capacity, nullptr);
    }
    const auto time_value_count = time_slot_count * lane_capacity;
    if (expanded_time_values.size() < time_value_count) {
      expanded_time_values.resize(time_value_count, 0.0);
      expanded_time_valid.resize(time_value_count, 0U);
    }
  }
};

struct CompiledMathDirectLeafPreparedOp {
  const CompiledMathSourceProductOp *op{nullptr};
  const CompiledMathSourceProductChannel *channel{nullptr};
  semantic::Index tile_time_plan_id{semantic::kInvalidIndex};
  bool time_is_bind{false};
  bool has_time_cap{false};
  bool cap_is_bind{false};
  bool conditioned_direct_leaf{false};
};

struct CompiledMathBatchIntegralFrame {
  std::vector<double> time_values;
  std::vector<std::uint8_t> time_valid;
  std::vector<CompiledMathWorkspace *> workspaces_by_lane;
  std::vector<ExactSourceChannels *> source_channels_by_lane;
  std::vector<CompiledSourceView *> parents_by_lane;
  std::vector<CompiledEvalWorkspace *> eval_workspaces_by_lane;
  std::vector<const std::vector<std::uint8_t> *> used_outcomes_by_lane;
  std::vector<semantic::Index> node_value_lane_map;
  std::vector<semantic::Index> active_lanes;
  std::vector<semantic::Index> parent_lanes;
  std::vector<semantic::Index> term_lanes_a;
  std::vector<semantic::Index> term_lanes_b;
  std::vector<semantic::Index> source_lanes_a;
  std::vector<semantic::Index> source_lanes_b;
  std::vector<double> weights;
  std::vector<double> values;
  std::vector<double> products;
  std::vector<double> factor_values;
  std::vector<double> source_values;
  std::vector<double> source_leaf_values;
  std::vector<double> bind_times;
  std::vector<const ExactLoadedLeafInput *> source_leaf_inputs_by_lane;
  std::vector<const ExactLoadedLeafInput *> source_leaf_inputs_by_op_parent;
  std::vector<ExactSourceChannels *> direct_source_channels_by_parent;
  std::vector<CompiledMathDirectLeafPreparedOp> direct_leaf_ops;
  std::vector<double> direct_parent_scale;
  std::vector<double> direct_parent_shift;
  std::vector<double> direct_source_time_by_op_parent;
  std::vector<double> direct_source_cap_by_op_parent;
  std::vector<double> direct_bound_lower_by_op_parent;
  std::vector<double> direct_bound_upper_by_op_parent;
  std::vector<double> direct_bound_exact_by_op_parent;
  std::vector<std::uint8_t> direct_bound_has_exact_by_op_parent;
  std::vector<double> direct_lower_cdf_by_op_parent;
  std::vector<double> direct_lower_survival_by_op_parent;
  std::vector<double> direct_upper_cdf_by_op_parent;
  std::vector<double> cached_integral_factor_values;
  std::vector<double> source_product_block_factor_values;
  std::vector<std::uint8_t> source_product_block_factor_valid;
  std::vector<CompiledSourceView *> condition_evaluators;
  std::vector<std::uint8_t> upper_found_by_lane;
  std::vector<double> upper_time_by_lane;
  std::vector<double> upper_normalizer_by_lane;
  std::vector<double> active_lower_bounds;
  std::vector<double> active_upper_bounds;
  std::vector<double> lower_by_lane;
  std::vector<double> upper_by_lane;
  std::vector<std::uint8_t> source_vector_valid_mask;
  std::vector<double> source_vector_pdf;
  std::vector<double> source_vector_cdf;
  std::vector<double> source_vector_survival;
  std::vector<semantic::Index> source_vector_cache_slot_by_program;
  std::vector<semantic::Index> cached_source_vector_programs;
  std::vector<CompiledMathBatchSourceProductScratch>
      source_product_batch_scratch;

  void ensure_lane_capacity(const std::size_t time_slot_count,
                            const std::size_t lane_capacity) {
    const auto time_value_count = time_slot_count * lane_capacity;
    if (time_values.size() < time_value_count) {
      time_values.resize(time_value_count, 0.0);
    }
    if (time_valid.size() < time_value_count) {
      time_valid.resize(time_value_count, 0U);
    }
    if (workspaces_by_lane.size() < lane_capacity) {
      workspaces_by_lane.resize(lane_capacity, nullptr);
    }
    if (source_channels_by_lane.size() < lane_capacity) {
      source_channels_by_lane.resize(lane_capacity, nullptr);
    }
    if (parents_by_lane.size() < lane_capacity) {
      parents_by_lane.resize(lane_capacity, nullptr);
    }
    if (eval_workspaces_by_lane.size() < lane_capacity) {
      eval_workspaces_by_lane.resize(lane_capacity, nullptr);
    }
    if (used_outcomes_by_lane.size() < lane_capacity) {
      used_outcomes_by_lane.resize(lane_capacity, nullptr);
    }
    if (node_value_lane_map.size() < lane_capacity) {
      node_value_lane_map.resize(lane_capacity);
    }
    if (active_lanes.size() < lane_capacity) {
      active_lanes.resize(lane_capacity);
    }
    if (parent_lanes.size() < lane_capacity) {
      parent_lanes.resize(lane_capacity);
    }
    if (term_lanes_a.size() < lane_capacity) {
      term_lanes_a.resize(lane_capacity);
    }
    if (term_lanes_b.size() < lane_capacity) {
      term_lanes_b.resize(lane_capacity);
    }
    if (source_lanes_a.size() < lane_capacity) {
      source_lanes_a.resize(lane_capacity);
    }
    if (source_lanes_b.size() < lane_capacity) {
      source_lanes_b.resize(lane_capacity);
    }
    if (weights.size() < lane_capacity) {
      weights.resize(lane_capacity, 0.0);
    }
    if (values.size() < lane_capacity) {
      values.resize(lane_capacity, 0.0);
    }
    if (products.size() < lane_capacity) {
      products.resize(lane_capacity, 0.0);
    }
    if (factor_values.size() < lane_capacity) {
      factor_values.resize(lane_capacity, 0.0);
    }
    if (source_values.size() < lane_capacity) {
      source_values.resize(lane_capacity, 0.0);
    }
    if (source_leaf_values.size() < lane_capacity) {
      source_leaf_values.resize(lane_capacity, 0.0);
    }
    if (bind_times.size() < lane_capacity) {
      bind_times.resize(lane_capacity, 0.0);
    }
    if (source_leaf_inputs_by_lane.size() < lane_capacity) {
      source_leaf_inputs_by_lane.resize(lane_capacity, nullptr);
    }
    if (direct_source_channels_by_parent.size() < lane_capacity) {
      direct_source_channels_by_parent.resize(lane_capacity, nullptr);
    }
    if (direct_parent_scale.size() < lane_capacity) {
      direct_parent_scale.resize(lane_capacity, 0.0);
    }
    if (direct_parent_shift.size() < lane_capacity) {
      direct_parent_shift.resize(lane_capacity, 0.0);
    }
    if (condition_evaluators.size() < lane_capacity) {
      condition_evaluators.resize(lane_capacity, nullptr);
    }
    if (upper_found_by_lane.size() < lane_capacity) {
      upper_found_by_lane.resize(lane_capacity, 0U);
    }
    if (upper_time_by_lane.size() < lane_capacity) {
      upper_time_by_lane.resize(
          lane_capacity,
          std::numeric_limits<double>::infinity());
    }
    if (upper_normalizer_by_lane.size() < lane_capacity) {
      upper_normalizer_by_lane.resize(lane_capacity, 0.0);
    }
    if (active_lower_bounds.size() < lane_capacity) {
      active_lower_bounds.resize(lane_capacity, 0.0);
    }
    if (active_upper_bounds.size() < lane_capacity) {
      active_upper_bounds.resize(lane_capacity, 0.0);
    }
    if (lower_by_lane.size() < lane_capacity) {
      lower_by_lane.resize(lane_capacity, 0.0);
    }
    if (upper_by_lane.size() < lane_capacity) {
      upper_by_lane.resize(lane_capacity, 0.0);
    }
  }

  void reset_source_vector_cache(
      const std::size_t program_count,
      const std::size_t lane_capacity,
      const std::vector<semantic::Index> &program_ids) {
    source_vector_cache_slot_by_program.assign(
        program_count,
        semantic::kInvalidIndex);
    std::size_t cache_slot_count = 0U;
    for (const auto program_id : program_ids) {
      if (program_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(program_id) >= program_count) {
        continue;
      }
      auto &slot =
          source_vector_cache_slot_by_program[
              static_cast<std::size_t>(program_id)];
      if (slot != semantic::kInvalidIndex) {
        continue;
      }
      slot = static_cast<semantic::Index>(cache_slot_count++);
    }
    const auto count = cache_slot_count * lane_capacity;
    if (source_vector_valid_mask.size() < count) {
      source_vector_valid_mask.resize(count, 0U);
      source_vector_pdf.resize(count, 0.0);
      source_vector_cdf.resize(count, 0.0);
      source_vector_survival.resize(count, 1.0);
    }
    std::fill(
        source_vector_valid_mask.begin(),
        source_vector_valid_mask.begin() + static_cast<std::ptrdiff_t>(count),
        0U);
  }

  semantic::Index source_vector_cache_slot(
      const semantic::Index program_id) const noexcept {
    if (program_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(program_id) >=
            source_vector_cache_slot_by_program.size()) {
      return semantic::kInvalidIndex;
    }
    return source_vector_cache_slot_by_program[
        static_cast<std::size_t>(program_id)];
  }

  CompiledMathBatchSourceProductScratch &source_product_scratch_layer(
      const std::size_t depth,
      const std::size_t lane_capacity,
      const std::size_t active_capacity) {
    constexpr std::size_t kMaxSourceProductScratchDepth = 32U;
    if (source_product_batch_scratch.size() < kMaxSourceProductScratchDepth) {
      source_product_batch_scratch.resize(kMaxSourceProductScratchDepth);
    } else if (source_product_batch_scratch.size() <= depth) {
      source_product_batch_scratch.resize(depth + 1U);
    }
    auto &scratch = source_product_batch_scratch[depth];
    scratch.ensure_size(lane_capacity, active_capacity);
    return scratch;
  }
};

struct CompiledMathBatchWorkspace {
  std::vector<CompiledSourceView *> condition_evaluators;
  std::vector<semantic::Index> active_lanes;
  std::vector<std::size_t> active_lane_indices;
  std::vector<std::size_t> compact_lane_indices;
  std::vector<CompiledMathWorkspace *> workspaces_by_lane;
  std::vector<ExactSourceChannels *> source_channels_by_lane;
  std::vector<CompiledSourceView *> parents_by_lane;
  std::vector<CompiledEvalWorkspace *> eval_workspaces_by_lane;
  std::vector<const std::vector<std::uint8_t> *> used_outcomes_by_lane;
  std::vector<semantic::Index> node_value_lane_map;
  std::vector<double> parent_time_values;
  std::vector<std::uint8_t> parent_time_valid;
  std::vector<ExactLoadedLeafInput> parent_leaf_inputs;
  std::vector<double> lower_by_lane;
  std::vector<double> upper_by_lane;
  std::vector<std::uint8_t> upper_found_by_lane;
  std::vector<double> upper_time_by_lane;
  std::vector<double> upper_normalizer_by_lane;
  std::vector<double> batch_values;
  std::vector<double> node_values;
  std::size_t node_value_lane_stride{0};
  std::vector<CompiledMathBatchIntegralFrame> integral_frames;
};

inline void compiled_math_batch_prepare_node_values(
    const CompiledMathProgram &program,
    const std::size_t lane_count,
    CompiledMathBatchWorkspace *workspace) {
  if (workspace == nullptr) {
    return;
  }
  workspace->node_value_lane_stride = lane_count;
  const auto value_count = program.nodes.size() * lane_count;
  if (workspace->node_values.size() < value_count) {
    workspace->node_values.resize(value_count, 0.0);
  } else {
    std::fill(
        workspace->node_values.begin(),
        workspace->node_values.begin() +
            static_cast<std::ptrdiff_t>(value_count),
        0.0);
  }
}

inline void compiled_math_batch_set_node_value(
    const CompiledMathNode &node,
    const semantic::Index lane,
    const double value,
    CompiledMathBatchWorkspace *workspace) {
  if (workspace == nullptr || workspace->node_value_lane_stride == 0U ||
      node.cache_slot == semantic::kInvalidIndex) {
    return;
  }
  const auto node_pos = static_cast<std::size_t>(node.cache_slot);
  const auto lane_pos = static_cast<std::size_t>(lane);
  workspace->node_values[
      node_pos * workspace->node_value_lane_stride + lane_pos] = value;
}

inline double compiled_math_batch_node_value(
    const semantic::Index node_id,
    const semantic::Index lane,
    const CompiledMathBatchWorkspace *workspace) {
  if (workspace == nullptr || workspace->node_value_lane_stride == 0U ||
      node_id == semantic::kInvalidIndex) {
    return 0.0;
  }
  const auto node_pos = static_cast<std::size_t>(node_id);
  const auto lane_pos = static_cast<std::size_t>(lane);
  return workspace->node_values[
      node_pos * workspace->node_value_lane_stride + lane_pos];
}

inline double *compiled_math_batch_node_value_array(
    const CompiledMathNode &node,
    CompiledMathBatchWorkspace *workspace) {
  if (workspace == nullptr || workspace->node_value_lane_stride == 0U ||
      node.cache_slot == semantic::kInvalidIndex) {
    return nullptr;
  }
  return workspace->node_values.data() +
         static_cast<std::size_t>(node.cache_slot) *
             workspace->node_value_lane_stride;
}

inline bool compiled_math_batch_product_is_live(
    const double value) noexcept {
  return std::isfinite(value) && value != 0.0;
}

inline bool compiled_math_batch_multiply_product_value_at(
    const std::size_t lane_pos,
    const double factor,
    double *products_by_lane) noexcept {
  const double current = products_by_lane[lane_pos];
  if (!compiled_math_batch_product_is_live(current)) {
    products_by_lane[lane_pos] = 0.0;
    return false;
  }
  const double product = current * factor;
  if (compiled_math_batch_product_is_live(product)) {
    products_by_lane[lane_pos] = product;
    return true;
  }
  products_by_lane[lane_pos] = 0.0;
  return false;
}

inline bool compiled_math_batch_lanes_are_contiguous(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::size_t *offset) {
  if (lane_count == 0U) {
    if (offset != nullptr) {
      *offset = 0U;
    }
    return true;
  }
  const auto first = static_cast<std::size_t>(lanes[0]);
  for (std::size_t i = 1; i < lane_count; ++i) {
    if (static_cast<std::size_t>(lanes[i]) != first + i) {
      return false;
    }
  }
  if (offset != nullptr) {
    *offset = first;
  }
  return true;
}

inline void compiled_math_batch_array_fill(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double value,
    double *out_by_lane) {
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    std::fill(out_by_lane + offset, out_by_lane + offset + lane_count, value);
    return;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    out_by_lane[static_cast<std::size_t>(lanes[i])] = value;
  }
}

inline std::size_t compiled_math_batch_array_compact_live(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *values_by_lane,
    semantic::Index *next_lanes) {
  std::size_t next_count = 0U;
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      if (values_by_lane[offset + i] != 0.0) {
        next_lanes[next_count++] =
            static_cast<semantic::Index>(offset + i);
      }
    }
    return next_count;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    if (values_by_lane[static_cast<std::size_t>(lane)] != 0.0) {
      next_lanes[next_count++] = lane;
    }
  }
  return next_count;
}

inline std::size_t compiled_math_batch_array_live_lanes_from_mask(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *products_by_lane,
    semantic::Index *live_lanes) {
  std::size_t live_count = 0U;
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      if (compiled_math_batch_product_is_live(products_by_lane[offset + i])) {
        live_lanes[live_count++] =
            static_cast<semantic::Index>(offset + i);
      }
    }
    return live_count;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    if (compiled_math_batch_product_is_live(
            products_by_lane[static_cast<std::size_t>(lane)])) {
      live_lanes[live_count++] = lane;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_live_lanes_from_dense_mask(
    const std::size_t lane_count,
    const double *products_by_lane,
    semantic::Index *live_lanes) {
  std::size_t live_count = 0U;
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    if (compiled_math_batch_product_is_live(products_by_lane[lane_pos])) {
      live_lanes[live_count++] = static_cast<semantic::Index>(lane_pos);
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_multiply_scalar_mask(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double scalar,
    double *products_by_lane) {
  if (scalar == 0.0 || !std::isfinite(scalar)) {
    compiled_math_batch_array_fill(lanes, lane_count, 0.0, products_by_lane);
    return 0U;
  }
  std::size_t live_count = 0U;
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    double *products = products_by_lane + offset;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const double current = products[i];
      if (!compiled_math_batch_product_is_live(current)) {
        products[i] = 0.0;
        continue;
      }
      const double product = current * scalar;
      if (compiled_math_batch_product_is_live(product)) {
        products[i] = product;
        ++live_count;
      } else {
        products[i] = 0.0;
      }
    }
    return live_count;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[i]);
    const double current = products_by_lane[lane_pos];
    if (!compiled_math_batch_product_is_live(current)) {
      products_by_lane[lane_pos] = 0.0;
      continue;
    }
    const double product = current * scalar;
    if (compiled_math_batch_product_is_live(product)) {
      products_by_lane[lane_pos] = product;
      ++live_count;
    } else {
      products_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_multiply_scalar_dense_compact(
    const std::size_t lane_count,
    const double scalar,
    double *products_by_lane,
    semantic::Index *next_lanes) {
  if (scalar == 0.0 || !std::isfinite(scalar)) {
    std::fill(products_by_lane, products_by_lane + lane_count, 0.0);
    return 0U;
  }
  std::size_t live_count = 0U;
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const double current = products_by_lane[lane_pos];
    if (!compiled_math_batch_product_is_live(current)) {
      products_by_lane[lane_pos] = 0.0;
      continue;
    }
    const double product = current * scalar;
    if (compiled_math_batch_product_is_live(product)) {
      products_by_lane[lane_pos] = product;
      next_lanes[live_count++] = static_cast<semantic::Index>(lane_pos);
    } else {
      products_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_multiply_values_dense_compact(
    const std::size_t lane_count,
    const double *factors_by_lane,
    double *products_by_lane,
    semantic::Index *next_lanes) {
  std::size_t live_count = 0U;
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const double current = products_by_lane[lane_pos];
    if (!compiled_math_batch_product_is_live(current)) {
      products_by_lane[lane_pos] = 0.0;
      continue;
    }
    const double product = current * factors_by_lane[lane_pos];
    if (compiled_math_batch_product_is_live(product)) {
      products_by_lane[lane_pos] = product;
      next_lanes[live_count++] = static_cast<semantic::Index>(lane_pos);
    } else {
      products_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_multiply_values_mask(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *factors_by_lane,
    double *products_by_lane) {
  std::size_t live_count = 0U;
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    const double *factors = factors_by_lane + offset;
    double *products = products_by_lane + offset;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const double current = products[i];
      if (!compiled_math_batch_product_is_live(current)) {
        products[i] = 0.0;
        continue;
      }
      const double product = current * factors[i];
      if (compiled_math_batch_product_is_live(product)) {
        products[i] = product;
        ++live_count;
      } else {
        products[i] = 0.0;
      }
    }
    return live_count;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[i]);
    const double current = products_by_lane[lane_pos];
    if (!compiled_math_batch_product_is_live(current)) {
      products_by_lane[lane_pos] = 0.0;
      continue;
    }
    const double product = current * factors_by_lane[lane_pos];
    if (compiled_math_batch_product_is_live(product)) {
      products_by_lane[lane_pos] = product;
      ++live_count;
    } else {
      products_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_store_live_values_dense_mask(
    const std::size_t lane_count,
    const double *values_by_lane,
    double *out_by_lane) {
  std::size_t live_count = 0U;
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const double value = values_by_lane[lane_pos];
    if (compiled_math_batch_product_is_live(value)) {
      out_by_lane[lane_pos] = value;
      ++live_count;
    } else {
      out_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_store_live_values_mask(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *values_by_lane,
    double *out_by_lane) {
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    return compiled_math_batch_array_store_live_values_dense_mask(
        lane_count,
        values_by_lane + offset,
        out_by_lane + offset);
  }
  std::size_t live_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[i]);
    const double value = values_by_lane[lane_pos];
    if (compiled_math_batch_product_is_live(value)) {
      out_by_lane[lane_pos] = value;
      ++live_count;
    } else {
      out_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_multiply_scalar_compact(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double scalar,
    double *products_by_lane,
    semantic::Index *next_lanes) {
  if (scalar == 0.0 || !std::isfinite(scalar)) {
    compiled_math_batch_array_fill(lanes, lane_count, 0.0, products_by_lane);
    return 0U;
  }
  std::size_t live_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double product = products_by_lane[lane_pos] * scalar;
    if (compiled_math_batch_product_is_live(product)) {
      products_by_lane[lane_pos] = product;
      next_lanes[live_count++] = lane;
    } else {
      products_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_multiply_values_compact(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *factors_by_lane,
    double *products_by_lane,
    semantic::Index *next_lanes) {
  std::size_t live_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double product =
        products_by_lane[lane_pos] * factors_by_lane[lane_pos];
    if (compiled_math_batch_product_is_live(product)) {
      products_by_lane[lane_pos] = product;
      next_lanes[live_count++] = lane;
    } else {
      products_by_lane[lane_pos] = 0.0;
    }
  }
  return live_count;
}

inline std::size_t compiled_math_batch_array_apply_expr_upper_op_compact(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const CompiledMathLaneTermOp &op,
    const CompiledMathLaneBatchState &state,
    const std::uint8_t *upper_found_by_lane,
    const double *upper_time_by_lane,
    const double *upper_normalizer_by_lane,
    double *products_by_lane,
    semantic::Index *next_lanes) {
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const auto lane = lanes[lane_pos];
    const auto lane_index = static_cast<std::size_t>(lane);
    const bool has_upper =
        upper_found_by_lane[lane_index] != 0U &&
        upper_normalizer_by_lane[lane_index] > 0.0;
    const double current_time = state.time(op.expr_upper_time_id, lane);
    if (op.expr_upper_mode == CompiledMathIntegralExprUpperMode::AfterOne) {
      if (!has_upper || current_time < upper_time_by_lane[lane_index]) {
        products_by_lane[lane_index] = 0.0;
      }
      continue;
    }
    if (!has_upper) {
      continue;
    }
    if (current_time >= upper_time_by_lane[lane_index]) {
      products_by_lane[lane_index] = 0.0;
      continue;
    }
    products_by_lane[lane_index] /= upper_normalizer_by_lane[lane_index];
    if (!compiled_math_batch_product_is_live(products_by_lane[lane_index])) {
      products_by_lane[lane_index] = 0.0;
    }
  }
  return compiled_math_batch_array_compact_live(
      lanes,
      lane_count,
      products_by_lane,
      next_lanes);
}

inline void compiled_math_batch_array_weighted_accumulate(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *weights_by_lane,
    const double *values_by_lane,
    double *out_by_lane) {
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    const double *weights = weights_by_lane + offset;
    const double *values = values_by_lane + offset;
    double *out = out_by_lane + offset;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const double value = values[i];
      if (compiled_math_batch_product_is_live(value)) {
        out[i] += weights[i] * value;
      }
    }
    return;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[i]);
    const double value = values_by_lane[lane_pos];
    if (compiled_math_batch_product_is_live(value)) {
      out_by_lane[lane_pos] += weights_by_lane[lane_pos] * value;
    }
  }
}

inline void compiled_math_batch_array_accumulate_product(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *left_by_lane,
    const double *right_by_lane,
    double *out_by_lane) {
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          lanes, lane_count, &offset)) {
    const double *left = left_by_lane + offset;
    const double *right = right_by_lane + offset;
    double *out = out_by_lane + offset;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const double product = left[i] * right[i];
      if (compiled_math_batch_product_is_live(product)) {
        out[i] += product;
      }
    }
    return;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[i]);
    const double product = left_by_lane[lane_pos] * right_by_lane[lane_pos];
    if (compiled_math_batch_product_is_live(product)) {
      out_by_lane[lane_pos] += product;
    }
  }
}

inline void compiled_math_batch_array_weighted_accumulate_parent(
    const semantic::Index *child_lanes,
    const semantic::Index *parent_lanes,
    const std::size_t lane_count,
    const double *weights_by_child_lane,
    const double *values_by_child_lane,
    double *out_by_parent_lane) {
  std::size_t offset = 0U;
  if (compiled_math_batch_lanes_are_contiguous(
          child_lanes,
          lane_count,
          &offset)) {
    const double *weights = weights_by_child_lane + offset;
    const double *values = values_by_child_lane + offset;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const double value = values[i];
      if (compiled_math_batch_product_is_live(value)) {
        out_by_parent_lane[static_cast<std::size_t>(parent_lanes[i])] +=
            weights[i] * value;
      }
    }
    return;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto child_lane = child_lanes[i];
    const auto child_pos = static_cast<std::size_t>(child_lane);
    const double value = values_by_child_lane[child_pos];
    if (compiled_math_batch_product_is_live(value)) {
      out_by_parent_lane[static_cast<std::size_t>(parent_lanes[i])] +=
          weights_by_child_lane[i] * value;
    }
  }
}

inline void compiled_math_batch_array_weighted_accumulate_parent_rows(
    const semantic::Index *parent_lanes,
    const std::size_t parent_count,
    const std::size_t child_width,
    const double *weights_by_child_lane,
    const double *values_by_child_lane,
    double *out_by_parent_lane) {
  for (std::size_t parent_pos = 0; parent_pos < parent_count; ++parent_pos) {
    const auto parent_lane = static_cast<std::size_t>(parent_lanes[parent_pos]);
    const auto row_offset = parent_pos * child_width;
    double sum = 0.0;
    for (std::size_t child_pos = 0; child_pos < child_width; ++child_pos) {
      const auto pos = row_offset + child_pos;
      const double value = values_by_child_lane[pos];
      if (compiled_math_batch_product_is_live(value)) {
        sum += weights_by_child_lane[pos] * value;
      }
    }
    out_by_parent_lane[parent_lane] += sum;
  }
}

inline std::size_t compiled_math_batch_active_workspace_lanes(
    const std::vector<CompiledMathBatchLane> &lanes,
    std::vector<semantic::Index> *active_lanes) {
  active_lanes->clear();
  active_lanes->reserve(lanes.size());
  for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
    if (lanes[lane_pos].workspace != nullptr) {
      active_lanes->push_back(static_cast<semantic::Index>(lane_pos));
    }
  }
  return active_lanes->size();
}

inline void compiled_math_batch_regular_constant(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double value,
    double *out_by_lane) {
  compiled_math_batch_array_fill(lanes, lane_count, value, out_by_lane);
}

inline void compiled_math_batch_regular_time_gate(
    const CompiledMathNode &node,
    const semantic::Index child_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *out_by_lane) {
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double raw = state.node_value(child_id, lane);
    out_by_lane[lane_pos] =
        state.time(node.time_id, lane) >= state.time(node.aux_id, lane)
            ? raw
            : 0.0;
  }
}

inline void compiled_math_batch_regular_product(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *out_by_lane) {
  for (std::size_t lane_pos_idx = 0; lane_pos_idx < lane_count;
       ++lane_pos_idx) {
    const auto lane = lanes[lane_pos_idx];
    const auto lane_pos = static_cast<std::size_t>(lane);
    double value = 1.0;
    for (semantic::Index i = 0; i < node.children.size; ++i) {
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset + i)];
      value *= state.node_value(child_id, lane);
      if (!compiled_math_batch_product_is_live(value)) {
        value = 0.0;
        break;
      }
    }
    out_by_lane[lane_pos] = value;
  }
}

inline void compiled_math_batch_regular_sum(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const bool clean_signed,
    double *out_by_lane) {
  for (std::size_t lane_pos_idx = 0; lane_pos_idx < lane_count;
       ++lane_pos_idx) {
    const auto lane = lanes[lane_pos_idx];
    double total = 0.0;
    for (semantic::Index i = 0; i < node.children.size; ++i) {
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset + i)];
      total += state.node_value(child_id, lane);
    }
    out_by_lane[static_cast<std::size_t>(lane)] =
        clean_signed ? clean_signed_value(total) : total;
  }
}

inline void compiled_math_batch_regular_unary(
    const CompiledMathNodeKind kind,
    const semantic::Index child_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *out_by_lane) {
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const double raw = state.node_value(child_id, lane);
    double value = raw;
    if (kind == CompiledMathNodeKind::ClampProbability) {
      value = clamp_probability(raw);
    } else if (kind == CompiledMathNodeKind::Complement) {
      value = clamp_probability(1.0 - raw);
    } else if (kind == CompiledMathNodeKind::Negate) {
      value = clean_signed_value(-raw);
    }
    out_by_lane[static_cast<std::size_t>(lane)] = value;
  }
}

inline CompiledMathBatchIntegralFrame &compiled_math_batch_integral_frame(
    CompiledMathBatchWorkspace *workspace,
    const std::size_t depth) {
  constexpr std::size_t kMaxBatchIntegralDepth = 18U;
  if (workspace->integral_frames.size() < kMaxBatchIntegralDepth) {
    workspace->integral_frames.resize(kMaxBatchIntegralDepth);
  } else if (workspace->integral_frames.size() <= depth) {
    workspace->integral_frames.resize(depth + 1U);
  }
  return workspace->integral_frames[depth];
}

inline std::size_t compiled_math_program_time_slot_count(
    const CompiledMathProgram &program) {
  semantic::Index max_time = static_cast<semantic::Index>(3);
  auto include_time = [&max_time](const semantic::Index time_id) {
    if (time_id != semantic::kInvalidIndex && time_id > max_time) {
      max_time = time_id;
    }
  };
  for (const auto &node : program.nodes) {
    include_time(node.time_id);
    if (node.kind == CompiledMathNodeKind::TimeGate ||
        compiled_math_is_source_value_node(node.kind)) {
      include_time(node.aux_id);
    }
  }
  for (const auto &kernel : program.integral_kernels) {
    include_time(kernel.bind_time_id);
  }
  for (const auto &channel : program.integral_kernel_source_product_channels) {
    include_time(channel.time_id);
    include_time(channel.time_cap_id);
  }
  for (const auto &term : program.source_condition_bound_terms) {
    include_time(term.time_id);
  }
  for (const auto &term : program.timed_upper_bound_terms) {
    include_time(term.time_id);
  }
  return static_cast<std::size_t>(max_time) + 1U;
}

inline bool compiled_math_source_product_ops_need_batch_cache(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops) {
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.cache_result) {
      return true;
    }
  }
  return false;
}

inline void compiled_math_add_cached_source_vector_program_id(
    const semantic::Index program_id,
    std::vector<semantic::Index> *program_ids) {
  if (program_id == semantic::kInvalidIndex || program_ids == nullptr) {
    return;
  }
  for (const auto existing : *program_ids) {
    if (existing == program_id) {
      return;
    }
  }
  program_ids->push_back(program_id);
}

inline void compiled_math_collect_cached_source_vector_programs_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    std::vector<semantic::Index> *program_ids) {
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.cache_result) {
      compiled_math_add_cached_source_vector_program_id(
          op.source_vector_program_id,
          program_ids);
    }
  }
}

inline void compiled_math_collect_kernel_cached_source_vector_programs(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    std::vector<semantic::Index> *program_ids) {
  if (program_ids == nullptr) {
    return;
  }
  program_ids->clear();
  const auto offset = static_cast<std::size_t>(
      kernel.source_product_block_cached_source_vector_programs.offset);
  const auto size = static_cast<std::size_t>(
      kernel.source_product_block_cached_source_vector_programs.size);
  if (offset + size >
      program.integral_kernel_source_product_block_cached_source_vector_programs
          .size()) {
    return;
  }
  for (std::size_t i = 0; i < size; ++i) {
    program_ids->push_back(
        program.integral_kernel_source_product_block_cached_source_vector_programs[
            offset + i]);
  }
}

inline bool compiled_math_kernel_needs_batch_source_cache(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel) {
  (void)program;
  return kernel.source_product_block_cached_source_vector_programs.size != 0;
}

inline void compiled_math_batch_write_impossible_fill_lane(
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  mask_by_lane[lane_pos] = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    pdf_by_lane[lane_pos] = 0.0;
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    cdf_by_lane[lane_pos] = 0.0;
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[lane_pos] = 1.0;
  }
}

inline void compiled_math_batch_write_certain_fill_lane(
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  mask_by_lane[lane_pos] = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    pdf_by_lane[lane_pos] = 0.0;
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    cdf_by_lane[lane_pos] = 1.0;
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[lane_pos] = 0.0;
  }
}

inline void compiled_math_batch_write_forced_fill_lane(
    const ExactRelation relation,
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  if (relation == ExactRelation::Before || relation == ExactRelation::At) {
    compiled_math_batch_write_certain_fill_lane(
        lane, fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
        survival_by_lane);
    return;
  }
  compiled_math_batch_write_impossible_fill_lane(
      lane, fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
      survival_by_lane);
  if (relation == ExactRelation::After &&
      (fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[static_cast<std::size_t>(lane)] = 1.0;
  }
}

inline void compiled_math_batch_write_fill_for_lanes(
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const std::uint8_t fill_mask,
    const bool certain,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  for (std::size_t i = 0; i < lane_count; ++i) {
    if (certain) {
      compiled_math_batch_write_certain_fill_lane(
          lanes[i], fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
          survival_by_lane);
    } else {
      compiled_math_batch_write_impossible_fill_lane(
          lanes[i], fill_mask, mask_by_lane, pdf_by_lane, cdf_by_lane,
          survival_by_lane);
    }
  }
}

inline void compiled_math_batch_copy_fill_lane(
    const semantic::Index lane,
    const std::uint8_t fill_mask,
    const std::uint8_t *src_mask_by_lane,
    const double *src_pdf_by_lane,
    const double *src_cdf_by_lane,
    const double *src_survival_by_lane,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  mask_by_lane[lane_pos] = src_mask_by_lane[lane_pos] | fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    pdf_by_lane[lane_pos] = src_pdf_by_lane[lane_pos];
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    cdf_by_lane[lane_pos] = src_cdf_by_lane[lane_pos];
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    survival_by_lane[lane_pos] = src_survival_by_lane[lane_pos];
  }
}

inline bool compiled_math_batch_lognormal_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      const auto *loaded =
          compiled_math_batch_leaf_input_for_lane(state, lane, leaf_index);
      if (loaded == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded->t0;
      const double m = loaded->params[0];
      const double s = loaded->params[1];
      if (!(x > 0.0) || !std::isfinite(m) || !std::isfinite(s) ||
          s <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded->q,
                      channel_mask);
        continue;
      }
      scratch->lanes_a[valid_count] = lane;
      scratch->lower[valid_count] = x;
      scratch->pdf_a[valid_count] = m;
      scratch->cdf_a[valid_count] = 1.0 / s;
      scratch->survival_a[valid_count] = loaded->q;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const double z =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
        scratch->upper[i] = -0.5 * z * z;
      }
      vvexp(scratch->upper.data(), scratch->upper.data(), &n);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double base_pdf =
            kInvSqrtTwoPi * scratch->upper[i] * scratch->cdf_a[i] /
            scratch->lower[i];
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                0.0,
                scratch->survival_a[i],
                channel_mask);
      }
    } else {
      for (std::size_t i = 0; i < valid_count; ++i) {
        scratch->upper[i] =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
      }
      compiled_math_batch_normal_cdf_from_finite_z(
          scratch->upper.data(),
          scratch->lower.data(),
          valid_count);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                scratch->upper[i],
                scratch->survival_a[i],
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return false;
}

inline void compiled_math_batch_gamma_regularized_lower_values(
    const double *shape_by_index,
    const double *scaled_x_by_index,
    const std::size_t count,
    double *out_by_index) {
  constexpr int kMaxIter = 1000;
  constexpr double kEps = 1e-14;
  constexpr double kTiny = 1e-300;
  for (std::size_t i = 0; i < count; ++i) {
    const double shape = shape_by_index[i];
    const double x = scaled_x_by_index[i];
    if (!std::isfinite(shape) || shape <= 0.0 || !std::isfinite(x) ||
        x <= 0.0) {
      out_by_index[i] = 0.0;
      continue;
    }
    const double log_scale = shape * std::log(x) - x - std::lgamma(shape);
    if (x < shape + 1.0) {
      double ap = shape;
      double delta = 1.0 / shape;
      double sum = delta;
      for (int n = 1; n <= kMaxIter; ++n) {
        ap += 1.0;
        delta *= x / ap;
        sum += delta;
        if (std::fabs(delta) <= std::fabs(sum) * kEps) {
          break;
        }
      }
      out_by_index[i] = clamp_probability(sum * std::exp(log_scale));
      continue;
    }

    double b = x + 1.0 - shape;
    double c = 1.0 / kTiny;
    double d = std::fabs(b) < kTiny ? kTiny : b;
    d = 1.0 / d;
    double h = d;
    for (int n = 1; n <= kMaxIter; ++n) {
      const double fn = static_cast<double>(n);
      const double an = -fn * (fn - shape);
      b += 2.0;
      d = an * d + b;
      if (std::fabs(d) < kTiny) {
        d = kTiny;
      }
      c = b + an / c;
      if (std::fabs(c) < kTiny) {
        c = kTiny;
      }
      d = 1.0 / d;
      const double delta = d * c;
      h *= delta;
      if (std::fabs(delta - 1.0) <= kEps) {
        break;
      }
    }
    const double upper = std::exp(log_scale) * h;
    out_by_index[i] = clamp_probability(1.0 - upper);
  }
}

inline bool compiled_math_batch_gamma_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      const auto *loaded =
          compiled_math_batch_leaf_input_for_lane(state, lane, leaf_index);
      if (loaded == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded->t0;
      const double shape = loaded->params[0];
      const double rate = loaded->params[1];
      if (!(x > 0.0) || !std::isfinite(shape) || shape <= 0.0 ||
          !std::isfinite(rate) || rate <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded->q,
                      channel_mask);
        continue;
      }
      scratch->lanes_a[valid_count] = lane;
      scratch->lower[valid_count] = x;
      scratch->pdf_a[valid_count] = shape;
      scratch->cdf_a[valid_count] = rate;
      scratch->survival_a[valid_count] = loaded->q;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    if (channel_mask != kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        scratch->upper[i] = scratch->cdf_a[i] * scratch->lower[i];
      }
      compiled_math_batch_gamma_regularized_lower_values(
          scratch->pdf_a.data(),
          scratch->upper.data(),
          valid_count,
          scratch->exact.data());
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto lane = scratch->lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                scratch->exact[i],
                scratch->survival_a[i],
                channel_mask);
      }
      return true;
    }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    for (std::size_t i = 0; i < valid_count; ++i) {
      const double shape = scratch->pdf_a[i];
      const double rate = scratch->cdf_a[i];
      scratch->upper[i] =
          shape * std::log(rate) + (shape - 1.0) * scratch->upper[i] -
          rate * scratch->lower[i] - std::lgamma(shape);
    }
    vvexp(scratch->upper.data(), scratch->upper.data(), &n);
#else
    for (std::size_t i = 0; i < valid_count; ++i) {
      scratch->upper[i] =
          std::exp(
              scratch->pdf_a[i] * std::log(scratch->cdf_a[i]) +
              (scratch->pdf_a[i] - 1.0) * std::log(scratch->lower[i]) -
              scratch->cdf_a[i] * scratch->lower[i] -
              std::lgamma(scratch->pdf_a[i]));
    }
#endif
    for (std::size_t i = 0; i < valid_count; ++i) {
      const auto lane = scratch->lanes_a[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      values_by_lane[lane_pos] =
          batch_source_product_finish_base_value(
              scratch->upper[i],
              0.0,
              scratch->survival_a[i],
              channel_mask);
    }
    return true;
  }
  return false;
}

inline bool compiled_math_batch_exgauss_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      const auto *loaded =
          compiled_math_batch_leaf_input_for_lane(state, lane, leaf_index);
      if (loaded == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded->t0;
      const double mu = loaded->params[0];
      const double sigma = loaded->params[1];
      const double tau = loaded->params[2];
      if (!(x > 0.0) || !std::isfinite(mu) || !std::isfinite(sigma) ||
          sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
        values_by_lane[lane_pos] =
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      loaded->q,
                      channel_mask);
        continue;
      }
      const double inv_tau = 1.0 / tau;
      const double sigma_over_tau = sigma * inv_tau;
      const double z0 = -mu / sigma;
      scratch->times[lane_pos] = x;
      scratch->pdf_b[lane_pos] = mu;
      scratch->cdf_b[lane_pos] = sigma;
      scratch->survival_b[lane_pos] = tau;
      scratch->survival_a[lane_pos] = loaded->q;
      scratch->lanes_a[valid_count] = lane;
      scratch->upper[valid_count] = z0;
      scratch->exact[valid_count] = z0 - sigma_over_tau;
      scratch->lower[valid_count] =
          sigma * sigma * inv_tau * inv_tau * 0.5 + mu * inv_tau;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }
    const auto n = static_cast<int>(valid_count);
    vvexp(scratch->lower.data(), scratch->lower.data(), &n);
    compiled_math_batch_normal_cdf_from_z(
        scratch->upper.data(),
        scratch->cdf_a.data(),
        valid_count);
    compiled_math_batch_normal_cdf_from_z(
        scratch->exact.data(),
        scratch->cdf_a.data(),
        valid_count);
    std::size_t live_count = 0U;
    for (std::size_t i = 0; i < valid_count; ++i) {
      const auto lane = scratch->lanes_a[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double lower_cdf =
          clamp_probability(scratch->upper[i] - scratch->lower[i] * scratch->exact[i]);
      const double lower_survival = 1.0 - lower_cdf;
      if (!(lower_survival > 0.0)) {
        values_by_lane[lane_pos] =
            batch_source_product_before_onset_value(channel_mask);
        continue;
      }
      scratch->cdf_c[lane_pos] = lower_cdf;
      scratch->survival_c[lane_pos] = lower_survival;
      scratch->lanes_b[live_count++] = lane;
    }
    if (live_count == 0U) {
      return true;
    }
    for (std::size_t i = 0; i < live_count; ++i) {
      const auto lane = scratch->lanes_b[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double x = scratch->times[lane_pos];
      const double mu = scratch->pdf_b[lane_pos];
      const double sigma = scratch->cdf_b[lane_pos];
      const double tau = scratch->survival_b[lane_pos];
      const double inv_tau = 1.0 / tau;
      const double sigma_over_tau = sigma * inv_tau;
      const double z = (x - mu) / sigma;
      scratch->pdf_a[i] = tau;
      scratch->pdf_c[i] = scratch->survival_a[lane_pos];
      scratch->onset_upper[i] = scratch->cdf_c[lane_pos];
      scratch->lower[i] =
          sigma * sigma * inv_tau * inv_tau * 0.5 -
          (x - mu) * inv_tau;
      scratch->exact[i] = z - sigma_over_tau;
      if (channel_mask != kLeafChannelPdf) {
        scratch->upper[i] = z;
      }
    }
    const auto live_n = static_cast<int>(live_count);
    vvexp(scratch->lower.data(), scratch->lower.data(), &live_n);
    compiled_math_batch_normal_cdf_from_z(
        scratch->exact.data(),
        scratch->cdf_a.data(),
        live_count);
    if (channel_mask != kLeafChannelPdf) {
      compiled_math_batch_normal_cdf_from_z(
          scratch->upper.data(),
          scratch->cdf_a.data(),
          live_count);
    }
    for (std::size_t i = 0; i < live_count; ++i) {
      const auto lane = scratch->lanes_b[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double lower_cdf = scratch->onset_upper[i];
      const double lower_survival = 1.0 - lower_cdf;
      if (channel_mask == kLeafChannelPdf) {
        const double base_pdf =
            (scratch->lower[i] * scratch->exact[i] / scratch->pdf_a[i]) /
            lower_survival;
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                0.0,
                scratch->pdf_c[i],
                channel_mask);
      } else {
        const double raw_cdf =
            clamp_probability(scratch->upper[i] -
                              scratch->lower[i] * scratch->exact[i]);
        const double base_cdf =
            (raw_cdf - lower_cdf) / lower_survival;
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                base_cdf,
                scratch->pdf_c[i],
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return false;
}

inline bool compiled_math_batch_lba_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t valid_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      const auto *loaded =
          compiled_math_batch_leaf_input_for_lane(state, lane, leaf_index);
      if (loaded == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded->t0;
      const double v = loaded->params[0];
      const double B = loaded->params[1];
      const double A = loaded->params[2];
      const double sv = loaded->params[3];
      if (!(x > 0.0)) {
        values_by_lane[lane_pos] =
            batch_source_product_before_onset_value(channel_mask);
        continue;
      }
      if (!std::isfinite(v) || !std::isfinite(B) ||
          !std::isfinite(A) || !std::isfinite(sv) || sv <= 0.0) {
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                0.0,
                loaded->q,
                channel_mask);
        continue;
      }
      scratch->lanes_a[valid_count] = lane;
      scratch->lower[valid_count] = x;
      scratch->upper[valid_count] = v;
      scratch->exact[valid_count] = B;
      scratch->pdf_a[valid_count] = A;
      scratch->cdf_a[valid_count] = sv;
      scratch->survival_a[valid_count] = loaded->q;
      scratch->pdf_b[valid_count] = v / sv;
      ++valid_count;
    }
    if (valid_count == 0U) {
      return true;
    }

    compiled_math_batch_normal_cdf_from_z(
        scratch->pdf_b.data(),
        scratch->survival_b.data(),
        valid_count);
    for (std::size_t i = 0; i < valid_count; ++i) {
      if (!std::isfinite(scratch->pdf_b[i]) || scratch->pdf_b[i] < 1e-10) {
        scratch->pdf_b[i] = 1e-10;
      }
    }

    std::size_t wide_count = 0U;
    std::size_t point_count = 0U;
    for (std::size_t i = 0; i < valid_count; ++i) {
      const double x = scratch->lower[i];
      const double v = scratch->upper[i];
      const double B = scratch->exact[i];
      const double A = scratch->pdf_a[i];
      const double sv = scratch->cdf_a[i];
      if (A > 1e-10) {
        const double zs = x * sv;
        if (!(zs > 0.0) || !std::isfinite(zs)) {
          const auto lane = scratch->lanes_a[i];
          values_by_lane[static_cast<std::size_t>(lane)] =
              batch_source_product_finish_base_value(
                  0.0,
                  0.0,
                  scratch->survival_a[i],
                  channel_mask);
          continue;
        }
        const double cmz = B - x * v;
        const double xx = cmz - A;
        const double cz = cmz / zs;
        const double cz_max = xx / zs;
        scratch->lanes_b[wide_count] = static_cast<semantic::Index>(i);
        scratch->cdf_b[wide_count] = cz;
        scratch->survival_b[wide_count] = cz_max;
        scratch->times[wide_count] = -0.5 * cz * cz;
        scratch->onset_upper[wide_count] = -0.5 * cz_max * cz_max;
        ++wide_count;
      } else {
        const double z = (B / x - v) / sv;
        scratch->lanes_c[point_count] = static_cast<semantic::Index>(i);
        scratch->cdf_c[point_count] = z;
        scratch->survival_c[point_count] = -0.5 * z * z;
        ++point_count;
      }
    }

    if (wide_count != 0U) {
      const auto n = static_cast<int>(wide_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->cdf_b.data(),
          scratch->mask_a.empty() ? nullptr : scratch->pdf_c.data(),
          wide_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->survival_b.data(),
          scratch->mask_a.empty() ? nullptr : scratch->pdf_c.data(),
          wide_count);
      vvexp(scratch->times.data(), scratch->times.data(), &n);
      vvexp(scratch->onset_upper.data(), scratch->onset_upper.data(), &n);
      for (std::size_t j = 0; j < wide_count; ++j) {
        const auto src = static_cast<std::size_t>(scratch->lanes_b[j]);
        const auto lane = scratch->lanes_a[src];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double x = scratch->lower[src];
        const double v = scratch->upper[src];
        const double B = scratch->exact[src];
        const double A = scratch->pdf_a[src];
        const double sv = scratch->cdf_a[src];
        const double denom = scratch->pdf_b[src];
        const double cmz = B - x * v;
        const double xx = cmz - A;
        const double phi_cz = kInvSqrtTwoPi * scratch->times[j];
        const double phi_cz_max = kInvSqrtTwoPi * scratch->onset_upper[j];
        double base_pdf = 0.0;
        double base_cdf = 0.0;
        if (channel_mask == kLeafChannelPdf) {
          base_pdf =
              (v * (scratch->cdf_b[j] - scratch->survival_b[j]) +
               sv * (phi_cz_max - phi_cz)) /
              (A * denom);
        } else {
          base_cdf =
              (1.0 + (x * sv * (phi_cz_max - phi_cz) +
                      xx * scratch->survival_b[j] -
                      cmz * scratch->cdf_b[j]) /
                         A) /
              denom;
        }
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                base_cdf,
                scratch->survival_a[src],
                channel_mask);
      }
    }

    if (point_count != 0U) {
      compiled_math_batch_normal_cdf_from_z(
          scratch->cdf_c.data(),
          scratch->pdf_c.data(),
          point_count);
      const auto n = static_cast<int>(point_count);
      vvexp(scratch->survival_c.data(), scratch->survival_c.data(), &n);
      for (std::size_t j = 0; j < point_count; ++j) {
        const auto src = static_cast<std::size_t>(scratch->lanes_c[j]);
        const auto lane = scratch->lanes_a[src];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double x = scratch->lower[src];
        const double B = scratch->exact[src];
        const double sv = scratch->cdf_a[src];
        const double denom = scratch->pdf_b[src];
        double base_pdf = 0.0;
        double base_cdf = 0.0;
        if (channel_mask == kLeafChannelPdf) {
          base_pdf =
              kInvSqrtTwoPi * scratch->survival_c[j] * B /
              (sv * x * x * denom);
        } else {
          base_cdf = clamp_probability(1.0 - scratch->cdf_c[j]) / denom;
        }
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                base_cdf,
                scratch->survival_a[src],
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return false;
}

inline bool compiled_math_batch_rdm_leaf_values_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || values_by_lane == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, lane_count);
    std::size_t common_count = 0U;
    std::size_t point_count = 0U;
    std::size_t small_drift_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      const auto *loaded =
          compiled_math_batch_leaf_input_for_lane(state, lane, leaf_index);
      if (loaded == nullptr) {
        return false;
      }
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double x =
          current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded->t0;
      const double v = loaded->params[0];
      const double B = loaded->params[1];
      const double A = loaded->params[2];
      const double s = loaded->params[3];
      if (!(x > 0.0)) {
        values_by_lane[lane_pos] =
            batch_source_product_before_onset_value(channel_mask);
        continue;
      }
      if (!std::isfinite(s) || s <= 0.0 || !std::isfinite(v) ||
          !std::isfinite(B) || !std::isfinite(A)) {
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                0.0,
                loaded->q,
                channel_mask);
        continue;
      }
      const double v_sc = v / s;
      const double B_sc = B / s;
      const double A_sc = A / s;
      const double a = 0.5 * A_sc;
      const double k = B_sc + a;
      const double l = v_sc;
      if (!std::isfinite(a) || !std::isfinite(k) || !std::isfinite(l)) {
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                0.0,
                loaded->q,
                channel_mask);
        continue;
      }
      if (a < 1e-10) {
        scratch->lanes_c[point_count] = lane;
        scratch->cdf_a[point_count] = x;
        scratch->cdf_b[point_count] = k;
        scratch->cdf_c[point_count] = l;
        scratch->survival_c[point_count] = loaded->q;
        ++point_count;
        continue;
      }
      if (l < 1e-10) {
        scratch->lanes_d[small_drift_count] = lane;
        scratch->pdf_b[small_drift_count] = x;
        scratch->survival_b[small_drift_count] = k;
        scratch->pdf_c[small_drift_count] = l;
        scratch->onset_upper[small_drift_count] = a;
        scratch->survival_c[small_drift_count] = loaded->q;
        ++small_drift_count;
        continue;
      }
      scratch->lanes_b[common_count] = lane;
      scratch->lower[common_count] = x;
      scratch->upper[common_count] = k;
      scratch->exact[common_count] = l;
      scratch->pdf_a[common_count] = a;
      scratch->survival_a[common_count] = loaded->q;
      ++common_count;
    }
    if (point_count != 0U) {
      if (channel_mask == kLeafChannelPdf) {
        for (std::size_t i = 0; i < point_count; ++i) {
          const double x = scratch->cdf_a[i];
          const double k = scratch->cdf_b[i];
          const double l = scratch->cdf_c[i];
          const double lambda = k * k;
          double exponent = 0.0;
          if (l == 0.0) {
            exponent = -0.5 * lambda / x;
          } else {
            const double mu = k / l;
            exponent =
                -(lambda / (2.0 * x)) *
                ((x * x) / (mu * mu) - 2.0 * x / mu + 1.0);
          }
          scratch->times[i] =
              exponent + 0.5 * std::log(lambda) -
              0.5 * std::log(2.0 * x * x * x * M_PI);
        }
        const auto n = static_cast<int>(point_count);
        vvexp(scratch->times.data(), scratch->times.data(), &n);
        for (std::size_t i = 0; i < point_count; ++i) {
          const auto lane = scratch->lanes_c[i];
          const auto lane_pos = static_cast<std::size_t>(lane);
          values_by_lane[lane_pos] =
              batch_source_product_finish_base_value(
                  scratch->times[i],
                  0.0,
                  scratch->survival_c[i],
                  channel_mask);
        }
      } else {
        for (std::size_t i = 0; i < point_count; ++i) {
          const double x = scratch->cdf_a[i];
          const double k = scratch->cdf_b[i];
          const double l = scratch->cdf_c[i];
          if (std::fabs(l) < 1e-12) {
            scratch->times[i] = k / std::sqrt(x);
            scratch->onset_upper[i] =
                std::numeric_limits<double>::quiet_NaN();
          } else {
            const double mu = k / l;
            const double z_base = std::sqrt((k * k) / x);
            scratch->times[i] = z_base * (1.0 + x / mu);
            scratch->onset_upper[i] = z_base * (1.0 - x / mu);
          }
        }
        compiled_math_batch_normal_cdf_from_z(
            scratch->times.data(),
            scratch->lower.data(),
            point_count);
        compiled_math_batch_normal_cdf_from_z(
            scratch->onset_upper.data(),
            scratch->lower.data(),
            point_count);
        for (std::size_t i = 0; i < point_count; ++i) {
          const auto lane = scratch->lanes_c[i];
          const auto lane_pos = static_cast<std::size_t>(lane);
          const double k = scratch->cdf_b[i];
          const double l = scratch->cdf_c[i];
          double base_cdf = 0.0;
          if (std::fabs(l) < 1e-12) {
            base_cdf = 2.0 * (1.0 - scratch->times[i]);
          } else {
            const double p1 = 1.0 - scratch->times[i];
            const double p2 = 1.0 - scratch->onset_upper[i];
            base_cdf =
                std::exp(
                    2.0 * k * l +
                    std::log(std::max(1e-300, p1))) +
                p2;
          }
          values_by_lane[lane_pos] =
              batch_source_product_finish_base_value(
                  0.0,
                  base_cdf,
                  scratch->survival_c[i],
                  channel_mask);
        }
      }
    }

    if (small_drift_count != 0U) {
      if (channel_mask == kLeafChannelPdf) {
        for (std::size_t i = 0; i < small_drift_count; ++i) {
          const double x = scratch->pdf_b[i];
          const double k = scratch->survival_b[i];
          const double a = scratch->onset_upper[i];
          scratch->times[i] = -((k - a) * (k - a)) / (2.0 * x);
          scratch->upper[i] = -((k + a) * (k + a)) / (2.0 * x);
        }
        const auto n = static_cast<int>(small_drift_count);
        vvexp(scratch->times.data(), scratch->times.data(), &n);
        vvexp(scratch->upper.data(), scratch->upper.data(), &n);
        for (std::size_t i = 0; i < small_drift_count; ++i) {
          const auto lane = scratch->lanes_d[i];
          const auto lane_pos = static_cast<std::size_t>(lane);
          const double x = scratch->pdf_b[i];
          const double a = scratch->onset_upper[i];
          const double term = scratch->times[i] - scratch->upper[i];
          const double base_pdf =
              std::exp(
                  -0.5 * (M_LN2 + kLogPi + std::log(x)) +
                  std::log(std::max(1e-300, term)) -
                  M_LN2 - std::log(a));
          values_by_lane[lane_pos] =
              batch_source_product_finish_base_value(
                  base_pdf,
                  0.0,
                  scratch->survival_c[i],
                  channel_mask);
        }
      } else {
        for (std::size_t i = 0; i < small_drift_count; ++i) {
          const double x = scratch->pdf_b[i];
          const double k = scratch->survival_b[i];
          const double a = scratch->onset_upper[i];
          const double sqt = std::sqrt(x);
          scratch->times[i] = (k + a) / sqt;
          scratch->upper[i] = (-k - a) / sqt;
          scratch->lower[i] =
              -0.5 * (((k + a) * (k + a)) / x - M_LN2 - kLogPi +
                      std::log(x)) -
              std::log(a);
          scratch->exact[i] =
              -0.5 * (((k - a) * (k - a)) / x - M_LN2 - kLogPi +
                      std::log(x)) -
              std::log(a);
        }
        compiled_math_batch_normal_cdf_from_z(
            scratch->times.data(),
            scratch->cdf_b.data(),
            small_drift_count);
        compiled_math_batch_normal_cdf_from_z(
            scratch->upper.data(),
            scratch->cdf_b.data(),
            small_drift_count);
        const auto n = static_cast<int>(small_drift_count);
        vvexp(scratch->lower.data(), scratch->lower.data(), &n);
        vvexp(scratch->exact.data(), scratch->exact.data(), &n);
        for (std::size_t i = 0; i < small_drift_count; ++i) {
          const auto lane = scratch->lanes_d[i];
          const auto lane_pos = static_cast<std::size_t>(lane);
          const double k = scratch->survival_b[i];
          const double a = scratch->onset_upper[i];
          const double t5a = 2.0 * scratch->times[i] - 1.0;
          const double t5b = 2.0 * scratch->upper[i] - 1.0;
          const double base_cdf =
              1.0 + scratch->lower[i] - scratch->exact[i] +
              ((-k + a) * t5a - (k - a) * t5b) / (2.0 * a);
          values_by_lane[lane_pos] =
              batch_source_product_finish_base_value(
                  0.0,
                  base_cdf,
                  scratch->survival_c[i],
                  channel_mask);
        }
      }
    }

    if (common_count == 0U) {
      return true;
    }

    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < common_count; ++i) {
        const double x = scratch->lower[i];
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double sqt = std::sqrt(x);
        scratch->cdf_b[i] = -std::pow(a - k + x * l, 2.0) / (2.0 * x);
        scratch->survival_b[i] =
            -std::pow(a + k - x * l, 2.0) / (2.0 * x);
        scratch->times[i] = (-k + a) / sqt + sqt * l;
        scratch->onset_upper[i] = (k + a) / sqt - sqt * l;
      }
      const auto n = static_cast<int>(common_count);
      vvexp(scratch->cdf_b.data(), scratch->cdf_b.data(), &n);
      vvexp(scratch->survival_b.data(), scratch->survival_b.data(), &n);
      compiled_math_batch_normal_cdf_from_z(
          scratch->times.data(),
          scratch->pdf_b.data(),
          common_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->onset_upper.data(),
          scratch->pdf_b.data(),
          common_count);
      for (std::size_t i = 0; i < common_count; ++i) {
        const auto lane = scratch->lanes_b[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double x = scratch->lower[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double sqt = std::sqrt(x);
        const double t1 =
            kInvSqrtTwo *
            (scratch->cdf_b[i] - scratch->survival_b[i]) /
            (std::sqrt(M_PI) * sqt);
        const double t2 =
            0.5 * l *
            ((2.0 * scratch->times[i] - 1.0) +
             (2.0 * scratch->onset_upper[i] - 1.0));
        const double base_pdf =
            std::exp(std::log(std::max(1e-300, t1 + t2)) -
                     M_LN2 - std::log(a));
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                base_pdf,
                0.0,
                scratch->survival_a[i],
                channel_mask);
      }
    } else {
      for (std::size_t i = 0; i < common_count; ++i) {
        const double x = scratch->lower[i];
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double sqt = std::sqrt(x);
        scratch->cdf_b[i] = -0.5 * std::pow(k - a - x * l, 2.0) / x;
        scratch->survival_b[i] =
            -0.5 * std::pow(a + k - x * l, 2.0) / x;
        scratch->pdf_c[i] = 0.5 * (std::log(x) - M_LN2 - kLogPi);
        scratch->times[i] = -(k - a + x * l) / sqt;
        scratch->onset_upper[i] = -(k + a + x * l) / sqt;
        scratch->pdf_b[i] = (k + a) / sqt - sqt * l;
        scratch->cdf_c[i] = (k - a) / sqt - sqt * l;
      }
      const auto n = static_cast<int>(common_count);
      vvexp(scratch->cdf_b.data(), scratch->cdf_b.data(), &n);
      vvexp(scratch->survival_b.data(), scratch->survival_b.data(), &n);
      vvexp(scratch->pdf_c.data(), scratch->pdf_c.data(), &n);
      compiled_math_batch_normal_log_cdf_from_z(
          scratch->times.data(),
          scratch->survival_c.data(),
          common_count);
      compiled_math_batch_normal_log_cdf_from_z(
          scratch->onset_upper.data(),
          scratch->survival_c.data(),
          common_count);
      for (std::size_t i = 0; i < common_count; ++i) {
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        scratch->times[i] += 2.0 * l * (k - a);
        scratch->onset_upper[i] += 2.0 * l * (k + a);
      }
      vvexp(scratch->times.data(), scratch->times.data(), &n);
      vvexp(scratch->onset_upper.data(), scratch->onset_upper.data(), &n);
      compiled_math_batch_normal_cdf_from_z(
          scratch->pdf_b.data(),
          scratch->survival_c.data(),
          common_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->cdf_c.data(),
          scratch->survival_c.data(),
          common_count);
      for (std::size_t i = 0; i < common_count; ++i) {
        const auto lane = scratch->lanes_b[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        const double x = scratch->lower[i];
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double t1 =
            scratch->pdf_c[i] *
            (scratch->cdf_b[i] - scratch->survival_b[i]);
        const double t2 =
            a + (scratch->onset_upper[i] - scratch->times[i]) /
                    (2.0 * l);
        const double t4a = 2.0 * scratch->pdf_b[i] - 1.0;
        const double t4b = 2.0 * scratch->cdf_c[i] - 1.0;
        const double t4 =
            0.5 * (x * l - a - k + 0.5 / l) * t4a +
            0.5 * (k - a - x * l - 0.5 / l) * t4b;
        const double base_cdf = 0.5 * (t4 + t2 + t1) / a;
        values_by_lane[lane_pos] =
            batch_source_product_finish_base_value(
                0.0,
                base_cdf,
                scratch->survival_a[i],
                channel_mask);
      }
    }
    return true;
  }
#else
  (void)scratch;
#endif
  return false;
}

inline bool compiled_math_batch_leaf_values_from_times(
    const std::uint8_t leaf_dist_kind,
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  switch (static_cast<leaf::DistKind>(leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    return compiled_math_batch_lognormal_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  case leaf::DistKind::Gamma:
    return compiled_math_batch_gamma_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  case leaf::DistKind::Exgauss:
    return compiled_math_batch_exgauss_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  case leaf::DistKind::LBA:
    return compiled_math_batch_lba_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  case leaf::DistKind::RDM:
    return compiled_math_batch_rdm_leaf_values_from_times(
        leaf_index,
        leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        channel_mask,
        values_by_lane,
        scratch);
  }
  return false;
}

inline bool compiled_math_batch_lognormal_leaf_fill_from_times(
    const semantic::Index leaf_index,
    const double leaf_onset_abs_value,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (leaf_index == semantic::kInvalidIndex ||
      current_time_by_lane == nullptr || scratch == nullptr) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf = (fill_mask & kLeafChannelCdf) != 0U;
  const bool need_survival = (fill_mask & kLeafChannelSurvival) != 0U;
  scratch->ensure_size(state.lane_count, lane_count);
  std::size_t valid_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const auto *loaded =
        compiled_math_batch_leaf_input_for_lane(state, lane, leaf_index);
    if (loaded == nullptr) {
      return false;
    }
    const double x =
        current_time_by_lane[lane_pos] - leaf_onset_abs_value - loaded->t0;
    const double m = loaded->params[0];
    const double s = loaded->params[1];
    if (!(x > 0.0) || !std::isfinite(m) || !std::isfinite(s) ||
        s <= 0.0) {
      if (need_pdf) {
        pdf_by_lane[lane_pos] = 0.0;
      }
      if (need_cdf) {
        cdf_by_lane[lane_pos] = 0.0;
      }
      if (need_survival) {
        survival_by_lane[lane_pos] = 1.0;
      }
      continue;
    }
    scratch->lanes_a[valid_count] = lane;
    scratch->lower[valid_count] = x;
    scratch->pdf_a[valid_count] = m;
    scratch->cdf_a[valid_count] = 1.0 / s;
    scratch->survival_a[valid_count] = loaded->q;
    ++valid_count;
  }
  if (valid_count == 0U) {
    return true;
  }

  const auto n = static_cast<int>(valid_count);
  vvlog(scratch->upper.data(), scratch->lower.data(), &n);
  for (std::size_t i = 0; i < valid_count; ++i) {
    scratch->upper[i] =
        (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
  }

  if (need_pdf) {
    for (std::size_t i = 0; i < valid_count; ++i) {
      const double z = scratch->upper[i];
      scratch->pdf_b[i] = -0.5 * z * z;
    }
    vvexp(scratch->pdf_b.data(), scratch->pdf_b.data(), &n);
    for (std::size_t i = 0; i < valid_count; ++i) {
      const auto lane = scratch->lanes_a[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double base_pdf =
          kInvSqrtTwoPi * scratch->pdf_b[i] * scratch->cdf_a[i] /
          scratch->lower[i];
      pdf_by_lane[lane_pos] =
          batch_source_product_finish_base_value(
              base_pdf,
              0.0,
              scratch->survival_a[i],
              kLeafChannelPdf);
    }
  }

  if (need_cdf || need_survival) {
    compiled_math_batch_normal_cdf_from_finite_z(
        scratch->upper.data(),
        scratch->cdf_b.data(),
        valid_count);
    for (std::size_t i = 0; i < valid_count; ++i) {
      const auto lane = scratch->lanes_a[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double cdf =
          batch_source_product_finish_base_value(
              0.0,
              scratch->upper[i],
              scratch->survival_a[i],
              kLeafChannelCdf);
      if (need_cdf) {
        cdf_by_lane[lane_pos] = cdf;
      }
      if (need_survival) {
        survival_by_lane[lane_pos] = clamp_probability(1.0 - cdf);
      }
    }
  }
  return true;
#else
  (void)fill_mask;
  (void)pdf_by_lane;
  (void)cdf_by_lane;
  (void)survival_by_lane;
  return false;
#endif
}

inline bool compiled_math_batch_resolve_source_bounds_for_source(
    const CompiledMathProgram &program,
    semantic::Index source_id,
    const CompiledSourceBoundPlan &bounds,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    std::size_t lane_count,
    double *lower_by_lane,
    double *upper_by_lane,
    double *exact_by_lane,
    std::uint8_t *has_exact_by_lane);

inline bool compiled_math_batch_resolve_source_exact_bounds_for_source(
    semantic::Index source_id,
    const CompiledSourceBoundPlan &bounds,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    std::size_t lane_count,
    double *exact_by_lane,
    std::uint8_t *has_exact_by_lane);

inline bool compiled_math_batch_source_vector_op_leaf_fill_values(
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    std::size_t lane_count,
    const double *current_time_by_lane,
    std::uint8_t fill_mask,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane,
    CompiledMathBatchSourceProductScratch *scratch);

inline bool compiled_math_batch_source_vector_op_leaf_fill(
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    std::size_t lane_count,
    const double *current_time_by_lane,
    std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane);

inline bool compiled_math_source_vector_slot_arrays(
    CompiledMathBatchSourceProductScratch *scratch,
    const std::size_t lane_stride,
    const semantic::Index slot_count,
    const semantic::Index mask_slot,
    const semantic::Index pdf_slot,
    const semantic::Index cdf_slot,
    const semantic::Index survival_slot,
    std::uint8_t **mask_by_lane,
    double **pdf_by_lane,
    double **cdf_by_lane,
    double **survival_by_lane) {
  if (scratch == nullptr || mask_slot == semantic::kInvalidIndex ||
      pdf_slot == semantic::kInvalidIndex ||
      cdf_slot == semantic::kInvalidIndex ||
      survival_slot == semantic::kInvalidIndex ||
      static_cast<std::size_t>(mask_slot) >=
          static_cast<std::size_t>(slot_count) ||
      static_cast<std::size_t>(pdf_slot) >=
          static_cast<std::size_t>(slot_count) ||
      static_cast<std::size_t>(cdf_slot) >=
          static_cast<std::size_t>(slot_count) ||
      static_cast<std::size_t>(survival_slot) >=
          static_cast<std::size_t>(slot_count)) {
    return false;
  }
  const auto mask_offset = static_cast<std::size_t>(mask_slot) * lane_stride;
  const auto pdf_offset = static_cast<std::size_t>(pdf_slot) * lane_stride;
  const auto cdf_offset = static_cast<std::size_t>(cdf_slot) * lane_stride;
  const auto survival_offset =
      static_cast<std::size_t>(survival_slot) * lane_stride;
  if (mask_offset + lane_stride > scratch->source_vector_mask_slots.size() ||
      pdf_offset + lane_stride > scratch->source_vector_pdf_slots.size() ||
      cdf_offset + lane_stride > scratch->source_vector_cdf_slots.size() ||
      survival_offset + lane_stride >
          scratch->source_vector_survival_slots.size()) {
    return false;
  }
  *mask_by_lane = scratch->source_vector_mask_slots.data() + mask_offset;
  *pdf_by_lane = scratch->source_vector_pdf_slots.data() + pdf_offset;
  *cdf_by_lane = scratch->source_vector_cdf_slots.data() + cdf_offset;
  *survival_by_lane =
      scratch->source_vector_survival_slots.data() + survival_offset;
  return true;
}

inline bool compiled_math_source_vector_const_slot_arrays(
    const CompiledMathBatchSourceProductScratch *scratch,
    const std::size_t lane_stride,
    const semantic::Index slot_count,
    const semantic::Index mask_slot,
    const semantic::Index pdf_slot,
    const semantic::Index cdf_slot,
    const semantic::Index survival_slot,
    const std::uint8_t **mask_by_lane,
    const double **pdf_by_lane,
    const double **cdf_by_lane,
    const double **survival_by_lane) {
  if (scratch == nullptr || mask_slot == semantic::kInvalidIndex ||
      pdf_slot == semantic::kInvalidIndex ||
      cdf_slot == semantic::kInvalidIndex ||
      survival_slot == semantic::kInvalidIndex ||
      static_cast<std::size_t>(mask_slot) >=
          static_cast<std::size_t>(slot_count) ||
      static_cast<std::size_t>(pdf_slot) >=
          static_cast<std::size_t>(slot_count) ||
      static_cast<std::size_t>(cdf_slot) >=
          static_cast<std::size_t>(slot_count) ||
      static_cast<std::size_t>(survival_slot) >=
          static_cast<std::size_t>(slot_count)) {
    return false;
  }
  const auto mask_offset = static_cast<std::size_t>(mask_slot) * lane_stride;
  const auto pdf_offset = static_cast<std::size_t>(pdf_slot) * lane_stride;
  const auto cdf_offset = static_cast<std::size_t>(cdf_slot) * lane_stride;
  const auto survival_offset =
      static_cast<std::size_t>(survival_slot) * lane_stride;
  if (mask_offset + lane_stride > scratch->source_vector_mask_slots.size() ||
      pdf_offset + lane_stride > scratch->source_vector_pdf_slots.size() ||
      cdf_offset + lane_stride > scratch->source_vector_cdf_slots.size() ||
      survival_offset + lane_stride >
          scratch->source_vector_survival_slots.size()) {
    return false;
  }
  *mask_by_lane = scratch->source_vector_mask_slots.data() + mask_offset;
  *pdf_by_lane = scratch->source_vector_pdf_slots.data() + pdf_offset;
  *cdf_by_lane = scratch->source_vector_cdf_slots.data() + cdf_offset;
  *survival_by_lane =
      scratch->source_vector_survival_slots.data() + survival_offset;
  return true;
}

inline bool compiled_math_batch_source_vector_exact_gate_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchSourceProductScratch *scratch,
    const std::uint8_t *child_mask_by_lane,
    const double *child_pdf_by_lane,
    const double *child_cdf_by_lane,
    const double *child_survival_by_lane,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  (void)child_mask_by_lane;
  const auto relation = compiled_math_source_vector_op_relation(op);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !compiled_math_source_product_bounds_have_overlay(
          op.bounds)) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      compiled_math_batch_write_forced_fill_lane(
          relation,
          lanes[i],
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
    }
    return true;
  }
  if (scratch == nullptr ||
      !compiled_math_batch_resolve_source_exact_bounds_for_source(
          op.source_id,
          op.bounds,
          state,
          lanes,
          lane_count,
          scratch->exact.data(),
          scratch->has_exact.data())) {
    return false;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    if (scratch->has_exact[lane_pos] != 0U) {
      if (current_time_by_lane[lane_pos] >= scratch->exact[lane_pos]) {
        compiled_math_batch_write_certain_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      } else {
        compiled_math_batch_write_impossible_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      }
      continue;
    }
    mask_by_lane[lane_pos] = fill_mask;
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      pdf_by_lane[lane_pos] = child_pdf_by_lane[lane_pos];
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      cdf_by_lane[lane_pos] = child_cdf_by_lane[lane_pos];
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      survival_by_lane[lane_pos] = child_survival_by_lane[lane_pos];
    }
  }
  return true;
}

inline bool compiled_math_execute_source_vector_ops(
    const CompiledMathProgram &program,
    const CompiledMathSourceVectorProgram &vector_program,
    const CompiledMathIndexSpan ops,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    std::size_t lane_count,
    const double *current_time_by_lane,
    std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t scratch_depth,
    CompiledMathBatchSourceProductScratch *scratch,
    std::size_t recursion_depth);

inline std::uint8_t compiled_math_source_vector_execution_fill_mask(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan ops,
    const std::uint8_t fill_mask) {
  std::uint8_t execution_fill_mask = fill_mask;
  for (semantic::Index op_idx = 0; op_idx < ops.size; ++op_idx) {
    const auto op_pos = static_cast<std::size_t>(ops.offset + op_idx);
    if (op_pos >= program.integral_kernel_source_vector_ops.size()) {
      continue;
    }
    const auto kind = program.integral_kernel_source_vector_ops[op_pos].kind;
    if (kind == CompiledMathSourceVectorOpKind::Conditioned &&
        (execution_fill_mask & kLeafChannelSurvival) != 0U) {
      execution_fill_mask |= kLeafChannelCdf;
    } else if (kind == CompiledMathSourceVectorOpKind::PoolKOfN) {
      execution_fill_mask |=
          static_cast<std::uint8_t>(kLeafChannelCdf | kLeafChannelSurvival);
    } else if (kind == CompiledMathSourceVectorOpKind::OnsetConvolution &&
               (execution_fill_mask & kLeafChannelSurvival) != 0U) {
      execution_fill_mask |= kLeafChannelCdf;
    }
  }
  return execution_fill_mask;
}

inline bool compiled_math_batch_source_vector_child_cdf_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceVectorProgram &vector_program,
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *time_by_lane,
    const std::uint8_t channel_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    double *cdf_out_by_lane,
    double *survival_out_by_lane,
    const std::size_t recursion_depth) {
  if (lane_count == 0U) {
    return true;
  }
  auto &child_scratch = frame->source_product_scratch_layer(
      scratch_depth, state.lane_count, lane_count);
  child_scratch.ensure_source_vector_slots(
      state.lane_count,
      static_cast<std::size_t>(vector_program.slot_count));
  if (!compiled_math_execute_source_vector_ops(
          program,
          vector_program,
          op.child_ops,
          state,
          lanes,
          lane_count,
          time_by_lane,
          channel_mask,
          frame,
          scratch_depth,
          &child_scratch,
          recursion_depth + 1U)) {
    return false;
  }
  const std::uint8_t *child_mask = nullptr;
  const double *child_pdf = nullptr;
  const double *child_cdf = nullptr;
  const double *child_survival = nullptr;
  if (!compiled_math_source_vector_const_slot_arrays(
          &child_scratch,
          state.lane_count,
          vector_program.slot_count,
          op.input_mask_slot,
          op.input_pdf_slot,
          op.input_cdf_slot,
          op.input_survival_slot,
          &child_mask,
          &child_pdf,
          &child_cdf,
          &child_survival)) {
    return false;
  }
  (void)child_mask;
  (void)child_pdf;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[i]);
    if ((channel_mask & kLeafChannelCdf) != 0U &&
        cdf_out_by_lane != nullptr) {
      cdf_out_by_lane[lane_pos] = child_cdf[lane_pos];
    }
    if ((channel_mask & kLeafChannelSurvival) != 0U &&
        survival_out_by_lane != nullptr) {
      survival_out_by_lane[lane_pos] = child_survival[lane_pos];
    }
  }
  return true;
}

inline bool compiled_math_batch_source_vector_conditioned_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceVectorProgram &vector_program,
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    CompiledMathBatchSourceProductScratch *scratch,
    const double *child_pdf_by_lane,
    const double *child_cdf_by_lane,
    const double *child_survival_by_lane,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane,
    const std::size_t recursion_depth) {
  if (scratch == nullptr || frame == nullptr) {
    return false;
  }
  const auto relation = compiled_math_source_vector_op_relation(op);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !compiled_math_source_product_bounds_have_overlay(
          op.bounds)) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      compiled_math_batch_write_forced_fill_lane(
          relation,
          lanes[i],
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
    }
    return true;
  }
  if (!compiled_math_batch_resolve_source_bounds_for_source(
          program,
          op.source_id,
          op.bounds,
          state,
          lanes,
          lane_count,
          scratch->lower.data(),
          scratch->upper.data(),
          scratch->exact.data(),
          scratch->has_exact.data())) {
    return false;
  }

  std::size_t value_count = 0U;
  std::size_t lower_cdf_count = 0U;
  std::size_t lower_survival_count = 0U;
  std::size_t upper_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower = scratch->lower[lane_pos];
    const double upper = scratch->upper[lane_pos];
    const double current_time = current_time_by_lane[lane_pos];
    if (std::isfinite(upper) && !(upper > lower)) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    if (scratch->has_exact[lane_pos] != 0U) {
      const double exact = scratch->exact[lane_pos];
      const bool exact_in_bounds =
          (!(lower > 0.0) || exact > lower) &&
          (!std::isfinite(upper) || exact < upper);
      if (exact_in_bounds && current_time >= exact) {
        compiled_math_batch_write_certain_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      } else {
        compiled_math_batch_write_impossible_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
      }
      continue;
    }
    if (std::isfinite(upper) && current_time >= upper) {
      compiled_math_batch_write_certain_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    if (current_time <= lower) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    scratch->lanes_a[value_count++] = lane;
    if (lower > 0.0) {
      scratch->times[lane_pos] = lower;
      if (std::isfinite(upper)) {
        scratch->lanes_b[lower_cdf_count++] = lane;
      } else {
        scratch->lanes_c[lower_survival_count++] = lane;
      }
    }
    if (std::isfinite(upper)) {
      scratch->lanes_d[upper_count++] = lane;
    }
  }
  if (value_count == 0U) {
    return true;
  }

  if (lower_cdf_count != 0U &&
      !compiled_math_batch_source_vector_child_cdf_fill(
          program,
          vector_program,
          op,
          state,
          scratch->lanes_b.data(),
          lower_cdf_count,
          scratch->times.data(),
          kLeafChannelCdf,
          frame,
          scratch_depth + 1U,
          scratch->cdf_b.data(),
          nullptr,
          recursion_depth)) {
    return false;
  }
  if (lower_survival_count != 0U &&
      !compiled_math_batch_source_vector_child_cdf_fill(
          program,
          vector_program,
          op,
          state,
          scratch->lanes_c.data(),
          lower_survival_count,
          scratch->times.data(),
          kLeafChannelSurvival,
          frame,
          scratch_depth + 1U,
          nullptr,
          scratch->survival_b.data(),
          recursion_depth)) {
    return false;
  }
  if (upper_count != 0U) {
    for (std::size_t i = 0; i < upper_count; ++i) {
      const auto lane = scratch->lanes_d[i];
      scratch->times[static_cast<std::size_t>(lane)] =
          scratch->upper[static_cast<std::size_t>(lane)];
    }
    if (!compiled_math_batch_source_vector_child_cdf_fill(
            program,
            vector_program,
            op,
            state,
            scratch->lanes_d.data(),
            upper_count,
            scratch->times.data(),
            kLeafChannelCdf,
            frame,
            scratch_depth + 1U,
            scratch->cdf_c.data(),
            nullptr,
            recursion_depth)) {
      return false;
    }
  }

  for (std::size_t i = 0; i < value_count; ++i) {
    const auto lane = scratch->lanes_a[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower = scratch->lower[lane_pos];
    const double upper = scratch->upper[lane_pos];
    mask_by_lane[lane_pos] = fill_mask;
    if (!(lower > 0.0) && !std::isfinite(upper)) {
      if ((fill_mask & kLeafChannelPdf) != 0U) {
        pdf_by_lane[lane_pos] = child_pdf_by_lane[lane_pos];
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        cdf_by_lane[lane_pos] = child_cdf_by_lane[lane_pos];
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        survival_by_lane[lane_pos] = child_survival_by_lane[lane_pos];
      }
      continue;
    }
    if (!std::isfinite(upper)) {
      const double survival_lower =
          lower > 0.0 ? scratch->survival_b[lane_pos] : 1.0;
      if (!std::isfinite(survival_lower) || !(survival_lower > 0.0)) {
        compiled_math_batch_write_impossible_fill_lane(
            lane,
            fill_mask,
            mask_by_lane,
            pdf_by_lane,
            cdf_by_lane,
            survival_by_lane);
        continue;
      }
      if ((fill_mask & kLeafChannelPdf) != 0U) {
        pdf_by_lane[lane_pos] =
            safe_density(child_pdf_by_lane[lane_pos] / survival_lower);
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        cdf_by_lane[lane_pos] =
            clamp_probability(
                (child_cdf_by_lane[lane_pos] + survival_lower - 1.0) /
                survival_lower);
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        survival_by_lane[lane_pos] =
            clamp_probability(
                child_survival_by_lane[lane_pos] / survival_lower);
      }
      continue;
    }
    const double lower_cdf = lower > 0.0 ? scratch->cdf_b[lane_pos] : 0.0;
    const double upper_cdf = scratch->cdf_c[lane_pos];
    const double mass = upper_cdf - lower_cdf;
    if (!std::isfinite(mass) || !(mass > 0.0)) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      pdf_by_lane[lane_pos] =
          safe_density(child_pdf_by_lane[lane_pos] / mass);
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      cdf_by_lane[lane_pos] =
          clamp_probability((child_cdf_by_lane[lane_pos] - lower_cdf) / mass);
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      survival_by_lane[lane_pos] =
          clamp_probability((upper_cdf - child_cdf_by_lane[lane_pos]) / mass);
    }
  }
  return true;
}

inline bool compiled_math_batch_source_vector_pool_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceVectorProgram &vector_program,
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const std::uint8_t fill_mask,
    CompiledMathBatchSourceProductScratch *scratch,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  const auto member_count = static_cast<std::size_t>(op.member_slots.size);
  if (scratch == nullptr || member_count == 0U ||
      op.member_slots.offset == semantic::kInvalidIndex ||
      static_cast<std::size_t>(op.member_slots.offset + op.member_slots.size) >
          program.integral_kernel_source_vector_slot_sets.size()) {
    return false;
  }
  const auto k = static_cast<std::size_t>(op.pool_k);
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const auto max_state = std::min(k, member_count);
  const auto width = max_state + 1U;
  const auto table_idx =
      [lane_stride = state.lane_count, width](const std::size_t row,
                                               const std::size_t col,
                                               const std::size_t lane_pos) {
        return (row * width + col) * lane_stride + lane_pos;
      };
  scratch->ensure_pool_size(state.lane_count, member_count);

  for (std::size_t member = 0U; member < member_count; ++member) {
    const auto &slots =
        program.integral_kernel_source_vector_slot_sets[
            static_cast<std::size_t>(op.member_slots.offset) + member];
    const std::uint8_t *member_mask = nullptr;
    const double *member_pdf = nullptr;
    const double *member_cdf = nullptr;
    const double *member_survival = nullptr;
    if (!compiled_math_source_vector_const_slot_arrays(
            scratch,
            state.lane_count,
            vector_program.slot_count,
            slots.mask,
            slots.pdf,
            slots.cdf,
            slots.survival,
            &member_mask,
            &member_pdf,
            &member_cdf,
            &member_survival)) {
      return false;
    }
    (void)member_mask;
    scratch->pool_member_pdf_ptrs[member] = member_pdf;
    scratch->pool_member_cdf_ptrs[member] = member_cdf;
    scratch->pool_member_survival_ptrs[member] = member_survival;
  }

  auto clear_pool_row = [&](std::vector<double> *table,
                            const std::size_t row) {
    for (std::size_t m = 0U; m < width; ++m) {
      auto *cell = table->data() + table_idx(row, m, 0U);
      for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
        cell[static_cast<std::size_t>(lanes[lane_i])] = 0.0;
      }
    }
  };

  clear_pool_row(&scratch->pool_prefix, 0U);
  if (need_pdf) {
    clear_pool_row(&scratch->pool_suffix, member_count);
  }

  for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[lane_i]);
    scratch->pool_prefix[table_idx(0U, 0U, lane_pos)] = 1.0;
    if (need_pdf) {
      scratch->pool_suffix[table_idx(member_count, 0U, lane_pos)] = 1.0;
    }
    mask_by_lane[lane_pos] = fill_mask;
    if (need_pdf) {
      scratch->pdf_c[lane_pos] = 0.0;
    }
    if (need_cdf) {
      scratch->survival_c[lane_pos] = 0.0;
    }
  }

  for (std::size_t member = 0U; member < member_count; ++member) {
    clear_pool_row(&scratch->pool_prefix, member + 1U);
    const auto *member_cdf = scratch->pool_member_cdf_ptrs[member];
    const auto *member_survival =
        scratch->pool_member_survival_ptrs[member];
    for (std::size_t m = 0U; m <= max_state; ++m) {
      const auto *prefix_prev =
          scratch->pool_prefix.data() + table_idx(member, m, 0U);
      auto *prefix_stay =
          scratch->pool_prefix.data() + table_idx(member + 1U, m, 0U);
      auto *prefix_hit =
          scratch->pool_prefix.data() +
          table_idx(member + 1U, std::min(m + 1U, max_state), 0U);
      for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
        const auto lane_pos = static_cast<std::size_t>(lanes[lane_i]);
        const double prev = prefix_prev[lane_pos];
        prefix_stay[lane_pos] += prev * member_survival[lane_pos];
        prefix_hit[lane_pos] += prev * member_cdf[lane_pos];
      }
    }
  }

  if (need_pdf && k > 0U) {
    for (std::size_t rev = 0U; rev < member_count; ++rev) {
      const std::size_t member = member_count - 1U - rev;
      const std::size_t suffix_member_count = member_count - member - 1U;
      clear_pool_row(&scratch->pool_suffix, member);
      const auto *member_cdf = scratch->pool_member_cdf_ptrs[member];
      const auto *member_survival =
          scratch->pool_member_survival_ptrs[member];
      for (std::size_t m = 0U; m <= std::min(max_state, suffix_member_count);
           ++m) {
        const auto *suffix_prev =
            scratch->pool_suffix.data() + table_idx(member + 1U, m, 0U);
        auto *suffix_stay =
            scratch->pool_suffix.data() + table_idx(member, m, 0U);
        auto *suffix_hit =
            scratch->pool_suffix.data() +
            table_idx(member, std::min(m + 1U, max_state), 0U);
        for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
          const auto lane_pos = static_cast<std::size_t>(lanes[lane_i]);
          const double prev = suffix_prev[lane_pos];
          suffix_stay[lane_pos] += prev * member_survival[lane_pos];
          suffix_hit[lane_pos] += prev * member_cdf[lane_pos];
        }
      }
    }

    for (std::size_t member = 0U; member < member_count; ++member) {
      const auto *member_pdf = scratch->pool_member_pdf_ptrs[member];
      for (std::size_t left = 0U; left < k; ++left) {
        if (left > member || left > max_state) {
          continue;
        }
        const auto right = k - 1U - left;
        if (right > member_count - member - 1U || right > max_state) {
          continue;
        }
        const auto *prefix =
            scratch->pool_prefix.data() + table_idx(member, left, 0U);
        const auto *suffix =
            scratch->pool_suffix.data() + table_idx(member + 1U, right, 0U);
        for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
          const auto lane_pos = static_cast<std::size_t>(lanes[lane_i]);
          scratch->pdf_c[lane_pos] +=
              member_pdf[lane_pos] * prefix[lane_pos] * suffix[lane_pos];
        }
      }
    }
  }

  if (need_pdf) {
    for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
      const auto lane_pos = static_cast<std::size_t>(lanes[lane_i]);
      pdf_by_lane[lane_pos] = safe_density(scratch->pdf_c[lane_pos]);
    }
  }
  if (need_cdf) {
    const auto survival_end = std::min(k, max_state + 1U);
    for (std::size_t m = 0U; m < survival_end; ++m) {
      const auto *prefix =
          scratch->pool_prefix.data() + table_idx(member_count, m, 0U);
      for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
        const auto lane_pos = static_cast<std::size_t>(lanes[lane_i]);
        scratch->survival_c[lane_pos] += prefix[lane_pos];
      }
    }
    for (std::size_t lane_i = 0U; lane_i < lane_count; ++lane_i) {
      const auto lane_pos = static_cast<std::size_t>(lanes[lane_i]);
      const double clamped_survival =
          clamp_probability(scratch->survival_c[lane_pos]);
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        survival_by_lane[lane_pos] = clamped_survival;
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        cdf_by_lane[lane_pos] = clamp_probability(1.0 - clamped_survival);
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_source_vector_onset_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceVectorProgram &vector_program,
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    CompiledMathBatchSourceProductScratch *scratch,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane,
    const std::size_t recursion_depth) {
  if (op.onset_source_id == semantic::kInvalidIndex ||
      op.leaf_index == semantic::kInvalidIndex ||
      frame == nullptr || scratch == nullptr || current_time_by_lane == nullptr) {
    return false;
  }

  std::size_t possible_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double upper =
        current_time_by_lane[lane_pos] - op.leaf_onset_lag;
    if (!(upper > 0.0)) {
      compiled_math_batch_write_impossible_fill_lane(
          lane,
          fill_mask,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane);
      continue;
    }
    scratch->onset_upper[lane_pos] = upper;
    scratch->lanes_a[possible_count++] = lane;
  }
  if (possible_count == 0U) {
    return true;
  }
  if (!compiled_math_batch_resolve_source_bounds_for_source(
          program,
          op.onset_source_id,
          op.onset_bounds,
          state,
          scratch->lanes_a.data(),
          possible_count,
          scratch->lower.data(),
          scratch->upper.data(),
          scratch->exact.data(),
          scratch->has_exact.data())) {
    return false;
  }

  std::size_t exact_count = 0U;
  std::size_t unresolved_count = 0U;
  for (std::size_t i = 0; i < possible_count; ++i) {
    const auto lane = scratch->lanes_a[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    if (scratch->has_exact[lane_pos] == 0U) {
      scratch->lanes_c[unresolved_count++] = lane;
      continue;
    }
    scratch->times[lane_pos] =
        current_time_by_lane[lane_pos] -
        scratch->exact[lane_pos] -
        op.leaf_onset_lag +
        op.leaf_onset_abs_value;
    scratch->lanes_b[exact_count++] = lane;
  }
  if (exact_count != 0U &&
      !compiled_math_batch_source_vector_op_leaf_fill(
          op,
          state,
          scratch->lanes_b.data(),
          exact_count,
          scratch->times.data(),
          fill_mask,
          frame,
          scratch_depth + 1U,
          mask_by_lane,
          pdf_by_lane,
          cdf_by_lane,
          survival_by_lane)) {
    return false;
  }
  if (unresolved_count == 0U) {
    return true;
  }

  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const std::uint8_t shifted_mask =
      static_cast<std::uint8_t>(
          (need_pdf ? kLeafChannelPdf : 0U) |
          (need_cdf ? kLeafChannelCdf : 0U));
  for (std::size_t i = 0; i < unresolved_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(scratch->lanes_c[i]);
    scratch->pdf_c[lane_pos] = 0.0;
    scratch->cdf_c[lane_pos] = 0.0;
  }

  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  const auto q_tile_order =
      std::max<std::size_t>(
          1U,
          std::min<std::size_t>(
              quadrature::kDefaultFiniteOrder,
              kSourceVectorOnsetQuadratureTargetExpandedLanes /
                  unresolved_count));
  for (std::size_t q_begin = 0U; q_begin < quadrature::kDefaultFiniteOrder;
       q_begin += q_tile_order) {
    const auto q_end =
        std::min(q_begin + q_tile_order,
                 static_cast<std::size_t>(quadrature::kDefaultFiniteOrder));
    const auto q_count = q_end - q_begin;
    const auto expanded_count = unresolved_count * q_count;
    scratch->ensure_expanded_lane_state(0U, expanded_count);
    std::size_t child_count = 0U;
    for (std::size_t sample = q_begin; sample < q_end; ++sample) {
      for (std::size_t i = 0U; i < unresolved_count; ++i) {
        const auto parent_lane = scratch->lanes_c[i];
        const auto parent_pos = static_cast<std::size_t>(parent_lane);
        const auto child_lane = static_cast<semantic::Index>(child_count);
        const double scale = 0.5 * scratch->onset_upper[parent_pos];
        const double source_time = scale * (rule.nodes[sample] + 1.0);
        scratch->expanded_lanes[child_count] = child_lane;
        scratch->expanded_parent_lanes[child_count] = parent_lane;
        scratch->expanded_weights[child_count] =
            scale * rule.weights[sample];
        scratch->expanded_source_times[child_count] = source_time;
        scratch->expanded_shifted_times[child_count] =
            current_time_by_lane[parent_pos] -
            source_time -
            op.leaf_onset_lag +
            op.leaf_onset_abs_value;
        ++child_count;
      }
    }

    const CompiledMathLaneBatchState expanded_state{
        expanded_count,
        state.time_slot_count,
        state.time_slots,
        state.workspaces_by_lane,
        state.source_channels_by_lane,
        state.parents_by_lane,
        state.eval_workspaces_by_lane,
        state.used_outcomes_by_lane,
        state.node_values,
        state.node_value_lane_stride,
        state.node_value_lane_map,
        state.leaf_inputs_by_lane,
        scratch->expanded_parent_lanes.data()};
    auto &child_scratch = frame->source_product_scratch_layer(
        scratch_depth + 1U, expanded_count, expanded_count);
    child_scratch.ensure_source_vector_slots(
        expanded_count,
        static_cast<std::size_t>(vector_program.slot_count));
    if (!compiled_math_execute_source_vector_ops(
            program,
            vector_program,
            op.child_ops,
            expanded_state,
            scratch->expanded_lanes.data(),
            expanded_count,
            scratch->expanded_source_times.data(),
            kLeafChannelPdf,
            frame,
            scratch_depth + 1U,
            &child_scratch,
            recursion_depth + 1U)) {
      return false;
    }
    const std::uint8_t *source_mask = nullptr;
    const double *source_pdf = nullptr;
    const double *source_cdf = nullptr;
    const double *source_survival = nullptr;
    if (!compiled_math_source_vector_const_slot_arrays(
            &child_scratch,
            expanded_count,
            vector_program.slot_count,
            op.input_mask_slot,
            op.input_pdf_slot,
            op.input_cdf_slot,
            op.input_survival_slot,
            &source_mask,
            &source_pdf,
            &source_cdf,
            &source_survival)) {
      return false;
    }
    (void)source_mask;
    (void)source_cdf;
    (void)source_survival;
    if (!compiled_math_batch_source_vector_op_leaf_fill_values(
            op,
            expanded_state,
            scratch->expanded_lanes.data(),
            expanded_count,
            scratch->expanded_shifted_times.data(),
            shifted_mask,
            child_scratch.pdf_b.data(),
            child_scratch.cdf_b.data(),
            child_scratch.survival_b.data(),
            &child_scratch)) {
      return false;
    }
    for (std::size_t i = 0U; i < expanded_count; ++i) {
      const auto child_lane = scratch->expanded_lanes[i];
      const auto child_pos = static_cast<std::size_t>(child_lane);
      const auto parent_lane = scratch->expanded_parent_lanes[i];
      const auto parent_pos = static_cast<std::size_t>(parent_lane);
      if (!(source_pdf[child_pos] > 0.0)) {
        continue;
      }
      const double weight =
          scratch->expanded_weights[i] * source_pdf[child_pos];
      if (need_pdf) {
        scratch->pdf_c[parent_pos] += weight * child_scratch.pdf_b[child_pos];
      }
      if (need_cdf) {
        scratch->cdf_c[parent_pos] += weight * child_scratch.cdf_b[child_pos];
      }
    }
  }
  for (std::size_t i = 0U; i < unresolved_count; ++i) {
    const auto lane = scratch->lanes_c[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    mask_by_lane[lane_pos] = fill_mask;
    if (need_pdf) {
      pdf_by_lane[lane_pos] = safe_density(scratch->pdf_c[lane_pos]);
    }
    if (need_cdf) {
      const double cdf = clamp_probability(scratch->cdf_c[lane_pos]);
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        cdf_by_lane[lane_pos] = cdf;
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        survival_by_lane[lane_pos] = clamp_probability(1.0 - cdf);
      }
    }
  }
  return true;
}

inline bool compiled_math_execute_source_vector_ops(
    const CompiledMathProgram &program,
    const CompiledMathSourceVectorProgram &vector_program,
    const CompiledMathIndexSpan ops,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    CompiledMathBatchSourceProductScratch *scratch,
    const std::size_t recursion_depth) {
  if (lane_count == 0U) {
    return true;
  }
  if (lanes == nullptr || current_time_by_lane == nullptr ||
      frame == nullptr || scratch == nullptr || recursion_depth > 32U) {
    return false;
  }
  scratch->ensure_source_vector_slots(
      state.lane_count,
      static_cast<std::size_t>(vector_program.slot_count));
  const auto execution_fill_mask =
      compiled_math_source_vector_execution_fill_mask(
          program,
          ops,
          fill_mask);
  for (semantic::Index op_idx = 0; op_idx < ops.size; ++op_idx) {
    const auto op_pos = static_cast<std::size_t>(ops.offset + op_idx);
    if (op_pos >= program.integral_kernel_source_vector_ops.size()) {
      return false;
    }
    const auto &op = program.integral_kernel_source_vector_ops[op_pos];
    std::uint8_t *op_mask = nullptr;
    double *op_pdf = nullptr;
    double *op_cdf = nullptr;
    double *op_survival = nullptr;
    if (!compiled_math_source_vector_slot_arrays(
            scratch,
            state.lane_count,
            vector_program.slot_count,
            op.output_mask_slot,
            op.output_pdf_slot,
            op.output_cdf_slot,
            op.output_survival_slot,
            &op_mask,
            &op_pdf,
            &op_cdf,
            &op_survival)) {
      return false;
    }
    switch (op.kind) {
    case CompiledMathSourceVectorOpKind::ConstantZero:
      compiled_math_batch_write_fill_for_lanes(
          lanes,
          lane_count,
          execution_fill_mask,
          false,
          op_mask,
          op_pdf,
          op_cdf,
          op_survival);
      break;
    case CompiledMathSourceVectorOpKind::ConstantOne:
      compiled_math_batch_write_fill_for_lanes(
          lanes,
          lane_count,
          execution_fill_mask,
          true,
          op_mask,
          op_pdf,
          op_cdf,
          op_survival);
      break;
    case CompiledMathSourceVectorOpKind::LeafAbsolute:
      if (!compiled_math_batch_source_vector_op_leaf_fill(
              op,
              state,
              lanes,
              lane_count,
              current_time_by_lane,
              execution_fill_mask,
              frame,
              scratch_depth,
              op_mask,
              op_pdf,
              op_cdf,
              op_survival)) {
        return false;
      }
      break;
    case CompiledMathSourceVectorOpKind::ExactGate: {
      const std::uint8_t *child_mask = nullptr;
      const double *child_pdf = nullptr;
      const double *child_cdf = nullptr;
      const double *child_survival = nullptr;
      if (!compiled_math_source_vector_const_slot_arrays(
              scratch,
              state.lane_count,
              vector_program.slot_count,
              op.input_mask_slot,
              op.input_pdf_slot,
              op.input_cdf_slot,
              op.input_survival_slot,
              &child_mask,
              &child_pdf,
              &child_cdf,
              &child_survival) ||
	          !compiled_math_batch_source_vector_exact_gate_fill(
	              program,
	              op,
	              state,
              lanes,
              lane_count,
              current_time_by_lane,
              execution_fill_mask,
              scratch,
              child_mask,
              child_pdf,
              child_cdf,
              child_survival,
              op_mask,
              op_pdf,
              op_cdf,
              op_survival)) {
        return false;
      }
      break;
    }
    case CompiledMathSourceVectorOpKind::Conditioned: {
      const std::uint8_t *child_mask = nullptr;
      const double *child_pdf = nullptr;
      const double *child_cdf = nullptr;
      const double *child_survival = nullptr;
      if (!compiled_math_source_vector_const_slot_arrays(
              scratch,
              state.lane_count,
              vector_program.slot_count,
              op.input_mask_slot,
              op.input_pdf_slot,
              op.input_cdf_slot,
              op.input_survival_slot,
              &child_mask,
              &child_pdf,
              &child_cdf,
              &child_survival)) {
        return false;
      }
      (void)child_mask;
	      if (!compiled_math_batch_source_vector_conditioned_fill(
	              program,
	              vector_program,
	              op,
	              state,
              lanes,
              lane_count,
              current_time_by_lane,
              execution_fill_mask,
              frame,
              scratch_depth,
              scratch,
              child_pdf,
              child_cdf,
              child_survival,
              op_mask,
              op_pdf,
              op_cdf,
              op_survival,
              recursion_depth)) {
        return false;
      }
      break;
    }
    case CompiledMathSourceVectorOpKind::OnsetConvolution:
      if (!compiled_math_batch_source_vector_onset_fill(
              program,
              vector_program,
              op,
              state,
              lanes,
              lane_count,
              current_time_by_lane,
              execution_fill_mask,
              frame,
              scratch_depth,
              scratch,
              op_mask,
              op_pdf,
              op_cdf,
              op_survival,
              recursion_depth)) {
        return false;
      }
      break;
    case CompiledMathSourceVectorOpKind::PoolKOfN:
      if (!compiled_math_batch_source_vector_pool_fill(
              program,
              vector_program,
              op,
              state,
              lanes,
              lane_count,
              execution_fill_mask,
              scratch,
              op_mask,
              op_pdf,
              op_cdf,
              op_survival)) {
        return false;
      }
      break;
    }
  }
  return true;
}

inline bool compiled_math_batch_source_vector_op_leaf_fill(
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  auto *scratch =
      frame == nullptr
          ? nullptr
          : &frame->source_product_scratch_layer(
                scratch_depth,
                state.lane_count,
                lane_count);
  for (std::size_t i = 0; i < lane_count; ++i) {
    mask_by_lane[static_cast<std::size_t>(lanes[i])] = fill_mask;
  }
  return compiled_math_batch_source_vector_op_leaf_fill_values(
      op,
      state,
      lanes,
      lane_count,
      current_time_by_lane,
      fill_mask,
      pdf_by_lane,
      cdf_by_lane,
      survival_by_lane,
      scratch);
}

inline bool compiled_math_batch_source_vector_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_vector_program_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane);

inline bool compiled_math_batch_source_vector_op_leaf_fill_values(
    const CompiledMathSourceVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (static_cast<leaf::DistKind>(op.leaf_dist_kind) ==
      leaf::DistKind::Lognormal) {
    return compiled_math_batch_lognormal_leaf_fill_from_times(
        op.leaf_index,
        op.leaf_onset_abs_value,
        state,
        lanes,
        lane_count,
        current_time_by_lane,
        fill_mask,
        pdf_by_lane,
        cdf_by_lane,
        survival_by_lane,
        scratch);
  }
  if ((fill_mask & kLeafChannelPdf) != 0U &&
      !compiled_math_batch_leaf_values_from_times(
          op.leaf_dist_kind,
          op.leaf_index,
          op.leaf_onset_abs_value,
          state,
          lanes,
          lane_count,
          current_time_by_lane,
          kLeafChannelPdf,
          pdf_by_lane,
          scratch)) {
    return false;
  }
  if ((fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    if (!compiled_math_batch_leaf_values_from_times(
            op.leaf_dist_kind,
            op.leaf_index,
            op.leaf_onset_abs_value,
            state,
            lanes,
            lane_count,
            current_time_by_lane,
            kLeafChannelCdf,
            cdf_by_lane,
            scratch)) {
      return false;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      for (std::size_t i = 0; i < lane_count; ++i) {
        const auto lane_pos = static_cast<std::size_t>(lanes[i]);
        survival_by_lane[lane_pos] =
            clamp_probability(1.0 - cdf_by_lane[lane_pos]);
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_resolve_source_bounds_for_source(
    const CompiledMathProgram &program,
    const semantic::Index source_id,
    const CompiledSourceBoundPlan &bounds,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *lower_by_lane,
    double *upper_by_lane,
    double *exact_by_lane,
    std::uint8_t *has_exact_by_lane) {
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr) {
      return false;
    }
    const auto lane_pos = static_cast<std::size_t>(lane);
    source_channels->source_product_sequence_bound_values(
        source_id,
        bounds,
        &lower_by_lane[lane_pos],
        &upper_by_lane[lane_pos],
        &exact_by_lane[lane_pos],
        &has_exact_by_lane[lane_pos]);
  }
  auto bound_term_time = [&](const CompiledSourceBoundTerm &term,
                             const semantic::Index lane) {
    if (term.time_id != semantic::kInvalidIndex &&
        state.time_slots.values != nullptr &&
        state.has_time(term.time_id, lane)) {
      return state.time(term.time_id, lane);
    }
    auto *source_channels = state.source_channels(lane);
    return source_channels == nullptr
               ? 0.0
               : source_channels->compiled_condition_fallback_time();
  };
  if (bounds.has_condition_lower) {
    for (semantic::Index term_idx = 0; term_idx < bounds.lower.size;
         ++term_idx) {
      const auto &term =
          program.source_condition_bound_terms[
              static_cast<std::size_t>(bounds.lower.offset + term_idx)];
      for (std::size_t i = 0; i < lane_count; ++i) {
        const auto lane = lanes[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        lower_by_lane[lane_pos] =
            std::max(lower_by_lane[lane_pos], bound_term_time(term, lane));
      }
    }
  }
  if (bounds.has_condition_upper) {
    for (semantic::Index term_idx = 0; term_idx < bounds.upper.size;
         ++term_idx) {
      const auto &term =
          program.source_condition_bound_terms[
              static_cast<std::size_t>(bounds.upper.offset + term_idx)];
      for (std::size_t i = 0; i < lane_count; ++i) {
        const auto lane = lanes[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        upper_by_lane[lane_pos] =
            std::min(upper_by_lane[lane_pos], bound_term_time(term, lane));
      }
    }
  }
  if (bounds.has_condition_exact) {
    const auto &term =
        program.source_condition_bound_terms[
            static_cast<std::size_t>(bounds.exact.offset)];
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      exact_by_lane[lane_pos] = bound_term_time(term, lane);
      has_exact_by_lane[lane_pos] = 1U;
    }
  }
  return true;
}

inline bool compiled_math_batch_resolve_source_exact_bounds_for_source(
    const semantic::Index source_id,
    const CompiledSourceBoundPlan &bounds,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    double *exact_by_lane,
    std::uint8_t *has_exact_by_lane) {
  if (exact_by_lane == nullptr || has_exact_by_lane == nullptr) {
    return false;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *source_channels = state.source_channels(lane);
    if (source_channels == nullptr) {
      return false;
    }
    const auto lane_pos = static_cast<std::size_t>(lane);
    source_channels->source_product_exact_bound_value(
        source_id,
        bounds,
        state.workspace(lane),
        &exact_by_lane[lane_pos],
        &has_exact_by_lane[lane_pos]);
  }
  return true;
}

inline bool compiled_math_batch_source_vector_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_vector_program_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t fill_mask,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t scratch_depth,
    std::uint8_t *mask_by_lane,
    double *pdf_by_lane,
    double *cdf_by_lane,
    double *survival_by_lane) {
  if (lane_count == 0U) {
    return true;
  }
  if (source_vector_program_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(source_vector_program_id) >=
          program.integral_kernel_source_vector_programs.size() ||
      lanes == nullptr || current_time_by_lane == nullptr ||
      frame == nullptr || mask_by_lane == nullptr || pdf_by_lane == nullptr ||
      cdf_by_lane == nullptr || survival_by_lane == nullptr) {
    return false;
  }
  const auto &vector_program =
      program.integral_kernel_source_vector_programs[
          static_cast<std::size_t>(source_vector_program_id)];
  auto &scratch = frame->source_product_scratch_layer(
      scratch_depth,
      state.lane_count,
      lane_count);
  scratch.ensure_source_vector_slots(
      state.lane_count,
      static_cast<std::size_t>(vector_program.slot_count));
  if (!compiled_math_execute_source_vector_ops(
          program,
          vector_program,
          vector_program.ops,
          state,
          lanes,
          lane_count,
          current_time_by_lane,
          fill_mask,
          frame,
          scratch_depth,
          &scratch,
          0U)) {
    return false;
  }
  const std::uint8_t *root_mask = nullptr;
  const double *root_pdf = nullptr;
  const double *root_cdf = nullptr;
  const double *root_survival = nullptr;
  if (!compiled_math_source_vector_const_slot_arrays(
          &scratch,
          state.lane_count,
          vector_program.slot_count,
          vector_program.output_mask_slot,
          vector_program.output_pdf_slot,
          vector_program.output_cdf_slot,
          vector_program.output_survival_slot,
          &root_mask,
          &root_pdf,
          &root_cdf,
          &root_survival)) {
    return false;
  }
  (void)root_mask;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane_pos = static_cast<std::size_t>(lanes[i]);
    mask_by_lane[lane_pos] = fill_mask;
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      pdf_by_lane[lane_pos] = root_pdf[lane_pos];
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      cdf_by_lane[lane_pos] = root_cdf[lane_pos];
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      survival_by_lane[lane_pos] = root_survival[lane_pos];
    }
  }
  return true;
}

inline void compiled_math_batch_store_source_vector_fill(
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t lane_stride,
    const semantic::Index cache_slot,
    const semantic::Index lane,
    const std::uint8_t mask,
    const double pdf,
    const double cdf,
    const double survival) {
  const auto pos =
      static_cast<std::size_t>(cache_slot) * lane_stride +
      static_cast<std::size_t>(lane);
  frame->source_vector_valid_mask[pos] |= mask;
  if ((mask & kLeafChannelPdf) != 0U) {
    frame->source_vector_pdf[pos] = pdf;
  }
  if ((mask & kLeafChannelCdf) != 0U) {
    frame->source_vector_cdf[pos] = cdf;
  }
  if ((mask & kLeafChannelSurvival) != 0U) {
    frame->source_vector_survival[pos] = survival;
  }
}

inline double compiled_math_batch_source_vector_cached_value(
    const CompiledMathBatchIntegralFrame &frame,
    const std::size_t lane_stride,
    const semantic::Index cache_slot,
    const semantic::Index lane,
    const std::uint8_t channel_mask) {
  const auto pos =
      static_cast<std::size_t>(cache_slot) * lane_stride +
      static_cast<std::size_t>(lane);
  if (channel_mask == kLeafChannelPdf) {
    return frame.source_vector_pdf[pos];
  }
  if (channel_mask == kLeafChannelCdf) {
    return frame.source_vector_cdf[pos];
  }
  if (channel_mask == kLeafChannelSurvival) {
    return frame.source_vector_survival[pos];
  }
  return 0.0;
}

inline double compiled_math_batch_source_product_array_value(
    const semantic::Index lane,
    const std::uint8_t channel_mask,
    const double *pdf_by_lane,
    const double *cdf_by_lane,
    const double *survival_by_lane) {
  const auto lane_pos = static_cast<std::size_t>(lane);
  if (channel_mask == kLeafChannelPdf) {
    return pdf_by_lane[lane_pos];
  }
  if (channel_mask == kLeafChannelCdf) {
    return cdf_by_lane[lane_pos];
  }
  if (channel_mask == kLeafChannelSurvival) {
    return survival_by_lane[lane_pos];
  }
  return 0.0;
}

inline double compiled_math_batch_forced_channel_value(
    const bool certain,
    const std::uint8_t channel_mask) noexcept {
  if (channel_mask == kLeafChannelPdf) {
    return 0.0;
  }
  if (channel_mask == kLeafChannelCdf) {
    return certain ? 1.0 : 0.0;
  }
  if (channel_mask == kLeafChannelSurvival) {
    return certain ? 0.0 : 1.0;
  }
  return 0.0;
}

inline bool compiled_math_batch_direct_leaf_values(
    const CompiledMathSourceProductChannel &channel,
    const semantic::Index source_product_channel_id,
    const CompiledMathLaneBatchState &state,
    const BatchActiveLaneSpan active,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (channel.leaf_index == semantic::kInvalidIndex) {
    return false;
  }
  (void)source_product_channel_id;
  if (scratch != nullptr) {
    scratch->ensure_size(state.lane_count, active.size);
    for (std::size_t i = 0; i < active.size; ++i) {
      const auto lane = active.lanes[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      scratch->times[lane_pos] =
          batch_source_product_channel_time(channel, state.time_slots, lane);
    }
    return compiled_math_batch_leaf_values_from_times(
        channel.leaf_dist_kind,
        channel.leaf_index,
        channel.leaf_onset_abs_value,
        state,
        active.lanes,
        active.size,
        scratch->times.data(),
        channel_mask,
        values_by_lane,
        scratch);
  }
  return false;
}

inline bool compiled_math_batch_conditioned_direct_leaf_values(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const double *current_time_by_lane,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    CompiledMathBatchSourceProductScratch *scratch) {
  if (lanes == nullptr || current_time_by_lane == nullptr ||
      values_by_lane == nullptr || scratch == nullptr ||
      channel.leaf_index == semantic::kInvalidIndex) {
    return false;
  }
  const auto relation = compiled_math_source_product_channel_relation(channel);
  if (compiled_math_source_product_relation_forces_fill(relation, channel_mask) &&
      !compiled_math_source_product_bounds_have_overlay(
          channel.bounds)) {
    const bool certain =
        relation == ExactRelation::Before || relation == ExactRelation::At;
    for (std::size_t i = 0; i < lane_count; ++i) {
      values_by_lane[static_cast<std::size_t>(lanes[i])] =
          compiled_math_batch_forced_channel_value(certain, channel_mask);
    }
    return true;
  }
  scratch->ensure_size(state.lane_count, lane_count);
  if (!compiled_math_batch_resolve_source_bounds_for_source(
          program,
          channel.source_id,
          channel.bounds,
          state,
          lanes,
          lane_count,
          scratch->lower.data(),
          scratch->upper.data(),
          scratch->exact.data(),
          scratch->has_exact.data())) {
    return false;
  }

  const auto value_mask =
      channel_mask == kLeafChannelPdf ? kLeafChannelPdf : kLeafChannelCdf;
  std::size_t value_count = 0U;
  std::size_t lower_count = 0U;
  std::size_t upper_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower = scratch->lower[lane_pos];
    const double upper = scratch->upper[lane_pos];
    const double current_time = current_time_by_lane[lane_pos];
    if (std::isfinite(upper) && !(upper > lower)) {
      values_by_lane[lane_pos] =
          compiled_math_batch_forced_channel_value(false, channel_mask);
      continue;
    }
    if (scratch->has_exact[lane_pos] != 0U) {
      const double exact = scratch->exact[lane_pos];
      const bool exact_in_bounds =
          (!(lower > 0.0) || exact > lower) &&
          (!std::isfinite(upper) || exact < upper);
      values_by_lane[lane_pos] =
          compiled_math_batch_forced_channel_value(
              exact_in_bounds && current_time >= exact,
              channel_mask);
      continue;
    }
    if (std::isfinite(upper) && current_time >= upper) {
      values_by_lane[lane_pos] =
          compiled_math_batch_forced_channel_value(true, channel_mask);
      continue;
    }
    if (current_time <= lower) {
      values_by_lane[lane_pos] =
          compiled_math_batch_forced_channel_value(false, channel_mask);
      continue;
    }
    scratch->times[lane_pos] = current_time;
    scratch->lanes_a[value_count++] = lane;
    if (lower > 0.0) {
      scratch->lanes_b[lower_count++] = lane;
    }
    if (std::isfinite(upper)) {
      scratch->lanes_c[upper_count++] = lane;
    }
  }
  if (value_count == 0U) {
    return true;
  }
  if (!compiled_math_batch_leaf_values_from_times(
          channel.leaf_dist_kind,
          channel.leaf_index,
          channel.leaf_onset_abs_value,
          state,
          scratch->lanes_a.data(),
          value_count,
          scratch->times.data(),
          value_mask,
          values_by_lane,
          scratch)) {
    return false;
  }
  if (lower_count != 0U) {
    for (std::size_t i = 0; i < lower_count; ++i) {
      const auto lane = scratch->lanes_b[i];
      scratch->times[static_cast<std::size_t>(lane)] =
          scratch->lower[static_cast<std::size_t>(lane)];
    }
    if (!compiled_math_batch_leaf_values_from_times(
            channel.leaf_dist_kind,
            channel.leaf_index,
            channel.leaf_onset_abs_value,
            state,
            scratch->lanes_b.data(),
            lower_count,
            scratch->times.data(),
            kLeafChannelCdf,
            scratch->cdf_b.data(),
            scratch)) {
      return false;
    }
  }
  if (upper_count != 0U) {
    for (std::size_t i = 0; i < upper_count; ++i) {
      const auto lane = scratch->lanes_c[i];
      scratch->times[static_cast<std::size_t>(lane)] =
          scratch->upper[static_cast<std::size_t>(lane)];
    }
    if (!compiled_math_batch_leaf_values_from_times(
            channel.leaf_dist_kind,
            channel.leaf_index,
            channel.leaf_onset_abs_value,
            state,
            scratch->lanes_c.data(),
            upper_count,
            scratch->times.data(),
            kLeafChannelCdf,
            scratch->cdf_c.data(),
            scratch)) {
      return false;
    }
  }
  for (std::size_t i = 0; i < value_count; ++i) {
    const auto lane = scratch->lanes_a[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower = scratch->lower[lane_pos];
    const double upper = scratch->upper[lane_pos];
    const double raw = values_by_lane[lane_pos];
    if (!(lower > 0.0) && !std::isfinite(upper)) {
      values_by_lane[lane_pos] =
          channel_mask == kLeafChannelSurvival
              ? clamp_probability(1.0 - raw)
              : raw;
      continue;
    }
    const double lower_cdf = lower > 0.0 ? scratch->cdf_b[lane_pos] : 0.0;
    const double lower_survival = clamp_probability(1.0 - lower_cdf);
    if (!std::isfinite(upper)) {
      if (!std::isfinite(lower_survival) || !(lower_survival > 0.0)) {
        values_by_lane[lane_pos] =
            compiled_math_batch_forced_channel_value(false, channel_mask);
      } else if (channel_mask == kLeafChannelPdf) {
        values_by_lane[lane_pos] = safe_density(raw / lower_survival);
      } else if (channel_mask == kLeafChannelCdf) {
        values_by_lane[lane_pos] =
            clamp_probability((raw - lower_cdf) / lower_survival);
      } else if (channel_mask == kLeafChannelSurvival) {
        values_by_lane[lane_pos] =
            clamp_probability((1.0 - raw) / lower_survival);
      }
      continue;
    }
    const double upper_cdf = scratch->cdf_c[lane_pos];
    const double mass = upper_cdf - lower_cdf;
    if (!std::isfinite(mass) || !(mass > 0.0)) {
      values_by_lane[lane_pos] =
          compiled_math_batch_forced_channel_value(false, channel_mask);
    } else if (channel_mask == kLeafChannelPdf) {
      values_by_lane[lane_pos] = safe_density(raw / mass);
    } else if (channel_mask == kLeafChannelCdf) {
      values_by_lane[lane_pos] =
          clamp_probability((raw - lower_cdf) / mass);
    } else if (channel_mask == kLeafChannelSurvival) {
      values_by_lane[lane_pos] =
          clamp_probability((upper_cdf - raw) / mass);
    }
  }
  return true;
}

inline bool compiled_math_batch_source_product_planned_values(
    const CompiledMathProgram &program,
    const semantic::Index source_product_channel_id,
    const semantic::Index source_vector_program_id,
    const CompiledMathSourceProductBlockFactorKind factor_kind,
    const std::uint8_t value_channel_mask,
    const std::uint8_t fill_channel_mask,
    const double constant_value,
    const bool cache_result,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    CompiledMathBatchIntegralFrame *frame,
    double *values_by_lane,
    const bool require_source_vector_program = false) {
  (void)require_source_vector_program;
  if (frame == nullptr || values_by_lane == nullptr) {
    return false;
  }
  if (lane_count == 0U) {
    return true;
  }
  if (lanes == nullptr) {
    return false;
  }
  const auto channel_mask = value_channel_mask;
  if (channel_mask == 0U) {
    compiled_math_batch_array_fill(
        lanes,
        lane_count,
        constant_value,
        values_by_lane);
    return true;
  }

  if (source_product_channel_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto channel_pos =
      static_cast<std::size_t>(source_product_channel_id);
  if (channel_pos >= program.integral_kernel_source_product_channels.size()) {
    return false;
  }
  const auto &channel =
      program.integral_kernel_source_product_channels[channel_pos];
  if (!cache_result &&
      factor_kind ==
          CompiledMathSourceProductBlockFactorKind::ConditionedDirectLeaf) {
    auto &scratch =
        frame->source_product_scratch_layer(0U, state.lane_count, lane_count);
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      frame->bind_times[static_cast<std::size_t>(lane)] =
          batch_source_product_channel_time(
              channel,
              state.time_slots,
              lane);
    }
    return compiled_math_batch_conditioned_direct_leaf_values(
        program,
        channel,
        state,
        lanes,
        lane_count,
        frame->bind_times.data(),
        channel_mask,
        values_by_lane,
        &scratch);
  }

  bool direct_leaf_for_all_lanes =
      factor_kind == CompiledMathSourceProductBlockFactorKind::DirectLeaf &&
      channel.leaf_index != semantic::kInvalidIndex;
  if (direct_leaf_for_all_lanes) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      auto *source_channels = state.source_channels(lanes[i]);
      if (source_channels == nullptr) {
        return false;
      }
      if (!source_channels->source_product_direct_leaf_available(
              source_product_channel_id)) {
        direct_leaf_for_all_lanes = false;
        break;
      }
    }
  }

  if (!cache_result && direct_leaf_for_all_lanes) {
    const BatchActiveLaneSpan active{lanes, lane_count};
    return compiled_math_batch_direct_leaf_values(
        channel,
        source_product_channel_id,
        state,
        active,
        channel_mask,
        values_by_lane,
        &frame->source_product_scratch_layer(
            0U,
            state.lane_count,
            lane_count));
  }

  auto &fill_scratch =
      frame->source_product_scratch_layer(
          0U,
          state.lane_count,
          lane_count);
  semantic::Index cache_slot = semantic::kInvalidIndex;
  if (!cache_result) {
    if (source_vector_program_id == semantic::kInvalidIndex) {
      return false;
    }
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      frame->bind_times[static_cast<std::size_t>(lane)] =
          batch_source_product_channel_time(
              channel,
              state.time_slots,
              lane);
    }
    const bool fill_ok =
        compiled_math_batch_source_vector_program_fill(
            program,
            source_vector_program_id,
            state,
            lanes,
            lane_count,
            frame->bind_times.data(),
            channel_mask,
            frame,
            1U,
            fill_scratch.mask_a.data(),
            fill_scratch.pdf_a.data(),
            fill_scratch.cdf_a.data(),
            fill_scratch.survival_a.data());
    if (!fill_ok) {
      return false;
    }
  } else {
    if (source_vector_program_id == semantic::kInvalidIndex) {
      return false;
    }
    cache_slot =
        frame->source_vector_cache_slot(source_vector_program_id);
    if (cache_slot == semantic::kInvalidIndex) {
      return false;
    }
    std::size_t missing_count = 0U;
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const auto program_pos =
          static_cast<std::size_t>(cache_slot) * state.lane_count +
          lane_pos;
      if ((frame->source_vector_valid_mask[program_pos] &
           channel_mask) != 0U) {
        continue;
      }
      frame->bind_times[lane_pos] =
          batch_source_product_channel_time(
              channel,
              state.time_slots,
              lane);
      fill_scratch.lanes_a[missing_count++] = lane;
    }
    if (missing_count != 0U) {
      const bool fill_ok =
          compiled_math_batch_source_vector_program_fill(
              program,
              source_vector_program_id,
              state,
              fill_scratch.lanes_a.data(),
              missing_count,
              frame->bind_times.data(),
              fill_channel_mask,
              frame,
              1U,
              fill_scratch.mask_a.data(),
              fill_scratch.pdf_a.data(),
              fill_scratch.cdf_a.data(),
              fill_scratch.survival_a.data());
      if (!fill_ok) {
        return false;
      }
      for (std::size_t i = 0; i < missing_count; ++i) {
        const auto lane = fill_scratch.lanes_a[i];
        const auto lane_pos = static_cast<std::size_t>(lane);
        compiled_math_batch_store_source_vector_fill(
            frame,
            state.lane_count,
            cache_slot,
            lane,
            fill_scratch.mask_a[lane_pos],
            fill_scratch.pdf_a[lane_pos],
            fill_scratch.cdf_a[lane_pos],
            fill_scratch.survival_a[lane_pos]);
      }
    }
  }

  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    values_by_lane[lane_pos] =
        cache_result
            ? compiled_math_batch_source_vector_cached_value(
                  *frame,
                  state.lane_count,
                  cache_slot,
                  lane,
                  channel_mask)
            : compiled_math_batch_source_product_array_value(
                  lane,
                  channel_mask,
                  fill_scratch.pdf_a.data(),
                  fill_scratch.cdf_a.data(),
                  fill_scratch.survival_a.data());
  }
  return true;
}

inline bool compiled_math_batch_source_product_op_values(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    CompiledMathBatchIntegralFrame *frame,
    double *values_by_lane) {
  return compiled_math_batch_source_product_planned_values(
      program,
      op.source_product_channel_id,
      op.source_vector_program_id,
      op.factor_kind,
      op.value_channel_mask,
      op.fill_channel_mask,
      op.constant_value,
      op.cache_result,
      state,
      lanes,
      lane_count,
      frame,
      values_by_lane);
}

inline bool compiled_math_batch_source_product_value_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    const CompiledMathLaneBatchState &state,
    const BatchActiveLaneSpan active,
    CompiledMathBatchIntegralFrame *frame,
    double *out_by_lane) {
  if (frame == nullptr || out_by_lane == nullptr) {
    return false;
  }
  compiled_math_batch_array_fill(
      active.lanes,
      active.size,
      1.0,
      out_by_lane);
  const auto *term_lanes = active.lanes;
  const std::size_t term_count = active.size;
  std::size_t live_count = term_count;

  for (semantic::Index op_idx = 0; op_idx < source_product_ops.size;
       ++op_idx) {
    if (live_count == 0U) {
      break;
    }
    const auto op_pos =
        static_cast<std::size_t>(source_product_ops.offset + op_idx);
    if (op_pos >= program.integral_kernel_source_product_ops.size()) {
      return false;
    }
    const auto &op = program.integral_kernel_source_product_ops[op_pos];
    if (op.value_channel_mask == 0U) {
      live_count =
          compiled_math_batch_array_multiply_scalar_mask(
              term_lanes,
              term_count,
              op.constant_value,
              out_by_lane);
      continue;
    }

    const semantic::Index *eval_lanes = term_lanes;
    std::size_t eval_count = term_count;
    if (live_count != term_count) {
      eval_count =
          compiled_math_batch_array_live_lanes_from_mask(
              term_lanes,
              term_count,
              out_by_lane,
              frame->source_lanes_a.data());
      eval_lanes = frame->source_lanes_a.data();
    }
    if (eval_count == 0U) {
      live_count = 0U;
      break;
    }
    if (!compiled_math_batch_source_product_planned_values(
            program,
            op.source_product_channel_id,
            op.source_vector_program_id,
            op.factor_kind,
            op.value_channel_mask,
            op.fill_channel_mask,
            op.constant_value,
            op.cache_result,
            state,
            eval_lanes,
            eval_count,
            frame,
            frame->source_leaf_values.data())) {
      return false;
    }
    live_count =
        compiled_math_batch_array_multiply_values_mask(
            eval_lanes,
            eval_count,
            frame->source_leaf_values.data(),
            out_by_lane);
  }
  return true;
}

inline bool compiled_math_batch_outcome_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan outcome_gate_nodes,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane) {
  const auto *used_outcomes = state.used_outcomes(lane);
  if (used_outcomes == nullptr) {
    return true;
  }
  auto *parent = state.parent(lane);
  if (parent == nullptr) {
    return false;
  }
  for (semantic::Index i = 0; i < outcome_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_outcome_gate_nodes[
            static_cast<std::size_t>(outcome_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (gate_node.kind != CompiledMathNodeKind::OutcomeSubsetUnused) {
      return false;
    }
    const auto offset = static_cast<std::size_t>(gate_node.subject_id);
    const auto size = static_cast<std::size_t>(gate_node.aux_id);
    const auto &indices = parent->plan().compiled_outcome_gate_indices;
    if (offset + size > indices.size()) {
      return false;
    }
    for (std::size_t j = 0; j < size; ++j) {
      const auto outcome_idx = indices[offset + j];
      if (outcome_idx != semantic::kInvalidIndex &&
          static_cast<std::size_t>(outcome_idx) < used_outcomes->size() &&
          (*used_outcomes)[static_cast<std::size_t>(outcome_idx)] != 0U) {
        return false;
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_time_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan time_gate_nodes,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane) {
  for (semantic::Index i = 0; i < time_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_time_gate_nodes[
            static_cast<std::size_t>(time_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (gate_node.kind != CompiledMathNodeKind::TimeGate) {
      return false;
    }
    if (state.time(gate_node.time_id, lane) <
        state.time(gate_node.aux_id, lane)) {
      return false;
    }
  }
  return true;
}

inline bool compiled_math_batch_upper_bounds_from_terms(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::uint8_t *found_by_lane,
    double *time_by_lane,
    double *normalizer_by_lane) {
  if (lanes == nullptr || found_by_lane == nullptr ||
      time_by_lane == nullptr || normalizer_by_lane == nullptr) {
    return false;
  }
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const auto lane = lanes[lane_pos];
    const auto lane_index = static_cast<std::size_t>(lane);
    found_by_lane[lane_index] = 0U;
    time_by_lane[lane_index] = std::numeric_limits<double>::infinity();
    normalizer_by_lane[lane_index] = 0.0;
  }
  if (node.aux_id == semantic::kInvalidIndex ||
      node.aux2_id == semantic::kInvalidIndex ||
      node.aux2_id == 0) {
    return true;
  }
  const auto offset = static_cast<std::size_t>(node.aux_id);
  const auto size = static_cast<std::size_t>(node.aux2_id);
  if (offset + size > program.timed_upper_bound_terms.size()) {
    return false;
  }
  for (std::size_t i = 0; i < size; ++i) {
    const auto &term = program.timed_upper_bound_terms[offset + i];
    if (term.normalizer_node_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(term.normalizer_node_id) >=
            program.nodes.size()) {
      return false;
    }
    for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
      const auto lane = lanes[lane_pos];
      double fact_time = state.time(node.time_id, lane);
      if (term.time_id != semantic::kInvalidIndex &&
          static_cast<std::size_t>(term.time_id) < state.time_slot_count &&
          state.has_time(term.time_id, lane)) {
        fact_time = state.time(term.time_id, lane);
      }
      const double normalizer =
          state.node_value(term.normalizer_node_id, lane);
      const auto lane_index = static_cast<std::size_t>(lane);
      if (std::isfinite(fact_time) && normalizer > 0.0 &&
          (found_by_lane[lane_index] == 0U ||
           fact_time < time_by_lane[lane_index])) {
        found_by_lane[lane_index] = 1U;
        time_by_lane[lane_index] = fact_time;
        normalizer_by_lane[lane_index] = normalizer;
      }
    }
  }
  return true;
}

inline bool compiled_math_batch_expr_upper_bounds_for_node(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    CompiledSourceView *const *condition_evaluators_by_lane,
    std::uint8_t *found_by_lane,
    double *time_by_lane,
    double *normalizer_by_lane) {
  if (!compiled_math_batch_upper_bounds_from_terms(
          program,
          node,
          state,
          lanes,
          lane_count,
          found_by_lane,
          time_by_lane,
          normalizer_by_lane)) {
    return false;
  }
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const auto lane = lanes[lane_pos];
    const auto lane_index = static_cast<std::size_t>(lane);
    auto *condition_evaluator =
        condition_evaluators_by_lane == nullptr
            ? nullptr
            : condition_evaluators_by_lane[lane_index];
    ExactTimedExprUpperBound sequence_upper;
    if (condition_evaluator != nullptr &&
        condition_evaluator->expr_upper_bound_for(
            node.subject_id,
            &sequence_upper) &&
        std::isfinite(sequence_upper.time) &&
        sequence_upper.normalizer > 0.0 &&
        (found_by_lane[lane_index] == 0U ||
         sequence_upper.time < time_by_lane[lane_index])) {
      found_by_lane[lane_index] = 1U;
      time_by_lane[lane_index] = sequence_upper.time;
      normalizer_by_lane[lane_index] = sequence_upper.normalizer;
    }
  }
  return true;
}

inline bool compiled_math_batch_node_evaluators_for_lanes(
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::vector<CompiledSourceView *> *evaluators_by_lane) {
  if (lanes == nullptr || evaluators_by_lane == nullptr) {
    return false;
  }
  evaluators_by_lane->assign(state.lane_count, nullptr);
  if (node.source_view_id == 0 ||
      node.source_view_id == semantic::kInvalidIndex) {
    for (std::size_t i = 0; i < lane_count; ++i) {
      const auto lane = lanes[i];
      (*evaluators_by_lane)[static_cast<std::size_t>(lane)] =
          state.parent(lane);
    }
    return true;
  }
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto lane = lanes[i];
    auto *eval_workspace = state.eval_workspace(lane);
    if (eval_workspace == nullptr) {
      return false;
    }
    (*evaluators_by_lane)[static_cast<std::size_t>(lane)] =
        eval_workspace->source_view_evaluator(
            node.source_view_id,
            state.parent(lane));
  }
  return true;
}

inline bool compiled_math_batch_node_evaluators_for_lanes(
    const CompiledMathNode &node,
    const std::vector<CompiledMathBatchLane> &lanes,
    std::vector<CompiledSourceView *> *evaluators_by_lane) {
  if (evaluators_by_lane == nullptr) {
    return false;
  }
  evaluators_by_lane->assign(lanes.size(), nullptr);
  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    auto *condition_evaluator = lanes[lane].parent;
    if (node.source_view_id != 0 &&
        node.source_view_id != semantic::kInvalidIndex) {
      if (lanes[lane].eval_workspace == nullptr) {
        return false;
      }
      condition_evaluator =
          lanes[lane].eval_workspace->source_view_evaluator(
              node.source_view_id,
              lanes[lane].parent);
    }
    (*evaluators_by_lane)[lane] = condition_evaluator;
  }
  return true;
}

inline bool compiled_math_batch_expr_upper_bounds_for_op(
    const CompiledMathProgram &program,
    const CompiledMathLaneTermOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::uint8_t *found_by_lane,
    double *time_by_lane,
    double *normalizer_by_lane) {
  if (lanes == nullptr || found_by_lane == nullptr ||
      time_by_lane == nullptr || normalizer_by_lane == nullptr ||
      op.expr_upper_time_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(op.expr_upper_time_id) >=
          state.time_slot_count) {
    return false;
  }
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const auto lane = lanes[lane_pos];
    const auto lane_index = static_cast<std::size_t>(lane);
    found_by_lane[lane_index] = 0U;
    time_by_lane[lane_index] = std::numeric_limits<double>::infinity();
    normalizer_by_lane[lane_index] = 0.0;
  }

  const auto upper_offset =
      static_cast<std::size_t>(op.expr_upper_timed_upper_bound_terms.offset);
  const auto upper_size =
      static_cast<std::size_t>(op.expr_upper_timed_upper_bound_terms.size);
  if (upper_offset + upper_size > program.timed_upper_bound_terms.size()) {
    return false;
  }
  for (std::size_t i = 0; i < upper_size; ++i) {
    const auto &term = program.timed_upper_bound_terms[upper_offset + i];
    if (term.normalizer_node_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(term.normalizer_node_id) >=
            program.nodes.size()) {
      return false;
    }
    for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
      const auto lane = lanes[lane_pos];
      double fact_time = state.time(op.expr_upper_time_id, lane);
      if (term.time_id != semantic::kInvalidIndex &&
          static_cast<std::size_t>(term.time_id) < state.time_slot_count &&
          state.has_time(term.time_id, lane)) {
        fact_time = state.time(term.time_id, lane);
      }
      const double normalizer =
          state.node_value(term.normalizer_node_id, lane);
      const auto lane_index = static_cast<std::size_t>(lane);
      if (std::isfinite(fact_time) && normalizer > 0.0 &&
          (found_by_lane[lane_index] == 0U ||
           fact_time < time_by_lane[lane_index])) {
        found_by_lane[lane_index] = 1U;
        time_by_lane[lane_index] = fact_time;
        normalizer_by_lane[lane_index] = normalizer;
      }
    }
  }

  const bool parent_source_view =
      op.expr_upper_source_view_id == 0 ||
      op.expr_upper_source_view_id == semantic::kInvalidIndex;
  for (std::size_t lane_pos = 0; lane_pos < lane_count; ++lane_pos) {
    const auto lane = lanes[lane_pos];
    auto *condition_evaluator = parent_source_view ? state.parent(lane) : nullptr;
    if (!parent_source_view) {
      auto *eval_workspace = state.eval_workspace(lane);
      if (eval_workspace == nullptr) {
        return false;
      }
      condition_evaluator =
          eval_workspace->source_view_evaluator(
              op.expr_upper_source_view_id,
              state.parent(lane));
    }
    ExactTimedExprUpperBound sequence_upper;
    const auto lane_index = static_cast<std::size_t>(lane);
    if (condition_evaluator != nullptr &&
        condition_evaluator->expr_upper_bound_for(
            op.expr_upper_subject_id,
            &sequence_upper) &&
        std::isfinite(sequence_upper.time) &&
        sequence_upper.normalizer > 0.0 &&
        (found_by_lane[lane_index] == 0U ||
         sequence_upper.time < time_by_lane[lane_index])) {
      found_by_lane[lane_index] = 1U;
      time_by_lane[lane_index] = sequence_upper.time;
      normalizer_by_lane[lane_index] = sequence_upper.normalizer;
    }
  }
  return true;
}

inline bool compiled_math_batch_apply_expr_upper_op_to_lanes(
    const CompiledMathProgram &program,
    const CompiledMathLaneTermOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::uint8_t *upper_found_by_lane,
    double *upper_time_by_lane,
    double *upper_normalizer_by_lane,
    double *products_by_lane,
    semantic::Index *next_lanes,
    std::size_t *next_count) {
  if (products_by_lane == nullptr || next_lanes == nullptr ||
      next_count == nullptr ||
      op.expr_upper_node_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(op.expr_upper_node_id) >= program.nodes.size()) {
    return false;
  }
  if (!compiled_math_batch_expr_upper_bounds_for_op(
          program,
          op,
          state,
          lanes,
          lane_count,
          upper_found_by_lane,
          upper_time_by_lane,
          upper_normalizer_by_lane)) {
    return false;
  }
  *next_count =
      compiled_math_batch_array_apply_expr_upper_op_compact(
          lanes,
          lane_count,
          op,
          state,
          upper_found_by_lane,
          upper_time_by_lane,
          upper_normalizer_by_lane,
          products_by_lane,
          next_lanes);
  return true;
}

inline bool compiled_math_direct_tile_time_for_id(
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathBatchIntegralFrame &frame,
    const std::size_t parent_pos,
    const semantic::Index parent_lane,
    const std::size_t q_abs,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const semantic::Index time_id,
    double *out) {
  if (out == nullptr || time_id == semantic::kInvalidIndex) {
    return false;
  }
  if (time_id == kernel.bind_time_id) {
    *out = frame.direct_parent_shift[parent_pos] +
           frame.direct_parent_scale[parent_pos] * rule.nodes[q_abs];
    return true;
  }
  if (static_cast<std::size_t>(time_id) >= parent_state.time_slot_count) {
    return false;
  }
  *out = parent_state.time_slots.get(time_id, parent_lane);
  return true;
}

inline bool compiled_math_direct_tile_time_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathIndexSpan time_gate_nodes,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathBatchIntegralFrame &frame,
    const std::size_t parent_pos,
    const semantic::Index parent_lane,
    const std::size_t q_abs,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    bool *open) {
  if (open == nullptr) {
    return false;
  }
  *open = true;
  for (semantic::Index i = 0; i < time_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_time_gate_nodes[
            static_cast<std::size_t>(time_gate_nodes.offset + i)];
    if (node_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(node_id) >= program.nodes.size()) {
      return false;
    }
    const auto &gate_node = program.nodes[static_cast<std::size_t>(node_id)];
    if (gate_node.kind != CompiledMathNodeKind::TimeGate) {
      return false;
    }
    double lhs = 0.0;
    double rhs = 0.0;
    if (!compiled_math_direct_tile_time_for_id(
            kernel,
            parent_state,
            frame,
            parent_pos,
            parent_lane,
            q_abs,
            rule,
            gate_node.time_id,
            &lhs) ||
        !compiled_math_direct_tile_time_for_id(
            kernel,
            parent_state,
            frame,
            parent_pos,
            parent_lane,
            q_abs,
            rule,
            gate_node.aux_id,
            &rhs)) {
      return false;
    }
    if (lhs < rhs) {
      *open = false;
      return true;
    }
  }
  return true;
}

inline bool compiled_math_direct_tile_expr_upper_bounds_for_op(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneTermOp &op,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathBatchIntegralFrame &frame,
    const std::size_t parent_pos,
    const semantic::Index parent_lane,
    const std::size_t q_abs,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    bool *found,
    double *upper_time,
    double *normalizer) {
  if (found == nullptr || upper_time == nullptr || normalizer == nullptr) {
    return false;
  }
  *found = false;
  *upper_time = std::numeric_limits<double>::infinity();
  *normalizer = 0.0;

  const auto upper_offset =
      static_cast<std::size_t>(op.expr_upper_timed_upper_bound_terms.offset);
  const auto upper_size =
      static_cast<std::size_t>(op.expr_upper_timed_upper_bound_terms.size);
  if (upper_offset + upper_size > program.timed_upper_bound_terms.size()) {
    return false;
  }
  double base_fact_time = 0.0;
  if (upper_size != 0U &&
      !compiled_math_direct_tile_time_for_id(
          kernel,
          parent_state,
          frame,
          parent_pos,
          parent_lane,
          q_abs,
          rule,
          op.expr_upper_time_id,
          &base_fact_time)) {
    return false;
  }
  for (std::size_t i = 0; i < upper_size; ++i) {
    const auto &term = program.timed_upper_bound_terms[upper_offset + i];
    if (term.normalizer_node_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(term.normalizer_node_id) >=
            program.nodes.size()) {
      return false;
    }
    double fact_time = base_fact_time;
    if (term.time_id != semantic::kInvalidIndex &&
        !compiled_math_direct_tile_time_for_id(
            kernel,
            parent_state,
            frame,
            parent_pos,
            parent_lane,
            q_abs,
            rule,
            term.time_id,
            &fact_time)) {
      return false;
    }
    const double term_normalizer =
        parent_state.node_value(term.normalizer_node_id, parent_lane);
    if (std::isfinite(fact_time) && term_normalizer > 0.0 &&
        (!*found || fact_time < *upper_time)) {
      *found = true;
      *upper_time = fact_time;
      *normalizer = term_normalizer;
    }
  }

  CompiledSourceView *condition_evaluator = nullptr;
  if (op.expr_upper_source_view_id == 0 ||
      op.expr_upper_source_view_id == semantic::kInvalidIndex) {
    condition_evaluator = parent_state.parent(parent_lane);
  } else {
    auto *eval_workspace = parent_state.eval_workspace(parent_lane);
    if (eval_workspace == nullptr) {
      return false;
    }
    condition_evaluator =
        eval_workspace->source_view_evaluator(
            op.expr_upper_source_view_id,
            parent_state.parent(parent_lane));
  }
  ExactTimedExprUpperBound sequence_upper;
  if (condition_evaluator != nullptr &&
      condition_evaluator->expr_upper_bound_for(
          op.expr_upper_subject_id,
          &sequence_upper) &&
      std::isfinite(sequence_upper.time) &&
      sequence_upper.normalizer > 0.0 &&
      (!*found || sequence_upper.time < *upper_time)) {
    *found = true;
    *upper_time = sequence_upper.time;
    *normalizer = sequence_upper.normalizer;
  }
  return true;
}

inline bool compiled_math_direct_tile_apply_expr_upper_op(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneTermOp &op,
    const CompiledMathLaneBatchState &parent_state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    const std::size_t child_width,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    CompiledMathBatchIntegralFrame *frame,
    semantic::Index *next_lanes,
    std::size_t *next_count) {
  if (lanes == nullptr || frame == nullptr || next_lanes == nullptr ||
      next_count == nullptr ||
      op.expr_upper_node_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(op.expr_upper_node_id) >= program.nodes.size()) {
    return false;
  }
  std::size_t out_count = 0U;
  for (std::size_t i = 0; i < lane_count; ++i) {
    const auto child_pos = static_cast<std::size_t>(lanes[i]);
    const auto parent_pos = child_pos / child_width;
    if (parent_pos >= eligible_count) {
      return false;
    }
    const auto parent_lane = frame->parent_lanes[parent_pos];
    const auto q_abs = q_begin + child_pos - parent_pos * child_width;
    double current_time = 0.0;
    if (!compiled_math_direct_tile_time_for_id(
            kernel,
            parent_state,
            *frame,
            parent_pos,
            parent_lane,
            q_abs,
            rule,
            op.expr_upper_time_id,
            &current_time)) {
      return false;
    }
    bool has_upper = false;
    double upper_time = std::numeric_limits<double>::infinity();
    double normalizer = 0.0;
    if (!compiled_math_direct_tile_expr_upper_bounds_for_op(
            program,
            kernel,
            op,
            parent_state,
            *frame,
            parent_pos,
            parent_lane,
            q_abs,
            rule,
            &has_upper,
            &upper_time,
            &normalizer)) {
      return false;
    }
    double &product = frame->products[child_pos];
    if (op.expr_upper_mode == CompiledMathIntegralExprUpperMode::AfterOne) {
      if (!has_upper || current_time < upper_time) {
        product = 0.0;
        continue;
      }
    } else if (has_upper) {
      if (current_time >= upper_time) {
        product = 0.0;
        continue;
      }
      product /= normalizer;
      if (!compiled_math_batch_product_is_live(product)) {
        product = 0.0;
        continue;
      }
    }
    if (compiled_math_batch_product_is_live(product)) {
      next_lanes[out_count++] = static_cast<semantic::Index>(child_pos);
    }
  }
  *next_count = out_count;
  return true;
}

inline void compiled_math_direct_tile_multiply_value_at(
    const std::size_t child_pos,
    const double value,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *live_count,
    double *factor_values_out = nullptr) {
  if (factor_values_out != nullptr) {
    if (compiled_math_batch_product_is_live(value)) {
      factor_values_out[child_pos] = value;
      ++(*live_count);
    } else {
      factor_values_out[child_pos] = 0.0;
    }
    return;
  }
  if (compiled_math_batch_multiply_product_value_at(
          child_pos,
          value,
          frame->products.data())) {
    ++(*live_count);
  }
}

inline void compiled_math_direct_tile_multiply_finished_value_at(
    const std::size_t child_pos,
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t channel_mask,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *live_count,
    double *factor_values_out = nullptr) {
  compiled_math_direct_tile_multiply_value_at(
      child_pos,
      batch_source_product_finish_base_value(
          base_pdf,
          base_cdf,
          q,
          channel_mask),
      frame,
      live_count,
      factor_values_out);
}

inline void compiled_math_direct_tile_multiply_live_value_at(
    const std::size_t child_pos,
    const double value,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *live_count,
    double *factor_values_out = nullptr) {
  if (factor_values_out != nullptr) {
    if (compiled_math_batch_product_is_live(value)) {
      factor_values_out[child_pos] = value;
      ++(*live_count);
    } else {
      factor_values_out[child_pos] = 0.0;
    }
    return;
  }
  const double product = frame->products[child_pos] * value;
  if (compiled_math_batch_product_is_live(product)) {
    frame->products[child_pos] = product;
    ++(*live_count);
  } else {
    frame->products[child_pos] = 0.0;
  }
}

inline void compiled_math_direct_tile_multiply_live_finished_value_at(
    const std::size_t child_pos,
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t channel_mask,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *live_count,
    double *factor_values_out = nullptr) {
  compiled_math_direct_tile_multiply_live_value_at(
      child_pos,
      batch_source_product_finish_base_value(
          base_pdf,
          base_cdf,
          q,
          channel_mask),
      frame,
      live_count,
      factor_values_out);
}

inline void compiled_math_direct_tile_multiply_live_pdf_at(
    const std::size_t child_pos,
    const double base_pdf,
    const double q,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *live_count,
    double *factor_values_out = nullptr) {
  compiled_math_direct_tile_multiply_live_value_at(
      child_pos,
      (1.0 - q) * safe_density(base_pdf),
      frame,
      live_count,
      factor_values_out);
}

inline void compiled_math_direct_tile_multiply_live_cdf_at(
    const std::size_t child_pos,
    const double base_cdf,
    const double q,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *live_count,
    double *factor_values_out = nullptr) {
  compiled_math_direct_tile_multiply_live_value_at(
      child_pos,
      clamp_probability((1.0 - q) * base_cdf),
      frame,
      live_count,
      factor_values_out);
}

inline void compiled_math_direct_tile_multiply_live_survival_at(
    const std::size_t child_pos,
    const double base_cdf,
    const double q,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *live_count,
    double *factor_values_out = nullptr) {
  const double cdf = clamp_probability((1.0 - q) * base_cdf);
  compiled_math_direct_tile_multiply_live_value_at(
      child_pos,
      clamp_probability(1.0 - cdf),
      frame,
      live_count,
      factor_values_out);
}

inline double compiled_math_direct_tile_forced_channel_value(
    const bool certain,
    const std::uint8_t channel_mask) noexcept {
  if (channel_mask == kLeafChannelPdf) {
    return 0.0;
  }
  if (channel_mask == kLeafChannelCdf) {
    return certain ? 1.0 : 0.0;
  }
  if (channel_mask == kLeafChannelSurvival) {
    return certain ? 0.0 : 1.0;
  }
  return 0.0;
}

inline double compiled_math_direct_tile_conditioned_continuous_value(
    const std::size_t op_parent_pos,
    const double raw_value,
    const std::uint8_t channel_mask,
    const CompiledMathBatchIntegralFrame &frame) noexcept {
  const double lower = frame.direct_bound_lower_by_op_parent[op_parent_pos];
  const double upper = frame.direct_bound_upper_by_op_parent[op_parent_pos];
  if (!(lower > 0.0) && !std::isfinite(upper)) {
    return channel_mask == kLeafChannelSurvival
               ? clamp_probability(1.0 - raw_value)
               : raw_value;
  }
  const double lower_cdf = frame.direct_lower_cdf_by_op_parent[op_parent_pos];
  const double lower_survival =
      frame.direct_lower_survival_by_op_parent[op_parent_pos];
  if (!std::isfinite(upper)) {
    if (!std::isfinite(lower_survival) || !(lower_survival > 0.0)) {
      return compiled_math_direct_tile_forced_channel_value(false, channel_mask);
    }
    if (channel_mask == kLeafChannelPdf) {
      return safe_density(raw_value / lower_survival);
    }
    if (channel_mask == kLeafChannelCdf) {
      return clamp_probability((raw_value - lower_cdf) / lower_survival);
    }
    if (channel_mask == kLeafChannelSurvival) {
      return clamp_probability((1.0 - raw_value) / lower_survival);
    }
    return 0.0;
  }
  const double upper_cdf = frame.direct_upper_cdf_by_op_parent[op_parent_pos];
  const double mass = upper_cdf - lower_cdf;
  if (!std::isfinite(mass) || !(mass > 0.0)) {
    return compiled_math_direct_tile_forced_channel_value(false, channel_mask);
  }
  if (channel_mask == kLeafChannelPdf) {
    return safe_density(raw_value / mass);
  }
  if (channel_mask == kLeafChannelCdf) {
    return clamp_probability((raw_value - lower_cdf) / mass);
  }
  if (channel_mask == kLeafChannelSurvival) {
    return clamp_probability((upper_cdf - raw_value) / mass);
  }
  return 0.0;
}

inline std::size_t compiled_math_direct_tile_child_pos(
    const semantic::Index *eval_lanes,
    const std::size_t eval_pos) noexcept {
  return eval_lanes == nullptr
             ? eval_pos
             : static_cast<std::size_t>(eval_lanes[eval_pos]);
}

inline double compiled_math_direct_tile_bind_time(
    const CompiledMathBatchIntegralFrame &frame,
    const std::size_t parent_pos,
    const std::size_t q_abs,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule) noexcept {
  return frame.direct_parent_shift[parent_pos] +
         frame.direct_parent_scale[parent_pos] * rule.nodes[q_abs];
}

inline bool compiled_math_bound_vector_plan_span(
    const CompiledMathProgram &program,
    const semantic::Index plan_id,
    CompiledMathIndexSpan *span_out) {
  if (span_out == nullptr ||
      plan_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(plan_id) >= program.bound_vector_plans.size()) {
    return false;
  }
  const auto span =
      program.bound_vector_plans[static_cast<std::size_t>(plan_id)].ops;
  if (static_cast<std::size_t>(span.offset + span.size) >
      program.bound_vector_ops.size()) {
    return false;
  }
  *span_out = span;
  return true;
}

struct CompiledMathBoundVectorEvalContext {
  const CompiledMathLaneBatchState *state{nullptr};
  const semantic::Index *lanes{nullptr};
  std::size_t lane_count{0U};
  const CompiledMathBatchIntegralFrame *frame{nullptr};
  std::size_t op_parent_offset{0U};
  std::size_t eligible_count{0U};
  std::size_t q_begin{0U};
  const quadrature::Rule<quadrature::kDefaultFiniteOrder> *rule{nullptr};
  std::size_t child_width{0U};
};

inline bool compiled_math_bound_vector_direct_child_parent(
    const CompiledMathBoundVectorEvalContext &context,
    const std::size_t i,
    std::size_t *child_pos_out,
    std::size_t *parent_pos_out) {
  if (context.frame == nullptr || context.rule == nullptr ||
      context.child_width == 0U || context.lanes == nullptr ||
      child_pos_out == nullptr || parent_pos_out == nullptr) {
    return false;
  }
  const auto child_pos = static_cast<std::size_t>(context.lanes[i]);
  const auto parent_pos = child_pos / context.child_width;
  if (parent_pos >= context.eligible_count) {
    return false;
  }
  *child_pos_out = child_pos;
  *parent_pos_out = parent_pos;
  return true;
}

inline bool compiled_math_eval_scalar_bound_vector_plan_compact(
    const CompiledMathProgram &program,
    const semantic::Index plan_id,
    const CompiledMathBoundVectorEvalContext &context,
    double *values) {
  if ((context.lane_count != 0U &&
       (context.lanes == nullptr || values == nullptr)) ||
      plan_id == semantic::kInvalidIndex) {
    return false;
  }
  CompiledMathIndexSpan span;
  if (!compiled_math_bound_vector_plan_span(program, plan_id, &span)) {
    return false;
  }
  if (span.size == 1 || span.size == 2) {
    const auto &base_op =
        program.bound_vector_ops[static_cast<std::size_t>(span.offset)];
    const CompiledMathBoundVectorOp *cap_op = nullptr;
    if (span.size == 2) {
      cap_op =
          &program.bound_vector_ops[
              static_cast<std::size_t>(span.offset + 1)];
    }
    const bool base_ok =
        base_op.kind == CompiledMathBoundVectorOpKind::DirectBindTime ||
        base_op.kind == CompiledMathBoundVectorOpKind::DirectPreparedTime;
    const bool cap_ok =
        cap_op == nullptr ||
        cap_op->kind == CompiledMathBoundVectorOpKind::MinDirectBindTime ||
        cap_op->kind == CompiledMathBoundVectorOpKind::MinDirectPreparedCap;
    if (base_ok && cap_ok) {
      if (context.frame == nullptr || context.rule == nullptr ||
          context.lanes == nullptr || context.child_width == 0U) {
        return false;
      }
      const bool base_is_bind =
          base_op.kind == CompiledMathBoundVectorOpKind::DirectBindTime;
      const bool cap_is_bind =
          cap_op != nullptr &&
          cap_op->kind == CompiledMathBoundVectorOpKind::MinDirectBindTime;
      const bool needs_bind = base_is_bind || cap_is_bind;
      auto cached_parent = std::numeric_limits<std::size_t>::max();
      double row_shift = 0.0;
      double row_scale = 0.0;
      double fixed_time = 0.0;
      double fixed_cap = 0.0;
      for (std::size_t i = 0; i < context.lane_count; ++i) {
        const auto child_pos = static_cast<std::size_t>(context.lanes[i]);
        const auto parent_pos = child_pos / context.child_width;
        if (parent_pos >= context.eligible_count) {
          return false;
        }
        if (parent_pos != cached_parent) {
          cached_parent = parent_pos;
          const auto op_parent_pos = context.op_parent_offset + parent_pos;
          if (needs_bind) {
            row_shift = context.frame->direct_parent_shift[parent_pos];
            row_scale = context.frame->direct_parent_scale[parent_pos];
          }
          if (!base_is_bind) {
            fixed_time =
                context.frame->direct_source_time_by_op_parent[op_parent_pos];
          }
          if (cap_op != nullptr && !cap_is_bind) {
            fixed_cap =
                context.frame->direct_source_cap_by_op_parent[op_parent_pos];
          }
        }
        const auto q_slot = child_pos - parent_pos * context.child_width;
        const double bind_time =
            needs_bind
                ? row_shift +
                      row_scale *
                          context.rule->nodes[context.q_begin + q_slot]
                : 0.0;
        double value = base_is_bind ? bind_time : fixed_time;
        if (cap_op != nullptr) {
          value = std::min(value, cap_is_bind ? bind_time : fixed_cap);
        }
        values[i] = value;
      }
      return true;
    }
  }
  bool wrote = false;
  for (semantic::Index op_idx = 0; op_idx < span.size; ++op_idx) {
    const auto &op =
        program.bound_vector_ops[
            static_cast<std::size_t>(span.offset + op_idx)];
    switch (op.kind) {
    case CompiledMathBoundVectorOpKind::ConstantZero:
      std::fill(
          values,
          values + static_cast<std::ptrdiff_t>(context.lane_count),
          0.0);
      wrote = true;
      break;
    case CompiledMathBoundVectorOpKind::StateTimeSlot:
      if (context.state == nullptr ||
          op.time_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(op.time_id) >=
              context.state->time_slot_count) {
        return false;
      }
      for (std::size_t i = 0; i < context.lane_count; ++i) {
        values[i] = context.state->time(op.time_id, context.lanes[i]);
      }
      wrote = true;
      break;
    case CompiledMathBoundVectorOpKind::DirectBindTime:
      for (std::size_t i = 0; i < context.lane_count; ++i) {
        std::size_t child_pos = 0U;
        std::size_t parent_pos = 0U;
        if (!compiled_math_bound_vector_direct_child_parent(
                context, i, &child_pos, &parent_pos)) {
          return false;
        }
        const auto q_slot = child_pos - parent_pos * context.child_width;
        values[i] =
            compiled_math_direct_tile_bind_time(
                *context.frame,
                parent_pos,
                context.q_begin + q_slot,
                *context.rule);
      }
      wrote = true;
      break;
    case CompiledMathBoundVectorOpKind::DirectPreparedTime:
    case CompiledMathBoundVectorOpKind::DirectPreparedCap:
      if (context.frame == nullptr || context.lanes == nullptr ||
          context.child_width == 0U) {
        return false;
      }
      for (std::size_t i = 0; i < context.lane_count; ++i) {
        std::size_t child_pos = 0U;
        std::size_t parent_pos = 0U;
        if (!compiled_math_bound_vector_direct_child_parent(
                context, i, &child_pos, &parent_pos)) {
          return false;
        }
        (void)child_pos;
        const auto op_parent_pos = context.op_parent_offset + parent_pos;
        values[i] =
            op.kind == CompiledMathBoundVectorOpKind::DirectPreparedTime
                ? context.frame->direct_source_time_by_op_parent[op_parent_pos]
                : context.frame->direct_source_cap_by_op_parent[op_parent_pos];
      }
      wrote = true;
      break;
    case CompiledMathBoundVectorOpKind::MinDirectBindTime:
      if (!wrote) {
        return false;
      }
      for (std::size_t i = 0; i < context.lane_count; ++i) {
        std::size_t child_pos = 0U;
        std::size_t parent_pos = 0U;
        if (!compiled_math_bound_vector_direct_child_parent(
                context, i, &child_pos, &parent_pos)) {
          return false;
        }
        const auto q_slot = child_pos - parent_pos * context.child_width;
        values[i] =
            std::min(
                values[i],
                compiled_math_direct_tile_bind_time(
                    *context.frame,
                    parent_pos,
                    context.q_begin + q_slot,
                    *context.rule));
      }
      break;
    case CompiledMathBoundVectorOpKind::MinDirectPreparedCap:
      if (!wrote) {
        return false;
      }
      for (std::size_t i = 0; i < context.lane_count; ++i) {
        std::size_t child_pos = 0U;
        std::size_t parent_pos = 0U;
        if (!compiled_math_bound_vector_direct_child_parent(
                context, i, &child_pos, &parent_pos)) {
          return false;
        }
        (void)child_pos;
        const auto op_parent_pos = context.op_parent_offset + parent_pos;
        values[i] =
            std::min(
                values[i],
                context.frame
                    ->direct_source_cap_by_op_parent[op_parent_pos]);
      }
      break;
    case CompiledMathBoundVectorOpKind::SourceSequenceBounds:
    case CompiledMathBoundVectorOpKind::ConditionLower:
    case CompiledMathBoundVectorOpKind::ConditionUpper:
    case CompiledMathBoundVectorOpKind::ConditionExact:
      return false;
    }
  }
  return wrote || context.lane_count == 0U;
}

inline bool compiled_math_bound_vector_parent_condition_time(
    const CompiledMathBoundVectorOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *parent_lanes,
    const std::size_t parent_pos,
    CompiledMathBatchIntegralFrame *frame,
    double *time_out) {
  if (parent_lanes == nullptr || frame == nullptr || time_out == nullptr) {
    return false;
  }
  auto *source_channels = frame->direct_source_channels_by_parent[parent_pos];
  if (source_channels == nullptr) {
    return false;
  }
  const auto parent_lane = parent_lanes[parent_pos];
  if (op.time_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(op.time_id) < state.time_slot_count &&
      state.has_time(op.time_id, parent_lane)) {
    *time_out = state.time(op.time_id, parent_lane);
    return true;
  }
  *time_out = source_channels->compiled_condition_fallback_time();
  return true;
}

inline bool compiled_math_eval_condition_bound_vector_plan_compact(
    const CompiledMathProgram &program,
    const semantic::Index plan_id,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *parent_lanes,
    const std::size_t parent_count,
    CompiledMathBatchIntegralFrame *frame,
    double *lower_out,
    double *upper_out,
    double *exact_out,
    std::uint8_t *has_exact_out) {
  if ((parent_count != 0U && parent_lanes == nullptr) || frame == nullptr ||
      lower_out == nullptr || upper_out == nullptr || exact_out == nullptr ||
      has_exact_out == nullptr) {
    return false;
  }
  CompiledMathIndexSpan span;
  if (!compiled_math_bound_vector_plan_span(program, plan_id, &span)) {
    return false;
  }
  for (semantic::Index op_idx = 0; op_idx < span.size; ++op_idx) {
    const auto &op =
        program.bound_vector_ops[
            static_cast<std::size_t>(span.offset + op_idx)];
    switch (op.kind) {
    case CompiledMathBoundVectorOpKind::SourceSequenceBounds:
      for (std::size_t parent_pos = 0; parent_pos < parent_count;
           ++parent_pos) {
        auto *source_channels =
            frame->direct_source_channels_by_parent[parent_pos];
        if (source_channels == nullptr) {
          return false;
        }
        lower_out[parent_pos] = 0.0;
        upper_out[parent_pos] = std::numeric_limits<double>::infinity();
        exact_out[parent_pos] = 0.0;
        has_exact_out[parent_pos] = 0U;
        source_channels->source_product_sequence_bound_values(
            op.source_id,
            op.bounds,
            &lower_out[parent_pos],
            &upper_out[parent_pos],
            &exact_out[parent_pos],
            &has_exact_out[parent_pos]);
      }
      break;
    case CompiledMathBoundVectorOpKind::ConditionLower:
      for (std::size_t parent_pos = 0; parent_pos < parent_count;
           ++parent_pos) {
        double bound_time = 0.0;
        if (!compiled_math_bound_vector_parent_condition_time(
                op, state, parent_lanes, parent_pos, frame, &bound_time)) {
          return false;
        }
        lower_out[parent_pos] = std::max(lower_out[parent_pos], bound_time);
      }
      break;
    case CompiledMathBoundVectorOpKind::ConditionUpper:
      for (std::size_t parent_pos = 0; parent_pos < parent_count;
           ++parent_pos) {
        double bound_time = 0.0;
        if (!compiled_math_bound_vector_parent_condition_time(
                op, state, parent_lanes, parent_pos, frame, &bound_time)) {
          return false;
        }
        upper_out[parent_pos] = std::min(upper_out[parent_pos], bound_time);
      }
      break;
    case CompiledMathBoundVectorOpKind::ConditionExact:
      for (std::size_t parent_pos = 0; parent_pos < parent_count;
           ++parent_pos) {
        double bound_time = 0.0;
        if (!compiled_math_bound_vector_parent_condition_time(
                op, state, parent_lanes, parent_pos, frame, &bound_time)) {
          return false;
        }
        exact_out[parent_pos] = bound_time;
        has_exact_out[parent_pos] = 1U;
      }
      break;
    case CompiledMathBoundVectorOpKind::ConstantZero:
    case CompiledMathBoundVectorOpKind::StateTimeSlot:
    case CompiledMathBoundVectorOpKind::DirectBindTime:
    case CompiledMathBoundVectorOpKind::DirectPreparedTime:
    case CompiledMathBoundVectorOpKind::DirectPreparedCap:
    case CompiledMathBoundVectorOpKind::MinDirectBindTime:
    case CompiledMathBoundVectorOpKind::MinDirectPreparedCap:
      return false;
    }
  }
  return true;
}

inline void compiled_math_direct_tile_weighted_accumulate_parent_rows(
    const semantic::Index *parent_lanes,
    const std::size_t parent_count,
    const std::size_t child_width,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const double *parent_scale,
    const double *values_by_child_lane,
    double *out_by_parent_lane) {
  for (std::size_t parent_pos = 0; parent_pos < parent_count; ++parent_pos) {
    const auto parent_lane = static_cast<std::size_t>(parent_lanes[parent_pos]);
    const auto row_offset = parent_pos * child_width;
    const double scale = parent_scale[parent_pos];
    double sum = 0.0;
    for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
      const auto pos = row_offset + q_slot;
      const double value = values_by_child_lane[pos];
      if (compiled_math_batch_product_is_live(value)) {
        sum += scale * rule.weights[q_begin + q_slot] * value;
      }
    }
    out_by_parent_lane[parent_lane] += sum;
  }
}

inline bool compiled_math_direct_tile_lognormal_values_dense(
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchSourceProductScratch *scratch,
    std::size_t *live_count_out,
    double *factor_values_out) {
  const auto channel_mask = prepared.op->value_channel_mask;
  const double onset = prepared.channel->leaf_onset_abs_value;
  std::size_t live_count = 0U;
  std::size_t valid_count = 0U;

  for (std::size_t parent_pos = 0; parent_pos < eligible_count; ++parent_pos) {
    const auto op_parent_pos = op_parent_offset + parent_pos;
    const auto *loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
    if (loaded == nullptr) {
      return false;
    }
    const double m = loaded->params[0];
    const double s = loaded->params[1];
    const double q = loaded->q;
    const double t0 = loaded->t0;
    const bool params_valid = std::isfinite(m) && std::isfinite(s) && s > 0.0;
    const double row_shift = frame->direct_parent_shift[parent_pos];
    const double row_scale = frame->direct_parent_scale[parent_pos];
    const double fixed_time =
        prepared.time_is_bind
            ? 0.0
            : frame->direct_source_time_by_op_parent[op_parent_pos];
    const double fixed_cap =
        prepared.has_time_cap && !prepared.cap_is_bind
            ? frame->direct_source_cap_by_op_parent[op_parent_pos]
            : 0.0;
    const auto row_offset = parent_pos * child_width;

    for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
      const auto child_pos = row_offset + q_slot;
      const double bind_time =
          row_shift + row_scale * rule.nodes[q_begin + q_slot];
      double current_time = prepared.time_is_bind ? bind_time : fixed_time;
      if (prepared.has_time_cap) {
        current_time =
            std::min(current_time, prepared.cap_is_bind ? bind_time : fixed_cap);
      }
      const double x = current_time - onset - t0;
      if (!(x > 0.0) || !params_valid) {
        compiled_math_direct_tile_multiply_live_value_at(
            child_pos,
            !(x > 0.0)
                ? batch_source_product_before_onset_value(channel_mask)
                : batch_source_product_finish_base_value(
                      0.0,
                      0.0,
                      q,
                      channel_mask),
            frame,
            &live_count,
            factor_values_out);
        continue;
      }
      scratch->lanes_a[valid_count] = static_cast<semantic::Index>(child_pos);
      scratch->lower[valid_count] = x;
      scratch->pdf_a[valid_count] = m;
      scratch->cdf_a[valid_count] = 1.0 / s;
      scratch->survival_a[valid_count] = q;
      ++valid_count;
    }
  }

  if (valid_count != 0U) {
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const double z =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
        scratch->upper[i] = -0.5 * z * z;
      }
      vvexp(scratch->upper.data(), scratch->upper.data(), &n);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto child_pos = static_cast<std::size_t>(scratch->lanes_a[i]);
        const double base_pdf =
            kInvSqrtTwoPi * scratch->upper[i] * scratch->cdf_a[i] /
            scratch->lower[i];
        compiled_math_direct_tile_multiply_live_pdf_at(
            child_pos,
            base_pdf,
            scratch->survival_a[i],
            frame,
            &live_count,
            factor_values_out);
      }
    } else {
      for (std::size_t i = 0; i < valid_count; ++i) {
        scratch->upper[i] =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
      }
      compiled_math_batch_normal_cdf_from_finite_z(
          scratch->upper.data(),
          scratch->lower.data(),
          valid_count);
      if (channel_mask == kLeafChannelCdf) {
        for (std::size_t i = 0; i < valid_count; ++i) {
          compiled_math_direct_tile_multiply_live_cdf_at(
              static_cast<std::size_t>(scratch->lanes_a[i]),
              scratch->upper[i],
              scratch->survival_a[i],
              frame,
              &live_count,
              factor_values_out);
        }
      } else if (channel_mask == kLeafChannelSurvival) {
        for (std::size_t i = 0; i < valid_count; ++i) {
          compiled_math_direct_tile_multiply_live_survival_at(
              static_cast<std::size_t>(scratch->lanes_a[i]),
              scratch->upper[i],
              scratch->survival_a[i],
              frame,
              &live_count,
              factor_values_out);
        }
      } else {
        for (std::size_t i = 0; i < valid_count; ++i) {
          compiled_math_direct_tile_multiply_live_finished_value_at(
              static_cast<std::size_t>(scratch->lanes_a[i]),
              0.0,
              scratch->upper[i],
              scratch->survival_a[i],
              channel_mask,
              frame,
              &live_count,
              factor_values_out);
        }
      }
    }
  }

  *live_count_out = live_count;
  return true;
}

inline bool compiled_math_direct_tile_lognormal_values_active(
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchSourceProductScratch *scratch,
    std::size_t *live_count_out,
    double *factor_values_out) {
  const auto channel_mask = prepared.op->value_channel_mask;
  const double onset = prepared.channel->leaf_onset_abs_value;
  std::size_t live_count = 0U;
  std::size_t valid_count = 0U;
  auto cached_parent = std::numeric_limits<std::size_t>::max();
  const ExactLoadedLeafInput *loaded = nullptr;
  double m = 0.0;
  double s = 0.0;
  double q = 0.0;
  double t0 = 0.0;
  double row_shift = 0.0;
  double row_scale = 0.0;
  double fixed_time = 0.0;
  double fixed_cap = 0.0;
  bool params_valid = false;

  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
    const auto child_pos = static_cast<std::size_t>(eval_lanes[eval_pos]);
    const auto parent_pos = child_pos / child_width;
    if (parent_pos >= eligible_count) {
      return false;
    }
    if (parent_pos != cached_parent) {
      cached_parent = parent_pos;
      const auto op_parent_pos = op_parent_offset + parent_pos;
      loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
      if (loaded == nullptr) {
        return false;
      }
      m = loaded->params[0];
      s = loaded->params[1];
      q = loaded->q;
      t0 = loaded->t0;
      params_valid = std::isfinite(m) && std::isfinite(s) && s > 0.0;
      row_shift = frame->direct_parent_shift[parent_pos];
      row_scale = frame->direct_parent_scale[parent_pos];
      fixed_time =
          prepared.time_is_bind
              ? 0.0
              : frame->direct_source_time_by_op_parent[op_parent_pos];
      fixed_cap =
          prepared.has_time_cap && !prepared.cap_is_bind
              ? frame->direct_source_cap_by_op_parent[op_parent_pos]
              : 0.0;
    }
    const auto q_slot = child_pos - parent_pos * child_width;
    const double bind_time =
        row_shift + row_scale * rule.nodes[q_begin + q_slot];
    double current_time = prepared.time_is_bind ? bind_time : fixed_time;
    if (prepared.has_time_cap) {
      current_time =
          std::min(current_time, prepared.cap_is_bind ? bind_time : fixed_cap);
    }
    const double x = current_time - onset - t0;
    if (!(x > 0.0) || !params_valid) {
      compiled_math_direct_tile_multiply_live_value_at(
          child_pos,
          !(x > 0.0)
              ? batch_source_product_before_onset_value(channel_mask)
              : batch_source_product_finish_base_value(
                    0.0,
                    0.0,
                    q,
                    channel_mask),
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    scratch->lanes_a[valid_count] = static_cast<semantic::Index>(child_pos);
    scratch->lower[valid_count] = x;
    scratch->pdf_a[valid_count] = m;
    scratch->cdf_a[valid_count] = 1.0 / s;
    scratch->survival_a[valid_count] = q;
    ++valid_count;
  }

  if (valid_count != 0U) {
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const double z =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
        scratch->upper[i] = -0.5 * z * z;
      }
      vvexp(scratch->upper.data(), scratch->upper.data(), &n);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto child_pos = static_cast<std::size_t>(scratch->lanes_a[i]);
        const double base_pdf =
            kInvSqrtTwoPi * scratch->upper[i] * scratch->cdf_a[i] /
            scratch->lower[i];
        compiled_math_direct_tile_multiply_live_pdf_at(
            child_pos,
            base_pdf,
            scratch->survival_a[i],
            frame,
            &live_count,
            factor_values_out);
      }
    } else {
      for (std::size_t i = 0; i < valid_count; ++i) {
        scratch->upper[i] =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
      }
      compiled_math_batch_normal_cdf_from_finite_z(
          scratch->upper.data(),
          scratch->lower.data(),
          valid_count);
      if (channel_mask == kLeafChannelCdf) {
        for (std::size_t i = 0; i < valid_count; ++i) {
          compiled_math_direct_tile_multiply_live_cdf_at(
              static_cast<std::size_t>(scratch->lanes_a[i]),
              scratch->upper[i],
              scratch->survival_a[i],
              frame,
              &live_count,
              factor_values_out);
        }
      } else if (channel_mask == kLeafChannelSurvival) {
        for (std::size_t i = 0; i < valid_count; ++i) {
          compiled_math_direct_tile_multiply_live_survival_at(
              static_cast<std::size_t>(scratch->lanes_a[i]),
              scratch->upper[i],
              scratch->survival_a[i],
              frame,
              &live_count,
              factor_values_out);
        }
      } else {
        for (std::size_t i = 0; i < valid_count; ++i) {
          compiled_math_direct_tile_multiply_live_finished_value_at(
              static_cast<std::size_t>(scratch->lanes_a[i]),
              0.0,
              scratch->upper[i],
              scratch->survival_a[i],
              channel_mask,
              frame,
              &live_count,
              factor_values_out);
        }
      }
    }
  }

  *live_count_out = live_count;
  return true;
}

inline bool compiled_math_direct_tile_lognormal_values(
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
	    CompiledMathBatchIntegralFrame *frame,
	    CompiledMathBatchSourceProductScratch *scratch,
	    std::size_t *live_count_out,
	    double *factor_values_out) {
	  if (prepared.op == nullptr || prepared.channel == nullptr ||
	      frame == nullptr || scratch == nullptr ||
	      live_count_out == nullptr || factor_values_out == nullptr ||
	      prepared.channel->leaf_index == semantic::kInvalidIndex ||
      static_cast<leaf::DistKind>(prepared.channel->leaf_dist_kind) !=
          leaf::DistKind::Lognormal) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  if (eval_lanes == nullptr && eval_count == tile_capacity) {
    return compiled_math_direct_tile_lognormal_values_dense(
        prepared,
        op_parent_offset,
        eligible_count,
        q_begin,
        rule,
        child_width,
        frame,
        scratch,
        live_count_out,
        factor_values_out);
  }
  if (eval_lanes != nullptr) {
    return compiled_math_direct_tile_lognormal_values_active(
        prepared,
        op_parent_offset,
        eligible_count,
        q_begin,
        rule,
        child_width,
        eval_lanes,
        eval_count,
        frame,
        scratch,
        live_count_out,
        factor_values_out);
  }
  std::size_t live_count = 0U;
  std::size_t valid_count = 0U;
  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
    const auto child_pos =
        compiled_math_direct_tile_child_pos(eval_lanes, eval_pos);
    const auto parent_pos = child_pos / child_width;
    if (parent_pos >= eligible_count) {
      return false;
    }
    const auto child_lane = static_cast<semantic::Index>(child_pos);
    const auto q_slot = child_pos - parent_pos * child_width;
    const auto op_parent_pos = op_parent_offset + parent_pos;
    const auto *loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
    if (loaded == nullptr) {
      return false;
    }
    const double bind_time =
        compiled_math_direct_tile_bind_time(
            *frame,
            parent_pos,
            q_begin + q_slot,
            rule);
    double current_time =
        prepared.time_is_bind
            ? bind_time
            : frame->direct_source_time_by_op_parent[op_parent_pos];
    if (prepared.has_time_cap) {
      const double cap =
          prepared.cap_is_bind
              ? bind_time
              : frame->direct_source_cap_by_op_parent[op_parent_pos];
      current_time = std::min(current_time, cap);
    }
    const double x =
        current_time - prepared.channel->leaf_onset_abs_value - loaded->t0;
    const double m = loaded->params[0];
    const double s = loaded->params[1];
    if (!(x > 0.0) || !std::isfinite(m) || !std::isfinite(s) ||
        s <= 0.0) {
      compiled_math_direct_tile_multiply_value_at(
          child_pos,
          !(x > 0.0)
              ? batch_source_product_before_onset_value(
                    prepared.op->value_channel_mask)
              : batch_source_product_finish_base_value(
                    0.0,
                    0.0,
                    loaded->q,
                    prepared.op->value_channel_mask),
	          frame,
	          &live_count,
	          factor_values_out);
      continue;
    }
    scratch->lanes_a[valid_count] = child_lane;
    scratch->lower[valid_count] = x;
    scratch->pdf_a[valid_count] = m;
    scratch->cdf_a[valid_count] = 1.0 / s;
    scratch->survival_a[valid_count] = loaded->q;
    ++valid_count;
  }

  if (valid_count != 0U) {
    const auto n = static_cast<int>(valid_count);
    vvlog(scratch->upper.data(), scratch->lower.data(), &n);
    if (prepared.op->value_channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < valid_count; ++i) {
        const double z =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
        scratch->upper[i] = -0.5 * z * z;
      }
      vvexp(scratch->upper.data(), scratch->upper.data(), &n);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto child_pos = static_cast<std::size_t>(scratch->lanes_a[i]);
        const double base_pdf =
            kInvSqrtTwoPi * scratch->upper[i] * scratch->cdf_a[i] /
            scratch->lower[i];
        compiled_math_direct_tile_multiply_finished_value_at(
            child_pos,
            base_pdf,
            0.0,
            scratch->survival_a[i],
            prepared.op->value_channel_mask,
	            frame,
	            &live_count,
	            factor_values_out);
      }
    } else {
      for (std::size_t i = 0; i < valid_count; ++i) {
        scratch->upper[i] =
            (scratch->upper[i] - scratch->pdf_a[i]) * scratch->cdf_a[i];
      }
      compiled_math_batch_normal_cdf_from_finite_z(
          scratch->upper.data(),
          scratch->lower.data(),
          valid_count);
      for (std::size_t i = 0; i < valid_count; ++i) {
        const auto child_pos = static_cast<std::size_t>(scratch->lanes_a[i]);
        compiled_math_direct_tile_multiply_finished_value_at(
            child_pos,
            0.0,
            scratch->upper[i],
            scratch->survival_a[i],
            prepared.op->value_channel_mask,
	            frame,
	            &live_count,
	            factor_values_out);
      }
    }
  }

  *live_count_out = live_count;
  return true;
#else
	  (void)op_parent_offset;
	  (void)eligible_count;
	  (void)q_begin;
	  (void)rule;
	  (void)child_width;
  (void)tile_capacity;
  (void)eval_count;
  return false;
#endif
}

inline bool compiled_math_direct_tile_gamma_update(
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchSourceProductScratch *scratch,
    std::size_t *live_count_out,
    double *factor_values_out = nullptr) {
  if (prepared.op == nullptr || prepared.channel == nullptr ||
      frame == nullptr || scratch == nullptr ||
      live_count_out == nullptr ||
      prepared.channel->leaf_index == semantic::kInvalidIndex ||
      static_cast<leaf::DistKind>(prepared.channel->leaf_dist_kind) !=
          leaf::DistKind::Gamma) {
    return false;
  }
  const auto channel_mask = prepared.op->value_channel_mask;
	  std::size_t live_count = 0U;
	  std::size_t valid_count = 0U;
	  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
	    const auto child_pos =
	        compiled_math_direct_tile_child_pos(eval_lanes, eval_pos);
	    const auto parent_pos = child_pos / child_width;
	    if (parent_pos >= eligible_count) {
	      return false;
	    }
	    const auto child_lane = static_cast<semantic::Index>(child_pos);
	    const auto q_slot = child_pos - parent_pos * child_width;
	    const auto op_parent_pos = op_parent_offset + parent_pos;
	    const auto *loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
	    if (loaded == nullptr) {
	      return false;
	    }
	    const double bind_time =
	        compiled_math_direct_tile_bind_time(
	            *frame,
	            parent_pos,
	            q_begin + q_slot,
	            rule);
    double current_time =
        prepared.time_is_bind
            ? bind_time
            : frame->direct_source_time_by_op_parent[op_parent_pos];
    if (prepared.has_time_cap) {
      const double cap =
          prepared.cap_is_bind
              ? bind_time
              : frame->direct_source_cap_by_op_parent[op_parent_pos];
      current_time = std::min(current_time, cap);
    }
    const double x =
        current_time - prepared.channel->leaf_onset_abs_value - loaded->t0;
    const double shape = loaded->params[0];
    const double rate = loaded->params[1];
    if (!(x > 0.0) || !std::isfinite(shape) || shape <= 0.0 ||
        !std::isfinite(rate) || rate <= 0.0) {
      compiled_math_direct_tile_multiply_value_at(
          child_pos,
          !(x > 0.0)
              ? batch_source_product_before_onset_value(channel_mask)
              : batch_source_product_finish_base_value(
                    0.0,
                    0.0,
                    loaded->q,
                    channel_mask),
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    scratch->lanes_a[valid_count] = child_lane;
    scratch->lower[valid_count] = x;
    scratch->pdf_a[valid_count] = shape;
    scratch->cdf_a[valid_count] = rate;
    scratch->survival_a[valid_count] = loaded->q;
    ++valid_count;
  }
  if (valid_count == 0U) {
    *live_count_out = live_count;
    return true;
  }
  if (channel_mask != kLeafChannelPdf) {
    for (std::size_t i = 0; i < valid_count; ++i) {
      scratch->upper[i] = scratch->cdf_a[i] * scratch->lower[i];
    }
    compiled_math_batch_gamma_regularized_lower_values(
        scratch->pdf_a.data(),
        scratch->upper.data(),
        valid_count,
        scratch->exact.data());
    for (std::size_t i = 0; i < valid_count; ++i) {
      compiled_math_direct_tile_multiply_finished_value_at(
          static_cast<std::size_t>(scratch->lanes_a[i]),
          0.0,
          scratch->exact[i],
          scratch->survival_a[i],
          channel_mask,
          frame,
          &live_count,
          factor_values_out);
    }
    *live_count_out = live_count;
    return true;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  const auto n = static_cast<int>(valid_count);
  vvlog(scratch->upper.data(), scratch->lower.data(), &n);
  for (std::size_t i = 0; i < valid_count; ++i) {
    const double shape = scratch->pdf_a[i];
    const double rate = scratch->cdf_a[i];
    scratch->upper[i] =
        shape * std::log(rate) + (shape - 1.0) * scratch->upper[i] -
        rate * scratch->lower[i] - std::lgamma(shape);
  }
  vvexp(scratch->upper.data(), scratch->upper.data(), &n);
#else
  for (std::size_t i = 0; i < valid_count; ++i) {
    scratch->upper[i] =
        std::exp(
            scratch->pdf_a[i] * std::log(scratch->cdf_a[i]) +
            (scratch->pdf_a[i] - 1.0) * std::log(scratch->lower[i]) -
            scratch->cdf_a[i] * scratch->lower[i] -
            std::lgamma(scratch->pdf_a[i]));
  }
#endif
  for (std::size_t i = 0; i < valid_count; ++i) {
    compiled_math_direct_tile_multiply_finished_value_at(
        static_cast<std::size_t>(scratch->lanes_a[i]),
        scratch->upper[i],
        0.0,
        scratch->survival_a[i],
        channel_mask,
        frame,
        &live_count,
        factor_values_out);
  }
  *live_count_out = live_count;
  return true;
}

inline bool compiled_math_direct_tile_exgauss_update(
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchSourceProductScratch *scratch,
    std::size_t *live_count_out,
    double *factor_values_out = nullptr) {
  if (prepared.op == nullptr || prepared.channel == nullptr ||
      frame == nullptr || scratch == nullptr ||
      live_count_out == nullptr ||
      prepared.channel->leaf_index == semantic::kInvalidIndex ||
      static_cast<leaf::DistKind>(prepared.channel->leaf_dist_kind) !=
          leaf::DistKind::Exgauss) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  const auto channel_mask = prepared.op->value_channel_mask;
	  std::size_t live_count = 0U;
	  std::size_t valid_count = 0U;
	  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
	    const auto child_pos =
	        compiled_math_direct_tile_child_pos(eval_lanes, eval_pos);
	    const auto parent_pos = child_pos / child_width;
	    if (parent_pos >= eligible_count) {
	      return false;
	    }
	    const auto child_lane = static_cast<semantic::Index>(child_pos);
	    const auto q_slot = child_pos - parent_pos * child_width;
	    const auto op_parent_pos = op_parent_offset + parent_pos;
	    const auto *loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
	    if (loaded == nullptr) {
	      return false;
	    }
	    const double bind_time =
	        compiled_math_direct_tile_bind_time(
	            *frame,
	            parent_pos,
	            q_begin + q_slot,
	            rule);
    double current_time =
        prepared.time_is_bind
            ? bind_time
            : frame->direct_source_time_by_op_parent[op_parent_pos];
    if (prepared.has_time_cap) {
      const double cap =
          prepared.cap_is_bind
              ? bind_time
              : frame->direct_source_cap_by_op_parent[op_parent_pos];
      current_time = std::min(current_time, cap);
    }
    const double x =
        current_time - prepared.channel->leaf_onset_abs_value - loaded->t0;
    const double mu = loaded->params[0];
    const double sigma = loaded->params[1];
    const double tau = loaded->params[2];
    if (!(x > 0.0) || !std::isfinite(mu) || !std::isfinite(sigma) ||
        sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
      compiled_math_direct_tile_multiply_value_at(
          child_pos,
          !(x > 0.0)
              ? batch_source_product_before_onset_value(channel_mask)
              : batch_source_product_finish_base_value(
                    0.0,
                    0.0,
                    loaded->q,
                    channel_mask),
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    const double inv_tau = 1.0 / tau;
    const double sigma_over_tau = sigma * inv_tau;
    const double z0 = -mu / sigma;
    scratch->times[child_pos] = x;
    scratch->pdf_b[child_pos] = mu;
    scratch->cdf_b[child_pos] = sigma;
    scratch->survival_b[child_pos] = tau;
    scratch->survival_a[child_pos] = loaded->q;
    scratch->lanes_a[valid_count] = child_lane;
    scratch->upper[valid_count] = z0;
    scratch->exact[valid_count] = z0 - sigma_over_tau;
    scratch->lower[valid_count] =
        sigma * sigma * inv_tau * inv_tau * 0.5 + mu * inv_tau;
    ++valid_count;
  }
  if (valid_count == 0U) {
    *live_count_out = live_count;
    return true;
  }
  const auto n = static_cast<int>(valid_count);
  vvexp(scratch->lower.data(), scratch->lower.data(), &n);
  compiled_math_batch_normal_cdf_from_z(
      scratch->upper.data(),
      scratch->cdf_a.data(),
      valid_count);
  compiled_math_batch_normal_cdf_from_z(
      scratch->exact.data(),
      scratch->cdf_a.data(),
      valid_count);
  std::size_t conditioned_count = 0U;
  for (std::size_t i = 0; i < valid_count; ++i) {
    const auto child_lane = scratch->lanes_a[i];
    const auto child_pos = static_cast<std::size_t>(child_lane);
    const double lower_cdf =
        clamp_probability(scratch->upper[i] - scratch->lower[i] * scratch->exact[i]);
    const double lower_survival = 1.0 - lower_cdf;
    if (!(lower_survival > 0.0)) {
      compiled_math_direct_tile_multiply_value_at(
          child_pos,
          batch_source_product_before_onset_value(channel_mask),
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    scratch->cdf_c[child_pos] = lower_cdf;
    scratch->survival_c[child_pos] = lower_survival;
    scratch->lanes_b[conditioned_count++] = child_lane;
  }
  if (conditioned_count == 0U) {
    *live_count_out = live_count;
    return true;
  }
  for (std::size_t i = 0; i < conditioned_count; ++i) {
    const auto child_lane = scratch->lanes_b[i];
    const auto child_pos = static_cast<std::size_t>(child_lane);
    const double x = scratch->times[child_pos];
    const double mu = scratch->pdf_b[child_pos];
    const double sigma = scratch->cdf_b[child_pos];
    const double tau = scratch->survival_b[child_pos];
    const double inv_tau = 1.0 / tau;
    const double sigma_over_tau = sigma * inv_tau;
    const double z = (x - mu) / sigma;
    scratch->pdf_a[i] = tau;
    scratch->pdf_c[i] = scratch->survival_a[child_pos];
    scratch->onset_upper[i] = scratch->cdf_c[child_pos];
    scratch->lower[i] =
        sigma * sigma * inv_tau * inv_tau * 0.5 -
        (x - mu) * inv_tau;
    scratch->exact[i] = z - sigma_over_tau;
    if (channel_mask != kLeafChannelPdf) {
      scratch->upper[i] = z;
    }
  }
  const auto conditioned_n = static_cast<int>(conditioned_count);
  vvexp(scratch->lower.data(), scratch->lower.data(), &conditioned_n);
  compiled_math_batch_normal_cdf_from_z(
      scratch->exact.data(),
      scratch->cdf_a.data(),
      conditioned_count);
  if (channel_mask != kLeafChannelPdf) {
    compiled_math_batch_normal_cdf_from_z(
        scratch->upper.data(),
        scratch->cdf_a.data(),
        conditioned_count);
  }
  for (std::size_t i = 0; i < conditioned_count; ++i) {
    const auto child_pos = static_cast<std::size_t>(scratch->lanes_b[i]);
    const double lower_cdf = scratch->onset_upper[i];
    const double lower_survival = 1.0 - lower_cdf;
    if (channel_mask == kLeafChannelPdf) {
      const double base_pdf =
          (scratch->lower[i] * scratch->exact[i] / scratch->pdf_a[i]) /
          lower_survival;
      compiled_math_direct_tile_multiply_finished_value_at(
          child_pos,
          base_pdf,
          0.0,
          scratch->pdf_c[i],
          channel_mask,
          frame,
          &live_count,
          factor_values_out);
    } else {
      const double raw_cdf =
          clamp_probability(scratch->upper[i] -
                            scratch->lower[i] * scratch->exact[i]);
      const double base_cdf = (raw_cdf - lower_cdf) / lower_survival;
      compiled_math_direct_tile_multiply_finished_value_at(
          child_pos,
          0.0,
          base_cdf,
          scratch->pdf_c[i],
          channel_mask,
          frame,
          &live_count,
          factor_values_out);
    }
  }
  *live_count_out = live_count;
  return true;
#else
	  (void)op_parent_offset;
	  (void)eligible_count;
	  (void)q_begin;
	  (void)rule;
	  (void)child_width;
  (void)tile_capacity;
  (void)eval_count;
  return false;
#endif
}

inline bool compiled_math_direct_tile_lba_update(
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchSourceProductScratch *scratch,
    std::size_t *live_count_out,
    double *factor_values_out = nullptr) {
  if (prepared.op == nullptr || prepared.channel == nullptr ||
      frame == nullptr || scratch == nullptr ||
      live_count_out == nullptr ||
      prepared.channel->leaf_index == semantic::kInvalidIndex ||
      static_cast<leaf::DistKind>(prepared.channel->leaf_dist_kind) !=
          leaf::DistKind::LBA) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  const auto channel_mask = prepared.op->value_channel_mask;
	  std::size_t live_count = 0U;
	  std::size_t valid_count = 0U;
	  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
	    const auto child_pos =
	        compiled_math_direct_tile_child_pos(eval_lanes, eval_pos);
	    const auto parent_pos = child_pos / child_width;
	    if (parent_pos >= eligible_count) {
	      return false;
	    }
	    const auto child_lane = static_cast<semantic::Index>(child_pos);
	    const auto q_slot = child_pos - parent_pos * child_width;
	    const auto op_parent_pos = op_parent_offset + parent_pos;
	    const auto *loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
	    if (loaded == nullptr) {
	      return false;
	    }
	    const double bind_time =
	        compiled_math_direct_tile_bind_time(
	            *frame,
	            parent_pos,
	            q_begin + q_slot,
	            rule);
    double current_time =
        prepared.time_is_bind
            ? bind_time
            : frame->direct_source_time_by_op_parent[op_parent_pos];
    if (prepared.has_time_cap) {
      const double cap =
          prepared.cap_is_bind
              ? bind_time
              : frame->direct_source_cap_by_op_parent[op_parent_pos];
      current_time = std::min(current_time, cap);
    }
    const double x =
        current_time - prepared.channel->leaf_onset_abs_value - loaded->t0;
    const double v = loaded->params[0];
    const double B = loaded->params[1];
    const double A = loaded->params[2];
    const double sv = loaded->params[3];
    if (!(x > 0.0)) {
      compiled_math_direct_tile_multiply_value_at(
          child_pos,
          batch_source_product_before_onset_value(channel_mask),
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    if (!std::isfinite(v) || !std::isfinite(B) ||
        !std::isfinite(A) || !std::isfinite(sv) || sv <= 0.0) {
      compiled_math_direct_tile_multiply_finished_value_at(
          child_pos,
          0.0,
          0.0,
	          loaded->q,
		          channel_mask,
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    scratch->lanes_a[valid_count] = child_lane;
    scratch->lower[valid_count] = x;
    scratch->upper[valid_count] = v;
    scratch->exact[valid_count] = B;
    scratch->pdf_a[valid_count] = A;
    scratch->cdf_a[valid_count] = sv;
    scratch->survival_a[valid_count] = loaded->q;
    scratch->pdf_b[valid_count] = v / sv;
    ++valid_count;
  }
  if (valid_count == 0U) {
    *live_count_out = live_count;
    return true;
  }

  compiled_math_batch_normal_cdf_from_z(
      scratch->pdf_b.data(),
      scratch->survival_b.data(),
      valid_count);
  for (std::size_t i = 0; i < valid_count; ++i) {
    if (!std::isfinite(scratch->pdf_b[i]) || scratch->pdf_b[i] < 1e-10) {
      scratch->pdf_b[i] = 1e-10;
    }
  }

  std::size_t wide_count = 0U;
  std::size_t point_count = 0U;
  for (std::size_t i = 0; i < valid_count; ++i) {
    const double x = scratch->lower[i];
    const double v = scratch->upper[i];
    const double B = scratch->exact[i];
    const double A = scratch->pdf_a[i];
    const double sv = scratch->cdf_a[i];
    if (A > 1e-10) {
      const double zs = x * sv;
      if (!(zs > 0.0) || !std::isfinite(zs)) {
        compiled_math_direct_tile_multiply_finished_value_at(
            static_cast<std::size_t>(scratch->lanes_a[i]),
            0.0,
            0.0,
            scratch->survival_a[i],
            channel_mask,
            frame,
            &live_count,
            factor_values_out);
        continue;
      }
      const double cmz = B - x * v;
      const double xx = cmz - A;
      const double cz = cmz / zs;
      const double cz_max = xx / zs;
      scratch->lanes_b[wide_count] = static_cast<semantic::Index>(i);
      scratch->cdf_b[wide_count] = cz;
      scratch->survival_b[wide_count] = cz_max;
      scratch->times[wide_count] = -0.5 * cz * cz;
      scratch->onset_upper[wide_count] = -0.5 * cz_max * cz_max;
      ++wide_count;
    } else {
      const double z = (B / x - v) / sv;
      scratch->lanes_c[point_count] = static_cast<semantic::Index>(i);
      scratch->cdf_c[point_count] = z;
      scratch->survival_c[point_count] = -0.5 * z * z;
      ++point_count;
    }
  }

  if (wide_count != 0U) {
    const auto n = static_cast<int>(wide_count);
    compiled_math_batch_normal_cdf_from_z(
        scratch->cdf_b.data(),
        scratch->mask_a.empty() ? nullptr : scratch->pdf_c.data(),
        wide_count);
    compiled_math_batch_normal_cdf_from_z(
        scratch->survival_b.data(),
        scratch->mask_a.empty() ? nullptr : scratch->pdf_c.data(),
        wide_count);
    vvexp(scratch->times.data(), scratch->times.data(), &n);
    vvexp(scratch->onset_upper.data(), scratch->onset_upper.data(), &n);
    for (std::size_t j = 0; j < wide_count; ++j) {
      const auto src = static_cast<std::size_t>(scratch->lanes_b[j]);
      const auto child_pos = static_cast<std::size_t>(scratch->lanes_a[src]);
      const double x = scratch->lower[src];
      const double v = scratch->upper[src];
      const double B = scratch->exact[src];
      const double A = scratch->pdf_a[src];
      const double sv = scratch->cdf_a[src];
      const double denom = scratch->pdf_b[src];
      const double cmz = B - x * v;
      const double xx = cmz - A;
      const double phi_cz = kInvSqrtTwoPi * scratch->times[j];
      const double phi_cz_max = kInvSqrtTwoPi * scratch->onset_upper[j];
      double base_pdf = 0.0;
      double base_cdf = 0.0;
      if (channel_mask == kLeafChannelPdf) {
        base_pdf =
            (v * (scratch->cdf_b[j] - scratch->survival_b[j]) +
             sv * (phi_cz_max - phi_cz)) /
            (A * denom);
      } else {
        base_cdf =
            (1.0 + (x * sv * (phi_cz_max - phi_cz) +
                    xx * scratch->survival_b[j] -
                    cmz * scratch->cdf_b[j]) /
                       A) /
            denom;
      }
      compiled_math_direct_tile_multiply_finished_value_at(
          child_pos,
          base_pdf,
          base_cdf,
          scratch->survival_a[src],
          channel_mask,
          frame,
          &live_count,
          factor_values_out);
    }
  }

  if (point_count != 0U) {
    compiled_math_batch_normal_cdf_from_z(
        scratch->cdf_c.data(),
        scratch->pdf_c.data(),
        point_count);
    const auto n = static_cast<int>(point_count);
    vvexp(scratch->survival_c.data(), scratch->survival_c.data(), &n);
    for (std::size_t j = 0; j < point_count; ++j) {
      const auto src = static_cast<std::size_t>(scratch->lanes_c[j]);
      const auto child_pos = static_cast<std::size_t>(scratch->lanes_a[src]);
      const double x = scratch->lower[src];
      const double B = scratch->exact[src];
      const double sv = scratch->cdf_a[src];
      const double denom = scratch->pdf_b[src];
      double base_pdf = 0.0;
      double base_cdf = 0.0;
      if (channel_mask == kLeafChannelPdf) {
        base_pdf =
            kInvSqrtTwoPi * scratch->survival_c[j] * B /
            (sv * x * x * denom);
      } else {
        base_cdf = clamp_probability(1.0 - scratch->cdf_c[j]) / denom;
      }
      compiled_math_direct_tile_multiply_finished_value_at(
          child_pos,
          base_pdf,
          base_cdf,
          scratch->survival_a[src],
          channel_mask,
          frame,
          &live_count,
          factor_values_out);
    }
  }
  *live_count_out = live_count;
  return true;
#else
	  (void)op_parent_offset;
	  (void)eligible_count;
	  (void)q_begin;
	  (void)rule;
	  (void)child_width;
  (void)tile_capacity;
  (void)eval_count;
  return false;
#endif
}

inline bool compiled_math_direct_tile_rdm_update(
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchSourceProductScratch *scratch,
    std::size_t *live_count_out,
    double *factor_values_out = nullptr) {
  if (prepared.op == nullptr || prepared.channel == nullptr ||
      frame == nullptr || scratch == nullptr ||
      live_count_out == nullptr ||
      prepared.channel->leaf_index == semantic::kInvalidIndex ||
      static_cast<leaf::DistKind>(prepared.channel->leaf_dist_kind) !=
          leaf::DistKind::RDM) {
    return false;
  }
#if defined(ACCUMULATR_USE_ACCELERATE_VFORCE)
  const auto channel_mask = prepared.op->value_channel_mask;
  std::size_t live_count = 0U;
	  std::size_t common_count = 0U;
	  std::size_t point_count = 0U;
	  std::size_t small_drift_count = 0U;
	  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
	    const auto child_pos =
	        compiled_math_direct_tile_child_pos(eval_lanes, eval_pos);
	    const auto parent_pos = child_pos / child_width;
	    if (parent_pos >= eligible_count) {
	      return false;
	    }
	    const auto child_lane = static_cast<semantic::Index>(child_pos);
	    const auto q_slot = child_pos - parent_pos * child_width;
	    const auto op_parent_pos = op_parent_offset + parent_pos;
	    const auto *loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
	    if (loaded == nullptr) {
	      return false;
	    }
	    const double bind_time =
	        compiled_math_direct_tile_bind_time(
	            *frame,
	            parent_pos,
	            q_begin + q_slot,
	            rule);
    double current_time =
        prepared.time_is_bind
            ? bind_time
            : frame->direct_source_time_by_op_parent[op_parent_pos];
    if (prepared.has_time_cap) {
      const double cap =
          prepared.cap_is_bind
              ? bind_time
              : frame->direct_source_cap_by_op_parent[op_parent_pos];
      current_time = std::min(current_time, cap);
    }
    const double x =
        current_time - prepared.channel->leaf_onset_abs_value - loaded->t0;
    const double v = loaded->params[0];
    const double B = loaded->params[1];
    const double A = loaded->params[2];
    const double s = loaded->params[3];
    if (!(x > 0.0)) {
      compiled_math_direct_tile_multiply_value_at(
          child_pos,
          batch_source_product_before_onset_value(channel_mask),
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    if (!std::isfinite(s) || s <= 0.0 || !std::isfinite(v) ||
        !std::isfinite(B) || !std::isfinite(A)) {
      compiled_math_direct_tile_multiply_finished_value_at(
          child_pos,
          0.0,
          0.0,
	          loaded->q,
		          channel_mask,
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    const double v_sc = v / s;
    const double B_sc = B / s;
    const double A_sc = A / s;
    const double a = 0.5 * A_sc;
    const double k = B_sc + a;
    const double l = v_sc;
    if (!std::isfinite(a) || !std::isfinite(k) || !std::isfinite(l)) {
      compiled_math_direct_tile_multiply_finished_value_at(
	          child_pos,
	          0.0,
	          0.0,
	          loaded->q,
	          channel_mask,
          frame,
          &live_count,
          factor_values_out);
      continue;
    }
    if (a < 1e-10) {
      scratch->lanes_c[point_count] = child_lane;
      scratch->cdf_a[point_count] = x;
      scratch->cdf_b[point_count] = k;
      scratch->cdf_c[point_count] = l;
      scratch->survival_c[point_count] = loaded->q;
      ++point_count;
      continue;
    }
    if (l < 1e-10) {
      scratch->lanes_d[small_drift_count] = child_lane;
      scratch->pdf_b[small_drift_count] = x;
      scratch->survival_b[small_drift_count] = k;
      scratch->pdf_c[small_drift_count] = l;
      scratch->onset_upper[small_drift_count] = a;
      scratch->survival_c[small_drift_count] = loaded->q;
      ++small_drift_count;
      continue;
    }
    scratch->lanes_b[common_count] = child_lane;
    scratch->lower[common_count] = x;
    scratch->upper[common_count] = k;
    scratch->exact[common_count] = l;
    scratch->pdf_a[common_count] = a;
    scratch->survival_a[common_count] = loaded->q;
    ++common_count;
  }

  if (point_count != 0U) {
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < point_count; ++i) {
        const double x = scratch->cdf_a[i];
        const double k = scratch->cdf_b[i];
        const double l = scratch->cdf_c[i];
        const double lambda = k * k;
        double exponent = 0.0;
        if (l == 0.0) {
          exponent = -0.5 * lambda / x;
        } else {
          const double mu = k / l;
          exponent =
              -(lambda / (2.0 * x)) *
              ((x * x) / (mu * mu) - 2.0 * x / mu + 1.0);
        }
        scratch->times[i] =
            exponent + 0.5 * std::log(lambda) -
            0.5 * std::log(2.0 * x * x * x * M_PI);
      }
      const auto n = static_cast<int>(point_count);
      vvexp(scratch->times.data(), scratch->times.data(), &n);
      for (std::size_t i = 0; i < point_count; ++i) {
        compiled_math_direct_tile_multiply_finished_value_at(
            static_cast<std::size_t>(scratch->lanes_c[i]),
            scratch->times[i],
            0.0,
            scratch->survival_c[i],
            channel_mask,
            frame,
            &live_count,
            factor_values_out);
      }
    } else {
      for (std::size_t i = 0; i < point_count; ++i) {
        const double x = scratch->cdf_a[i];
        const double k = scratch->cdf_b[i];
        const double l = scratch->cdf_c[i];
        if (std::fabs(l) < 1e-12) {
          scratch->times[i] = k / std::sqrt(x);
          scratch->onset_upper[i] =
              std::numeric_limits<double>::quiet_NaN();
        } else {
          const double mu = k / l;
          const double z_base = std::sqrt((k * k) / x);
          scratch->times[i] = z_base * (1.0 + x / mu);
          scratch->onset_upper[i] = z_base * (1.0 - x / mu);
        }
      }
      compiled_math_batch_normal_cdf_from_z(
          scratch->times.data(),
          scratch->lower.data(),
          point_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->onset_upper.data(),
          scratch->lower.data(),
          point_count);
      for (std::size_t i = 0; i < point_count; ++i) {
        const double k = scratch->cdf_b[i];
        const double l = scratch->cdf_c[i];
        double base_cdf = 0.0;
        if (std::fabs(l) < 1e-12) {
          base_cdf = 2.0 * (1.0 - scratch->times[i]);
        } else {
          const double p1 = 1.0 - scratch->times[i];
          const double p2 = 1.0 - scratch->onset_upper[i];
          base_cdf =
              std::exp(
                  2.0 * k * l +
                  std::log(std::max(1e-300, p1))) +
              p2;
        }
        compiled_math_direct_tile_multiply_finished_value_at(
            static_cast<std::size_t>(scratch->lanes_c[i]),
            0.0,
            base_cdf,
            scratch->survival_c[i],
            channel_mask,
            frame,
            &live_count,
            factor_values_out);
      }
    }
  }

  if (small_drift_count != 0U) {
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < small_drift_count; ++i) {
        const double x = scratch->pdf_b[i];
        const double k = scratch->survival_b[i];
        const double a = scratch->onset_upper[i];
        scratch->times[i] = -((k - a) * (k - a)) / (2.0 * x);
        scratch->upper[i] = -((k + a) * (k + a)) / (2.0 * x);
      }
      const auto n = static_cast<int>(small_drift_count);
      vvexp(scratch->times.data(), scratch->times.data(), &n);
      vvexp(scratch->upper.data(), scratch->upper.data(), &n);
      for (std::size_t i = 0; i < small_drift_count; ++i) {
        const double x = scratch->pdf_b[i];
        const double a = scratch->onset_upper[i];
        const double term = scratch->times[i] - scratch->upper[i];
        const double base_pdf =
            std::exp(
                -0.5 * (M_LN2 + kLogPi + std::log(x)) +
                std::log(std::max(1e-300, term)) -
                M_LN2 - std::log(a));
        compiled_math_direct_tile_multiply_finished_value_at(
            static_cast<std::size_t>(scratch->lanes_d[i]),
            base_pdf,
            0.0,
            scratch->survival_c[i],
            channel_mask,
            frame,
            &live_count,
            factor_values_out);
      }
    } else {
      for (std::size_t i = 0; i < small_drift_count; ++i) {
        const double x = scratch->pdf_b[i];
        const double k = scratch->survival_b[i];
        const double a = scratch->onset_upper[i];
        const double sqt = std::sqrt(x);
        scratch->times[i] = (k + a) / sqt;
        scratch->upper[i] = (-k - a) / sqt;
        scratch->lower[i] =
            -0.5 * (((k + a) * (k + a)) / x - M_LN2 - kLogPi +
                    std::log(x)) -
            std::log(a);
        scratch->exact[i] =
            -0.5 * (((k - a) * (k - a)) / x - M_LN2 - kLogPi +
                    std::log(x)) -
            std::log(a);
      }
      compiled_math_batch_normal_cdf_from_z(
          scratch->times.data(),
          scratch->cdf_b.data(),
          small_drift_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->upper.data(),
          scratch->cdf_b.data(),
          small_drift_count);
      const auto n = static_cast<int>(small_drift_count);
      vvexp(scratch->lower.data(), scratch->lower.data(), &n);
      vvexp(scratch->exact.data(), scratch->exact.data(), &n);
      for (std::size_t i = 0; i < small_drift_count; ++i) {
        const double k = scratch->survival_b[i];
        const double a = scratch->onset_upper[i];
        const double t5a = 2.0 * scratch->times[i] - 1.0;
        const double t5b = 2.0 * scratch->upper[i] - 1.0;
        const double base_cdf =
            1.0 + scratch->lower[i] - scratch->exact[i] +
            ((-k + a) * t5a - (k - a) * t5b) / (2.0 * a);
        compiled_math_direct_tile_multiply_finished_value_at(
            static_cast<std::size_t>(scratch->lanes_d[i]),
            0.0,
            base_cdf,
            scratch->survival_c[i],
            channel_mask,
            frame,
            &live_count,
            factor_values_out);
      }
    }
  }

  if (common_count != 0U) {
    if (channel_mask == kLeafChannelPdf) {
      for (std::size_t i = 0; i < common_count; ++i) {
        const double x = scratch->lower[i];
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double sqt = std::sqrt(x);
        scratch->cdf_b[i] = -std::pow(a - k + x * l, 2.0) / (2.0 * x);
        scratch->survival_b[i] =
            -std::pow(a + k - x * l, 2.0) / (2.0 * x);
        scratch->times[i] = (-k + a) / sqt + sqt * l;
        scratch->onset_upper[i] = (k + a) / sqt - sqt * l;
      }
      const auto n = static_cast<int>(common_count);
      vvexp(scratch->cdf_b.data(), scratch->cdf_b.data(), &n);
      vvexp(scratch->survival_b.data(), scratch->survival_b.data(), &n);
      compiled_math_batch_normal_cdf_from_z(
          scratch->times.data(),
          scratch->pdf_b.data(),
          common_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->onset_upper.data(),
          scratch->pdf_b.data(),
          common_count);
      for (std::size_t i = 0; i < common_count; ++i) {
        const double x = scratch->lower[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double sqt = std::sqrt(x);
        const double t1 =
            kInvSqrtTwo *
            (scratch->cdf_b[i] - scratch->survival_b[i]) /
            (std::sqrt(M_PI) * sqt);
        const double t2 =
            0.5 * l *
            ((2.0 * scratch->times[i] - 1.0) +
             (2.0 * scratch->onset_upper[i] - 1.0));
        const double base_pdf =
            std::exp(std::log(std::max(1e-300, t1 + t2)) -
                     M_LN2 - std::log(a));
        compiled_math_direct_tile_multiply_finished_value_at(
            static_cast<std::size_t>(scratch->lanes_b[i]),
            base_pdf,
            0.0,
            scratch->survival_a[i],
            channel_mask,
            frame,
            &live_count,
            factor_values_out);
      }
    } else {
      for (std::size_t i = 0; i < common_count; ++i) {
        const double x = scratch->lower[i];
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double sqt = std::sqrt(x);
        scratch->cdf_b[i] = -0.5 * std::pow(k - a - x * l, 2.0) / x;
        scratch->survival_b[i] =
            -0.5 * std::pow(a + k - x * l, 2.0) / x;
        scratch->pdf_c[i] = 0.5 * (std::log(x) - M_LN2 - kLogPi);
        scratch->times[i] = -(k - a + x * l) / sqt;
        scratch->onset_upper[i] = -(k + a + x * l) / sqt;
        scratch->pdf_b[i] = (k + a) / sqt - sqt * l;
        scratch->cdf_c[i] = (k - a) / sqt - sqt * l;
      }
      const auto n = static_cast<int>(common_count);
      vvexp(scratch->cdf_b.data(), scratch->cdf_b.data(), &n);
      vvexp(scratch->survival_b.data(), scratch->survival_b.data(), &n);
      vvexp(scratch->pdf_c.data(), scratch->pdf_c.data(), &n);
      compiled_math_batch_normal_log_cdf_from_z(
          scratch->times.data(),
          scratch->survival_c.data(),
          common_count);
      compiled_math_batch_normal_log_cdf_from_z(
          scratch->onset_upper.data(),
          scratch->survival_c.data(),
          common_count);
      for (std::size_t i = 0; i < common_count; ++i) {
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        scratch->times[i] += 2.0 * l * (k - a);
        scratch->onset_upper[i] += 2.0 * l * (k + a);
      }
      vvexp(scratch->times.data(), scratch->times.data(), &n);
      vvexp(scratch->onset_upper.data(), scratch->onset_upper.data(), &n);
      compiled_math_batch_normal_cdf_from_z(
          scratch->pdf_b.data(),
          scratch->survival_c.data(),
          common_count);
      compiled_math_batch_normal_cdf_from_z(
          scratch->cdf_c.data(),
          scratch->survival_c.data(),
          common_count);
      for (std::size_t i = 0; i < common_count; ++i) {
        const double x = scratch->lower[i];
        const double k = scratch->upper[i];
        const double l = scratch->exact[i];
        const double a = scratch->pdf_a[i];
        const double t1 =
            scratch->pdf_c[i] *
            (scratch->cdf_b[i] - scratch->survival_b[i]);
        const double t2 =
            a + (scratch->onset_upper[i] - scratch->times[i]) /
                    (2.0 * l);
        const double t4a = 2.0 * scratch->pdf_b[i] - 1.0;
        const double t4b = 2.0 * scratch->cdf_c[i] - 1.0;
        const double t4 =
            0.5 * (x * l - a - k + 0.5 / l) * t4a +
            0.5 * (k - a - x * l - 0.5 / l) * t4b;
        const double base_cdf = 0.5 * (t4 + t2 + t1) / a;
        compiled_math_direct_tile_multiply_finished_value_at(
            static_cast<std::size_t>(scratch->lanes_b[i]),
            0.0,
            base_cdf,
            scratch->survival_a[i],
            channel_mask,
            frame,
            &live_count,
            factor_values_out);
      }
    }
  }
  *live_count_out = live_count;
  return true;
#else
	  (void)op_parent_offset;
	  (void)eligible_count;
	  (void)q_begin;
	  (void)rule;
	  (void)child_width;
  (void)tile_capacity;
  (void)eval_count;
  return false;
#endif
}

inline bool compiled_math_direct_tile_prepare_leaf_time_inputs(
    const CompiledMathProgram &program,
    const semantic::Index plan_id,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchSourceProductScratch *scratch,
    double *time_by_lane,
    const semantic::Index **active_lanes_out) {
  if (frame == nullptr || scratch == nullptr || time_by_lane == nullptr ||
      active_lanes_out == nullptr) {
    return false;
  }
  if (eval_lanes == nullptr && frame->source_lanes_a.size() < eval_count) {
    return false;
  }
  const semantic::Index *active_lanes = eval_lanes;
  if (active_lanes == nullptr) {
    for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
      frame->source_lanes_a[eval_pos] =
          static_cast<semantic::Index>(eval_pos);
    }
    active_lanes = frame->source_lanes_a.data();
  }
  const CompiledMathBoundVectorEvalContext time_context{
      nullptr,
      active_lanes,
      eval_count,
      frame,
      op_parent_offset,
      eligible_count,
      q_begin,
      &rule,
      child_width};
  if (!compiled_math_eval_scalar_bound_vector_plan_compact(
          program,
          plan_id,
          time_context,
          scratch->pdf_a.data())) {
    return false;
  }
  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
    const auto child_pos = static_cast<std::size_t>(active_lanes[eval_pos]);
    const auto parent_pos = child_pos / child_width;
    if (parent_pos >= eligible_count) {
      return false;
    }
    const auto op_parent_pos = op_parent_offset + parent_pos;
    const auto *loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
    if (loaded == nullptr) {
      return false;
    }
    time_by_lane[child_pos] = scratch->pdf_a[eval_pos];
    frame->source_leaf_inputs_by_lane[child_pos] = loaded;
  }
  *active_lanes_out = active_lanes;
  return true;
}

inline bool compiled_math_direct_tile_value_leaf_update(
    const CompiledMathProgram &program,
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
	    CompiledMathBatchIntegralFrame *frame,
	    CompiledMathBatchSourceProductScratch *scratch,
	    std::size_t *live_count_out,
	    double *factor_values_out) {
	  if (prepared.op == nullptr || prepared.channel == nullptr ||
	      frame == nullptr || scratch == nullptr || live_count_out == nullptr ||
	      factor_values_out == nullptr ||
	      prepared.channel->leaf_index == semantic::kInvalidIndex ||
	      prepared.tile_time_plan_id == semantic::kInvalidIndex) {
	    return false;
	  }
	  if (static_cast<leaf::DistKind>(prepared.channel->leaf_dist_kind) ==
	      leaf::DistKind::Lognormal) {
	    std::size_t value_live_count = 0U;
	    if (!compiled_math_direct_tile_lognormal_values(
	            prepared,
	            op_parent_offset,
	            eligible_count,
	            q_begin,
	            rule,
	            child_width,
	            tile_capacity,
	            eval_lanes,
	            eval_count,
	            frame,
	            scratch,
	            &value_live_count,
	            factor_values_out)) {
	      return false;
	    }
	    *live_count_out = value_live_count;
	    return true;
	  }
	  const semantic::Index *active_lanes = nullptr;
	  if (!compiled_math_direct_tile_prepare_leaf_time_inputs(
	          program,
	          prepared.tile_time_plan_id,
	          op_parent_offset,
	          eligible_count,
	          q_begin,
          rule,
	          child_width,
	          eval_lanes,
		          eval_count,
		          frame,
		          scratch,
		          scratch->times.data(),
		          &active_lanes) ||
	      active_lanes == nullptr) {
	    return false;
	  }
	  const CompiledMathLaneBatchState tile_leaf_state{
	      tile_capacity,
	      0U,
      BatchTimeSlotView{},
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      nullptr,
	      0U,
	      nullptr,
	      &frame->source_leaf_inputs_by_lane};
	  if (!compiled_math_batch_leaf_values_from_times(
	          prepared.channel->leaf_dist_kind,
	          prepared.channel->leaf_index,
          prepared.channel->leaf_onset_abs_value,
          tile_leaf_state,
	          active_lanes,
	          eval_count,
	          scratch->times.data(),
	          prepared.op->value_channel_mask,
	          factor_values_out,
	          scratch)) {
	    return false;
	  }
	  *live_count_out =
	      eval_lanes == nullptr
	          ? compiled_math_batch_array_store_live_values_dense_mask(
	                eval_count,
	                factor_values_out,
	                factor_values_out)
	          : compiled_math_batch_array_store_live_values_mask(
	                active_lanes,
	                eval_count,
	                factor_values_out,
	                factor_values_out);
	  return true;
	}

inline bool compiled_math_direct_tile_classify_conditioned_leaf_lanes(
    const CompiledMathProgram &program,
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
    CompiledMathBatchIntegralFrame *frame,
	    CompiledMathBatchSourceProductScratch *scratch,
	    std::size_t *live_count,
	    std::size_t *value_count,
	    double *factor_values_out) {
		  if (prepared.op == nullptr || frame == nullptr || scratch == nullptr ||
		      live_count == nullptr || value_count == nullptr ||
		      factor_values_out == nullptr) {
		    return false;
		  }
		  const auto channel_mask = prepared.op->value_channel_mask;
		  *live_count = 0U;
		  *value_count = 0U;
		  (void)program;
		  const semantic::Index *active_lanes = eval_lanes;
		  if (active_lanes == nullptr) {
		    if (frame->source_lanes_a.size() < eval_count) {
		      return false;
		    }
		    for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
		      frame->source_lanes_a[eval_pos] =
		          static_cast<semantic::Index>(eval_pos);
		    }
		    active_lanes = frame->source_lanes_a.data();
		  }
		  const bool needs_bind_time =
		      prepared.time_is_bind ||
		      (prepared.has_time_cap && prepared.cap_is_bind);
		  auto cached_parent = std::numeric_limits<std::size_t>::max();
		  const ExactLoadedLeafInput *loaded = nullptr;
		  double row_shift = 0.0;
		  double row_scale = 0.0;
		  double fixed_time = 0.0;
		  double fixed_cap = 0.0;
		  double lower = 0.0;
		  double upper = std::numeric_limits<double>::infinity();
		  double exact = 0.0;
		  std::uint8_t has_exact = 0U;
		  for (std::size_t eval_pos = 0; eval_pos < eval_count; ++eval_pos) {
	    const auto child_pos =
	        static_cast<std::size_t>(active_lanes[eval_pos]);
	    const auto parent_pos = child_pos / child_width;
	    if (parent_pos >= eligible_count) {
	      return false;
	    }
		    if (parent_pos != cached_parent) {
		      cached_parent = parent_pos;
		      const auto op_parent_pos = op_parent_offset + parent_pos;
		      loaded = frame->source_leaf_inputs_by_op_parent[op_parent_pos];
		      if (loaded == nullptr) {
		        return false;
		      }
		      if (needs_bind_time) {
		        row_shift = frame->direct_parent_shift[parent_pos];
		        row_scale = frame->direct_parent_scale[parent_pos];
		      }
		      fixed_time =
		          prepared.time_is_bind
		              ? 0.0
		              : frame->direct_source_time_by_op_parent[op_parent_pos];
		      fixed_cap =
		          prepared.has_time_cap && !prepared.cap_is_bind
		              ? frame->direct_source_cap_by_op_parent[op_parent_pos]
		              : 0.0;
		      lower = frame->direct_bound_lower_by_op_parent[op_parent_pos];
		      upper = frame->direct_bound_upper_by_op_parent[op_parent_pos];
		      exact = frame->direct_bound_exact_by_op_parent[op_parent_pos];
		      has_exact = frame->direct_bound_has_exact_by_op_parent[op_parent_pos];
		    }
		    const auto q_slot = child_pos - parent_pos * child_width;
		    const double bind_time =
		        needs_bind_time
		            ? row_shift + row_scale * rule.nodes[q_begin + q_slot]
		            : 0.0;
		    double current_time =
		        prepared.time_is_bind ? bind_time : fixed_time;
		    if (prepared.has_time_cap) {
		      current_time =
		          std::min(
		              current_time,
		              prepared.cap_is_bind ? bind_time : fixed_cap);
		    }
	    if (has_exact != 0U) {
	      const bool exact_in_bounds =
	          (!(lower > 0.0) || exact > lower) &&
	          (!std::isfinite(upper) || exact < upper);
	      compiled_math_direct_tile_multiply_live_value_at(
	          child_pos,
	          compiled_math_direct_tile_forced_channel_value(
	              exact_in_bounds && current_time >= exact,
	              channel_mask),
	          frame,
	          live_count,
	          factor_values_out);
	      continue;
	    }
	    if (std::isfinite(upper) && current_time >= upper) {
	      compiled_math_direct_tile_multiply_live_value_at(
	          child_pos,
          compiled_math_direct_tile_forced_channel_value(true, channel_mask),
          frame,
          live_count,
          factor_values_out);
      continue;
    }
    if (current_time <= lower) {
      compiled_math_direct_tile_multiply_live_value_at(
          child_pos,
          compiled_math_direct_tile_forced_channel_value(false, channel_mask),
          frame,
          live_count,
          factor_values_out);
	      continue;
	    }

		    scratch->times[child_pos] = current_time;
	    frame->source_leaf_inputs_by_lane[child_pos] = loaded;
	    scratch->lanes_a[(*value_count)++] =
	        static_cast<semantic::Index>(child_pos);
  }
  return true;
}

inline void compiled_math_direct_tile_finalize_conditioned_leaf_values(
    const std::size_t op_parent_offset,
    const std::size_t child_width,
    const std::uint8_t channel_mask,
    const std::size_t value_count,
	    CompiledMathBatchIntegralFrame *frame,
	    CompiledMathBatchSourceProductScratch *scratch,
	    std::size_t *live_count,
	    double *factor_values_out) {
  for (std::size_t i = 0; i < value_count; ++i) {
    const auto child_pos = static_cast<std::size_t>(scratch->lanes_a[i]);
    const auto parent_pos = child_pos / child_width;
    const auto op_parent_pos = op_parent_offset + parent_pos;
    compiled_math_direct_tile_multiply_live_value_at(
        child_pos,
        compiled_math_direct_tile_conditioned_continuous_value(
            op_parent_pos,
            frame->source_leaf_values[child_pos],
            channel_mask,
            *frame),
        frame,
        live_count,
        factor_values_out);
  }
}

inline bool compiled_math_direct_tile_conditioned_leaf_update(
    const CompiledMathProgram &program,
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
	    CompiledMathBatchIntegralFrame *frame,
	    CompiledMathBatchSourceProductScratch *scratch,
	    std::size_t *live_count_out,
	    double *factor_values_out) {
  if (prepared.op == nullptr || prepared.channel == nullptr ||
      frame == nullptr ||
      scratch == nullptr || live_count_out == nullptr ||
	      factor_values_out == nullptr ||
	      !prepared.conditioned_direct_leaf ||
      prepared.channel->leaf_index == semantic::kInvalidIndex ||
      prepared.tile_time_plan_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto channel_mask = prepared.op->value_channel_mask;
  const auto value_mask =
      channel_mask == kLeafChannelPdf ? kLeafChannelPdf : kLeafChannelCdf;
  std::size_t live_count = 0U;
  std::size_t value_count = 0U;
  if (!compiled_math_direct_tile_classify_conditioned_leaf_lanes(
          program,
          prepared,
          op_parent_offset,
          eligible_count,
          q_begin,
          rule,
          child_width,
          eval_lanes,
          eval_count,
          frame,
          scratch,
          &live_count,
          &value_count,
          factor_values_out)) {
    return false;
  }
  if (value_count == 0U) {
    *live_count_out = live_count;
    return true;
  }

  const CompiledMathLaneBatchState tile_leaf_state{
      tile_capacity,
      0U,
      BatchTimeSlotView{},
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      0U,
      nullptr,
      &frame->source_leaf_inputs_by_lane};
  if (!compiled_math_batch_leaf_values_from_times(
          prepared.channel->leaf_dist_kind,
          prepared.channel->leaf_index,
          prepared.channel->leaf_onset_abs_value,
          tile_leaf_state,
          scratch->lanes_a.data(),
          value_count,
          scratch->times.data(),
          value_mask,
          frame->source_leaf_values.data(),
          scratch)) {
    return false;
  }
  compiled_math_direct_tile_finalize_conditioned_leaf_values(
      op_parent_offset,
      child_width,
      channel_mask,
      value_count,
      frame,
      scratch,
      &live_count,
      factor_values_out);
  *live_count_out = live_count;
  return true;
}

inline bool compiled_math_direct_tile_leaf_update(
    const CompiledMathProgram &program,
    const CompiledMathDirectLeafPreparedOp &prepared,
    const std::size_t op_parent_offset,
    const std::size_t eligible_count,
    const std::size_t q_begin,
    const quadrature::Rule<quadrature::kDefaultFiniteOrder> &rule,
    const std::size_t child_width,
    const std::size_t tile_capacity,
    const semantic::Index *eval_lanes,
    const std::size_t eval_count,
	    CompiledMathBatchIntegralFrame *frame,
	    CompiledMathBatchSourceProductScratch *scratch,
	    std::size_t *live_count_out,
	    double *factor_values_out) {
	  if (prepared.channel == nullptr) {
    return false;
  }
  if (prepared.conditioned_direct_leaf) {
    return compiled_math_direct_tile_conditioned_leaf_update(
        program,
        prepared,
        op_parent_offset,
        eligible_count,
        q_begin,
        rule,
        child_width,
        tile_capacity,
        eval_lanes,
        eval_count,
        frame,
        scratch,
        live_count_out,
        factor_values_out);
  }
  return compiled_math_direct_tile_value_leaf_update(
      program,
      prepared,
      op_parent_offset,
      eligible_count,
      q_begin,
      rule,
      child_width,
      tile_capacity,
      eval_lanes,
      eval_count,
      frame,
      scratch,
      live_count_out,
      factor_values_out);
}

inline bool evaluate_compiled_integral_kernel_lane_batch_pdf_antiderivative(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane);

inline bool evaluate_compiled_integral_kernel_lane_batch_direct_leaf(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane);

inline bool evaluate_compiled_integral_kernel_lane_batch(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane);

inline bool compiled_math_lane_term_source_product_factor_values(
    const CompiledMathProgram &program,
    const CompiledMathLaneTermOp &op,
    const CompiledMathLaneBatchState &state,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    CompiledMathBatchIntegralFrame *frame,
    double *values_out) {
  if (frame == nullptr || values_out == nullptr) {
    return false;
  }
  if (lane_count == 0U) {
    return true;
  }
  if (lanes == nullptr || op.kind == CompiledMathLaneTermOpKind::IntegralFactor) {
    return false;
  }
  return compiled_math_batch_source_product_planned_values(
      program,
      op.source_product_channel_id,
      op.source_vector_program_id,
      op.factor_kind,
      op.value_channel_mask,
      op.fill_channel_mask,
      op.constant_value,
      op.cache_result,
      state,
      lanes,
      lane_count,
      frame,
      values_out,
      op.kind == CompiledMathLaneTermOpKind::SourceProgramFactor);
}

inline bool compiled_math_ensure_direct_leaf_parent_capacity(
    const std::size_t entry_count,
    const std::size_t parent_count,
    CompiledMathBatchIntegralFrame *frame) {
  if (frame == nullptr) {
    return false;
  }
  if (frame->direct_leaf_ops.size() < entry_count) {
    frame->direct_leaf_ops.resize(entry_count);
  }
  const auto entry_parent_count = entry_count * parent_count;
  if (frame->source_leaf_inputs_by_op_parent.size() < entry_parent_count) {
    frame->source_leaf_inputs_by_op_parent.resize(entry_parent_count, nullptr);
  }
  if (frame->direct_source_time_by_op_parent.size() < entry_parent_count) {
    frame->direct_source_time_by_op_parent.resize(entry_parent_count, 0.0);
  }
  if (frame->direct_source_cap_by_op_parent.size() < entry_parent_count) {
    frame->direct_source_cap_by_op_parent.resize(entry_parent_count, 0.0);
  }
  if (frame->direct_bound_lower_by_op_parent.size() < entry_parent_count) {
    frame->direct_bound_lower_by_op_parent.resize(entry_parent_count, 0.0);
  }
  if (frame->direct_bound_upper_by_op_parent.size() < entry_parent_count) {
    frame->direct_bound_upper_by_op_parent.resize(
        entry_parent_count,
        std::numeric_limits<double>::infinity());
  }
  if (frame->direct_bound_exact_by_op_parent.size() < entry_parent_count) {
    frame->direct_bound_exact_by_op_parent.resize(entry_parent_count, 0.0);
  }
  if (frame->direct_bound_has_exact_by_op_parent.size() <
      entry_parent_count) {
    frame->direct_bound_has_exact_by_op_parent.resize(entry_parent_count, 0U);
  }
  if (frame->direct_lower_cdf_by_op_parent.size() < entry_parent_count) {
    frame->direct_lower_cdf_by_op_parent.resize(entry_parent_count, 0.0);
  }
  if (frame->direct_lower_survival_by_op_parent.size() <
      entry_parent_count) {
    frame->direct_lower_survival_by_op_parent.resize(
        entry_parent_count,
        1.0);
  }
  if (frame->direct_upper_cdf_by_op_parent.size() < entry_parent_count) {
    frame->direct_upper_cdf_by_op_parent.resize(entry_parent_count, 0.0);
  }
  return true;
}

template <typename DirectLeafPlan>
inline bool compiled_math_prepare_direct_leaf_parent_values(
    const CompiledMathProgram &program,
    const DirectLeafPlan &plan,
    const CompiledMathLaneBatchState &parent_state,
    const semantic::Index *parent_lanes,
    const std::size_t parent_count,
    const std::size_t entry_parent_offset,
    const std::size_t entry_parent_capacity,
    CompiledMathBatchIntegralFrame *frame) {
  if (frame == nullptr || (parent_count != 0U && parent_lanes == nullptr) ||
      plan.direct_leaf_index == semantic::kInvalidIndex ||
      (!plan.direct_leaf_time_is_bind &&
       plan.direct_leaf_parent_time_plan_id == semantic::kInvalidIndex) ||
      (plan.direct_leaf_has_time_cap && !plan.direct_leaf_cap_is_bind &&
       plan.direct_leaf_parent_cap_plan_id == semantic::kInvalidIndex)) {
    return false;
  }
  const bool conditioned = plan.direct_leaf_conditioned;
  if (conditioned &&
      plan.direct_leaf_source_id == semantic::kInvalidIndex) {
    return false;
  }
  if ((!plan.direct_leaf_time_is_bind &&
       plan.direct_leaf_parent_time_plan_id == semantic::kInvalidIndex) ||
      (plan.direct_leaf_has_time_cap &&
       !plan.direct_leaf_cap_is_bind &&
       plan.direct_leaf_parent_cap_plan_id == semantic::kInvalidIndex) ||
      (conditioned &&
       plan.direct_leaf_condition_bounds_plan_id == semantic::kInvalidIndex)) {
    return false;
  }
  std::size_t lower_bound_count = 0U;
  std::size_t upper_bound_count = 0U;
  semantic::Index parent_time_id = semantic::kInvalidIndex;
  semantic::Index parent_cap_id = semantic::kInvalidIndex;
  if (!plan.direct_leaf_time_is_bind) {
    CompiledMathIndexSpan parent_time_span;
    if (!compiled_math_bound_vector_plan_span(
            program,
            plan.direct_leaf_parent_time_plan_id,
            &parent_time_span) ||
        parent_time_span.size != 1) {
      return false;
    }
    const auto &time_op =
        program.bound_vector_ops[
            static_cast<std::size_t>(parent_time_span.offset)];
    if (time_op.kind != CompiledMathBoundVectorOpKind::StateTimeSlot ||
        time_op.time_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(time_op.time_id) >=
            parent_state.time_slot_count) {
      return false;
    }
    parent_time_id = time_op.time_id;
  }
  if (plan.direct_leaf_has_time_cap && !plan.direct_leaf_cap_is_bind) {
    CompiledMathIndexSpan parent_cap_span;
    if (!compiled_math_bound_vector_plan_span(
            program,
            plan.direct_leaf_parent_cap_plan_id,
            &parent_cap_span) ||
        parent_cap_span.size != 1) {
      return false;
    }
    const auto &cap_op =
        program.bound_vector_ops[
            static_cast<std::size_t>(parent_cap_span.offset)];
    if (cap_op.kind != CompiledMathBoundVectorOpKind::StateTimeSlot ||
        cap_op.time_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(cap_op.time_id) >=
            parent_state.time_slot_count) {
      return false;
    }
    parent_cap_id = cap_op.time_id;
  }
  for (std::size_t parent_pos = 0; parent_pos < parent_count; ++parent_pos) {
    auto *source_channels =
        frame->direct_source_channels_by_parent[parent_pos];
    if (source_channels == nullptr) {
      return false;
    }
    const auto entry_parent_pos = entry_parent_offset + parent_pos;
    if (!conditioned &&
        !source_channels->source_product_direct_leaf_available(
            plan.source_product_channel_id)) {
      return false;
    }
    frame->source_leaf_inputs_by_op_parent[entry_parent_pos] =
        &source_channels->source_product_leaf_input(plan.direct_leaf_index);
    const auto parent_lane = parent_lanes[parent_pos];
    frame->direct_source_time_by_op_parent[entry_parent_pos] =
        parent_time_id == semantic::kInvalidIndex
            ? 0.0
            : parent_state.time(parent_time_id, parent_lane);
    frame->direct_source_cap_by_op_parent[entry_parent_pos] =
        parent_cap_id == semantic::kInvalidIndex
            ? 0.0
            : parent_state.time(parent_cap_id, parent_lane);
    if (!conditioned) {
      continue;
    }
    frame->direct_lower_cdf_by_op_parent[entry_parent_pos] = 0.0;
    frame->direct_lower_survival_by_op_parent[entry_parent_pos] = 1.0;
    frame->direct_upper_cdf_by_op_parent[entry_parent_pos] = 0.0;
  }
  if (!conditioned) {
    return true;
  }
  if (!compiled_math_eval_condition_bound_vector_plan_compact(
          program,
          plan.direct_leaf_condition_bounds_plan_id,
          parent_state,
          parent_lanes,
          parent_count,
          frame,
          frame->direct_bound_lower_by_op_parent.data() + entry_parent_offset,
          frame->direct_bound_upper_by_op_parent.data() + entry_parent_offset,
          frame->direct_bound_exact_by_op_parent.data() + entry_parent_offset,
          frame->direct_bound_has_exact_by_op_parent.data() +
              entry_parent_offset)) {
    return false;
  }

  for (std::size_t parent_pos = 0; parent_pos < parent_count; ++parent_pos) {
    const auto entry_parent_pos = entry_parent_offset + parent_pos;
    frame->source_leaf_inputs_by_lane[entry_parent_pos] =
        frame->source_leaf_inputs_by_op_parent[entry_parent_pos];
    const double lower =
        frame->direct_bound_lower_by_op_parent[entry_parent_pos];
    if (lower > 0.0) {
      frame->lower_by_lane[entry_parent_pos] = lower;
      frame->source_lanes_a[lower_bound_count++] =
          static_cast<semantic::Index>(entry_parent_pos);
    }
    const double upper =
        frame->direct_bound_upper_by_op_parent[entry_parent_pos];
    if (std::isfinite(upper)) {
      frame->upper_by_lane[entry_parent_pos] = upper;
      frame->source_lanes_b[upper_bound_count++] =
          static_cast<semantic::Index>(entry_parent_pos);
    }
  }

  if (lower_bound_count == 0U && upper_bound_count == 0U) {
    return true;
  }
  auto &bound_scratch =
      frame->source_product_scratch_layer(
          0U,
          entry_parent_capacity,
          parent_count);
  const CompiledMathLaneBatchState bound_leaf_state{
      entry_parent_capacity,
      0U,
      BatchTimeSlotView{},
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      0U,
      nullptr,
      &frame->source_leaf_inputs_by_lane};
  if (lower_bound_count != 0U &&
      !compiled_math_batch_leaf_values_from_times(
          plan.direct_leaf_dist_kind,
          plan.direct_leaf_index,
          plan.direct_leaf_onset_abs_value,
          bound_leaf_state,
          frame->source_lanes_a.data(),
          lower_bound_count,
          frame->lower_by_lane.data(),
          kLeafChannelCdf,
          frame->direct_lower_cdf_by_op_parent.data(),
          &bound_scratch)) {
    return false;
  }
  for (std::size_t i = 0; i < lower_bound_count; ++i) {
    const auto entry_parent_pos =
        static_cast<std::size_t>(frame->source_lanes_a[i]);
    frame->direct_lower_survival_by_op_parent[entry_parent_pos] =
        clamp_probability(
            1.0 - frame->direct_lower_cdf_by_op_parent[entry_parent_pos]);
  }
  if (upper_bound_count != 0U &&
      !compiled_math_batch_leaf_values_from_times(
          plan.direct_leaf_dist_kind,
          plan.direct_leaf_index,
          plan.direct_leaf_onset_abs_value,
          bound_leaf_state,
          frame->source_lanes_b.data(),
          upper_bound_count,
          frame->upper_by_lane.data(),
          kLeafChannelCdf,
          frame->direct_upper_cdf_by_op_parent.data(),
          &bound_scratch)) {
    return false;
  }
  return true;
}

inline bool compiled_math_prepare_direct_tile_block_factors(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const std::size_t eligible_count,
    CompiledMathBatchIntegralFrame *frame,
    std::size_t *factor_count_out) {
  if (frame == nullptr || factor_count_out == nullptr) {
    return false;
  }
  const auto factor_count =
      static_cast<std::size_t>(kernel.source_product_block_factors.size);
  *factor_count_out = factor_count;
  const auto factor_parent_count = factor_count * eligible_count;
  if (!compiled_math_ensure_direct_leaf_parent_capacity(
          factor_count,
          eligible_count,
          frame)) {
    return false;
  }

  const auto factor_offset =
      static_cast<std::size_t>(kernel.source_product_block_factors.offset);
  for (std::size_t factor_pos = 0; factor_pos < factor_count; ++factor_pos) {
    const auto &factor =
        program.integral_kernel_source_product_block_factors[
            factor_offset + factor_pos];
    auto &prepared = frame->direct_leaf_ops[factor_pos];
    prepared.op = nullptr;
    prepared.channel = nullptr;
    prepared.tile_time_plan_id = semantic::kInvalidIndex;
    prepared.time_is_bind = false;
    prepared.has_time_cap = false;
    prepared.cap_is_bind = false;
    prepared.conditioned_direct_leaf = false;
    if (factor.factor_kind ==
        CompiledMathSourceProductBlockFactorKind::IntegralFactor) {
      continue;
    }
    if (factor.factor_kind ==
        CompiledMathSourceProductBlockFactorKind::SourceProgram) {
      continue;
    }
    if (factor.value_channel_mask == 0U) {
      if (factor.op_id != semantic::kInvalidIndex &&
          static_cast<std::size_t>(factor.op_id) <
              program.integral_kernel_source_product_ops.size()) {
        prepared.op =
            &program.integral_kernel_source_product_ops[
                static_cast<std::size_t>(factor.op_id)];
      }
      continue;
    }
    if (factor.op_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(factor.op_id) >=
            program.integral_kernel_source_product_ops.size() ||
        factor.source_product_channel_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(factor.source_product_channel_id) >=
            program.integral_kernel_source_product_channels.size()) {
      return false;
    }
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(factor.op_id)];
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(factor.source_product_channel_id)];
    prepared.op = &op;
    prepared.channel = &channel;
    prepared.conditioned_direct_leaf = factor.direct_leaf_conditioned;
    if (factor.direct_leaf_index == semantic::kInvalidIndex) {
      return false;
    }
    prepared.time_is_bind = factor.direct_leaf_time_is_bind;
    prepared.has_time_cap = factor.direct_leaf_has_time_cap;
    prepared.cap_is_bind = factor.direct_leaf_cap_is_bind;
    prepared.tile_time_plan_id = factor.direct_leaf_tile_time_plan_id;
    if (!compiled_math_prepare_direct_leaf_parent_values(
            program,
            factor,
            parent_state,
            frame->parent_lanes.data(),
            eligible_count,
            factor_pos * eligible_count,
            factor_parent_count,
            frame)) {
      return false;
    }
  }
  return true;
}

inline void compiled_math_reset_block_factor_cache(
    const std::size_t factor_count,
    const std::size_t tile_capacity,
    CompiledMathBatchIntegralFrame *frame) {
  if (frame == nullptr) {
    return;
  }
  const auto cache_size = factor_count * tile_capacity;
  if (frame->source_product_block_factor_values.size() < cache_size) {
    frame->source_product_block_factor_values.resize(cache_size, 0.0);
  }
  if (frame->source_product_block_factor_valid.size() < cache_size) {
    frame->source_product_block_factor_valid.resize(cache_size, 0U);
  }
  std::fill(
      frame->source_product_block_factor_valid.begin(),
      frame->source_product_block_factor_valid.begin() +
          static_cast<std::ptrdiff_t>(cache_size),
      0U);
}

struct CompiledMathLaneTermRunState {
  semantic::Index *current_lanes{nullptr};
  semantic::Index *next_lanes{nullptr};
  std::size_t term_count{0U};
};

inline bool compiled_math_apply_lane_term_constant_factor(
    const CompiledMathLaneTermOp &op,
    CompiledMathLaneTermRunState *state,
    CompiledMathBatchIntegralFrame *frame) {
  if (state == nullptr || state->current_lanes == nullptr ||
      state->next_lanes == nullptr || frame == nullptr) {
    return false;
  }
  const auto next_count =
      compiled_math_batch_array_multiply_scalar_compact(
          state->current_lanes,
          state->term_count,
          op.constant_value,
          frame->products.data(),
          state->next_lanes);
  std::swap(state->current_lanes, state->next_lanes);
  state->term_count = next_count;
  return true;
}

template <typename FillValues>
inline bool compiled_math_apply_lane_term_cached_factor_values(
    const CompiledMathLaneTermOp &op,
    CompiledMathLaneTermRunState *state,
    CompiledMathBatchIntegralFrame *frame,
    const std::size_t factor_count,
    const std::size_t factor_capacity,
    FillValues fill_values) {
  if (frame == nullptr || state == nullptr ||
      state->current_lanes == nullptr || state->next_lanes == nullptr ||
      op.factor_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(op.factor_id) >= factor_count) {
    return false;
  }
  const auto factor_pos = static_cast<std::size_t>(op.factor_id);
  if (op.use_count > 1) {
    auto *factor_values =
        frame->source_product_block_factor_values.data() +
        factor_pos * factor_capacity;
    auto *factor_valid =
        frame->source_product_block_factor_valid.data() +
        factor_pos * factor_capacity;
    std::size_t missing_count = 0U;
    for (std::size_t i = 0; i < state->term_count; ++i) {
      const auto lane_pos = static_cast<std::size_t>(state->current_lanes[i]);
      if (factor_valid[lane_pos] == 0U) {
        frame->source_lanes_a[missing_count++] =
            static_cast<semantic::Index>(lane_pos);
      }
    }
    if (missing_count != 0U) {
      if (!fill_values(
              frame->source_lanes_a.data(),
              missing_count,
              factor_values)) {
        return false;
      }
      for (std::size_t i = 0; i < missing_count; ++i) {
        factor_valid[static_cast<std::size_t>(frame->source_lanes_a[i])] = 1U;
      }
    }
    const auto next_count =
        compiled_math_batch_array_multiply_values_compact(
            state->current_lanes,
            state->term_count,
            factor_values,
            frame->products.data(),
            state->next_lanes);
    std::swap(state->current_lanes, state->next_lanes);
    state->term_count = next_count;
    return true;
  }
  if (!fill_values(
          state->current_lanes,
          state->term_count,
          frame->factor_values.data())) {
    return false;
  }
  const auto next_count =
      compiled_math_batch_array_multiply_values_compact(
          state->current_lanes,
          state->term_count,
          frame->factor_values.data(),
          frame->products.data(),
          state->next_lanes);
  std::swap(state->current_lanes, state->next_lanes);
  state->term_count = next_count;
  return true;
}

inline bool compiled_math_apply_lane_term_cached_integral_factor_op(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneTermOp &op,
    const CompiledMathLaneBatchState &child_state,
    CompiledMathLaneTermRunState *state,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t factor_count,
    const std::size_t factor_capacity,
    const std::size_t depth) {
  if (frame == nullptr || batch_workspace == nullptr || state == nullptr ||
      state->current_lanes == nullptr || state->next_lanes == nullptr ||
      op.factor_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(op.factor_id) >= factor_count) {
    return false;
  }
  const auto factor_pos = static_cast<std::size_t>(op.factor_id);
  auto materialize_values =
      [&](const semantic::Index *lanes,
          const std::size_t lane_count,
          double *values_out) -> bool {
    if (values_out == nullptr) {
      return false;
    }
    if (lane_count == 0U) {
      return true;
    }
    if (lanes == nullptr ||
        op.kind != CompiledMathLaneTermOpKind::IntegralFactor ||
        op.integral_factor_cache_slot == semantic::kInvalidIndex ||
        static_cast<std::size_t>(op.integral_factor_cache_slot) >=
            static_cast<std::size_t>(
                kernel.cached_integral_factor_plans.size)) {
      return false;
    }
    const auto plan_pos =
        static_cast<std::size_t>(kernel.cached_integral_factor_plans.offset) +
        static_cast<std::size_t>(op.integral_factor_cache_slot);
    if (plan_pos >= program.integral_kernel_cached_integral_factor_plans.size()) {
      return false;
    }
    const auto &factor_plan =
        program.integral_kernel_cached_integral_factor_plans[plan_pos];
    if (factor_plan.node_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(factor_plan.node_id) >= program.nodes.size() ||
        factor_plan.kernel_slot == semantic::kInvalidIndex ||
        static_cast<std::size_t>(factor_plan.kernel_slot) >=
            program.integral_kernels.size() ||
        factor_plan.lower_bound_plan_id == semantic::kInvalidIndex ||
        factor_plan.upper_bound_plan_id == semantic::kInvalidIndex ||
        frame->active_lower_bounds.size() < lane_count ||
        frame->active_upper_bounds.size() < lane_count) {
      return false;
    }
    const CompiledMathBoundVectorEvalContext bound_context{
        &child_state,
        lanes,
        lane_count,
        nullptr,
        0U,
        0U,
        0U,
        nullptr,
        0U};
    if (!compiled_math_eval_scalar_bound_vector_plan_compact(
            program,
            factor_plan.lower_bound_plan_id,
            bound_context,
            frame->active_lower_bounds.data()) ||
        !compiled_math_eval_scalar_bound_vector_plan_compact(
            program,
            factor_plan.upper_bound_plan_id,
            bound_context,
            frame->active_upper_bounds.data())) {
      return false;
    }
    const CompiledMathActiveBoundBatch factor_bounds{
        lanes,
        lane_count,
        frame->active_lower_bounds.data(),
        frame->active_upper_bounds.data()};
    if (!evaluate_compiled_integral_kernel_lane_batch(
            program,
            program.nodes[static_cast<std::size_t>(factor_plan.node_id)],
            program.integral_kernels[
                static_cast<std::size_t>(factor_plan.kernel_slot)],
            child_state,
            factor_bounds,
            batch_workspace,
            depth + 1U,
            values_out)) {
      return false;
    }
    if (factor_plan.clamp_probability) {
      for (std::size_t i = 0; i < lane_count; ++i) {
        const auto lane_pos = static_cast<std::size_t>(lanes[i]);
        values_out[lane_pos] = clamp_probability(values_out[lane_pos]);
      }
    }
    return true;
  };
  if (op.use_count > 1) {
    auto *factor_values =
        frame->source_product_block_factor_values.data() +
        factor_pos * factor_capacity;
    auto *factor_valid =
        frame->source_product_block_factor_valid.data() +
        factor_pos * factor_capacity;
    std::size_t missing_count = 0U;
    for (std::size_t i = 0; i < state->term_count; ++i) {
      const auto lane_pos = static_cast<std::size_t>(state->current_lanes[i]);
      if (factor_valid[lane_pos] == 0U) {
        frame->source_lanes_a[missing_count++] =
            static_cast<semantic::Index>(lane_pos);
      }
    }
    if (missing_count != 0U) {
      if (!materialize_values(
              frame->source_lanes_a.data(),
              missing_count,
              factor_values)) {
        return false;
      }
      for (std::size_t i = 0; i < missing_count; ++i) {
        factor_valid[static_cast<std::size_t>(frame->source_lanes_a[i])] = 1U;
      }
    }
    const auto next_count =
        compiled_math_batch_array_multiply_values_compact(
            state->current_lanes,
            state->term_count,
            factor_values,
            frame->products.data(),
            state->next_lanes);
    std::swap(state->current_lanes, state->next_lanes);
    state->term_count = next_count;
    return true;
  }
  if (!materialize_values(
          state->current_lanes,
          state->term_count,
          frame->factor_values.data())) {
    return false;
  }
  const auto next_count =
      compiled_math_batch_array_multiply_values_compact(
          state->current_lanes,
          state->term_count,
          frame->factor_values.data(),
          frame->products.data(),
          state->next_lanes);
  std::swap(state->current_lanes, state->next_lanes);
  state->term_count = next_count;
  return true;
}

template <typename Domain>
inline bool compiled_math_execute_lane_term_schedules(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan schedule_span,
    Domain *domain) {
  if (domain == nullptr) {
    return false;
  }
  for (semantic::Index schedule_idx = 0;
       schedule_idx < schedule_span.size;
       ++schedule_idx) {
    const auto schedule_pos =
        static_cast<std::size_t>(schedule_span.offset + schedule_idx);
    if (schedule_pos >= program.integral_kernel_lane_term_schedules.size()) {
      return false;
    }
    const auto &schedule =
        program.integral_kernel_lane_term_schedules[schedule_pos];
    CompiledMathLaneTermRunState state;
    if (!domain->seed_schedule(schedule, &state)) {
      return false;
    }
    bool factors_started = false;
    for (semantic::Index op_idx = 0; op_idx < schedule.ops.size; ++op_idx) {
      if (state.term_count == 0U) {
        break;
      }
      const auto op_pos =
          static_cast<std::size_t>(schedule.ops.offset + op_idx);
      if (op_pos >= program.integral_kernel_lane_term_ops.size()) {
        return false;
      }
      const auto &op = program.integral_kernel_lane_term_ops[op_pos];
      switch (op.kind) {
      case CompiledMathLaneTermOpKind::ExprUpper:
        if (!domain->apply_expr_upper_op(op, &state)) {
          return false;
        }
        break;
      case CompiledMathLaneTermOpKind::ConstantFactor:
        if (!factors_started) {
          if (!domain->begin_schedule(schedule)) {
            return false;
          }
          factors_started = true;
        }
        if (!domain->apply_constant_factor_op(op, &state)) {
          return false;
        }
        break;
      case CompiledMathLaneTermOpKind::DirectLeafFactor:
      case CompiledMathLaneTermOpKind::ConditionedDirectLeafFactor:
        if (!factors_started) {
          if (!domain->begin_schedule(schedule)) {
            return false;
          }
          factors_started = true;
        }
        if (!domain->apply_direct_leaf_factor_op(op, &state)) {
          return false;
        }
        break;
      case CompiledMathLaneTermOpKind::IntegralFactor:
        if (!factors_started) {
          if (!domain->begin_schedule(schedule)) {
            return false;
          }
          factors_started = true;
        }
        if (!domain->apply_integral_factor_op(op, &state)) {
          return false;
        }
        break;
      case CompiledMathLaneTermOpKind::SourceProgramFactor:
        if (!factors_started) {
          if (!domain->begin_schedule(schedule)) {
            return false;
          }
          factors_started = true;
        }
        if (!domain->apply_source_product_factor_op(op, &state)) {
          return false;
        }
        break;
      }
    }
    if (state.term_count != 0U &&
        !domain->accumulate_term(state.current_lanes, state.term_count)) {
      return false;
    }
  }
  return domain->finish_terms();
}

struct CompiledMathDirectTileBlockTermDomain {
  const CompiledMathProgram *program{nullptr};
  const CompiledMathIntegralKernel *kernel{nullptr};
  const CompiledMathLaneBatchState *parent_state{nullptr};
  const CompiledMathLaneBatchState *integral_child_state{nullptr};
  CompiledMathBatchIntegralFrame *frame{nullptr};
  CompiledMathBatchWorkspace *batch_workspace{nullptr};
  const quadrature::Rule<quadrature::kDefaultFiniteOrder> *rule{nullptr};
  std::size_t q_begin{0U};
  std::size_t child_width{0U};
  std::size_t eligible_count{0U};
  std::size_t tile_capacity{0U};
  std::size_t prepared_factor_count{0U};
  std::size_t depth{0U};
  double *out_by_parent_lane{nullptr};
  CompiledMathBatchSourceProductScratch *direct_scratch{nullptr};
  bool has_accumulated_terms{false};

  bool seed_schedule(
      const CompiledMathLaneTermSchedule &schedule,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || kernel == nullptr || parent_state == nullptr ||
        frame == nullptr || rule == nullptr || state == nullptr) {
      return false;
    }
    std::size_t count = 0U;
    switch (schedule.gate_kind) {
    case CompiledMathLaneTermGateKind::None:
      for (std::size_t parent_pos = 0; parent_pos < eligible_count;
           ++parent_pos) {
        for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
          const auto child_pos = parent_pos * child_width + q_slot;
          frame->products[child_pos] = schedule.sign;
          frame->term_lanes_a[count++] =
              static_cast<semantic::Index>(child_pos);
        }
      }
      break;
    case CompiledMathLaneTermGateKind::Outcome:
      for (std::size_t parent_pos = 0; parent_pos < eligible_count;
           ++parent_pos) {
        const auto parent_lane = frame->parent_lanes[parent_pos];
        if (!compiled_math_batch_outcome_gates_open(
                *program,
                schedule.outcome_gate_nodes,
                *parent_state,
                parent_lane)) {
          for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
            frame->products[parent_pos * child_width + q_slot] = 0.0;
          }
          continue;
        }
        for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
          const auto child_pos = parent_pos * child_width + q_slot;
          frame->products[child_pos] = schedule.sign;
          frame->term_lanes_a[count++] =
              static_cast<semantic::Index>(child_pos);
        }
      }
      break;
    case CompiledMathLaneTermGateKind::Time:
      for (std::size_t parent_pos = 0; parent_pos < eligible_count;
           ++parent_pos) {
        const auto parent_lane = frame->parent_lanes[parent_pos];
        for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
          const auto child_pos = parent_pos * child_width + q_slot;
          bool gates_open = true;
          if (!compiled_math_direct_tile_time_gates_open(
                  *program,
                  *kernel,
                  schedule.time_gate_nodes,
                  *parent_state,
                  *frame,
                  parent_pos,
                  parent_lane,
                  q_begin + q_slot,
                  *rule,
                  &gates_open)) {
            return false;
          }
          if (!gates_open) {
            frame->products[child_pos] = 0.0;
            continue;
          }
          frame->products[child_pos] = schedule.sign;
          frame->term_lanes_a[count++] =
              static_cast<semantic::Index>(child_pos);
        }
      }
      break;
    case CompiledMathLaneTermGateKind::OutcomeAndTime:
      for (std::size_t parent_pos = 0; parent_pos < eligible_count;
           ++parent_pos) {
        const auto parent_lane = frame->parent_lanes[parent_pos];
        if (!compiled_math_batch_outcome_gates_open(
                *program,
                schedule.outcome_gate_nodes,
                *parent_state,
                parent_lane)) {
          for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
            frame->products[parent_pos * child_width + q_slot] = 0.0;
          }
          continue;
        }
        for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
          const auto child_pos = parent_pos * child_width + q_slot;
          bool gates_open = true;
          if (!compiled_math_direct_tile_time_gates_open(
                  *program,
                  *kernel,
                  schedule.time_gate_nodes,
                  *parent_state,
                  *frame,
                  parent_pos,
                  parent_lane,
                  q_begin + q_slot,
                  *rule,
                  &gates_open)) {
            return false;
          }
          if (!gates_open) {
            frame->products[child_pos] = 0.0;
            continue;
          }
          frame->products[child_pos] = schedule.sign;
          frame->term_lanes_a[count++] =
              static_cast<semantic::Index>(child_pos);
        }
      }
      break;
    }
    state->current_lanes = frame->term_lanes_a.data();
    state->next_lanes = frame->term_lanes_b.data();
    state->term_count = count;
    return true;
  }

  bool begin_schedule(const CompiledMathLaneTermSchedule &) {
    return frame != nullptr;
  }

  bool apply_expr_upper_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || kernel == nullptr || parent_state == nullptr ||
        frame == nullptr || rule == nullptr || state == nullptr ||
        state->current_lanes == nullptr || state->next_lanes == nullptr) {
      return false;
    }
    std::size_t next_count = 0U;
    if (!compiled_math_direct_tile_apply_expr_upper_op(
            *program,
            *kernel,
            op,
            *parent_state,
            state->current_lanes,
            state->term_count,
            child_width,
            eligible_count,
            q_begin,
            *rule,
            frame,
            state->next_lanes,
            &next_count)) {
      return false;
    }
    std::swap(state->current_lanes, state->next_lanes);
    state->term_count = next_count;
    return true;
  }

  bool apply_constant_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    return compiled_math_apply_lane_term_constant_factor(op, state, frame);
  }

  bool apply_integral_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || integral_child_state == nullptr ||
        frame == nullptr || batch_workspace == nullptr) {
      return false;
    }
    return compiled_math_apply_lane_term_cached_integral_factor_op(
        *program,
        *kernel,
        op,
        *integral_child_state,
        state,
        frame,
        batch_workspace,
        prepared_factor_count,
        tile_capacity,
        depth);
  }

  bool apply_source_product_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || integral_child_state == nullptr ||
        frame == nullptr) {
      return false;
    }
    return compiled_math_apply_lane_term_cached_factor_values(
        op,
        state,
        frame,
        prepared_factor_count,
        tile_capacity,
        [&](const semantic::Index *lanes,
            const std::size_t lane_count,
            double *values_out) {
          return compiled_math_lane_term_source_product_factor_values(
              *program,
              op,
              *integral_child_state,
              lanes,
              lane_count,
              frame,
              values_out);
        });
  }

  bool apply_direct_leaf_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    if (op.factor_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(op.factor_id) >= prepared_factor_count) {
      return false;
    }
    const auto factor_pos = static_cast<std::size_t>(op.factor_id);
    const auto &prepared = frame->direct_leaf_ops[factor_pos];
    const auto *prepared_op = prepared.op;
    if (prepared_op == nullptr || program == nullptr || kernel == nullptr ||
        rule == nullptr ||
        state == nullptr || state->current_lanes == nullptr ||
        state->next_lanes == nullptr) {
      return false;
    }
	    if (direct_scratch == nullptr) {
	      direct_scratch =
	          &frame->source_product_scratch_layer(0U, tile_capacity, tile_capacity);
	    }
	    return compiled_math_apply_lane_term_cached_factor_values(
	        op,
	        state,
	        frame,
	        prepared_factor_count,
	        tile_capacity,
	        [&](const semantic::Index *lanes,
	            const std::size_t lane_count,
	            double *values_out) {
	          std::size_t factor_live_count = 0U;
	          const bool ok = compiled_math_direct_tile_leaf_update(
	              *program,
	              prepared,
	              factor_pos * eligible_count,
	              eligible_count,
	              q_begin,
	              *rule,
	              child_width,
	              tile_capacity,
	              lanes,
	              lane_count,
	              frame,
	              direct_scratch,
	              &factor_live_count,
	              values_out);
	          (void)factor_live_count;
	          return ok;
	        });
	  }

  bool accumulate_term(
      const semantic::Index *current_lanes,
      const std::size_t term_count) {
    if (kernel == nullptr || frame == nullptr || rule == nullptr ||
        out_by_parent_lane == nullptr || current_lanes == nullptr) {
      return false;
    }
    if (kernel->clean_signed_source_sum) {
      has_accumulated_terms = true;
      for (std::size_t i = 0; i < term_count; ++i) {
        const auto child_pos = static_cast<std::size_t>(current_lanes[i]);
        frame->values[child_pos] += frame->products[child_pos];
      }
    } else {
      compiled_math_direct_tile_weighted_accumulate_parent_rows(
          frame->parent_lanes.data(),
          eligible_count,
          child_width,
          q_begin,
          *rule,
          frame->direct_parent_scale.data(),
          frame->products.data(),
          out_by_parent_lane);
    }
    return true;
  }

  bool finish_terms() {
    if (kernel == nullptr || frame == nullptr || rule == nullptr ||
        out_by_parent_lane == nullptr) {
      return false;
    }
    if (!kernel->clean_signed_source_sum) {
      return true;
    }
    if (!has_accumulated_terms) {
      return true;
    }
    for (std::size_t child_pos = 0; child_pos < tile_capacity; ++child_pos) {
      frame->values[child_pos] = clean_signed_value(frame->values[child_pos]);
    }
    compiled_math_direct_tile_weighted_accumulate_parent_rows(
        frame->parent_lanes.data(),
        eligible_count,
        child_width,
        q_begin,
        *rule,
        frame->direct_parent_scale.data(),
        frame->values.data(),
        out_by_parent_lane);
    return true;
  }
};

struct CompiledMathChildLaneTermDomain {
  const CompiledMathProgram *program{nullptr};
  const CompiledMathIntegralKernel *kernel{nullptr};
  const CompiledMathLaneBatchState *child_state{nullptr};
  CompiledMathBatchIntegralFrame *frame{nullptr};
  CompiledMathBatchWorkspace *batch_workspace{nullptr};
  std::size_t child_count{0U};
  std::size_t child_capacity{0U};
  std::size_t depth{0U};
  bool has_accumulated_terms{false};

  bool seed_schedule(
      const CompiledMathLaneTermSchedule &schedule,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || child_state == nullptr || frame == nullptr ||
        state == nullptr) {
      return false;
    }
    std::size_t count = 0U;
    switch (schedule.gate_kind) {
    case CompiledMathLaneTermGateKind::None:
      for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
        const auto lane = frame->active_lanes[child_pos];
        frame->term_lanes_a[count++] = lane;
        frame->products[static_cast<std::size_t>(lane)] = schedule.sign;
      }
      break;
    case CompiledMathLaneTermGateKind::Outcome:
      for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
        const auto lane = frame->active_lanes[child_pos];
        if (!compiled_math_batch_outcome_gates_open(
                *program,
                schedule.outcome_gate_nodes,
                *child_state,
                lane)) {
          continue;
        }
        frame->term_lanes_a[count++] = lane;
        frame->products[static_cast<std::size_t>(lane)] = schedule.sign;
      }
      break;
    case CompiledMathLaneTermGateKind::Time:
      for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
        const auto lane = frame->active_lanes[child_pos];
        if (!compiled_math_batch_time_gates_open(
                *program,
                schedule.time_gate_nodes,
                *child_state,
                lane)) {
          continue;
        }
        frame->term_lanes_a[count++] = lane;
        frame->products[static_cast<std::size_t>(lane)] = schedule.sign;
      }
      break;
    case CompiledMathLaneTermGateKind::OutcomeAndTime:
      for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
        const auto lane = frame->active_lanes[child_pos];
        if (!compiled_math_batch_outcome_gates_open(
                *program,
                schedule.outcome_gate_nodes,
                *child_state,
                lane) ||
            !compiled_math_batch_time_gates_open(
                *program,
                schedule.time_gate_nodes,
                *child_state,
                lane)) {
          continue;
        }
        frame->term_lanes_a[count++] = lane;
        frame->products[static_cast<std::size_t>(lane)] = schedule.sign;
      }
      break;
    }
    state->current_lanes = frame->term_lanes_a.data();
    state->next_lanes = frame->term_lanes_b.data();
    state->term_count = count;
    return true;
  }

  bool begin_schedule(const CompiledMathLaneTermSchedule &) {
    return true;
  }

  bool apply_expr_upper_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || child_state == nullptr || frame == nullptr ||
        state == nullptr || state->current_lanes == nullptr ||
        state->next_lanes == nullptr) {
      return false;
    }
    std::size_t next_count = 0U;
    if (!compiled_math_batch_apply_expr_upper_op_to_lanes(
            *program,
            op,
            *child_state,
            state->current_lanes,
            state->term_count,
            frame->upper_found_by_lane.data(),
            frame->upper_time_by_lane.data(),
            frame->upper_normalizer_by_lane.data(),
            frame->products.data(),
            state->next_lanes,
            &next_count)) {
      return false;
    }
    std::swap(state->current_lanes, state->next_lanes);
    state->term_count = next_count;
    return true;
  }

  bool apply_constant_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    return compiled_math_apply_lane_term_constant_factor(op, state, frame);
  }

  bool apply_integral_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || kernel == nullptr || child_state == nullptr ||
        frame == nullptr || batch_workspace == nullptr) {
      return false;
    }
    return compiled_math_apply_lane_term_cached_integral_factor_op(
        *program,
        *kernel,
        op,
        *child_state,
        state,
        frame,
        batch_workspace,
        static_cast<std::size_t>(kernel->source_product_block_factors.size),
        child_capacity,
        depth);
  }

  bool apply_source_product_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    if (program == nullptr || kernel == nullptr || child_state == nullptr ||
        frame == nullptr) {
      return false;
    }
    return compiled_math_apply_lane_term_cached_factor_values(
        op,
        state,
        frame,
        static_cast<std::size_t>(kernel->source_product_block_factors.size),
        child_capacity,
        [&](const semantic::Index *lanes,
            const std::size_t lane_count,
            double *values_out) {
          return compiled_math_lane_term_source_product_factor_values(
              *program,
              op,
              *child_state,
              lanes,
              lane_count,
              frame,
              values_out);
        });
  }

  bool apply_direct_leaf_factor_op(
      const CompiledMathLaneTermOp &op,
      CompiledMathLaneTermRunState *state) {
    return apply_source_product_factor_op(op, state);
  }

  bool accumulate_term(
      const semantic::Index *current_lanes,
      const std::size_t term_count) {
    if (frame == nullptr || current_lanes == nullptr) {
      return false;
    }
    for (std::size_t i = 0; i < term_count; ++i) {
      const auto lane_pos = static_cast<std::size_t>(current_lanes[i]);
      frame->values[lane_pos] += frame->products[lane_pos];
    }
    has_accumulated_terms = true;
    return true;
  }

  bool finish_terms() {
    if (kernel == nullptr || frame == nullptr) {
      return false;
    }
    if (!kernel->clean_signed_source_sum) {
      return true;
    }
    if (!has_accumulated_terms) {
      return true;
    }
    for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
      const auto lane = frame->active_lanes[child_pos];
      const auto lane_index = static_cast<std::size_t>(lane);
      frame->values[lane_index] = clean_signed_value(frame->values[lane_index]);
    }
    return true;
  }
};

inline bool compiled_math_evaluate_child_lane_source_product_sum_block(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &child_state,
    const std::size_t child_count,
    const std::size_t child_capacity,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth) {
  if (frame == nullptr || batch_workspace == nullptr ||
      kernel.source_product_block_terms.size !=
          kernel.source_product_terms.size) {
    return false;
  }
  const auto block_factor_count =
      static_cast<std::size_t>(kernel.source_product_block_factors.size);
  if (kernel.source_product_block_reuses_factors) {
    compiled_math_reset_block_factor_cache(
        block_factor_count,
        child_capacity,
        frame);
  }
  for (std::size_t i = 0; i < child_count; ++i) {
    frame->values[static_cast<std::size_t>(frame->active_lanes[i])] = 0.0;
  }
  CompiledMathChildLaneTermDomain domain{
      &program,
      &kernel,
      &child_state,
      frame,
      batch_workspace,
      child_count,
      child_capacity,
      depth};
  return compiled_math_execute_lane_term_schedules(
      program,
      kernel.lane_term_schedules,
      &domain);
}

inline bool evaluate_compiled_integral_kernel_lane_batch_direct_source_product_sum(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane) {
  if (batch_workspace == nullptr || out_by_parent_lane == nullptr ||
      (bounds.count != 0U &&
       (bounds.lanes == nullptr || bounds.lower == nullptr ||
        bounds.upper == nullptr)) ||
      kernel.eval_kind !=
          CompiledMathIntegralKernelEvalKind::SourceProductBlock) {
    return false;
  }
  const auto block_factor_count =
      static_cast<std::size_t>(kernel.source_product_block_factors.size);
  auto &frame = compiled_math_batch_integral_frame(batch_workspace, depth);
  const auto child_capacity =
      bounds.count * kFiniteDirectLeafQuadratureTileOrder;
  const auto factor_parent_capacity = bounds.count * block_factor_count;
  const auto time_slot_count =
      std::max(
          parent_state.time_slot_count,
          static_cast<std::size_t>(kernel.bind_time_id) + 1U);
  frame.ensure_lane_capacity(
      time_slot_count,
      std::max(
          std::max(parent_state.lane_count, child_capacity),
          factor_parent_capacity));

  std::size_t eligible_count = 0U;
  for (std::size_t i = 0; i < bounds.count; ++i) {
    const auto lane = bounds.lanes[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    out_by_parent_lane[lane_pos] = 0.0;
    const double lower = bounds.lower[i];
    const double upper = bounds.upper[i];
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      continue;
    }
    frame.parent_lanes[eligible_count] = lane;
    frame.direct_parent_scale[eligible_count] = 0.5 * (upper - lower);
    frame.direct_parent_shift[eligible_count] = 0.5 * (upper + lower);
    frame.direct_source_channels_by_parent[eligible_count] =
        parent_state.source_channels(lane);
    ++eligible_count;
  }
  if (eligible_count == 0U) {
    return true;
  }

  std::size_t prepared_factor_count = block_factor_count;
  if (kernel.source_product_block_needs_direct_leaf_prepare) {
    if (!compiled_math_prepare_direct_tile_block_factors(
            program,
            kernel,
            parent_state,
            eligible_count,
            &frame,
            &prepared_factor_count) ||
        prepared_factor_count != block_factor_count) {
      return false;
    }
  }

  const bool needs_child_state = kernel.source_product_block_needs_child_state;
  const bool cache_reused_block_factors =
      kernel.source_product_block_reuses_factors;
  const bool cache_source_programs =
      kernel.source_product_block_needs_batch_source_cache;
  if (cache_source_programs) {
    compiled_math_collect_kernel_cached_source_vector_programs(
        program,
        kernel,
        &frame.cached_source_vector_programs);
  }
  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  for (std::size_t q_begin = 0; q_begin < quadrature::kDefaultFiniteOrder;
       q_begin += kFiniteDirectLeafQuadratureTileOrder) {
    const auto q_end =
        std::min(
            q_begin + kFiniteDirectLeafQuadratureTileOrder,
            static_cast<std::size_t>(quadrature::kDefaultFiniteOrder));
    const auto child_width = q_end - q_begin;
    const auto tile_capacity = eligible_count * child_width;
    CompiledMathLaneBatchState child_state{};
    const CompiledMathLaneBatchState *integral_child_state = nullptr;
    if (needs_child_state) {
      std::size_t child_count = 0U;
      for (std::size_t parent_pos = 0; parent_pos < eligible_count;
           ++parent_pos) {
        const auto parent_lane = frame.parent_lanes[parent_pos];
        for (std::size_t q_slot = 0; q_slot < child_width; ++q_slot) {
          const auto child_lane = static_cast<semantic::Index>(child_count);
          frame.active_lanes[child_count] = child_lane;
          frame.node_value_lane_map[child_count] =
              parent_state.node_value_lane_map == nullptr
                  ? parent_lane
                  : parent_state.node_value_lane_map[
                        static_cast<std::size_t>(parent_lane)];
          frame.workspaces_by_lane[child_count] =
              parent_state.workspace(parent_lane);
          frame.source_channels_by_lane[child_count] =
              parent_state.source_channels(parent_lane);
          frame.parents_by_lane[child_count] = parent_state.parent(parent_lane);
          frame.eval_workspaces_by_lane[child_count] =
              parent_state.eval_workspace(parent_lane);
          frame.used_outcomes_by_lane[child_count] =
              parent_state.used_outcomes(parent_lane);
          for (std::size_t time_slot = 0; time_slot < time_slot_count;
               ++time_slot) {
            double value = 0.0;
            if (time_slot < parent_state.time_slot_count) {
              value =
                  parent_state.time_slots.get(
                      static_cast<semantic::Index>(time_slot),
                      parent_lane);
            }
            frame.time_values[time_slot * tile_capacity + child_count] =
                value;
            frame.time_valid[time_slot * tile_capacity + child_count] =
                time_slot < parent_state.time_slot_count &&
                        parent_state.time_slots.has(
                            static_cast<semantic::Index>(time_slot),
                            parent_lane)
                    ? 1U
                    : 0U;
          }
          frame.time_values[
              static_cast<std::size_t>(kernel.bind_time_id) * tile_capacity +
              child_count] =
              frame.direct_parent_shift[parent_pos] +
              frame.direct_parent_scale[parent_pos] * rule.nodes[q_begin + q_slot];
          frame.time_valid[
              static_cast<std::size_t>(kernel.bind_time_id) * tile_capacity +
              child_count] = 1U;
          ++child_count;
        }
      }
      child_state = CompiledMathLaneBatchState{
          tile_capacity,
          time_slot_count,
          BatchTimeSlotView{
              frame.time_values.data(),
              tile_capacity,
              frame.time_valid.data()},
          &frame.workspaces_by_lane,
          &frame.source_channels_by_lane,
          &frame.parents_by_lane,
          &frame.eval_workspaces_by_lane,
          &frame.used_outcomes_by_lane,
          parent_state.node_values,
          parent_state.node_value_lane_stride,
          frame.node_value_lane_map.data()};
      (void)child_count;
      integral_child_state = &child_state;
    }
    if (cache_source_programs) {
      frame.reset_source_vector_cache(
          program.integral_kernel_source_vector_programs.size(),
          tile_capacity,
          frame.cached_source_vector_programs);
    }
    if (cache_reused_block_factors) {
      compiled_math_reset_block_factor_cache(
          block_factor_count,
          tile_capacity,
          &frame);
    }
    if (kernel.clean_signed_source_sum) {
      std::fill(
          frame.values.begin(),
          frame.values.begin() + static_cast<std::ptrdiff_t>(tile_capacity),
          0.0);
    }

    CompiledMathDirectTileBlockTermDomain domain{
        &program,
        &kernel,
        &parent_state,
        integral_child_state,
        &frame,
        batch_workspace,
        &rule,
        q_begin,
        child_width,
        eligible_count,
        tile_capacity,
        prepared_factor_count,
        depth,
        out_by_parent_lane,
        nullptr};
    if (!compiled_math_execute_lane_term_schedules(
            program,
            kernel.lane_term_schedules,
            &domain)) {
      return false;
    }
  }
  return true;
}

inline bool evaluate_compiled_integral_kernel_lane_batch_pdf_antiderivative(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane) {
  if (batch_workspace == nullptr || out_by_parent_lane == nullptr ||
      (bounds.count != 0U &&
       (bounds.lanes == nullptr || bounds.lower == nullptr ||
        bounds.upper == nullptr)) ||
      kernel.eval_kind !=
          CompiledMathIntegralKernelEvalKind::PdfAntiderivative) {
    return false;
  }

  const CompiledMathSourceProductOp *pdf_op = nullptr;
  const CompiledMathSourceProductChannel *pdf_channel = nullptr;
  double constant_scale = 1.0;
  for (semantic::Index op_idx = 0; op_idx < kernel.source_product_ops.size;
       ++op_idx) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(
                kernel.source_product_ops.offset + op_idx)];
    if (op.value_channel_mask == 0U) {
      constant_scale *= op.constant_value;
      continue;
    }
    if (pdf_op != nullptr || op.value_channel_mask != kLeafChannelPdf ||
        op.source_product_channel_id == semantic::kInvalidIndex) {
      return false;
    }
    const auto channel_pos =
        static_cast<std::size_t>(op.source_product_channel_id);
    if (channel_pos >= program.integral_kernel_source_product_channels.size()) {
      return false;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[channel_pos];
    if (!channel.direct_leaf_absolute_candidate ||
        channel.leaf_index == semantic::kInvalidIndex ||
        channel.time_id != kernel.bind_time_id ||
        (channel.time_cap_id != semantic::kInvalidIndex &&
         channel.time_cap_id != kernel.bind_time_id)) {
      return false;
    }
    pdf_op = &op;
    pdf_channel = &channel;
  }
  if (pdf_op == nullptr || pdf_channel == nullptr) {
    return false;
  }
  for (std::size_t i = 0; i < bounds.count; ++i) {
    out_by_parent_lane[static_cast<std::size_t>(bounds.lanes[i])] = 0.0;
  }
  if (bounds.count == 0U || constant_scale == 0.0) {
    return true;
  }

  auto &frame = compiled_math_batch_integral_frame(batch_workspace, depth);
  const auto lane_capacity =
      std::max(parent_state.lane_count, bounds.count);
  frame.ensure_lane_capacity(parent_state.time_slot_count, lane_capacity);
  std::size_t eval_count = 0U;
  for (std::size_t active_pos = 0; active_pos < bounds.count; ++active_pos) {
    const auto lane = bounds.lanes[active_pos];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double lower = bounds.lower[active_pos];
    const double upper = bounds.upper[active_pos];
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      continue;
    }
    frame.source_lanes_a[eval_count++] = lane;
    frame.lower_by_lane[lane_pos] = lower;
    frame.upper_by_lane[lane_pos] = upper;
  }
  if (eval_count == 0U) {
    return true;
  }

  auto &scratch =
      frame.source_product_scratch_layer(0U, parent_state.lane_count, eval_count);
  if (!compiled_math_batch_leaf_values_from_times(
          pdf_channel->leaf_dist_kind,
          pdf_channel->leaf_index,
          pdf_channel->leaf_onset_abs_value,
          parent_state,
          frame.source_lanes_a.data(),
          eval_count,
          frame.upper_by_lane.data(),
          kLeafChannelCdf,
          frame.source_leaf_values.data(),
          &scratch)) {
    return false;
  }
  if (!compiled_math_batch_leaf_values_from_times(
          pdf_channel->leaf_dist_kind,
          pdf_channel->leaf_index,
          pdf_channel->leaf_onset_abs_value,
          parent_state,
          frame.source_lanes_a.data(),
          eval_count,
          frame.lower_by_lane.data(),
          kLeafChannelCdf,
          frame.source_values.data(),
          &scratch)) {
    return false;
  }
  for (std::size_t i = 0; i < eval_count; ++i) {
    const auto lane = frame.source_lanes_a[i];
    const auto lane_pos = static_cast<std::size_t>(lane);
    const double value =
        constant_scale *
        (frame.source_leaf_values[lane_pos] - frame.source_values[lane_pos]);
    out_by_parent_lane[lane_pos] = std::isfinite(value) ? value : 0.0;
  }
  return true;
}

inline bool evaluate_compiled_integral_kernel_lane_batch_direct_leaf(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane) {
  if (batch_workspace == nullptr || out_by_parent_lane == nullptr ||
      (bounds.count != 0U &&
       (bounds.lanes == nullptr || bounds.lower == nullptr ||
        bounds.upper == nullptr)) ||
      kernel.eval_kind != CompiledMathIntegralKernelEvalKind::DirectLeaf) {
    return false;
  }
  const auto op_count =
      static_cast<std::size_t>(kernel.source_product_ops.size);
  auto &frame = compiled_math_batch_integral_frame(batch_workspace, depth);
  const auto time_slot_count =
      std::max(
          parent_state.time_slot_count,
          static_cast<std::size_t>(kernel.bind_time_id) + 1U);
  const auto max_child_capacity =
      bounds.count * kFiniteDirectLeafQuadratureTileOrder;
  const auto max_op_parent_capacity = bounds.count * op_count;
  frame.ensure_lane_capacity(
      time_slot_count,
      std::max(
          parent_state.lane_count,
          std::max(max_child_capacity, max_op_parent_capacity)));
  std::size_t eligible_count = 0U;
  for (std::size_t active_pos = 0; active_pos < bounds.count; ++active_pos) {
    const auto lane = bounds.lanes[active_pos];
    const auto lane_pos = static_cast<std::size_t>(lane);
    out_by_parent_lane[lane_pos] = 0.0;
    const double lower = bounds.lower[active_pos];
    const double upper = bounds.upper[active_pos];
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      continue;
    }
    frame.direct_parent_scale[eligible_count] = 0.5 * (upper - lower);
    frame.direct_parent_shift[eligible_count] = 0.5 * (upper + lower);
    frame.direct_source_channels_by_parent[eligible_count] =
        parent_state.source_channels(lane);
    frame.parent_lanes[eligible_count++] = lane;
  }
  if (eligible_count == 0U) {
    return true;
  }
  const auto op_parent_count = op_count * eligible_count;
  if (!compiled_math_ensure_direct_leaf_parent_capacity(
          op_count,
          eligible_count,
          &frame)) {
    return false;
  }

  for (std::size_t op_pos = 0; op_pos < op_count; ++op_pos) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(
                kernel.source_product_ops.offset +
                static_cast<semantic::Index>(op_pos))];
    auto &prepared = frame.direct_leaf_ops[op_pos];
    prepared.op = &op;
    prepared.channel = nullptr;
    prepared.tile_time_plan_id = semantic::kInvalidIndex;
    prepared.time_is_bind = false;
    prepared.has_time_cap = false;
    prepared.cap_is_bind = false;
    prepared.conditioned_direct_leaf = false;
    if (op.value_channel_mask == 0U) {
      continue;
    }
    if (op.source_product_channel_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(op.source_product_channel_id) >=
            program.integral_kernel_source_product_channels.size() ||
        op.direct_leaf_plan_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(op.direct_leaf_plan_id) >=
            program.integral_kernel_direct_leaf_op_plans.size()) {
      return false;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];
    const auto &direct_leaf_plan =
        program.integral_kernel_direct_leaf_op_plans[
            static_cast<std::size_t>(op.direct_leaf_plan_id)];
    prepared.channel = &channel;
    prepared.conditioned_direct_leaf =
        direct_leaf_plan.direct_leaf_conditioned;
    if (direct_leaf_plan.direct_leaf_index == semantic::kInvalidIndex) {
      return false;
    }
    prepared.time_is_bind = direct_leaf_plan.direct_leaf_time_is_bind;
    prepared.has_time_cap = direct_leaf_plan.direct_leaf_has_time_cap;
    prepared.cap_is_bind = direct_leaf_plan.direct_leaf_cap_is_bind;
    prepared.tile_time_plan_id =
        direct_leaf_plan.direct_leaf_tile_time_plan_id;
    if (!compiled_math_prepare_direct_leaf_parent_values(
            program,
            direct_leaf_plan,
            parent_state,
            frame.parent_lanes.data(),
            eligible_count,
            op_pos * eligible_count,
            op_parent_count,
            &frame)) {
      return false;
    }
  }

  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  for (std::size_t q_begin = 0; q_begin < quadrature::kDefaultFiniteOrder;
       q_begin += kFiniteDirectLeafQuadratureTileOrder) {
    const auto q_end =
        std::min(
            q_begin + kFiniteDirectLeafQuadratureTileOrder,
            static_cast<std::size_t>(quadrature::kDefaultFiniteOrder));
    const auto child_width = q_end - q_begin;
	    const auto tile_capacity = eligible_count * child_width;
	    std::fill(
	        frame.products.begin(),
	        frame.products.begin() + static_cast<std::ptrdiff_t>(tile_capacity),
	        1.0);
	    const std::size_t term_count = tile_capacity;
	    std::size_t live_count = term_count;
	    std::size_t active_count = term_count;
	    semantic::Index *active_lanes = nullptr;
	    semantic::Index *next_lanes = frame.source_lanes_a.data();
	    auto &direct_scratch =
	        frame.source_product_scratch_layer(
	            0U,
            tile_capacity,
            tile_capacity);

    for (std::size_t op_pos = 0; op_pos < op_count; ++op_pos) {
      if (live_count == 0U) {
        break;
      }
      const auto &prepared = frame.direct_leaf_ops[op_pos];
	      const auto *op = prepared.op;
	      if (op == nullptr) {
	        return false;
	      }
	      if (op->value_channel_mask == 0U) {
	        if (active_lanes == nullptr) {
	          active_count =
	              compiled_math_batch_array_multiply_scalar_dense_compact(
	                  term_count,
	                  op->constant_value,
	                  frame.products.data(),
	                  next_lanes);
	          if (active_count != term_count) {
	            active_lanes = next_lanes;
	            next_lanes =
	                active_lanes == frame.source_lanes_a.data()
	                    ? frame.source_lanes_b.data()
	                    : frame.source_lanes_a.data();
	          }
	        } else {
	          active_count =
	              compiled_math_batch_array_multiply_scalar_compact(
	                  active_lanes,
	                  active_count,
	                  op->constant_value,
	                  frame.products.data(),
	                  next_lanes);
	          std::swap(active_lanes, next_lanes);
	        }
	        live_count = active_count;
	        continue;
	      }
	      const auto *channel = prepared.channel;
	      if (channel == nullptr) {
	        return false;
	      }
	      const auto op_parent_offset = op_pos * eligible_count;
	      if (active_count == 0U) {
	        live_count = 0U;
	        break;
	      }
	      std::size_t factor_live_count = 0U;
	      auto *factor_values = frame.factor_values.data();
	      if (!compiled_math_direct_tile_leaf_update(
	              program,
	              prepared,
              op_parent_offset,
              eligible_count,
              q_begin,
	              rule,
	              child_width,
	              tile_capacity,
	              active_lanes,
	              active_count,
	              &frame,
	              &direct_scratch,
	              &factor_live_count,
	              factor_values)) {
	        return false;
	      }
	      (void)factor_live_count;
	      if (active_lanes == nullptr) {
	        active_count =
	            compiled_math_batch_array_multiply_values_dense_compact(
	                term_count,
	                factor_values,
	                frame.products.data(),
	                next_lanes);
	        if (active_count != term_count) {
	          active_lanes = next_lanes;
	          next_lanes =
	              active_lanes == frame.source_lanes_a.data()
	                  ? frame.source_lanes_b.data()
	                  : frame.source_lanes_a.data();
	        }
	      } else {
	        active_count =
	            compiled_math_batch_array_multiply_values_compact(
	                active_lanes,
	                active_count,
	                factor_values,
	                frame.products.data(),
	                next_lanes);
	        std::swap(active_lanes, next_lanes);
	      }
	      live_count = active_count;
	    }

    compiled_math_direct_tile_weighted_accumulate_parent_rows(
        frame.parent_lanes.data(),
        eligible_count,
        q_end - q_begin,
        q_begin,
        rule,
        frame.direct_parent_scale.data(),
        frame.products.data(),
        out_by_parent_lane);
  }
  return true;
}

struct CompiledMathQuadratureChildBatch {
  CompiledMathLaneBatchState state{};
  std::size_t child_count{0U};
  std::size_t child_capacity{0U};
};

inline bool compiled_math_build_quadrature_child_batch(
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathQuadratureChildBatch *child_batch) {
  if (frame == nullptr || child_batch == nullptr ||
      (bounds.count != 0U &&
       (bounds.lanes == nullptr || bounds.lower == nullptr ||
        bounds.upper == nullptr))) {
    return false;
  }
  const auto time_slot_count =
      std::max(
          parent_state.time_slot_count,
          static_cast<std::size_t>(kernel.bind_time_id) + 1U);
  const auto child_capacity =
      bounds.count * quadrature::kDefaultFiniteOrder;
  frame->ensure_lane_capacity(time_slot_count, child_capacity);

  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  std::size_t child_count = 0U;
  for (std::size_t parent_pos = 0; parent_pos < bounds.count; ++parent_pos) {
    const auto parent_lane = bounds.lanes[parent_pos];
    const double lower = bounds.lower[parent_pos];
    const double upper = bounds.upper[parent_pos];
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      continue;
    }
    const double scale = 0.5 * (upper - lower);
    const double shift = 0.5 * (upper + lower);
    for (std::size_t q = 0; q < quadrature::kDefaultFiniteOrder; ++q) {
      const auto child_lane = static_cast<semantic::Index>(child_count);
      frame->active_lanes[child_count] = child_lane;
      frame->parent_lanes[child_count] = parent_lane;
      frame->node_value_lane_map[child_count] =
          parent_state.node_value_lane_map == nullptr
              ? parent_lane
              : parent_state.node_value_lane_map[
                    static_cast<std::size_t>(parent_lane)];
      frame->weights[child_count] = scale * rule.weights[q];
      frame->workspaces_by_lane[child_count] =
          parent_state.workspace(parent_lane);
      frame->source_channels_by_lane[child_count] =
          parent_state.source_channels(parent_lane);
      frame->parents_by_lane[child_count] = parent_state.parent(parent_lane);
      frame->eval_workspaces_by_lane[child_count] =
          parent_state.eval_workspace(parent_lane);
      frame->used_outcomes_by_lane[child_count] =
          parent_state.used_outcomes(parent_lane);
      for (std::size_t time_slot = 0; time_slot < time_slot_count;
           ++time_slot) {
        double value = 0.0;
        if (time_slot < parent_state.time_slot_count) {
          value =
              parent_state.time_slots.get(
                  static_cast<semantic::Index>(time_slot),
                  parent_lane);
        }
        frame->time_values[time_slot * child_capacity + child_count] = value;
        frame->time_valid[time_slot * child_capacity + child_count] =
            time_slot < parent_state.time_slot_count &&
                    parent_state.time_slots.has(
                        static_cast<semantic::Index>(time_slot),
                        parent_lane)
                ? 1U
                : 0U;
      }
      frame->time_values[
          static_cast<std::size_t>(kernel.bind_time_id) * child_capacity +
          child_count] = shift + scale * rule.nodes[q];
      frame->time_valid[
          static_cast<std::size_t>(kernel.bind_time_id) * child_capacity +
          child_count] = 1U;
      ++child_count;
    }
  }

  child_batch->child_count = child_count;
  child_batch->child_capacity = child_capacity;
  child_batch->state = CompiledMathLaneBatchState{
      child_capacity,
      time_slot_count,
      BatchTimeSlotView{
          frame->time_values.data(),
          child_capacity,
          frame->time_valid.data()},
      &frame->workspaces_by_lane,
      &frame->source_channels_by_lane,
      &frame->parents_by_lane,
      &frame->eval_workspaces_by_lane,
      &frame->used_outcomes_by_lane,
      parent_state.node_values,
      parent_state.node_value_lane_stride,
      frame->node_value_lane_map.data()};
  return true;
}

inline bool compiled_math_evaluate_child_lane_kernel_values(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathQuadratureChildBatch &child_batch,
    CompiledMathBatchIntegralFrame *frame,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth) {
  if (frame == nullptr || batch_workspace == nullptr) {
    return false;
  }
  const BatchActiveLaneSpan child_active{
      frame->active_lanes.data(),
      child_batch.child_count};
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
    return compiled_math_batch_source_product_value_for_ops(
        program,
        kernel.source_product_ops,
        child_batch.state,
        child_active,
        frame,
        frame->values.data());
  }
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProductSum) {
    return compiled_math_evaluate_child_lane_source_product_sum_block(
        program,
        kernel,
        child_batch.state,
        child_batch.child_count,
        child_batch.child_capacity,
        frame,
        batch_workspace,
        depth);
  }
  return false;
}

inline bool evaluate_compiled_integral_kernel_lane_batch_quadrature_child_program(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane) {
  if (batch_workspace == nullptr || out_by_parent_lane == nullptr ||
      (bounds.count != 0U &&
       (bounds.lanes == nullptr || bounds.lower == nullptr ||
        bounds.upper == nullptr)) ||
      (kernel.eval_kind != CompiledMathIntegralKernelEvalKind::FlatQuadrature &&
       kernel.eval_kind !=
           CompiledMathIntegralKernelEvalKind::GenericQuadrature)) {
    return false;
  }

  auto &frame = compiled_math_batch_integral_frame(batch_workspace, depth);
  CompiledMathQuadratureChildBatch child_batch{};
  if (!compiled_math_build_quadrature_child_batch(
          kernel,
          parent_state,
          bounds,
          &frame,
          &child_batch)) {
    return false;
  }
  if (child_batch.child_count == 0U) {
    return true;
  }

  if (kernel.source_product_block_needs_batch_source_cache) {
    compiled_math_collect_kernel_cached_source_vector_programs(
        program,
        kernel,
        &frame.cached_source_vector_programs);
    frame.reset_source_vector_cache(
        program.integral_kernel_source_vector_programs.size(),
        child_batch.child_capacity,
        frame.cached_source_vector_programs);
  }
  if (!compiled_math_evaluate_child_lane_kernel_values(
          program,
          kernel,
          child_batch,
          &frame,
          batch_workspace,
          depth)) {
    return false;
  }

  compiled_math_batch_array_weighted_accumulate_parent(
      frame.active_lanes.data(),
      frame.parent_lanes.data(),
      child_batch.child_count,
      frame.weights.data(),
      frame.values.data(),
      out_by_parent_lane);
  return true;
}

inline bool evaluate_compiled_integral_kernel_lane_batch(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathIntegralKernel &kernel,
    const CompiledMathLaneBatchState &parent_state,
    const CompiledMathActiveBoundBatch &bounds,
    CompiledMathBatchWorkspace *batch_workspace,
    const std::size_t depth,
    double *out_by_parent_lane) {
  if (batch_workspace == nullptr || out_by_parent_lane == nullptr ||
      (bounds.count != 0U &&
       (bounds.lanes == nullptr || bounds.lower == nullptr ||
        bounds.upper == nullptr))) {
    return false;
  }
  (void)node;
  for (std::size_t i = 0; i < bounds.count; ++i) {
    out_by_parent_lane[static_cast<std::size_t>(bounds.lanes[i])] =
        0.0;
  }
  if (bounds.count == 0U) {
    return true;
  }
  switch (kernel.eval_kind) {
  case CompiledMathIntegralKernelEvalKind::PdfAntiderivative:
    return evaluate_compiled_integral_kernel_lane_batch_pdf_antiderivative(
        program,
        kernel,
        parent_state,
        bounds,
        batch_workspace,
        depth,
        out_by_parent_lane);
  case CompiledMathIntegralKernelEvalKind::DirectLeaf:
    return evaluate_compiled_integral_kernel_lane_batch_direct_leaf(
        program,
        kernel,
        parent_state,
        bounds,
        batch_workspace,
        depth,
        out_by_parent_lane);
  case CompiledMathIntegralKernelEvalKind::SourceProductBlock:
    return evaluate_compiled_integral_kernel_lane_batch_direct_source_product_sum(
        program,
        kernel,
        parent_state,
        bounds,
        batch_workspace,
        depth,
        out_by_parent_lane);
  case CompiledMathIntegralKernelEvalKind::FlatQuadrature:
  case CompiledMathIntegralKernelEvalKind::GenericQuadrature:
    return evaluate_compiled_integral_kernel_lane_batch_quadrature_child_program(
        program,
        kernel,
        parent_state,
        bounds,
        batch_workspace,
        depth,
        out_by_parent_lane);
  }
  return false;
}

inline bool evaluate_compiled_integral_node_batch_lane_native(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathIntegralKernel &kernel,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *node_values) {
  if (batch_workspace == nullptr || node_values == nullptr) {
    return false;
  }

  const std::size_t time_slot_count =
      compiled_math_program_time_slot_count(program);
  auto &handled_lane_indices = batch_workspace->active_lane_indices;
  auto &compact_lane_indices = batch_workspace->compact_lane_indices;
  auto &active_lanes = batch_workspace->active_lanes;
  auto &lower_by_lane = batch_workspace->lower_by_lane;
  auto &upper_by_lane = batch_workspace->upper_by_lane;
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &source_channels_by_lane = batch_workspace->source_channels_by_lane;
  auto &parents_by_lane = batch_workspace->parents_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;
  auto &node_value_lane_map = batch_workspace->node_value_lane_map;
  auto &batch_values = batch_workspace->batch_values;

  handled_lane_indices.clear();
  compact_lane_indices.clear();
  active_lanes.clear();
  lower_by_lane.clear();
  upper_by_lane.clear();

  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    auto *workspace = lanes[lane].workspace;
    auto *parent = lanes[lane].parent;
    if (workspace == nullptr || lane >= condition_evaluators.size() ||
        condition_evaluators[lane] == nullptr || parent == nullptr ||
        parent->source_channels() == nullptr) {
      handled_lane_indices.push_back(lane);
      (*node_values)[lane] = 0.0;
      continue;
    }
    const auto compact_lane =
        static_cast<semantic::Index>(compact_lane_indices.size());
    handled_lane_indices.push_back(lane);
    compact_lane_indices.push_back(lane);
    const double upper = compiled_math_node_time(node, *workspace);
    if (std::isfinite(upper) && upper > 0.0) {
      active_lanes.push_back(compact_lane);
      lower_by_lane.push_back(0.0);
      upper_by_lane.push_back(upper);
    }
  }

  if (handled_lane_indices.empty()) {
    return false;
  }
  if (compact_lane_indices.empty()) {
    return true;
  }
  const std::size_t compact_count = compact_lane_indices.size();
  parent_time_values.assign(time_slot_count * compact_count, 0.0);
  parent_time_valid.assign(time_slot_count * compact_count, 0U);
  workspaces_by_lane.assign(compact_count, nullptr);
  source_channels_by_lane.assign(compact_count, nullptr);
  parents_by_lane.assign(compact_count, nullptr);
  eval_workspaces_by_lane.assign(compact_count, nullptr);
  used_outcomes_by_lane.assign(compact_count, nullptr);
  node_value_lane_map.assign(compact_count, semantic::kInvalidIndex);
  batch_values.assign(compact_count, 0.0);

  for (std::size_t compact_lane = 0; compact_lane < compact_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    auto *workspace = lanes[lane].workspace;
    for (std::size_t time_slot = 0; time_slot < time_slot_count;
         ++time_slot) {
      const auto time_pos = time_slot * compact_count + compact_lane;
      if (time_slot < workspace->time_values.size()) {
        parent_time_values[time_pos] = workspace->time_values[time_slot];
      }
      parent_time_valid[time_pos] =
          workspace->has_time(static_cast<semantic::Index>(time_slot))
              ? 1U
              : 0U;
    }
    workspaces_by_lane[compact_lane] = workspace;
    source_channels_by_lane[compact_lane] =
        lanes[lane].parent->source_channels();
    parents_by_lane[compact_lane] = lanes[lane].parent;
    eval_workspaces_by_lane[compact_lane] = lanes[lane].eval_workspace;
    used_outcomes_by_lane[compact_lane] = workspace->used_outcomes;
    node_value_lane_map[compact_lane] = static_cast<semantic::Index>(lane);
  }

  const CompiledMathLaneBatchState parent_state{
      compact_count,
      time_slot_count,
      BatchTimeSlotView{
          parent_time_values.data(),
          compact_count,
          parent_time_valid.data()},
      &workspaces_by_lane,
      &source_channels_by_lane,
      &parents_by_lane,
      &eval_workspaces_by_lane,
      &used_outcomes_by_lane,
      batch_workspace->node_values.data(),
      batch_workspace->node_value_lane_stride,
      node_value_lane_map.data()};
  const CompiledMathActiveBoundBatch active_bounds{
      active_lanes.data(),
      active_lanes.size(),
      lower_by_lane.data(),
      upper_by_lane.data()};
  if (!evaluate_compiled_integral_kernel_lane_batch(
          program,
          node,
          kernel,
          parent_state,
          active_bounds,
          batch_workspace,
          0U,
          batch_values.data())) {
    return false;
  }

  for (std::size_t compact_lane = 0; compact_lane < compact_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    double value = batch_values[compact_lane];
    if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent) {
      value = clamp_probability(value);
    } else if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw) {
      value = std::isfinite(value) ? clean_signed_value(value) : 0.0;
    } else {
      value = std::isfinite(value) ? value : 0.0;
    }
    (*node_values)[lane] = value;
  }
  return true;
}

inline void evaluate_compiled_integral_node_batch(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *node_values) {
  node_values->assign(lanes.size(), 0.0);
  if (node.integral_kernel_slot == semantic::kInvalidIndex) {
    return;
  }
  const auto kernel_pos =
      static_cast<std::size_t>(node.integral_kernel_slot);
  if (kernel_pos >= program.integral_kernels.size()) {
    throw std::runtime_error(
        "compiled integral node points outside integral kernels");
  }
  const auto &kernel = program.integral_kernels[kernel_pos];
  std::vector<std::uint8_t> batched(lanes.size(), 0U);
  if (evaluate_compiled_integral_node_batch_lane_native(
          program,
          node,
          kernel,
          lanes,
          condition_evaluators,
          batch_workspace,
          node_values)) {
    for (const auto lane : batch_workspace->active_lane_indices) {
      batched[lane] = 1U;
    }
  }

  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    if (batched[lane] != 0U) {
      continue;
    }
    auto *workspace = lanes[lane].workspace;
    if (workspace == nullptr || lane >= condition_evaluators.size() ||
        condition_evaluators[lane] == nullptr ||
        node.integral_kernel_slot == semantic::kInvalidIndex) {
      (*node_values)[lane] = 0.0;
      continue;
    }
    const double current_time = compiled_math_node_time(node, *workspace);
    if (!(current_time > 0.0)) {
      (*node_values)[lane] = 0.0;
      continue;
    }
    throw std::runtime_error(
        "compiled batch integral node has no batch implementation");
  }
}

inline bool prepare_compiled_source_node_batch_base(
    const CompiledMathProgram &program,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  if (batch_workspace == nullptr) {
    return false;
  }
  auto &compact_lane_indices = batch_workspace->compact_lane_indices;
  auto &active_lanes = batch_workspace->active_lanes;
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;
  auto &node_value_lane_map = batch_workspace->node_value_lane_map;

  compact_lane_indices.clear();
  active_lanes.clear();
  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    if (lanes[lane].workspace == nullptr) {
      return false;
    }
    compact_lane_indices.push_back(lane);
    active_lanes.push_back(static_cast<semantic::Index>(lane));
  }
  const std::size_t lane_count = compact_lane_indices.size();
  if (lane_count == 0U) {
    return true;
  }

  parent_time_values.resize(time_slot_count * lane_count);
  parent_time_valid.resize(time_slot_count * lane_count);
  workspaces_by_lane.resize(lane_count);
  eval_workspaces_by_lane.resize(lane_count);
  used_outcomes_by_lane.resize(lane_count);
  node_value_lane_map.resize(lane_count);

  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    auto *workspace = lanes[lane].workspace;
    for (std::size_t time_slot = 0; time_slot < time_slot_count;
         ++time_slot) {
      const auto time_pos = time_slot * lane_count + compact_lane;
      if (time_slot < workspace->time_values.size()) {
        parent_time_values[time_pos] = workspace->time_values[time_slot];
      }
      parent_time_valid[time_pos] =
          workspace->has_time(static_cast<semantic::Index>(time_slot))
              ? 1U
              : 0U;
    }
    workspaces_by_lane[compact_lane] = workspace;
    eval_workspaces_by_lane[compact_lane] = lanes[lane].eval_workspace;
    used_outcomes_by_lane[compact_lane] = workspace->used_outcomes;
    node_value_lane_map[compact_lane] = static_cast<semantic::Index>(lane);
  }
  (void)program;
  return true;
}

inline bool evaluate_compiled_source_node_batch_prepared(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  if (batch_workspace == nullptr ||
      node.source_vector_program_id == semantic::kInvalidIndex) {
    return false;
  }

  const auto value_mask = compiled_math_source_factor_channel_mask(node.kind);
  const auto fill_mask = compiled_math_source_node_fill_mask(node.kind);
  if (value_mask == 0U || fill_mask == 0U) {
    return false;
  }

  auto &compact_lane_indices = batch_workspace->compact_lane_indices;
  auto &active_lanes = batch_workspace->active_lanes;
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &source_channels_by_lane = batch_workspace->source_channels_by_lane;
  auto &parents_by_lane = batch_workspace->parents_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;
  auto &node_value_lane_map = batch_workspace->node_value_lane_map;
  const std::size_t lane_count = compact_lane_indices.size();
  if (lane_count == 0U) {
    return true;
  }

  source_channels_by_lane.resize(lane_count);
  parents_by_lane.resize(lane_count);
  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    if (lane >= condition_evaluators.size() ||
        condition_evaluators[lane] == nullptr ||
        condition_evaluators[lane]->source_channels() == nullptr) {
      return false;
    }
    source_channels_by_lane[compact_lane] =
        condition_evaluators[lane]->source_channels();
    parents_by_lane[compact_lane] = condition_evaluators[lane];
  }

  auto &frame = compiled_math_batch_integral_frame(batch_workspace, 0U);
  frame.ensure_lane_capacity(time_slot_count, lane_count);
  auto &scratch =
      frame.source_product_scratch_layer(0U, lane_count, lane_count);
  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    scratch.times[compact_lane] =
        compiled_math_source_node_time(node, *lanes[lane].workspace);
    scratch.mask_a[compact_lane] = 0U;
    scratch.pdf_a[compact_lane] = 0.0;
    scratch.cdf_a[compact_lane] = 0.0;
    scratch.survival_a[compact_lane] = 1.0;
  }

  const CompiledMathLaneBatchState state{
      lane_count,
      time_slot_count,
      BatchTimeSlotView{
          parent_time_values.data(),
          lane_count,
          parent_time_valid.data()},
      &workspaces_by_lane,
      &source_channels_by_lane,
      &parents_by_lane,
      &eval_workspaces_by_lane,
      &used_outcomes_by_lane,
      batch_workspace->node_values.data(),
      batch_workspace->node_value_lane_stride,
      node_value_lane_map.data()};
  const auto source_vector_program_id = node.source_vector_program_id;
  if (source_vector_program_id == semantic::kInvalidIndex) {
    return false;
  }
  const bool supported =
      compiled_math_batch_source_vector_program_fill(
          program,
          source_vector_program_id,
          state,
          active_lanes.data(),
          lane_count,
          scratch.times.data(),
          fill_mask,
          &frame,
          1U,
          scratch.mask_a.data(),
          scratch.pdf_a.data(),
          scratch.cdf_a.data(),
          scratch.survival_a.data());
  if (!supported) {
    return false;
  }

  for (std::size_t compact_lane = 0; compact_lane < lane_count;
       ++compact_lane) {
    const auto lane = compact_lane_indices[compact_lane];
    const double raw =
        value_mask == kLeafChannelPdf
            ? safe_density(scratch.pdf_a[compact_lane])
            : (value_mask == kLeafChannelCdf
                   ? clamp_probability(scratch.cdf_a[compact_lane])
                   : clamp_probability(scratch.survival_a[compact_lane]));
    compiled_math_batch_set_node_value(
        node,
        static_cast<semantic::Index>(lane),
        raw,
        batch_workspace);
  }
  return true;
}

inline CompiledMathLaneBatchState prepare_compiled_regular_node_batch_state(
    const CompiledMathProgram &program,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  auto &parent_time_values = batch_workspace->parent_time_values;
  auto &parent_time_valid = batch_workspace->parent_time_valid;
  auto &workspaces_by_lane = batch_workspace->workspaces_by_lane;
  auto &parents_by_lane = batch_workspace->parents_by_lane;
  auto &eval_workspaces_by_lane = batch_workspace->eval_workspaces_by_lane;
  auto &used_outcomes_by_lane = batch_workspace->used_outcomes_by_lane;

  const std::size_t lane_count = lanes.size();
  parent_time_values.assign(time_slot_count * lane_count, 0.0);
  parent_time_valid.assign(time_slot_count * lane_count, 0U);
  workspaces_by_lane.assign(lane_count, nullptr);
  parents_by_lane.assign(lane_count, nullptr);
  eval_workspaces_by_lane.assign(lane_count, nullptr);
  used_outcomes_by_lane.assign(lane_count, nullptr);

  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    auto *workspace = lanes[lane].workspace;
    if (workspace != nullptr) {
      for (std::size_t time_slot = 0; time_slot < time_slot_count;
           ++time_slot) {
        const auto time_pos = time_slot * lane_count + lane;
        if (time_slot < workspace->time_values.size()) {
          parent_time_values[time_pos] = workspace->time_values[time_slot];
        }
        parent_time_valid[time_pos] =
            workspace->has_time(static_cast<semantic::Index>(time_slot))
                ? 1U
                : 0U;
      }
      used_outcomes_by_lane[lane] = workspace->used_outcomes;
    }
    workspaces_by_lane[lane] = workspace;
    parents_by_lane[lane] = lanes[lane].parent;
    eval_workspaces_by_lane[lane] = lanes[lane].eval_workspace;
  }

  return CompiledMathLaneBatchState{
      lane_count,
      time_slot_count,
      BatchTimeSlotView{
          parent_time_values.data(),
          lane_count,
          parent_time_valid.data()},
      &workspaces_by_lane,
      nullptr,
      &parents_by_lane,
      &eval_workspaces_by_lane,
      &used_outcomes_by_lane,
      batch_workspace->node_values.data(),
      batch_workspace->node_value_lane_stride,
      nullptr};
}

inline double compiled_math_batch_union_kernel_density_for_lane(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane,
    const ExactVariantPlan &plan,
    const bool multi_subset_only) {
  const auto expr_pos = static_cast<std::size_t>(node.subject_id);
  if (expr_pos >= plan.expr_kernels.size()) {
    return 0.0;
  }
  const auto &kernel = plan.expr_kernels[expr_pos];
  if (kernel.children.empty()) {
    return 0.0;
  }
  const auto child_count = kernel.children.size;
  auto child_value = [&](const semantic::Index child_index) {
    if (child_index >= node.children.size) {
      return 0.0;
    }
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + child_index)];
    return state.node_value(child_id, lane);
  };
  if (!multi_subset_only && child_count == 1U) {
    return child_value(0);
  }

  double total = 0.0;
  for (semantic::Index subset_idx = 0;
       subset_idx < kernel.union_subset_span.size;
       ++subset_idx) {
    const auto &subset = plan.expr_union_subsets[
        static_cast<std::size_t>(
            kernel.union_subset_span.offset + subset_idx)];
    if (multi_subset_only && subset.child_positions.size <= 1U) {
      continue;
    }
    if (subset.child_positions.empty()) {
      continue;
    }
    double subset_value = 0.0;
    if (subset.child_positions.size == 1U) {
      const auto active_pos = plan.expr_union_subset_child_positions[
          static_cast<std::size_t>(subset.child_positions.offset)];
      subset_value = child_value(active_pos);
    } else {
      for (semantic::Index i = 0; i < subset.child_positions.size; ++i) {
        const auto active_pos = plan.expr_union_subset_child_positions[
            static_cast<std::size_t>(subset.child_positions.offset + i)];
        double term = child_value(active_pos);
        if (!std::isfinite(term) || term == 0.0) {
          continue;
        }
        for (semantic::Index j = 0; j < subset.child_positions.size; ++j) {
          if (i == j) {
            continue;
          }
          const auto other_pos = plan.expr_union_subset_child_positions[
              static_cast<std::size_t>(subset.child_positions.offset + j)];
          term *= child_value(child_count + other_pos);
          if (!std::isfinite(term) || term == 0.0) {
            break;
          }
        }
        subset_value += term;
      }
      subset_value = clean_signed_value(subset_value);
    }
    total += static_cast<double>(subset.sign) * subset_value;
  }
  return clean_signed_value(total);
}

inline double compiled_math_batch_union_kernel_cdf_for_lane(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledMathLaneBatchState &state,
    const semantic::Index lane) {
  double total = 0.0;
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + i)];
    total += state.node_value(child_id, lane);
  }
  return clamp_probability(total);
}

inline void evaluate_compiled_regular_node_batch(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const std::vector<CompiledMathBatchLane> &lanes,
    const std::vector<CompiledSourceView *> &condition_evaluators,
    const std::size_t time_slot_count,
    CompiledMathBatchWorkspace *batch_workspace) {
  if (compiled_math_is_source_value_node(node.kind)) {
    throw std::runtime_error(
        "compiled batch source node reached regular-node evaluator");
  }
  const auto state =
      prepare_compiled_regular_node_batch_state(
          program,
          lanes,
          time_slot_count,
          batch_workspace);
  const auto write_value = [&](const std::size_t lane_pos,
                               const double value) {
    compiled_math_batch_set_node_value(
        node,
        static_cast<semantic::Index>(lane_pos),
        value,
        batch_workspace);
  };
  const auto first_child_id = [&]() {
    return program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
  };
  double *node_out =
      compiled_math_batch_node_value_array(node, batch_workspace);
  if (node_out == nullptr) {
    throw std::runtime_error("compiled batch node has no batch value storage");
  }
  auto &regular_active_lanes = batch_workspace->active_lanes;
  const auto regular_active_count =
      compiled_math_batch_active_workspace_lanes(
          lanes,
          &regular_active_lanes);
  const auto prepare_upper_lane_scratch = [&]() {
    auto &active_lanes = batch_workspace->active_lanes;
    active_lanes.clear();
    active_lanes.reserve(lanes.size());
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      if (lanes[lane_pos].workspace != nullptr) {
        active_lanes.push_back(static_cast<semantic::Index>(lane_pos));
      }
    }
    const auto lane_count = lanes.size();
    batch_workspace->upper_found_by_lane.assign(lane_count, 0U);
    batch_workspace->upper_time_by_lane.assign(
        lane_count,
        std::numeric_limits<double>::infinity());
    batch_workspace->upper_normalizer_by_lane.assign(lane_count, 0.0);
    return active_lanes.size();
  };

  switch (node.kind) {
  case CompiledMathNodeKind::Constant:
    compiled_math_batch_regular_constant(
        regular_active_lanes.data(),
        regular_active_count,
        node.constant,
        node_out);
    return;
  case CompiledMathNodeKind::SourcePdf:
  case CompiledMathNodeKind::SourceCdf:
  case CompiledMathNodeKind::SourceSurvival:
    throw std::runtime_error(
        "compiled batch source node reached regular-node evaluator");
  case CompiledMathNodeKind::ExprDensity:
  case CompiledMathNodeKind::ExprCdf:
  case CompiledMathNodeKind::ExprSurvival:
    throw std::runtime_error(
        "compiled math root contains an interpreter expression node");
  case CompiledMathNodeKind::UnionKernelMultiSubsetCdf:
  case CompiledMathNodeKind::IntegralZeroToCurrent:
  case CompiledMathNodeKind::IntegralZeroToCurrentRaw:
    throw std::runtime_error("integral node was not handled by batch path");
  case CompiledMathNodeKind::TimeGate: {
    const auto child_id = first_child_id();
    compiled_math_batch_regular_time_gate(
        node,
        child_id,
        state,
        regular_active_lanes.data(),
        regular_active_count,
        node_out);
    return;
  }
  case CompiledMathNodeKind::SimpleGuardDensity:
  case CompiledMathNodeKind::SimpleGuardCdf: {
    const auto child_id = first_child_id();
    const auto active_count = prepare_upper_lane_scratch();
    if (active_count != 0U &&
        !compiled_math_batch_upper_bounds_from_terms(
            program,
            node,
            state,
            batch_workspace->active_lanes.data(),
            active_count,
            batch_workspace->upper_found_by_lane.data(),
            batch_workspace->upper_time_by_lane.data(),
            batch_workspace->upper_normalizer_by_lane.data())) {
      throw std::runtime_error(
          "compiled batch simple guard upper bound is unsupported");
    }
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      auto *parent = lanes[lane_pos].parent;
      auto *condition_evaluator = condition_evaluators[lane_pos];
      if (workspace == nullptr) {
        continue;
      }
      if (condition_evaluator == nullptr || parent == nullptr) {
        write_value(lane_pos, 0.0);
        continue;
      }
      const auto &kernel = parent->plan().expr_kernels[
          static_cast<std::size_t>(node.subject_id)];
      if (!kernel.simple_event_guard || kernel.has_unless) {
        throw std::runtime_error(
            "compiled simple guard reached a non-simple guard");
      }
      const auto lane = static_cast<semantic::Index>(lane_pos);
      const double raw = state.node_value(child_id, lane);
      if (condition_evaluator->source_order_known_before(
              kernel.guard_blocker_source_id,
              kernel.guard_ref_source_id,
              node.condition_id,
              workspace)) {
        write_value(lane_pos, 0.0);
        continue;
      }
      const auto lane_index = static_cast<std::size_t>(lane);
      if (batch_workspace->upper_found_by_lane[lane_index] == 0U) {
        write_value(
            lane_pos,
            node.kind == CompiledMathNodeKind::SimpleGuardCdf
                ? clamp_probability(raw)
                : clean_signed_value(raw));
        continue;
      }
      const double current_time = state.time(node.time_id, lane);
      double value = 0.0;
      const double upper_time = batch_workspace->upper_time_by_lane[lane_index];
      const double normalizer =
          batch_workspace->upper_normalizer_by_lane[lane_index];
      if (node.kind == CompiledMathNodeKind::SimpleGuardCdf) {
        if (!(normalizer > 0.0)) {
          value = 0.0;
        } else if (current_time >= upper_time) {
          value = 1.0;
        } else {
          value = clamp_probability(raw / normalizer);
        }
      } else {
        if (!(normalizer > 0.0) || current_time >= upper_time) {
          value = 0.0;
        } else {
          value = clean_signed_value(raw / normalizer);
        }
      }
      write_value(lane_pos, value);
    }
    return;
  }
  case CompiledMathNodeKind::UnionKernelDensity:
  case CompiledMathNodeKind::UnionKernelMultiSubsetDensity: {
    const bool multi_subset_only =
        node.kind == CompiledMathNodeKind::UnionKernelMultiSubsetDensity;
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      auto *parent = lanes[lane_pos].parent;
      if (workspace == nullptr) {
        continue;
      }
      write_value(
          lane_pos,
          parent == nullptr
              ? 0.0
              : compiled_math_batch_union_kernel_density_for_lane(
                    program,
                    node,
                    state,
                    static_cast<semantic::Index>(lane_pos),
                    parent->plan(),
                    multi_subset_only));
    }
    return;
  }
  case CompiledMathNodeKind::UnionKernelCdf:
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace != nullptr) {
        write_value(
            lane_pos,
            compiled_math_batch_union_kernel_cdf_for_lane(
                program,
                node,
                state,
                static_cast<semantic::Index>(lane_pos)));
      }
    }
    return;
  case CompiledMathNodeKind::OutcomeSubsetUnused:
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      auto *parent = lanes[lane_pos].parent;
      if (workspace == nullptr) {
        continue;
      }
      double value = 1.0;
      if (workspace->used_outcomes != nullptr) {
        if (parent == nullptr) {
          throw std::runtime_error(
              "compiled outcome-used gate has no parent source view");
        }
        const auto offset = static_cast<std::size_t>(node.subject_id);
        const auto size = static_cast<std::size_t>(node.aux_id);
        const auto &indices = parent->plan().compiled_outcome_gate_indices;
        if (offset + size > indices.size()) {
          throw std::runtime_error(
              "compiled outcome-used gate points outside the plan");
        }
        for (std::size_t i = 0; i < size; ++i) {
          const auto outcome_idx = indices[offset + i];
          if (outcome_idx != semantic::kInvalidIndex &&
              static_cast<std::size_t>(outcome_idx) <
                  workspace->used_outcomes->size() &&
              (*workspace->used_outcomes)[
                  static_cast<std::size_t>(outcome_idx)] != 0U) {
            value = 0.0;
            break;
          }
        }
      }
      write_value(lane_pos, value);
    }
    return;
  case CompiledMathNodeKind::ExprUpperBoundDensity:
  case CompiledMathNodeKind::ExprUpperBoundCdf: {
    const auto child_id = first_child_id();
    const auto active_count = prepare_upper_lane_scratch();
    if (active_count != 0U &&
        !compiled_math_batch_expr_upper_bounds_for_node(
            program,
            node,
            state,
            batch_workspace->active_lanes.data(),
            active_count,
            condition_evaluators.data(),
            batch_workspace->upper_found_by_lane.data(),
            batch_workspace->upper_time_by_lane.data(),
            batch_workspace->upper_normalizer_by_lane.data())) {
      throw std::runtime_error(
          "compiled batch expr upper bound is unsupported");
    }
    for (std::size_t lane_pos = 0; lane_pos < lanes.size(); ++lane_pos) {
      auto *workspace = lanes[lane_pos].workspace;
      if (workspace == nullptr) {
        continue;
      }
      const auto lane = static_cast<semantic::Index>(lane_pos);
      const double raw = state.node_value(child_id, lane);
      const auto lane_index = static_cast<std::size_t>(lane);
      if (batch_workspace->upper_found_by_lane[lane_index] == 0U) {
        write_value(lane_pos, raw);
        continue;
      }
      const double current_time = state.time(node.time_id, lane);
      double value = 0.0;
      const double upper_time = batch_workspace->upper_time_by_lane[lane_index];
      const double normalizer =
          batch_workspace->upper_normalizer_by_lane[lane_index];
      if (node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
        if (!(normalizer > 0.0)) {
          value = 0.0;
        } else if (current_time >= upper_time) {
          value = 1.0;
        } else {
          value = clamp_probability(raw / normalizer);
        }
      } else {
        if (!(normalizer > 0.0) ||
            current_time >= upper_time) {
          value = 0.0;
        } else {
          value = safe_density(raw / normalizer);
        }
      }
      write_value(lane_pos, value);
    }
    return;
  }
  case CompiledMathNodeKind::Product:
    compiled_math_batch_regular_product(
        program,
        node,
        state,
        regular_active_lanes.data(),
        regular_active_count,
        node_out);
    return;
  case CompiledMathNodeKind::Sum:
  case CompiledMathNodeKind::CleanSignedSum:
    compiled_math_batch_regular_sum(
        program,
        node,
        state,
        regular_active_lanes.data(),
        regular_active_count,
        node.kind == CompiledMathNodeKind::CleanSignedSum,
        node_out);
    return;
  case CompiledMathNodeKind::ClampProbability: {
    const auto child_id = first_child_id();
    compiled_math_batch_regular_unary(
        node.kind,
        child_id,
        state,
        regular_active_lanes.data(),
        regular_active_count,
        node_out);
    return;
  }
  case CompiledMathNodeKind::Complement: {
    const auto child_id = first_child_id();
    compiled_math_batch_regular_unary(
        node.kind,
        child_id,
        state,
        regular_active_lanes.data(),
        regular_active_count,
        node_out);
    return;
  }
  case CompiledMathNodeKind::Negate: {
    const auto child_id = first_child_id();
    compiled_math_batch_regular_unary(
        node.kind,
        child_id,
        state,
        regular_active_lanes.data(),
        regular_active_count,
        node_out);
    return;
  }
  }
}

inline void evaluate_compiled_node_span_batch(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan schedule,
    const semantic::Index result_node_id,
    const std::vector<CompiledMathBatchLane> &lanes,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane,
    const std::vector<semantic::Index> *schedule_nodes = nullptr,
    const bool workspace_prepared = false) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (result_node_id == semantic::kInvalidIndex || lanes.empty()) {
    return;
  }
  if (batch_workspace == nullptr) {
    throw std::runtime_error("compiled math batch evaluation requires a workspace");
  }
  if (!workspace_prepared) {
    for (const auto &lane : lanes) {
      if (lane.workspace != nullptr) {
        lane.workspace->ensure_size(program);
      }
    }
  }
  compiled_math_batch_prepare_node_values(program, lanes.size(), batch_workspace);
  auto &condition_evaluators = batch_workspace->condition_evaluators;
  condition_evaluators.assign(lanes.size(), nullptr);
  std::vector<double> integral_values;
  const std::size_t time_slot_count =
      compiled_math_program_time_slot_count(program);
  const auto &nodes =
      schedule_nodes == nullptr ? program.root_schedule_nodes : *schedule_nodes;
  bool source_node_batch_base_prepared = false;
  bool source_node_batch_base_supported = false;
  for (semantic::Index schedule_idx = 0; schedule_idx < schedule.size;
       ++schedule_idx) {
    const auto node_id = nodes[
        static_cast<std::size_t>(schedule.offset + schedule_idx)];
    const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
    if (!compiled_math_batch_node_evaluators_for_lanes(
            node,
            lanes,
            &condition_evaluators)) {
      throw std::runtime_error(
          "compiled node requires a planned source view but no workspace was supplied");
    }

    if (compiled_math_is_source_value_node(node.kind) &&
        node.source_vector_program_id != semantic::kInvalidIndex) {
      if (!source_node_batch_base_prepared) {
        source_node_batch_base_supported =
            prepare_compiled_source_node_batch_base(
                program,
                lanes,
                time_slot_count,
                batch_workspace);
        source_node_batch_base_prepared = true;
      }
      if (source_node_batch_base_supported &&
          evaluate_compiled_source_node_batch_prepared(
              program,
              node,
              lanes,
              condition_evaluators,
              time_slot_count,
              batch_workspace)) {
        continue;
      }
      throw std::runtime_error(
          "compiled batch source node has no prepared source-program batch implementation");
    }

    if (compiled_math_is_integral_node(node.kind)) {
      evaluate_compiled_integral_node_batch(
          program,
          node,
          lanes,
          condition_evaluators,
          batch_workspace,
          &integral_values);
      for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
        compiled_math_batch_set_node_value(
            node,
            static_cast<semantic::Index>(lane),
            integral_values[lane],
            batch_workspace);
      }
      continue;
    }

    evaluate_compiled_regular_node_batch(
        program,
        node,
        lanes,
        condition_evaluators,
        time_slot_count,
        batch_workspace);
  }

  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    (*out_by_lane)[lane] =
        compiled_math_batch_node_value(
            result_node_id,
            static_cast<semantic::Index>(lane),
            batch_workspace);
  }
}

inline void evaluate_compiled_math_root_batch(
    const CompiledMathProgram &program,
    const semantic::Index root_id,
    const std::vector<CompiledMathBatchLane> &lanes,
    CompiledMathBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (root_id == semantic::kInvalidIndex || lanes.empty()) {
    return;
  }
  const auto &root = program.roots[static_cast<std::size_t>(root_id)];
  evaluate_compiled_node_span_batch(
      program,
      root.schedule,
      root.node_id,
      lanes,
      batch_workspace,
      out_by_lane,
      nullptr,
      true);
}

} // namespace detail
} // namespace accumulatr::eval
