#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "eval_query.hpp"
#include "exact_sequence.hpp"
#include "observed_plan.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

inline double logsumexp_records(const std::vector<double> &values) {
  if (values.empty()) {
    return R_NegInf;
  }
  double anchor = R_NegInf;
  for (const auto value : values) {
    if (std::isfinite(value) && value > anchor) {
      anchor = value;
    }
  }
  if (!std::isfinite(anchor)) {
    return R_NegInf;
  }
  double sum = 0.0;
  for (const auto value : values) {
    if (std::isfinite(value)) {
      sum += std::exp(value - anchor);
    }
  }
  return sum > 0.0 ? anchor + std::log(sum) : R_NegInf;
}

inline semantic::Index resolve_variant_index_by_component_code(
    const semantic::Index component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code);

struct TrustedParamMatrix {
  const double *base{nullptr};
  int nrow{0};
  int component_weight_start{0};

  TrustedParamMatrix(SEXP paramsSEXP, const int component_weight_param_count)
      : base(REAL(paramsSEXP)),
        nrow(Rf_nrows(paramsSEXP)),
        component_weight_start(
            Rf_ncols(paramsSEXP) - component_weight_param_count) {}

  inline double component_weight(const int row,
                                 const int weight_param_index) const {
    return base[
        static_cast<R_xlen_t>(component_weight_start + weight_param_index) *
            nrow +
        row];
  }

};

inline bool trials_have_observed_components(
    const PreparedTrialLayout &layout,
    const int *component,
    const int *ok = nullptr) {
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    if (!trial_is_selected(ok, trial_index)) {
      continue;
    }
    if (integer_cell_is_na(
            component,
            static_cast<R_xlen_t>(layout.spans[trial_index].start_row))) {
      return false;
    }
  }
  return true;
}

inline std::vector<double> resolve_component_weights(
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &component_codes,
    const TrustedParamMatrix &params,
    const int row) {
  std::vector<double> weights(component_codes.size(), 0.0);
  if (component_codes.empty()) {
    return weights;
  }

  if (model.component_mode != "sample") {
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      weights[i] = model.components[component_index].weight;
    }
  } else {
    const auto reference_id =
        !model.component_reference.empty() && std::any_of(
                                                 component_codes.begin(),
                                                 component_codes.end(),
                                                 [&](const auto code) {
                                                   return model.components
                                                              [static_cast<std::size_t>(code - 1)]
                                                                  .id ==
                                                          model.component_reference;
                                                 })
            ? model.component_reference
            : model.components[static_cast<std::size_t>(component_codes.front() - 1)]
                  .id;
    double sum_nonref = 0.0;
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      const auto &component = model.components[component_index];
      if (component.id == reference_id) {
        continue;
      }
      double weight = component.weight;
      if (component.weight_param_index >= 0) {
        weight = params.component_weight(row, component.weight_param_index);
      }
      weights[i] = weight;
      sum_nonref += weight;
    }
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      const auto &component = model.components[component_index];
      if (component.id != reference_id) {
        continue;
      }
      double ref_weight = 1.0 - sum_nonref;
      weights[i] = ref_weight;
      break;
    }
  }

  double total = 0.0;
  for (const auto weight : weights) {
    total += weight;
  }
  const double inv_total = 1.0 / total;
  for (auto &weight : weights) {
    weight *= inv_total;
  }
  return weights;
}

struct TrialComponentChoice {
  semantic::Index component_code{semantic::kInvalidIndex};
  double weight{0.0};
};

inline bool param_leaf_blocks_equal(SEXP paramsSEXP,
                                    const int lhs_first_row,
                                    const int rhs_first_row,
                                    const std::size_t row_count) {
  if (lhs_first_row == rhs_first_row) {
    return true;
  }
  if (lhs_first_row < 0 || rhs_first_row < 0) {
    return false;
  }
  const int nrow = Rf_nrows(paramsSEXP);
  const int ncol = Rf_ncols(paramsSEXP);
  if (lhs_first_row + static_cast<int>(row_count) > nrow ||
      rhs_first_row + static_cast<int>(row_count) > nrow) {
    return false;
  }
  const double *base = REAL(paramsSEXP);
  for (int col = 0; col < ncol; ++col) {
    const auto col_offset = static_cast<R_xlen_t>(col) * nrow;
    for (std::size_t i = 0; i < row_count; ++i) {
      const auto row_delta = static_cast<R_xlen_t>(i);
      const double lhs =
          base[col_offset + static_cast<R_xlen_t>(lhs_first_row) + row_delta];
      const double rhs =
          base[col_offset + static_cast<R_xlen_t>(rhs_first_row) + row_delta];
      if (lhs != rhs) {
        return false;
      }
    }
  }
  return true;
}

inline std::uint64_t param_leaf_block_hash(SEXP paramsSEXP,
                                           const int first_row,
                                           const std::size_t row_count) {
  if (first_row < 0) {
    return 0U;
  }
  const int nrow = Rf_nrows(paramsSEXP);
  const int ncol = Rf_ncols(paramsSEXP);
  if (first_row + static_cast<int>(row_count) > nrow) {
    return 0U;
  }
  constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
  constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
  std::uint64_t hash = kFnvOffset;
  const double *base = REAL(paramsSEXP);
  for (int col = 0; col < ncol; ++col) {
    const auto col_offset = static_cast<R_xlen_t>(col) * nrow;
    for (std::size_t i = 0; i < row_count; ++i) {
      double value =
          base[col_offset + static_cast<R_xlen_t>(first_row) +
               static_cast<R_xlen_t>(i)];
      if (value == 0.0) {
        value = 0.0;
      }
      std::uint64_t bits = 0U;
      std::memcpy(&bits, &value, sizeof(bits));
      hash ^= bits;
      hash *= kFnvPrime;
    }
  }
  return hash;
}

struct RtFreeObservationPlanCacheKey {
  semantic::Index component_code{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index state_code{semantic::kInvalidIndex};
  std::size_t leaf_count{0U};
  std::uint64_t param_hash{0U};

  bool operator==(const RtFreeObservationPlanCacheKey &other) const noexcept {
    return component_code == other.component_code &&
           variant_index == other.variant_index &&
           state_code == other.state_code &&
           leaf_count == other.leaf_count &&
           param_hash == other.param_hash;
  }
};

struct RtFreeObservationPlanCacheKeyHash {
  std::size_t operator()(const RtFreeObservationPlanCacheKey &key) const noexcept {
    std::uint64_t hash = 1469598103934665603ULL;
    auto mix = [&](const std::uint64_t value) {
      hash ^= value;
      hash *= 1099511628211ULL;
    };
    mix(static_cast<std::uint64_t>(key.component_code));
    mix(static_cast<std::uint64_t>(key.variant_index));
    mix(static_cast<std::uint64_t>(key.state_code));
    mix(static_cast<std::uint64_t>(key.leaf_count));
    mix(key.param_hash);
    return static_cast<std::size_t>(hash);
  }
};

struct RtFreeObservationPlanCacheEntry {
  int first_param_row{0};
  double value{0.0};
};

inline bool rt_free_observation_cache_lookup(
    const std::unordered_map<
        RtFreeObservationPlanCacheKey,
        std::vector<RtFreeObservationPlanCacheEntry>,
        RtFreeObservationPlanCacheKeyHash> &cache,
    SEXP paramsSEXP,
    const RtFreeObservationPlanCacheKey &key,
    const int first_param_row,
    double *value) {
  const auto found = cache.find(key);
  if (found == cache.end()) {
    return false;
  }
  for (const auto &entry : found->second) {
    if (!param_leaf_blocks_equal(
            paramsSEXP,
            entry.first_param_row,
            first_param_row,
            key.leaf_count)) {
      continue;
    }
    *value = entry.value;
    return true;
  }
  return false;
}

inline std::vector<TrialComponentChoice> resolve_trial_components(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const semantic::SemanticModel &model,
    const PreparedTrialLayout &layout,
    const int *component,
    const TrustedParamMatrix &params,
    const std::size_t trial_index) {
  const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
  std::vector<TrialComponentChoice> resolved;

  if (!integer_cell_is_na(component, row)) {
    resolved.push_back(TrialComponentChoice{
        static_cast<semantic::Index>(component[row]),
        1.0});
    return resolved;
  }

  std::vector<semantic::Index> component_codes;
  const auto max_component_code = std::min(
      component_plans_by_code.size(),
      model.components.size() + 1U);
  for (std::size_t component_code = 1; component_code < max_component_code;
       ++component_code) {
    if (component_plans_by_code[component_code].present) {
      component_codes.push_back(static_cast<semantic::Index>(component_code));
    }
  }
  if (component_codes.empty()) {
    return resolved;
  }

  const auto weights = resolve_component_weights(
      model,
      component_codes,
      params,
      static_cast<int>(row));
  resolved.reserve(component_codes.size());
  for (std::size_t i = 0; i < component_codes.size(); ++i) {
    if (!(weights[i] > 0.0)) {
      continue;
    }
    resolved.push_back(TrialComponentChoice{component_codes[i], weights[i]});
  }
  return resolved;
}

inline semantic::Index resolve_variant_index_by_component_code(
    const semantic::Index component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code) {
  const auto idx = static_cast<std::size_t>(component_code);
  return exact_variant_index_by_component_code[idx];
}

inline bool identity_trials_are_all_finite(
    const PreparedTrialLayout &layout,
    const int *label,
    const double *rt,
    const int *ok = nullptr) {
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    if (!trial_is_selected(ok, trial_index)) {
      continue;
    }
    const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
    if (integer_cell_is_na(label, row) || Rcpp::NumericVector::is_na(rt[row])) {
      return false;
    }
  }
  return true;
}

inline double evaluate_observation_plan_direct(
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const ObservationProbabilityPlan &obs_plan,
    const semantic::Index trial_index,
    const semantic::Index variant_index,
    const double observed_rt,
    const double min_ll,
    const int *row_map,
    const int row_offset,
    ExactStepWorkspacePool *workspace_pool,
    std::vector<double> *values) {
  if (values == nullptr || workspace_pool == nullptr) {
    return min_ll;
  }
  values->assign(obs_plan.ops.size(), 0.0);
  const auto &exact_plan = exact_plans[static_cast<std::size_t>(variant_index)];
  ParamView params(paramsSEXP, row_map, row_offset);
  const int first_param_row =
      row_map == nullptr
          ? static_cast<int>(
                layout.spans[static_cast<std::size_t>(trial_index)].start_row)
          : 0;
  auto &workspace = workspace_pool->get(exact_plans, variant_index);

  for (std::size_t op_index = 0; op_index < obs_plan.ops.size(); ++op_index) {
    const auto &op = obs_plan.ops[op_index];
    double value = 0.0;
    switch (op.kind) {
    case ObservationPlanOpKind::Constant:
      value = op.constant;
      break;
    case ObservationPlanOpKind::LogDensity: {
      const auto target_idx =
          exact_plan.outcome_index_by_code[static_cast<std::size_t>(
              op.semantic_code)];
      const double density = exact_unranked_target_density(
          exact_plan,
          params,
          first_param_row,
          target_idx,
          observed_rt,
          &workspace);
      value =
          std::isfinite(density) && density > 0.0 && op.weight > 0.0
              ? std::log(op.weight) + std::log(density)
              : min_ll;
      break;
    }
    case ObservationPlanOpKind::FiniteOutcomeProbability: {
      const auto target_idx =
          exact_plan.outcome_index_by_code[static_cast<std::size_t>(
              op.semantic_code)];
      const double probability = exact_finite_probability_program_value(
          exact_plan,
          params,
          first_param_row,
          target_idx,
          &workspace);
      value =
          std::isfinite(probability) && probability > 0.0 && op.weight > 0.0
              ? op.weight * probability
              : 0.0;
      break;
    }
    case ObservationPlanOpKind::NoResponseProbability: {
      const double probability = exact_no_response_program_value(
          exact_plan,
          params,
          first_param_row,
          &workspace);
      value =
          std::isfinite(probability) && probability > 0.0 && op.weight > 0.0
              ? op.weight * probability
              : 0.0;
      break;
    }
    case ObservationPlanOpKind::WeightedSum:
      if (op.value_kind == ObservationPlanValueKind::Log) {
        double anchor = R_NegInf;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = obs_plan.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double child_value =
              (*values)[static_cast<std::size_t>(child)];
          if (std::isfinite(child_value) && child_value > anchor) {
            anchor = child_value;
          }
        }
        if (!std::isfinite(anchor)) {
          value = min_ll;
          break;
        }
        double sum = 0.0;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = obs_plan.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double child_value =
              (*values)[static_cast<std::size_t>(child)];
          if (std::isfinite(child_value)) {
            sum += std::exp(child_value - anchor);
          }
        }
        value = sum > 0.0 ? anchor + std::log(sum) : min_ll;
      } else {
        double sum = 0.0;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = obs_plan.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double child_value =
              (*values)[static_cast<std::size_t>(child)];
          if (std::isfinite(child_value) && child_value > 0.0) {
            sum += child_value;
          }
        }
        value = sum;
      }
      break;
    case ObservationPlanOpKind::Complement: {
      double sum = 0.0;
      for (semantic::Index i = 0; i < op.children.size; ++i) {
        const auto child = obs_plan.child_ops[
            static_cast<std::size_t>(op.children.offset + i)];
        if (child == semantic::kInvalidIndex) {
          continue;
        }
        const double child_value =
            (*values)[static_cast<std::size_t>(child)];
        if (std::isfinite(child_value)) {
          sum += child_value;
        }
      }
      value = std::max(0.0, 1.0 - sum);
      break;
    }
    case ObservationPlanOpKind::Log: {
      double probability = 0.0;
      if (op.children.size > 0) {
        const auto child =
            obs_plan.child_ops[static_cast<std::size_t>(op.children.offset)];
        if (child != semantic::kInvalidIndex) {
          probability = (*values)[static_cast<std::size_t>(child)];
        }
      }
      value =
          std::isfinite(probability) && probability > 0.0
              ? std::log(probability)
              : min_ll;
      break;
    }
    }
    (*values)[op_index] = value;
  }

  return obs_plan.root == semantic::kInvalidIndex
             ? min_ll
             : (*values)[static_cast<std::size_t>(obs_plan.root)];
}

struct ObservationBatchLane {
  std::size_t trial_index{0};
  double component_weight{1.0};
  double observed_rt{NA_REAL};
  const int *row_map{nullptr};
  int row_offset{0};
  int first_param_row{0};
};

struct ObservationBatchGroup {
  semantic::Index variant_index{semantic::kInvalidIndex};
  const ObservationProbabilityPlan *plan{nullptr};
  std::vector<ObservationBatchLane> lanes;
};

struct ObservationBatchWorkspace {
  semantic::Index workspace_variant_index{semantic::kInvalidIndex};
  std::vector<std::unique_ptr<ExactStepWorkspace>> step_workspaces;
  std::vector<ExactProbabilityBatchLane> exact_lanes;
  ExactProbabilityBatchWorkspace exact_workspace;
  std::vector<double> exact_values;
  std::vector<double> unique_exact_values;
  std::vector<double> op_values;
  std::vector<double> child_values;
};

inline ObservationBatchGroup &find_or_add_observation_batch_group(
    std::vector<ObservationBatchGroup> *groups,
    const semantic::Index variant_index,
    const ObservationProbabilityPlan &plan) {
  for (auto &group : *groups) {
    if (group.variant_index == variant_index && group.plan == &plan) {
      return group;
    }
  }
  groups->push_back(ObservationBatchGroup{variant_index, &plan, {}});
  return groups->back();
}

inline void observation_batch_workspace_ensure_step_workspaces(
    ObservationBatchWorkspace *workspace,
    const ExactVariantPlan &plan,
    const semantic::Index variant_index,
    const std::size_t lane_count) {
  if (workspace->workspace_variant_index != variant_index) {
    workspace->step_workspaces.clear();
    workspace->workspace_variant_index = variant_index;
  }
  while (workspace->step_workspaces.size() < lane_count) {
    workspace->step_workspaces.push_back(
        std::make_unique<ExactStepWorkspace>(plan));
  }
}

inline void evaluate_exact_program_for_observation_batch(
    const ExactVariantPlan &exact_plan,
    const ExactProbabilityProgram &program,
    SEXP paramsSEXP,
    const ObservationBatchGroup &group,
    const double observed_time,
    const bool use_lane_observed_time,
    ObservationBatchWorkspace *workspace,
    std::vector<double> *out_by_lane) {
  const auto lane_count = group.lanes.size();
  out_by_lane->assign(lane_count, 0.0);
  const auto leaf_count =
      static_cast<std::size_t>(exact_plan.lowered.program.layout.n_leaves);
  if (!use_lane_observed_time && !std::isfinite(observed_time) &&
      lane_count > 1U) {
    bool can_deduplicate = true;
    std::vector<std::size_t> unique_lanes;
    std::vector<std::uint64_t> unique_hashes;
    std::vector<std::size_t> unique_index_by_lane;
    unique_lanes.reserve(lane_count);
    unique_hashes.reserve(lane_count);
    unique_index_by_lane.reserve(lane_count);
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const auto &batch_lane = group.lanes[lane];
      if (batch_lane.row_map != nullptr) {
        can_deduplicate = false;
        break;
      }
      const auto hash =
          param_leaf_block_hash(paramsSEXP, batch_lane.first_param_row, leaf_count);
      std::size_t unique_index = unique_lanes.size();
      for (std::size_t i = 0; i < unique_lanes.size(); ++i) {
        if (unique_hashes[i] != hash) {
          continue;
        }
        const auto &unique_lane = group.lanes[unique_lanes[i]];
        if (param_leaf_blocks_equal(
                paramsSEXP,
                unique_lane.first_param_row,
                batch_lane.first_param_row,
                leaf_count)) {
          unique_index = i;
          break;
        }
      }
      if (unique_index == unique_lanes.size()) {
        unique_lanes.push_back(lane);
        unique_hashes.push_back(hash);
      }
      unique_index_by_lane.push_back(unique_index);
    }
    if (can_deduplicate && unique_lanes.size() < lane_count) {
      workspace->exact_lanes.clear();
      workspace->exact_lanes.reserve(unique_lanes.size());
      for (const auto lane : unique_lanes) {
        const auto &batch_lane = group.lanes[lane];
        workspace->exact_lanes.emplace_back(
            paramsSEXP,
            nullptr,
            0,
            batch_lane.first_param_row,
            observed_time,
            workspace->step_workspaces[lane].get());
      }
      evaluate_exact_probability_program_batch(
          exact_plan,
          program,
          workspace->exact_lanes,
          &workspace->exact_workspace,
          &workspace->unique_exact_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        (*out_by_lane)[lane] =
            workspace->unique_exact_values[unique_index_by_lane[lane]];
      }
      return;
    }
  }
  workspace->exact_lanes.clear();
  workspace->exact_lanes.reserve(lane_count);
  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    const auto &batch_lane = group.lanes[lane];
    workspace->exact_lanes.emplace_back(
        paramsSEXP,
        batch_lane.row_map,
        batch_lane.row_offset,
        batch_lane.first_param_row,
        use_lane_observed_time ? batch_lane.observed_rt : observed_time,
        workspace->step_workspaces[lane].get());
  }
  evaluate_exact_probability_program_batch(
      exact_plan,
      program,
      workspace->exact_lanes,
      &workspace->exact_workspace,
      out_by_lane);
}

inline void evaluate_observation_batch_group(
    const std::vector<ExactVariantPlan> &exact_plans,
    SEXP paramsSEXP,
    const double min_ll,
    ObservationBatchGroup *group,
    ObservationBatchWorkspace *workspace,
    std::vector<double> *out_by_lane) {
  const auto lane_count = group->lanes.size();
  out_by_lane->assign(lane_count, min_ll);
  if (group->plan == nullptr ||
      group->variant_index == semantic::kInvalidIndex ||
      lane_count == 0U) {
    return;
  }
  const auto &obs_plan = *group->plan;
  if (obs_plan.empty()) {
    return;
  }
  const auto &exact_plan =
      exact_plans[static_cast<std::size_t>(group->variant_index)];
  observation_batch_workspace_ensure_step_workspaces(
      workspace,
      exact_plan,
      group->variant_index,
      lane_count);
  auto &op_values = workspace->op_values;
  op_values.assign(obs_plan.ops.size() * lane_count, 0.0);
  auto op_value = [&](const std::size_t op_index,
                      const std::size_t lane) -> double & {
    return op_values[op_index * lane_count + lane];
  };

  for (std::size_t op_index = 0; op_index < obs_plan.ops.size(); ++op_index) {
    const auto &op = obs_plan.ops[op_index];
    switch (op.kind) {
    case ObservationPlanOpKind::Constant:
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        op_value(op_index, lane) = op.constant;
      }
      break;

    case ObservationPlanOpKind::LogDensity: {
      const auto target_idx =
          exact_plan.outcome_index_by_code[static_cast<std::size_t>(
              op.semantic_code)];
      if (target_idx == semantic::kInvalidIndex || !(op.weight > 0.0)) {
        for (std::size_t lane = 0; lane < lane_count; ++lane) {
          op_value(op_index, lane) = min_ll;
        }
        break;
      }
      const auto program_index =
          exact_plan.probability_programs.density_by_outcome[
              static_cast<std::size_t>(target_idx)];
      const auto &program =
          exact_plan.probability_programs.programs[
              static_cast<std::size_t>(program_index)];
      evaluate_exact_program_for_observation_batch(
          exact_plan,
          program,
          paramsSEXP,
          *group,
          NA_REAL,
          true,
          workspace,
          &workspace->exact_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        const double rt = group->lanes[lane].observed_rt;
        const double density = workspace->exact_values[lane];
        op_value(op_index, lane) =
            std::isfinite(rt) && rt > 0.0 &&
                    std::isfinite(density) && density > 0.0
                ? std::log(op.weight) + std::log(density)
                : min_ll;
      }
      break;
    }

    case ObservationPlanOpKind::FiniteOutcomeProbability: {
      const auto target_idx =
          exact_plan.outcome_index_by_code[static_cast<std::size_t>(
              op.semantic_code)];
      if (target_idx == semantic::kInvalidIndex || !(op.weight > 0.0)) {
        break;
      }
      const auto program_index =
          exact_plan.probability_programs.finite_probability_by_outcome[
              static_cast<std::size_t>(target_idx)];
      const auto &program =
          exact_plan.probability_programs.programs[
              static_cast<std::size_t>(program_index)];
      evaluate_exact_program_for_observation_batch(
          exact_plan,
          program,
          paramsSEXP,
          *group,
          NA_REAL,
          false,
          workspace,
          &workspace->exact_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        const double probability = workspace->exact_values[lane];
        op_value(op_index, lane) =
            std::isfinite(probability) && probability > 0.0
                ? op.weight * probability
                : 0.0;
      }
      break;
    }

    case ObservationPlanOpKind::NoResponseProbability: {
      if (!(op.weight > 0.0)) {
        break;
      }
      const auto &program =
          exact_plan.probability_programs.programs[
              static_cast<std::size_t>(
                  exact_plan.probability_programs.no_response_probability)];
      evaluate_exact_program_for_observation_batch(
          exact_plan,
          program,
          paramsSEXP,
          *group,
          NA_REAL,
          false,
          workspace,
          &workspace->exact_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        const double probability = workspace->exact_values[lane];
        op_value(op_index, lane) =
            std::isfinite(probability) && probability > 0.0
                ? op.weight * probability
                : 0.0;
      }
      break;
    }

    case ObservationPlanOpKind::WeightedSum:
      if (op.value_kind == ObservationPlanValueKind::Log) {
        for (std::size_t lane = 0; lane < lane_count; ++lane) {
          double anchor = R_NegInf;
          for (semantic::Index i = 0; i < op.children.size; ++i) {
            const auto child = obs_plan.child_ops[
                static_cast<std::size_t>(op.children.offset + i)];
            if (child == semantic::kInvalidIndex) {
              continue;
            }
            const double child_value =
                op_value(static_cast<std::size_t>(child), lane);
            if (std::isfinite(child_value) && child_value > anchor) {
              anchor = child_value;
            }
          }
          if (!std::isfinite(anchor)) {
            op_value(op_index, lane) = min_ll;
            continue;
          }
          double sum = 0.0;
          for (semantic::Index i = 0; i < op.children.size; ++i) {
            const auto child = obs_plan.child_ops[
                static_cast<std::size_t>(op.children.offset + i)];
            if (child == semantic::kInvalidIndex) {
              continue;
            }
            const double child_value =
                op_value(static_cast<std::size_t>(child), lane);
            if (std::isfinite(child_value)) {
              sum += std::exp(child_value - anchor);
            }
          }
          op_value(op_index, lane) =
              sum > 0.0 ? anchor + std::log(sum) : min_ll;
        }
      } else {
        for (std::size_t lane = 0; lane < lane_count; ++lane) {
          double sum = 0.0;
          for (semantic::Index i = 0; i < op.children.size; ++i) {
            const auto child = obs_plan.child_ops[
                static_cast<std::size_t>(op.children.offset + i)];
            if (child == semantic::kInvalidIndex) {
              continue;
            }
            const double child_value =
                op_value(static_cast<std::size_t>(child), lane);
            if (std::isfinite(child_value) && child_value > 0.0) {
              sum += child_value;
            }
          }
          op_value(op_index, lane) = sum;
        }
      }
      break;

    case ObservationPlanOpKind::Complement:
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        double sum = 0.0;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = obs_plan.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double child_value =
              op_value(static_cast<std::size_t>(child), lane);
          if (std::isfinite(child_value)) {
            sum += child_value;
          }
        }
        op_value(op_index, lane) = std::max(0.0, 1.0 - sum);
      }
      break;

    case ObservationPlanOpKind::Log:
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        double probability = 0.0;
        if (op.children.size > 0) {
          const auto child =
              obs_plan.child_ops[static_cast<std::size_t>(op.children.offset)];
          if (child != semantic::kInvalidIndex) {
            probability = op_value(static_cast<std::size_t>(child), lane);
          }
        }
        op_value(op_index, lane) =
            std::isfinite(probability) && probability > 0.0
                ? std::log(probability)
                : min_ll;
      }
      break;
    }
  }

  if (obs_plan.root == semantic::kInvalidIndex) {
    return;
  }
  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    (*out_by_lane)[lane] =
        op_value(static_cast<std::size_t>(obs_plan.root), lane);
  }
}

inline Rcpp::List evaluate_outcome_queries_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP) {
  const auto table = read_prepared_data_view(dataSEXP, layout);
  const int *label =
      INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector probability(n_trials, 0.0);
  ExactStepWorkspacePool workspace_pool(exact_plans.size());
  std::vector<double> plan_values;

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
    if (integer_cell_is_na(table.component, row) || integer_cell_is_na(label, row)) {
      continue;
    }

    const auto component_code =
        static_cast<semantic::Index>(table.component[row]);
    if (component_code <= 0 ||
        component_code >=
            static_cast<semantic::Index>(component_plans_by_code.size())) {
      continue;
    }

    const auto &component_plan =
        component_plans_by_code[static_cast<std::size_t>(component_code)];
    if (!component_plan.present) {
      continue;
    }

    const auto variant_index = resolve_variant_index_by_component_code(
        component_code,
        exact_variant_index_by_component_code);
    if (variant_index == semantic::kInvalidIndex) {
      continue;
    }

    const auto observed_code =
        static_cast<semantic::Index>(label[row]);
    if (observed_code <= 0) {
      continue;
    }
    const auto state_code =
        missing_rt_observation_state_code(component_plan, observed_code);
    const auto &plan =
        observation_probability_plan_for_state(component_plan, state_code);
    const double value = evaluate_observation_plan_direct(
        exact_plans,
        layout,
        paramsSEXP,
        plan,
        static_cast<semantic::Index>(trial_index),
        variant_index,
        NA_REAL,
        std::log(1e-10),
        nullptr,
        0,
        &workspace_pool,
        &plan_values);
    if (std::isfinite(value) && value > 0.0) {
      probability[static_cast<R_xlen_t>(trial_index)] = value;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = probability,
      Rcpp::Named("n_trials") = static_cast<int>(n_trials));
}

inline SEXP evaluate_observed_trials_cached_impl(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const bool observed_identity,
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    SEXP expandSEXP,
    const int *ok = nullptr) {
  BatchSourceProductDebugScope batch_source_product_debug_scope;
  BatchFiniteIntegralDebugScope batch_finite_integral_debug_scope;
  const auto table = read_prepared_data_view(dataSEXP, layout);

  if (observed_identity) {
    const int *label =
        INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
    const double *rt =
        REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
    if (trials_have_observed_components(layout, table.component, ok)) {
      if (identity_trials_are_all_finite(layout, label, rt, ok)) {
        SEXP result = evaluate_exact_trials_cached(
            exact_variant_index_by_component_code,
            exact_plans,
            layout,
            paramsSEXP,
            dataSEXP,
            min_ll,
            expandSEXP,
            ok);
        return result;
      }
    }
  }

  const int *label =
      INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
  const double *rt =
      REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
  const TrustedParamMatrix trusted_params(
      paramsSEXP,
      model.component_weight_param_count);

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector loglik(n_trials, min_ll);
  std::vector<std::vector<double>> trial_values(n_trials);
  std::vector<ObservationBatchGroup> observation_groups;

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto &span = layout.spans[trial_index];
    const auto row = static_cast<R_xlen_t>(span.start_row);
    if (!trial_is_selected(ok, trial_index)) {
      continue;
    }
    const auto observed_label_code =
        integer_cell_is_na(label, row)
            ? semantic::kInvalidIndex
            : static_cast<semantic::Index>(label[row]);
    const double observed_rt = rt[row];

    const auto components = resolve_trial_components(
        component_plans_by_code,
        model,
        layout,
        table.component,
        trusted_params,
        trial_index);
    const auto latent_trial = integer_cell_is_na(table.component, row);
    for (const auto &choice : components) {
      if (choice.component_code <= 0 ||
          choice.component_code >=
              static_cast<semantic::Index>(component_plans_by_code.size())) {
        continue;
      }
      const auto &component_plan =
          component_plans_by_code[static_cast<std::size_t>(choice.component_code)];
      if (!component_plan.present) {
        continue;
      }
      const auto variant_index = resolve_variant_index_by_component_code(
          choice.component_code,
          exact_variant_index_by_component_code);
      if (variant_index == semantic::kInvalidIndex) {
        continue;
      }

      const auto state_code = observation_state_code(
          component_plan,
          observed_label_code,
          observed_rt);
      if (state_code == semantic::kInvalidIndex) {
        continue;
      }
      const auto &plan =
          observation_log_plan_for_state(component_plan, state_code);
      const double plan_rt =
          observation_state_uses_rt(component_plan, state_code)
              ? observed_rt
              : NA_REAL;

      const int *row_map_ptr = nullptr;
      int row_offset = 0;
      const auto &exact_plan =
          exact_plans[static_cast<std::size_t>(variant_index)];
      (void)exact_plan;
      if (latent_trial) {
        const auto &leaf_offsets =
            exact_plan.leaf_row_offsets;
        row_map_ptr = leaf_offsets.data();
        row_offset = static_cast<int>(span.start_row);
      }
      const int first_param_row =
          row_map_ptr == nullptr ? static_cast<int>(span.start_row) : 0;
      auto &group = find_or_add_observation_batch_group(
          &observation_groups,
          variant_index,
          plan);
      group.lanes.push_back(ObservationBatchLane{
          trial_index,
          choice.weight,
          plan_rt,
          row_map_ptr,
          row_offset,
          first_param_row});
    }
  }

  ObservationBatchWorkspace batch_workspace;
  std::vector<double> group_values;
  for (auto &group : observation_groups) {
    evaluate_observation_batch_group(
        exact_plans,
        paramsSEXP,
        min_ll,
        &group,
        &batch_workspace,
        &group_values);
    for (std::size_t lane = 0; lane < group.lanes.size(); ++lane) {
      const auto &batch_lane = group.lanes[lane];
      const double value = group_values[lane];
      if (std::isfinite(value) && batch_lane.component_weight > 0.0) {
        trial_values[batch_lane.trial_index].push_back(
            std::log(batch_lane.component_weight) + value);
      }
    }
  }

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto value = logsumexp_records(trial_values[trial_index]);
    loglik[static_cast<R_xlen_t>(trial_index)] =
        std::isfinite(value) ? value : min_ll;
  }

  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    if (!std::isfinite(loglik[i])) {
      loglik[i] = min_ll;
    }
  }

  const double total_loglik = aggregate_trial_loglik(loglik, expandSEXP);

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("n_trials") = static_cast<int>(n_trials));
}

inline SEXP evaluate_observed_trials_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const bool observed_identity,
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    SEXP expandSEXP,
    const int *ok = nullptr) {
  return evaluate_observed_trials_cached_impl(
      component_plans_by_code,
      observed_identity,
      model,
      exact_variant_index_by_component_code,
      exact_plans,
      layout,
      paramsSEXP,
      dataSEXP,
      min_ll,
      expandSEXP,
      ok);
}

} // namespace detail
} // namespace accumulatr::eval
