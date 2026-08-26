#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "exact_sequence.hpp"
#include "lane_math.hpp"
#include "observation_model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ObservationLaneGroup {
  ObservationLaneGroup(const semantic::Index variant_index_,
                       const ObservationProbabilityPlan *plan_)
      : variant_index(variant_index_), plan(plan_) {}

  void clear() {
    lanes.clear();
    destination_indices.clear();
    component_weights.clear();
  }

  void reserve(const std::size_t size) {
    lanes.reserve(size);
    destination_indices.reserve(size);
    component_weights.reserve(size);
  }

  semantic::Index variant_index{semantic::kInvalidIndex};
  const ObservationProbabilityPlan *plan{nullptr};
  ObservationLaneBatch lanes;
  std::vector<std::size_t> destination_indices;
  std::vector<double> component_weights;
};

struct ObservationLaneWorkspace {
  explicit ObservationLaneWorkspace(const std::size_t plan_count)
      : exact_workspaces(plan_count) {}

  void begin() {
    for (auto &group : groups) {
      group.clear();
    }
  }

  ExactStepLaneWorkspacePool exact_workspaces;
  std::vector<double> exact_values;
  std::vector<double> op_values;
  std::vector<double> reduction_values;
  std::vector<ObservationLaneGroup> groups;
  std::vector<double> group_values;
  std::vector<semantic::Index> component_codes;
  std::vector<double> component_weights;
};

inline ObservationLaneGroup &find_observation_lane_group(
    std::vector<ObservationLaneGroup> *groups,
    const semantic::Index variant_index,
    const ObservationProbabilityPlan &plan) {
  for (auto &group : *groups) {
    if (group.variant_index == variant_index && group.plan == &plan) {
      return group;
    }
  }
  groups->emplace_back(variant_index, &plan);
  groups->back().reserve(kExactLaneTileSize);
  return groups->back();
}

inline void finish_log_density_lanes(
    const double *density,
    const std::size_t lane_count,
    const double weight,
    const double min_ll,
    double *out) {
  log_lanes(density, out, lane_count);
  const double log_weight = weight == 1.0 ? 0.0 : std::log(weight);
  for (std::size_t lane = 0U; lane < lane_count; ++lane) {
    out[lane] = std::isfinite(out[lane])
                    ? log_weight + out[lane]
                    : min_ll;
  }
}

inline void evaluate_observation_lanes(
    const std::vector<ExactVariantPlan> &exact_plans,
    const double min_ll,
    const semantic::Index variant_index,
    const ObservationProbabilityPlan &observation,
    const ObservationLaneBatchView lanes,
    ObservationLaneWorkspace *workspace,
    std::vector<double> *out) {
  const auto lane_count = lanes.size;
  if (observation.empty() || variant_index == semantic::kInvalidIndex ||
      lane_count == 0U) {
    out->assign(lane_count, min_ll);
    return;
  }
  out->resize(lane_count);
  const auto &exact_plan =
      exact_plans[static_cast<std::size_t>(variant_index)];
  auto &exact_workspace = workspace->exact_workspaces.get(
      exact_plans, variant_index);
  const auto root_index = static_cast<std::size_t>(observation.root);
  workspace->op_values.resize((observation.ops.size() - 1U) * lane_count);
  const auto op_values = [&](const std::size_t op_index) {
    if (op_index == root_index) {
      return out->data();
    }
    const auto storage_index = op_index < root_index
                                   ? op_index
                                   : op_index - 1U;
    return workspace->op_values.data() + storage_index * lane_count;
  };

  for (std::size_t op_index = 0;
       op_index < observation.ops.size();
       ++op_index) {
    const auto &op = observation.ops[op_index];
    double *values = op_values(op_index);
    switch (op.kind) {
    case ObservationPlanOpKind::Constant:
      std::fill_n(values, lane_count, op.constant);
      break;
    case ObservationPlanOpKind::LogDensity: {
      const auto target = exact_plan.outcome_index_by_code[
          static_cast<std::size_t>(op.semantic_code)];
      exact_unranked_target_density_lanes(
          exact_plan,
          lanes,
          target,
          &exact_workspace,
          &workspace->exact_values);
      finish_log_density_lanes(
          workspace->exact_values.data(),
          lane_count,
          op.weight,
          min_ll,
          values);
      break;
    }
    case ObservationPlanOpKind::FiniteOutcomeProbability: {
      const auto target = exact_plan.outcome_index_by_code[
          static_cast<std::size_t>(op.semantic_code)];
      exact_finite_outcome_probability_lanes(
          exact_plan,
          lanes,
          target,
          &exact_workspace,
          &workspace->exact_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        const double probability = workspace->exact_values[lane];
        values[lane] = std::isfinite(probability) && probability > 0.0
                           ? op.weight * probability
                           : 0.0;
      }
      break;
    }
    case ObservationPlanOpKind::NoResponseProbability:
      exact_terminal_no_response_probability_lanes(
          exact_plan,
          lanes,
          &exact_workspace,
          &workspace->exact_values);
      std::copy(
          workspace->exact_values.begin(),
          workspace->exact_values.end(),
          values);
      break;
    case ObservationPlanOpKind::WeightedSum:
      if (op.value_kind == ObservationPlanValueKind::Log) {
        std::fill_n(values, lane_count, R_NegInf);
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = observation.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double *child_values =
              op_values(static_cast<std::size_t>(child));
          for (std::size_t lane = 0; lane < lane_count; ++lane) {
            if (std::isfinite(child_values[lane]) &&
                child_values[lane] > values[lane]) {
              values[lane] = child_values[lane];
            }
          }
        }
        workspace->reduction_values.assign(lane_count, 0.0);
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = observation.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double *child_values =
              op_values(static_cast<std::size_t>(child));
          for (std::size_t lane = 0; lane < lane_count; ++lane) {
            if (std::isfinite(values[lane]) &&
                std::isfinite(child_values[lane])) {
              workspace->reduction_values[lane] +=
                  std::exp(child_values[lane] - values[lane]);
            }
          }
        }
        for (std::size_t lane = 0; lane < lane_count; ++lane) {
          const double sum = workspace->reduction_values[lane];
          values[lane] = std::isfinite(values[lane]) && sum > 0.0
                             ? values[lane] + std::log(sum)
                             : min_ll;
        }
      } else {
        std::fill_n(values, lane_count, 0.0);
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = observation.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double *child_values =
              op_values(static_cast<std::size_t>(child));
          for (std::size_t lane = 0; lane < lane_count; ++lane) {
            if (std::isfinite(child_values[lane]) &&
                child_values[lane] > 0.0) {
              values[lane] += child_values[lane];
            }
          }
        }
      }
      break;
    case ObservationPlanOpKind::Complement:
      std::fill_n(values, lane_count, 1.0);
      for (semantic::Index i = 0; i < op.children.size; ++i) {
        const auto child = observation.child_ops[
            static_cast<std::size_t>(op.children.offset + i)];
        if (child == semantic::kInvalidIndex) {
          continue;
        }
        const double *child_values =
            op_values(static_cast<std::size_t>(child));
        for (std::size_t lane = 0; lane < lane_count; ++lane) {
          if (std::isfinite(child_values[lane])) {
            values[lane] -= child_values[lane];
          }
        }
      }
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        values[lane] = std::max(0.0, values[lane]);
      }
      break;
    case ObservationPlanOpKind::Log: {
      const double *probability = nullptr;
      if (op.children.size > 0) {
        const auto child = observation.child_ops[
            static_cast<std::size_t>(op.children.offset)];
        if (child != semantic::kInvalidIndex) {
          probability = op_values(static_cast<std::size_t>(child));
        }
      }
      if (probability == nullptr) {
        std::fill_n(values, lane_count, min_ll);
        break;
      }
      log_lanes(probability, values, lane_count);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        values[lane] = std::isfinite(values[lane])
                           ? values[lane]
                           : min_ll;
      }
      break;
    }
    }
  }
}

inline void evaluate_observation_lane_group(
    const std::vector<ExactVariantPlan> &exact_plans,
    const double min_ll,
    const ObservationLaneGroup &group,
    const ParamMatrixView &parameter_matrix,
    ObservationLaneWorkspace *workspace,
    std::vector<double> *out) {
  if (group.plan == nullptr) {
    out->assign(group.lanes.size(), min_ll);
    return;
  }
  evaluate_observation_lanes(
      exact_plans,
      min_ll,
      group.variant_index,
      *group.plan,
      group.lanes.view(parameter_matrix),
      workspace,
      out);
}

} // namespace detail
} // namespace accumulatr::eval
