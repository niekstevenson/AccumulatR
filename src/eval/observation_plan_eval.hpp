#pragma once

#include <Rcpp.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "exact_sequence.hpp"
#include "observation_model.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ObservationInterval {
  double lower{0.0};
  double upper{R_PosInf};
  bool include_terminal_no_response{false};
};

inline bool param_leaf_blocks_equal(SEXP paramsSEXP,
                                    const double *onset,
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
  if (onset != nullptr) {
    for (std::size_t i = 0; i < row_count; ++i) {
      const auto row_delta = static_cast<R_xlen_t>(i);
      if (onset[static_cast<R_xlen_t>(lhs_first_row) + row_delta] !=
          onset[static_cast<R_xlen_t>(rhs_first_row) + row_delta]) {
        return false;
      }
    }
  }
  return true;
}

inline std::uint64_t param_leaf_block_hash(SEXP paramsSEXP,
                                           const double *onset,
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
  if (onset != nullptr) {
    for (std::size_t i = 0; i < row_count; ++i) {
      double value =
          onset[static_cast<R_xlen_t>(first_row) + static_cast<R_xlen_t>(i)];
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
    const double *onset,
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
            onset,
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

inline double evaluate_observed_label_interval_probability(
    const std::vector<ExactVariantPlan> &exact_plans,
    SEXP paramsSEXP,
    const double *onset,
    const ComponentObservationPlan &component_plan,
    const semantic::Index observed_code,
    const semantic::Index variant_index,
    const ObservationInterval &interval,
    const int *row_map,
    const int row_offset,
    const int first_param_row,
    ExactStepWorkspacePool *workspace_pool) {
  if (workspace_pool == nullptr ||
      observed_code == semantic::kInvalidIndex ||
      static_cast<std::size_t>(observed_code) >=
          component_plan.interval_by_code.size()) {
    return 0.0;
  }
  const auto &branches =
      component_plan.interval_by_code[static_cast<std::size_t>(observed_code)];
  if (branches.empty()) {
    return 0.0;
  }
  const auto &exact_plan = exact_plans[static_cast<std::size_t>(variant_index)];
  ParamView params(paramsSEXP, onset, row_map, row_offset);
  auto &workspace = workspace_pool->get(exact_plans, variant_index);
  double total = 0.0;
  for (const auto &branch : branches) {
    if (!(branch.weight > 0.0) ||
        branch.semantic_code == semantic::kInvalidIndex ||
        static_cast<std::size_t>(branch.semantic_code) >=
            exact_plan.outcome_index_by_code.size()) {
      continue;
    }
    const auto target_idx =
        exact_plan.outcome_index_by_code[
            static_cast<std::size_t>(branch.semantic_code)];
    total += branch.weight *
             exact_outcome_probability_between(
                 exact_plan,
                 params,
                 first_param_row,
                 target_idx,
                 interval.lower,
                 interval.upper,
                 &workspace);
  }
  return std::isfinite(total) ? clamp_probability(total) : 0.0;
}

inline double evaluate_unknown_label_interval_probability(
    const std::vector<ExactVariantPlan> &exact_plans,
    SEXP paramsSEXP,
    const double *onset,
    const ComponentObservationPlan &component_plan,
    const semantic::Index variant_index,
    const ObservationInterval &interval,
    const int *row_map,
    const int row_offset,
    const int first_param_row,
    ExactStepWorkspacePool *workspace_pool) {
  if (workspace_pool == nullptr) {
    return 0.0;
  }
  const auto &exact_plan = exact_plans[static_cast<std::size_t>(variant_index)];
  ParamView params(paramsSEXP, onset, row_map, row_offset);
  auto &workspace = workspace_pool->get(exact_plans, variant_index);
  double total = exact_all_outcome_probability_between(
      exact_plan,
      params,
      first_param_row,
      component_plan.finite_outcome_codes,
      interval.lower,
      interval.upper,
      &workspace);
  if (interval.include_terminal_no_response) {
    total += exact_terminal_no_response_probability(
        exact_plan,
        params,
        first_param_row);
  }
  return std::isfinite(total) ? clamp_probability(total) : 0.0;
}

inline double evaluate_observation_plan_at_row(
    const std::vector<ExactVariantPlan> &exact_plans,
    SEXP paramsSEXP,
    const double *onset,
    const ObservationProbabilityPlan &obs_plan,
    const semantic::Index variant_index,
    const double observed_rt,
    const double min_ll,
    const int *row_map,
    const int row_offset,
    const int first_param_row,
    ExactStepWorkspacePool *workspace_pool,
    std::vector<double> *values) {
  if (values == nullptr || workspace_pool == nullptr) {
    return min_ll;
  }
  values->assign(obs_plan.ops.size(), 0.0);
  const auto &exact_plan = exact_plans[static_cast<std::size_t>(variant_index)];
  ParamView params(paramsSEXP, onset, row_map, row_offset);
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
      const double probability = exact_finite_outcome_probability(
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
    case ObservationPlanOpKind::NoResponseProbability:
      value = exact_terminal_no_response_probability(
          exact_plan,
          params,
          first_param_row);
      break;
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

inline double evaluate_observation_plan_direct(
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    const double *onset,
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
  const int first_param_row =
      row_map == nullptr
          ? static_cast<int>(
                layout.trials[static_cast<std::size_t>(trial_index)].start_row)
          : 0;
  return evaluate_observation_plan_at_row(
      exact_plans,
      paramsSEXP,
      onset,
      obs_plan,
      variant_index,
      observed_rt,
      min_ll,
      row_map,
      row_offset,
      first_param_row,
      workspace_pool,
      values);
}

} // namespace detail
} // namespace accumulatr::eval
