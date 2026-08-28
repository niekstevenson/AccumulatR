#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "exact_adaptive.hpp"
#include "exact_sequence.hpp"

namespace accumulatr::eval {
namespace detail {

inline constexpr double kExactResponseTimeLimit = 30.0;

struct ExactOutcomeTerm {
  semantic::Index target{semantic::kInvalidIndex};
  double weight{1.0};
};

enum class ExactResponseMeasure : std::uint8_t {
  SelectedResponse = 0,
  ObservableResponses = 1
};

struct ExactIntervalLaneWorkspace {
  AdaptiveLaneWorkspace adaptive;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> endpoint_times;
  std::vector<std::size_t> endpoint_lanes;
  std::vector<std::size_t> endpoint_destinations;
  std::vector<double> endpoint_values;
  std::vector<double> mapped_times;
  std::vector<double> terminal_values;
  std::vector<double> root_values;
  std::vector<double> state_values;
  std::vector<double> term_values;
  std::vector<double> trigger_weights;
};

inline void exact_interval_trigger_weights(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    ExactIntervalLaneWorkspace *workspace) {
  const auto state_count = plan.trigger_state_table.states.size();
  workspace->trigger_weights.resize(state_count * lanes.size);
  for (std::size_t state = 0U; state < state_count; ++state) {
    const auto &compiled = plan.trigger_state_table.states[state];
    auto *destination = workspace->trigger_weights.data() + state * lanes.size;
    if (compiled.weight_terms.empty()) {
      std::fill_n(destination, lanes.size, compiled.fixed_weight);
      continue;
    }
    exact_compiled_trigger_state_weights_lanes(
        plan, lanes, lanes.size, compiled, &workspace->root_values);
    std::copy_n(workspace->root_values.data(), lanes.size, destination);
  }
}

inline void exact_response_survival_endpoints(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    ExactStepLaneWorkspace *exact_workspace,
    ExactIntervalLaneWorkspace *workspace) {
  const auto request_count = workspace->endpoint_times.size();
  if (request_count == 0U) {
    workspace->endpoint_values.clear();
    return;
  }
  if (exact_has_single_unweighted_trigger_state(plan)) {
    const auto &compiled = plan.trigger_state_table.states.front();
    exact_workspace->bind_initial_sources(
        lanes, exact_compiled_trigger_shared_started(plan, compiled));
    workspace->endpoint_values.resize(request_count);
    for (std::size_t begin = 0U;
         begin < request_count;
         begin += kExactLaneTileSize) {
      const auto count = std::min(kExactLaneTileSize, request_count - begin);
      auto &frame = exact_workspace->prepare_mapped_times(
          workspace->endpoint_times.data() + begin,
          workspace->endpoint_lanes.data() + begin,
          count);
      evaluate_compiled_lane_root(
          plan,
          plan.finite_response_survival_root_id,
          &exact_workspace->executor,
          &frame,
          count,
          &workspace->root_values);
      for (std::size_t position = 0U; position < count; ++position) {
        const double value = workspace->root_values[position];
        workspace->endpoint_values[begin + position] =
            std::isfinite(value) && value > 0.0 ? value : 0.0;
      }
    }
    return;
  }
  workspace->endpoint_values.assign(request_count, 0.0);
  exact_interval_trigger_weights(plan, lanes, workspace);
  for (std::size_t state = 0U;
       state < plan.trigger_state_table.states.size();
       ++state) {
    const auto &compiled = plan.trigger_state_table.states[state];
    exact_workspace->bind_initial_sources(
        lanes, exact_compiled_trigger_shared_started(plan, compiled));
    for (std::size_t begin = 0U;
         begin < request_count;
         begin += kExactLaneTileSize) {
      const auto count = std::min(kExactLaneTileSize, request_count - begin);
      auto &frame = exact_workspace->prepare_mapped_times(
          workspace->endpoint_times.data() + begin,
          workspace->endpoint_lanes.data() + begin,
          count);
      evaluate_compiled_lane_root(
          plan,
          plan.finite_response_survival_root_id,
          &exact_workspace->executor,
          &frame,
          count,
          &workspace->root_values);
      for (std::size_t position = 0U; position < count; ++position) {
        const auto request = begin + position;
        const auto source = workspace->endpoint_lanes[request];
        const double value = workspace->root_values[position];
        const double weight = workspace->trigger_weights[
            state * lanes.size + source];
        if (weight > 0.0 && std::isfinite(value) && value > 0.0) {
          workspace->endpoint_values[request] += weight * value;
        }
      }
    }
  }
}

inline void exact_add_response_endpoint(
    const ObservationLaneBatchView lanes,
    const std::size_t lane,
    const std::size_t destination,
    const double time,
    ExactIntervalLaneWorkspace *workspace) {
  for (std::size_t request = workspace->endpoint_times.size();
       request > 0U;
       --request) {
    const auto previous_lane = workspace->endpoint_lanes[request - 1U];
    if (lanes.row_maps[previous_lane] != lanes.row_maps[lane] ||
        lanes.row_offsets[previous_lane] != lanes.row_offsets[lane]) {
      break;
    }
    if (workspace->endpoint_times[request - 1U] == time) {
      workspace->endpoint_destinations[destination] = request - 1U;
      return;
    }
  }
  workspace->endpoint_destinations[destination] =
      workspace->endpoint_times.size();
  workspace->endpoint_times.push_back(time);
  workspace->endpoint_lanes.push_back(lane);
}

inline void exact_finite_response_density_mapped_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    const double *times,
    const std::size_t *source_lanes,
    const std::size_t count,
    ExactStepLaneWorkspace *exact_workspace,
    ExactIntervalLaneWorkspace *workspace,
    double *out) {
  if (exact_has_single_unweighted_trigger_state(plan)) {
    const auto &compiled = plan.trigger_state_table.states.front();
    exact_workspace->bind_initial_sources(
        lanes, exact_compiled_trigger_shared_started(plan, compiled));
    auto &frame = exact_workspace->prepare_mapped_times(
        times, source_lanes, count);
    evaluate_compiled_lane_root(
        plan,
        plan.finite_response_density_root_id,
        &exact_workspace->executor,
        &frame,
        count,
        &workspace->root_values);
    for (std::size_t position = 0U; position < count; ++position) {
      const double value = workspace->root_values[position];
      out[position] =
          std::isfinite(value) && value > 0.0 ? value : 0.0;
    }
    return;
  }
  std::fill_n(out, count, 0.0);
  for (std::size_t state = 0U;
       state < plan.trigger_state_table.states.size();
       ++state) {
    const auto &compiled = plan.trigger_state_table.states[state];
    exact_workspace->bind_initial_sources(
        lanes, exact_compiled_trigger_shared_started(plan, compiled));
    auto &frame = exact_workspace->prepare_mapped_times(
        times, source_lanes, count);
    evaluate_compiled_lane_root(
        plan,
        plan.finite_response_density_root_id,
        &exact_workspace->executor,
        &frame,
        count,
        &workspace->root_values);
    for (std::size_t position = 0U; position < count; ++position) {
      const auto source = source_lanes[position];
      const double value = workspace->root_values[position];
      const double weight = workspace->trigger_weights[
          state * lanes.size + source];
      if (weight > 0.0 && std::isfinite(value) && value > 0.0) {
        out[position] += weight * value;
      }
    }
  }
}

inline void exact_finite_response_probability_between_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    const double *lower,
    const double *upper,
    ExactStepLaneWorkspace *exact_workspace,
    ExactIntervalLaneWorkspace *workspace,
    std::vector<double> *out) {
  workspace->lower.resize(lanes.size);
  workspace->upper.resize(lanes.size);
  bool has_infinite_upper = false;
  for (std::size_t lane = 0U; lane < lanes.size; ++lane) {
    const double lo = std::isfinite(lower[lane])
                          ? std::max(0.0, lower[lane])
                          : 0.0;
    if (upper[lane] == std::numeric_limits<double>::infinity() &&
        upper[lane] > lower[lane]) {
      has_infinite_upper = true;
      workspace->lower[lane] = 0.0;
      workspace->upper[lane] = 1.0;
    } else if (std::isfinite(upper[lane]) && upper[lane] > lo) {
      workspace->lower[lane] = lo;
      workspace->upper[lane] = upper[lane];
    } else {
      workspace->lower[lane] = 0.0;
      workspace->upper[lane] = 0.0;
    }
  }
  if (!exact_has_single_unweighted_trigger_state(plan)) {
    exact_interval_trigger_weights(plan, lanes, workspace);
  }
  const auto evaluate = [&](const double *integration_times,
                            const std::size_t *source_lanes,
                            const std::size_t count,
                            double *values) {
    const double *evaluation_times = integration_times;
    if (has_infinite_upper) {
      workspace->mapped_times.resize(count);
      for (std::size_t position = 0U; position < count; ++position) {
        const auto source = source_lanes[position];
        if (upper[source] == std::numeric_limits<double>::infinity()) {
          const double unit_time = integration_times[position];
          const double remaining = 1.0 - unit_time;
          const double origin = std::isfinite(lower[source])
                                    ? std::max(0.0, lower[source])
                                    : 0.0;
          workspace->mapped_times[position] =
              origin + unit_time / remaining;
        } else {
          workspace->mapped_times[position] = integration_times[position];
        }
      }
      evaluation_times = workspace->mapped_times.data();
    }
    exact_finite_response_density_mapped_lanes(
        plan,
        lanes,
        evaluation_times,
        source_lanes,
        count,
        exact_workspace,
        workspace,
        values);
    if (has_infinite_upper) {
      for (std::size_t position = 0U; position < count; ++position) {
        const auto source = source_lanes[position];
        if (upper[source] == std::numeric_limits<double>::infinity()) {
          const double remaining = 1.0 - integration_times[position];
          values[position] /= remaining * remaining;
        }
        if (!std::isfinite(values[position]) || values[position] <= 0.0) {
          values[position] = 0.0;
        }
      }
    }
  };
  adaptive_integrate_lane_batch(
      lanes.size,
      workspace->lower.data(),
      workspace->upper.data(),
      evaluate,
      &workspace->adaptive,
      out);
  for (double &value : *out) {
    value = clamp_probability(value);
  }
}

inline void exact_any_response_probability_between_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    const double *lower,
    const double *upper,
    ExactStepLaneWorkspace *exact_workspace,
    ExactIntervalLaneWorkspace *workspace,
    std::vector<double> *out) {
  const bool direct_terminal =
      plan.no_response.direct_leaf_failure_product &&
      !plan.no_response.leaf_indices.empty() &&
      plan.finite_response_survival_root_id != semantic::kInvalidIndex;
  if (!direct_terminal) {
    exact_finite_response_probability_between_lanes(
        plan, lanes, lower, upper, exact_workspace, workspace, out);
    return;
  }
  out->assign(lanes.size, 0.0);
  workspace->endpoint_times.clear();
  workspace->endpoint_lanes.clear();
  workspace->endpoint_destinations.assign(2U * lanes.size,
                                          std::numeric_limits<std::size_t>::max());
  bool needs_terminal = false;
  for (std::size_t lane = 0U; lane < lanes.size; ++lane) {
    if (std::isnan(lower[lane]) || std::isnan(upper[lane]) ||
        !(upper[lane] > lower[lane])) {
      continue;
    }
    const bool terminal =
        upper[lane] == std::numeric_limits<double>::infinity();
    needs_terminal = needs_terminal || terminal;
    if (lower[lane] != -std::numeric_limits<double>::infinity()) {
      exact_add_response_endpoint(
          lanes,
          lane,
          2U * lane,
          std::max(0.0, lower[lane]),
          workspace);
    }
    if (std::isfinite(upper[lane])) {
      exact_add_response_endpoint(
          lanes, lane, 2U * lane + 1U, upper[lane], workspace);
    }
  }
  exact_response_survival_endpoints(
      plan, lanes, exact_workspace, workspace);
  if (needs_terminal) {
    exact_terminal_no_response_probability_lanes(
        plan, lanes, exact_workspace, &workspace->terminal_values);
  }
  const auto absent = std::numeric_limits<std::size_t>::max();
  for (std::size_t lane = 0U; lane < lanes.size; ++lane) {
    if (std::isnan(lower[lane]) || std::isnan(upper[lane]) ||
        !(upper[lane] > lower[lane])) {
      continue;
    }
    const bool terminal =
        upper[lane] == std::numeric_limits<double>::infinity();
    const auto upper_index = workspace->endpoint_destinations[2U * lane + 1U];
    if (!terminal && upper_index == absent) {
      continue;
    }
    const auto lower_index = workspace->endpoint_destinations[2U * lane];
    const double lower_survival =
        lower_index == absent ? 1.0 : workspace->endpoint_values[lower_index];
    const double upper_survival =
        terminal
            ? workspace->terminal_values[lane]
            : workspace->endpoint_values[upper_index];
    (*out)[lane] = clamp_probability(
        std::max(0.0, lower_survival - upper_survival));
  }
}

inline void exact_weighted_outcome_density_mapped_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    const std::vector<ExactOutcomeTerm> &terms,
    const double *times,
    const std::size_t *source_lanes,
    const std::size_t count,
    ExactStepLaneWorkspace *exact_workspace,
    ExactIntervalLaneWorkspace *workspace,
    double *out) {
  if (exact_has_single_unweighted_trigger_state(plan)) {
    const auto &compiled = plan.trigger_state_table.states.front();
    exact_workspace->bind_initial_sources(
        lanes, exact_compiled_trigger_shared_started(plan, compiled));
    auto &frame = exact_workspace->prepare_mapped_times(
        times, source_lanes, count);
    std::fill_n(out, count, 0.0);
    for (const auto &term : terms) {
      evaluate_exact_step_distribution_prepared_lanes(
          plan,
          &frame,
          count,
          term.target,
          false,
          exact_workspace,
          &workspace->term_values);
      for (std::size_t position = 0U; position < count; ++position) {
        const double value = workspace->term_values[position];
        if (std::isfinite(value) && value > 0.0) {
          out[position] += term.weight * value;
        }
      }
    }
    return;
  }
  std::fill_n(out, count, 0.0);
  for (std::size_t state = 0U;
       state < plan.trigger_state_table.states.size();
       ++state) {
    const auto &compiled = plan.trigger_state_table.states[state];
    exact_workspace->bind_initial_sources(
        lanes, exact_compiled_trigger_shared_started(plan, compiled));
    auto &frame = exact_workspace->prepare_mapped_times(
        times, source_lanes, count);
    workspace->state_values.assign(count, 0.0);
    for (const auto &term : terms) {
      evaluate_exact_step_distribution_prepared_lanes(
          plan,
          &frame,
          count,
          term.target,
          false,
          exact_workspace,
          &workspace->term_values);
      for (std::size_t position = 0U; position < count; ++position) {
        const double value = workspace->term_values[position];
        if (std::isfinite(value) && value > 0.0) {
          workspace->state_values[position] += term.weight * value;
        }
      }
    }
    for (std::size_t position = 0U; position < count; ++position) {
      const auto source = source_lanes[position];
      const double weight = workspace->trigger_weights[
          state * lanes.size + source];
      if (weight > 0.0) {
        out[position] += weight * workspace->state_values[position];
      }
    }
  }
}

inline void exact_integrated_outcome_probability_between_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    const std::vector<ExactOutcomeTerm> &terms,
    const double *lower,
    const double *upper,
    const bool transform_infinite_tail,
    ExactStepLaneWorkspace *exact_workspace,
    ExactIntervalLaneWorkspace *workspace,
    std::vector<double> *out) {
  if (terms.empty()) {
    out->assign(lanes.size, 0.0);
    return;
  }
  workspace->lower.resize(lanes.size);
  workspace->upper.resize(lanes.size);
  bool has_infinite_upper = false;
  for (std::size_t lane = 0U; lane < lanes.size; ++lane) {
    const double lo = std::isfinite(lower[lane]) ? std::max(0.0, lower[lane])
                                                 : 0.0;
    if (transform_infinite_tail &&
        upper[lane] == std::numeric_limits<double>::infinity() &&
        upper[lane] > lower[lane]) {
      has_infinite_upper = true;
      workspace->lower[lane] = 0.0;
      workspace->upper[lane] = 1.0;
    } else {
      const double hi = std::isfinite(upper[lane])
                            ? upper[lane]
                            : kExactResponseTimeLimit;
      workspace->lower[lane] = lo;
      workspace->upper[lane] = hi > lo ? hi : lo;
    }
  }
  if (!exact_has_single_unweighted_trigger_state(plan)) {
    exact_interval_trigger_weights(plan, lanes, workspace);
  }
  const auto evaluate = [&](const double *times,
                            const std::size_t *source_lanes,
                            const std::size_t count,
                            double *values) {
    const double *evaluation_times = times;
    if (has_infinite_upper) {
      workspace->mapped_times.resize(count);
      for (std::size_t position = 0U; position < count; ++position) {
        const auto source = source_lanes[position];
        if (upper[source] == std::numeric_limits<double>::infinity()) {
          const double remaining = 1.0 - times[position];
          const double origin = std::isfinite(lower[source])
                                    ? std::max(0.0, lower[source])
                                    : 0.0;
          workspace->mapped_times[position] =
              origin + times[position] / remaining;
        } else {
          workspace->mapped_times[position] = times[position];
        }
      }
      evaluation_times = workspace->mapped_times.data();
    }
    exact_weighted_outcome_density_mapped_lanes(
        plan,
        lanes,
        terms,
        evaluation_times,
        source_lanes,
        count,
        exact_workspace,
        workspace,
        values);
    if (has_infinite_upper) {
      for (std::size_t position = 0U; position < count; ++position) {
        const auto source = source_lanes[position];
        if (upper[source] == std::numeric_limits<double>::infinity()) {
          const double remaining = 1.0 - times[position];
          values[position] /= remaining * remaining;
        }
        if (!std::isfinite(values[position]) || values[position] <= 0.0) {
          values[position] = 0.0;
        }
      }
    }
  };
  adaptive_integrate_lane_batch(
      lanes.size,
      workspace->lower.data(),
      workspace->upper.data(),
      evaluate,
      &workspace->adaptive,
      out);
  for (double &value : *out) {
    value = clamp_probability(value);
  }
}

inline void exact_response_probability_between_lanes(
    const ExactVariantPlan &plan,
    const ExactResponseMeasure measure,
    const bool complete_outcome_partition,
    const std::vector<ExactOutcomeTerm> &terms,
    const ObservationLaneBatchView lanes,
    const double *lower,
    const double *upper,
    ExactStepLaneWorkspace *exact_workspace,
    ExactIntervalLaneWorkspace *workspace,
    std::vector<double> *out) {
  if (measure == ExactResponseMeasure::ObservableResponses &&
      complete_outcome_partition) {
    exact_any_response_probability_between_lanes(
        plan, lanes, lower, upper, exact_workspace, workspace, out);
    return;
  }
  exact_integrated_outcome_probability_between_lanes(
      plan,
      lanes,
      terms,
      lower,
      upper,
      measure == ExactResponseMeasure::ObservableResponses,
      exact_workspace,
      workspace,
      out);
}

} // namespace detail
} // namespace accumulatr::eval
