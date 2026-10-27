#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace accumulatr::eval {
namespace detail {

inline constexpr std::size_t kAdaptiveNodeBatchSize = 512U;
inline constexpr double kAdaptiveAbsoluteTolerance = 1e-12;
inline constexpr double kAdaptiveRelativeTolerance = 1e-3;
inline constexpr std::size_t kAdaptiveMaximumEvaluations = 6000U;

struct AdaptiveLanePanel {
  std::size_t lane{0U};
  double lower{0.0};
  double upper{0.0};
  double value{0.0};
  double error{0.0};
};

struct AdaptiveLaneState {
  std::vector<AdaptiveLanePanel> heap;
  double value{0.0};
  double error{0.0};
  std::size_t evaluations{0U};
};

struct AdaptiveLaneWorkspace {
  std::vector<AdaptiveLaneState> lanes;
  std::vector<std::size_t> active;
  std::vector<std::size_t> next_active;
  std::vector<AdaptiveLanePanel> parents;
  std::vector<AdaptiveLanePanel> panels;
  std::vector<double> mapped_times;
  std::vector<std::size_t> mapped_lanes;
  std::vector<double> values;
};

inline bool adaptive_panel_less(
    const AdaptiveLanePanel &lhs,
    const AdaptiveLanePanel &rhs) noexcept {
  return lhs.error < rhs.error;
}

inline void adaptive_heap_push(
    AdaptiveLaneState *state,
    const AdaptiveLanePanel &panel) {
  state->heap.push_back(panel);
  std::push_heap(
      state->heap.begin(), state->heap.end(), adaptive_panel_less);
}

inline AdaptiveLanePanel adaptive_heap_pop(AdaptiveLaneState *state) {
  std::pop_heap(
      state->heap.begin(), state->heap.end(), adaptive_panel_less);
  const auto panel = state->heap.back();
  state->heap.pop_back();
  return panel;
}

inline bool adaptive_integral_converged(
    const double value,
    const double error,
    const double absolute_tolerance,
    const double relative_tolerance) noexcept {
  return error <= std::max(relative_tolerance * std::fabs(value), absolute_tolerance);
}

template <std::size_t Order>
inline void adaptive_map_panel(
    const AdaptiveLanePanel &panel,
    const std::size_t panel_index,
    AdaptiveLaneWorkspace *workspace) noexcept {
  constexpr auto half = (Order - 1U) / 2U;
  const auto &nodes = boost::math::quadrature::gauss_kronrod<double, Order>::abscissa();
  const double center = 0.5 * (panel.lower + panel.upper);
  const double delta = 0.5 * (panel.upper - panel.lower);
  const auto offset = panel_index * Order;
  for (std::size_t node = 0U; node < half; ++node) {
    const double distance = delta * nodes[half - node];
    workspace->mapped_times[offset + 2U * node] = center - distance;
    workspace->mapped_times[offset + 2U * node + 1U] = center + distance;
    workspace->mapped_lanes[offset + 2U * node] = panel.lane;
    workspace->mapped_lanes[offset + 2U * node + 1U] = panel.lane;
  }
  workspace->mapped_times[offset + Order - 1U] = center;
  workspace->mapped_lanes[offset + Order - 1U] = panel.lane;
}

template <std::size_t Order>
inline void adaptive_finish_panel(AdaptiveLanePanel *panel, const double *values) noexcept {
  constexpr auto half = (Order - 1U) / 2U;
  const auto &kronrod_weights = boost::math::quadrature::gauss_kronrod<double, Order>::weights();
  const auto &gauss_weights = boost::math::quadrature::gauss<double, half>::weights();
  double kronrod = values[Order - 1U] * kronrod_weights[0];
  double gauss = half % 2U ? values[Order - 1U] * gauss_weights[0] : 0.0;
  for (std::size_t node = 0U; node < half; ++node) {
    const double pair = values[2U * node] + values[2U * node + 1U];
    kronrod += pair * kronrod_weights[half - node];
    if ((node & 1U) != 0U) gauss += pair * gauss_weights[(half - node) / 2U];
  }
  const double scale = std::fabs(0.5 * (panel->upper - panel->lower));
  kronrod *= scale;
  panel->value = kronrod;
  panel->error = std::fabs(kronrod - gauss * scale);
}

template <std::size_t Order, typename Evaluate>
inline void adaptive_evaluate_panels(
    Evaluate &evaluate,
    AdaptiveLaneWorkspace *workspace) {
  const auto node_count =
      workspace->panels.size() * Order;
  workspace->mapped_times.resize(node_count);
  workspace->mapped_lanes.resize(node_count);
  workspace->values.resize(node_count);
  for (std::size_t panel = 0U;
       panel < workspace->panels.size();
       ++panel) {
    adaptive_map_panel<Order>(
        workspace->panels[panel], panel, workspace);
  }
  for (std::size_t begin = 0U;
       begin < node_count;
       begin += kAdaptiveNodeBatchSize) {
    const auto count =
        std::min(kAdaptiveNodeBatchSize, node_count - begin);
    evaluate(
        workspace->mapped_times.data() + begin,
        workspace->mapped_lanes.data() + begin,
        count,
        workspace->values.data() + begin);
  }
  for (std::size_t panel = 0U;
       panel < workspace->panels.size();
       ++panel) {
    adaptive_finish_panel<Order>(
        &workspace->panels[panel],
        workspace->values.data() + panel * Order);
  }
}

struct AdaptiveNoRefinement {
  void operator()(AdaptiveLaneWorkspace *) const noexcept {}
};

template <std::size_t Order = 15U, std::size_t InitialPanels = 1U,
          typename Evaluate, typename Refine = AdaptiveNoRefinement>
inline void adaptive_integrate_lane_batch(
    const std::size_t lane_count,
    const double *lower,
    const double *upper,
    Evaluate &&evaluate,
    AdaptiveLaneWorkspace *workspace,
    std::vector<double> *out,
    const double absolute_tolerance = kAdaptiveAbsoluteTolerance,
    const double relative_tolerance = kAdaptiveRelativeTolerance,
    const double *absolute_scales = nullptr,
    const double *tolerance_divisors = nullptr,
    Refine refine = {}) {
  const auto converged = [&](std::size_t lane) {
    const auto &state = workspace->lanes[lane];
    const double absolute = absolute_scales == nullptr ? absolute_tolerance
        : absolute_tolerance / absolute_scales[lane];
    const double divisor = tolerance_divisors == nullptr ? 1.0 : tolerance_divisors[lane];
    return adaptive_integral_converged(state.value, state.error,
        absolute / divisor, relative_tolerance / divisor);
  };
  if (workspace->lanes.size() < lane_count) {
    workspace->lanes.resize(lane_count);
  }
  workspace->active.clear();
  workspace->panels.clear();
  out->assign(lane_count, 0.0);

  for (std::size_t lane = 0U; lane < lane_count; ++lane) {
    auto &state = workspace->lanes[lane];
    state.heap.clear();
    state.value = 0.0;
    state.error = 0.0;
    state.evaluations = 0U;
    if (std::isfinite(lower[lane]) &&
        std::isfinite(upper[lane]) &&
        upper[lane] > lower[lane]) {
      const double width = (upper[lane] - lower[lane]) / InitialPanels;
      for (std::size_t i = 0; i < InitialPanels; ++i)
        workspace->panels.push_back(AdaptiveLanePanel{lane,
            lower[lane] + i * width, i + 1U == InitialPanels ? upper[lane]
                : lower[lane] + (i + 1U) * width, 0.0, 0.0});
    }
  }

  adaptive_evaluate_panels<Order>(evaluate, workspace);
  refine(workspace);
  for (const auto &panel : workspace->panels) {
    auto &state = workspace->lanes[panel.lane];
    state.value += panel.value;
    state.error += panel.error;
    state.evaluations += Order;
    adaptive_heap_push(&state, panel);
  }
  for (std::size_t lane = 0; lane < lane_count; ++lane)
    if (!converged(lane)) workspace->active.push_back(lane);

  while (!workspace->active.empty()) {
    const auto active_count = workspace->active.size();
    workspace->parents.resize(active_count);
    workspace->panels.resize(2U * active_count);
    for (std::size_t position = 0U;
         position < active_count;
         ++position) {
      const auto lane = workspace->active[position];
      const auto parent = adaptive_heap_pop(&workspace->lanes[lane]);
      workspace->parents[position] = parent;
      const double width = 0.5 * (parent.upper - parent.lower);
      workspace->panels[2U * position] = AdaptiveLanePanel{
          lane,
          parent.lower + width,
          parent.upper,
          0.0,
          0.0};
      workspace->panels[2U * position + 1U] = AdaptiveLanePanel{
          lane,
          parent.lower,
          parent.upper - width,
          0.0,
          0.0};
    }

    adaptive_evaluate_panels<Order>(evaluate, workspace);
    refine(workspace);
    workspace->next_active.clear();
    for (std::size_t position = 0U;
         position < active_count;
         ++position) {
      const auto lane = workspace->active[position];
      auto &state = workspace->lanes[lane];
      const auto &parent = workspace->parents[position];
      const auto &right = workspace->panels[2U * position];
      const auto &left = workspace->panels[2U * position + 1U];
      adaptive_heap_push(&state, right);
      adaptive_heap_push(&state, left);
      state.value += right.value + left.value - parent.value;
      state.error += right.error + left.error - parent.error;
      state.evaluations += 2U * Order;
      const bool at_limit =
          state.evaluations >= kAdaptiveMaximumEvaluations;
      if (!converged(lane) &&
          !at_limit &&
          std::isfinite(state.value)) {
        workspace->next_active.push_back(lane);
      }
    }
    workspace->active.swap(workspace->next_active);
  }

  for (std::size_t lane = 0U; lane < lane_count; ++lane) {
    for (const auto &panel : workspace->lanes[lane].heap) (*out)[lane] += panel.value;
  }
}

} // namespace detail
} // namespace accumulatr::eval
