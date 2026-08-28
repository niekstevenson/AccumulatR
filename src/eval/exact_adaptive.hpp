#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace accumulatr::eval {
namespace detail {

inline constexpr std::size_t kAdaptiveNodeBatchSize = 512U;
inline constexpr std::size_t kKronrod15NodeCount = 15U;
inline constexpr double kAdaptiveAbsoluteTolerance = 1e-8;
inline constexpr double kAdaptiveRelativeTolerance = 1e-6;
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

inline bool adaptive_lane_converged(
    const AdaptiveLaneState &state) noexcept {
  return state.error <=
         std::max(
             kAdaptiveRelativeTolerance * std::fabs(state.value),
             kAdaptiveAbsoluteTolerance);
}

inline void adaptive_map_kronrod15_panel(
    const AdaptiveLanePanel &panel,
    const std::size_t panel_index,
    AdaptiveLaneWorkspace *workspace) noexcept {
  static constexpr std::array<double, 7> nodes{
      -0.99145537112081263920685469752598,
      -0.94910791234275852452618968404809,
      -0.86486442335976907278971278864098,
      -0.74153118559939443986386477328079,
      -0.58608723546769113029414483825842,
      -0.40584515137739716690660641207696,
      -0.20778495500789846760068940377324};
  const double center = 0.5 * (panel.lower + panel.upper);
  const double delta = 0.5 * (panel.upper - panel.lower);
  const auto offset = panel_index * kKronrod15NodeCount;
  for (std::size_t node = 0U; node < 7U; ++node) {
    const double distance = delta * nodes[node];
    workspace->mapped_times[offset + 2U * node] = center + distance;
    workspace->mapped_times[offset + 2U * node + 1U] =
        center - distance;
    workspace->mapped_lanes[offset + 2U * node] = panel.lane;
    workspace->mapped_lanes[offset + 2U * node + 1U] = panel.lane;
  }
  workspace->mapped_times[offset + 14U] = center;
  workspace->mapped_lanes[offset + 14U] = panel.lane;
}

inline void adaptive_finish_kronrod15_panel(
    AdaptiveLanePanel *panel,
    const double *values) noexcept {
  static constexpr std::array<double, 8> kronrod_weights{
      0.022935322010529224963732008059913,
      0.063092092629978553290700663189093,
      0.10479001032225018383987632254189,
      0.14065325971552591874518959051021,
      0.16900472663926790282619833697952,
      0.19035057806478540991325640242055,
      0.20443294007529889241416199923466,
      0.20948214108472782801299917489173};
  static constexpr std::array<double, 4> gauss_weights{
      0.12948496616886969327061143267787,
      0.27970539148927666790146777142378,
      0.38183005050511894495036977548818,
      0.41795918367346938775510204081633};

  double kronrod = values[14U] * kronrod_weights[7U];
  double gauss = values[14U] * gauss_weights[3U];
  for (std::size_t node = 0U; node < 7U; ++node) {
    const double pair = values[2U * node] + values[2U * node + 1U];
    kronrod += pair * kronrod_weights[node];
    if ((node & 1U) != 0U) {
      gauss += pair * gauss_weights[node / 2U];
    }
  }
  const double scale = std::fabs(0.5 * (panel->upper - panel->lower));
  kronrod *= scale;
  panel->value = kronrod;
  panel->error = std::fabs(kronrod - gauss * scale);
}

template <typename Evaluate>
inline void adaptive_evaluate_kronrod15_panels(
    Evaluate &evaluate,
    AdaptiveLaneWorkspace *workspace) {
  const auto node_count =
      workspace->panels.size() * kKronrod15NodeCount;
  workspace->mapped_times.resize(node_count);
  workspace->mapped_lanes.resize(node_count);
  workspace->values.resize(node_count);
  for (std::size_t panel = 0U;
       panel < workspace->panels.size();
       ++panel) {
    adaptive_map_kronrod15_panel(
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
    adaptive_finish_kronrod15_panel(
        &workspace->panels[panel],
        workspace->values.data() + panel * kKronrod15NodeCount);
  }
}

inline double adaptive_drain_lane(AdaptiveLaneState *state) {
  double value = 0.0;
  while (!state->heap.empty()) {
    value += adaptive_heap_pop(state).value;
  }
  return value;
}

template <typename Evaluate>
inline void adaptive_integrate_lane_batch(
    const std::size_t lane_count,
    const double *lower,
    const double *upper,
    Evaluate &&evaluate,
    AdaptiveLaneWorkspace *workspace,
    std::vector<double> *out) {
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
      workspace->panels.push_back(
          AdaptiveLanePanel{lane, lower[lane], upper[lane], 0.0, 0.0});
    }
  }

  adaptive_evaluate_kronrod15_panels(evaluate, workspace);
  for (const auto &panel : workspace->panels) {
    auto &state = workspace->lanes[panel.lane];
    state.value = panel.value;
    state.error = panel.error;
    state.evaluations = kKronrod15NodeCount;
    adaptive_heap_push(&state, panel);
    const bool at_limit =
        state.evaluations >= kAdaptiveMaximumEvaluations;
    if (!adaptive_lane_converged(state) && !at_limit) {
      workspace->active.push_back(panel.lane);
    }
  }

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

    adaptive_evaluate_kronrod15_panels(evaluate, workspace);
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
      state.evaluations += 2U * kKronrod15NodeCount;
      const bool at_limit =
          state.evaluations >= kAdaptiveMaximumEvaluations;
      if (!adaptive_lane_converged(state) &&
          !at_limit &&
          std::isfinite(state.value)) {
        workspace->next_active.push_back(lane);
      }
    }
    workspace->active.swap(workspace->next_active);
  }

  for (std::size_t lane = 0U; lane < lane_count; ++lane) {
    (*out)[lane] = adaptive_drain_lane(&workspace->lanes[lane]);
  }
}

} // namespace detail
} // namespace accumulatr::eval
