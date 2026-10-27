#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "exact_source_lane_eval.hpp"
#include "exact_adaptive.hpp"
#include "cumulative_lane.hpp"

namespace accumulatr::eval {
namespace detail {

struct CompiledLaneScratch {
  void ensure(const std::size_t count) {
    lanes.resize(count);
    positions.resize(count);
    parents.resize(count);
    times.resize(count);
    values.resize(count);
    child_values.resize(count);
    products.resize(count);
  }

  std::vector<semantic::Index> lanes;
  std::vector<std::size_t> positions;
  std::vector<std::size_t> parents;
  std::vector<double> times;
  std::vector<double> integral_lower;
  std::vector<double> integral_upper;
  std::vector<double> integral_origins;
  std::vector<double> density_values;
  std::vector<std::uint8_t> infinite_integrals;
  AdaptiveLaneWorkspace adaptive;
  std::vector<double> values;
  std::vector<double> child_values;
  std::vector<double> products;
  SourceLaneFill source_values;
  std::vector<double> integral_scales;
  std::vector<double> factor_scales;
  CumulativeLaneWorkspace cumulative;
};

struct CompiledLaneExecutor {
  explicit CompiledLaneExecutor(const CompiledMathProgram &program)
      : lanes(program) {}

  CompiledLaneScratch &scratch(const std::size_t depth) {
    while (scratch_layers.size() <= depth) {
      scratch_layers.push_back(std::make_unique<CompiledLaneScratch>());
    }
    return *scratch_layers[depth];
  }

  double terminal_relative_tolerance{1e-9};
  CompiledLaneWorkspace lanes;
  SourceLaneWorkspace sources;
  std::vector<std::unique_ptr<CompiledLaneScratch>> scratch_layers;
};

inline double compiled_lane_node_time(const CompiledMathNode &node,
                                      const CompiledLaneFrame &frame,
                                      const std::size_t lane) {
  return frame.time(node.time_id, lane);
}

inline double compiled_lane_source_node_time(
    const CompiledMathNode &node,
    const CompiledLaneFrame &frame,
    const std::size_t lane) {
  const double time = frame.time(node.time_id, lane);
  return node.aux_id == semantic::kInvalidIndex
             ? time
             : std::min(time, frame.time(node.aux_id, lane));
}

inline bool compiled_lane_outcome_used(
    const CompiledMathNode &node,
    const ExactVariantPlan &plan,
    const CompiledLaneFrame &frame,
    const std::size_t lane) {
  if (!frame.has_sequence_history) {
    return false;
  }
  const auto *used = frame.used_outcomes[lane];
  bool any_used = false;
  if (used != nullptr) {
    const auto offset = static_cast<std::size_t>(node.subject_id);
    const auto size = static_cast<std::size_t>(node.aux_id);
    for (std::size_t i = 0; i < size; ++i) {
      const auto outcome = plan.compiled_outcome_gate_indices[offset + i];
      if ((*used)[static_cast<std::size_t>(outcome)] != 0U) {
        any_used = true;
        break;
      }
    }
  }
  return any_used;
}

struct CompiledLaneUpperBound {
  bool found{false};
  double time{std::numeric_limits<double>::infinity()};
  double normalizer{0.0};
};

inline CompiledLaneUpperBound compiled_lane_expr_upper_bound(
    const CompiledMathNode &node,
    CompiledLaneExecutor *executor,
    const CompiledLaneFrame &frame,
    const std::size_t lane) {
  ExactTimedExprUpperBound sequence_upper;
  const auto source_lane = executor->lanes.source_lane(frame, lane);
  if (executor->lanes.source_state().expr_upper_bound(
          source_lane, node.subject_id, &sequence_upper)) {
    return CompiledLaneUpperBound{
        true, sequence_upper.time, sequence_upper.normalizer};
  }
  return {};
}

inline void evaluate_compiled_lane_root(
    const ExactVariantPlan &plan, semantic::Index root_id,
    CompiledLaneExecutor *executor, CompiledLaneFrame *frame,
    std::size_t lane_count, std::vector<double> *out,
    std::size_t integral_depth);

inline void evaluate_compiled_lane_schedule(
    const ExactVariantPlan &plan,
    const CompiledMathIndexSpan schedule,
    const semantic::Index result_node_id,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *frame,
    std::size_t lane_count,
    std::vector<double> *out,
    std::size_t integral_depth = 0U);

inline void evaluate_compiled_lane_integral_node(
    const ExactVariantPlan &plan,
    const CompiledMathNode &node,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *parent,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::vector<double> *out,
    const std::size_t integral_depth,
    const double *absolute_scales = nullptr);

inline void evaluate_compiled_lane_execution(
    const ExactVariantPlan &plan,
    const CompiledMathRoot &root,
    const CompiledMathExecutionPlan &execution,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *frame,
    std::size_t lane_count,
    std::vector<double> *out,
    std::size_t integral_depth);

inline void evaluate_source_product_ops_lanes(
    const ExactVariantPlan &plan,
    const CompiledMathSourceProductOps ops,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *frame,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::vector<double> *out,
    const std::size_t scratch_depth) {
  const auto &program = plan.compiled_math;
  auto &scratch = executor->scratch(scratch_depth);
  if (lane_count == 0U) {
    out->clear();
    return;
  }
  if (ops.empty()) {
    out->assign(lane_count, 1.0);
    return;
  }
  for (semantic::Index op_index = 0;
       op_index < ops.size;
       ++op_index) {
    const auto &op = program.source_product_ops[
        static_cast<std::size_t>(ops.offset + op_index)];
    if (op.value_channel_mask == 0U) {
      out->assign(lane_count, 0.0);
      return;
    }
    const double *times = nullptr;
    if (lanes == nullptr && op.time_cap_id == semantic::kInvalidIndex) {
      times = frame->time_values.data() +
              static_cast<std::size_t>(op.time_id) * frame->stride;
    } else {
      scratch.times.resize(lane_count);
      for (std::size_t lane_index = 0; lane_index < lane_count; ++lane_index) {
        const auto frame_lane = compiled_frame_lane(lanes, lane_index);
        const auto lane = static_cast<std::size_t>(frame_lane);
        const double time = frame->time(op.time_id, lane);
        scratch.times[lane_index] =
            op.time_cap_id == semantic::kInvalidIndex
                ? time
                : std::min(time, frame->time(op.time_cap_id, lane));
      }
      times = scratch.times.data();
    }
    auto &source_values = evaluate_source_program_cached_lanes(
        program,
        op.source_product_program_id,
        &executor->lanes,
        frame,
        lanes,
        times,
        lane_count,
        op.fill_channel_mask,
        op.cache_result,
        &executor->sources,
        &scratch.source_values);
    if (op_index == 0) {
      out->swap(source_lane_fill_storage(
          &source_values, op.value_channel_mask));
      continue;
    }
    const double *factors = source_lane_fill_values(
        source_values, op.value_channel_mask);
    const bool final_factor = op_index + 1 == ops.size;
    if (final_factor && ops.can_overflow) {
      for (std::size_t lane_index = 0;
           lane_index < lane_count;
           ++lane_index) {
        double &product = (*out)[lane_index];
        product *= factors[lane_index];
        if (!std::isfinite(product)) {
          product = 0.0;
        }
      }
    } else {
      for (std::size_t lane_index = 0;
           lane_index < lane_count;
           ++lane_index) {
        (*out)[lane_index] *= factors[lane_index];
      }
    }
  }
}

inline void evaluate_source_product_terms(
    const ExactVariantPlan &plan,
    const CompiledMathExecutionPlan &execution,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *frame,
    const std::size_t lane_count,
    std::vector<double> *out,
    const std::size_t integral_depth) {
  const auto &program = plan.compiled_math;
  auto &scratch = executor->scratch(integral_depth * 4U + 1U);
  scratch.ensure(lane_count);
  const bool has_common_sources = !execution.source_product_ops.empty();
  if (has_common_sources) {
    evaluate_source_product_ops_lanes(
        plan,
        execution.source_product_ops,
        executor,
        frame,
        nullptr,
        lane_count,
        &scratch.child_values,
        integral_depth * 4U + 2U);
  }
  out->assign(lane_count, 0.0);
  for (semantic::Index term_index = 0;
       term_index < execution.source_product_terms.size;
       ++term_index) {
    const auto &term = program.source_product_terms[
        static_cast<std::size_t>(execution.source_product_terms.offset +
                                 term_index)];
    std::size_t active_count = lane_count;
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      scratch.positions[lane] = lane;
    }
    for (semantic::Index i = 0;
         i < term.time_gate_nodes.size && active_count > 0U;
         ++i) {
      const auto node_id = program.time_gate_nodes[
          static_cast<std::size_t>(term.time_gate_nodes.offset + i)];
      const auto &gate = program.nodes[static_cast<std::size_t>(node_id)];
      std::size_t next_count = 0U;
      for (std::size_t j = 0; j < active_count; ++j) {
        const auto lane = scratch.positions[j];
        const double current = frame->time(gate.time_id, lane);
        const double bound = frame->time(gate.aux_id, lane);
        const bool open = gate.kind == CompiledMathNodeKind::StrictTimeGate
                              ? current > bound
                              : current >= bound;
        if (open) {
          scratch.positions[next_count++] = lane;
        }
      }
      active_count = next_count;
    }
    for (std::size_t j = 0; j < active_count; ++j) {
      scratch.lanes[j] =
          static_cast<semantic::Index>(scratch.positions[j]);
    }
    if (active_count == 0U) {
      continue;
    }

    evaluate_source_product_ops_lanes(
        plan, term.source_product_ops, executor, frame, scratch.lanes.data(),
        active_count, &scratch.values, integral_depth * 4U + 2U);
    std::size_t nonzero_count = 0U;
    for (std::size_t j = 0; j < active_count; ++j) {
      const auto lane = scratch.positions[j];
      const double product = term.sign * scratch.values[j];
      if (product != 0.0 && std::isfinite(product) &&
          (!has_common_sources || scratch.child_values[lane] != 0.0)) {
        scratch.products[nonzero_count] = product;
        scratch.lanes[nonzero_count] = scratch.lanes[j];
        scratch.positions[nonzero_count++] = lane;
      }
    }
    active_count = nonzero_count;

    for (semantic::Index factor_index = 0;
         factor_index < term.integral_factor_nodes.size && active_count > 0U;
         ++factor_index) {
      const auto node_id = program.integral_factor_nodes[
          static_cast<std::size_t>(term.integral_factor_nodes.offset +
                                   factor_index)];
      const auto &factor_node =
          program.nodes[static_cast<std::size_t>(node_id)];
      scratch.factor_scales.resize(active_count);
      for (std::size_t j = 0; j < active_count; ++j) {
        scratch.factor_scales[j] = std::fabs(scratch.products[j] *
            (has_common_sources ? scratch.child_values[scratch.positions[j]] : 1.0));
      }
      evaluate_compiled_lane_integral_node(
          plan,
          factor_node,
          executor,
          frame,
          scratch.lanes.data(),
          active_count,
          &scratch.values,
          integral_depth + 1U, scratch.factor_scales.data());
      std::size_t next_count = 0U;
      for (std::size_t j = 0; j < active_count; ++j) {
        const double product = scratch.products[j] * scratch.values[j];
        if (std::isfinite(product) && product != 0.0) {
          scratch.products[next_count] = product;
          scratch.lanes[next_count] = scratch.lanes[j];
          scratch.positions[next_count] = scratch.positions[j];
          ++next_count;
        }
      }
      active_count = next_count;
    }
    for (semantic::Index factor_index = 0;
         factor_index < term.expr_upper_factors.size && active_count > 0U;
         ++factor_index) {
      const auto &factor = program.expr_upper_factors[
          static_cast<std::size_t>(term.expr_upper_factors.offset +
                                   factor_index)];
      const auto &factor_node =
          program.nodes[static_cast<std::size_t>(factor.node_id)];
      std::size_t next_count = 0U;
      for (std::size_t j = 0; j < active_count; ++j) {
        const auto lane = static_cast<std::size_t>(scratch.lanes[j]);
        double product = scratch.products[j];
        const auto upper = compiled_lane_expr_upper_bound(
            factor_node, executor, *frame, lane);
        const bool has_upper = upper.found;
        const double time = compiled_lane_node_time(
            factor_node, *frame, lane);
        if (factor.mode == CompiledMathIntegralExprUpperMode::AfterOne) {
          if (!has_upper || time < upper.time) {
            continue;
          }
        } else if (has_upper) {
          if (time >= upper.time) {
            continue;
          }
          product /= upper.normalizer;
        }
        if (std::isfinite(product) && product != 0.0) {
          scratch.products[next_count] = product;
          scratch.lanes[next_count] = scratch.lanes[j];
          scratch.positions[next_count] = scratch.positions[j];
          ++next_count;
        }
      }
      active_count = next_count;
    }
    if (active_count == 0U) {
      continue;
    }
    for (std::size_t j = 0; j < active_count; ++j) {
      (*out)[scratch.positions[j]] += scratch.products[j];
    }
  }
  if (has_common_sources) {
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const double value = (*out)[lane] * scratch.child_values[lane];
      (*out)[lane] = std::isfinite(value) ? value : 0.0;
    }
  }
}

inline void evaluate_compiled_lane_integral_node(
    const ExactVariantPlan &plan,
    const CompiledMathNode &node,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *parent,
    const semantic::Index *lanes,
    const std::size_t lane_count,
    std::vector<double> *out,
    const std::size_t integral_depth,
    const double *absolute_scales) {
  const auto &program = plan.compiled_math;
  const auto &kernel = program.integral_kernels[
      static_cast<std::size_t>(node.integral_kernel_slot)];
  auto &scratch = executor->scratch(integral_depth * 4U);
  scratch.ensure(kAdaptiveNodeBatchSize);
  scratch.integral_lower.assign(lane_count, 0.0);
  scratch.integral_upper.resize(lane_count);
  const auto cached = [&](std::size_t request, std::size_t pos) {
    if (parent->cache_epoch[pos] != parent->cache_current_epoch) return false;
    const double scale = absolute_scales == nullptr ? 1.0 : absolute_scales[request];
    const auto &estimate = parent->cache_estimates[pos];
    const auto lane = static_cast<std::size_t>(compiled_frame_lane(lanes, request));
    const double relative = compiled_lane_node_time(node, *parent, lane) ==
            std::numeric_limits<double>::infinity()
        ? executor->terminal_relative_tolerance : kAdaptiveRelativeTolerance;
    return adaptive_integral_converged(estimate.value, estimate.error,
        kAdaptiveAbsoluteTolerance / scale, relative);
  };
  bool has_infinite_upper = false;
  for (std::size_t request = 0U; request < lane_count; ++request) {
    const auto lane = static_cast<std::size_t>(compiled_frame_lane(lanes, request));
    const auto cache_pos = parent->integral_cache_pos(node.integral_kernel_slot, lane);
    scratch.integral_upper[request] =
        cached(request, cache_pos)
            ? 0.0 : compiled_lane_node_time(node, *parent, lane);
    has_infinite_upper = has_infinite_upper || scratch.integral_upper[request] ==
        std::numeric_limits<double>::infinity();
  }
  // PDF/CDF factors are exactly zero before their source's support begins.
  if (!parent->has_sequence_history) {
    for (const auto leaf : kernel.support_leaves) {
      const auto input = executor->lanes.source_state().leaf_batch(leaf);
      for (std::size_t i = 0; i < lane_count; ++i) {
        const auto lane = static_cast<std::size_t>(compiled_frame_lane(lanes, i));
        const auto source = executor->lanes.source_lane(*parent, lane);
        const auto row = input.physical_row(source);
        scratch.integral_lower[i] = std::max(scratch.integral_lower[i],
            input.onset(row) + input.t0(row));
      }
    }
  }
  if (has_infinite_upper) {
    scratch.infinite_integrals.assign(lane_count, 0U);
    scratch.integral_origins = scratch.integral_lower;
    for (std::size_t i = 0; i < lane_count; ++i) {
      if (scratch.integral_upper[i] == std::numeric_limits<double>::infinity()) {
        scratch.infinite_integrals[i] = 1U;
        scratch.integral_lower[i] = 0.0;
        scratch.integral_upper[i] = 1.0;
      }
    }
  }
  const bool measured = has_infinite_upper && !parent->has_sequence_history && kernel.density_measure >= 0;
  if (measured) {
    // A monotone probability factor that is zero at its maximum makes the
    // entire integral zero, regardless of how narrow the PDF is.
    for (std::size_t begin = 0; begin < lane_count; begin += kAdaptiveNodeBatchSize) {
      const auto count = std::min(kAdaptiveNodeBatchSize, lane_count - begin);
      for (std::size_t i = 0; i < count; ++i)
        scratch.parents[i] = compiled_frame_lane(lanes, begin + i);
      auto &child = executor->lanes.integral_frame(integral_depth, count);
      child.copy_lanes_from(*parent, scratch.parents.data(), count);
      for (semantic::Index k = 0; k < kernel.density_factors.size; ++k) {
        const auto &op = program.source_product_ops[kernel.density_factors.offset + k];
        for (std::size_t i = 0; i < count; ++i)
          scratch.times[i] = op.value_channel_mask == kLeafChannelCdf
              ? (scratch.infinite_integrals[begin + i] ? std::numeric_limits<double>::infinity()
                  : scratch.integral_upper[begin + i]) : scratch.integral_origins[begin + i];
        child.set_time_plane(kernel.bind_time_id, scratch.times.data(), count);
        for (std::size_t i = 0; i < count; ++i) {
          const double t = child.time(op.time_id, i);
          scratch.times[i] = op.time_cap_id < 0 ? t : std::min(t, child.time(op.time_cap_id, i));
        }
        const auto &fill = evaluate_source_program_lanes(program, op.source_product_program_id,
            &executor->lanes, child, nullptr, scratch.times.data(), count,
            op.value_channel_mask, &executor->sources, 0);
        const auto values = source_lane_fill_values(fill, op.value_channel_mask);
        bool active = false;
        for (std::size_t i = 0; i < count; ++i) {
          if (values[i] == 0.0) scratch.integral_upper[begin + i] = scratch.integral_lower[begin + i];
          active = active || scratch.integral_upper[begin + i] > scratch.integral_lower[begin + i];
        }
        if (!active) break;
      }
    }
  }
  auto evaluate = [&](const double *times, const std::size_t *requests,
                      const std::size_t child_count, double *values) {
    for (std::size_t i = 0U; i < child_count; ++i) {
      scratch.parents[i] = compiled_frame_lane(lanes, requests[i]);
    }
    if (has_infinite_upper) {
      for (std::size_t i = 0; i < child_count; ++i) {
        scratch.times[i] = scratch.infinite_integrals[requests[i]]
            ? scratch.integral_origins[requests[i]] + times[i] / (1.0-times[i]) : times[i];
      }
    }
    auto &child = executor->lanes.integral_frame(integral_depth, child_count);
    child.copy_lanes_from(*parent, scratch.parents.data(), child_count);
    child.set_time_plane(kernel.bind_time_id,
        has_infinite_upper ? scratch.times.data() : times, child_count);
    auto &child_values = scratch.child_values;
    if (measured) {
      evaluate_source_product_ops_lanes(plan, kernel.density_factors, executor, &child,
          nullptr, child_count, &child_values, (integral_depth + 1U) * 4U);
      const auto &density = evaluate_source_program_lanes(program, kernel.density_measure,
          &executor->lanes, child, nullptr, scratch.times.data(), child_count,
          kLeafChannelPdf, &executor->sources, 0);
      scratch.density_values.resize(scratch.adaptive.values.size());
      const auto offset = values - scratch.adaptive.values.data();
      for (std::size_t i = 0; i < child_count; ++i) {
        const double v = 1.0 - times[i];
        const double mass = density.pdf[i] / (scratch.infinite_integrals[requests[i]] ? v * v : 1.0);
        scratch.density_values[offset + i] = mass;
        values[i] = mass * child_values[i];
      }
      return;
    }
    const auto &execution = child.has_sequence_history
                                ? kernel.execution
                                : kernel.initial_execution;
    evaluate_compiled_lane_execution(plan, program.roots[kernel.root_id], execution,
        executor, &child, child_count, &child_values, integral_depth + 1U);

    for (std::size_t i = 0U; i < child_count; ++i) {
      values[i] = std::isfinite(child_values[i]) ? child_values[i] : 0.0;
    }
    if (has_infinite_upper) {
      for (std::size_t i = 0; i < child_count; ++i) {
        if (scratch.infinite_integrals[requests[i]]) values[i] /= (1.0-times[i])*(1.0-times[i]);
      }
    }
  };
  const bool shared = !parent->has_sequence_history &&
      !has_infinite_upper &&
          lane_count >= kCumulativeMinimumRequests &&
          !kernel.cumulative_leaves.empty()
      ? prepare_cumulative_lanes(plan, kernel, executor->lanes.source_state(),
          *parent, lanes, lane_count, scratch.integral_lower.data(),
          scratch.integral_upper.data(), &scratch.cumulative)
      : false;
  const double *scales = absolute_scales;
  if (shared && absolute_scales != nullptr) {
    scratch.integral_scales.assign(absolute_scales, absolute_scales + lane_count);
    for (auto i = scratch.cumulative.order.rbegin(); i != scratch.cumulative.order.rend(); ++i) {
      const auto previous = scratch.cumulative.previous[*i];
      if (previous != lane_count) scratch.integral_scales[previous] = std::max(
          scratch.integral_scales[previous], scratch.integral_scales[*i]);
    }
    scales = scratch.integral_scales.data();
  }
  if (has_infinite_upper) {
    const auto refine = [&](AdaptiveLaneWorkspace *work) {
      if (!measured) return;
      // The embedded rules can both miss a narrow PDF. Compare its quadrature
      // mass with the CDF over each panel, using the same distribution interface.
      const auto map_time = [&](double u, std::size_t request) {
        return scratch.infinite_integrals[request]
            ? scratch.integral_origins[request] + u / (1.0 - u) : u;
      };
      for (std::size_t begin = 0; begin < work->panels.size(); begin += kAdaptiveNodeBatchSize / 2) {
        const auto count = std::min(kAdaptiveNodeBatchSize / 2, work->panels.size() - begin);
        for (std::size_t i = 0; i < count; ++i) {
          const auto &panel = work->panels[begin + i];
          scratch.lanes[2*i] = scratch.lanes[2*i+1] = compiled_frame_lane(lanes, panel.lane);
          scratch.times[2*i] = map_time(panel.lower, panel.lane);
          scratch.times[2*i+1] = map_time(panel.upper, panel.lane);
        }
        const auto &bounds = evaluate_source_program_lanes(program, kernel.density_measure,
            &executor->lanes, *parent, scratch.lanes.data(), scratch.times.data(), 2 * count,
            kLeafChannelCdf | kLeafChannelSurvival, &executor->sources, 0);
        for (std::size_t i = 0; i < count; ++i) {
          auto &panel = work->panels[begin + i];
          auto mass = panel;
          adaptive_finish_panel<21>(&mass, scratch.density_values.data() + 21 * (begin + i));
          const double exact = bounds.cdf[2*i] < 0.5
              ? bounds.cdf[2*i+1] - bounds.cdf[2*i]
              : bounds.survival[2*i] - bounds.survival[2*i+1];
          panel.error = std::max(panel.error, std::fabs(mass.value - exact));
        }
      }
    };
    adaptive_integrate_lane_batch<21, 2>(lane_count, scratch.integral_lower.data(),
        scratch.integral_upper.data(), evaluate, &scratch.adaptive, out,
        kAdaptiveAbsoluteTolerance, executor->terminal_relative_tolerance, scales, nullptr, refine);
  } else {
    adaptive_integrate_lane_batch(lane_count, scratch.integral_lower.data(),
        scratch.integral_upper.data(), evaluate, &scratch.adaptive, out,
        kAdaptiveAbsoluteTolerance, kAdaptiveRelativeTolerance, scales,
        shared ? scratch.cumulative.tolerance_divisors.data() : nullptr);
  }
  if (shared) {
    finish_cumulative_lanes(lane_count, scratch.integral_upper.data(), evaluate, &scratch.adaptive,
        &scratch.cumulative, out, absolute_scales);
  }
  for (std::size_t request = 0; request < lane_count; ++request) {
    const auto lane = static_cast<std::size_t>(
        compiled_frame_lane(lanes, request));
    const auto cache_pos = parent->integral_cache_pos(node.integral_kernel_slot, lane);
    if (scratch.integral_upper[request] == 0.0 && cached(request, cache_pos)) {
      (*out)[request] = parent->cache_estimates[cache_pos].value;
      continue;
    }
    if (!(compiled_lane_node_time(node, *parent, lane) > 0.0)) {
      continue;
    }
    double value = (*out)[request];
    // A raw integral may still be multiplied by a large external density.
    value = node.kind == CompiledMathNodeKind::IntegralZeroToCurrent
                ? clamp_probability(value)
                : (std::isfinite(value) ? value : 0.0);
    (*out)[request] = value;
    parent->cache_epoch[cache_pos] = parent->cache_current_epoch;
    parent->cache_estimates[cache_pos] = {value, shared
        ? scratch.cumulative.errors[request] : scratch.adaptive.lanes[request].error};
  }
}

inline void evaluate_compiled_lane_schedule(
    const ExactVariantPlan &plan,
    const CompiledMathIndexSpan schedule,
    const semantic::Index result_node_id,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *frame,
    const std::size_t lane_count,
    std::vector<double> *out,
    const std::size_t integral_depth) {
  if (result_node_id == semantic::kInvalidIndex || lane_count == 0U) {
    out->assign(lane_count, 0.0);
    return;
  }
  out->resize(lane_count);
  const auto &program = plan.compiled_math;
  auto &scratch = executor->scratch(integral_depth * 4U);
  scratch.ensure(lane_count);
  for (semantic::Index schedule_index = 0;
       schedule_index < schedule.size;
       ++schedule_index) {
    const auto node_id = program.root_schedule_nodes[
        static_cast<std::size_t>(schedule.offset + schedule_index)];
    const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
    double *node_out = frame->values_for(node_id);
    switch (node.kind) {
    case CompiledMathNodeKind::Constant:
      std::fill_n(node_out, lane_count, node.constant);
      break;
    case CompiledMathNodeKind::SourcePdf:
    case CompiledMathNodeKind::SourceCdf:
    case CompiledMathNodeKind::SourceSurvival: {
      for (std::size_t i = 0; i < lane_count; ++i) {
        scratch.times[i] = compiled_lane_source_node_time(node, *frame, i);
      }
      const auto value_mask = compiled_math_source_factor_channel_mask(node.kind);
      const auto fill_mask = exact_source_node_fill_mask(node.kind);
      const auto &source_values = evaluate_source_program_cached_lanes(
          program,
          node.source_program_id,
          &executor->lanes,
          frame,
          nullptr,
          scratch.times.data(),
          lane_count,
          fill_mask,
          node.cache_source_program,
          &executor->sources,
          &scratch.source_values);
      const double *values = source_lane_fill_values(
          source_values, value_mask);
      std::copy_n(values, lane_count, node_out);
      break;
    }
    case CompiledMathNodeKind::TimeGate:
    case CompiledMathNodeKind::StrictTimeGate: {
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset)];
      const double *child = frame->values_for(child_id);
      for (std::size_t i = 0; i < lane_count; ++i) {
        const bool open = node.kind == CompiledMathNodeKind::StrictTimeGate
                              ? frame->time(node.time_id, i) >
                                    frame->time(node.aux_id, i)
                              : frame->time(node.time_id, i) >=
                                    frame->time(node.aux_id, i);
        node_out[i] = open ? child[i] : 0.0;
      }
      break;
    }
    case CompiledMathNodeKind::OutcomeSelect: {
      // Branch roots have their own schedules; evaluate only selected lanes.
      for (std::size_t branch = 0; branch < 2U; ++branch) {
        std::size_t count = 0;
        for (std::size_t lane = 0; lane < lane_count; ++lane) {
          if (compiled_lane_outcome_used(node, plan, *frame, lane) == (branch == 1U)) {
            scratch.positions[count++] = lane;
          }
        }
        if (count == 0U) continue;
        const auto &branch_node = program.nodes[
            program.child_nodes[node.children.offset + branch]];
        if (branch_node.kind == CompiledMathNodeKind::Constant) {
          for (std::size_t i = 0; i < count; ++i) {
            node_out[scratch.positions[i]] = branch_node.constant;
          }
          continue;
        }
        auto &child = executor->lanes.integral_frame(integral_depth, count);
        child.copy_lanes_from(*frame, scratch.positions.data(), count);
        evaluate_compiled_lane_root(plan, node.branch_roots[branch], executor,
            &child, count, &scratch.values, integral_depth + 1U);
        for (std::size_t i = 0; i < count; ++i) {
          node_out[scratch.positions[i]] = scratch.values[i];
        }
      }
      break;
    }
    case CompiledMathNodeKind::IntegralZeroToCurrent:
    case CompiledMathNodeKind::IntegralZeroToCurrentRaw:
      evaluate_compiled_lane_integral_node(
          plan,
          node,
          executor,
          frame,
          nullptr,
          lane_count,
          &scratch.values,
          integral_depth);
      std::copy_n(scratch.values.data(), lane_count, node_out);
      break;
    case CompiledMathNodeKind::ExprUpperBoundDensity:
    case CompiledMathNodeKind::ExprUpperBoundCdf: {
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset)];
      const double *child = frame->values_for(child_id);
      for (std::size_t i = 0; i < lane_count; ++i) {
        const double raw = child[i];
        const auto upper = compiled_lane_expr_upper_bound(
            node, executor, *frame, i);
        if (!upper.found) {
          node_out[i] = raw;
          continue;
        }
        const double time = compiled_lane_node_time(node, *frame, i);
        if (node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
          node_out[i] =
              time >= upper.time
                  ? 1.0
                  : clamp_probability(raw / upper.normalizer);
        } else {
          node_out[i] =
              time >= upper.time
                  ? 0.0
                  : safe_density(raw / upper.normalizer);
        }
      }
      break;
    }
    case CompiledMathNodeKind::Product:
      std::fill_n(node_out, lane_count, 1.0);
      for (semantic::Index child_index = 0;
           child_index < node.children.size;
           ++child_index) {
        const auto child_id = program.child_nodes[
            static_cast<std::size_t>(node.children.offset + child_index)];
        const double *child = frame->values_for(child_id);
        for (std::size_t i = 0; i < lane_count; ++i) {
          const double value = node_out[i] * child[i];
          node_out[i] = std::isfinite(value) ? value : 0.0;
        }
      }
      break;
    case CompiledMathNodeKind::Sum:
    case CompiledMathNodeKind::CleanSignedSum:
      std::fill_n(node_out, lane_count, 0.0);
      for (semantic::Index child_index = 0;
           child_index < node.children.size;
           ++child_index) {
        const auto child_id = program.child_nodes[
            static_cast<std::size_t>(node.children.offset + child_index)];
        const double *child = frame->values_for(child_id);
        for (std::size_t i = 0; i < lane_count; ++i) {
          node_out[i] += child[i];
        }
      }
      if (node.kind == CompiledMathNodeKind::CleanSignedSum) {
        for (std::size_t i = 0; i < lane_count; ++i) {
          node_out[i] = clean_signed_value(node_out[i]);
        }
      }
      break;
    case CompiledMathNodeKind::ClampProbability:
    case CompiledMathNodeKind::Complement:
    case CompiledMathNodeKind::Negate: {
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset)];
      const double *child = frame->values_for(child_id);
      for (std::size_t i = 0; i < lane_count; ++i) {
        node_out[i] =
            node.kind == CompiledMathNodeKind::ClampProbability
                ? clamp_probability(child[i])
                : (node.kind == CompiledMathNodeKind::Complement
                       ? clamp_probability(1.0 - child[i])
                       : clean_signed_value(-child[i]));
      }
      break;
    }
    }
  }
  std::copy_n(frame->values_for(result_node_id), lane_count, out->begin());
}

inline void evaluate_compiled_lane_execution(
    const ExactVariantPlan &plan,
    const CompiledMathRoot &root,
    const CompiledMathExecutionPlan &execution,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *frame,
    const std::size_t lane_count,
    std::vector<double> *out,
    const std::size_t integral_depth) {
  if (execution.kind == CompiledMathExecutionKind::SourceProduct) {
    evaluate_source_product_ops_lanes(
        plan,
        execution.source_product_ops,
        executor,
        frame,
        nullptr,
        lane_count,
        out,
        integral_depth * 4U);
  } else if (execution.kind == CompiledMathExecutionKind::SourceProductSum) {
    evaluate_source_product_terms(
        plan,
        execution,
        executor,
        frame,
        lane_count,
        out,
        integral_depth);
  } else {
    evaluate_compiled_lane_schedule(
        plan,
        root.schedule,
        root.node_id,
        executor,
        frame,
        lane_count,
        out,
        integral_depth);
  }
  if (execution.clean_signed_source_sum) {
    for (auto &value : *out) value = clean_signed_value(value);
  }
}

inline void evaluate_compiled_lane_root(
    const ExactVariantPlan &plan,
    const semantic::Index root_id,
    CompiledLaneExecutor *executor,
    CompiledLaneFrame *frame,
    const std::size_t lane_count,
    std::vector<double> *out,
    const std::size_t integral_depth = 0U) {
  if (root_id == semantic::kInvalidIndex || lane_count == 0U) {
    out->assign(lane_count, 0.0);
    return;
  }
  const auto &root = plan.compiled_math.roots[root_id];
  evaluate_compiled_lane_execution(plan, root,
      frame->has_sequence_history ? root.execution : root.initial_execution,
      executor, frame, lane_count, out, integral_depth);
}

} // namespace detail
} // namespace accumulatr::eval
