#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <vector>

#include "exact_source_lane_eval.hpp"

namespace accumulatr::eval {
namespace detail {

struct CompiledLaneScratch {
  void ensure(const std::size_t count) {
    lanes.resize(count);
    positions.resize(count);
    parents.resize(count);
    times.resize(count);
    weights.resize(count);
    values.resize(count);
    child_values.resize(count);
    products.resize(count);
  }

  std::vector<semantic::Index> lanes;
  std::vector<std::size_t> positions;
  std::vector<std::size_t> parents;
  std::vector<double> times;
  std::vector<double> weights;
  std::vector<double> values;
  std::vector<double> child_values;
  std::vector<double> products;
  SourceLaneFill source_values;
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

inline bool compiled_lane_outcome_gate_open(
    const CompiledMathNode &node,
    const ExactVariantPlan &plan,
    const CompiledLaneFrame &frame,
    const std::size_t lane) {
  if (!frame.has_sequence_history) {
    return node.kind == CompiledMathNodeKind::OutcomeSubsetUnused;
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
  return node.kind == CompiledMathNodeKind::OutcomeSubsetUsed
             ? any_used
             : !any_used;
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
    const std::size_t integral_depth);

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
         i < term.outcome_gate_nodes.size && active_count > 0U;
         ++i) {
      const auto node_id = program.outcome_gate_nodes[
          static_cast<std::size_t>(term.outcome_gate_nodes.offset + i)];
      const auto &gate = program.nodes[static_cast<std::size_t>(node_id)];
      if (!frame->has_sequence_history) {
        if (gate.kind == CompiledMathNodeKind::OutcomeSubsetUsed) {
          active_count = 0U;
        }
        continue;
      }
      std::size_t next_count = 0U;
      for (std::size_t j = 0; j < active_count; ++j) {
        const auto lane = scratch.positions[j];
        if (compiled_lane_outcome_gate_open(gate, plan, *frame, lane)) {
          scratch.positions[next_count++] = lane;
        }
      }
      active_count = next_count;
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
      scratch.products[j] = term.sign;
    }
    if (active_count == 0U) {
      continue;
    }

    for (semantic::Index factor_index = 0;
         factor_index < term.integral_factor_nodes.size && active_count > 0U;
         ++factor_index) {
      const auto node_id = program.integral_factor_nodes[
          static_cast<std::size_t>(term.integral_factor_nodes.offset +
                                   factor_index)];
      const auto &factor_node =
          program.nodes[static_cast<std::size_t>(node_id)];
      evaluate_compiled_lane_integral_node(
          plan,
          factor_node,
          executor,
          frame,
          scratch.lanes.data(),
          active_count,
          &scratch.values,
          integral_depth + 1U);
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
    evaluate_source_product_ops_lanes(
        plan,
        term.source_product_ops,
        executor,
        frame,
        scratch.lanes.data(),
        active_count,
        &scratch.values,
        integral_depth * 4U + 2U);
    for (std::size_t j = 0; j < active_count; ++j) {
      (*out)[scratch.positions[j]] +=
          scratch.products[j] * scratch.values[j];
    }
  }
  if (has_common_sources) {
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const double value = (*out)[lane] * scratch.child_values[lane];
      (*out)[lane] =
          std::isfinite(value) && value != 0.0 ? value : 0.0;
    }
  }
  if (execution.clean_signed_source_sum) {
    for (auto &value : *out) {
      value = clean_signed_value(value);
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
  const std::size_t integral_depth) {
  const auto &program = plan.compiled_math;
  const auto &kernel = program.integral_kernels[
      static_cast<std::size_t>(node.integral_kernel_slot)];
  auto &scratch = executor->scratch(integral_depth * 4U);
  constexpr auto parent_tile =
      exact_expanded_parent_lane_tile_size<
          quadrature::kGenericFiniteOrder>();
  scratch.ensure(kExactExpandedLaneTileSize);
  out->assign(lane_count, 0.0);
  bool any_children = false;
  for (std::size_t request_begin = 0U;
       request_begin < lane_count;
       request_begin += parent_tile) {
    const auto request_end =
        std::min(lane_count, request_begin + parent_tile);
    std::size_t child_count = 0U;
    for (std::size_t request = request_begin;
         request < request_end;
         ++request) {
      const auto lane = static_cast<std::size_t>(
          compiled_frame_lane(lanes, request));
      const double upper = compiled_lane_node_time(node, *parent, lane);
      if (!(upper > 0.0)) {
        continue;
      }
      const auto cache_pos = parent->integral_cache_pos(
          node.integral_kernel_slot, lane);
      if (parent->cache_epoch[cache_pos] == parent->cache_current_epoch) {
        (*out)[request] = parent->cache_values[cache_pos];
        continue;
      }
      const auto nodes = quadrature::map_rule_to_finite_interval<
          quadrature::kGenericFiniteOrder>(0.0, upper);
      for (std::size_t q = 0; q < quadrature::kGenericFiniteOrder; ++q) {
        scratch.parents[child_count] = lane;
        scratch.positions[child_count] = request;
        scratch.times[child_count] = nodes.nodes[q];
        scratch.weights[child_count] = nodes.weights[q];
        ++child_count;
      }
    }
    if (child_count == 0U) {
      continue;
    }
    any_children = true;

    auto &child = executor->lanes.integral_frame(integral_depth, child_count);
    child.copy_lanes_from(*parent, scratch.parents.data(), child_count);
    child.set_time_plane(
        kernel.bind_time_id, scratch.times.data(), child_count);
    auto &child_values = scratch.child_values;
    const auto &execution = child.has_sequence_history
                                ? kernel.execution
                                : kernel.initial_execution;
    if (execution.kind == CompiledMathExecutionKind::Schedule) {
      const auto &root =
          program.roots[static_cast<std::size_t>(kernel.root_id)];
      evaluate_compiled_lane_schedule(
          plan,
          root.schedule,
          root.node_id,
          executor,
          &child,
          child_count,
          &child_values,
          integral_depth + 1U);
    } else if (execution.kind == CompiledMathExecutionKind::SourceProduct) {
      evaluate_source_product_ops_lanes(
          plan,
          execution.source_product_ops,
          executor,
          &child,
          nullptr,
          child_count,
          &child_values,
          integral_depth * 4U + 3U);
      if (execution.clean_signed_source_sum) {
        for (auto &value : child_values) {
          value = clean_signed_value(value);
        }
      }
    } else {
      evaluate_source_product_terms(
          plan,
          execution,
          executor,
          &child,
          child_count,
          &child_values,
          integral_depth + 1U);
    }

    for (std::size_t child_index = 0;
         child_index < child_count;
         ++child_index) {
      const auto request = scratch.positions[child_index];
      const double value = child_values[child_index];
      if (std::isfinite(value) && value != 0.0) {
        (*out)[request] += scratch.weights[child_index] * value;
      }
    }
  }
  if (!any_children) {
    return;
  }
  for (std::size_t request = 0; request < lane_count; ++request) {
    const auto lane = static_cast<std::size_t>(
        compiled_frame_lane(lanes, request));
    if (!(compiled_lane_node_time(node, *parent, lane) > 0.0)) {
      continue;
    }
    double value = (*out)[request];
    value = node.kind == CompiledMathNodeKind::IntegralZeroToCurrent
                ? clamp_probability(value)
                : (std::isfinite(value) ? clean_signed_value(value) : 0.0);
    (*out)[request] = value;
    const auto cache_pos = parent->integral_cache_pos(
        node.integral_kernel_slot, lane);
    parent->cache_epoch[cache_pos] = parent->cache_current_epoch;
    parent->cache_values[cache_pos] = value;
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
    case CompiledMathNodeKind::ExprDensity:
    case CompiledMathNodeKind::ExprCdf:
    case CompiledMathNodeKind::ExprSurvival:
      std::terminate();
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
    case CompiledMathNodeKind::OutcomeSubsetUnused:
    case CompiledMathNodeKind::OutcomeSubsetUsed:
      for (std::size_t i = 0; i < lane_count; ++i) {
        node_out[i] = compiled_lane_outcome_gate_open(
                             node, plan, *frame, i)
                             ? 1.0
                             : 0.0;
      }
      break;
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
          node_out[i] =
              std::isfinite(value) && value != 0.0 ? value : 0.0;
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
  const auto &root =
      plan.compiled_math.roots[static_cast<std::size_t>(root_id)];
  const auto &execution = frame->has_sequence_history
                              ? root.execution
                              : root.initial_execution;
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
    if (execution.clean_signed_source_sum) {
      for (auto &value : *out) {
        value = clean_signed_value(value);
      }
    }
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
}

} // namespace detail
} // namespace accumulatr::eval
