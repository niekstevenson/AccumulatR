#pragma once

#include <cmath>
#include <stdexcept>

#include "exact_source_product_eval.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool compiled_math_outcome_gate_open_for_node(
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const CompiledSourceView *parent) {
  if (node.kind != CompiledMathNodeKind::OutcomeSubsetUnused &&
      node.kind != CompiledMathNodeKind::OutcomeSubsetUsed) {
    throw std::runtime_error("integral outcome gate contains a non-gate node");
  }
  if (workspace.used_outcomes == nullptr) {
    return node.kind == CompiledMathNodeKind::OutcomeSubsetUnused;
  }
  if (parent == nullptr) {
    throw std::runtime_error(
        "integral outcome gate requires a compiled source view");
  }
  const auto offset = static_cast<std::size_t>(node.subject_id);
  const auto size = static_cast<std::size_t>(node.aux_id);
  const auto &indices = parent->plan().compiled_outcome_gate_indices;
  if (offset + size > indices.size()) {
    throw std::runtime_error(
        "compiled outcome-used gate points outside the plan");
  }
  bool any_used = false;
  for (std::size_t i = 0; i < size; ++i) {
    const auto outcome_idx = indices[offset + i];
    if (outcome_idx != semantic::kInvalidIndex &&
        static_cast<std::size_t>(outcome_idx) <
            workspace.used_outcomes->size() &&
        (*workspace.used_outcomes)[static_cast<std::size_t>(outcome_idx)] !=
            0U) {
      any_used = true;
      break;
    }
  }
  return node.kind == CompiledMathNodeKind::OutcomeSubsetUsed
             ? any_used
             : !any_used;
}

inline bool compiled_math_integral_outcome_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan outcome_gate_nodes,
    const CompiledMathWorkspace &workspace,
    const CompiledSourceView *parent) {
  for (semantic::Index i = 0; i < outcome_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_outcome_gate_nodes[
            static_cast<std::size_t>(outcome_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (!compiled_math_outcome_gate_open_for_node(
            gate_node,
            workspace,
            parent)) {
      return false;
    }
  }
  return true;
}

inline bool compiled_math_integral_time_gates_open(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan time_gate_nodes,
    const CompiledMathWorkspace &workspace) {
  for (semantic::Index i = 0; i < time_gate_nodes.size; ++i) {
    const auto node_id =
        program.integral_kernel_time_gate_nodes[
            static_cast<std::size_t>(time_gate_nodes.offset + i)];
    const auto &gate_node =
        program.nodes[static_cast<std::size_t>(node_id)];
    if (gate_node.kind != CompiledMathNodeKind::TimeGate &&
        gate_node.kind != CompiledMathNodeKind::StrictTimeGate) {
      throw std::runtime_error("integral time gate contains a non-gate node");
    }
    const auto current =
        workspace.time_values[static_cast<std::size_t>(gate_node.time_id)];
    const auto gate =
        workspace.time_values[static_cast<std::size_t>(gate_node.aux_id)];
    const bool open =
        gate_node.kind == CompiledMathNodeKind::StrictTimeGate
            ? current > gate
            : current >= gate;
    if (!open) {
      return false;
    }
  }
  return true;
}

inline double compiled_math_integral_factor_value_for_node(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace) {
  if (!compiled_math_is_integral_node(node.kind)) {
    throw std::runtime_error("integral factor contains a non-integral node");
  }
  const double current_time = compiled_math_node_time(node, *workspace);
  if (!(current_time > 0.0)) {
    return 0.0;
  }
  double value =
      evaluate_compiled_integral_kernel(
          program,
          node,
          workspace,
          parent,
          scenario,
          eval_workspace,
          0.0,
          current_time);
  if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent) {
    value = clamp_probability(value);
  }
  return value;
}

inline bool compiled_math_apply_expr_upper_factor(
    const CompiledMathProgram &program,
    const CompiledMathIntegralExprUpperFactor &factor,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *eval_workspace,
    double *product) {
  const auto &node =
      program.nodes[static_cast<std::size_t>(factor.node_id)];
  if (node.kind != CompiledMathNodeKind::ExprUpperBoundDensity &&
      node.kind != CompiledMathNodeKind::ExprUpperBoundCdf) {
    throw std::runtime_error(
        "integral expression upper factor contains a non-upper node");
  }
  const auto upper =
      compiled_math_expr_upper_bound_for_node(
          program,
          node,
          *workspace,
          compiled_math_node_evaluator(node, parent, eval_workspace));
  const bool has_upper = upper.found && upper.normalizer > 0.0;
  const double current_time = compiled_math_node_time(node, *workspace);
  if (factor.mode == CompiledMathIntegralExprUpperMode::AfterOne) {
    if (!has_upper || current_time < upper.time) {
      *product = 0.0;
      return false;
    }
    return true;
  }
  if (!has_upper) {
    return true;
  }
  if (current_time >= upper.time) {
    *product = 0.0;
    return false;
  }
  *product /= upper.normalizer;
  return *product != 0.0;
}

inline double evaluate_compiled_source_product_integral_kernel(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledEvalWorkspace *eval_workspace,
    const double lower,
    const double upper) {
  if (!(upper > lower)) {
    return 0.0;
  }
  auto *source_channels = parent == nullptr ? nullptr : parent->source_channels();
  if (source_channels == nullptr) {
    return 0.0;
  }
  auto &integral_workspace = workspace->integral_workspace_for(program);
  const std::vector<std::uint8_t> *term_outcome_open_ptr = nullptr;
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProductSum) {
    integral_workspace.integral_term_open.assign(
        static_cast<std::size_t>(kernel.source_product_terms.size), 1U);
    for (semantic::Index term_idx = 0;
         term_idx < kernel.source_product_terms.size;
         ++term_idx) {
      const auto &term =
          program.integral_kernel_source_product_terms[
              static_cast<std::size_t>(
                  kernel.source_product_terms.offset + term_idx)];
      integral_workspace.integral_term_open[static_cast<std::size_t>(term_idx)] =
          compiled_math_integral_outcome_gates_open(
              program,
              term.outcome_gate_nodes,
              integral_workspace,
              parent)
              ? 1U
              : 0U;
    }
    term_outcome_open_ptr = &integral_workspace.integral_term_open;
  }
  double sum = 0.0;
  CompiledMathWorkspace::RebindableTimeBinding bind(
      &integral_workspace,
      kernel.bind_time_id);
  auto evaluate_sample = [&]() {
    double value = 0.0;
    if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
      value = compiled_math_source_product_value_for_ops(
          program,
          kernel.source_product_ops,
          &integral_workspace,
          source_channels);
    } else {
      double total = 0.0;
      for (semantic::Index term_idx = 0;
           term_idx < kernel.source_product_terms.size;
           ++term_idx) {
        const auto &term =
            program.integral_kernel_source_product_terms[
                static_cast<std::size_t>(
                    kernel.source_product_terms.offset + term_idx)];
        const auto term_pos = static_cast<std::size_t>(term_idx);
        if (term_outcome_open_ptr != nullptr &&
            term_pos < term_outcome_open_ptr->size() &&
            (*term_outcome_open_ptr)[term_pos] == 0U) {
          continue;
        }
        if (!compiled_math_integral_time_gates_open(
                program, term.time_gate_nodes, integral_workspace)) {
          continue;
        }
        double product = term.sign;
        for (semantic::Index i = 0; i < term.integral_factor_nodes.size; ++i) {
          const auto node_id =
              program.integral_kernel_integral_factor_nodes[
                  static_cast<std::size_t>(
                      term.integral_factor_nodes.offset + i)];
          const auto &integral_node =
              program.nodes[static_cast<std::size_t>(node_id)];
          product *= compiled_math_integral_factor_value_for_node(
              program,
              integral_node,
              &integral_workspace,
              parent,
              nullptr,
              eval_workspace);
          if (product == 0.0) {
            product = 0.0;
            break;
          }
        }
        if (product == 0.0) {
          continue;
        }
        for (semantic::Index i = 0; i < term.expr_upper_factors.size; ++i) {
          const auto &factor =
              program.integral_kernel_expr_upper_factors[
                  static_cast<std::size_t>(
                      term.expr_upper_factors.offset + i)];
          if (!compiled_math_apply_expr_upper_factor(
                  program,
                  factor,
                  &integral_workspace,
                  parent,
                  eval_workspace,
                  &product)) {
            break;
          }
        }
        if (product == 0.0) {
          continue;
        }
        product *= compiled_math_source_product_value_for_ops(
            program,
            term.source_product_ops,
            &integral_workspace,
            source_channels);
        if (product == 0.0) {
          continue;
        }
        total += product;
      }
      value =
          kernel.clean_signed_source_sum ? clean_signed_value(total) : total;
    }
    return value;
  };
  auto integrate_nodes = [&](const auto &nodes) {
    double local_sum = 0.0;
    for (std::size_t sample_idx = 0;
         sample_idx < nodes.nodes.size();
         ++sample_idx) {
      bind.set(nodes.nodes[sample_idx]);
      integral_workspace.reset_source_product_program_cache();
      const double value = evaluate_sample();
      if (value == 0.0) {
        continue;
      }
      local_sum += nodes.weights[sample_idx] * value;
    }
    return local_sum;
  };
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProductSum) {
    const auto nodes =
        quadrature::map_rule_to_finite_interval<
            quadrature::kGenericFiniteOrder>(lower, upper);
    sum = integrate_nodes(nodes);
  } else {
    const auto batch = quadrature::build_finite_batch(lower, upper);
    sum = integrate_nodes(batch.nodes);
  }
  return sum;
}

inline double evaluate_compiled_generic_integral_kernel(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace,
    const double lower,
    const double upper) {
  if (!(upper > lower)) {
    return 0.0;
  }
  const auto root_pos = static_cast<std::size_t>(kernel.root_id);
  if (root_pos >= program.roots.size()) {
    throw std::runtime_error(
        "generic integral kernel references an unplanned integrand root");
  }
  auto &integral_workspace = workspace->integral_workspace_for(program);
  integral_workspace.ensure_size(program);
  const auto &root = program.roots[root_pos];
  const auto nodes =
      quadrature::map_rule_to_finite_interval<
          quadrature::kGenericFiniteOrder>(lower, upper);
  double sum = 0.0;
  CompiledMathWorkspace::RebindableTimeBinding bind(
      &integral_workspace,
      kernel.bind_time_id);
  for (std::size_t sample_idx = 0;
       sample_idx < quadrature::kGenericFiniteOrder;
       ++sample_idx) {
    bind.set(nodes.nodes[sample_idx]);
    integral_workspace.reset_cache();
    const double value = evaluate_compiled_node_span(
        program,
        root.schedule,
        root.node_id,
        &integral_workspace,
        parent,
        scenario,
        eval_workspace,
        nullptr,
        true);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    sum += nodes.weights[sample_idx] * value;
  }
  return std::isfinite(sum) ? sum : 0.0;
}

inline double evaluate_compiled_integral_kernel(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledMathWorkspace *workspace,
    CompiledSourceView *parent,
    CompiledSourceView *scenario,
    CompiledEvalWorkspace *eval_workspace,
    const double lower,
    const double upper) {
  if (node.integral_kernel_slot == semantic::kInvalidIndex) {
    throw std::runtime_error(
        "compiled integral node has no planned integral kernel");
  }
  const auto kernel_pos =
      static_cast<std::size_t>(node.integral_kernel_slot);
  if (kernel_pos >= program.integral_kernels.size()) {
    throw std::runtime_error(
        "compiled integral node points outside integral kernels");
  }
  const auto &kernel = program.integral_kernels[kernel_pos];
  if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct ||
      kernel.kind == CompiledMathIntegralKernelKind::SourceProductSum) {
    return evaluate_compiled_source_product_integral_kernel(
        program,
        kernel,
        workspace,
        parent,
        eval_workspace,
        lower,
        upper);
  }
  if (kernel.kind == CompiledMathIntegralKernelKind::Generic) {
    return evaluate_compiled_generic_integral_kernel(
        program,
        kernel,
        workspace,
        parent,
        scenario,
        eval_workspace,
        lower,
        upper);
  }
  throw std::runtime_error("compiled integral kernel has no planned evaluator");
}

} // namespace detail
} // namespace accumulatr::eval
