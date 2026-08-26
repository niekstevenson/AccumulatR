#pragma once

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "compiled_math_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool compiled_math_is_source_value_node(
    const CompiledMathNodeKind kind) noexcept {
  return kind == CompiledMathNodeKind::SourcePdf ||
         kind == CompiledMathNodeKind::SourceCdf ||
         kind == CompiledMathNodeKind::SourceSurvival;
}

inline bool compiled_math_is_integral_node(
    const CompiledMathNodeKind kind) noexcept {
  return kind == CompiledMathNodeKind::IntegralZeroToCurrent ||
         kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw;
}

inline semantic::Index compiled_math_integral_root_id(
    const CompiledMathNode &node) noexcept {
  return node.subject_id;
}

inline semantic::Index compiled_math_integral_bind_time_id(
    const CompiledMathNode &node) noexcept {
  return node.aux2_id == semantic::kInvalidIndex ? node.time_id : node.aux2_id;
}

inline std::uint8_t compiled_math_source_factor_channel_mask(
    const CompiledMathNodeKind kind) noexcept;

struct CompiledMathSourceProductTermBuild {
  std::vector<semantic::Index> source_value_nodes;
  std::vector<semantic::Index> outcome_gate_nodes;
  std::vector<semantic::Index> time_gate_nodes;
  std::vector<semantic::Index> integral_factor_nodes;
  std::vector<CompiledMathExprUpperFactor> expr_upper_factors;
  double sign{1.0};
};

inline bool compiled_math_expand_source_product_terms(
    const CompiledMathProgram &program,
    const semantic::Index node_id,
    const double sign,
    const bool initial_state,
    std::vector<CompiledMathSourceProductTermBuild> *terms,
    bool *clean_signed) {
  if (node_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &node = program.nodes[static_cast<std::size_t>(node_id)];
  if (compiled_math_is_source_value_node(node.kind)) {
    terms->push_back(
        CompiledMathSourceProductTermBuild{{node_id}, {}, {}, {}, {}, sign});
    return true;
  }
  if (node.kind == CompiledMathNodeKind::Constant) {
    const double term_sign = sign * node.constant;
    if (term_sign != 0.0) {
      terms->push_back(
          CompiledMathSourceProductTermBuild{{}, {}, {}, {}, {}, term_sign});
    }
    return true;
  }
  if (node.kind == CompiledMathNodeKind::OutcomeSubsetUnused ||
      node.kind == CompiledMathNodeKind::OutcomeSubsetUsed) {
    if (initial_state) {
      if (node.kind == CompiledMathNodeKind::OutcomeSubsetUnused) {
        terms->push_back(
            CompiledMathSourceProductTermBuild{{}, {}, {}, {}, {}, sign});
      }
      return true;
    }
    terms->push_back(
        CompiledMathSourceProductTermBuild{{}, {node_id}, {}, {}, {}, sign});
    return true;
  }
  if (compiled_math_is_integral_node(node.kind)) {
    if (node.integral_kernel_slot == semantic::kInvalidIndex) {
      return false;
    }
    terms->push_back(
        CompiledMathSourceProductTermBuild{{}, {}, {}, {node_id}, {}, sign});
    return true;
  }
  if ((node.kind == CompiledMathNodeKind::ExprUpperBoundDensity ||
       node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) &&
      node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    std::vector<CompiledMathSourceProductTermBuild> bounded_terms;
    if (!compiled_math_expand_source_product_terms(
            program,
            child_id,
            sign,
            initial_state,
            &bounded_terms,
            clean_signed)) {
      return false;
    }
    for (auto &term : bounded_terms) {
      term.expr_upper_factors.push_back(
          CompiledMathExprUpperFactor{
              node_id,
              CompiledMathIntegralExprUpperMode::BeforeScale});
    }
    terms->insert(
        terms->end(),
        std::make_move_iterator(bounded_terms.begin()),
        std::make_move_iterator(bounded_terms.end()));
    if (node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
      terms->push_back(
          CompiledMathSourceProductTermBuild{
              {},
              {},
              {},
              {},
              {CompiledMathExprUpperFactor{
                  node_id,
                  CompiledMathIntegralExprUpperMode::AfterOne}},
              sign});
    }
    return true;
  }
  if ((node.kind == CompiledMathNodeKind::TimeGate ||
       node.kind == CompiledMathNodeKind::StrictTimeGate) &&
      node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    std::vector<CompiledMathSourceProductTermBuild> gated_terms;
    if (!compiled_math_expand_source_product_terms(
            program,
            child_id,
            sign,
            initial_state,
            &gated_terms,
            clean_signed)) {
      return false;
    }
    for (auto &term : gated_terms) {
      term.time_gate_nodes.push_back(node_id);
    }
    terms->insert(
        terms->end(),
        std::make_move_iterator(gated_terms.begin()),
        std::make_move_iterator(gated_terms.end()));
    return true;
  }
  if (node.kind == CompiledMathNodeKind::Sum ||
      node.kind == CompiledMathNodeKind::CleanSignedSum) {
    if (node.kind == CompiledMathNodeKind::CleanSignedSum &&
        clean_signed != nullptr) {
      *clean_signed = true;
    }
    if (node.children.empty()) {
      return false;
    }
    for (semantic::Index i = 0; i < node.children.size; ++i) {
      const auto child_id = program.child_nodes[
          static_cast<std::size_t>(node.children.offset + i)];
      if (!compiled_math_expand_source_product_terms(
              program,
              child_id,
              sign,
              initial_state,
              terms,
              clean_signed)) {
        return false;
      }
    }
    return true;
  }
  if (node.kind == CompiledMathNodeKind::Negate && node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    return compiled_math_expand_source_product_terms(
        program, child_id, -sign, initial_state, terms, clean_signed);
  }
  if (node.kind == CompiledMathNodeKind::ClampProbability &&
      node.children.size == 1U) {
    if (clean_signed != nullptr) {
      *clean_signed = true;
    }
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    return compiled_math_expand_source_product_terms(
        program, child_id, sign, initial_state, terms, clean_signed);
  }
  if (node.kind == CompiledMathNodeKind::Complement &&
      node.children.size == 1U) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset)];
    terms->push_back(CompiledMathSourceProductTermBuild{{}, {}, {}, {}, {}, sign});
    return compiled_math_expand_source_product_terms(
        program, child_id, -sign, initial_state, terms, clean_signed);
  }
  if (node.kind != CompiledMathNodeKind::Product) {
    return false;
  }

  std::vector<CompiledMathSourceProductTermBuild> product_terms(1);
  product_terms.front().sign = sign;
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + i)];
    std::vector<CompiledMathSourceProductTermBuild> child_terms;
    if (!compiled_math_expand_source_product_terms(
            program,
            child_id,
            1.0,
            initial_state,
            &child_terms,
            clean_signed)) {
      return false;
    }
    std::vector<CompiledMathSourceProductTermBuild> next_terms;
    next_terms.reserve(product_terms.size() * child_terms.size());
    for (const auto &lhs : product_terms) {
      for (const auto &rhs : child_terms) {
        CompiledMathSourceProductTermBuild combined;
        combined.sign = lhs.sign * rhs.sign;
        combined.source_value_nodes.reserve(
            lhs.source_value_nodes.size() + rhs.source_value_nodes.size());
        combined.source_value_nodes.insert(
            combined.source_value_nodes.end(),
            lhs.source_value_nodes.begin(),
            lhs.source_value_nodes.end());
        combined.source_value_nodes.insert(
            combined.source_value_nodes.end(),
            rhs.source_value_nodes.begin(),
            rhs.source_value_nodes.end());
        combined.outcome_gate_nodes.reserve(
            lhs.outcome_gate_nodes.size() + rhs.outcome_gate_nodes.size());
        combined.outcome_gate_nodes.insert(
            combined.outcome_gate_nodes.end(),
            lhs.outcome_gate_nodes.begin(),
            lhs.outcome_gate_nodes.end());
        combined.outcome_gate_nodes.insert(
            combined.outcome_gate_nodes.end(),
            rhs.outcome_gate_nodes.begin(),
            rhs.outcome_gate_nodes.end());
        combined.time_gate_nodes.reserve(
            lhs.time_gate_nodes.size() + rhs.time_gate_nodes.size());
        combined.time_gate_nodes.insert(
            combined.time_gate_nodes.end(),
            lhs.time_gate_nodes.begin(),
            lhs.time_gate_nodes.end());
        combined.time_gate_nodes.insert(
            combined.time_gate_nodes.end(),
            rhs.time_gate_nodes.begin(),
            rhs.time_gate_nodes.end());
        combined.integral_factor_nodes.reserve(
            lhs.integral_factor_nodes.size() +
            rhs.integral_factor_nodes.size());
        combined.integral_factor_nodes.insert(
            combined.integral_factor_nodes.end(),
            lhs.integral_factor_nodes.begin(),
            lhs.integral_factor_nodes.end());
        combined.integral_factor_nodes.insert(
            combined.integral_factor_nodes.end(),
            rhs.integral_factor_nodes.begin(),
            rhs.integral_factor_nodes.end());
        combined.expr_upper_factors.reserve(
            lhs.expr_upper_factors.size() + rhs.expr_upper_factors.size());
        combined.expr_upper_factors.insert(
            combined.expr_upper_factors.end(),
            lhs.expr_upper_factors.begin(),
            lhs.expr_upper_factors.end());
        combined.expr_upper_factors.insert(
            combined.expr_upper_factors.end(),
            rhs.expr_upper_factors.begin(),
            rhs.expr_upper_factors.end());
        if (combined.sign != 0.0) {
          next_terms.push_back(std::move(combined));
        }
      }
    }
    product_terms = std::move(next_terms);
  }
  terms->insert(terms->end(), product_terms.begin(), product_terms.end());
  return true;
}

inline bool compiled_math_collect_source_product_terms(
    CompiledMathProgram *program,
    const semantic::Index node_id,
    const double sign,
    const bool initial_state,
    std::vector<CompiledMathSourceProductTerm> *terms,
    bool *clean_signed) {
  std::vector<CompiledMathSourceProductTermBuild> built_terms;
  if (!compiled_math_expand_source_product_terms(
      *program,
      node_id,
      sign,
      initial_state,
      &built_terms,
      clean_signed)) {
    return false;
  }
  if (built_terms.empty()) {
    CompiledMathSourceProductTerm term;
    term.sign = 0.0;
    terms->push_back(term);
    return true;
  }
  for (const auto &built : built_terms) {
    const auto factor_offset = static_cast<semantic::Index>(
        program->source_value_factors.size());
    for (const auto source_node_id : built.source_value_nodes) {
      const auto &source_node =
          program->nodes[static_cast<std::size_t>(source_node_id)];
      if (!compiled_math_is_source_value_node(source_node.kind)) {
        throw std::runtime_error(
            "source-product integral factor is not a planned source value");
      }
      program->source_value_factors.push_back(
          CompiledMathSourceValueFactor{
              source_node.subject_id,
              source_node.condition_id,
              source_node.source_view_id,
              source_node.time_id,
              source_node.aux_id,
              source_node.kind});
    }
    const auto gate_offset = static_cast<semantic::Index>(
        program->outcome_gate_nodes.size());
    program->outcome_gate_nodes.insert(
        program->outcome_gate_nodes.end(),
        built.outcome_gate_nodes.begin(),
        built.outcome_gate_nodes.end());
    const auto time_gate_offset = static_cast<semantic::Index>(
        program->time_gate_nodes.size());
    program->time_gate_nodes.insert(
        program->time_gate_nodes.end(),
        built.time_gate_nodes.begin(),
        built.time_gate_nodes.end());
    const auto integral_factor_offset = static_cast<semantic::Index>(
        program->integral_factor_nodes.size());
    program->integral_factor_nodes.insert(
        program->integral_factor_nodes.end(),
        built.integral_factor_nodes.begin(),
        built.integral_factor_nodes.end());
    const auto expr_upper_offset = static_cast<semantic::Index>(
        program->expr_upper_factors.size());
    program->expr_upper_factors.insert(
        program->expr_upper_factors.end(),
        built.expr_upper_factors.begin(),
        built.expr_upper_factors.end());
    CompiledMathSourceProductTerm term;
    term.source_value_factors = CompiledMathIndexSpan{
        factor_offset,
        static_cast<semantic::Index>(built.source_value_nodes.size())};
    term.outcome_gate_nodes = CompiledMathIndexSpan{
        gate_offset,
        static_cast<semantic::Index>(built.outcome_gate_nodes.size())};
    term.time_gate_nodes = CompiledMathIndexSpan{
        time_gate_offset,
        static_cast<semantic::Index>(built.time_gate_nodes.size())};
    term.integral_factor_nodes = CompiledMathIndexSpan{
        integral_factor_offset,
        static_cast<semantic::Index>(built.integral_factor_nodes.size())};
    term.expr_upper_factors = CompiledMathIndexSpan{
        expr_upper_offset,
        static_cast<semantic::Index>(built.expr_upper_factors.size())};
    term.sign = built.sign;
    terms->push_back(term);
  }
  return true;
}

inline bool compiled_math_same_source_product_channel(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathSourceValueFactor &factor,
    const semantic::Index effective_source_view_id) noexcept {
  return channel.source_id == factor.source_id &&
         channel.condition_id == factor.condition_id &&
         channel.source_view_id == effective_source_view_id &&
         channel.time_id == factor.time_id &&
         channel.time_cap_id == factor.time_cap_id;
}

inline std::uint8_t compiled_math_source_factor_channel_mask(
    const CompiledMathNodeKind kind) noexcept {
  switch (kind) {
  case CompiledMathNodeKind::SourcePdf:
    return 1U;
  case CompiledMathNodeKind::SourceCdf:
    return 2U;
  case CompiledMathNodeKind::SourceSurvival:
    return 4U;
  default:
    break;
  }
  return 0U;
}

inline void compiled_math_compile_source_product_channels(
    CompiledMathProgram *program,
    const semantic::Index kernel_source_view_id,
    const std::size_t source_value_factor_mark) {
  const auto channel_offset = static_cast<semantic::Index>(
      program->source_product_channels.size());
  for (std::size_t factor_pos = source_value_factor_mark;
       factor_pos < program->source_value_factors.size();
       ++factor_pos) {
    auto &factor = program->source_value_factors[factor_pos];
    const auto required_channels =
        compiled_math_source_factor_channel_mask(factor.kind);
    const auto factor_source_view_id =
        factor.source_view_id == semantic::kInvalidIndex ? 0 : factor.source_view_id;
    const auto effective_source_view_id =
        factor_source_view_id == 0 ? kernel_source_view_id : factor_source_view_id;
    semantic::Index channel_id = semantic::kInvalidIndex;
    for (std::size_t channel_pos = static_cast<std::size_t>(channel_offset);
         channel_pos < program->source_product_channels.size();
         ++channel_pos) {
      if (compiled_math_same_source_product_channel(
              program->source_product_channels[channel_pos],
              factor,
              effective_source_view_id)) {
        channel_id = static_cast<semantic::Index>(channel_pos);
        break;
      }
    }
    if (channel_id == semantic::kInvalidIndex) {
      channel_id = static_cast<semantic::Index>(
          program->source_product_channels.size());
      program->source_product_channels.push_back(
          CompiledMathSourceProductChannel{
              factor.source_id,
              factor.condition_id,
              effective_source_view_id,
              factor.time_id,
              factor.time_cap_id,
              required_channels});
    } else {
      program->source_product_channels[
          static_cast<std::size_t>(channel_id)].required_channels |=
          required_channels;
    }
    factor.source_product_channel_id = channel_id;
  }
}

inline bool compiled_math_same_source_factor(
    const CompiledMathSourceValueFactor &lhs,
    const CompiledMathSourceValueFactor &rhs) noexcept {
  return lhs.source_id == rhs.source_id &&
         lhs.condition_id == rhs.condition_id &&
         lhs.source_view_id == rhs.source_view_id &&
         lhs.time_id == rhs.time_id &&
         lhs.time_cap_id == rhs.time_cap_id &&
         lhs.kind == rhs.kind;
}

inline std::vector<CompiledMathSourceValueFactor>
compiled_math_source_factors(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan span) {
  const auto begin = program.source_value_factors.begin() + span.offset;
  return std::vector<CompiledMathSourceValueFactor>(begin, begin + span.size);
}

inline void compiled_math_extract_common_source_factors(
    CompiledMathProgram *program,
    std::vector<CompiledMathSourceProductTerm> *terms,
    const std::size_t factor_mark,
    CompiledMathExecutionPlan *execution) {
  if (terms->size() < 2U) {
    return;
  }
  std::vector<std::vector<CompiledMathSourceValueFactor>> factors;
  factors.reserve(terms->size());
  for (const auto &term : *terms) {
    factors.push_back(compiled_math_source_factors(
        *program, term.source_value_factors));
  }
  auto common = factors.front();
  for (std::size_t term_index = 1U;
       term_index < factors.size() && !common.empty();
       ++term_index) {
    std::vector<std::uint8_t> matched(factors[term_index].size(), 0U);
    auto keep = common.begin();
    for (auto candidate = common.begin(); candidate != common.end(); ++candidate) {
      bool found = false;
      for (std::size_t i = 0; i < factors[term_index].size(); ++i) {
        if (matched[i] == 0U && compiled_math_same_source_factor(
                                    *candidate, factors[term_index][i])) {
          matched[i] = 1U;
          found = true;
          break;
        }
      }
      if (found) {
        *keep++ = *candidate;
      }
    }
    common.erase(keep, common.end());
  }
  if (common.empty()) {
    return;
  }

  program->source_value_factors.resize(factor_mark);
  const auto common_offset = static_cast<semantic::Index>(
      program->source_value_factors.size());
  program->source_value_factors.insert(
      program->source_value_factors.end(), common.begin(), common.end());
  execution->source_value_factors = CompiledMathIndexSpan{
      common_offset, static_cast<semantic::Index>(common.size())};

  for (std::size_t term_index = 0; term_index < terms->size(); ++term_index) {
    auto &remaining = factors[term_index];
    for (const auto &factor : common) {
      const auto found = std::find_if(
          remaining.begin(),
          remaining.end(),
          [&](const CompiledMathSourceValueFactor &candidate) {
            return compiled_math_same_source_factor(candidate, factor);
          });
      if (found != remaining.end()) {
        remaining.erase(found);
      }
    }
    const auto offset = static_cast<semantic::Index>(
        program->source_value_factors.size());
    program->source_value_factors.insert(
        program->source_value_factors.end(), remaining.begin(), remaining.end());
    (*terms)[term_index].source_value_factors = CompiledMathIndexSpan{
        offset, static_cast<semantic::Index>(remaining.size())};
  }
}

inline bool compiled_math_plan_execution(
    CompiledMathProgram *program,
    const semantic::Index node_id,
    const semantic::Index source_view_id,
    const bool initial_state,
    CompiledMathExecutionPlan *execution) {
  std::vector<CompiledMathSourceProductTerm> source_product_terms;
  bool clean_signed_source_sum = false;
  const auto source_value_factor_mark =
      program->source_value_factors.size();
  const auto source_product_channel_mark =
      program->source_product_channels.size();
  const auto outcome_gate_node_mark =
      program->outcome_gate_nodes.size();
  const auto time_gate_node_mark =
      program->time_gate_nodes.size();
  const auto integral_factor_node_mark =
      program->integral_factor_nodes.size();
  const auto expr_upper_factor_mark =
      program->expr_upper_factors.size();
  if (compiled_math_collect_source_product_terms(
          program,
          node_id,
          1.0,
          initial_state,
          &source_product_terms,
          &clean_signed_source_sum)) {
    if (source_product_terms.size() == 1U &&
        source_product_terms.front().sign == 1.0 &&
        source_product_terms.front().outcome_gate_nodes.empty() &&
        source_product_terms.front().time_gate_nodes.empty() &&
        source_product_terms.front().integral_factor_nodes.empty() &&
        source_product_terms.front().expr_upper_factors.empty()) {
      execution->kind = CompiledMathExecutionKind::SourceProduct;
      execution->source_value_factors =
          source_product_terms.front().source_value_factors;
      execution->clean_signed_source_sum = clean_signed_source_sum;
    } else {
      execution->kind = CompiledMathExecutionKind::SourceProductSum;
      compiled_math_extract_common_source_factors(
          program,
          &source_product_terms,
          source_value_factor_mark,
          execution);
      execution->source_product_terms = CompiledMathIndexSpan{
          static_cast<semantic::Index>(
              program->source_product_terms.size()),
          static_cast<semantic::Index>(source_product_terms.size())};
      execution->clean_signed_source_sum = clean_signed_source_sum;
      program->source_product_terms.insert(
          program->source_product_terms.end(),
          source_product_terms.begin(),
          source_product_terms.end());
    }
    compiled_math_compile_source_product_channels(
        program,
        source_view_id,
        source_value_factor_mark);
    return true;
  }
  program->source_value_factors.resize(source_value_factor_mark);
  program->source_product_channels.resize(source_product_channel_mark);
  program->outcome_gate_nodes.resize(outcome_gate_node_mark);
  program->time_gate_nodes.resize(time_gate_node_mark);
  program->integral_factor_nodes.resize(integral_factor_node_mark);
  program->expr_upper_factors.resize(expr_upper_factor_mark);
  *execution = CompiledMathExecutionPlan{};
  return false;
}

inline void compiled_math_plan_root_executions(
    CompiledMathProgram *program) {
  for (auto &root : program->roots) {
    compiled_math_plan_execution(
        program, root.node_id, 0, false, &root.execution);
    compiled_math_plan_execution(
        program, root.node_id, 0, true, &root.initial_execution);
  }
}

inline semantic::Index compiled_math_integral_kernel_slot(
    CompiledMathProgram *program,
    const CompiledMathNode &node) {
  if (!compiled_math_is_integral_node(node.kind)) {
    return semantic::kInvalidIndex;
  }
  const auto root_id = compiled_math_integral_root_id(node);
  if (root_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto root_pos = static_cast<std::size_t>(root_id);
  if (root_pos >= program->roots.size()) {
    throw std::runtime_error(
        "compiled integral node references an unplanned integrand root");
  }
  const auto slot =
      static_cast<semantic::Index>(program->integral_kernels.size());
  CompiledMathIntegralKernel kernel;
  kernel.root_id = root_id;
  kernel.bind_time_id = compiled_math_integral_bind_time_id(node);
  const auto source_view_id =
      node.source_view_id == semantic::kInvalidIndex ? 0 : node.source_view_id;
  compiled_math_plan_execution(
      program,
      program->roots[root_pos].node_id,
      source_view_id,
      false,
      &kernel.execution);
  compiled_math_plan_execution(
      program,
      program->roots[root_pos].node_id,
      source_view_id,
      true,
      &kernel.initial_execution);
  program->integral_kernels.push_back(kernel);
  return slot;
}

} // namespace detail
} // namespace accumulatr::eval
