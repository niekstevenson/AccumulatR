#pragma once

#include "compiled_math_builders.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool probability_integrand_matches(const CompiledMathProgram &program,
    const CompiledMathSourceProductTermBuild &term, const semantic::Index primitive,
    const semantic::Index removed) {
  auto factors = term.source_value_nodes;
  const auto found = std::find(factors.begin(), factors.end(), removed);
  if (found != factors.end()) {
    factors.erase(found);
    if (!term.integral_factor_nodes.empty()) return false;
  } else if (term.integral_factor_nodes.size() != 1 ||
             term.integral_factor_nodes[0] != removed) return false;
  const auto &integral = program.nodes[primitive];
  const auto observed = static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  if (integral.time_id != observed || integral.source_view_id != 0) return false;
  if (integral.kind == CompiledMathNodeKind::SourceCdf) {
    if (factors.size() != 1 || integral.aux_id != semantic::kInvalidIndex) return false;
    const auto &density = program.nodes[factors[0]];
    return density.kind == CompiledMathNodeKind::SourcePdf &&
        density.subject_id == integral.subject_id && density.time_id == observed &&
        density.source_view_id == 0 && density.aux_id == semantic::kInvalidIndex;
  }
  const auto root = program.roots[integral.subject_id].node_id;
  const auto bind = compiled_math_integral_bind_time_id(integral);
  if (bind != observed && compiled_math_depends_on_integral_scope(program, root, observed)) return false;
  std::vector<CompiledMathSourceProductTermBuild> terms;
  bool clean = false;
  if (!compiled_math_expand_source_product_terms(program, root, 1.0, true, &terms, &clean) ||
      terms.size() != 1 || terms[0].sign != 1.0 ||
      !terms[0].integral_factor_nodes.empty() || !terms[0].time_gate_nodes.empty() ||
      !terms[0].expr_upper_factors.empty() ||
      terms[0].source_value_nodes.size() != factors.size()) return false;
  std::vector<bool> matched(factors.size(), false);
  for (const auto id : terms[0].source_value_nodes) {
    const auto &a = program.nodes[id];
    bool found = false;
    for (std::size_t i = 0; i < factors.size(); ++i) {
      const auto &b = program.nodes[factors[i]];
      if (!matched[i] && a.kind == b.kind && a.subject_id == b.subject_id &&
          a.time_id == bind && b.time_id == observed &&
          a.source_view_id == 0 && b.source_view_id == 0 &&
          a.aux_id == semantic::kInvalidIndex && b.aux_id == semantic::kInvalidIndex) {
        matched[i] = found = true; break;
      }
    }
    if (!found) return false;
  }
  return true;
}

// Recognize one ordering integral, with no outside time dependence.
inline std::vector<semantic::Index> probability_pair_sources(
    const CompiledMathProgram &program, const semantic::Index id) {
  const auto &node = program.nodes[id];
  if (!compiled_math_is_integral_node(node.kind) || node.source_view_id != 0 ||
      node.time_id != static_cast<semantic::Index>(CompiledMathTimeSlot::Observed)) return {};
  std::vector<CompiledMathSourceProductTermBuild> terms;
  bool clean = false;
  if (!compiled_math_expand_source_product_terms(program,
          program.roots[node.subject_id].node_id, 1.0, true, &terms, &clean) ||
      terms.size() != 1 || terms[0].sign != 1.0 ||
      terms[0].source_value_nodes.size() != 2 ||
      !terms[0].integral_factor_nodes.empty() || !terms[0].time_gate_nodes.empty() ||
      !terms[0].expr_upper_factors.empty()) return {};
  auto pair = terms[0].source_value_nodes;
  if (program.nodes[pair[1]].kind == CompiledMathNodeKind::SourcePdf) std::swap(pair[0], pair[1]);
  if (program.nodes[pair[0]].kind != CompiledMathNodeKind::SourcePdf ||
      (program.nodes[pair[1]].kind != CompiledMathNodeKind::SourceCdf &&
       program.nodes[pair[1]].kind != CompiledMathNodeKind::SourceSurvival)) return {};
  for (const auto source : pair) {
    const auto &value = program.nodes[source];
    if (value.source_view_id != 0 || value.aux_id != semantic::kInvalidIndex ||
        value.time_id != compiled_math_integral_bind_time_id(node)) return {};
  }
  return pair;
}

inline semantic::Index probability_complement_node(CompiledMathProgram *program,
    const semantic::Index id) {
  const auto node = program->nodes[id];
  if (node.kind == CompiledMathNodeKind::SourceCdf) return compiled_math_source_node(
      program, CompiledMathNodeKind::SourceSurvival, node.subject_id,
      node.time_id, node.source_view_id, node.aux_id);
  const auto pair = probability_pair_sources(*program, id);
  if (!pair.empty()) {
    const auto a = program->nodes[pair[0]].subject_id;
    const auto b = program->nodes[pair[1]].subject_id;
    const auto other = program->nodes[pair[1]].kind;
    const auto sa = compiled_math_source_node(program, CompiledMathNodeKind::SourceSurvival, a);
    const auto sb = compiled_math_source_node(program, CompiledMathNodeKind::SourceSurvival, b);
    const auto fb = compiled_math_source_node(program, CompiledMathNodeKind::SourcePdf, b);
    const auto boundary = other == CompiledMathNodeKind::SourceSurvival
        ? compiled_math_algebra_node(program, CompiledMathNodeKind::Product, {sa, sb})
        : compiled_math_algebra_node(program, CompiledMathNodeKind::Sum,
            {sa,compiled_math_algebra_node(program, CompiledMathNodeKind::Product,
                {compiled_math_source_node(program, CompiledMathNodeKind::SourceCdf, a), sb})});
    const auto integrand = compiled_math_algebra_node(program, CompiledMathNodeKind::Product,
        {fb, compiled_math_source_node(program, other, a)});
    const auto reverse = compiled_math_integral_node(program,
        CompiledMathNodeKind::IntegralZeroToCurrentRaw,
        compiled_math_make_root(program, integrand), node.time_id, 0, node.time_id);
    // Integrate the complementary ordering directly, rather than subtracting
    // a probability near one. The boundary retains source failure probability.
    return compiled_math_algebra_node(program, CompiledMathNodeKind::Sum, {boundary, reverse});
  }
  return compiled_math_unary_node(program, CompiledMathNodeKind::Complement, id);
}

inline semantic::Index probability_sum_complement(CompiledMathProgram *program,
    const std::vector<semantic::Index> &nodes) {
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const auto c = probability_pair_sources(*program, nodes[i]);
    if (c.empty()) continue;
    for (std::size_t k = 0; k < nodes.size(); ++k) {
      if (i == k) continue;
      const auto product = program->nodes[nodes[k]];
      if (product.kind != CompiledMathNodeKind::Product || product.children.size != 2) continue;
      for (std::size_t j = 0; j < 2; ++j) {
        const auto a_id = program->child_nodes[product.children.offset + j];
        const auto b_id = program->child_nodes[product.children.offset + 1 - j];
        const auto a = probability_pair_sources(*program, a_id);
        if (a.empty() || program->nodes[a[0]].subject_id != program->nodes[c[0]].subject_id ||
            program->nodes[a[1]].subject_id != program->nodes[c[1]].subject_id ||
            program->nodes[a[1]].kind == program->nodes[c[1]].kind) continue;
        // C + A = F when their integrands differ only by CDF versus survival.
        // Thus 1 - (C + A B) = (1 - F) + A (1 - B).
        const auto survival = compiled_math_source_node(program, CompiledMathNodeKind::SourceSurvival,
            program->nodes[a[0]].subject_id);
        const auto residual = compiled_math_algebra_node(program, CompiledMathNodeKind::Product,
            {a_id, probability_complement_node(program, b_id)});
        std::vector<semantic::Index> result{survival, residual};
        for (std::size_t r = 0; r < nodes.size(); ++r) if (r != i && r != k)
          result.push_back(compiled_math_unary_node(program, CompiledMathNodeKind::Negate, nodes[r]));
        return compiled_math_algebra_node(program, CompiledMathNodeKind::Sum, std::move(result));
      }
    }
  }
  return probability_complement_node(program, compiled_math_algebra_node(program,
      CompiledMathNodeKind::CleanSignedSum,nodes));
}

// Zero means the expression contains subtraction or an unknown sign.
inline int probability_expression_sign(const CompiledMathProgram &program,
    const semantic::Index id) {
  const auto &node = program.nodes[id];
  if (compiled_math_is_source_value_node(node.kind)) return 1;
  if (node.kind == CompiledMathNodeKind::Constant) return node.constant < 0 ? -1 : 1;
  if (compiled_math_is_integral_node(node.kind))
    return probability_expression_sign(program, program.roots[node.subject_id].node_id);
  if (node.kind == CompiledMathNodeKind::Complement) return 0;
  int sign = 1;
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child = probability_expression_sign(program, program.child_nodes[node.children.offset + i]);
    if (child == 0) return 0;
    if (node.kind == CompiledMathNodeKind::Product) sign *= child;
    else if (i == 0) sign = child;
    else if (sign != child) return 0;
  }
  return node.kind == CompiledMathNodeKind::Negate ? -sign : sign;
}

inline void compile_outcome_terminal_roots(CompiledMathProgram *program,
    const semantic::Index density_root, semantic::Index *probability, semantic::Index *survival) {
  const auto observed = static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto integrate = [&](semantic::Index root) {
    const auto node = program->nodes[program->roots[root].node_id];
    if (node.kind == CompiledMathNodeKind::SourcePdf && node.time_id == observed &&
        node.aux_id == semantic::kInvalidIndex && node.source_view_id == 0)
      return compiled_math_source_node(program, CompiledMathNodeKind::SourceCdf, node.subject_id);
    return compiled_math_integral_node(program, CompiledMathNodeKind::IntegralZeroToCurrentRaw,
                                      root, observed, 0, observed);
  };
  std::vector<CompiledMathSourceProductTermBuild> terms;
  bool clean = false;
  if (!compiled_math_expand_source_product_terms(*program,
          program->roots[density_root].node_id, 1.0, true, &terms, &clean) ||
      std::any_of(terms.begin(), terms.end(), [](const auto &term) {
        return !term.expr_upper_factors.empty();
      })) {
    const auto node = integrate(density_root);
    *probability = compiled_math_make_root(program, node);
    *survival = compiled_math_make_root(program, probability_complement_node(program, node));
    return;
  }
  std::vector<semantic::Index> nodes;
  std::vector<bool> used(terms.size(), false);
  const auto primitives = [&](const CompiledMathSourceProductTermBuild &term) {
    auto ids = term.integral_factor_nodes;
    for (const auto id : term.source_value_nodes)
      if (program->nodes[id].kind == CompiledMathNodeKind::SourceCdf) ids.push_back(id);
    return ids;
  };
  for (std::size_t i = 0; i < terms.size(); ++i) {
    if (used[i]) continue;
    const auto &a = terms[i];
    semantic::Index node = semantic::kInvalidIndex;
    if (a.time_gate_nodes.empty()) {
      for (const auto fa : primitives(a)) {
        for (std::size_t j = i + 1; j < terms.size(); ++j) {
          const auto &b = terms[j];
          if (used[j] || a.sign != b.sign || !b.time_gate_nodes.empty()) continue;
          for (const auto fb : primitives(b)) {
            // At least one integral vanishes at zero, including with shifted sources.
            if (!compiled_math_is_integral_node(program->nodes[fa].kind) &&
                !compiled_math_is_integral_node(program->nodes[fb].kind)) continue;
            if (probability_integrand_matches(*program, a, fb, fa) &&
                probability_integrand_matches(*program, b, fa, fb)) {
              // Integrate f G + F g as F G, including primitive source CDFs.
              node = compiled_math_algebra_node(program, CompiledMathNodeKind::Product, {fa, fb});
              used[j] = true;
              break;
            }
          }
          if (node != semantic::kInvalidIndex) break;
        }
        if (node != semantic::kInvalidIndex) break;
      }
    }
    if (node == semantic::kInvalidIndex) {
      auto factors = a.source_value_nodes;
      factors.insert(factors.end(), a.integral_factor_nodes.begin(), a.integral_factor_nodes.end());
      auto integrand = compiled_math_algebra_node(program,
          CompiledMathNodeKind::Product, std::move(factors));
      for (const auto gate_id : a.time_gate_nodes) {
        const auto gate = program->nodes[gate_id];
        integrand = gate.kind == CompiledMathNodeKind::StrictTimeGate
            ? compiled_math_strict_time_gate_node(program, integrand, gate.time_id, gate.aux_id)
            : compiled_math_time_gate_node(program, integrand, gate.time_id, gate.aux_id);
      }
      node = integrate(compiled_math_make_root(program, integrand));
    }
    if (a.sign != 1.0) node = compiled_math_algebra_node(program,
        CompiledMathNodeKind::Product, {compiled_math_constant(program, a.sign), node});
    nodes.push_back(node);
  }
  *probability = compiled_math_make_root(program,
      compiled_math_algebra_node(program, CompiledMathNodeKind::CleanSignedSum, nodes));
  *survival = compiled_math_make_root(program, probability_sum_complement(program, nodes));
}
} // namespace detail
} // namespace accumulatr::eval
