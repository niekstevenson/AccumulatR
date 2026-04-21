#pragma once

#include <Rcpp.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "layout.hpp"
#include "../compile/project_semantic.hpp"

namespace accumulatr::runtime {

struct ExactProgram {
  RuntimeLayout layout{};

  std::vector<LeafRuntimeDescriptor> leaf_descriptors;
  std::vector<std::uint8_t> leaf_dist_kind;

  std::vector<std::uint8_t> onset_kind;
  std::vector<std::uint8_t> onset_source_kind;
  std::vector<semantic::Index> onset_source_index;
  std::vector<double> onset_lag;
  std::vector<double> onset_abs_value;

  std::vector<semantic::Index> leaf_trigger_index;

  std::vector<std::uint8_t> trigger_kind;
  std::vector<double> trigger_fixed_q;
  std::vector<std::uint8_t> trigger_has_fixed_q;
  std::vector<semantic::Index> trigger_member_offsets;
  std::vector<semantic::Index> trigger_member_indices;

  std::vector<semantic::Index> pool_k;
  std::vector<semantic::Index> pool_member_offsets;
  std::vector<semantic::Index> pool_member_indices;
  std::vector<std::uint8_t> pool_member_kind;

  std::vector<std::uint8_t> expr_kind;
  std::vector<semantic::Index> expr_arg_offsets;
  std::vector<semantic::Index> expr_args;
  std::vector<semantic::Index> expr_ref_child;
  std::vector<semantic::Index> expr_blocker_child;
  std::vector<semantic::Index> expr_source_index;
  std::vector<std::uint8_t> expr_source_kind;
  std::vector<int> expr_event_k;

  std::vector<semantic::Index> outcome_expr_root;
  std::vector<semantic::Index> observed_label_index;

  ParameterLayout parameter_layout{};
};

struct LoweredExactVariant {
  std::string component_id;
  double weight{1.0};
  std::string weight_name;
  ExactProgram program{};
  std::vector<std::string> param_keys;
  std::vector<std::string> leaf_ids;
  std::vector<std::string> pool_ids;
  std::vector<std::string> outcome_labels;
};

struct LoweredExactModel {
  std::string component_mode{"fixed"};
  std::string component_reference;
  std::vector<LoweredExactVariant> variants;
};

namespace detail {

inline void require_exact_variant(const compile::CompiledVariant &variant) {
  if (variant.backend == compile::BackendKind::Exact) {
    return;
  }
  std::string message = "cannot lower direct variant";
  if (!variant.component_id.empty()) {
    message += " '" + variant.component_id + "'";
  }
  throw std::runtime_error(message);
}

template <typename T>
Rcpp::IntegerVector as_integer_vector_exact(const std::vector<T> &values) {
  Rcpp::IntegerVector out(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = static_cast<int>(values[i]);
  }
  return out;
}

inline Rcpp::NumericVector as_numeric_vector_exact(
    const std::vector<double> &values) {
  Rcpp::NumericVector out(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = values[i];
  }
  return out;
}

class ExactSlotAllocator {
public:
  semantic::Index slot_for(const std::string &key) {
    if (key.empty()) {
      return semantic::kInvalidIndex;
    }
    const auto it = slot_by_key_.find(key);
    if (it != slot_by_key_.end()) {
      return it->second;
    }
    const auto slot = static_cast<semantic::Index>(keys_.size());
    keys_.push_back(key);
    slot_by_key_.emplace(key, slot);
    return slot;
  }

  const std::vector<std::string> &keys() const noexcept {
    return keys_;
  }

private:
  std::unordered_map<std::string, semantic::Index> slot_by_key_;
  std::vector<std::string> keys_;
};

inline Rcpp::List to_r_list(const ExactProgram &program) {
  return Rcpp::List::create(
      Rcpp::Named("layout") = Rcpp::List::create(
          Rcpp::Named("n_leaves") = program.layout.n_leaves,
          Rcpp::Named("n_pools") = program.layout.n_pools,
          Rcpp::Named("n_outcomes") = program.layout.n_outcomes,
          Rcpp::Named("n_params") = program.layout.n_params,
          Rcpp::Named("n_triggers") = program.layout.n_triggers),
      Rcpp::Named("leaf_dist_kind") =
          as_integer_vector_exact(program.leaf_dist_kind),
      Rcpp::Named("onset_kind") = as_integer_vector_exact(program.onset_kind),
      Rcpp::Named("onset_source_kind") =
          as_integer_vector_exact(program.onset_source_kind),
      Rcpp::Named("onset_source_index") =
          as_integer_vector_exact(program.onset_source_index),
      Rcpp::Named("onset_lag") = as_numeric_vector_exact(program.onset_lag),
      Rcpp::Named("onset_abs_value") =
          as_numeric_vector_exact(program.onset_abs_value),
      Rcpp::Named("leaf_trigger_index") =
          as_integer_vector_exact(program.leaf_trigger_index),
      Rcpp::Named("trigger_kind") =
          as_integer_vector_exact(program.trigger_kind),
      Rcpp::Named("trigger_fixed_q") =
          as_numeric_vector_exact(program.trigger_fixed_q),
      Rcpp::Named("trigger_has_fixed_q") =
          as_integer_vector_exact(program.trigger_has_fixed_q),
      Rcpp::Named("trigger_member_offsets") =
          as_integer_vector_exact(program.trigger_member_offsets),
      Rcpp::Named("trigger_member_indices") =
          as_integer_vector_exact(program.trigger_member_indices),
      Rcpp::Named("pool_k") = as_integer_vector_exact(program.pool_k),
      Rcpp::Named("pool_member_offsets") =
          as_integer_vector_exact(program.pool_member_offsets),
      Rcpp::Named("pool_member_indices") =
          as_integer_vector_exact(program.pool_member_indices),
      Rcpp::Named("pool_member_kind") =
          as_integer_vector_exact(program.pool_member_kind),
      Rcpp::Named("expr_kind") = as_integer_vector_exact(program.expr_kind),
      Rcpp::Named("expr_arg_offsets") =
          as_integer_vector_exact(program.expr_arg_offsets),
      Rcpp::Named("expr_args") = as_integer_vector_exact(program.expr_args),
      Rcpp::Named("expr_ref_child") =
          as_integer_vector_exact(program.expr_ref_child),
      Rcpp::Named("expr_blocker_child") =
          as_integer_vector_exact(program.expr_blocker_child),
      Rcpp::Named("expr_source_index") =
          as_integer_vector_exact(program.expr_source_index),
      Rcpp::Named("expr_source_kind") =
          as_integer_vector_exact(program.expr_source_kind),
      Rcpp::Named("expr_event_k") =
          as_integer_vector_exact(program.expr_event_k),
      Rcpp::Named("outcome_expr_root") =
          as_integer_vector_exact(program.outcome_expr_root),
      Rcpp::Named("observed_label_index") =
          as_integer_vector_exact(program.observed_label_index),
      Rcpp::Named("parameter_layout") = Rcpp::List::create(
          Rcpp::Named("leaf_param_offsets") =
              as_integer_vector_exact(program.parameter_layout.leaf_param_offsets),
          Rcpp::Named("leaf_param_slots") =
              as_integer_vector_exact(program.parameter_layout.leaf_param_slots),
          Rcpp::Named("leaf_q_slots") =
              as_integer_vector_exact(program.parameter_layout.leaf_q_slots),
          Rcpp::Named("leaf_t0_slots") =
              as_integer_vector_exact(program.parameter_layout.leaf_t0_slots)));
}

} // namespace detail

inline LoweredExactVariant lower_exact_variant(
    const compile::CompiledVariant &variant) {
  detail::require_exact_variant(variant);

  LoweredExactVariant lowered;
  lowered.component_id = variant.component_id;
  lowered.weight = variant.weight;
  lowered.weight_name = variant.weight_name;

  const auto &model = variant.model;
  auto &program = lowered.program;

  program.layout.n_leaves = static_cast<int>(model.leaves.size());
  program.layout.n_pools = static_cast<int>(model.pools.size());
  program.layout.n_outcomes = static_cast<int>(model.outcomes.size());
  program.layout.n_triggers = static_cast<int>(model.triggers.size());

  lowered.leaf_ids.reserve(model.leaves.size());
  lowered.pool_ids.reserve(model.pools.size());
  lowered.outcome_labels.reserve(model.outcomes.size());

  program.leaf_dist_kind.reserve(model.leaves.size());
  program.leaf_descriptors.reserve(model.leaves.size());
  program.onset_kind.reserve(model.leaves.size());
  program.onset_source_kind.reserve(model.leaves.size());
  program.onset_source_index.reserve(model.leaves.size());
  program.onset_lag.reserve(model.leaves.size());
  program.onset_abs_value.reserve(model.leaves.size());
  program.leaf_trigger_index.reserve(model.leaves.size());
  program.parameter_layout.leaf_param_offsets.reserve(model.leaves.size() + 1U);
  program.parameter_layout.leaf_q_slots.reserve(model.leaves.size());
  program.parameter_layout.leaf_t0_slots.reserve(model.leaves.size());

  detail::ExactSlotAllocator slots;
  program.parameter_layout.leaf_param_offsets.push_back(0);
  for (const auto &leaf : model.leaves) {
    const auto param_offset = static_cast<semantic::Index>(
        program.parameter_layout.leaf_param_slots.size());
    lowered.leaf_ids.push_back(leaf.id);
    program.leaf_dist_kind.push_back(static_cast<std::uint8_t>(leaf.dist));
    program.onset_kind.push_back(static_cast<std::uint8_t>(leaf.onset.kind));
    program.onset_source_kind.push_back(
        static_cast<std::uint8_t>(leaf.onset.source.kind));
    program.onset_source_index.push_back(leaf.onset.source.index);
    program.onset_lag.push_back(leaf.onset.lag);
    program.onset_abs_value.push_back(leaf.onset.absolute_value);
    program.leaf_trigger_index.push_back(leaf.trigger_index);

    for (const auto &param_name : leaf.params.dist_param_names) {
      program.parameter_layout.leaf_param_slots.push_back(
          slots.slot_for(param_name));
    }
    const auto param_end = static_cast<semantic::Index>(
        program.parameter_layout.leaf_param_slots.size());
    program.parameter_layout.leaf_param_offsets.push_back(param_end);
    program.parameter_layout.leaf_q_slots.push_back(
        slots.slot_for(leaf.params.q_name));
    program.parameter_layout.leaf_t0_slots.push_back(
        slots.slot_for(leaf.params.t0_name));
    program.leaf_descriptors.push_back(LeafRuntimeDescriptor{
        static_cast<std::uint8_t>(leaf.dist),
        static_cast<std::uint8_t>(leaf.onset.kind),
        static_cast<std::uint8_t>(leaf.onset.source.kind),
        leaf.onset.source.index,
        leaf.onset.lag,
        leaf.onset.absolute_value,
        leaf.trigger_index,
        param_offset,
        static_cast<int>(param_end - param_offset)});
  }

  program.trigger_kind.reserve(model.triggers.size());
  program.trigger_fixed_q.reserve(model.triggers.size());
  program.trigger_has_fixed_q.reserve(model.triggers.size());
  program.trigger_member_offsets.reserve(model.triggers.size() + 1U);
  program.trigger_member_offsets.push_back(0);
  for (const auto &trigger : model.triggers) {
    program.trigger_kind.push_back(static_cast<std::uint8_t>(trigger.kind));
    program.trigger_fixed_q.push_back(trigger.fixed_q);
    program.trigger_has_fixed_q.push_back(
        trigger.has_fixed_q ? static_cast<std::uint8_t>(1)
                            : static_cast<std::uint8_t>(0));
    for (const auto leaf_index : trigger.leaf_indices) {
      program.trigger_member_indices.push_back(leaf_index);
    }
    program.trigger_member_offsets.push_back(
        static_cast<semantic::Index>(program.trigger_member_indices.size()));
  }

  program.pool_k.reserve(model.pools.size());
  program.pool_member_offsets.reserve(model.pools.size() + 1U);
  program.pool_member_offsets.push_back(0);
  for (const auto &pool : model.pools) {
    lowered.pool_ids.push_back(pool.id);
    program.pool_k.push_back(pool.k);
    for (const auto &member : pool.members) {
      program.pool_member_indices.push_back(member.index);
      program.pool_member_kind.push_back(
          static_cast<std::uint8_t>(member.kind));
    }
    program.pool_member_offsets.push_back(
        static_cast<semantic::Index>(program.pool_member_indices.size()));
  }

  program.expr_kind.reserve(model.expr_nodes.size());
  program.expr_arg_offsets.reserve(model.expr_nodes.size() + 1U);
  program.expr_ref_child.reserve(model.expr_nodes.size());
  program.expr_blocker_child.reserve(model.expr_nodes.size());
  program.expr_source_index.reserve(model.expr_nodes.size());
  program.expr_source_kind.reserve(model.expr_nodes.size());
  program.expr_event_k.reserve(model.expr_nodes.size());
  program.expr_arg_offsets.push_back(0);
  for (const auto &expr : model.expr_nodes) {
    program.expr_kind.push_back(static_cast<std::uint8_t>(expr.kind));
    for (const auto child : expr.children) {
      program.expr_args.push_back(child);
    }
    for (const auto child : expr.unless_children) {
      program.expr_args.push_back(child);
    }
    program.expr_arg_offsets.push_back(
        static_cast<semantic::Index>(program.expr_args.size()));
    program.expr_ref_child.push_back(expr.reference_child);
    program.expr_blocker_child.push_back(expr.blocker_child);
    program.expr_source_index.push_back(expr.source.index);
    program.expr_source_kind.push_back(
        static_cast<std::uint8_t>(expr.source.kind));
    program.expr_event_k.push_back(expr.event_k);
  }

  program.outcome_expr_root.reserve(model.outcomes.size());
  program.observed_label_index.reserve(model.outcomes.size());
  for (std::size_t i = 0; i < model.outcomes.size(); ++i) {
    const auto &outcome = model.outcomes[i];
    lowered.outcome_labels.push_back(outcome.label);
    program.outcome_expr_root.push_back(outcome.expr_root);
    program.observed_label_index.push_back(static_cast<semantic::Index>(i));
  }

  lowered.param_keys = slots.keys();
  program.layout.n_params = static_cast<int>(lowered.param_keys.size());

  return lowered;
}

inline LoweredExactModel lower_exact_model(
    const compile::CompiledModel &compiled) {
  LoweredExactModel lowered;
  lowered.component_mode = compiled.component_mode;
  lowered.component_reference = compiled.component_reference;
  lowered.variants.reserve(compiled.variants.size());
  for (const auto &variant : compiled.variants) {
    lowered.variants.push_back(lower_exact_variant(variant));
  }
  return lowered;
}

inline Rcpp::List to_r_list(const LoweredExactModel &lowered) {
  Rcpp::List variants(lowered.variants.size());
  for (std::size_t i = 0; i < lowered.variants.size(); ++i) {
    const auto &variant = lowered.variants[i];
    variants[i] = Rcpp::List::create(
        Rcpp::Named("component_id") = variant.component_id,
        Rcpp::Named("weight") = variant.weight,
        Rcpp::Named("weight_name") = variant.weight_name,
        Rcpp::Named("param_keys") = variant.param_keys,
        Rcpp::Named("leaf_ids") = variant.leaf_ids,
        Rcpp::Named("pool_ids") = variant.pool_ids,
        Rcpp::Named("outcome_labels") = variant.outcome_labels,
        Rcpp::Named("program") = detail::to_r_list(variant.program));
  }

  return Rcpp::List::create(
      Rcpp::Named("component_mode") = lowered.component_mode,
      Rcpp::Named("component_reference") = lowered.component_reference,
      Rcpp::Named("variants") = variants);
}

} // namespace accumulatr::runtime
