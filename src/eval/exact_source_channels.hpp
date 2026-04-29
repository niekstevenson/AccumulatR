#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "eval_query.hpp"
#include "exact_planner.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactLoadedLeafInput {
  std::array<double, 8> params{};
  double q{0.0};
  double t0{0.0};
};

class CompiledSourceChannels {
public:
  struct SourceProductScalarFill {
    std::uint8_t mask{0U};
    double pdf{0.0};
    double cdf{0.0};
    double survival{1.0};
  };

  struct ResolvedSourceBounds {
    const double *exact_time{nullptr};
    double lower{0.0};
    double upper{std::numeric_limits<double>::infinity()};
    bool has_condition_lower{false};
  };

  explicit CompiledSourceChannels(const ExactVariantPlan &plan)
      : plan_(plan),
        program_(plan.lowered.program),
        leaf_inputs_(static_cast<std::size_t>(program_.layout.n_leaves)),
        source_product_direct_available_(
            plan.compiled_math.integral_kernel_source_product_channels.size(),
            0U) {}

  void reset(const ParamView &params,
             const int first_param_row,
             const ExactTriggerState &trigger_state,
             const ExactSequenceState &sequence_state,
             const double top_time) {
    params_ = &params;
    first_param_row_ = first_param_row;
    trigger_state_ = &trigger_state;
    sequence_state_ = &sequence_state;
    conditional_time_ = top_time;
    for (int i = 0; i < program_.layout.n_leaves; ++i) {
      const auto pos = static_cast<std::size_t>(i);
      const auto &desc = program_.leaf_descriptors[pos];
      const int row = first_param_row_ + i;
      auto &loaded = leaf_inputs_[pos];
      const int n_local = std::min<int>(desc.param_count, 8);
      for (int j = 0; j < n_local; ++j) {
        loaded.params[static_cast<std::size_t>(j)] = params_->p(row, j);
      }
      for (int j = n_local; j < 8; ++j) {
        loaded.params[static_cast<std::size_t>(j)] = 0.0;
      }
      loaded.q = leaf_q(i, row);
      loaded.t0 = params_->t0(row);
    }
    const auto &source_product_channels =
        plan_.compiled_math.integral_kernel_source_product_channels;
    if (source_product_direct_available_.size() <
        source_product_channels.size()) {
      source_product_direct_available_.resize(
          source_product_channels.size(),
          0U);
    }
    for (std::size_t i = 0; i < source_product_channels.size(); ++i) {
      const auto &channel = source_product_channels[i];
      source_product_direct_available_[i] =
          channel.direct_leaf_absolute_candidate &&
                  sequence_bounds_inactive_for(channel.source_id, channel.bounds)
              ? 1U
              : 0U;
    }
  }

  bool source_product_direct_leaf_available(
      const semantic::Index source_product_channel_id) const {
    const auto pos = static_cast<std::size_t>(source_product_channel_id);
    return pos < source_product_direct_available_.size() &&
           source_product_direct_available_[pos] != 0U;
  }

  const ExactLoadedLeafInput &source_product_leaf_input(
      const semantic::Index leaf_index) const {
    return leaf_inputs_[static_cast<std::size_t>(leaf_index)];
  }

  const double *exact_time_for_source(
      const semantic::Index source_id,
      const semantic::Index condition_id = 0,
      const CompiledMathWorkspace *workspace = nullptr) const {
    return exact_time_for(source_id, condition_id, workspace);
  }

  bool source_order_known_before(
      const semantic::Index before_source_id,
      const semantic::Index after_source_id,
      const semantic::Index condition_id = 0,
      const CompiledMathWorkspace *workspace = nullptr) const {
    if (before_source_id == semantic::kInvalidIndex ||
        after_source_id == semantic::kInvalidIndex ||
        before_source_id == after_source_id) {
      return false;
    }
    return condition_source_order_known_before(
        before_source_id, after_source_id, condition_id, workspace);
  }

  bool expr_upper_bound_for(
      const semantic::Index expr_id,
      ExactTimedExprUpperBound *out) const {
    if (expr_id == semantic::kInvalidIndex || sequence_state_ == nullptr) {
      return false;
    }
    const auto pos = static_cast<std::size_t>(expr_id);
    if (pos >= sequence_state_->expr_upper_bounds.size() ||
        pos >= sequence_state_->expr_upper_normalizers.size()) {
      return false;
    }
    const double time = sequence_state_->expr_upper_bounds[pos];
    const double normalizer = sequence_state_->expr_upper_normalizers[pos];
    if (!std::isfinite(time) || !(normalizer > 0.0)) {
      return false;
    }
    if (out != nullptr) {
      out->expr_id = expr_id;
      out->time = time;
      out->normalizer = normalizer;
    }
    return true;
  }

  ResolvedSourceBounds source_product_resolved_bounds(
      const semantic::Index source_id,
      const CompiledSourceBoundPlan &bounds,
      const CompiledMathWorkspace *workspace) const {
    return compiled_source_bounds(source_id, bounds, workspace);
  }

private:
  const ExactVariantPlan &plan_;
  const runtime::ExactProgram &program_;
  const ParamView *params_{nullptr};
  int first_param_row_;
  const ExactTriggerState *trigger_state_{nullptr};
  const ExactSequenceState *sequence_state_{nullptr};
  double conditional_time_;
  std::vector<ExactLoadedLeafInput> leaf_inputs_;
  std::vector<std::uint8_t> source_product_direct_available_;
  mutable double compiled_condition_time_{
      std::numeric_limits<double>::quiet_NaN()};

  CompiledSourceBoundPlan empty_source_bound_plan_{};

  bool sequence_bounds_inactive_for(
      const semantic::Index source_id,
      const CompiledSourceBoundPlan &bounds) const {
    if (sequence_state_ == nullptr) {
      return true;
    }
    if (bounds.use_sequence_lower && sequence_state_->lower_bound > 0.0) {
      return false;
    }
    const auto source_pos = static_cast<std::size_t>(source_id);
    if (bounds.use_sequence_exact &&
        source_pos < sequence_state_->exact_times.size() &&
        std::isfinite(sequence_state_->exact_times[source_pos])) {
      return false;
    }
    if (bounds.use_sequence_upper &&
        source_pos < sequence_state_->upper_bounds.size() &&
        std::isfinite(sequence_state_->upper_bounds[source_pos])) {
      return false;
    }
    return true;
  }

  const CompiledSourceBoundPlan &source_bound_plan(
      const semantic::Index condition_id,
      const semantic::Index source_id) const {
    const auto slot =
        compiled_source_bound_plan_slot(
            plan_.compiled_math, condition_id, source_id);
    if (slot == semantic::kInvalidIndex) {
      return empty_source_bound_plan_;
    }
    return plan_.compiled_math.source_condition_bound_plans[
        static_cast<std::size_t>(slot)];
  }

  double compiled_bound_term_time(
      const CompiledMathWorkspace *workspace,
      const CompiledSourceBoundTerm &term) const {
    const auto time_id = term.time_id;
    if (workspace != nullptr && workspace->has_time(time_id)) {
      return workspace->time(time_id);
    }
    return conditional_time_;
  }

  const double *exact_time_for(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    return compiled_source_bounds(
        source_id,
        source_bound_plan(condition_id, source_id),
        workspace).exact_time;
  }

  ResolvedSourceBounds compiled_source_bounds(
      const semantic::Index source_id,
      const CompiledSourceBoundPlan &bounds,
      const CompiledMathWorkspace *workspace) const {
    ResolvedSourceBounds resolved;
    if (bounds.use_sequence_lower && sequence_state_ != nullptr) {
      resolved.lower = sequence_state_->lower_bound;
    }
    if (source_id == semantic::kInvalidIndex) {
      return resolved;
    }
    const auto source_pos = static_cast<std::size_t>(source_id);
    const bool sequence_exact =
        bounds.use_sequence_exact &&
        sequence_state_ != nullptr &&
        source_pos < sequence_state_->exact_times.size() &&
        std::isfinite(sequence_state_->exact_times[source_pos]);
    const bool sequence_upper =
        bounds.use_sequence_upper &&
        sequence_state_ != nullptr &&
        source_pos < sequence_state_->upper_bounds.size() &&
        std::isfinite(sequence_state_->upper_bounds[source_pos]);
    if (sequence_exact || sequence_upper) {
      resolved.lower = 0.0;
    }
    if (sequence_upper) {
      resolved.upper = sequence_state_->upper_bounds[source_pos];
    }

    if (bounds.has_condition_lower) {
      for (semantic::Index i = 0; i < bounds.lower.size; ++i) {
        const auto term_pos =
            static_cast<std::size_t>(bounds.lower.offset + i);
        resolved.lower = std::max(
            resolved.lower,
            compiled_bound_term_time(
                workspace,
                plan_.compiled_math.source_condition_bound_terms[term_pos]));
      }
      resolved.has_condition_lower = true;
    }
    if (bounds.has_condition_upper) {
      for (semantic::Index i = 0; i < bounds.upper.size; ++i) {
        const auto term_pos =
            static_cast<std::size_t>(bounds.upper.offset + i);
        resolved.upper = std::min(
            resolved.upper,
            compiled_bound_term_time(
                workspace,
                plan_.compiled_math.source_condition_bound_terms[term_pos]));
      }
    }
    if (bounds.has_condition_exact) {
      compiled_condition_time_ =
          compiled_bound_term_time(
              workspace,
              plan_.compiled_math.source_condition_bound_terms[
                  static_cast<std::size_t>(bounds.exact.offset)]);
      resolved.exact_time = &compiled_condition_time_;
      return resolved;
    }
    if (sequence_exact) {
      resolved.exact_time = &sequence_state_->exact_times[source_pos];
    }
    return resolved;
  }

  bool condition_source_order_known_before(
      const semantic::Index before_source_id,
      const semantic::Index after_source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    (void)workspace;
    if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
      return false;
    }
    const auto source_count = static_cast<std::size_t>(
        plan_.compiled_math.condition_source_relation_source_count);
    const auto condition_slot = static_cast<std::size_t>(condition_id);
    const auto before_pos = static_cast<std::size_t>(before_source_id);
    const auto after_pos = static_cast<std::size_t>(after_source_id);
    if (source_count == 0U ||
        before_pos >= source_count ||
        after_pos >= source_count) {
      return false;
    }
    const auto relation_offset = condition_slot * source_count;
    if (relation_offset + after_pos >=
        plan_.compiled_math.condition_source_relations.size()) {
      return false;
    }
    const auto before_relation =
        static_cast<ExactRelation>(
            plan_.compiled_math.condition_source_relations[
                relation_offset + before_pos]);
    const auto after_relation =
        static_cast<ExactRelation>(
            plan_.compiled_math.condition_source_relations[
                relation_offset + after_pos]);
    if (before_relation != ExactRelation::Unknown &&
        after_relation != ExactRelation::Unknown &&
        before_relation < after_relation) {
      return true;
    }
    return false;
  }

  double leaf_q(const semantic::Index leaf_index, const int row) const {
    const auto trigger_index =
        program_.leaf_trigger_index[static_cast<std::size_t>(leaf_index)];
    if (trigger_index != semantic::kInvalidIndex &&
        static_cast<semantic::TriggerKind>(
            program_.trigger_kind[static_cast<std::size_t>(trigger_index)]) ==
            semantic::TriggerKind::Shared &&
        trigger_state_->shared_started[static_cast<std::size_t>(trigger_index)] <= 1U) {
      return trigger_state_->shared_started[static_cast<std::size_t>(trigger_index)] == 1U
                 ? 0.0
                 : 1.0;
    }
    return params_->q(row);
  }
};

using ExactSourceChannels = CompiledSourceChannels;

inline double exact_compiled_trigger_state_weight(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactCompiledTriggerState &compiled_state) {
  double weight = compiled_state.fixed_weight;
  const auto &table = plan.trigger_state_table;
  for (semantic::Index i = 0; i < compiled_state.weight_terms.size; ++i) {
    const auto &term =
        table.weight_terms[
            static_cast<std::size_t>(
                compiled_state.weight_terms.offset + i)];
    const double q =
        clamp_probability(params.q(first_param_row + term.leaf_index));
    weight *= term.shared_started == 0U ? q : (1.0 - q);
    if (!(weight > 0.0)) {
      return 0.0;
    }
  }
  return weight;
}

inline ExactTriggerState exact_compiled_trigger_state_view(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactCompiledTriggerState &compiled_state) {
  const auto weight =
      exact_compiled_trigger_state_weight(
          plan, params, first_param_row, compiled_state);
  const auto shared_offset =
      static_cast<std::size_t>(compiled_state.shared_started_offset);
  return ExactTriggerState{
      weight,
      shared_offset < plan.trigger_state_table.shared_started_values.size()
          ? plan.trigger_state_table.shared_started_values.data() + shared_offset
          : nullptr};
}

} // namespace detail
} // namespace accumulatr::eval
