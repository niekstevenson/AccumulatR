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

inline leaf::EventChannels forced_channels(const ExactRelation relation) {
  switch (relation) {
  case ExactRelation::Before:
  case ExactRelation::At:
    return leaf::EventChannels::certain();
  case ExactRelation::After:
    return leaf::EventChannels::impossible();
  case ExactRelation::Unknown:
    break;
  }
  return leaf::EventChannels::impossible();
}

class CompiledSourceChannels {
public:
  struct TimedCache {
    explicit TimedCache(const std::size_t size = 0U)
        : channels(size), epoch(size, 0U) {}

    void invalidate() {
      ++current_epoch;
      if (current_epoch == 0U) {
        current_epoch = 1U;
        std::fill(epoch.begin(), epoch.end(), 0U);
      }
    }

    void reset(const double t) {
      time = t;
      condition_id = 0;
      invalidate();
    }

    void set_time(const double t) {
      set_context(t, 0);
    }

    void set_context(const double t, const semantic::Index frame_id) {
      if (std::fabs(t - time) <= 1e-12 && condition_id == frame_id) {
        return;
      }
      time = t;
      condition_id = frame_id;
      invalidate();
    }

    double time{0.0};
    semantic::Index condition_id{0};
    std::vector<leaf::EventChannels> channels;
    std::vector<std::uint32_t> epoch;
    std::uint32_t current_epoch{1U};
  };

  class BaseTimeGuard {
  public:
    BaseTimeGuard(CompiledSourceChannels *channels,
                  semantic::Index recursive_index,
                  TimedCache *cache,
                  const double time)
        : channels_(channels), recursive_index(recursive_index) {
      if (channels_ == nullptr || cache == nullptr) {
        return;
      }
      channels_->push_base_cache(cache, time);
      active_ = true;
    }

    BaseTimeGuard(const BaseTimeGuard &) = delete;
    BaseTimeGuard &operator=(const BaseTimeGuard &) = delete;

    ~BaseTimeGuard() {
      if (channels_ != nullptr && active_) {
        channels_->pop_base_cache(recursive_index);
      }
    }

  private:
    CompiledSourceChannels *channels_{nullptr};
    semantic::Index recursive_index{-1};
    bool active_{false};
  };


  explicit CompiledSourceChannels(const ExactVariantPlan &plan)
      : plan_(plan),
        program_(plan.lowered.program),
        leaf_inputs_(static_cast<std::size_t>(program_.layout.n_leaves)),
        conditional_cache_(static_cast<std::size_t>(plan.source_count)),
        base_current_cache_(static_cast<std::size_t>(plan.source_count)),
        base_lower_cache_(static_cast<std::size_t>(plan.source_count)) {
    recursive_base_caches_.reserve(
        static_cast<std::size_t>(std::max(program_.layout.n_leaves, 1)));
    base_cache_stack_.reserve(
        static_cast<std::size_t>(std::max(program_.layout.n_leaves, 1)));
    build_source_kernel_evaluators();
  }

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
    recursive_cache_depth_ = 0;
    base_cache_stack_.clear();
    conditional_cache_.reset(top_time);
    base_current_cache_.reset(top_time);
    base_lower_cache_.reset(sequence_state.lower_bound);
    for (auto &cache : recursive_base_caches_) {
      cache.invalidate();
    }
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
  }

  [[nodiscard]] BaseTimeGuard base_time_guard(const double time) {
    const auto recursive_index = acquire_recursive_cache();
    return BaseTimeGuard(
        this,
        recursive_index,
        &recursive_base_caches_[static_cast<std::size_t>(recursive_index)],
        time);
  }

  leaf::EventChannels evaluate_source_channel_plan_at(
      const CompiledSourceChannelPlan &channel_plan,
      const double time,
      const CompiledMathWorkspace *workspace = nullptr) {
    if (can_evaluate_direct_leaf_absolute(channel_plan)) {
      return evaluate_direct_leaf_absolute_at(channel_plan, time, workspace);
    }
    const auto &request = channel_plan.request;
    conditional_cache_.set_context(
        time, compiled_condition_cache_id(request.condition_id, workspace));
    return load_conditional_source(
        &conditional_cache_, time, channel_plan, workspace);
  }

  bool source_channel_condition_invariant(
      const CompiledSourceChannelPlan &channel_plan) const {
    return channel_plan.kernel == CompiledSourceChannelKernelKind::LeafAbsolute &&
           !channel_plan.has_source_condition_overlay;
  }

  bool has_exact_time_overlay(
      const semantic::Index source_id,
      const semantic::Index condition_id = 0,
      const CompiledMathWorkspace *workspace = nullptr) const {
    if (source_id == semantic::kInvalidIndex) {
      return false;
    }
    return has_condition_exact_time(source_id, condition_id, workspace);
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

  bool has_guard_upper_bound_overlay(
      const semantic::Index condition_id = 0) const {
    return condition_has_fact(
        condition_id, CompiledMathConditionFactKind::GuardUpperBound);
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

  bool has_source_order_overlay(
      const semantic::Index condition_id = 0) const {
    return condition_has_fact(
        condition_id, CompiledMathConditionFactKind::SourceOrder);
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
  TimedCache conditional_cache_;
  TimedCache base_current_cache_;
  TimedCache base_lower_cache_;
  std::vector<TimedCache> recursive_base_caches_;
  std::vector<TimedCache *> base_cache_stack_;
  using SourceKernelEvaluator = leaf::EventChannels (CompiledSourceChannels::*)(
      const ExactSourceKernel &,
      semantic::Index,
      const CompiledMathWorkspace *);
  std::vector<SourceKernelEvaluator> source_kernel_evaluators_;
  mutable double compiled_condition_time_{
      std::numeric_limits<double>::quiet_NaN()};
  semantic::Index recursive_cache_depth_{0};

  CompiledSourceBoundPlan empty_source_bound_plan_{};

  struct ResolvedSourceBounds {
    const double *exact_time{nullptr};
    double lower{0.0};
    double upper{std::numeric_limits<double>::infinity()};
    bool has_condition_lower{false};
  };

  static void hash_condition_part(std::size_t *seed,
                                  const std::size_t value) noexcept {
    *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
  }

  static semantic::Index compiled_condition_cache_id(
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) {
    if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
      return 0;
    }
    std::size_t seed = 0x517cc1b727220a95ULL;
    hash_condition_part(&seed, static_cast<std::size_t>(condition_id));
    if (workspace != nullptr) {
      const auto *time_dependencies =
          workspace->condition_time_dependencies(condition_id);
      if (time_dependencies != nullptr) {
        for (const auto time_id : *time_dependencies) {
          const auto time_pos = static_cast<std::size_t>(time_id);
          if (time_pos >= workspace->time_valid.size() ||
              workspace->time_valid[time_pos] == 0U) {
            continue;
          }
          hash_condition_part(
              &seed, static_cast<std::size_t>(time_id) + 1U);
          hash_condition_part(
              &seed,
              std::hash<double>{}(workspace->time_values[time_pos]));
        }
      } else {
        for (std::size_t i = 0; i < workspace->time_values.size(); ++i) {
          if (i >= workspace->time_valid.size() ||
              workspace->time_valid[i] == 0U) {
            continue;
          }
          hash_condition_part(&seed, i + 1U);
          hash_condition_part(
              &seed,
              std::hash<double>{}(workspace->time_values[i]));
        }
      }
    }
    return static_cast<semantic::Index>(
        (seed & static_cast<std::size_t>(0x3fffffffU)) + 1U);
  }

  void build_source_kernel_evaluators() {
    source_kernel_evaluators_.assign(
        plan_.source_kernels.size(),
        nullptr);
    for (std::size_t i = 0; i < plan_.source_kernels.size(); ++i) {
      switch (plan_.source_kernels[i].kind) {
      case CompiledSourceChannelKernelKind::LeafAbsolute:
        source_kernel_evaluators_[i] =
            &CompiledSourceChannels::evaluate_leaf_absolute_kernel;
        break;
      case CompiledSourceChannelKernelKind::LeafOnsetConvolution:
        source_kernel_evaluators_[i] =
            &CompiledSourceChannels::evaluate_leaf_onset_convolution_kernel;
        break;
      case CompiledSourceChannelKernelKind::PoolKOfN:
        source_kernel_evaluators_[i] =
            &CompiledSourceChannels::evaluate_pool_kofn_kernel;
        break;
      case CompiledSourceChannelKernelKind::Invalid:
        break;
      }
    }
  }

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

  bool can_evaluate_direct_leaf_absolute(
      const CompiledSourceChannelPlan &channel_plan) const {
    const auto &request = channel_plan.request;
    if (request.source_id == semantic::kInvalidIndex ||
        channel_plan.source_kernel_slot == semantic::kInvalidIndex ||
        !source_channel_condition_invariant(channel_plan)) {
      return false;
    }
    return sequence_bounds_inactive_for(request.source_id, channel_plan.bounds);
  }

  leaf::EventChannels evaluate_direct_leaf_absolute_at(
      const CompiledSourceChannelPlan &channel_plan,
      const double time,
      const CompiledMathWorkspace *workspace) {
    const auto source_id = channel_plan.request.source_id;
    const auto *kernel = source_kernel_for(channel_plan.source_kernel_slot);
    if (kernel == nullptr ||
        kernel->source_id != source_id ||
        kernel->kind != CompiledSourceChannelKernelKind::LeafAbsolute) {
      return impossible_channels();
    }
    base_current_cache_.set_context(time, 0);
    const auto pos = static_cast<std::size_t>(source_id);
    if (base_current_cache_.epoch[pos] == base_current_cache_.current_epoch) {
      return base_current_cache_.channels[pos];
    }
    const bool needs_frame = active_base_cache() != &base_current_cache_;
    if (needs_frame) {
      base_cache_stack_.push_back(&base_current_cache_);
    }
    base_current_cache_.channels[pos] =
        evaluate_leaf_absolute_kernel(*kernel, 0, workspace);
    base_current_cache_.epoch[pos] = base_current_cache_.current_epoch;
    if (needs_frame) {
      base_cache_stack_.pop_back();
    }
    return base_current_cache_.channels[pos];
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

  const CompiledMathCondition *condition_by_id(
      const semantic::Index condition_id) const {
    if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
      return nullptr;
    }
    const auto pos = static_cast<std::size_t>(condition_id - 1U);
    if (pos >= plan_.compiled_math.conditions.size()) {
      return nullptr;
    }
    return &plan_.compiled_math.conditions[pos];
  }

  static const CompiledMathPairFactEntry *condition_pair_entry_for(
      const CompiledMathPairFactLookup &lookup,
      const semantic::Index first,
      const semantic::Index second) {
    if (first == semantic::kInvalidIndex ||
        second == semantic::kInvalidIndex ||
        lookup.entries.empty()) {
      return nullptr;
    }
    const auto found = std::lower_bound(
        lookup.entries.begin(),
        lookup.entries.end(),
        std::pair<semantic::Index, semantic::Index>{first, second},
        [](const CompiledMathPairFactEntry &entry, const auto &target) {
          return entry.first < target.first ||
                 (entry.first == target.first &&
                  entry.second < target.second);
        });
    if (found == lookup.entries.end() ||
        found->first != first ||
        found->second != second) {
      return nullptr;
    }
    return &*found;
  }

  double condition_fact_time(
      const CompiledMathWorkspace *workspace,
      const CompiledMathCondition &condition,
      const std::size_t fact_index) const {
    const auto time_id =
        fact_index < condition.fact_time_ids.size()
            ? condition.fact_time_ids[fact_index]
            : static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
    if (workspace != nullptr && workspace->has_time(time_id)) {
      return workspace->time(time_id);
    }
    return conditional_time_;
  }

  bool condition_has_fact(
      const semantic::Index condition_id,
      const CompiledMathConditionFactKind target_kind) const {
    const auto *condition = condition_by_id(condition_id);
    if (condition == nullptr) {
      return false;
    }
    switch (target_kind) {
    case CompiledMathConditionFactKind::SourceExact:
      for (const auto &span : condition->source_exact_fact_spans) {
        if (!span.empty()) {
          return true;
        }
      }
      break;
    case CompiledMathConditionFactKind::SourceLowerBound:
      for (const auto &span : condition->source_lower_fact_spans) {
        if (!span.empty()) {
          return true;
        }
      }
      break;
    case CompiledMathConditionFactKind::SourceUpperBound:
      for (const auto &span : condition->source_upper_fact_spans) {
        if (!span.empty()) {
          return true;
        }
      }
      break;
    case CompiledMathConditionFactKind::ExprUpperBound:
      for (const auto &span : condition->expr_upper_fact_spans) {
        if (!span.empty()) {
          return true;
        }
      }
      break;
    case CompiledMathConditionFactKind::SourceOrder:
      return !condition->source_order_fact_lookup.entries.empty();
    case CompiledMathConditionFactKind::GuardUpperBound:
      return !condition->guard_upper_fact_lookup.entries.empty();
    }
    return false;
  }

  bool has_condition_exact_time(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    (void)workspace;
    if (source_id == semantic::kInvalidIndex) {
      return false;
    }
    return source_bound_plan(condition_id, source_id).has_condition_exact;
  }

  const double *exact_time_for(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    return resolve_source_bounds(
        source_id,
        condition_id,
        source_bound_plan(condition_id, source_id),
        workspace).exact_time;
  }

  ResolvedSourceBounds resolve_source_bounds(
      const semantic::Index source_id,
      const semantic::Index condition_id,
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

    const auto *condition = condition_by_id(condition_id);
    if (condition != nullptr && bounds.has_condition_lower) {
      for (semantic::Index i = 0; i < bounds.lower.size; ++i) {
        const auto fact_pos = static_cast<std::size_t>(
            condition->source_lower_fact_indices[
                static_cast<std::size_t>(bounds.lower.offset + i)]);
        resolved.lower = std::max(
            resolved.lower,
            condition_fact_time(workspace, *condition, fact_pos));
      }
      resolved.has_condition_lower = true;
    }
    if (condition != nullptr && bounds.has_condition_upper) {
      for (semantic::Index i = 0; i < bounds.upper.size; ++i) {
        const auto fact_pos = static_cast<std::size_t>(
            condition->source_upper_fact_indices[
                static_cast<std::size_t>(bounds.upper.offset + i)]);
        resolved.upper = std::min(
            resolved.upper,
            condition_fact_time(workspace, *condition, fact_pos));
      }
    }
    if (condition != nullptr && bounds.has_condition_exact) {
      const auto fact_pos = static_cast<std::size_t>(
          condition->source_exact_fact_indices[
              static_cast<std::size_t>(bounds.exact.offset)]);
      compiled_condition_time_ =
          condition_fact_time(workspace, *condition, fact_pos);
      resolved.exact_time = &compiled_condition_time_;
      return resolved;
    }
    if (sequence_exact) {
      resolved.exact_time = &sequence_state_->exact_times[source_pos];
    }
    return resolved;
  }

  ResolvedSourceBounds resolve_source_bounds(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    return resolve_source_bounds(
        source_id,
        condition_id,
        source_bound_plan(condition_id, source_id),
        workspace);
  }

  bool condition_source_order_known_before(
      const semantic::Index before_source_id,
      const semantic::Index after_source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    (void)workspace;
    const auto *condition = condition_by_id(condition_id);
    if (condition == nullptr) {
      return false;
    }
    if (condition_pair_entry_for(
            condition->source_order_fact_lookup,
            before_source_id,
            after_source_id) != nullptr) {
      return true;
    }
    bool have_before_relation = false;
    bool have_after_relation = false;
    ExactRelation before_relation{ExactRelation::Unknown};
    ExactRelation after_relation{ExactRelation::Unknown};
    for (std::size_t i = 0; i < condition->source_ids.size(); ++i) {
      if (condition->source_ids[i] == before_source_id) {
        before_relation = static_cast<ExactRelation>(condition->relations[i]);
        have_before_relation = true;
      } else if (condition->source_ids[i] == after_source_id) {
        after_relation = static_cast<ExactRelation>(condition->relations[i]);
        have_after_relation = true;
      }
    }
    if (have_before_relation && have_after_relation &&
        before_relation != ExactRelation::Unknown &&
        after_relation != ExactRelation::Unknown &&
        before_relation < after_relation) {
      return true;
    }
    return false;
  }

  semantic::Index acquire_recursive_cache() {
    if (recursive_cache_depth_ ==
        static_cast<semantic::Index>(recursive_base_caches_.size())) {
      recursive_base_caches_.emplace_back(
          static_cast<std::size_t>(plan_.source_count));
    }
    return recursive_cache_depth_++;
  }

  void push_base_cache(TimedCache *cache, const double t) {
    cache->set_time(t);
    base_cache_stack_.push_back(cache);
  }

  void pop_base_cache(const semantic::Index recursive_index) {
    base_cache_stack_.pop_back();
    if (recursive_index != -1) {
      --recursive_cache_depth_;
    }
  }

  void invalidate_all_timed_caches() {
    conditional_cache_.invalidate();
    base_current_cache_.invalidate();
    base_lower_cache_.invalidate();
    for (auto &cache : recursive_base_caches_) {
      cache.invalidate();
    }
  }

  TimedCache *active_base_cache() {
    if (base_cache_stack_.empty()) {
      return nullptr;
    }
    return base_cache_stack_.back();
  }

  double active_base_time() const {
    return base_cache_stack_.empty() ? conditional_time_
                                     : base_cache_stack_.back()->time;
  }

  leaf::EventChannels conditionalize(const leaf::EventChannels uncond,
                                     const leaf::EventChannels lower) const {
    const double surv_lb = lower.survival;
    if (!std::isfinite(surv_lb) || !(surv_lb > 0.0)) {
      return impossible_channels();
    }
    leaf::EventChannels out;
    out.pdf = safe_density(uncond.pdf / surv_lb);
    out.cdf = clamp_probability((uncond.cdf - lower.cdf) / surv_lb);
    out.survival = clamp_probability(uncond.survival / surv_lb);
    return out;
  }

  leaf::EventChannels conditionalize_between(
      const leaf::EventChannels uncond,
      const leaf::EventChannels lower,
      const leaf::EventChannels upper) const {
    const double mass = upper.cdf - lower.cdf;
    if (!std::isfinite(mass) || !(mass > 0.0)) {
      return impossible_channels();
    }
    leaf::EventChannels out;
    out.pdf = safe_density(uncond.pdf / mass);
    out.cdf = clamp_probability((uncond.cdf - lower.cdf) / mass);
    out.survival = clamp_probability((upper.cdf - uncond.cdf) / mass);
    return out;
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

  leaf::EventChannels evaluate_conditioned_source_kernel(
      const CompiledSourceChannelPlan &channel_plan,
      const double current_time,
      const CompiledMathWorkspace *workspace) {
    const auto &request = channel_plan.request;
    const auto source_id = request.source_id;
    const auto condition_id = request.condition_id;
    const auto *planned_kernel =
        source_kernel_for(channel_plan.source_kernel_slot);
    if (planned_kernel == nullptr ||
        planned_kernel->source_id != source_id ||
        planned_kernel->kind != channel_plan.kernel) {
      return impossible_channels();
    }
    const auto bounds = resolve_source_bounds(
        source_id, condition_id, channel_plan.bounds, workspace);
    const double lower_bound = bounds.lower;
    const double upper_bound = bounds.upper;
    if (std::isfinite(upper_bound) && !(upper_bound > lower_bound)) {
      return impossible_channels();
    }
    if (const double *exact_time = bounds.exact_time) {
      if (lower_bound > 0.0 && !(*exact_time > lower_bound)) {
        return impossible_channels();
      }
      if (std::isfinite(upper_bound) && !(*exact_time < upper_bound)) {
        return impossible_channels();
      }
      if (!(current_time >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }

    const auto uncond =
        load_base_source(
            &base_current_cache_, current_time, source_id, condition_id, workspace);
    if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
      return uncond;
    }

    leaf::EventChannels lower = leaf::EventChannels::impossible();
    if (!bounds.has_condition_lower &&
        std::fabs(lower_bound - sequence_state_->lower_bound) <= 1e-12) {
      lower = load_base_source(
          &base_lower_cache_, lower_bound, source_id, condition_id, workspace);
    } else if (lower_bound > 0.0) {
      const auto frame = base_time_guard(lower_bound);
      lower = base_source(source_id, condition_id, workspace);
    }
    if (!std::isfinite(upper_bound)) {
      return conditionalize(uncond, lower);
    }
    const auto frame = base_time_guard(upper_bound);
    const auto upper = base_source(source_id, condition_id, workspace);
    if (!std::isfinite(upper.cdf - lower.cdf) ||
        !(upper.cdf - lower.cdf > 0.0)) {
      return impossible_channels();
    }
    if (current_time >= upper_bound) {
      return leaf::EventChannels::certain();
    }
    if (current_time <= lower_bound) {
      return impossible_channels();
    }
    return conditionalize_between(uncond, lower, upper);
  }

  leaf::EventChannels base_source(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) {
    if (auto *cache = active_base_cache(); cache != nullptr) {
      return load_base_source(
          cache, cache->time, source_id, condition_id, workspace);
    }
    return load_base_source(
        &base_current_cache_, conditional_time_, source_id, condition_id, workspace);
  }

  leaf::EventChannels compute_base_source(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) {
    const auto bounds = resolve_source_bounds(source_id, condition_id, workspace);
    if (const double *exact_time = bounds.exact_time) {
      if (!(active_base_time() >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }
    const auto *kernel = source_kernel_for(source_id);
    if (kernel == nullptr) {
      return impossible_channels();
    }
    const auto source_pos = static_cast<std::size_t>(source_id);
    if (source_pos >= source_kernel_evaluators_.size() ||
        source_kernel_evaluators_[source_pos] == nullptr) {
      return impossible_channels();
    }
    return (this->*source_kernel_evaluators_[source_pos])(
        *kernel,
        condition_id,
        workspace);
  }

  const ExactSourceKernel *source_kernel_for(
      const semantic::Index source_id) const {
    if (source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(source_id) >= plan_.source_kernels.size()) {
      return nullptr;
    }
    return &plan_.source_kernels[static_cast<std::size_t>(source_id)];
  }

  leaf::EventChannels evaluate_leaf_absolute_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) {
    (void)condition_id;
    (void)workspace;
    const auto pos = static_cast<std::size_t>(kernel.leaf_index);
    if (kernel.leaf_index == semantic::kInvalidIndex ||
        pos >= program_.leaf_descriptors.size() ||
        pos >= leaf_inputs_.size()) {
      return impossible_channels();
    }
    const auto &desc = program_.leaf_descriptors[pos];
    const auto &loaded = leaf_inputs_[pos];
    const double t = active_base_time();
    return standard_leaf_channels(
        desc.dist_kind,
        loaded.params.data(),
        std::min(desc.param_count, 8),
        loaded.q,
        loaded.t0,
        t - desc.onset_abs_value);
  }

  leaf::EventChannels evaluate_leaf_onset_convolution_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) {
    const auto pos = static_cast<std::size_t>(kernel.leaf_index);
    if (kernel.leaf_index == semantic::kInvalidIndex ||
        pos >= program_.leaf_descriptors.size() ||
        pos >= leaf_inputs_.size()) {
      return impossible_channels();
    }
    const auto &desc = program_.leaf_descriptors[pos];
    const auto &loaded = leaf_inputs_[pos];
    const double t = active_base_time();

    const double lag = desc.onset_lag;
    const double upper = t - lag;
    if (!(upper > 0.0)) {
      return impossible_channels();
    }

    auto shifted_channels = [&](const double source_time) {
      return standard_leaf_channels(
          desc.dist_kind,
          loaded.params.data(),
          std::min(desc.param_count, 8),
          loaded.q,
          loaded.t0,
          t - source_time - lag);
    };

    const auto onset_source_id = kernel.onset_source_id;
    const auto onset_bounds =
        resolve_source_bounds(onset_source_id, condition_id, workspace);
    if (const double *exact_time = onset_bounds.exact_time) {
      return shifted_channels(*exact_time);
    }

    const auto batch = quadrature::build_finite_batch(0.0, upper);
    double pdf = 0.0;
    double cdf = 0.0;
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double u = batch.nodes.nodes[i];
      const auto frame = base_time_guard(u);
      const auto source = base_source(onset_source_id, condition_id, workspace);
      if (!(source.pdf > 0.0)) {
        continue;
      }
      const auto shifted = shifted_channels(u);
      const double weight = batch.nodes.weights[i] * source.pdf;
      pdf += weight * shifted.pdf;
      cdf += weight * shifted.cdf;
    }

    leaf::EventChannels out;
    out.pdf = safe_density(pdf);
    out.cdf = clamp_probability(cdf);
    out.survival = clamp_probability(1.0 - out.cdf);
    return out;
  }

  leaf::EventChannels evaluate_pool_kofn_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) {
    const auto begin = kernel.pool_member_offset;
    const auto end = begin + kernel.pool_member_count;
    const auto k = kernel.pool_k;
    if (end > static_cast<semantic::Index>(
                  program_.pool_member_source_ids.size())) {
      return impossible_channels();
    }
    const auto n_members = static_cast<int>(end - begin);

    std::vector<leaf::EventChannels> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(base_source(
          program_.pool_member_source_ids[static_cast<std::size_t>(i)],
          condition_id,
          workspace));
    }

    const auto width = static_cast<std::size_t>(n_members + 1);
    const auto table_size = width * width;
    std::vector<double> prefix(table_size, 0.0);
    std::vector<double> suffix(table_size, 0.0);
    const auto idx = [width](const int row, const int col) {
      return static_cast<std::size_t>(row) * width +
             static_cast<std::size_t>(col);
    };

    prefix[idx(0, 0)] = 1.0;
    for (int i = 0; i < n_members; ++i) {
      for (int m = 0; m <= i; ++m) {
        prefix[idx(i + 1, m)] +=
            prefix[idx(i, m)] *
            members[static_cast<std::size_t>(i)].survival;
        prefix[idx(i + 1, m + 1)] +=
            prefix[idx(i, m)] *
            members[static_cast<std::size_t>(i)].cdf;
      }
    }
    suffix[idx(n_members, 0)] = 1.0;
    for (int i = n_members - 1; i >= 0; --i) {
      const int count = n_members - i - 1;
      for (int m = 0; m <= count; ++m) {
        suffix[idx(i, m)] +=
            suffix[idx(i + 1, m)] *
            members[static_cast<std::size_t>(i)].survival;
        suffix[idx(i, m + 1)] +=
            suffix[idx(i + 1, m)] *
            members[static_cast<std::size_t>(i)].cdf;
      }
    }

    double surv = 0.0;
    for (int m = 0; m < k; ++m) {
      surv += prefix[idx(n_members, m)];
    }
    double pdf = 0.0;
    for (int i = 0; i < n_members; ++i) {
      double others_exact = 0.0;
      for (int left = 0; left < k; ++left) {
        const int right = (k - 1) - left;
        if (right < 0 || right > (n_members - i - 1)) {
          continue;
        }
        others_exact +=
            prefix[idx(i, left)] * suffix[idx(i + 1, right)];
      }
      pdf += members[static_cast<std::size_t>(i)].pdf * others_exact;
    }

    leaf::EventChannels out;
    out.pdf = safe_density(pdf);
    out.survival = clamp_probability(surv);
    out.cdf = clamp_probability(1.0 - out.survival);
    return out;
  }

  leaf::EventChannels load_conditional_source(
      TimedCache *cache,
      const double time,
      const CompiledSourceChannelPlan &channel_plan,
      const CompiledMathWorkspace *workspace) {
    const auto source_id = channel_plan.request.source_id;
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->epoch[pos] == cache->current_epoch) {
      return cache->channels[pos];
    }
    cache->channels[pos] =
        evaluate_conditioned_source_kernel(
            channel_plan, time, workspace);
    cache->epoch[pos] = cache->current_epoch;
    return cache->channels[pos];
  }

  leaf::EventChannels load_base_source(TimedCache *cache,
                                       const double time,
                                       const semantic::Index source_id,
                                       const semantic::Index condition_id,
                                       const CompiledMathWorkspace *workspace) {
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    cache->set_context(time, compiled_condition_cache_id(condition_id, workspace));
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->epoch[pos] == cache->current_epoch) {
      return cache->channels[pos];
    }
    const bool needs_frame = active_base_cache() != cache;
    if (needs_frame) {
      base_cache_stack_.push_back(cache);
    }
    cache->channels[pos] =
        compute_base_source(source_id, condition_id, workspace);
    cache->epoch[pos] = cache->current_epoch;
    if (needs_frame) {
      base_cache_stack_.pop_back();
    }
    return cache->channels[pos];
  }
};

using ExactSourceChannels = CompiledSourceChannels;

inline void enumerate_trigger_states_into(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    std::vector<ExactTriggerState> *states,
    std::vector<ExactTriggerState> *buffer = nullptr) {
  if (states == nullptr) {
    return;
  }
  states->clear();
  states->resize(1);
  states->front().weight = 1.0;
  states->front().shared_started.assign(
      static_cast<std::size_t>(plan.lowered.program.layout.n_triggers), 2U);

  for (const auto trigger_index : plan.shared_trigger_indices) {
    const auto pos = static_cast<std::size_t>(trigger_index);
    double q = 0.0;
    if (plan.lowered.program.trigger_has_fixed_q[pos] != 0U) {
      q = plan.lowered.program.trigger_fixed_q[pos];
    } else {
      const auto member_begin = plan.lowered.program.trigger_member_offsets[pos];
      if (member_begin != plan.lowered.program.trigger_member_offsets[pos + 1U]) {
        q = params.q(first_param_row +
                     plan.lowered.program.trigger_member_indices[static_cast<std::size_t>(
                         member_begin)]);
      }
    }
    q = clamp_probability(q);

    std::vector<ExactTriggerState> local_buffer;
    auto *next = buffer == nullptr ? &local_buffer : buffer;
    next->clear();
    next->reserve(states->size() * 2U);
    for (const auto &state : *states) {
      if (q > 0.0) {
        auto fail = state;
        fail.weight *= q;
        fail.shared_started[pos] = 0U;
        next->push_back(std::move(fail));
      }
      if (q < 1.0) {
        auto start = state;
        start.weight *= (1.0 - q);
        start.shared_started[pos] = 1U;
        next->push_back(std::move(start));
      }
    }
    states->swap(*next);
  }
}

inline std::vector<ExactTriggerState> enumerate_trigger_states(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row) {
  std::vector<ExactTriggerState> states;
  enumerate_trigger_states_into(plan, params, first_param_row, &states);
  return states;
}

} // namespace detail
} // namespace accumulatr::eval
