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
  struct SourceProductScalarFill {
    std::uint8_t mask{0U};
    double pdf{0.0};
    double cdf{0.0};
    double survival{1.0};
  };

  struct TimedCache {
    explicit TimedCache(const std::size_t size = 0U)
        : channels(size), epoch(size, 0U), valid_mask(size, 0U) {}

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

    bool has(const std::size_t pos, const std::uint8_t mask) const {
      return epoch[pos] == current_epoch && (valid_mask[pos] & mask) == mask;
    }

    leaf::EventChannels get(const std::size_t pos) const {
      return channels[pos];
    }

    void store(const std::size_t pos,
               const leaf::EventChannels value,
               const std::uint8_t mask) {
      if (mask == kLeafChannelAll) {
        channels[pos] = value;
        valid_mask[pos] = kLeafChannelAll;
        epoch[pos] = current_epoch;
        return;
      }
      if (epoch[pos] != current_epoch) {
        channels[pos] = leaf::EventChannels::impossible();
        valid_mask[pos] = 0U;
        epoch[pos] = current_epoch;
      }
      if ((mask & kLeafChannelPdf) != 0U) {
        channels[pos].pdf = value.pdf;
      }
      if ((mask & kLeafChannelCdf) != 0U) {
        channels[pos].cdf = value.cdf;
      }
      if ((mask & kLeafChannelSurvival) != 0U) {
        channels[pos].survival = value.survival;
      }
      valid_mask[pos] |= mask;
    }

    double time{0.0};
    semantic::Index condition_id{0};
    std::vector<leaf::EventChannels> channels;
    std::vector<std::uint32_t> epoch;
    std::vector<std::uint8_t> valid_mask;
    std::uint32_t current_epoch{1U};
  };

  struct SourceProductTimedCache {
    explicit SourceProductTimedCache(const std::size_t size = 0U)
        : fills(size), epoch(size, 0U), valid_mask(size, 0U) {}

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

    bool has(const std::size_t pos, const std::uint8_t mask) const {
      return epoch[pos] == current_epoch && (valid_mask[pos] & mask) == mask;
    }

    SourceProductScalarFill get(const std::size_t pos) const {
      return fills[pos];
    }

    void store(const std::size_t pos,
               const SourceProductScalarFill value,
               const std::uint8_t mask) {
      if (mask == kLeafChannelAll) {
        fills[pos] = value;
        valid_mask[pos] = kLeafChannelAll;
        epoch[pos] = current_epoch;
        return;
      }
      if (epoch[pos] != current_epoch) {
        fills[pos] = SourceProductScalarFill{};
        valid_mask[pos] = 0U;
        epoch[pos] = current_epoch;
      }
      fills[pos].mask |= mask;
      if ((mask & kLeafChannelPdf) != 0U) {
        fills[pos].pdf = value.pdf;
      }
      if ((mask & kLeafChannelCdf) != 0U) {
        fills[pos].cdf = value.cdf;
      }
      if ((mask & kLeafChannelSurvival) != 0U) {
        fills[pos].survival = value.survival;
      }
      valid_mask[pos] |= mask;
    }

    double time{0.0};
    semantic::Index condition_id{0};
    std::vector<SourceProductScalarFill> fills;
    std::vector<std::uint32_t> epoch;
    std::vector<std::uint8_t> valid_mask;
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

  class SourceProductBaseTimeGuard {
  public:
    SourceProductBaseTimeGuard(CompiledSourceChannels *channels,
                               semantic::Index recursive_index,
                               SourceProductTimedCache *cache,
                               const double time)
        : channels_(channels), recursive_index_(recursive_index) {
      if (channels_ == nullptr || cache == nullptr) {
        return;
      }
      channels_->push_source_product_base_cache(cache, time);
      active_ = true;
    }

    SourceProductBaseTimeGuard(const SourceProductBaseTimeGuard &) = delete;
    SourceProductBaseTimeGuard &operator=(
        const SourceProductBaseTimeGuard &) = delete;

    ~SourceProductBaseTimeGuard() {
      if (channels_ != nullptr && active_) {
        channels_->pop_source_product_base_cache(recursive_index_);
      }
    }

  private:
    CompiledSourceChannels *channels_{nullptr};
    semantic::Index recursive_index_{-1};
    bool active_{false};
  };


  explicit CompiledSourceChannels(const ExactVariantPlan &plan)
      : plan_(plan),
        program_(plan.lowered.program),
        leaf_inputs_(static_cast<std::size_t>(program_.layout.n_leaves)),
        source_product_direct_available_(
            plan.compiled_math.integral_kernel_source_product_channels.size(),
            0U),
        conditional_cache_(static_cast<std::size_t>(plan.source_count)),
        base_current_cache_(static_cast<std::size_t>(plan.source_count)),
        base_lower_cache_(static_cast<std::size_t>(plan.source_count)),
        source_product_conditional_cache_(
            static_cast<std::size_t>(plan.source_count)),
        source_product_base_current_cache_(
            static_cast<std::size_t>(plan.source_count)),
        source_product_base_lower_cache_(
            static_cast<std::size_t>(plan.source_count)) {
    recursive_base_caches_.reserve(
        static_cast<std::size_t>(std::max(program_.layout.n_leaves, 1)));
    base_cache_stack_.reserve(
        static_cast<std::size_t>(std::max(program_.layout.n_leaves, 1)));
    source_product_recursive_base_caches_.reserve(
        static_cast<std::size_t>(std::max(program_.layout.n_leaves, 1)));
    source_product_base_cache_stack_.reserve(
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
    source_product_recursive_cache_depth_ = 0;
    base_cache_stack_.clear();
    source_product_base_cache_stack_.clear();
    conditional_cache_.reset(top_time);
    base_current_cache_.reset(top_time);
    base_lower_cache_.reset(sequence_state.lower_bound);
    source_product_conditional_cache_.reset(top_time);
    source_product_base_current_cache_.reset(top_time);
    source_product_base_lower_cache_.reset(sequence_state.lower_bound);
    for (auto &cache : recursive_base_caches_) {
      cache.invalidate();
    }
    for (auto &cache : source_product_recursive_base_caches_) {
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

  [[nodiscard]] BaseTimeGuard base_time_guard(const double time) {
    const auto recursive_index = acquire_recursive_cache();
    return BaseTimeGuard(
        this,
        recursive_index,
        &recursive_base_caches_[static_cast<std::size_t>(recursive_index)],
        time);
  }

  [[nodiscard]] SourceProductBaseTimeGuard source_product_base_time_guard(
      const double time) {
    const auto recursive_index = acquire_source_product_recursive_cache();
    return SourceProductBaseTimeGuard(
        this,
        recursive_index,
        &source_product_recursive_base_caches_[
            static_cast<std::size_t>(recursive_index)],
        time);
  }

  leaf::EventChannels evaluate_source_channel_plan_at(
      const CompiledSourceChannelPlan &channel_plan,
      const double time,
      const CompiledMathWorkspace *workspace = nullptr,
      const std::uint8_t channel_mask = kLeafChannelAll) {
    if (channel_plan.direct_leaf_absolute_candidate &&
        sequence_bounds_inactive_for(
            channel_plan.request.source_id, channel_plan.bounds)) {
      return evaluate_direct_leaf_absolute_at(
          channel_plan,
          time,
          workspace,
          channel_mask);
    }
    conditional_cache_.set_context(
        time, condition_cache_id(channel_plan, workspace));
    return load_conditional_source(
        &conditional_cache_, time, channel_plan, workspace, channel_mask);
  }

  SourceProductScalarFill evaluate_source_product_channel_fill(
      const CompiledMathSourceProductChannel &channel,
      const std::uint8_t requested_mask,
      const double time,
      const semantic::Index parent_source_view_id,
      const CompiledMathWorkspace *workspace) {
    const auto effective_source_view_id =
        channel.source_view_id == 0 ||
                channel.source_view_id == semantic::kInvalidIndex
            ? parent_source_view_id
            : channel.source_view_id;
    const auto relation =
        channel.has_static_source_view_relation
            ? static_cast<ExactRelation>(channel.static_source_view_relation)
            : exact_compiled_source_view_relation(
                  plan_,
                  effective_source_view_id,
                  channel.source_id);
    if (relation != ExactRelation::Unknown &&
        !channel.has_source_condition_overlay) {
      return source_product_forced_fill(relation, requested_mask);
    }
    if (channel.direct_leaf_absolute_candidate &&
        sequence_bounds_inactive_for(channel.source_id, channel.bounds)) {
      return evaluate_source_product_leaf_absolute_fill(
          channel,
          time,
          requested_mask);
    }
    source_product_conditional_cache_.set_context(
        time, source_product_condition_cache_id(channel, workspace));
    return load_source_product_conditioned_fill(
        &source_product_conditional_cache_,
        time,
        channel,
        workspace,
        requested_mask);
  }

  bool source_product_direct_leaf_available(
      const CompiledMathSourceProductChannel &channel) const {
    return channel.direct_leaf_absolute_candidate &&
           sequence_bounds_inactive_for(channel.source_id, channel.bounds);
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
  TimedCache conditional_cache_;
  TimedCache base_current_cache_;
  TimedCache base_lower_cache_;
  std::vector<TimedCache> recursive_base_caches_;
  std::vector<TimedCache *> base_cache_stack_;
  SourceProductTimedCache source_product_conditional_cache_;
  SourceProductTimedCache source_product_base_current_cache_;
  SourceProductTimedCache source_product_base_lower_cache_;
  std::vector<SourceProductTimedCache> source_product_recursive_base_caches_;
  std::vector<SourceProductTimedCache *> source_product_base_cache_stack_;
  using SourceKernelEvaluator = leaf::EventChannels (CompiledSourceChannels::*)(
      const ExactSourceKernel &,
      semantic::Index,
      const CompiledMathWorkspace *,
      std::uint8_t);
  std::vector<SourceKernelEvaluator> source_kernel_evaluators_;
  mutable double compiled_condition_time_{
      std::numeric_limits<double>::quiet_NaN()};
  semantic::Index recursive_cache_depth_{0};
  semantic::Index source_product_recursive_cache_depth_{0};

  CompiledSourceBoundPlan empty_source_bound_plan_{};

  struct ResolvedSourceBounds {
    const double *exact_time{nullptr};
    double lower{0.0};
    double upper{std::numeric_limits<double>::infinity()};
    bool has_condition_lower{false};
  };

  semantic::Index dynamic_condition_cache_id(
      const semantic::Index condition_id,
      const semantic::Index static_cache_id,
      const CompiledMathIndexSpan time_dependencies,
      const CompiledMathWorkspace *workspace) const {
    if (time_dependencies.size == 0 || workspace == nullptr) {
      return static_cache_id;
    }
    std::size_t seed = 0x517cc1b727220a95ULL;
    compiled_math_condition_cache_hash_part(
        &seed, static_cast<std::size_t>(condition_id));
    for (semantic::Index i = 0; i < time_dependencies.size; ++i) {
      const auto time_id =
          plan_.compiled_math.condition_cache_time_dependencies[
              static_cast<std::size_t>(time_dependencies.offset + i)];
      const auto time_pos = static_cast<std::size_t>(time_id);
      compiled_math_condition_cache_hash_part(
          &seed, static_cast<std::size_t>(time_id) + 1U);
      compiled_math_condition_cache_hash_part(
          &seed,
          std::hash<double>{}(workspace->time_values[time_pos]));
    }
    return compiled_math_condition_cache_index(seed);
  }

  semantic::Index condition_cache_id(
      const CompiledSourceChannelPlan &channel_plan,
      const CompiledMathWorkspace *workspace) const {
    if (!channel_plan.condition_cache_dynamic) {
      return channel_plan.condition_static_cache_id;
    }
    return dynamic_condition_cache_id(
        channel_plan.request.condition_id,
        channel_plan.condition_static_cache_id,
        channel_plan.condition_time_dependencies,
        workspace);
  }

  semantic::Index source_product_condition_cache_id(
      const CompiledMathSourceProductChannel &channel,
      const CompiledMathWorkspace *workspace) const {
    if (!channel.condition_cache_dynamic) {
      return channel.condition_static_cache_id;
    }
    return dynamic_condition_cache_id(
        channel.condition_id,
        channel.condition_static_cache_id,
        channel.condition_time_dependencies,
        workspace);
  }

  semantic::Index condition_cache_id(
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace) const {
    if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
      return 0;
    }
    const auto condition_pos = static_cast<std::size_t>(condition_id);
    if (condition_pos >= plan_.compiled_math.condition_cache_plans.size()) {
      return 0;
    }
    const auto &cache_plan =
        plan_.compiled_math.condition_cache_plans[condition_pos];
    if (!cache_plan.dynamic) {
      return cache_plan.static_cache_id;
    }
    return dynamic_condition_cache_id(
        condition_id,
        cache_plan.static_cache_id,
        cache_plan.time_dependencies,
        workspace);
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

  leaf::EventChannels evaluate_direct_leaf_absolute_at(
      const CompiledSourceChannelPlan &channel_plan,
      const double time,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask) {
    const auto source_id = channel_plan.request.source_id;
    const auto &kernel =
        plan_.source_kernels[
            static_cast<std::size_t>(channel_plan.source_kernel_slot)];
    if (channel_mask != kLeafChannelAll) {
      return evaluate_leaf_absolute_kernel_at(kernel, time, channel_mask);
    }
    base_current_cache_.set_context(time, 0);
    const auto pos = static_cast<std::size_t>(source_id);
    if (base_current_cache_.has(pos, kLeafChannelAll)) {
      return base_current_cache_.channels[pos];
    }
    const bool needs_frame = active_base_cache() != &base_current_cache_;
    if (needs_frame) {
      base_cache_stack_.push_back(&base_current_cache_);
    }
    base_current_cache_.channels[pos] =
        evaluate_leaf_absolute_kernel(kernel, 0, workspace, kLeafChannelAll);
    base_current_cache_.epoch[pos] = base_current_cache_.current_epoch;
    base_current_cache_.valid_mask[pos] = kLeafChannelAll;
    if (needs_frame) {
      base_cache_stack_.pop_back();
    }
    return base_current_cache_.channels[pos];
  }

  static SourceProductScalarFill source_product_impossible_fill(
      const std::uint8_t mask) {
    SourceProductScalarFill fill;
    fill.mask = mask;
    return fill;
  }

  static SourceProductScalarFill source_product_certain_fill(
      const std::uint8_t mask) {
    SourceProductScalarFill fill;
    fill.mask = mask;
    if ((mask & kLeafChannelCdf) != 0U) {
      fill.cdf = 1.0;
    }
    if ((mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 0.0;
    }
    return fill;
  }

  static SourceProductScalarFill source_product_forced_fill(
      const ExactRelation relation,
      const std::uint8_t mask) {
    SourceProductScalarFill fill;
    fill.mask = mask;
    if (relation == ExactRelation::Before || relation == ExactRelation::At) {
      if ((mask & kLeafChannelCdf) != 0U) {
        fill.cdf = 1.0;
      }
      if ((mask & kLeafChannelSurvival) != 0U) {
        fill.survival = 0.0;
      }
      return fill;
    }
    if (relation == ExactRelation::After) {
      if ((mask & kLeafChannelSurvival) != 0U) {
        fill.survival = 1.0;
      }
      return fill;
    }
    return fill;
  }

  SourceProductScalarFill evaluate_source_product_leaf_fill_at(
      const semantic::Index leaf_index,
      const double leaf_time,
      const std::uint8_t fill_mask) const {
    const auto leaf_pos = static_cast<std::size_t>(leaf_index);
    if (leaf_index == semantic::kInvalidIndex ||
        leaf_pos >= leaf_inputs_.size() ||
        leaf_pos >= program_.leaf_descriptors.size()) {
      return source_product_impossible_fill(fill_mask);
    }
    const auto &desc = program_.leaf_descriptors[leaf_pos];
    const auto &loaded = leaf_inputs_[leaf_pos];
    const double x = leaf_time - loaded.t0;
    if (!(x > 0.0)) {
      return source_product_impossible_fill(fill_mask);
    }
    switch (static_cast<leaf::DistKind>(desc.dist_kind)) {
    case leaf::DistKind::Lognormal:
      return source_product_lognormal_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::Gamma:
      return source_product_gamma_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::Exgauss:
      return source_product_exgauss_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::LBA:
      return source_product_lba_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::RDM:
      return source_product_rdm_leaf_fill(loaded, x, fill_mask);
    }
    return source_product_impossible_fill(fill_mask);
  }

  SourceProductScalarFill evaluate_source_product_leaf_absolute_fill(
      const CompiledMathSourceProductChannel &channel,
      const double time,
      const std::uint8_t fill_mask) const {
    const auto leaf_pos = static_cast<std::size_t>(channel.leaf_index);
    if (channel.leaf_index == semantic::kInvalidIndex ||
        leaf_pos >= leaf_inputs_.size()) {
      SourceProductScalarFill fill;
      fill.mask = fill_mask;
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        fill.survival = 1.0;
      }
      return fill;
    }
    const auto &loaded = leaf_inputs_[leaf_pos];
    const double x = time - channel.leaf_onset_abs_value - loaded.t0;
    if (!(x > 0.0)) {
      SourceProductScalarFill fill;
      fill.mask = fill_mask;
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        fill.survival = 1.0;
      }
      return fill;
    }
    switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
    case leaf::DistKind::Lognormal:
      return source_product_lognormal_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::Gamma:
      return source_product_gamma_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::Exgauss:
      return source_product_exgauss_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::LBA:
      return source_product_lba_leaf_fill(loaded, x, fill_mask);
    case leaf::DistKind::RDM:
      return source_product_rdm_leaf_fill(loaded, x, fill_mask);
    }
    return SourceProductScalarFill{};
  }

  static SourceProductScalarFill source_product_finish_base_fill(
      const double base_pdf,
      const double base_cdf,
      const double q,
      const std::uint8_t fill_mask) {
    const double start_prob = 1.0 - q;
    SourceProductScalarFill fill;
    fill.mask = fill_mask;
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      fill.pdf = start_prob * safe_density(base_pdf);
    }
    if ((fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
      const double cdf =
          clamp_probability(start_prob * clamp_probability(base_cdf));
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        fill.cdf = cdf;
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        fill.survival = clamp_probability(1.0 - cdf);
      }
    }
    return fill;
  }

  static SourceProductScalarFill source_product_lognormal_leaf_fill(
      const ExactLoadedLeafInput &loaded,
      const double x,
      const std::uint8_t fill_mask) {
    const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    const double m = loaded.params[0];
    const double s = loaded.params[1];
    return source_product_finish_base_fill(
        need_pdf ? R::dlnorm(x, m, s, 0) : 0.0,
        need_cdf ? R::plnorm(x, m, s, 1, 0) : 0.0,
        loaded.q,
        fill_mask);
  }

  static SourceProductScalarFill source_product_gamma_leaf_fill(
      const ExactLoadedLeafInput &loaded,
      const double x,
      const std::uint8_t fill_mask) {
    const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    const double shape = loaded.params[0];
    const double scale = 1.0 / loaded.params[1];
    return source_product_finish_base_fill(
        need_pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
        need_cdf ? R::pgamma(x, shape, scale, 1, 0) : 0.0,
        loaded.q,
        fill_mask);
  }

  static SourceProductScalarFill source_product_exgauss_leaf_fill(
      const ExactLoadedLeafInput &loaded,
      const double x,
      const std::uint8_t fill_mask) {
    const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    const double mu = loaded.params[0];
    const double sigma = loaded.params[1];
    const double tau = loaded.params[2];
    const double lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
    const double lower_survival = 1.0 - lower_cdf;
    if (!(lower_survival > 0.0)) {
      SourceProductScalarFill fill;
      fill.mask = fill_mask;
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        fill.survival = 1.0;
      }
      return fill;
    }
    return source_product_finish_base_fill(
        need_pdf ? exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival : 0.0,
        need_cdf
            ? (exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) /
                  lower_survival
            : 0.0,
        loaded.q,
        fill_mask);
  }

  static SourceProductScalarFill source_product_lba_leaf_fill(
      const ExactLoadedLeafInput &loaded,
      const double x,
      const std::uint8_t fill_mask) {
    const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    return source_product_finish_base_fill(
        need_pdf
            ? lba_pdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3])
            : 0.0,
        need_cdf
            ? lba_cdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3])
            : 0.0,
        loaded.q,
        fill_mask);
  }

  static SourceProductScalarFill source_product_rdm_leaf_fill(
      const ExactLoadedLeafInput &loaded,
      const double x,
      const std::uint8_t fill_mask) {
    const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    return source_product_finish_base_fill(
        need_pdf
            ? rdm_pdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3])
            : 0.0,
        need_cdf
            ? rdm_cdf_fast(x, loaded.params[0], loaded.params[1],
                           loaded.params[2], loaded.params[3])
            : 0.0,
        loaded.q,
        fill_mask);
  }

  leaf::EventChannels evaluate_leaf_absolute_kernel_at(
      const ExactSourceKernel &kernel,
      const double time,
      const std::uint8_t channel_mask) const {
    const auto pos = static_cast<std::size_t>(kernel.leaf_index);
    if (kernel.leaf_index == semantic::kInvalidIndex ||
        pos >= program_.leaf_descriptors.size() ||
        pos >= leaf_inputs_.size()) {
      return impossible_channels();
    }
    const auto &desc = program_.leaf_descriptors[pos];
    const auto &loaded = leaf_inputs_[pos];
    return standard_leaf_channels_mask(
        desc.dist_kind,
        loaded.params.data(),
        std::min(desc.param_count, 8),
        loaded.q,
        loaded.t0,
        time - desc.onset_abs_value,
        channel_mask);
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

  semantic::Index acquire_recursive_cache() {
    if (recursive_cache_depth_ ==
        static_cast<semantic::Index>(recursive_base_caches_.size())) {
      recursive_base_caches_.emplace_back(
          static_cast<std::size_t>(plan_.source_count));
    }
    return recursive_cache_depth_++;
  }

  semantic::Index acquire_source_product_recursive_cache() {
    if (source_product_recursive_cache_depth_ ==
        static_cast<semantic::Index>(
            source_product_recursive_base_caches_.size())) {
      source_product_recursive_base_caches_.emplace_back(
          static_cast<std::size_t>(plan_.source_count));
    }
    return source_product_recursive_cache_depth_++;
  }

  void push_base_cache(TimedCache *cache, const double t) {
    cache->set_time(t);
    base_cache_stack_.push_back(cache);
  }

  void push_source_product_base_cache(SourceProductTimedCache *cache,
                                      const double t) {
    cache->set_time(t);
    source_product_base_cache_stack_.push_back(cache);
  }

  void pop_base_cache(const semantic::Index recursive_index) {
    base_cache_stack_.pop_back();
    if (recursive_index != -1) {
      --recursive_cache_depth_;
    }
  }

  void pop_source_product_base_cache(const semantic::Index recursive_index) {
    source_product_base_cache_stack_.pop_back();
    if (recursive_index != -1) {
      --source_product_recursive_cache_depth_;
    }
  }

  void invalidate_all_timed_caches() {
    conditional_cache_.invalidate();
    base_current_cache_.invalidate();
    base_lower_cache_.invalidate();
    for (auto &cache : recursive_base_caches_) {
      cache.invalidate();
    }
    source_product_conditional_cache_.invalidate();
    source_product_base_current_cache_.invalidate();
    source_product_base_lower_cache_.invalidate();
    for (auto &cache : source_product_recursive_base_caches_) {
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

  SourceProductTimedCache *active_source_product_base_cache() {
    if (source_product_base_cache_stack_.empty()) {
      return nullptr;
    }
    return source_product_base_cache_stack_.back();
  }

  double active_source_product_base_time() const {
    return source_product_base_cache_stack_.empty()
               ? conditional_time_
               : source_product_base_cache_stack_.back()->time;
  }

  leaf::EventChannels conditionalize(const leaf::EventChannels uncond,
                                     const leaf::EventChannels lower,
                                     const std::uint8_t channel_mask) const {
    const double surv_lb = lower.survival;
    if (!std::isfinite(surv_lb) || !(surv_lb > 0.0)) {
      return impossible_channels();
    }
    leaf::EventChannels out;
    if ((channel_mask & kLeafChannelPdf) != 0U) {
      out.pdf = safe_density(uncond.pdf / surv_lb);
    }
    if ((channel_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamp_probability(
          (uncond.cdf + surv_lb - 1.0) / surv_lb);
    }
    if ((channel_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability(uncond.survival / surv_lb);
    }
    return out;
  }

  leaf::EventChannels conditionalize_between(
      const leaf::EventChannels uncond,
      const leaf::EventChannels lower,
      const leaf::EventChannels upper,
      const std::uint8_t channel_mask) const {
    const double mass = upper.cdf - lower.cdf;
    if (!std::isfinite(mass) || !(mass > 0.0)) {
      return impossible_channels();
    }
    leaf::EventChannels out;
    if ((channel_mask & kLeafChannelPdf) != 0U) {
      out.pdf = safe_density(uncond.pdf / mass);
    }
    if ((channel_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamp_probability((uncond.cdf - lower.cdf) / mass);
    }
    if ((channel_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability((upper.cdf - uncond.cdf) / mass);
    }
    return out;
  }

  SourceProductScalarFill source_product_conditionalize(
      const SourceProductScalarFill uncond,
      const SourceProductScalarFill lower,
      const std::uint8_t fill_mask) const {
    const double surv_lb = lower.survival;
    if (!std::isfinite(surv_lb) || !(surv_lb > 0.0)) {
      return source_product_impossible_fill(fill_mask);
    }
    SourceProductScalarFill out;
    out.mask = fill_mask;
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      out.pdf = safe_density(uncond.pdf / surv_lb);
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamp_probability(
          (uncond.cdf + surv_lb - 1.0) / surv_lb);
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability(uncond.survival / surv_lb);
    }
    return out;
  }

  SourceProductScalarFill source_product_conditionalize_between(
      const SourceProductScalarFill uncond,
      const SourceProductScalarFill lower,
      const SourceProductScalarFill upper,
      const std::uint8_t fill_mask) const {
    const double mass = upper.cdf - lower.cdf;
    if (!std::isfinite(mass) || !(mass > 0.0)) {
      return source_product_impossible_fill(fill_mask);
    }
    SourceProductScalarFill out;
    out.mask = fill_mask;
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      out.pdf = safe_density(uncond.pdf / mass);
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamp_probability((uncond.cdf - lower.cdf) / mass);
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability((upper.cdf - uncond.cdf) / mass);
    }
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
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask) {
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
    const auto bounds =
        compiled_source_bounds(source_id, channel_plan.bounds, workspace);
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
            &base_current_cache_,
            current_time,
            source_id,
            condition_id,
            workspace,
            channel_mask);
    if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
      return uncond;
    }

    leaf::EventChannels lower = leaf::EventChannels::impossible();
    std::uint8_t lower_mask = kLeafChannelSurvival;
    if (std::isfinite(upper_bound)) {
      lower_mask = kLeafChannelCdf;
    }
    if (!bounds.has_condition_lower &&
        std::fabs(lower_bound - sequence_state_->lower_bound) <= 1e-12) {
      lower = load_base_source(
          &base_lower_cache_,
          lower_bound,
          source_id,
          condition_id,
          workspace,
          lower_mask);
    } else if (lower_bound > 0.0) {
      const auto frame = base_time_guard(lower_bound);
      lower = base_source(source_id, condition_id, workspace, lower_mask);
    }
    if (!std::isfinite(upper_bound)) {
      return conditionalize(uncond, lower, channel_mask);
    }
    const auto frame = base_time_guard(upper_bound);
    const auto upper =
        base_source(source_id, condition_id, workspace, kLeafChannelCdf);
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
    return conditionalize_between(uncond, lower, upper, channel_mask);
  }

  SourceProductScalarFill evaluate_source_product_conditioned_fill(
      const CompiledMathSourceProductChannel &channel,
      const double current_time,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask) {
    const auto source_id = channel.source_id;
    const auto condition_id = channel.condition_id;
    const auto bounds =
        compiled_source_bounds(source_id, channel.bounds, workspace);
    const double lower_bound = bounds.lower;
    const double upper_bound = bounds.upper;
    if (std::isfinite(upper_bound) && !(upper_bound > lower_bound)) {
      return source_product_impossible_fill(fill_mask);
    }
    if (const double *exact_time = bounds.exact_time) {
      if (lower_bound > 0.0 && !(*exact_time > lower_bound)) {
        return source_product_impossible_fill(fill_mask);
      }
      if (std::isfinite(upper_bound) && !(*exact_time < upper_bound)) {
        return source_product_impossible_fill(fill_mask);
      }
      if (!(current_time >= *exact_time)) {
        return source_product_impossible_fill(fill_mask);
      }
      return source_product_certain_fill(fill_mask);
    }

    std::uint8_t uncond_mask = fill_mask;
    if (std::isfinite(upper_bound) &&
        (fill_mask & kLeafChannelSurvival) != 0U) {
      uncond_mask |= kLeafChannelCdf;
    }
    const auto uncond =
        load_source_product_base_fill(
            &source_product_base_current_cache_,
            current_time,
            source_id,
            condition_id,
            workspace,
            uncond_mask);
    if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
      return uncond;
    }

    SourceProductScalarFill lower = source_product_impossible_fill(
        kLeafChannelSurvival);
    std::uint8_t lower_mask = kLeafChannelSurvival;
    if (std::isfinite(upper_bound)) {
      lower_mask = kLeafChannelCdf;
    }
    if (!bounds.has_condition_lower &&
        sequence_state_ != nullptr &&
        std::fabs(lower_bound - sequence_state_->lower_bound) <= 1e-12) {
      lower = load_source_product_base_fill(
          &source_product_base_lower_cache_,
          lower_bound,
          source_id,
          condition_id,
          workspace,
          lower_mask);
    } else if (lower_bound > 0.0) {
      const auto frame = source_product_base_time_guard(lower_bound);
      lower = source_product_base_fill(
          source_id, condition_id, workspace, lower_mask);
    }
    if (!std::isfinite(upper_bound)) {
      return source_product_conditionalize(uncond, lower, fill_mask);
    }
    const auto frame = source_product_base_time_guard(upper_bound);
    const auto upper =
        source_product_base_fill(
            source_id, condition_id, workspace, kLeafChannelCdf);
    if (!std::isfinite(upper.cdf - lower.cdf) ||
        !(upper.cdf - lower.cdf > 0.0)) {
      return source_product_impossible_fill(fill_mask);
    }
    if (current_time >= upper_bound) {
      return source_product_certain_fill(fill_mask);
    }
    if (current_time <= lower_bound) {
      return source_product_impossible_fill(fill_mask);
    }
    return source_product_conditionalize_between(
        uncond, lower, upper, fill_mask);
  }

  SourceProductScalarFill source_product_base_fill(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask = kLeafChannelAll) {
    if (auto *cache = active_source_product_base_cache(); cache != nullptr) {
      return load_source_product_base_fill(
          cache, cache->time, source_id, condition_id, workspace, fill_mask);
    }
    return load_source_product_base_fill(
        &source_product_base_current_cache_,
        conditional_time_,
        source_id,
        condition_id,
        workspace,
        fill_mask);
  }

  SourceProductScalarFill compute_source_product_base_fill(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask) {
    const auto bounds =
        compiled_source_bounds(
            source_id, source_bound_plan(condition_id, source_id), workspace);
    if (const double *exact_time = bounds.exact_time) {
      if (!(active_source_product_base_time() >= *exact_time)) {
        return source_product_impossible_fill(fill_mask);
      }
      return source_product_certain_fill(fill_mask);
    }
    const auto *kernel = source_kernel_for(source_id);
    if (kernel == nullptr) {
      return source_product_impossible_fill(fill_mask);
    }
    switch (kernel->kind) {
    case CompiledSourceChannelKernelKind::LeafAbsolute:
      return evaluate_source_product_leaf_absolute_kernel(
          *kernel, condition_id, workspace, fill_mask);
    case CompiledSourceChannelKernelKind::LeafOnsetConvolution:
      return evaluate_source_product_leaf_onset_convolution_kernel(
          *kernel, condition_id, workspace, fill_mask);
    case CompiledSourceChannelKernelKind::PoolKOfN:
      return evaluate_source_product_pool_kofn_kernel(
          *kernel, condition_id, workspace, fill_mask);
    case CompiledSourceChannelKernelKind::Invalid:
      break;
    }
    return source_product_impossible_fill(fill_mask);
  }

  leaf::EventChannels base_source(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask = kLeafChannelAll) {
    if (auto *cache = active_base_cache(); cache != nullptr) {
      return load_base_source(
          cache, cache->time, source_id, condition_id, workspace, channel_mask);
    }
    return load_base_source(
        &base_current_cache_,
        conditional_time_,
        source_id,
        condition_id,
        workspace,
        channel_mask);
  }

  leaf::EventChannels compute_base_source(
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask) {
    const auto bounds =
        compiled_source_bounds(
            source_id, source_bound_plan(condition_id, source_id), workspace);
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
        workspace,
        channel_mask);
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
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask) {
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
    return standard_leaf_channels_mask(
        desc.dist_kind,
        loaded.params.data(),
        std::min(desc.param_count, 8),
        loaded.q,
        loaded.t0,
        t - desc.onset_abs_value,
        channel_mask);
  }

  leaf::EventChannels evaluate_leaf_onset_convolution_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask) {
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

    const bool need_pdf = (channel_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (channel_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    const std::uint8_t shifted_mask =
        (need_pdf ? kLeafChannelPdf : 0U) |
        (need_cdf ? kLeafChannelCdf : 0U);

    auto shifted_channels = [&](const double source_time) {
      return standard_leaf_channels_mask(
          desc.dist_kind,
          loaded.params.data(),
          std::min(desc.param_count, 8),
          loaded.q,
          loaded.t0,
          t - source_time - lag,
          shifted_mask);
    };

    const auto onset_source_id = kernel.onset_source_id;
    const auto onset_bounds =
        compiled_source_bounds(
            onset_source_id,
            source_bound_plan(condition_id, onset_source_id),
            workspace);
    if (const double *exact_time = onset_bounds.exact_time) {
      const auto shifted = shifted_channels(*exact_time);
      leaf::EventChannels out;
      if (need_pdf) {
        out.pdf = shifted.pdf;
      }
      if ((channel_mask & kLeafChannelCdf) != 0U) {
        out.cdf = shifted.cdf;
      }
      if ((channel_mask & kLeafChannelSurvival) != 0U) {
        out.survival = clamp_probability(1.0 - shifted.cdf);
      }
      return out;
    }

    const auto batch = quadrature::build_finite_batch(0.0, upper);
    double pdf = 0.0;
    double cdf = 0.0;
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double u = batch.nodes.nodes[i];
      const auto frame = base_time_guard(u);
      const auto source =
          base_source(
              onset_source_id,
              condition_id,
              workspace,
              kLeafChannelPdf);
      if (!(source.pdf > 0.0)) {
        continue;
      }
      const auto shifted = shifted_channels(u);
      const double weight = batch.nodes.weights[i] * source.pdf;
      if (need_pdf) {
        pdf += weight * shifted.pdf;
      }
      if (need_cdf) {
        cdf += weight * shifted.cdf;
      }
    }

    leaf::EventChannels out;
    if (need_pdf) {
      out.pdf = safe_density(pdf);
    }
    if (need_cdf) {
      const double clamped = clamp_probability(cdf);
      if ((channel_mask & kLeafChannelCdf) != 0U) {
        out.cdf = clamped;
      }
      if ((channel_mask & kLeafChannelSurvival) != 0U) {
        out.survival = clamp_probability(1.0 - clamped);
      }
    }
    return out;
  }

  leaf::EventChannels evaluate_pool_kofn_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask) {
    const auto begin = kernel.pool_member_offset;
    const auto end = begin + kernel.pool_member_count;
    const auto k = kernel.pool_k;
    if (end > static_cast<semantic::Index>(
                  program_.pool_member_source_ids.size())) {
      return impossible_channels();
    }
    const auto n_members = static_cast<int>(end - begin);
    const bool need_pdf = (channel_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (channel_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    const std::uint8_t member_mask =
        kLeafChannelCdf | kLeafChannelSurvival |
        (need_pdf ? kLeafChannelPdf : 0U);

    std::vector<leaf::EventChannels> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(base_source(
          program_.pool_member_source_ids[static_cast<std::size_t>(i)],
          condition_id,
          workspace,
          member_mask));
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

    double pdf = 0.0;
    if (need_pdf) {
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
    }
    double surv = 0.0;
    if (need_cdf) {
      for (int m = 0; m < k; ++m) {
        surv += prefix[idx(n_members, m)];
      }
    }

    leaf::EventChannels out;
    if (need_pdf) {
      out.pdf = safe_density(pdf);
    }
    if (need_cdf) {
      const double clamped_survival = clamp_probability(surv);
      if ((channel_mask & kLeafChannelSurvival) != 0U) {
        out.survival = clamped_survival;
      }
      if ((channel_mask & kLeafChannelCdf) != 0U) {
        out.cdf = clamp_probability(1.0 - clamped_survival);
      }
    }
    return out;
  }

  SourceProductScalarFill evaluate_source_product_leaf_absolute_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask) {
    (void)condition_id;
    (void)workspace;
    const auto pos = static_cast<std::size_t>(kernel.leaf_index);
    if (kernel.leaf_index == semantic::kInvalidIndex ||
        pos >= program_.leaf_descriptors.size() ||
        pos >= leaf_inputs_.size()) {
      return source_product_impossible_fill(fill_mask);
    }
    const auto &desc = program_.leaf_descriptors[pos];
    return evaluate_source_product_leaf_fill_at(
        kernel.leaf_index,
        active_source_product_base_time() - desc.onset_abs_value,
        fill_mask);
  }

  SourceProductScalarFill evaluate_source_product_leaf_onset_convolution_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask) {
    const auto pos = static_cast<std::size_t>(kernel.leaf_index);
    if (kernel.leaf_index == semantic::kInvalidIndex ||
        pos >= program_.leaf_descriptors.size() ||
        pos >= leaf_inputs_.size()) {
      return source_product_impossible_fill(fill_mask);
    }
    const auto &desc = program_.leaf_descriptors[pos];
    const double t = active_source_product_base_time();
    const double lag = desc.onset_lag;
    const double upper = t - lag;
    if (!(upper > 0.0)) {
      return source_product_impossible_fill(fill_mask);
    }

    const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    const std::uint8_t shifted_mask =
        (need_pdf ? kLeafChannelPdf : 0U) |
        (need_cdf ? kLeafChannelCdf : 0U);
    auto shifted_fill = [&](const double source_time) {
      return evaluate_source_product_leaf_fill_at(
          kernel.leaf_index,
          t - source_time - lag,
          shifted_mask);
    };

    const auto onset_source_id = kernel.onset_source_id;
    const auto onset_bounds =
        compiled_source_bounds(
            onset_source_id,
            source_bound_plan(condition_id, onset_source_id),
            workspace);
    if (const double *exact_time = onset_bounds.exact_time) {
      const auto shifted = shifted_fill(*exact_time);
      SourceProductScalarFill out;
      out.mask = fill_mask;
      if (need_pdf) {
        out.pdf = shifted.pdf;
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        out.cdf = shifted.cdf;
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        out.survival = clamp_probability(1.0 - shifted.cdf);
      }
      return out;
    }

    const auto batch = quadrature::build_finite_batch(0.0, upper);
    double pdf = 0.0;
    double cdf = 0.0;
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double u = batch.nodes.nodes[i];
      const auto frame = source_product_base_time_guard(u);
      const auto source =
          source_product_base_fill(
              onset_source_id,
              condition_id,
              workspace,
              kLeafChannelPdf);
      if (!(source.pdf > 0.0)) {
        continue;
      }
      const auto shifted = shifted_fill(u);
      const double weight = batch.nodes.weights[i] * source.pdf;
      if (need_pdf) {
        pdf += weight * shifted.pdf;
      }
      if (need_cdf) {
        cdf += weight * shifted.cdf;
      }
    }

    SourceProductScalarFill out;
    out.mask = fill_mask;
    if (need_pdf) {
      out.pdf = safe_density(pdf);
    }
    if (need_cdf) {
      const double clamped = clamp_probability(cdf);
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        out.cdf = clamped;
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        out.survival = clamp_probability(1.0 - clamped);
      }
    }
    return out;
  }

  SourceProductScalarFill evaluate_source_product_pool_kofn_kernel(
      const ExactSourceKernel &kernel,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask) {
    const auto begin = kernel.pool_member_offset;
    const auto end = begin + kernel.pool_member_count;
    const auto k = kernel.pool_k;
    if (end > static_cast<semantic::Index>(
                  program_.pool_member_source_ids.size())) {
      return source_product_impossible_fill(fill_mask);
    }
    const auto n_members = static_cast<int>(end - begin);
    const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
    const bool need_cdf =
        (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
    const std::uint8_t member_mask =
        kLeafChannelCdf | kLeafChannelSurvival |
        (need_pdf ? kLeafChannelPdf : 0U);

    std::vector<SourceProductScalarFill> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(source_product_base_fill(
          program_.pool_member_source_ids[static_cast<std::size_t>(i)],
          condition_id,
          workspace,
          member_mask));
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

    double pdf = 0.0;
    if (need_pdf) {
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
    }
    double surv = 0.0;
    if (need_cdf) {
      for (int m = 0; m < k; ++m) {
        surv += prefix[idx(n_members, m)];
      }
    }

    SourceProductScalarFill out;
    out.mask = fill_mask;
    if (need_pdf) {
      out.pdf = safe_density(pdf);
    }
    if (need_cdf) {
      const double clamped_survival = clamp_probability(surv);
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        out.survival = clamped_survival;
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        out.cdf = clamp_probability(1.0 - clamped_survival);
      }
    }
    return out;
  }

  leaf::EventChannels load_conditional_source(
      TimedCache *cache,
      const double time,
      const CompiledSourceChannelPlan &channel_plan,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t channel_mask) {
    const auto source_id = channel_plan.request.source_id;
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->has(pos, channel_mask)) {
      return cache->get(pos);
    }
    const auto channels =
        evaluate_conditioned_source_kernel(
            channel_plan, time, workspace, channel_mask);
    cache->store(pos, channels, channel_mask);
    return cache->get(pos);
  }

  SourceProductScalarFill load_source_product_conditioned_fill(
      SourceProductTimedCache *cache,
      const double time,
      const CompiledMathSourceProductChannel &channel,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask) {
    const auto source_id = channel.source_id;
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return source_product_impossible_fill(fill_mask);
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->has(pos, fill_mask)) {
      return cache->get(pos);
    }
    const auto fill =
        evaluate_source_product_conditioned_fill(
            channel, time, workspace, fill_mask);
    cache->store(pos, fill, fill_mask);
    return cache->get(pos);
  }

  SourceProductScalarFill load_source_product_base_fill(
      SourceProductTimedCache *cache,
      const double time,
      const semantic::Index source_id,
      const semantic::Index condition_id,
      const CompiledMathWorkspace *workspace,
      const std::uint8_t fill_mask) {
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return source_product_impossible_fill(fill_mask);
    }
    cache->set_context(time, condition_cache_id(condition_id, workspace));
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->has(pos, fill_mask)) {
      return cache->get(pos);
    }
    const bool needs_frame = active_source_product_base_cache() != cache;
    if (needs_frame) {
      source_product_base_cache_stack_.push_back(cache);
    }
    const auto fill =
        compute_source_product_base_fill(
            source_id, condition_id, workspace, fill_mask);
    cache->store(pos, fill, fill_mask);
    if (needs_frame) {
      source_product_base_cache_stack_.pop_back();
    }
    return cache->get(pos);
  }

  leaf::EventChannels load_base_source(TimedCache *cache,
                                       const double time,
                                       const semantic::Index source_id,
                                       const semantic::Index condition_id,
                                       const CompiledMathWorkspace *workspace,
                                       const std::uint8_t channel_mask) {
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    cache->set_context(time, condition_cache_id(condition_id, workspace));
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->has(pos, channel_mask)) {
      return cache->get(pos);
    }
    const bool needs_frame = active_base_cache() != cache;
    if (needs_frame) {
      base_cache_stack_.push_back(cache);
    }
    const auto channels =
        compute_base_source(source_id, condition_id, workspace, channel_mask);
    cache->store(pos, channels, channel_mask);
    if (needs_frame) {
      base_cache_stack_.pop_back();
    }
    return cache->get(pos);
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
