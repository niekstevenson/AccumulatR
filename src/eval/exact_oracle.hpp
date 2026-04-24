#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

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

class ExactOracleScratch {
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
      invalidate();
    }

    void set_time(const double t) {
      if (std::fabs(t - time) <= 1e-12) {
        return;
      }
      reset(t);
    }

    double time{0.0};
    std::vector<leaf::EventChannels> channels;
    std::vector<std::uint32_t> epoch;
    std::uint32_t current_epoch{1U};
  };

  class ConditionalTimeGuard {
  public:
    ConditionalTimeGuard(ExactOracleScratch *oracle, const double time)
        : oracle_(oracle) {
      if (oracle_ == nullptr) {
        return;
      }
      previous_time_ = oracle_->conditional_time_;
      changed_ = oracle_->set_conditional_time(time);
    }

    ConditionalTimeGuard(const ConditionalTimeGuard &) = delete;
    ConditionalTimeGuard &operator=(const ConditionalTimeGuard &) = delete;

    ~ConditionalTimeGuard() {
      if (oracle_ != nullptr && changed_) {
        oracle_->set_conditional_time(previous_time_);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    double previous_time_{0.0};
    bool changed_{false};
  };

  class BaseTimeGuard {
  public:
    BaseTimeGuard(ExactOracleScratch *oracle,
                  semantic::Index recursive_index,
                  TimedCache *cache,
                  const double time)
        : oracle_(oracle), recursive_index(recursive_index) {
      if (oracle_ == nullptr || cache == nullptr) {
        return;
      }
      oracle_->push_base_cache(cache, time);
      active_ = true;
    }

    BaseTimeGuard(const BaseTimeGuard &) = delete;
    BaseTimeGuard &operator=(const BaseTimeGuard &) = delete;

    ~BaseTimeGuard() {
      if (oracle_ != nullptr && active_) {
        oracle_->pop_base_cache(recursive_index);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    semantic::Index recursive_index{-1};
    bool active_{false};
  };

  explicit ExactOracleScratch(const ExactVariantPlan &plan)
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

  [[nodiscard]] ConditionalTimeGuard conditional_time_guard(const double time) {
    return ConditionalTimeGuard(this, time);
  }

  [[nodiscard]] BaseTimeGuard base_time_guard(const double time) {
    const auto recursive_index = acquire_recursive_cache();
    return BaseTimeGuard(
        this,
        recursive_index,
        &recursive_base_caches_[static_cast<std::size_t>(recursive_index)],
        time);
  }

  double conditional_time() const noexcept {
    return conditional_time_;
  }

  leaf::EventChannels conditional_source(const semantic::Index source_id) {
    conditional_cache_.set_time(conditional_time_);
    return load_conditional_source(&conditional_cache_, source_id);
  }

  leaf::EventChannels source_channels(const semantic::Index source_id,
                                      const double t) {
    const auto guard = conditional_time_guard(t);
    return conditional_source(source_id);
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
  semantic::Index recursive_cache_depth_{0};

  bool set_conditional_time(const double t) {
    if (std::fabs(t - conditional_time_) <= 1e-12) {
      return false;
    }
    conditional_time_ = t;
    conditional_cache_.set_time(t);
    return true;
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

  const double *exact_time_for(const semantic::Index source_id) const {
    if (source_id == semantic::kInvalidIndex) {
      return nullptr;
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos >= sequence_state_->exact_times.size()) {
      return nullptr;
    }
    const auto value = sequence_state_->exact_times[pos];
    if (!std::isfinite(value)) {
      return nullptr;
    }
    return &sequence_state_->exact_times[pos];
  }

  leaf::EventChannels conditionalize(const leaf::EventChannels uncond,
                                     const leaf::EventChannels lower) const {
    if (!(sequence_state_->lower_bound > 0.0)) {
      return uncond;
    }
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

  leaf::EventChannels compute_conditional_source(const semantic::Index source_id) {
    if (const double *exact_time = exact_time_for(source_id)) {
      if (!(conditional_time_ >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }

    const auto uncond =
        load_base_source(&base_current_cache_, conditional_time_, source_id);
    if (!(sequence_state_->lower_bound > 0.0)) {
      return uncond;
    }

    const auto lower = load_base_source(
        &base_lower_cache_, sequence_state_->lower_bound, source_id);
    return conditionalize(uncond, lower);
  }

  leaf::EventChannels base_source(const semantic::Index source_id) {
    if (auto *cache = active_base_cache(); cache != nullptr) {
      return load_base_source(cache, cache->time, source_id);
    }
    return load_base_source(&base_current_cache_, conditional_time_, source_id);
  }

  leaf::EventChannels compute_base_source(const semantic::Index source_id) {
    if (const double *exact_time = exact_time_for(source_id)) {
      if (!(active_base_time() >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }
    if (source_id < program_.layout.n_leaves) {
      return base_leaf_channels(source_id);
    }
    return base_pool_channels(source_id - program_.layout.n_leaves);
  }

  leaf::EventChannels base_leaf_channels(const semantic::Index index) {
    const auto pos = static_cast<std::size_t>(index);
    const auto &desc = program_.leaf_descriptors[pos];
    const auto &loaded = leaf_inputs_[pos];
    const auto onset_kind = static_cast<semantic::OnsetKind>(desc.onset_kind);
    const double t = active_base_time();
    if (onset_kind == semantic::OnsetKind::Absolute) {
      return standard_leaf_channels(
          desc.dist_kind,
          loaded.params.data(),
          std::min(desc.param_count, 8),
          loaded.q,
          loaded.t0,
          t - desc.onset_abs_value);
    }

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

    const auto onset_source_id = desc.onset_source_id;
    if (const double *exact_time = exact_time_for(onset_source_id)) {
      return shifted_channels(*exact_time);
    }

    const auto batch = quadrature::build_finite_batch(0.0, upper);
    double pdf = 0.0;
    double cdf = 0.0;
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double u = batch.nodes.nodes[i];
      const auto frame = base_time_guard(u);
      const auto source = base_source(onset_source_id);
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

  leaf::EventChannels base_pool_channels(const semantic::Index index) {
    const auto pos = static_cast<std::size_t>(index);
    const auto begin = program_.pool_member_offsets[pos];
    const auto end = program_.pool_member_offsets[pos + 1U];
    const auto k = program_.pool_k[pos];
    const auto n_members = static_cast<int>(end - begin);

    std::vector<leaf::EventChannels> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(base_source(
          program_.pool_member_source_ids[static_cast<std::size_t>(i)]));
    }

    std::vector<std::vector<double>> prefix(
        static_cast<std::size_t>(n_members + 1),
        std::vector<double>(static_cast<std::size_t>(n_members + 1), 0.0));
    std::vector<std::vector<double>> suffix(
        static_cast<std::size_t>(n_members + 1),
        std::vector<double>(static_cast<std::size_t>(n_members + 1), 0.0));
    prefix[0][0] = 1.0;
    for (int i = 0; i < n_members; ++i) {
      for (int m = 0; m <= i; ++m) {
        prefix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].survival;
        prefix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m + 1)] +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].cdf;
      }
    }
    suffix[static_cast<std::size_t>(n_members)][0] = 1.0;
    for (int i = n_members - 1; i >= 0; --i) {
      const int count = n_members - i - 1;
      for (int m = 0; m <= count; ++m) {
        suffix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] +=
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].survival;
        suffix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m + 1)] +=
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].cdf;
      }
    }

    double surv = 0.0;
    for (int m = 0; m < k; ++m) {
      surv += prefix[static_cast<std::size_t>(n_members)][static_cast<std::size_t>(m)];
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
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(left)] *
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(right)];
      }
      pdf += members[static_cast<std::size_t>(i)].pdf * others_exact;
    }

    leaf::EventChannels out;
    out.pdf = safe_density(pdf);
    out.survival = clamp_probability(surv);
    out.cdf = clamp_probability(1.0 - out.survival);
    return out;
  }

  leaf::EventChannels load_conditional_source(TimedCache *cache,
                                              const semantic::Index source_id) {
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->epoch[pos] == cache->current_epoch) {
      return cache->channels[pos];
    }
    cache->channels[pos] = compute_conditional_source(source_id);
    cache->epoch[pos] = cache->current_epoch;
    return cache->channels[pos];
  }

  leaf::EventChannels load_base_source(TimedCache *cache,
                                       const double time,
                                       const semantic::Index source_id) {
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    cache->set_time(time);
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->epoch[pos] == cache->current_epoch) {
      return cache->channels[pos];
    }
    const bool needs_frame = active_base_cache() != cache;
    if (needs_frame) {
      base_cache_stack_.push_back(cache);
    }
    cache->channels[pos] = compute_base_source(source_id);
    cache->epoch[pos] = cache->current_epoch;
    if (needs_frame) {
      base_cache_stack_.pop_back();
    }
    return cache->channels[pos];
  }
};

using ExactSourceOracle = ExactOracleScratch;

inline std::vector<ExactTriggerState> enumerate_trigger_states(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row) {
  std::vector<ExactTriggerState> states(1);
  states.front().weight = 1.0;
  states.front().shared_started.assign(
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

    std::vector<ExactTriggerState> next;
    next.reserve(states.size() * 2U);
    for (const auto &state : states) {
      if (q > 0.0) {
        auto fail = state;
        fail.weight *= q;
        fail.shared_started[pos] = 0U;
        next.push_back(std::move(fail));
      }
      if (q < 1.0) {
        auto start = state;
        start.weight *= (1.0 - q);
        start.shared_started[pos] = 1U;
        next.push_back(std::move(start));
      }
    }
    states.swap(next);
  }

  return states;
}

inline double source_cdf_at(ExactSourceOracle *oracle,
                            const semantic::Index source_id,
                            const double t) {
  const auto guard = oracle->conditional_time_guard(t);
  return clamp_probability(oracle->conditional_source(source_id).cdf);
}

inline double source_pdf_at(ExactSourceOracle *oracle,
                            const semantic::Index source_id,
                            const double t) {
  const auto guard = oracle->conditional_time_guard(t);
  return safe_density(oracle->conditional_source(source_id).pdf);
}

inline double source_survival_at(ExactSourceOracle *oracle,
                                 const semantic::Index source_id,
                                 const double t) {
  const auto guard = oracle->conditional_time_guard(t);
  return clamp_probability(oracle->conditional_source(source_id).survival);
}

} // namespace detail
} // namespace accumulatr::eval
