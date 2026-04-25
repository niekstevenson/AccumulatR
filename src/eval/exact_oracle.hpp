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

  class ExactTimeOverlayGuard {
  public:
    ExactTimeOverlayGuard(ExactOracleScratch *oracle,
                          const semantic::Index source_id,
                          const double time)
        : oracle_(oracle), source_id_(source_id) {
      if (oracle_ == nullptr || source_id == semantic::kInvalidIndex ||
          !std::isfinite(time)) {
        return;
      }
      previous_time_ = oracle_->set_exact_time_overlay(source_id, time);
      active_ = true;
    }

    ExactTimeOverlayGuard(const ExactTimeOverlayGuard &) = delete;
    ExactTimeOverlayGuard &operator=(const ExactTimeOverlayGuard &) = delete;

    ~ExactTimeOverlayGuard() {
      if (oracle_ != nullptr && active_) {
        oracle_->set_exact_time_overlay(source_id_, previous_time_);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    semantic::Index source_id_{semantic::kInvalidIndex};
    double previous_time_{std::numeric_limits<double>::quiet_NaN()};
    bool active_{false};
  };

  class SourceLowerBoundOverlayGuard {
  public:
    SourceLowerBoundOverlayGuard(ExactOracleScratch *oracle,
                                 const semantic::Index source_id,
                                 const double time)
        : oracle_(oracle), source_id_(source_id) {
      if (oracle_ == nullptr || source_id == semantic::kInvalidIndex ||
          !std::isfinite(time)) {
        return;
      }
      previous_time_ = oracle_->set_source_lower_bound_overlay(source_id, time);
      active_ = true;
    }

    SourceLowerBoundOverlayGuard(const SourceLowerBoundOverlayGuard &) = delete;
    SourceLowerBoundOverlayGuard &operator=(
        const SourceLowerBoundOverlayGuard &) = delete;

    ~SourceLowerBoundOverlayGuard() {
      if (oracle_ != nullptr && active_) {
        oracle_->set_source_lower_bound_overlay(source_id_, previous_time_);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    semantic::Index source_id_{semantic::kInvalidIndex};
    double previous_time_{std::numeric_limits<double>::quiet_NaN()};
    bool active_{false};
  };

  class SourceUpperBoundOverlayGuard {
  public:
    SourceUpperBoundOverlayGuard(ExactOracleScratch *oracle,
                                 const semantic::Index source_id,
                                 const double time)
        : oracle_(oracle), source_id_(source_id) {
      if (oracle_ == nullptr || source_id == semantic::kInvalidIndex ||
          !std::isfinite(time)) {
        return;
      }
      previous_time_ = oracle_->set_source_upper_bound_overlay(source_id, time);
      active_ = true;
    }

    SourceUpperBoundOverlayGuard(const SourceUpperBoundOverlayGuard &) = delete;
    SourceUpperBoundOverlayGuard &operator=(
        const SourceUpperBoundOverlayGuard &) = delete;

    ~SourceUpperBoundOverlayGuard() {
      if (oracle_ != nullptr && active_) {
        oracle_->set_source_upper_bound_overlay(source_id_, previous_time_);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    semantic::Index source_id_{semantic::kInvalidIndex};
    double previous_time_{std::numeric_limits<double>::quiet_NaN()};
    bool active_{false};
  };

  struct ExprUpperBoundOverlayState {
    double time{std::numeric_limits<double>::quiet_NaN()};
    double normalizer{std::numeric_limits<double>::quiet_NaN()};
  };

  class ExprUpperBoundOverlayGuard {
  public:
    ExprUpperBoundOverlayGuard(ExactOracleScratch *oracle,
                               const semantic::Index expr_id,
                               const double time,
                               const double normalizer)
        : oracle_(oracle), expr_id_(expr_id) {
      if (oracle_ == nullptr || expr_id == semantic::kInvalidIndex ||
          !std::isfinite(time) || !(normalizer > 0.0)) {
        return;
      }
      previous_ =
          oracle_->set_expr_upper_bound_overlay(expr_id, time, normalizer);
      active_ = true;
    }

    ExprUpperBoundOverlayGuard(const ExprUpperBoundOverlayGuard &) = delete;
    ExprUpperBoundOverlayGuard &operator=(
        const ExprUpperBoundOverlayGuard &) = delete;

    ~ExprUpperBoundOverlayGuard() {
      if (oracle_ != nullptr && active_) {
        oracle_->set_expr_upper_bound_overlay(
            expr_id_, previous_.time, previous_.normalizer);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    semantic::Index expr_id_{semantic::kInvalidIndex};
    ExprUpperBoundOverlayState previous_{};
    bool active_{false};
  };

  class SourceOrderOverlayGuard {
  public:
    SourceOrderOverlayGuard(ExactOracleScratch *oracle,
                            const std::vector<ExactSourceOrderFact> &facts)
        : oracle_(oracle) {
      if (oracle_ == nullptr || facts.empty()) {
        return;
      }
      previous_size_ = oracle_->source_order_overlays_.size();
      for (const auto &fact : facts) {
        if (fact.before_source_id == semantic::kInvalidIndex ||
            fact.after_source_id == semantic::kInvalidIndex ||
            fact.before_source_id == fact.after_source_id) {
          continue;
        }
        oracle_->source_order_overlays_.push_back(fact);
      }
      active_ = oracle_->source_order_overlays_.size() != previous_size_;
    }

    SourceOrderOverlayGuard(const SourceOrderOverlayGuard &) = delete;
    SourceOrderOverlayGuard &operator=(const SourceOrderOverlayGuard &) = delete;

    ~SourceOrderOverlayGuard() {
      if (oracle_ != nullptr && active_) {
        oracle_->source_order_overlays_.resize(previous_size_);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    std::size_t previous_size_{0U};
    bool active_{false};
  };

  struct GuardUpperBoundOverlay {
    semantic::Index ref_source_id{semantic::kInvalidIndex};
    semantic::Index blocker_source_id{semantic::kInvalidIndex};
    double time{std::numeric_limits<double>::quiet_NaN()};
    double normalizer{0.0};
  };

  class GuardUpperBoundOverlayGuard {
  public:
    GuardUpperBoundOverlayGuard(ExactOracleScratch *oracle,
                                const semantic::Index ref_source_id,
                                const semantic::Index blocker_source_id,
                                const double time,
                                const double normalizer)
        : oracle_(oracle) {
      if (oracle_ == nullptr || ref_source_id == semantic::kInvalidIndex ||
          blocker_source_id == semantic::kInvalidIndex || !std::isfinite(time) ||
          !(normalizer > 0.0)) {
        return;
      }
      previous_size_ = oracle_->guard_upper_bound_overlays_.size();
      oracle_->guard_upper_bound_overlays_.push_back(GuardUpperBoundOverlay{
          ref_source_id, blocker_source_id, time, normalizer});
      active_ = true;
    }

    GuardUpperBoundOverlayGuard(const GuardUpperBoundOverlayGuard &) = delete;
    GuardUpperBoundOverlayGuard &operator=(
        const GuardUpperBoundOverlayGuard &) = delete;

    ~GuardUpperBoundOverlayGuard() {
      if (oracle_ != nullptr && active_) {
        oracle_->guard_upper_bound_overlays_.resize(previous_size_);
      }
    }

  private:
    ExactOracleScratch *oracle_{nullptr};
    std::size_t previous_size_{0U};
    bool active_{false};
  };

  explicit ExactOracleScratch(const ExactVariantPlan &plan)
      : plan_(plan),
        program_(plan.lowered.program),
        leaf_inputs_(static_cast<std::size_t>(program_.layout.n_leaves)),
        conditional_cache_(static_cast<std::size_t>(plan.source_count)),
        base_current_cache_(static_cast<std::size_t>(plan.source_count)),
        base_lower_cache_(static_cast<std::size_t>(plan.source_count)),
        exact_time_overlays_(
            static_cast<std::size_t>(plan.source_count),
            std::numeric_limits<double>::quiet_NaN()),
        source_lower_bound_overlays_(
            static_cast<std::size_t>(plan.source_count),
            std::numeric_limits<double>::quiet_NaN()),
        source_upper_bound_overlays_(
            static_cast<std::size_t>(plan.source_count),
            std::numeric_limits<double>::quiet_NaN()),
        expr_upper_bound_overlays_(program_.expr_kind.size(),
                                   std::numeric_limits<double>::quiet_NaN()),
        expr_upper_bound_normalizers_(program_.expr_kind.size(),
            std::numeric_limits<double>::quiet_NaN()) {
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
    std::fill(
        exact_time_overlays_.begin(),
        exact_time_overlays_.end(),
        std::numeric_limits<double>::quiet_NaN());
    std::fill(
        source_lower_bound_overlays_.begin(),
        source_lower_bound_overlays_.end(),
        std::numeric_limits<double>::quiet_NaN());
    std::fill(
        source_upper_bound_overlays_.begin(),
        source_upper_bound_overlays_.end(),
        std::numeric_limits<double>::quiet_NaN());
    std::fill(
        expr_upper_bound_overlays_.begin(),
        expr_upper_bound_overlays_.end(),
        std::numeric_limits<double>::quiet_NaN());
    std::fill(
        expr_upper_bound_normalizers_.begin(),
        expr_upper_bound_normalizers_.end(),
        std::numeric_limits<double>::quiet_NaN());
    source_order_overlays_.clear();
    guard_upper_bound_overlays_.clear();
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

  leaf::EventChannels conditional_source(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) {
    conditional_cache_.set_context(
        conditional_time_, condition_frame_id(condition_frame));
    return load_conditional_source(
        &conditional_cache_, source_id, condition_frame);
  }

  leaf::EventChannels source_channels(const semantic::Index source_id,
                                      const double t,
                                      const ExactConditionFrame *condition_frame = nullptr) {
    const auto guard = conditional_time_guard(t);
    return conditional_source(source_id, condition_frame);
  }

  [[nodiscard]] ExactTimeOverlayGuard exact_time_overlay(
      const semantic::Index source_id,
      const double time) {
    return ExactTimeOverlayGuard(this, source_id, time);
  }

  [[nodiscard]] SourceLowerBoundOverlayGuard source_lower_bound_overlay(
      const semantic::Index source_id,
      const double time) {
    return SourceLowerBoundOverlayGuard(this, source_id, time);
  }

  [[nodiscard]] SourceUpperBoundOverlayGuard source_upper_bound_overlay(
      const semantic::Index source_id,
      const double time) {
    return SourceUpperBoundOverlayGuard(this, source_id, time);
  }

  [[nodiscard]] ExprUpperBoundOverlayGuard expr_upper_bound_overlay(
      const semantic::Index expr_id,
      const double time,
      const double normalizer) {
    return ExprUpperBoundOverlayGuard(this, expr_id, time, normalizer);
  }

  [[nodiscard]] SourceOrderOverlayGuard source_order_overlay(
      const std::vector<ExactSourceOrderFact> &facts) {
    return SourceOrderOverlayGuard(this, facts);
  }

  [[nodiscard]] GuardUpperBoundOverlayGuard guard_upper_bound_overlay(
      const semantic::Index ref_source_id,
      const semantic::Index blocker_source_id,
      const double time,
      const double normalizer) {
    return GuardUpperBoundOverlayGuard(
        this, ref_source_id, blocker_source_id, time, normalizer);
  }

  bool has_exact_time_overlay(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (source_id == semantic::kInvalidIndex) {
      return false;
    }
    if (condition_frame != nullptr) {
      const auto pos = static_cast<std::size_t>(source_id);
      if (pos < condition_frame->source_exact_times.size() &&
          std::isfinite(condition_frame->source_exact_times[pos])) {
        return true;
      }
    }
    const auto pos = static_cast<std::size_t>(source_id);
    return pos < exact_time_overlays_.size() &&
           std::isfinite(exact_time_overlays_[pos]);
  }

  bool has_source_condition_overlay(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (source_id == semantic::kInvalidIndex) {
      return false;
    }
    if (condition_frame != nullptr) {
      const auto pos = static_cast<std::size_t>(source_id);
      if ((pos < condition_frame->source_exact_times.size() &&
           std::isfinite(condition_frame->source_exact_times[pos])) ||
          (pos < condition_frame->source_lower_bounds.size() &&
           std::isfinite(condition_frame->source_lower_bounds[pos])) ||
          (pos < condition_frame->source_upper_bounds.size() &&
           std::isfinite(condition_frame->source_upper_bounds[pos]))) {
        return true;
      }
    }
    const auto pos = static_cast<std::size_t>(source_id);
    return (pos < exact_time_overlays_.size() &&
            std::isfinite(exact_time_overlays_[pos])) ||
           (pos < source_lower_bound_overlays_.size() &&
            std::isfinite(source_lower_bound_overlays_[pos])) ||
           (pos < source_upper_bound_overlays_.size() &&
            std::isfinite(source_upper_bound_overlays_[pos]));
  }

  const double *exact_time_for_source(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    return exact_time_for(source_id, condition_frame);
  }

  const ExactTimedExprUpperBound *expr_upper_bound_for(
      const semantic::Index expr_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (expr_id == semantic::kInvalidIndex) {
      return nullptr;
    }
    if (condition_frame != nullptr && condition_frame->has_expr_upper_bounds) {
      const ExactTimedExprUpperBound *best = nullptr;
      for (const auto *frame = condition_frame; frame != nullptr;
           frame = frame->parent) {
        for (const auto &fact : frame->expr_upper_bounds) {
          if (fact.expr_id == expr_id && std::isfinite(fact.time) &&
              fact.normalizer > 0.0 &&
              (best == nullptr || fact.time < best->time)) {
            best = &fact;
          }
        }
      }
      if (best != nullptr) {
        return best;
      }
    }
    const auto pos = static_cast<std::size_t>(expr_id);
    if (pos >= expr_upper_bound_overlays_.size() ||
        !std::isfinite(expr_upper_bound_overlays_[pos])) {
      return nullptr;
    }
    legacy_expr_upper_bound_ = ExactTimedExprUpperBound{
        expr_id,
        expr_upper_bound_overlays_[pos],
        expr_upper_bound_normalizers_[pos]};
    return &legacy_expr_upper_bound_;
  }

  double expr_upper_bound_normalizer(const semantic::Index expr_id) const {
    if (expr_id == semantic::kInvalidIndex) {
      return 0.0;
    }
    const auto pos = static_cast<std::size_t>(expr_id);
    if (pos >= expr_upper_bound_normalizers_.size()) {
      return 0.0;
    }
    return expr_upper_bound_normalizers_[pos];
  }

  bool source_order_known_before(
      const semantic::Index before_source_id,
      const semantic::Index after_source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (before_source_id == semantic::kInvalidIndex ||
        after_source_id == semantic::kInvalidIndex ||
        before_source_id == after_source_id) {
      return false;
    }
    if (condition_frame != nullptr && condition_frame->has_source_order_facts) {
      for (const auto *frame = condition_frame; frame != nullptr;
           frame = frame->parent) {
        for (const auto &fact : frame->source_order_facts) {
          if (fact.before_source_id == before_source_id &&
              fact.after_source_id == after_source_id) {
            return true;
          }
        }
      }
    }
    for (const auto &fact : source_order_overlays_) {
      if (fact.before_source_id == before_source_id &&
          fact.after_source_id == after_source_id) {
        return true;
      }
    }
    return false;
  }

  const ExactTimedGuardUpperBound *guard_upper_bound_for(
      const semantic::Index ref_source_id,
      const semantic::Index blocker_source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (condition_frame != nullptr && condition_frame->has_guard_upper_bounds) {
      const ExactTimedGuardUpperBound *best = nullptr;
      for (const auto *frame = condition_frame; frame != nullptr;
           frame = frame->parent) {
        for (const auto &fact : frame->guard_upper_bounds) {
          if (fact.ref_source_id == ref_source_id &&
              fact.blocker_source_id == blocker_source_id &&
              std::isfinite(fact.time) && fact.normalizer > 0.0 &&
              (best == nullptr || fact.time < best->time)) {
            best = &fact;
          }
        }
      }
      if (best != nullptr) {
        return best;
      }
    }
    for (auto it = guard_upper_bound_overlays_.rbegin();
         it != guard_upper_bound_overlays_.rend();
         ++it) {
      if (it->ref_source_id == ref_source_id &&
          it->blocker_source_id == blocker_source_id) {
        legacy_guard_upper_bound_ = ExactTimedGuardUpperBound{
            it->ref_source_id,
            it->blocker_source_id,
            it->time,
            it->normalizer};
        return &legacy_guard_upper_bound_;
      }
    }
    return nullptr;
  }

  bool has_guard_upper_bound_overlay(
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (condition_frame != nullptr &&
        condition_frame->has_guard_upper_bounds) {
      return true;
    }
    return !guard_upper_bound_overlays_.empty();
  }

  bool has_source_order_overlay(
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (condition_frame != nullptr &&
        condition_frame->has_source_order_facts) {
      return true;
    }
    return !source_order_overlays_.empty();
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
  std::vector<double> exact_time_overlays_;
  std::vector<double> source_lower_bound_overlays_;
  std::vector<double> source_upper_bound_overlays_;
  std::vector<double> expr_upper_bound_overlays_;
  std::vector<double> expr_upper_bound_normalizers_;
  std::vector<ExactSourceOrderFact> source_order_overlays_;
  std::vector<GuardUpperBoundOverlay> guard_upper_bound_overlays_;
  mutable ExactTimedExprUpperBound legacy_expr_upper_bound_;
  mutable ExactTimedGuardUpperBound legacy_guard_upper_bound_;
  semantic::Index recursive_cache_depth_{0};

  static semantic::Index condition_frame_id(
      const ExactConditionFrame *condition_frame) {
    return condition_frame == nullptr ? 0 : condition_frame->id;
  }

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

  void invalidate_all_timed_caches() {
    conditional_cache_.invalidate();
    base_current_cache_.invalidate();
    base_lower_cache_.invalidate();
    for (auto &cache : recursive_base_caches_) {
      cache.invalidate();
    }
  }

  double set_exact_time_overlay(const semantic::Index source_id,
                                const double time) {
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos >= exact_time_overlays_.size()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double previous = exact_time_overlays_[pos];
    if ((std::isfinite(previous) && std::isfinite(time) &&
         std::fabs(previous - time) <= 1e-12) ||
        (!std::isfinite(previous) && !std::isfinite(time))) {
      return previous;
    }
    exact_time_overlays_[pos] = time;
    invalidate_all_timed_caches();
    return previous;
  }

  double set_source_lower_bound_overlay(const semantic::Index source_id,
                                        const double time) {
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos >= source_lower_bound_overlays_.size()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double previous = source_lower_bound_overlays_[pos];
    if ((std::isfinite(previous) && std::isfinite(time) &&
         std::fabs(previous - time) <= 1e-12) ||
        (!std::isfinite(previous) && !std::isfinite(time))) {
      return previous;
    }
    source_lower_bound_overlays_[pos] = time;
    invalidate_all_timed_caches();
    return previous;
  }

  double set_source_upper_bound_overlay(const semantic::Index source_id,
                                        const double time) {
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos >= source_upper_bound_overlays_.size()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double previous = source_upper_bound_overlays_[pos];
    if ((std::isfinite(previous) && std::isfinite(time) &&
         std::fabs(previous - time) <= 1e-12) ||
        (!std::isfinite(previous) && !std::isfinite(time))) {
      return previous;
    }
    source_upper_bound_overlays_[pos] = time;
    invalidate_all_timed_caches();
    return previous;
  }

  ExprUpperBoundOverlayState set_expr_upper_bound_overlay(
      const semantic::Index expr_id,
      const double time,
      const double normalizer) {
    const auto pos = static_cast<std::size_t>(expr_id);
    if (pos >= expr_upper_bound_overlays_.size()) {
      return ExprUpperBoundOverlayState{};
    }
    ExprUpperBoundOverlayState previous{
        expr_upper_bound_overlays_[pos],
        expr_upper_bound_normalizers_[pos]};
    expr_upper_bound_overlays_[pos] = time;
    expr_upper_bound_normalizers_[pos] = normalizer;
    return previous;
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

  const double *exact_time_for(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (source_id == semantic::kInvalidIndex) {
      return nullptr;
    }
    if (condition_frame != nullptr) {
      const auto pos = static_cast<std::size_t>(source_id);
      if (pos < condition_frame->source_exact_times.size() &&
          std::isfinite(condition_frame->source_exact_times[pos])) {
        return &condition_frame->source_exact_times[pos];
      }
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos < exact_time_overlays_.size()) {
      const auto overlay = exact_time_overlays_[pos];
      if (std::isfinite(overlay)) {
        return &exact_time_overlays_[pos];
      }
    }
    if (pos >= sequence_state_->exact_times.size()) {
      return nullptr;
    }
    const auto value = sequence_state_->exact_times[pos];
    if (!std::isfinite(value)) {
      return nullptr;
    }
    return &sequence_state_->exact_times[pos];
  }

  bool has_source_lower_bound_overlay(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (source_id == semantic::kInvalidIndex) {
      return false;
    }
    if (condition_frame != nullptr) {
      const auto pos = static_cast<std::size_t>(source_id);
      if (pos < condition_frame->source_lower_bounds.size() &&
          std::isfinite(condition_frame->source_lower_bounds[pos])) {
        return true;
      }
    }
    const auto pos = static_cast<std::size_t>(source_id);
    return pos < source_lower_bound_overlays_.size() &&
           std::isfinite(source_lower_bound_overlays_[pos]);
  }

  bool has_source_upper_bound_overlay(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (source_id == semantic::kInvalidIndex) {
      return false;
    }
    if (condition_frame != nullptr) {
      const auto pos = static_cast<std::size_t>(source_id);
      if (pos < condition_frame->source_upper_bounds.size() &&
          std::isfinite(condition_frame->source_upper_bounds[pos])) {
        return true;
      }
    }
    const auto pos = static_cast<std::size_t>(source_id);
    return pos < source_upper_bound_overlays_.size() &&
           std::isfinite(source_upper_bound_overlays_[pos]);
  }

  double lower_bound_for(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    double lower_bound = sequence_state_->lower_bound;
    if (source_id != semantic::kInvalidIndex) {
      if (condition_frame != nullptr) {
        const auto pos = static_cast<std::size_t>(source_id);
        if (pos < condition_frame->source_lower_bounds.size() &&
            std::isfinite(condition_frame->source_lower_bounds[pos])) {
          lower_bound = std::max(
              lower_bound, condition_frame->source_lower_bounds[pos]);
        }
      }
      const auto pos = static_cast<std::size_t>(source_id);
      if (pos < source_lower_bound_overlays_.size() &&
          std::isfinite(source_lower_bound_overlays_[pos])) {
        lower_bound = std::max(lower_bound, source_lower_bound_overlays_[pos]);
      }
    }
    return lower_bound;
  }

  double upper_bound_for(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) const {
    if (source_id == semantic::kInvalidIndex) {
      return std::numeric_limits<double>::infinity();
    }
    double upper_bound = std::numeric_limits<double>::infinity();
    if (condition_frame != nullptr) {
      const auto pos = static_cast<std::size_t>(source_id);
      if (pos < condition_frame->source_upper_bounds.size() &&
          std::isfinite(condition_frame->source_upper_bounds[pos])) {
        upper_bound =
            std::min(upper_bound, condition_frame->source_upper_bounds[pos]);
      }
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos < source_upper_bound_overlays_.size() &&
        std::isfinite(source_upper_bound_overlays_[pos])) {
      upper_bound = std::min(upper_bound, source_upper_bound_overlays_[pos]);
    }
    return upper_bound;
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

  leaf::EventChannels compute_conditional_source(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) {
    const double lower_bound = lower_bound_for(source_id, condition_frame);
    const double upper_bound = upper_bound_for(source_id, condition_frame);
    if (std::isfinite(upper_bound) && !(upper_bound > lower_bound)) {
      return impossible_channels();
    }
    if (const double *exact_time = exact_time_for(source_id, condition_frame)) {
      if (lower_bound > 0.0 && !(*exact_time > lower_bound)) {
        return impossible_channels();
      }
      if (std::isfinite(upper_bound) && !(*exact_time < upper_bound)) {
        return impossible_channels();
      }
      if (!(conditional_time_ >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }

    const auto uncond =
        load_base_source(
            &base_current_cache_, conditional_time_, source_id, condition_frame);
    if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
      return uncond;
    }

    leaf::EventChannels lower = leaf::EventChannels::impossible();
    if (!has_source_lower_bound_overlay(source_id, condition_frame) &&
        std::fabs(lower_bound - sequence_state_->lower_bound) <= 1e-12) {
      lower = load_base_source(
          &base_lower_cache_, lower_bound, source_id, condition_frame);
    } else if (lower_bound > 0.0) {
      const auto frame = base_time_guard(lower_bound);
      lower = base_source(source_id, condition_frame);
    }
    if (!std::isfinite(upper_bound)) {
      return conditionalize(uncond, lower);
    }
    const auto frame = base_time_guard(upper_bound);
    const auto upper = base_source(source_id, condition_frame);
    if (!std::isfinite(upper.cdf - lower.cdf) ||
        !(upper.cdf - lower.cdf > 0.0)) {
      return impossible_channels();
    }
    if (conditional_time_ >= upper_bound) {
      return leaf::EventChannels::certain();
    }
    if (conditional_time_ <= lower_bound) {
      return impossible_channels();
    }
    return conditionalize_between(uncond, lower, upper);
  }

  leaf::EventChannels base_source(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) {
    if (auto *cache = active_base_cache(); cache != nullptr) {
      return load_base_source(cache, cache->time, source_id, condition_frame);
    }
    return load_base_source(
        &base_current_cache_, conditional_time_, source_id, condition_frame);
  }

  leaf::EventChannels compute_base_source(
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) {
    if (const double *exact_time = exact_time_for(source_id, condition_frame)) {
      if (!(active_base_time() >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }
    if (source_id < program_.layout.n_leaves) {
      return base_leaf_channels(source_id, condition_frame);
    }
    return base_pool_channels(
        source_id - program_.layout.n_leaves, condition_frame);
  }

  leaf::EventChannels base_leaf_channels(
      const semantic::Index index,
      const ExactConditionFrame *condition_frame = nullptr) {
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
    if (const double *exact_time =
            exact_time_for(onset_source_id, condition_frame)) {
      return shifted_channels(*exact_time);
    }

    const auto batch = quadrature::build_finite_batch(0.0, upper);
    double pdf = 0.0;
    double cdf = 0.0;
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double u = batch.nodes.nodes[i];
      const auto frame = base_time_guard(u);
      const auto source = base_source(onset_source_id, condition_frame);
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

  leaf::EventChannels base_pool_channels(
      const semantic::Index index,
      const ExactConditionFrame *condition_frame = nullptr) {
    const auto pos = static_cast<std::size_t>(index);
    const auto begin = program_.pool_member_offsets[pos];
    const auto end = program_.pool_member_offsets[pos + 1U];
    const auto k = program_.pool_k[pos];
    const auto n_members = static_cast<int>(end - begin);

    std::vector<leaf::EventChannels> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(base_source(
          program_.pool_member_source_ids[static_cast<std::size_t>(i)],
          condition_frame));
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

  leaf::EventChannels load_conditional_source(
      TimedCache *cache,
      const semantic::Index source_id,
      const ExactConditionFrame *condition_frame = nullptr) {
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->epoch[pos] == cache->current_epoch) {
      return cache->channels[pos];
    }
    cache->channels[pos] = compute_conditional_source(source_id, condition_frame);
    cache->epoch[pos] = cache->current_epoch;
    return cache->channels[pos];
  }

  leaf::EventChannels load_base_source(TimedCache *cache,
                                       const double time,
                                       const semantic::Index source_id,
                                       const ExactConditionFrame *condition_frame = nullptr) {
    if (cache == nullptr || source_id == semantic::kInvalidIndex) {
      return impossible_channels();
    }
    cache->set_context(time, condition_frame_id(condition_frame));
    const auto pos = static_cast<std::size_t>(source_id);
    if (cache->epoch[pos] == cache->current_epoch) {
      return cache->channels[pos];
    }
    const bool needs_frame = active_base_cache() != cache;
    if (needs_frame) {
      base_cache_stack_.push_back(cache);
    }
    cache->channels[pos] = compute_base_source(source_id, condition_frame);
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
                            const double t,
                            const ExactConditionFrame *condition_frame = nullptr) {
  const auto guard = oracle->conditional_time_guard(t);
  return clamp_probability(
      oracle->conditional_source(source_id, condition_frame).cdf);
}

inline double source_pdf_at(ExactSourceOracle *oracle,
                            const semantic::Index source_id,
                            const double t,
                            const ExactConditionFrame *condition_frame = nullptr) {
  const auto guard = oracle->conditional_time_guard(t);
  return safe_density(
      oracle->conditional_source(source_id, condition_frame).pdf);
}

inline double source_survival_at(ExactSourceOracle *oracle,
                                 const semantic::Index source_id,
                                 const double t,
                                 const ExactConditionFrame *condition_frame = nullptr) {
  const auto guard = oracle->conditional_time_guard(t);
  return clamp_probability(
      oracle->conditional_source(source_id, condition_frame).survival);
}

} // namespace detail
} // namespace accumulatr::eval
