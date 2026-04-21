#pragma once

#include <array>

#include "exact_planner.hpp"
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

class ExactSourceOracle {
public:
  ExactSourceOracle(const ExactVariantPlan &plan,
                    const ParamView &params,
                    const int first_param_row,
                    const ExactTriggerState &trigger_state,
                    const ExactSequenceState &sequence_state,
                    const double top_time)
      : plan_(plan),
        program_(plan.lowered.program),
        params_(params),
        first_param_row_(first_param_row),
        trigger_state_(trigger_state),
        sequence_state_(sequence_state),
        current_time_(top_time),
        leaf_inputs_(static_cast<std::size_t>(program_.layout.n_leaves)),
        source_cache_(static_cast<std::size_t>(plan.source_count)),
        source_epoch_(static_cast<std::size_t>(plan.source_count), 0U) {
    for (int i = 0; i < program_.layout.n_leaves; ++i) {
      const auto pos = static_cast<std::size_t>(i);
      const auto &desc = program_.leaf_descriptors[pos];
      const int row = first_param_row_ + i;
      auto &loaded = leaf_inputs_[pos];
      const int n_local = std::min<int>(desc.param_count, 8);
      for (int j = 0; j < n_local; ++j) {
        loaded.params[static_cast<std::size_t>(j)] = params_.p(row, j);
      }
      for (int j = n_local; j < 8; ++j) {
        loaded.params[static_cast<std::size_t>(j)] = 0.0;
      }
      loaded.q = leaf_q(i, row);
      loaded.t0 = params_.t0(row);
    }
  }

  leaf::EventChannels source_channels(const semantic::SourceKind kind,
                                      const semantic::Index index,
                                      const double t) {
    set_time(t);
    return source_channels(ExactSourceKey{kind, index});
  }

private:
  const ExactVariantPlan &plan_;
  const runtime::ExactProgram &program_;
  const ParamView &params_;
  int first_param_row_;
  const ExactTriggerState &trigger_state_;
  const ExactSequenceState &sequence_state_;
  double current_time_;
  std::vector<ExactLoadedLeafInput> leaf_inputs_;
  std::vector<leaf::EventChannels> source_cache_;
  std::vector<std::uint32_t> source_epoch_;
  std::uint32_t current_epoch_{1U};

  const double *exact_time_for(const ExactSourceKey key) const {
    const auto source_id = source_ordinal(plan_, key);
    if (source_id == semantic::kInvalidIndex) {
      return nullptr;
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos >= sequence_state_.exact_times.size()) {
      return nullptr;
    }
    const auto value = sequence_state_.exact_times[pos];
    if (!std::isfinite(value)) {
      return nullptr;
    }
    return &sequence_state_.exact_times[pos];
  }

  leaf::EventChannels conditionalize(const leaf::EventChannels uncond,
                                     const leaf::EventChannels lower) const {
    if (!(sequence_state_.lower_bound > 0.0)) {
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
        trigger_state_.shared_started[static_cast<std::size_t>(trigger_index)] <= 1U) {
      return trigger_state_.shared_started[static_cast<std::size_t>(trigger_index)] == 1U
                 ? 0.0
                 : 1.0;
    }
    return params_.q(row);
  }

  void set_time(const double t) {
    if (std::fabs(t - current_time_) <= 1e-12) {
      return;
    }
    current_time_ = t;
    ++current_epoch_;
    if (current_epoch_ == 0U) {
      current_epoch_ = 1U;
      std::fill(source_epoch_.begin(), source_epoch_.end(), 0U);
    }
  }

  leaf::EventChannels source_channels(const ExactSourceKey key) {
    const auto source_id = source_ordinal(plan_, key);
    if (source_id == semantic::kInvalidIndex) {
      return compute_source_channels(key.kind, key.index, current_time_);
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (source_epoch_[pos] == current_epoch_) {
      return source_cache_[pos];
    }
    source_cache_[pos] = compute_source_channels(key.kind, key.index, current_time_);
    source_epoch_[pos] = current_epoch_;
    return source_cache_[pos];
  }

  leaf::EventChannels compute_source_channels(const semantic::SourceKind kind,
                                              const semantic::Index index,
                                              const double t) {
    const ExactSourceKey key{kind, index};
    if (const double *exact_time = exact_time_for(key)) {
      if (!(t >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }
    const auto uncond = base_source_channels(kind, index, t);
    if (!(sequence_state_.lower_bound > 0.0)) {
      return uncond;
    }
    const auto lower =
        base_source_channels(kind, index, sequence_state_.lower_bound);
    return conditionalize(uncond, lower);
  }

  leaf::EventChannels base_source_channels(const semantic::SourceKind kind,
                                           const semantic::Index index,
                                           const double t) {
    const ExactSourceKey key{kind, index};
    if (const double *exact_time = exact_time_for(key)) {
      if (!(t >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }
    if (kind == semantic::SourceKind::Leaf) {
      return base_leaf_channels(index, t);
    }
    return base_pool_channels(index, t);
  }

  leaf::EventChannels base_leaf_channels(const semantic::Index index,
                                         const double t) {
    const auto pos = static_cast<std::size_t>(index);
    const auto &desc = program_.leaf_descriptors[pos];
    const auto &loaded = leaf_inputs_[pos];
    const auto onset_kind = static_cast<semantic::OnsetKind>(desc.onset_kind);
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
    const auto source_kind = static_cast<semantic::SourceKind>(desc.onset_source_kind);
    const auto source_index = desc.onset_source_index;
    const ExactSourceKey onset_key{source_kind, source_index};
    auto shifted_channels = [&](const double source_time) {
      return standard_leaf_channels(
          desc.dist_kind,
          loaded.params.data(),
          std::min(desc.param_count, 8),
          loaded.q,
          loaded.t0,
          t - source_time - lag);
    };
    if (const double *exact_time = exact_time_for(onset_key)) {
      return shifted_channels(*exact_time);
    }
    const auto batch = quadrature::build_finite_batch(0.0, upper);
    double pdf = 0.0;
    double cdf = 0.0;
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double u = batch.nodes.nodes[i];
      const auto source = base_source_channels(source_kind, source_index, u);
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

  leaf::EventChannels base_pool_channels(const semantic::Index index,
                                         const double t) {
    const auto pos = static_cast<std::size_t>(index);
    const auto begin = program_.pool_member_offsets[pos];
    const auto end = program_.pool_member_offsets[pos + 1U];
    const auto k = program_.pool_k[pos];
    const auto n_members = static_cast<int>(end - begin);

    std::vector<leaf::EventChannels> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(base_source_channels(
          static_cast<semantic::SourceKind>(
              program_.pool_member_kind[static_cast<std::size_t>(i)]),
          program_.pool_member_indices[static_cast<std::size_t>(i)],
          t));
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
};

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
                            const ExactSourceKey key,
                            const double t) {
  return clamp_probability(oracle->source_channels(key.kind, key.index, t).cdf);
}

inline double source_pdf_at(ExactSourceOracle *oracle,
                            const ExactSourceKey key,
                            const double t) {
  return safe_density(oracle->source_channels(key.kind, key.index, t).pdf);
}

inline double source_survival_at(ExactSourceOracle *oracle,
                                 const ExactSourceKey key,
                                 const double t) {
  return clamp_probability(oracle->source_channels(key.kind, key.index, t).survival);
}

} // namespace detail
} // namespace accumulatr::eval
