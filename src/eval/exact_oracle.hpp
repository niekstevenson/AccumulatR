#pragma once

#include "exact_planner.hpp"

namespace accumulatr::eval {
namespace detail {

template <typename Fn>
double simpson_integrate(Fn &&fn, double lower, double upper, int n = 256) {
  if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
    return 0.0;
  }
  if (n < 2) {
    n = 2;
  }
  if (n % 2 != 0) {
    ++n;
  }
  const double h = (upper - lower) / static_cast<double>(n);
  double sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double x = lower + h * static_cast<double>(i);
    const double fx = fn(x);
    const double w = (i == 0 || i == n) ? 1.0 : ((i % 2 == 0) ? 2.0 : 4.0);
    sum += w * (std::isfinite(fx) ? fx : 0.0);
  }
  return sum * h / 3.0;
}

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
      : program_(plan.lowered.program),
        params_(params),
        first_param_row_(first_param_row),
        trigger_state_(trigger_state),
        sequence_state_(sequence_state),
        top_time_(top_time),
        leaf_cache_(static_cast<std::size_t>(program_.layout.n_leaves)),
        leaf_ready_(static_cast<std::size_t>(program_.layout.n_leaves), 0U),
        pool_cache_(static_cast<std::size_t>(program_.layout.n_pools)),
        pool_ready_(static_cast<std::size_t>(program_.layout.n_pools), 0U) {}

  leaf::EventChannels source_channels(const semantic::SourceKind kind,
                                      const semantic::Index index,
                                      const double t) {
    if (std::fabs(t - top_time_) <= 1e-12) {
      if (kind == semantic::SourceKind::Leaf) {
        return cached_leaf(index, t);
      }
      if (kind == semantic::SourceKind::Pool) {
        return cached_pool(index, t);
      }
      throw std::runtime_error("exact kernel does not support special sources");
    }
    return compute_source_channels(kind, index, t);
  }

private:
  const runtime::ExactProgram &program_;
  const ParamView &params_;
  int first_param_row_;
  const ExactTriggerState &trigger_state_;
  const ExactSequenceState &sequence_state_;
  double top_time_;
  std::vector<leaf::EventChannels> leaf_cache_;
  std::vector<std::uint8_t> leaf_ready_;
  std::vector<leaf::EventChannels> pool_cache_;
  std::vector<std::uint8_t> pool_ready_;

  const double *exact_time_for(const ExactSourceKey key) const {
    const auto it = sequence_state_.exact_times.find(key);
    if (it == sequence_state_.exact_times.end()) {
      return nullptr;
    }
    return &it->second;
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

  leaf::EventChannels cached_leaf(const semantic::Index index, const double t) {
    const auto pos = static_cast<std::size_t>(index);
    if (leaf_ready_[pos] != 0U) {
      return leaf_cache_[pos];
    }
    leaf_ready_[pos] = 1U;
    leaf_cache_[pos] = compute_source_channels(semantic::SourceKind::Leaf, index, t);
    return leaf_cache_[pos];
  }

  leaf::EventChannels cached_pool(const semantic::Index index, const double t) {
    const auto pos = static_cast<std::size_t>(index);
    if (pool_ready_[pos] != 0U) {
      return pool_cache_[pos];
    }
    pool_ready_[pos] = 1U;
    pool_cache_[pos] = compute_source_channels(semantic::SourceKind::Pool, index, t);
    return pool_cache_[pos];
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
    if (kind == semantic::SourceKind::Pool) {
      return base_pool_channels(index, t);
    }
    throw std::runtime_error("exact kernel does not support special sources");
  }

  leaf::EventChannels base_leaf_channels(const semantic::Index index,
                                         const double t) {
    const auto pos = static_cast<std::size_t>(index);
    const int row = first_param_row_ + index;
    const auto begin = program_.parameter_layout.leaf_param_offsets[pos];
    const auto end = program_.parameter_layout.leaf_param_offsets[pos + 1U];
    double local_params[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const int n_local = std::min<int>(end - begin, 8);
    for (int j = 0; j < n_local; ++j) {
      local_params[j] = params_.p(row, j);
    }

    const double q = leaf_q(index, row);
    const double t0 = params_.t0(row);
    const auto onset_kind =
        static_cast<semantic::OnsetKind>(program_.onset_kind[pos]);
    if (onset_kind == semantic::OnsetKind::Absolute) {
      return standard_leaf_channels(
          program_.leaf_dist_kind[pos],
          local_params,
          end - begin,
          q,
          t0,
          t - program_.onset_abs_value[pos]);
    }

    const double lag = program_.onset_lag[pos];
    const double upper = t - lag;
    if (!(upper > 0.0)) {
      return impossible_channels();
    }
    const auto source_kind =
        static_cast<semantic::SourceKind>(program_.onset_source_kind[pos]);
    const auto source_index = program_.onset_source_index[pos];
    const ExactSourceKey onset_key{source_kind, source_index};
    auto shifted_channels = [&](const double source_time) {
      return standard_leaf_channels(
          program_.leaf_dist_kind[pos],
          local_params,
          end - begin,
          q,
          t0,
          t - source_time - lag);
    };
    if (const double *exact_time = exact_time_for(onset_key)) {
      return shifted_channels(*exact_time);
    }
    const auto pdf = simpson_integrate(
        [&](const double u) {
          const auto source = base_source_channels(source_kind, source_index, u);
          return source.pdf * shifted_channels(u).pdf;
        },
        0.0,
        upper);
    const auto cdf = simpson_integrate(
        [&](const double u) {
          const auto source = base_source_channels(source_kind, source_index, u);
          return source.pdf * shifted_channels(u).cdf;
        },
        0.0,
        upper);

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
    if (n_members <= 0 || k < 1 || k > n_members) {
      return impossible_channels();
    }

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
