#include "trial_params.h"

#include <cstring>
#include <unordered_map>
#include <unordered_set>

#include "proto.h"

namespace {

inline void na_key_mix_u64(std::uint64_t &h1, std::uint64_t &h2,
                           std::uint64_t value) {
  h1 = mix_hash64(h1 ^ (value + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2)));
  h2 = mix_hash64(h2 ^ (value + 0xbf58476d1ce4e5b9ULL + (h2 << 6) + (h2 >> 2)));
}

inline void na_key_mix_i32(std::uint64_t &h1, std::uint64_t &h2, int value) {
  const std::uint64_t v = static_cast<std::uint64_t>(
      static_cast<std::uint32_t>(static_cast<std::int32_t>(value)));
  na_key_mix_u64(h1, h2, v);
}

inline void na_key_mix_u32(std::uint64_t &h1, std::uint64_t &h2,
                           std::uint32_t value) {
  na_key_mix_u64(h1, h2, static_cast<std::uint64_t>(value));
}

inline void na_key_mix_double(std::uint64_t &h1, std::uint64_t &h2,
                              double value) {
  na_key_mix_u64(h1, h2, canonical_double_bits(value));
}

inline void key_mix_bool(uuber::NAMapCacheKey &key, bool value) {
  na_key_mix_u64(key.hash1, key.hash2, value ? 1ULL : 0ULL);
}

inline void key_mix_i32(uuber::NAMapCacheKey &key, int value) {
  na_key_mix_i32(key.hash1, key.hash2, value);
}

inline void key_mix_u64(uuber::NAMapCacheKey &key, std::uint64_t value) {
  na_key_mix_u64(key.hash1, key.hash2, value);
}

inline void key_mix_double(uuber::NAMapCacheKey &key, double value) {
  na_key_mix_double(key.hash1, key.hash2, value);
}

inline void key_mix_string(uuber::NAMapCacheKey &key, const std::string &value) {
  key_mix_u64(key, static_cast<std::uint64_t>(value.size()));
  for (unsigned char ch : value) {
    key_mix_u64(key, static_cast<std::uint64_t>(ch));
  }
}

inline std::uint64_t acc_signature_except_q(const uuber::TrialParamsSoA &soa,
                                            std::size_t idx) {
  uuber::NAMapCacheKey key{0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL};
  key_mix_i32(key, idx < soa.dist_code.size() ? soa.dist_code[idx] : -1);
  key_mix_double(key, idx < soa.onset.size() ? soa.onset[idx]
                                             : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.t0.size() ? soa.t0[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p1.size() ? soa.p1[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p2.size() ? soa.p2[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p3.size() ? soa.p3[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p4.size() ? soa.p4[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p5.size() ? soa.p5[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p6.size() ? soa.p6[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p7.size() ? soa.p7[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  key_mix_double(key, idx < soa.p8.size() ? soa.p8[idx]
                                          : std::numeric_limits<double>::quiet_NaN());
  return key.hash1 ^ (key.hash2 + 0x9e3779b97f4a7c15ULL + (key.hash1 << 6) +
                      (key.hash1 >> 2));
}

inline void refresh_trial_params_soa_signatures(uuber::TrialParamsSoA &soa) {
  if (!soa.valid || soa.n_acc <= 0) {
    soa.acc_signature_except_q.clear();
    return;
  }
  const std::size_t n_acc = static_cast<std::size_t>(soa.n_acc);
  soa.acc_signature_except_q.resize(n_acc);
  for (std::size_t acc_idx = 0; acc_idx < n_acc; ++acc_idx) {
    soa.acc_signature_except_q[acc_idx] = acc_signature_except_q(soa, acc_idx);
  }
}

inline void build_base_trial_params_soa_from_context(
    const uuber::NativeContext &ctx, uuber::TrialParamsSoA &out) {
  out = uuber::TrialParamsSoA{};
  out.n_acc = static_cast<int>(ctx.accumulators.size());
  out.dist_code.resize(ctx.accumulators.size(), 0);
  out.onset.resize(ctx.accumulators.size(), 0.0);
  out.q.resize(ctx.accumulators.size(), 0.0);
  out.t0.resize(ctx.accumulators.size(), 0.0);
  out.p1.resize(ctx.accumulators.size(), 0.0);
  out.p2.resize(ctx.accumulators.size(), 0.0);
  out.p3.resize(ctx.accumulators.size(), 0.0);
  out.p4.resize(ctx.accumulators.size(), 0.0);
  out.p5.resize(ctx.accumulators.size(), 0.0);
  out.p6.resize(ctx.accumulators.size(), 0.0);
  out.p7.resize(ctx.accumulators.size(), 0.0);
  out.p8.resize(ctx.accumulators.size(), 0.0);
  for (std::size_t i = 0; i < ctx.accumulators.size(); ++i) {
    const uuber::NativeAccumulator &acc = ctx.accumulators[i];
    out.dist_code[i] = acc.dist_cfg.code;
    out.onset[i] = acc.onset;
    out.q[i] = acc.q;
    out.t0[i] = acc.dist_cfg.t0;
    out.p1[i] = acc.dist_cfg.p1;
    out.p2[i] = acc.dist_cfg.p2;
    out.p3[i] = acc.dist_cfg.p3;
    out.p4[i] = acc.dist_cfg.p4;
    out.p5[i] = acc.dist_cfg.p5;
    out.p6[i] = acc.dist_cfg.p6;
    out.p7[i] = acc.dist_cfg.p7;
    out.p8[i] = acc.dist_cfg.p8;
  }
  out.valid = true;
  refresh_trial_params_soa_signatures(out);
}

struct TrialSharedTriggerGroup {
  std::vector<int> acc_indices;
  double q{0.0};
};

inline bool resolve_trial_shared_trigger_q_value(
    const uuber::NativeContext &ctx, const TrialParamSet &params, int acc_idx,
    double &q_out) {
  if (acc_idx >= 0 && acc_idx < static_cast<int>(params.acc_params.size())) {
    const double shared_q =
        params.acc_params[static_cast<std::size_t>(acc_idx)].shared_q;
    if (std::isfinite(shared_q)) {
      q_out = clamp_probability(shared_q);
      return true;
    }
  }
  if (acc_idx >= 0 && acc_idx < static_cast<int>(ctx.accumulators.size())) {
    q_out = clamp_probability(
        ctx.accumulators[static_cast<std::size_t>(acc_idx)].q);
    return true;
  }
  return false;
}

inline std::vector<TrialSharedTriggerGroup> derive_shared_trigger_groups(
    const uuber::NativeContext &ctx, const TrialParamSet &params) {
  std::vector<TrialSharedTriggerGroup> groups;
  std::unordered_map<std::string, int> group_index;
  const int n_acc = std::min(static_cast<int>(params.acc_params.size()),
                             static_cast<int>(ctx.accumulators.size()));
  for (int acc_idx = 0; acc_idx < n_acc; ++acc_idx) {
    const TrialAccumulatorParams &acc =
        params.acc_params[static_cast<std::size_t>(acc_idx)];
    if (acc.shared_trigger_id.empty()) {
      continue;
    }
    auto it = group_index.find(acc.shared_trigger_id);
    if (it == group_index.end()) {
      double q_val = 0.0;
      if (!resolve_trial_shared_trigger_q_value(ctx, params, acc_idx, q_val)) {
        continue;
      }
      TrialSharedTriggerGroup group;
      group.q = q_val;
      groups.push_back(std::move(group));
      const int group_idx = static_cast<int>(groups.size()) - 1;
      group_index.emplace(acc.shared_trigger_id, group_idx);
      it = group_index.find(acc.shared_trigger_id);
    } else {
      double q_val = 0.0;
      if (resolve_trial_shared_trigger_q_value(ctx, params, acc_idx, q_val) &&
          std::fabs(q_val - groups[static_cast<std::size_t>(it->second)].q) >
              1e-12) {
        Rcpp::stop(
            "Inconsistent shared trigger q values for trigger '%s'",
            acc.shared_trigger_id.c_str());
      }
    }
    groups[static_cast<std::size_t>(it->second)].acc_indices.push_back(acc_idx);
  }
  return groups;
}

inline void initialize_trigger_q_states_inplace(TrialParamSet &params,
                                                const uuber::NativeContext &ctx,
                                                const SharedTriggerPlan &plan,
                                                bool fail) {
  for (int i = 0; i < static_cast<int>(plan.trigger_q.size()); ++i) {
    apply_trigger_q_soa_inplace_plan(ctx, params, plan, i, fail);
  }
}

inline void compute_trial_param_fingerprints(
    const TrialParamSet &params, uuber::NAMapCacheKey &source_key_out,
    uuber::NAMapCacheKey &value_key_out) {
  uuber::NAMapCacheKey source_key{0x243f6a8885a308d3ULL,
                                  0x13198a2e03707344ULL};
  uuber::NAMapCacheKey value_key{0xa4093822299f31d0ULL,
                                 0x082efa98ec4e6c89ULL};
  key_mix_bool(source_key, params.shared_trigger_layout_matches_context);
  key_mix_bool(value_key, params.shared_trigger_layout_matches_context);
  const std::uint64_t acc_count =
      static_cast<std::uint64_t>(params.acc_params.size());
  key_mix_u64(source_key, acc_count);
  key_mix_u64(value_key, acc_count);
  std::uint64_t shared_gate_count = 0ULL;
  for (const TrialAccumulatorParams &acc : params.acc_params) {
    const bool has_shared = !acc.shared_trigger_id.empty();
    key_mix_bool(source_key, has_shared);
    if (has_shared) {
      ++shared_gate_count;
      key_mix_double(source_key, acc.shared_q);
      key_mix_string(source_key, acc.shared_trigger_id);
    }
    key_mix_double(value_key, acc.onset);
    key_mix_i32(value_key, acc.onset_kind);
    key_mix_i32(value_key, acc.onset_source_acc_idx);
    key_mix_i32(value_key, acc.onset_source_pool_idx);
    key_mix_double(value_key, acc.onset_lag);
    key_mix_double(value_key, acc.q);
    key_mix_double(value_key, acc.shared_q);
    key_mix_i32(value_key, acc.dist_cfg.code);
    key_mix_double(value_key, acc.dist_cfg.t0);
    key_mix_double(value_key, acc.dist_cfg.p1);
    key_mix_double(value_key, acc.dist_cfg.p2);
    key_mix_double(value_key, acc.dist_cfg.p3);
    key_mix_double(value_key, acc.dist_cfg.p4);
    key_mix_double(value_key, acc.dist_cfg.p5);
    key_mix_double(value_key, acc.dist_cfg.p6);
    key_mix_double(value_key, acc.dist_cfg.p7);
    key_mix_double(value_key, acc.dist_cfg.p8);
    key_mix_string(value_key, acc.shared_trigger_id);
    key_mix_u64(value_key,
                static_cast<std::uint64_t>(acc.component_indices.size()));
    for (int comp_idx : acc.component_indices) {
      key_mix_i32(value_key, comp_idx);
    }
  }
  key_mix_u64(source_key, shared_gate_count);
  source_key_out = source_key;
  value_key_out = value_key;
}

} // namespace

std::uint64_t compute_trial_param_fingerprint(const TrialParamSet &params) {
  uuber::NAMapCacheKey source_key;
  uuber::NAMapCacheKey value_key;
  compute_trial_param_fingerprints(params, source_key, value_key);
  return source_key.hash1;
}

void refresh_trial_param_fingerprint(TrialParamSet &params) {
  compute_trial_param_fingerprints(params, params.shared_trigger_source_key,
                                   params.value_key);
  params.shared_trigger_source_fingerprint =
      params.shared_trigger_source_key.hash1;
  params.value_fingerprint = params.value_key.hash1;
  params.value_fingerprint_valid = true;
}

bool trial_paramsets_equivalent(const TrialParamSet &a,
                                const TrialParamSet &b) {
  if (a.shared_trigger_layout_matches_context !=
          b.shared_trigger_layout_matches_context ||
      a.acc_params.size() != b.acc_params.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.acc_params.size(); ++i) {
    const TrialAccumulatorParams &x = a.acc_params[i];
    const TrialAccumulatorParams &y = b.acc_params[i];
    const bool x_shared = !x.shared_trigger_id.empty();
    const bool y_shared = !y.shared_trigger_id.empty();
    if (x_shared != y_shared) {
      return false;
    }
    if (!x_shared) {
      continue;
    }
    if (x.shared_q != y.shared_q) {
      return false;
    }
  }
  return true;
}

std::uint64_t compute_trial_param_value_fingerprint(const TrialParamSet &params) {
  if (params.value_fingerprint_valid) {
    return params.value_fingerprint;
  }
  uuber::NAMapCacheKey source_key;
  uuber::NAMapCacheKey value_key;
  compute_trial_param_fingerprints(params, source_key, value_key);
  return value_key.hash1;
}

bool trial_paramsets_value_equivalent(const TrialParamSet &a,
                                      const TrialParamSet &b) {
  if (a.shared_trigger_layout_matches_context !=
          b.shared_trigger_layout_matches_context ||
      a.acc_params.size() != b.acc_params.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.acc_params.size(); ++i) {
    const TrialAccumulatorParams &x = a.acc_params[i];
    const TrialAccumulatorParams &y = b.acc_params[i];
    if (canonical_double_bits(x.onset) != canonical_double_bits(y.onset) ||
        x.onset_kind != y.onset_kind ||
        x.onset_source_acc_idx != y.onset_source_acc_idx ||
        x.onset_source_pool_idx != y.onset_source_pool_idx ||
        canonical_double_bits(x.onset_lag) != canonical_double_bits(y.onset_lag) ||
        canonical_double_bits(x.q) != canonical_double_bits(y.q) ||
        canonical_double_bits(x.shared_q) != canonical_double_bits(y.shared_q) ||
        x.dist_cfg.code != y.dist_cfg.code ||
        canonical_double_bits(x.dist_cfg.t0) != canonical_double_bits(y.dist_cfg.t0) ||
        canonical_double_bits(x.dist_cfg.p1) != canonical_double_bits(y.dist_cfg.p1) ||
        canonical_double_bits(x.dist_cfg.p2) != canonical_double_bits(y.dist_cfg.p2) ||
        canonical_double_bits(x.dist_cfg.p3) != canonical_double_bits(y.dist_cfg.p3) ||
        canonical_double_bits(x.dist_cfg.p4) != canonical_double_bits(y.dist_cfg.p4) ||
        canonical_double_bits(x.dist_cfg.p5) != canonical_double_bits(y.dist_cfg.p5) ||
        canonical_double_bits(x.dist_cfg.p6) != canonical_double_bits(y.dist_cfg.p6) ||
        canonical_double_bits(x.dist_cfg.p7) != canonical_double_bits(y.dist_cfg.p7) ||
        canonical_double_bits(x.dist_cfg.p8) != canonical_double_bits(y.dist_cfg.p8) ||
        x.shared_trigger_id != y.shared_trigger_id ||
        x.component_indices != y.component_indices) {
      return false;
    }
  }
  return true;
}

uuber::NAMapCacheKey
build_na_map_cache_key_idx(const uuber::OutcomeContextInfo &info,
                           const TrialParamSet *params_ptr,
                           const std::vector<int> &component_indices,
                           const std::vector<double> &comp_weights,
                           int outcome_idx_context) {
  uuber::NAMapCacheKey key;
  if (!params_ptr) {
    return key;
  }
  std::uint64_t h1 = 0x243f6a8885a308d3ULL;
  std::uint64_t h2 = 0x13198a2e03707344ULL;
  na_key_mix_i32(h1, h2, info.node_id);
  na_key_mix_i32(h1, h2, outcome_idx_context);
  na_key_mix_u32(h1, h2,
                 static_cast<std::uint32_t>(params_ptr->acc_params.size()));
  const bool has_soa_q =
      params_ptr->shared_trigger_mask_valid &&
      params_ptr->soa_cache_valid && params_ptr->soa_cache.valid &&
      params_ptr->soa_cache.q.size() >= params_ptr->acc_params.size();
  for (std::size_t acc_idx = 0; acc_idx < params_ptr->acc_params.size();
       ++acc_idx) {
    const TrialAccumulatorParams &acc = params_ptr->acc_params[acc_idx];
    const double effective_q =
        has_soa_q ? params_ptr->soa_cache.q[acc_idx] : acc.q;
    na_key_mix_double(h1, h2, effective_q);
    na_key_mix_double(h1, h2, acc.shared_q);
    na_key_mix_double(h1, h2, acc.onset);
    na_key_mix_i32(h1, h2, acc.onset_kind);
    na_key_mix_i32(h1, h2, acc.onset_source_acc_idx);
    na_key_mix_i32(h1, h2, acc.onset_source_pool_idx);
    na_key_mix_double(h1, h2, acc.onset_lag);
    na_key_mix_double(h1, h2, acc.dist_cfg.t0);
    na_key_mix_double(h1, h2, acc.dist_cfg.p1);
    na_key_mix_double(h1, h2, acc.dist_cfg.p2);
    na_key_mix_double(h1, h2, acc.dist_cfg.p3);
    na_key_mix_double(h1, h2, acc.dist_cfg.p4);
    na_key_mix_double(h1, h2, acc.dist_cfg.p5);
    na_key_mix_double(h1, h2, acc.dist_cfg.p6);
    na_key_mix_double(h1, h2, acc.dist_cfg.p7);
    na_key_mix_double(h1, h2, acc.dist_cfg.p8);
  }
  na_key_mix_u32(h1, h2, static_cast<std::uint32_t>(component_indices.size()));
  for (std::size_t i = 0; i < component_indices.size(); ++i) {
    na_key_mix_i32(h1, h2, component_indices[i]);
    const double w = (i < comp_weights.size()) ? comp_weights[i] : 0.0;
    na_key_mix_double(h1, h2, w);
  }
  key.hash1 = h1;
  key.hash2 = h2;
  return key;
}

TrialAccumulatorParams base_params(const uuber::NativeAccumulator &base) {
  TrialAccumulatorParams p;
  p.onset = base.onset;
  p.onset_kind = base.onset_kind;
  p.onset_source_acc_idx = base.onset_source_acc_idx;
  p.onset_source_pool_idx = base.onset_source_pool_idx;
  p.onset_lag = base.onset_lag;
  p.q = base.q;
  p.dist_cfg = base.dist_cfg;
  p.shared_q = base.q;
  p.shared_trigger_id = base.shared_trigger_id;
  p.has_components = !base.components.empty();
  p.components = base.components;
  p.component_indices = base.component_indices;
  p.has_override = false;
  return p;
}

TrialParamSet build_base_paramset(const uuber::NativeContext &ctx) {
  TrialParamSet ps;
  ps.acc_params.reserve(ctx.accumulators.size());
  for (const auto &acc : ctx.accumulators) {
    ps.acc_params.push_back(base_params(acc));
  }
  refresh_trial_param_fingerprint(ps);
  return ps;
}

void materialize_trial_params_soa(const uuber::NativeContext &ctx,
                                  const TrialParamSet *trial_params,
                                  uuber::TrialParamsSoA &out) {
  if (ctx.base_params_soa.valid &&
      ctx.base_params_soa.n_acc == static_cast<int>(ctx.accumulators.size())) {
    out = ctx.base_params_soa;
  } else {
    build_base_trial_params_soa_from_context(ctx, out);
  }
  if (!trial_params || !trial_params->has_any_override) {
    return;
  }
  const int n_acc =
      std::min(static_cast<int>(trial_params->acc_params.size()), out.n_acc);
  for (int acc_idx = 0; acc_idx < n_acc; ++acc_idx) {
    const TrialAccumulatorParams &entry =
        trial_params->acc_params[static_cast<std::size_t>(acc_idx)];
    if (!entry.has_override) {
      continue;
    }
    out.dist_code[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.code;
    out.onset[static_cast<std::size_t>(acc_idx)] = entry.onset;
    const double effective_q =
        (!entry.shared_trigger_id.empty() && std::isfinite(entry.shared_q))
            ? entry.shared_q
            : entry.q;
    out.q[static_cast<std::size_t>(acc_idx)] =
        clamp_probability(effective_q);
    out.t0[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.t0;
    out.p1[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p1;
    out.p2[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p2;
    out.p3[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p3;
    out.p4[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p4;
    out.p5[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p5;
    out.p6[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p6;
    out.p7[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p7;
    out.p8[static_cast<std::size_t>(acc_idx)] = entry.dist_cfg.p8;
  }
  out.valid = true;
  refresh_trial_params_soa_signatures(out);
}

std::size_t shared_trigger_count(const SharedTriggerPlan &plan) {
  return plan.trigger_q.size();
}

double shared_trigger_mask_weight(const SharedTriggerPlan &plan,
                                  std::uint64_t mask) {
  double weight = 1.0;
  for (std::size_t i = 0; i < plan.trigger_q.size(); ++i) {
    const double q = clamp_probability(plan.trigger_q[i]);
    const bool fail = ((mask >> i) & 1ULL) != 0ULL;
    weight *= fail ? q : (1.0 - q);
    if (!std::isfinite(weight) || weight <= 0.0) {
      return 0.0;
    }
  }
  return weight;
}

SharedTriggerPlan build_shared_trigger_plan(const uuber::NativeContext &ctx,
                                            const TrialParamSet *params_ptr) {
  SharedTriggerPlan plan;
  if (!params_ptr) {
    return plan;
  }
  const std::vector<TrialSharedTriggerGroup> groups =
      derive_shared_trigger_groups(ctx, *params_ptr);
  for (const TrialSharedTriggerGroup &group : groups) {
    if (group.acc_indices.empty()) {
      continue;
    }
    plan.trigger_acc_begin.push_back(
        static_cast<int>(plan.trigger_acc_indices.size()));
    plan.trigger_acc_count.push_back(
        static_cast<int>(group.acc_indices.size()));
    plan.trigger_acc_indices.insert(plan.trigger_acc_indices.end(),
                                    group.acc_indices.begin(),
                                    group.acc_indices.end());
    plan.trigger_q.push_back(group.q);
  }

  const std::size_t trigger_count = shared_trigger_count(plan);
  if (trigger_count == 0) {
    return plan;
  }
  if (trigger_count >= 63) {
    Rcpp::stop("Too many shared triggers for density evaluation");
  }
  const std::uint64_t n_states = 1ULL << trigger_count;
  constexpr std::uint64_t kMaxPrecomputedWeights = 1ULL << 20;
  if (n_states <= kMaxPrecomputedWeights) {
    plan.mask_weights.resize(n_states);
    for (std::uint64_t mask = 0; mask < n_states; ++mask) {
      plan.mask_weights[static_cast<std::size_t>(mask)] =
          shared_trigger_mask_weight(plan, mask);
    }
  }
  return plan;
}

void apply_trigger_state_inplace(TrialParamSet &params, int acc_idx, bool fail) {
  if (acc_idx < 0 || acc_idx >= static_cast<int>(params.acc_params.size())) {
    return;
  }
  const double q_val = fail ? 1.0 : 0.0;
  TrialAccumulatorParams &param =
      params.acc_params[static_cast<std::size_t>(acc_idx)];
  param.q = q_val;
  param.has_override = true;
  params.has_any_override = true;
  params.value_fingerprint_valid = false;
  if (params.soa_cache_valid && params.soa_cache.valid &&
      acc_idx < params.soa_cache.n_acc &&
      static_cast<std::size_t>(acc_idx) < params.soa_cache.q.size()) {
    params.soa_cache.q[static_cast<std::size_t>(acc_idx)] = q_val;
  } else {
    params.soa_cache_valid = false;
  }
}

void ensure_trial_params_soa(const uuber::NativeContext &ctx,
                             TrialParamSet &params) {
  if (params.soa_cache_valid && params.soa_cache.valid &&
      params.soa_cache.n_acc == static_cast<int>(ctx.accumulators.size()) &&
      params.soa_cache.q.size() == ctx.accumulators.size()) {
    return;
  }
  materialize_trial_params_soa(ctx, &params, params.soa_cache);
  params.soa_cache_valid = true;
}

const uuber::TrialParamsSoA *
resolve_trial_params_soa(const uuber::NativeContext &ctx,
                         const TrialParamSet *trial_params) {
  if (trial_params && trial_params->soa_cache_valid &&
      trial_params->soa_cache.valid &&
      trial_params->soa_cache.n_acc == static_cast<int>(ctx.accumulators.size()) &&
      trial_params->soa_cache.q.size() == ctx.accumulators.size()) {
    return &trial_params->soa_cache;
  }
  if (trial_params && trial_params->has_any_override) {
    materialize_trial_params_soa(ctx, trial_params, trial_params->soa_cache);
    trial_params->soa_cache_valid = true;
    if (trial_params->soa_cache.valid &&
        trial_params->soa_cache.n_acc ==
            static_cast<int>(ctx.accumulators.size()) &&
        trial_params->soa_cache.q.size() == ctx.accumulators.size()) {
      return &trial_params->soa_cache;
    }
  }
  if (ctx.base_params_soa.valid &&
      ctx.base_params_soa.n_acc == static_cast<int>(ctx.accumulators.size()) &&
      ctx.base_params_soa.q.size() == ctx.accumulators.size()) {
    return &ctx.base_params_soa;
  }
  return nullptr;
}

void apply_trigger_q_soa_inplace(TrialParamSet &params, int acc_idx, bool fail) {
  if (acc_idx < 0) {
    return;
  }
  if (!(params.soa_cache_valid && params.soa_cache.valid) ||
      acc_idx >= params.soa_cache.n_acc ||
      static_cast<std::size_t>(acc_idx) >= params.soa_cache.q.size()) {
    apply_trigger_state_inplace(params, acc_idx, fail);
    return;
  }
  params.soa_cache.q[static_cast<std::size_t>(acc_idx)] = fail ? 1.0 : 0.0;
}

void apply_trigger_q_soa_inplace_plan(const uuber::NativeContext &ctx,
                                      TrialParamSet &params,
                                      const SharedTriggerPlan &plan,
                                      int trigger_bit_idx, bool fail) {
  (void)ctx;
  if (trigger_bit_idx < 0 ||
      trigger_bit_idx >= static_cast<int>(plan.trigger_acc_begin.size()) ||
      trigger_bit_idx >= static_cast<int>(plan.trigger_acc_count.size())) {
    return;
  }
  const int begin =
      plan.trigger_acc_begin[static_cast<std::size_t>(trigger_bit_idx)];
  const int count =
      plan.trigger_acc_count[static_cast<std::size_t>(trigger_bit_idx)];
  if (begin < 0 || count <= 0) {
    return;
  }
  const int end = std::min(begin + count,
                           static_cast<int>(plan.trigger_acc_indices.size()));
  for (int i = begin; i < end; ++i) {
    apply_trigger_q_soa_inplace(
        params, plan.trigger_acc_indices[static_cast<std::size_t>(i)], fail);
  }
}

namespace {

inline void apply_trigger_q_soa_inplace_plan(
    const uuber::NativeContext &ctx, uuber::TrialParamsSoA &params_soa,
    const SharedTriggerPlan &plan, int trigger_bit_idx, bool fail) {
  (void)ctx;
  if (trigger_bit_idx < 0 ||
      trigger_bit_idx >= static_cast<int>(plan.trigger_acc_begin.size()) ||
      trigger_bit_idx >= static_cast<int>(plan.trigger_acc_count.size())) {
    return;
  }
  const int begin =
      plan.trigger_acc_begin[static_cast<std::size_t>(trigger_bit_idx)];
  const int count =
      plan.trigger_acc_count[static_cast<std::size_t>(trigger_bit_idx)];
  if (begin < 0 || count <= 0) {
    return;
  }
  const int end = std::min(begin + count,
                           static_cast<int>(plan.trigger_acc_indices.size()));
  const double q_value = fail ? 1.0 : 0.0;
  for (int i = begin; i < end; ++i) {
    const int acc_idx = plan.trigger_acc_indices[static_cast<std::size_t>(i)];
    if (acc_idx < 0 || acc_idx >= params_soa.n_acc ||
        static_cast<std::size_t>(acc_idx) >= params_soa.q.size()) {
      continue;
    }
    params_soa.q[static_cast<std::size_t>(acc_idx)] = q_value;
  }
}

} // namespace

void prepare_trigger_scratch(const uuber::NativeContext &ctx,
                             const TrialParamSet *base_params,
                             const SharedTriggerPlan &plan,
                             TrialParamSet &scratch) {
  if (!base_params) {
    return;
  }
  const std::uint64_t base_fingerprint =
      base_params->shared_trigger_source_fingerprint != 0ULL
          ? base_params->shared_trigger_source_fingerprint
          : compute_trial_param_fingerprint(*base_params);
  const bool needs_copy =
      !scratch.shared_trigger_scratch_ready ||
      scratch.acc_params.size() != base_params->acc_params.size() ||
      scratch.shared_trigger_source_fingerprint != base_fingerprint ||
      scratch.shared_trigger_source_params != base_params;
  if (needs_copy) {
    scratch = *base_params;
    scratch.shared_trigger_source_fingerprint = base_fingerprint;
    scratch.shared_trigger_source_params = base_params;
    scratch.shared_trigger_scratch_ready = true;
    scratch.shared_trigger_plan_identity = nullptr;
    scratch.shared_trigger_current_mask = 0ULL;
    scratch.shared_trigger_mask_valid = false;
  } else if (scratch.shared_trigger_source_params != base_params) {
    scratch.shared_trigger_source_params = base_params;
  }
  ensure_trial_params_soa(ctx, scratch);
  const void *plan_identity = static_cast<const void *>(&plan);
  if (!scratch.shared_trigger_mask_valid ||
      scratch.shared_trigger_plan_identity != plan_identity) {
    initialize_trigger_q_states_inplace(scratch, ctx, plan, false);
    scratch.shared_trigger_plan_identity = plan_identity;
    scratch.shared_trigger_current_mask = 0ULL;
    scratch.shared_trigger_mask_valid = true;
    return;
  }
  if (scratch.shared_trigger_current_mask == 0ULL) {
    return;
  }
  std::uint64_t diff = scratch.shared_trigger_current_mask;
  while (diff != 0ULL) {
    const int bit = static_cast<int>(__builtin_ctzll(diff));
    apply_trigger_q_soa_inplace_plan(ctx, scratch, plan, bit, false);
    diff &= (diff - 1ULL);
  }
  scratch.shared_trigger_current_mask = 0ULL;
  scratch.shared_trigger_mask_valid = true;
}

bool build_shared_trigger_mask_soa_batch(const uuber::NativeContext &ctx,
                                         const TrialParamSet *base_params,
                                         const SharedTriggerPlan &plan,
                                         SharedTriggerMaskSoABatch &out) {
  out.mask_params.clear();
  out.mask_param_ptrs.clear();
  out.mask_weights.clear();
  if (!base_params) {
    return false;
  }

  const uuber::TrialParamsSoA *base_soa = resolve_trial_params_soa(ctx, base_params);
  if (!base_soa || !base_soa->valid) {
    return false;
  }

  uuber::TrialParamsSoA success_soa = *base_soa;
  for (int acc_idx : plan.trigger_acc_indices) {
    if (acc_idx < 0 || acc_idx >= success_soa.n_acc ||
        static_cast<std::size_t>(acc_idx) >= success_soa.q.size()) {
      continue;
    }
    success_soa.q[static_cast<std::size_t>(acc_idx)] = 0.0;
  }

  const std::size_t trigger_count = shared_trigger_count(plan);
  if (trigger_count == 0u) {
    out.mask_params.push_back(success_soa);
    out.mask_weights.push_back(1.0);
  } else {
    if (trigger_count >= 63u) {
      Rcpp::stop("Too many shared triggers for density evaluation");
    }
    const std::uint64_t n_states = 1ULL << trigger_count;
    for (std::uint64_t mask = 0ULL; mask < n_states; ++mask) {
      const double weight = plan.mask_weights.empty()
                                ? shared_trigger_mask_weight(plan, mask)
                                : plan.mask_weights[static_cast<std::size_t>(mask)];
      if (!(std::isfinite(weight) && weight > 0.0)) {
        continue;
      }
      out.mask_params.push_back(success_soa);
      uuber::TrialParamsSoA &mask_soa = out.mask_params.back();
      for (std::size_t bit_idx = 0; bit_idx < trigger_count; ++bit_idx) {
        if (((mask >> bit_idx) & 1ULL) == 0ULL) {
          continue;
        }
        apply_trigger_q_soa_inplace_plan(
            ctx, mask_soa, plan, static_cast<int>(bit_idx), true);
      }
      mask_soa.valid = true;
      out.mask_weights.push_back(weight);
    }
  }

  out.mask_param_ptrs.reserve(out.mask_params.size());
  for (const uuber::TrialParamsSoA &mask_params : out.mask_params) {
    out.mask_param_ptrs.push_back(&mask_params);
  }
  return !out.mask_param_ptrs.empty() &&
         out.mask_param_ptrs.size() == out.mask_weights.size();
}

bool prepare_trial_params_runtime(const uuber::NativeContext &ctx,
                                  TrialParamSet *params_ptr,
                                  PreparedTrialParamsRuntime &out) {
  static_cast<void>(ctx);
  out = PreparedTrialParamsRuntime{};
  if (!params_ptr) {
    return false;
  }
  if (!params_ptr->value_fingerprint_valid) {
    refresh_trial_param_fingerprint(*params_ptr);
  }
  out.params = params_ptr;
  out.shared_trigger_source_key = params_ptr->shared_trigger_source_key;
  out.value_key = params_ptr->value_key;
  out.shared_trigger_source_fingerprint =
      params_ptr->shared_trigger_source_fingerprint;
  out.value_fingerprint = params_ptr->value_fingerprint;
  out.valid = true;
  return true;
}

bool ensure_prepared_trial_params_soa(const uuber::NativeContext &ctx,
                                      PreparedTrialParamsRuntime &runtime) {
  if (!runtime.valid || runtime.params == nullptr) {
    return false;
  }
  if (runtime.soa != nullptr && runtime.soa->valid) {
    return true;
  }
  runtime.soa = resolve_trial_params_soa(ctx, runtime.params);
  return runtime.soa != nullptr && runtime.soa->valid;
}

bool ensure_prepared_trial_params_trigger_plan(
    const uuber::NativeContext &ctx, PreparedTrialParamsRuntime &runtime) {
  if (!runtime.valid || runtime.params == nullptr) {
    return false;
  }
  if (runtime.trigger_plan_ready) {
    return true;
  }
  runtime.trigger_plan = build_shared_trigger_plan(ctx, runtime.params);
  runtime.trigger_plan_ready = true;
  return true;
}

bool ensure_prepared_trial_params_mask_batch(
    const uuber::NativeContext &ctx, PreparedTrialParamsRuntime &runtime) {
  if (!runtime.valid || runtime.params == nullptr) {
    return false;
  }
  if (!ensure_prepared_trial_params_trigger_plan(ctx, runtime)) {
    return false;
  }
  if (runtime.mask_batch_ready) {
    return !runtime.mask_batch.mask_param_ptrs.empty();
  }
  runtime.mask_batch_ready = true;
  if (shared_trigger_count(runtime.trigger_plan) == 0u) {
    runtime.mask_batch = SharedTriggerMaskSoABatch{};
    return false;
  }
  return build_shared_trigger_mask_soa_batch(ctx, runtime.params,
                                             runtime.trigger_plan,
                                             runtime.mask_batch);
}

namespace {

std::vector<uuber::ProtoParamEntry> params_from_rcpp(const Rcpp::List &params) {
  Rcpp::CharacterVector names = params.names();
  std::vector<uuber::ProtoParamEntry> out;
  out.reserve(params.size());
  for (R_xlen_t i = 0; i < params.size(); ++i) {
    if (params[i] == R_NilValue) {
      continue;
    }
    uuber::ProtoParamEntry entry;
    entry.name = Rcpp::as<std::string>(names[i]);
    SEXP val = params[i];
    if (Rf_isLogical(val)) {
      Rcpp::LogicalVector lv(val);
      if (lv.size() == 1) {
        entry.tag = uuber::ParamValueTag::LogicalScalar;
        entry.logical_scalar = lv[0];
      } else {
        entry.tag = uuber::ParamValueTag::LogicalVector;
        entry.logical_values.assign(lv.begin(), lv.end());
      }
    } else {
      Rcpp::NumericVector nv(val);
      if (nv.size() == 1) {
        entry.tag = uuber::ParamValueTag::NumericScalar;
        entry.numeric_scalar = nv[0];
      } else {
        entry.tag = uuber::ParamValueTag::NumericVector;
        entry.numeric_values.assign(nv.begin(), nv.end());
      }
    }
    out.push_back(std::move(entry));
  }
  return out;
}

std::vector<std::string> string_vector_from_entry(SEXP entry) {
  if (Rf_isNull(entry)) {
    return {};
  }
  Rcpp::CharacterVector vec(entry);
  std::vector<std::string> out;
  out.reserve(vec.size());
  for (R_xlen_t i = 0; i < vec.size(); ++i) {
    if (vec[i] == NA_STRING) {
      continue;
    }
    out.push_back(Rcpp::as<std::string>(vec[i]));
  }
  return out;
}

std::vector<int> component_indices_from_labels(
    const uuber::NativeContext &ctx, const std::vector<std::string> &labels) {
  std::vector<int> out;
  out.reserve(labels.size());
  for (const auto &label : labels) {
    auto it = ctx.component_index.find(label);
    if (it == ctx.component_index.end()) {
      Rcpp::stop("Unknown component label '%s' in trial override",
                 label.c_str());
    }
    out.push_back(it->second);
  }
  return out;
}

} // namespace

std::unique_ptr<TrialParamSet>
build_trial_params_from_df(const uuber::NativeContext &ctx,
                           const Rcpp::Nullable<Rcpp::DataFrame> &rows_opt) {
  if (rows_opt.isNull()) {
    return nullptr;
  }
  Rcpp::DataFrame rows(rows_opt.get());
  const R_xlen_t n = rows.nrows();
  if (n == 0) {
    return nullptr;
  }

  const bool has_acc = rows.containsElementNamed("accumulator");
  Rcpp::RObject acc_col_obj;
  if (has_acc) {
    acc_col_obj = rows["accumulator"];
  }

  Rcpp::CharacterVector dist_col = rows.containsElementNamed("dist")
                                       ? Rcpp::CharacterVector(rows["dist"])
                                       : Rcpp::CharacterVector();
  Rcpp::NumericVector onset_col = rows.containsElementNamed("onset")
                                      ? Rcpp::NumericVector(rows["onset"])
                                      : Rcpp::NumericVector();
  Rcpp::NumericVector q_col = rows.containsElementNamed("q")
                                  ? Rcpp::NumericVector(rows["q"])
                                  : Rcpp::NumericVector();
  Rcpp::NumericVector shared_q_col =
      rows.containsElementNamed("shared_trigger_q")
          ? Rcpp::NumericVector(rows["shared_trigger_q"])
          : Rcpp::NumericVector();
  Rcpp::CharacterVector shared_col =
      rows.containsElementNamed("shared_trigger_id")
          ? Rcpp::CharacterVector(rows["shared_trigger_id"])
          : Rcpp::CharacterVector();
  Rcpp::List comps_col = rows.containsElementNamed("components")
                             ? Rcpp::List(rows["components"])
                             : Rcpp::List();
  Rcpp::List params_list_col = rows.containsElementNamed("params")
                                   ? Rcpp::List(rows["params"])
                                   : Rcpp::List();

  const std::unordered_set<std::string> base_cols = {
      "trial",   "component", "accumulator",      "type",
      "outcome", "rt",        "params",           "onset",
      "q",       "condition", "component_weight", "shared_trigger_id",
      "dist",    "components"};
  std::vector<std::string> param_cols;
  std::vector<SEXP> param_col_data;
  std::vector<int> param_col_types;
  Rcpp::CharacterVector df_names = rows.names();
  if (!df_names.isNULL()) {
    for (R_xlen_t i = 0; i < df_names.size(); ++i) {
      if (df_names[i] == NA_STRING) {
        continue;
      }
      const std::string col_name = Rcpp::as<std::string>(df_names[i]);
      if (base_cols.find(col_name) != base_cols.end()) {
        continue;
      }
      SEXP column = rows[col_name];
      param_cols.push_back(col_name);
      param_col_data.push_back(column);
      param_col_types.push_back(TYPEOF(column));
    }
  }

  auto params_set = std::make_unique<TrialParamSet>();
  params_set->acc_params.assign(ctx.accumulators.size(),
                                TrialAccumulatorParams{});

  auto upsert_param_entry = [](std::vector<uuber::ProtoParamEntry> &entries,
                               uuber::ProtoParamEntry entry) {
    for (auto &existing : entries) {
      if (existing.name == entry.name) {
        existing = std::move(entry);
        return;
      }
    }
    entries.push_back(std::move(entry));
  };

  auto append_scalar_entry = [&](std::vector<uuber::ProtoParamEntry> &entries,
                                 const std::string &name, double value,
                                 bool logical_scalar = false) {
    if (!std::isfinite(value)) {
      return;
    }
    uuber::ProtoParamEntry entry;
    entry.name = name;
    if (logical_scalar) {
      entry.tag = uuber::ParamValueTag::LogicalScalar;
      entry.logical_scalar = static_cast<int>(value != 0.0);
    } else {
      entry.tag = uuber::ParamValueTag::NumericScalar;
      entry.numeric_scalar = value;
    }
    upsert_param_entry(entries, std::move(entry));
  };

  auto append_numeric_vector_entry =
      [&](std::vector<uuber::ProtoParamEntry> &entries, const std::string &name,
          const Rcpp::NumericVector &vec) {
        if (vec.size() == 0) {
          return;
        }
        uuber::ProtoParamEntry entry;
        entry.name = name;
        entry.tag = uuber::ParamValueTag::NumericVector;
        entry.numeric_values.reserve(vec.size());
        for (double val : vec) {
          entry.numeric_values.push_back(val);
        }
        upsert_param_entry(entries, std::move(entry));
      };

  auto append_logical_vector_entry =
      [&](std::vector<uuber::ProtoParamEntry> &entries, const std::string &name,
          const Rcpp::LogicalVector &vec) {
        if (vec.size() == 0) {
          return;
        }
        uuber::ProtoParamEntry entry;
        entry.name = name;
        entry.tag = uuber::ParamValueTag::LogicalVector;
        entry.logical_values.reserve(vec.size());
        for (int val : vec) {
          entry.logical_values.push_back(val);
        }
        upsert_param_entry(entries, std::move(entry));
      };

  for (R_xlen_t i = 0; i < n; ++i) {
    int acc_idx = -1;
    if (has_acc) {
      if (!Rf_isNumeric(acc_col_obj)) {
        Rcpp::stop("Column 'accumulator' must be numeric (1-based index)");
      }
      Rcpp::NumericVector acc_num(acc_col_obj);
      if (i >= acc_num.size()) {
        continue;
      }
      const double raw_val = acc_num[i];
      if (Rcpp::NumericVector::is_na(raw_val)) {
        continue;
      }
      const int idx_val = static_cast<int>(std::llround(raw_val)) - 1;
      if (idx_val < 0 || idx_val >= static_cast<int>(ctx.accumulators.size())) {
        continue;
      }
      acc_idx = idx_val;
    }
    if (acc_idx < 0 || acc_idx >= static_cast<int>(ctx.accumulators.size())) {
      continue;
    }
    const uuber::NativeAccumulator &base = ctx.accumulators[acc_idx];

    TrialAccumulatorParams override;
    override.onset_kind = base.onset_kind;
    override.onset_source_acc_idx = base.onset_source_acc_idx;
    override.onset_source_pool_idx = base.onset_source_pool_idx;
    override.onset_lag = base.onset_lag;
    override.onset =
        (i < onset_col.size() && !Rcpp::NumericVector::is_na(onset_col[i]))
            ? static_cast<double>(onset_col[i])
            : base.onset;
    bool has_shared = false;
    if (i < shared_col.size() && shared_col[i] != NA_STRING) {
      override.shared_trigger_id = Rcpp::as<std::string>(shared_col[i]);
      has_shared = true;
    } else {
      override.shared_trigger_id = base.shared_trigger_id;
      has_shared = !override.shared_trigger_id.empty();
    }
    if (override.shared_trigger_id != base.shared_trigger_id) {
      params_set->shared_trigger_layout_matches_context = false;
    }
    if (has_shared) {
      double gate_q = base.q;
      if (i < shared_q_col.size() &&
          !Rcpp::NumericVector::is_na(shared_q_col[i])) {
        gate_q = static_cast<double>(shared_q_col[i]);
      }
      gate_q = clamp_probability(gate_q);
      override.shared_q = gate_q;
      override.q = 0.0;
    } else {
      if (i < q_col.size() && !Rcpp::NumericVector::is_na(q_col[i])) {
        override.q = clamp_probability(static_cast<double>(q_col[i]));
      } else {
        override.q = base.q;
      }
      override.shared_q = override.q;
    }
    SEXP comp_entry = (i < comps_col.size()) ? comps_col[i] : R_NilValue;
    if (comp_entry != R_NilValue) {
      override.components = string_vector_from_entry(comp_entry);
      override.has_components = !override.components.empty();
      if (override.has_components) {
        override.component_indices =
            component_indices_from_labels(ctx, override.components);
      }
    } else {
      override.has_components = false;
    }

    std::string dist_name = base.dist;
    if (i < dist_col.size() && dist_col[i] != NA_STRING) {
      dist_name = Rcpp::as<std::string>(dist_col[i]);
    }
    override.dist_cfg = base.dist_cfg;
    std::vector<uuber::ProtoParamEntry> param_entries;
    if (params_list_col.size() > 0) {
      SEXP param_entry =
          (i < params_list_col.size()) ? params_list_col[i] : R_NilValue;
      if (param_entry != R_NilValue) {
        Rcpp::List param_list(param_entry);
        std::vector<uuber::ProtoParamEntry> entries =
            params_from_rcpp(param_list);
        for (auto &entry : entries) {
          upsert_param_entry(param_entries, entry);
        }
      }
    }
    for (std::size_t pc = 0; pc < param_cols.size(); ++pc) {
      SEXP column = param_col_data[pc];
      const int type = param_col_types[pc];
      const std::string &col_name = param_cols[pc];
      switch (type) {
      case REALSXP: {
        Rcpp::NumericVector vec(column);
        if (i < vec.size()) {
          const double val = vec[i];
          if (!Rcpp::NumericVector::is_na(val)) {
            append_scalar_entry(param_entries, col_name, val, false);
          }
        }
        break;
      }
      case INTSXP: {
        Rcpp::IntegerVector vec(column);
        if (i < vec.size()) {
          const int val = vec[i];
          if (val != NA_INTEGER) {
            append_scalar_entry(param_entries, col_name,
                                static_cast<double>(val), false);
          }
        }
        break;
      }
      case LGLSXP: {
        Rcpp::LogicalVector vec(column);
        if (i < vec.size()) {
          const int val = vec[i];
          if (val != NA_LOGICAL) {
            append_scalar_entry(param_entries, col_name,
                                static_cast<double>(val), true);
          }
        }
        break;
      }
      case STRSXP:
        break;
      case VECSXP: {
        Rcpp::List vec(column);
        if (i >= vec.size()) {
          break;
        }
        SEXP cell = vec[i];
        if (cell == R_NilValue) {
          break;
        }
        if (Rf_isReal(cell)) {
          Rcpp::NumericVector nv(cell);
          if (nv.size() == 1) {
            const double val = nv[0];
            if (!Rcpp::NumericVector::is_na(val)) {
              append_scalar_entry(param_entries, col_name, val, false);
            }
          } else {
            append_numeric_vector_entry(param_entries, col_name, nv);
          }
        } else if (Rf_isInteger(cell)) {
          Rcpp::IntegerVector iv(cell);
          if (iv.size() == 1) {
            const int val = iv[0];
            if (val != NA_INTEGER) {
              append_scalar_entry(param_entries, col_name,
                                  static_cast<double>(val), false);
            }
          } else {
            Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(iv);
            append_numeric_vector_entry(param_entries, col_name, nv);
          }
        } else if (Rf_isLogical(cell)) {
          Rcpp::LogicalVector lv(cell);
          if (lv.size() == 1) {
            const int val = lv[0];
            if (val != NA_LOGICAL) {
              append_scalar_entry(param_entries, col_name,
                                  static_cast<double>(val), true);
            }
          } else {
            append_logical_vector_entry(param_entries, col_name, lv);
          }
        }
        break;
      }
      default:
        break;
      }
    }
    if (!param_entries.empty() || dist_name != base.dist) {
      if (param_entries.empty()) {
        continue;
      }
      override.dist_cfg = uuber::resolve_acc_params_entries(dist_name, param_entries);
    }

    override.has_override = true;
    params_set->acc_params[acc_idx] = std::move(override);
    params_set->has_any_override = true;
    params_set->soa_cache_valid = false;
  }

  refresh_trial_param_fingerprint(*params_set);
  return params_set;
}
