#include "pool_math.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace {

inline double clamp_probability(double value) {
  if (!std::isfinite(value)) {
    return 0.0;
  }
  if (value < 0.0) {
    return 0.0;
  }
  if (value > 1.0) {
    return 1.0;
  }
  return value;
}

inline double clamp_unit(double value) { return clamp_probability(value); }

struct PoolPolyScratch {
  std::vector<std::vector<double>> prefix;
  std::vector<std::vector<double>> suffix;
  std::vector<double> surv;
  std::vector<double> fail;
};

inline void ensure_prefix_shape(std::vector<std::vector<double>> &prefix,
                                std::size_t n) {
  if (prefix.size() != n + 1) {
    prefix.resize(n + 1);
  }
  for (std::size_t i = 0; i <= n; ++i) {
    if (prefix[i].size() != i + 1) {
      prefix[i].assign(i + 1, 0.0);
    }
  }
}

inline void ensure_suffix_shape(std::vector<std::vector<double>> &suffix,
                                std::size_t n) {
  if (suffix.size() != n + 1) {
    suffix.resize(n + 1);
  }
  for (std::size_t i = 0; i <= n; ++i) {
    const std::size_t len = (n - i) + 1;
    if (suffix[i].size() != len) {
      suffix[i].assign(len, 0.0);
    }
  }
}

inline PoolPolyScratch &fetch_pool_poly_scratch(std::size_t n,
                                                bool need_suffix) {
  thread_local std::unordered_map<std::size_t, PoolPolyScratch> scratch_map;
  PoolPolyScratch &scratch = scratch_map[n];
  ensure_prefix_shape(scratch.prefix, n);
  if (need_suffix) {
    ensure_suffix_shape(scratch.suffix, n);
  }
  return scratch;
}

inline void ensure_prob_buffers(PoolPolyScratch &scratch, std::size_t n) {
  if (scratch.surv.size() != n) {
    scratch.surv.resize(n);
  }
  if (scratch.fail.size() != n) {
    scratch.fail.resize(n);
  }
}

inline void fill_prefix_buffers(std::vector<std::vector<double>> &prefix,
                                const std::vector<double> &surv,
                                const std::vector<double> &fail) {
  const std::size_t n = surv.size();
  std::fill(prefix[0].begin(), prefix[0].end(), 0.0);
  prefix[0][0] = 1.0;
  for (std::size_t i = 0; i < n; ++i) {
    const std::vector<double> &prev = prefix[i];
    std::vector<double> &out = prefix[i + 1];
    std::fill(out.begin(), out.end(), 0.0);
    const double surv_i = surv[i];
    const double fail_i = fail[i];
    for (std::size_t j = 0; j < prev.size(); ++j) {
      const double base = prev[j];
      if (base == 0.0) {
        continue;
      }
      out[j] += base * surv_i;
      out[j + 1] += base * fail_i;
    }
  }
}

inline void fill_suffix_buffers(std::vector<std::vector<double>> &suffix,
                                const std::vector<double> &surv,
                                const std::vector<double> &fail) {
  const std::size_t n = surv.size();
  std::fill(suffix[n].begin(), suffix[n].end(), 0.0);
  suffix[n][0] = 1.0;
  for (std::size_t idx = n; idx-- > 0;) {
    const std::vector<double> &next = suffix[idx + 1];
    std::vector<double> &out = suffix[idx];
    std::fill(out.begin(), out.end(), 0.0);
    const double surv_i = surv[idx];
    const double fail_i = fail[idx];
    for (std::size_t j = 0; j < next.size(); ++j) {
      const double base = next[j];
      if (base == 0.0) {
        continue;
      }
      out[j] += base * surv_i;
      out[j + 1] += base * fail_i;
    }
  }
}

inline double coefficient_for_order(const std::vector<double> &pref,
                                    const std::vector<double> &suff,
                                    int order) {
  if (order < 0) {
    return 0.0;
  }
  double total = 0.0;
  const int max_pref = static_cast<int>(pref.size()) - 1;
  const int max_index = std::min(order, max_pref);
  for (int i = 0; i <= max_index; ++i) {
    const int suffix_idx = order - i;
    if (suffix_idx < 0 || suffix_idx >= static_cast<int>(suff.size())) {
      continue;
    }
    total += pref[static_cast<std::size_t>(i)] *
             suff[static_cast<std::size_t>(suffix_idx)];
  }
  return total;
}

inline void fill_surv_fail_buffers(PoolPolyScratch &scratch,
                                   const double *survival, std::size_t n) {
  ensure_prob_buffers(scratch, n);
  std::vector<double> &surv = scratch.surv;
  std::vector<double> &fail = scratch.fail;
  for (std::size_t i = 0; i < n; ++i) {
    const double s = clamp_unit(survival[i]);
    surv[i] = s;
    fail[i] = clamp_unit(1.0 - s);
  }
}

inline double pool_density_fast_impl(const double *density,
                                     const double *survival, std::size_t n,
                                     int k) {
  if (n == 0 || k < 1 || k > static_cast<int>(n)) {
    return 0.0;
  }

  PoolPolyScratch &scratch = fetch_pool_poly_scratch(n, true);
  fill_surv_fail_buffers(scratch, survival, n);
  fill_prefix_buffers(scratch.prefix, scratch.surv, scratch.fail);
  fill_suffix_buffers(scratch.suffix, scratch.surv, scratch.fail);

  double total = 0.0;
  for (std::size_t idx = 0; idx < n; ++idx) {
    const double density_value = density[idx];
    if (!std::isfinite(density_value) || density_value <= 0.0) {
      continue;
    }
    const double coeff = coefficient_for_order(
        scratch.prefix[idx], scratch.suffix[idx + 1], k - 1);
    if (!std::isfinite(coeff) || coeff <= 0.0) {
      continue;
    }
    total += density_value * coeff;
  }
  if (!std::isfinite(total) || total < 0.0) {
    return 0.0;
  }
  return total;
}

inline double pool_survival_fast_impl(const double *survival, std::size_t n,
                                      int k) {
  if (n == 0) {
    return 1.0;
  }
  if (k < 1) {
    return 0.0;
  }

  PoolPolyScratch &scratch = fetch_pool_poly_scratch(n, false);
  fill_surv_fail_buffers(scratch, survival, n);
  fill_prefix_buffers(scratch.prefix, scratch.surv, scratch.fail);

  const std::vector<double> &poly = scratch.prefix.back();
  const int upto = std::min<int>(k, static_cast<int>(poly.size()));
  if (upto <= 0) {
    return 0.0;
  }

  double total = 0.0;
  for (int i = 0; i < upto; ++i) {
    total += poly[static_cast<std::size_t>(i)];
  }
  return clamp_probability(total);
}

} // namespace

double pool_density_fast(const std::vector<double> &density,
                         const std::vector<double> &survival, int k) {
  if (density.size() != survival.size()) {
    return 0.0;
  }
  return pool_density_fast_impl(density.data(), survival.data(), density.size(),
                                k);
}

double pool_survival_fast(const std::vector<double> &survival, int k) {
  return pool_survival_fast_impl(survival.data(), survival.size(), k);
}
