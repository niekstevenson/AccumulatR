#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

constexpr std::uint64_t kFNV64Offset = 1469598103934665603ULL;
constexpr std::uint64_t kFNV64Prime = 1099511628211ULL;
constexpr std::uint64_t kCanonicalQuietNaNBits = 0x7ff8000000000000ULL;

template <typename T> inline T clamp(T val, T lo, T hi) {
  if (!std::isfinite(val)) {
    return val;
  }
  if (val < lo) {
    return lo;
  }
  if (val > hi) {
    return hi;
  }
  return val;
}

inline void sort_unique(std::vector<int> &vec) {
  if (vec.empty()) {
    return;
  }
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

inline std::vector<int> integer_vector_to_std(const Rcpp::IntegerVector &src,
                                              bool subtract_one = false) {
  std::vector<int> out;
  out.reserve(src.size());
  for (int value : src) {
    if (value == NA_INTEGER) {
      continue;
    }
    out.push_back(subtract_one ? value - 1 : value);
  }
  sort_unique(out);
  return out;
}

inline bool is_invalid_positive(double value) {
  return !std::isfinite(value) || value <= 0.0;
}

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

inline double safe_density(double value) {
  if (!std::isfinite(value) || value <= 0.0) {
    return 0.0;
  }
  return value;
}

inline std::uint64_t canonical_double_bits(double value) {
  if (value == 0.0) {
    return 0ULL;
  }
  if (std::isnan(value)) {
    return kCanonicalQuietNaNBits;
  }
  std::uint64_t bits = 0ULL;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline void hash_append_bytes(std::uint64_t &hash, const void *data,
                              std::size_t n_bytes) {
  const auto *bytes = static_cast<const unsigned char *>(data);
  for (std::size_t i = 0; i < n_bytes; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= kFNV64Prime;
  }
}

inline void hash_append_u64(std::uint64_t &hash, std::uint64_t value) {
  hash_append_bytes(hash, &value, sizeof(value));
}

inline void hash_append_bool(std::uint64_t &hash, bool value) {
  const std::uint8_t v = value ? 1U : 0U;
  hash_append_bytes(hash, &v, sizeof(v));
}

inline void hash_append_double(std::uint64_t &hash, double value) {
  hash_append_u64(hash, canonical_double_bits(value));
}

inline std::uint64_t mix_hash64(std::uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}
