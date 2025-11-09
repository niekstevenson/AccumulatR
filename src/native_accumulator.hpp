#pragma once

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <vector>

#include "native_proto.hpp"

namespace uuber {

enum AccDistKind {
  ACC_DIST_LOGNORMAL = 1,
  ACC_DIST_GAMMA = 2,
  ACC_DIST_EXGAUSS = 3
};

struct AccDistParams {
  int code{ACC_DIST_LOGNORMAL};
  double p1{};
  double p2{};
  double p3{};
};

inline std::string normalize_dist_name(const std::string& dist) {
  std::string out(dist);
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return out;
}

inline const ProtoParamEntry* find_param(const std::vector<ProtoParamEntry>& entries,
                                         const std::string& name) {
  for (const auto& entry : entries) {
    if (entry.name == name) {
      return &entry;
    }
  }
  return nullptr;
}

inline double extract_param_numeric(const ProtoParamEntry* entry,
                                    const std::string& dist_name,
                                    const std::string& param_name) {
  if (entry == nullptr) {
    throw std::runtime_error("Accumulator distribution '" + dist_name +
                             "' missing parameter '" + param_name + "'");
  }
  if (entry->tag == ParamValueTag::NumericScalar) {
    return entry->numeric_scalar;
  }
  if (entry->tag == ParamValueTag::NumericVector && !entry->numeric_values.empty()) {
    return entry->numeric_values.front();
  }
  throw std::runtime_error("Accumulator distribution '" + dist_name +
                           "' parameter '" + param_name + "' must be numeric");
}

inline AccDistParams resolve_acc_params_entries(const std::string& dist,
                                                const std::vector<ProtoParamEntry>& entries) {
  if (entries.empty()) {
    throw std::runtime_error("Accumulator distribution '" + dist + "' missing parameter list");
  }
  AccDistParams cfg{};
  std::string dist_name = normalize_dist_name(dist);
  if (dist_name == "lognormal") {
    cfg.code = ACC_DIST_LOGNORMAL;
    cfg.p1 = extract_param_numeric(find_param(entries, "meanlog"), dist_name, "meanlog");
    cfg.p2 = extract_param_numeric(find_param(entries, "sdlog"), dist_name, "sdlog");
    cfg.p3 = 0.0;
  } else if (dist_name == "gamma") {
    cfg.code = ACC_DIST_GAMMA;
    cfg.p1 = extract_param_numeric(find_param(entries, "shape"), dist_name, "shape");
    cfg.p2 = extract_param_numeric(find_param(entries, "rate"), dist_name, "rate");
    cfg.p3 = 0.0;
  } else if (dist_name == "exgauss") {
    cfg.code = ACC_DIST_EXGAUSS;
    cfg.p1 = extract_param_numeric(find_param(entries, "mu"), dist_name, "mu");
    cfg.p2 = extract_param_numeric(find_param(entries, "sigma"), dist_name, "sigma");
    cfg.p3 = extract_param_numeric(find_param(entries, "tau"), dist_name, "tau");
  } else {
    throw std::runtime_error("Unsupported accumulator distribution '" + dist_name + "'");
  }
  return cfg;
}

} // namespace uuber
