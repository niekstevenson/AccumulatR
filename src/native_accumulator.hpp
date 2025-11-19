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
  double t0{0.0};
};

inline std::string normalize_dist_name(const std::string& dist) {
  std::string out(dist);
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return out;
}

inline double numeric_param(const std::vector<ProtoParamEntry>& entries,
                            const std::string& name,
                            double default_value = 0.0) {
  for (const auto& entry : entries) {
    if (entry.name == name) {
      if (entry.tag == ParamValueTag::NumericScalar) {
        return entry.numeric_scalar;
      }
      if (!entry.numeric_values.empty()) {
        return entry.numeric_values.front();
      }
    }
  }
  return default_value;
}

inline AccDistParams resolve_acc_params_entries(const std::string& dist,
                                                const std::vector<ProtoParamEntry>& entries) {
  AccDistParams cfg{};
  std::string dist_name = normalize_dist_name(dist);
  if (dist_name == "lognormal") {
    cfg.code = ACC_DIST_LOGNORMAL;
    cfg.p1 = numeric_param(entries, "meanlog");
    cfg.p2 = numeric_param(entries, "sdlog");
  } else if (dist_name == "gamma") {
    cfg.code = ACC_DIST_GAMMA;
    cfg.p1 = numeric_param(entries, "shape");
    cfg.p2 = numeric_param(entries, "rate");
  } else if (dist_name == "exgauss") {
    cfg.code = ACC_DIST_EXGAUSS;
    cfg.p1 = numeric_param(entries, "mu");
    cfg.p2 = numeric_param(entries, "sigma");
    cfg.p3 = numeric_param(entries, "tau");
  }
  cfg.t0 = numeric_param(entries, "t0");
  return cfg;
}

} // namespace uuber
