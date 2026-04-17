#pragma once

#include <cctype>
#include <cstdint>
#include <string_view>

namespace accumulatr::leaf {

enum class DistKind : std::uint8_t {
  Lognormal = 0,
  Gamma = 1,
  Exgauss = 2,
  LBA = 3,
  RDM = 4
};

constexpr std::string_view to_string(DistKind kind) noexcept {
  switch (kind) {
  case DistKind::Lognormal:
    return "lognormal";
  case DistKind::Gamma:
    return "gamma";
  case DistKind::Exgauss:
    return "exgauss";
  case DistKind::LBA:
    return "LBA";
  case DistKind::RDM:
    return "RDM";
  }
  return "unknown";
}

constexpr int dist_param_count(DistKind kind) noexcept {
  switch (kind) {
  case DistKind::Lognormal:
    return 2;
  case DistKind::Gamma:
    return 2;
  case DistKind::Exgauss:
    return 3;
  case DistKind::LBA:
    return 4;
  case DistKind::RDM:
    return 4;
  }
  return 0;
}

inline bool try_parse_dist_kind(std::string_view text, DistKind *out) noexcept {
  auto eq_lower = [text](std::string_view ref) {
    if (text.size() != ref.size()) {
      return false;
    }
    for (std::size_t i = 0; i < text.size(); ++i) {
      const unsigned char lhs = static_cast<unsigned char>(text[i]);
      const unsigned char rhs = static_cast<unsigned char>(ref[i]);
      if (std::tolower(lhs) != std::tolower(rhs)) {
        return false;
      }
    }
    return true;
  };

  DistKind parsed{};
  if (eq_lower("lognormal")) {
    parsed = DistKind::Lognormal;
  } else if (eq_lower("gamma")) {
    parsed = DistKind::Gamma;
  } else if (eq_lower("exgauss")) {
    parsed = DistKind::Exgauss;
  } else if (eq_lower("lba")) {
    parsed = DistKind::LBA;
  } else if (eq_lower("rdm")) {
    parsed = DistKind::RDM;
  } else {
    return false;
  }

  if (out != nullptr) {
    *out = parsed;
  }
  return true;
}

} // namespace accumulatr::leaf
