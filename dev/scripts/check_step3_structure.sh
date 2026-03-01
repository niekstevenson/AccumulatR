#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"
cd "$ROOT_DIR"

fail() {
  echo "[FAIL] $1" >&2
  exit 1
}

if rg -n "KernelGuardEvalMode::GeneralCompiled" src/context.h src/prep_builder.cpp src/distributions.cpp >/dev/null; then
  fail "GeneralCompiled guard mode is still referenced."
fi

if rg -n "node_contains_guard\(" src/distributions.cpp >/dev/null; then
  fail "Runtime recursive guard-structure probe (node_contains_guard) is still present."
fi

guard_block="$(awk '
/double guard_cdf_internal\(const GuardEvalInput &input, double t,/ { ++n }
n >= 2 { print }
n >= 2 && /std::vector<int> gather_blocker_sources\(/ { exit }
' src/distributions.cpp)"

if [[ -z "$guard_block" ]]; then
  fail "Unable to locate guard_cdf_internal definition block."
fi

if printf "%s" "$guard_block" | rg -n "integrate_boost_fn|guard_density_internal\(input, u" >/dev/null; then
  fail "guard_cdf_internal still contains generic quadrature fallback logic."
fi

if ! printf "%s" "$guard_block" | rg -n "LinearChainODE" >/dev/null; then
  fail "guard_cdf_internal is not enforcing LinearChainODE mode."
fi

echo "[OK] Step 3 structural checks passed."
