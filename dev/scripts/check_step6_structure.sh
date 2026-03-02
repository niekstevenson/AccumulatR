#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"
cd "$ROOT_DIR"

fail() {
  echo "[FAIL] $1" >&2
  exit 1
}

if rg -n "if \(acc_specialized\)|eval_unforced = \[\&\]" src/distributions.cpp >/dev/null; then
  fail "Specialized-vs-generic fallback branches are still present in coupling evaluators."
fi

if rg -n "uses_coupling_program" src/distributions.cpp >/dev/null; then
  fail "Outcome probability path still branches on uses_coupling_program."
fi

if ! rg -n "evaluate_outcome_coupling_unified\(" src/distributions.cpp >/dev/null; then
  fail "Unified coupling evaluator entrypoint is missing."
fi

if ! rg -n "GenericNodeIntegral" src/context.h src/prep_builder.cpp src/distributions.cpp >/dev/null; then
  fail "GenericNodeIntegral coupling kind is not wired through context/prep/runtime."
fi

prob_block="$(awk '
/double native_outcome_probability_bits_impl_idx\(/ { ++n }
n >= 1 { print }
n >= 1 && /^double native_outcome_probability_bits_impl_idx\(/ && ++m >= 2 { exit }
' src/distributions.cpp)"

if [[ -z "$prob_block" ]]; then
  fail "Unable to locate native_outcome_probability_bits_impl_idx definition block."
fi

if printf "%s" "$prob_block" | rg -n "integrate_boost_fn_0_to_upper" >/dev/null; then
  fail "Outcome probability function still contains scalar adaptive integration fallback."
fi

if ! printf "%s" "$prob_block" | rg -n "evaluate_outcome_coupling_unified\(" >/dev/null; then
  fail "Outcome probability function is not routing competitor handling via unified coupling evaluator."
fi

echo "[OK] Step 6 structural checks passed."
