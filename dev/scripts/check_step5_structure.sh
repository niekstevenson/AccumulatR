#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"
cd "$ROOT_DIR"

fail() {
  echo "[FAIL] $1" >&2
  exit 1
}

if ! rg -n "build_time_batch_0_to_upper\\(" src/quadrature_batch.h src/quadrature_batch.cpp >/dev/null; then
  fail "Missing build_time_batch_0_to_upper finite/infinite batch builder."
fi

if ! rg -n "evaluate_coupling_acc_batch\\(" src/distributions.cpp >/dev/null; then
  fail "Missing vectorized accumulator batch evaluator for coupling paths."
fi

if ! rg -n "uuber::eval_pdf_vec\\(|uuber::eval_cdf_vec\\(" src/distributions.cpp >/dev/null; then
  fail "Vectorized distribution primitives are not referenced in coupling path."
fi

unified_block="$(awk '
/double evaluate_outcome_coupling_unified\(/ { in_fn=1 }
in_fn { print }
in_fn && /^double evaluate_coupling_mass_nway\(/ { exit }
' src/distributions.cpp)"
if [[ -z "$unified_block" ]]; then
  fail "Unable to locate evaluate_outcome_coupling_unified definition block."
fi

if printf "%s" "$unified_block" | rg -n "integrate_boost_fn_0_to_upper" >/dev/null; then
  fail "Unified coupling evaluator still uses scalar adaptive integration."
fi

if ! printf "%s" "$unified_block" | rg -n "evaluate_labelref_batch\\(" >/dev/null; then
  fail "Unified coupling evaluator is not using shared batch label evaluation."
fi

if ! printf "%s" "$unified_block" | rg -n "integrate_node_density_batch\\(" >/dev/null; then
  fail "Unified coupling evaluator is missing generic node-integral batch path."
fi

echo "[OK] Step 5 structural checks passed."
