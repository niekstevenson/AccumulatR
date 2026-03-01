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

mass_nway_special="$(awk '
/double evaluate_coupling_mass_nway\(/ { in_fn=1 }
in_fn && /if \(acc_specialized\) \{/ { in_block=1 }
in_block { print }
in_block && /auto eval_unforced = \[\&\]\(const LabelRef &ref, double t\)/ { exit }
' src/distributions.cpp)"
if [[ -z "$mass_nway_special" ]]; then
  fail "Unable to locate acc-specialized block for evaluate_coupling_mass_nway."
fi
if printf "%s" "$mass_nway_special" | rg -n "integrate_boost_fn_0_to_upper" >/dev/null; then
  fail "evaluate_coupling_mass_nway acc-specialized block still uses scalar adaptive integration."
fi

mass_pair_special="$(awk '
/double evaluate_coupling_mass_pair\(/ { in_fn=1 }
in_fn && /if \(acc_specialized\) \{/ { in_block=1 }
in_block { print }
in_block && /auto eval_unforced = \[\&\]\(const LabelRef &ref, double t\)/ { exit }
' src/distributions.cpp)"
if [[ -z "$mass_pair_special" ]]; then
  fail "Unable to locate acc-specialized block for evaluate_coupling_mass_pair."
fi
if printf "%s" "$mass_pair_special" | rg -n "integrate_boost_fn_0_to_upper|auto total_integrand" >/dev/null; then
  fail "evaluate_coupling_mass_pair acc-specialized block still uses scalar adaptive integration."
fi

density_nway_special="$(awk '
/double evaluate_coupling_density_nway\(/ { in_fn=1 }
in_fn && /if \(acc_specialized\) \{/ { in_block=1 }
in_block { print }
in_block && /auto eval_unforced = \[\&\]\(const LabelRef &ref, double s\)/ { exit }
' src/distributions.cpp)"
if [[ -z "$density_nway_special" ]]; then
  fail "Unable to locate acc-specialized block for evaluate_coupling_density_nway."
fi
if printf "%s" "$density_nway_special" | rg -n "integrate_boost_fn_0_to_upper|auto order_integrand" >/dev/null; then
  fail "evaluate_coupling_density_nway acc-specialized block still uses scalar adaptive integration."
fi

density_pair_special="$(awk '
/double evaluate_coupling_density_pair\(/ { in_fn=1 }
in_fn && /if \(acc_specialized\) \{/ { in_block=1 }
in_block { print }
in_block && /auto eval_unforced = \[\&\]\(const LabelRef &ref, double s\)/ { exit }
' src/distributions.cpp)"
if [[ -z "$density_pair_special" ]]; then
  fail "Unable to locate acc-specialized block for evaluate_coupling_density_pair."
fi
if printf "%s" "$density_pair_special" | rg -n "integrate_boost_fn_0_to_upper|auto i_integrand" >/dev/null; then
  fail "evaluate_coupling_density_pair acc-specialized block still uses scalar adaptive integration."
fi

echo "[OK] Step 5 structural checks passed."
