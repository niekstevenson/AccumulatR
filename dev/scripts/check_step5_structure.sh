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

if ! rg -n "struct KernelBatchRuntimeState|KernelEventBatchEvalFn|KernelGuardBatchEvalFn|eval_kernel_node_batch_incremental\\(|eval_kernel_nodes_batch_incremental\\(" src/kernel_executor.h >/dev/null; then
  fail "Missing kernel-level batch executor contract in kernel_executor.h."
fi

if ! rg -n "reset_kernel_batch_runtime\\(|invalidate_kernel_batch_runtime_from_slot\\(|eval_kernel_node_batch_incremental\\(|eval_kernel_nodes_batch_incremental\\(" src/kernel_executor.cpp >/dev/null; then
  fail "Missing kernel-level batch executor implementation in kernel_executor.cpp."
fi

if ! rg -n "eval_kernel_competitor_product_batch_incremental\\(" src/kernel_executor.h src/kernel_executor.cpp >/dev/null; then
  fail "Missing executor-side compiled competitor batch product contract."
fi

unified_block="$(awk '
/^double evaluate_outcome_coupling_unified\(/ { in_fn=1 }
in_fn {
  print
  open_count += gsub(/\{/, "{")
  close_count += gsub(/\}/, "}")
  if (open_count > 0 && open_count == close_count) {
    exit
  }
}
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

if printf "%s" "$unified_block" | rg -n "eval_node_with_forced_state_view\\(" >/dev/null; then
  fail "Unified coupling evaluator still uses scalar direct NodeRef helper."
fi

if ! printf "%s" "$unified_block" | rg -n "eval_node_with_forced_state_view_batch\\(|evaluate_generic_direct_cdf_resolved_provider\\(|provider_runtime\\.eval_direct_cdf\\(" >/dev/null; then
  fail "Unified coupling evaluator is not routing direct NodeRef through the batch executor adapter."
fi

if printf "%s" "$unified_block" | rg -n "evaluate_generic_direct_cdf_resolved_provider\\(|provider_runtime\\.eval_direct_cdf\\(" >/dev/null; then
  noderef_cdf_block="$(awk '
/^inline double evaluate_generic_direct_cdf_resolved_provider\(/ { in_fn=1 }
in_fn {
  print
  open_count += gsub(/\{/, "{")
  close_count += gsub(/\}/, "}")
  if (open_count > 0 && open_count == close_count) {
    exit
  }
}
' src/distributions.cpp)"
  if [[ -z "$noderef_cdf_block" ]]; then
    fail "Unable to locate evaluate_generic_direct_cdf_noderef definition block."
  fi
  if printf "%s" "$noderef_cdf_block" | rg -n "eval_node_with_forced_state_view\\(" >/dev/null; then
    fail "Resolved NodeRef direct-CDF path still uses scalar direct helper."
  fi
  if ! printf "%s" "$noderef_cdf_block" | rg -n "eval_node_with_forced_state_view_batch\\(" >/dev/null; then
    fail "Resolved NodeRef direct-CDF path is not using the batch executor adapter."
  fi
fi

competitor_start_line="$(rg -n "^bool node_density_with_competitors_batch_internal\\(" src/distributions.cpp | tail -n 1 | cut -d: -f1)"
competitor_end_line="$(awk "NR>${competitor_start_line} && /^double node_density_with_competitors_internal\\(/ { print NR - 1; exit }" src/distributions.cpp)"
competitor_block="$(sed -n "${competitor_start_line},${competitor_end_line}p" src/distributions.cpp)"
if [[ -z "$competitor_block" ]]; then
  fail "Unable to locate node_density_with_competitors_batch_internal definition block."
fi

if ! printf "%s" "$competitor_block" | rg -n "eval_node_batch_with_state_dense\\(" >/dev/null; then
  fail "Competitor-bearing batch path is not using batched base-node execution."
fi

if ! printf "%s" "$competitor_block" | rg -n "uuber::eval_kernel_competitor_product_batch_incremental\\(" >/dev/null; then
  fail "Competitor-bearing batch path is not delegating compiled competitor ops to the executor layer."
fi

if printf "%s" "$competitor_block" | rg -n "eval_node_recursive\\(" >/dev/null; then
  fail "Competitor-bearing batch path still calls scalar eval_node_recursive inside the time-batch helper."
fi

if printf "%s" "$competitor_block" | rg -n "run_competitor_compiled_ops_product\\(" >/dev/null; then
  fail "Competitor-bearing batch path still calls scalar competitor compiled-op execution inside the time-batch helper."
fi

if printf "%s" "$competitor_block" | rg -n "run_competitor_compiled_ops_product_batch\\(|run_competitor_batch_survival_op_batch\\(|eval_nodes_batch_with_state_dense\\(" >/dev/null; then
  fail "Deleted helper-local batch orchestration is still present inside the competitor batch helper."
fi

if rg -n "^bool run_competitor_batch_survival_op_batch\\(|^bool run_competitor_compiled_ops_product_batch\\(|^bool eval_nodes_batch_with_state_dense\\(" src/distributions.cpp >/dev/null; then
  fail "Deleted helper-local batch orchestration symbols still exist in distributions.cpp."
fi

echo "[OK] Step 5 structural checks passed."
