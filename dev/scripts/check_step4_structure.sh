#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"
cd "$ROOT_DIR"

fail() {
  echo "[FAIL] $1" >&2
  exit 1
}

if rg -n "enforce_ranked_sequence|needs_ranked_scratch|ranked_trigger_scratch" src/distributions.cpp >/dev/null; then
  fail "Ranked-only execution/scratch flags are still referenced in hot paths."
fi

seq_block="$(awk '
/double evaluate_sequence_density_kernel_idx\(/ { ++n }
n >= 1 { print }
n >= 1 && /^double evaluate_trial_contribution_kernel_idx\(/ { exit }
' src/distributions.cpp)"

if [[ -z "$seq_block" ]]; then
  fail "Unable to locate evaluate_sequence_density_kernel_idx definition block."
fi

if printf "%s" "$seq_block" | rg -n "single_outcome_idx|single_node_idx|single_time|ranked_outcome_indices|ranked_node_indices|ranked_times_ptr" >/dev/null; then
  fail "Legacy single/ranked split fields are still referenced in sequence evaluator."
fi

if ! printf "%s" "$seq_block" | rg -n "sequence_length == 1" >/dev/null; then
  fail "Unified sequence evaluator missing single-step dispatch by sequence_length."
fi

if ! printf "%s" "$seq_block" | rg -n "sequence_prefix_density_resolved" >/dev/null; then
  fail "Unified sequence evaluator is not routing multi-step sequence density through the shared transition path."
fi

key_block="$(awk '
/struct SequenceStateKey \{/ { ++n }
n >= 1 { print }
n >= 1 && /^};$/ { exit }
' src/distributions.cpp)"

if [[ -z "$key_block" ]]; then
  fail "Unable to locate SequenceStateKey definition block."
fi

if printf "%s" "$key_block" | rg -n "std::string|string" >/dev/null; then
  fail "SequenceStateKey contains string fields; hot-path keying must be hash/bitset based."
fi

if ! rg -n "sequence_state_key\(" src/distributions.cpp >/dev/null; then
  fail "sequence_state_key helper missing."
fi

echo "[OK] Step 4 structural checks passed."
