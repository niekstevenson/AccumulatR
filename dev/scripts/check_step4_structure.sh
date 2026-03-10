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

if rg -n "sequence_times_ptr|sequence_times_storage|sequence_times_vec" src/distributions.cpp >/dev/null; then
  fail "Legacy sequence-time payload/storage split is still present."
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

if printf "%s" "$seq_block" | rg -n "sequence_length == 1" >/dev/null; then
  fail "Unified sequence evaluator still dispatches by raw sequence_length."
fi

if ! printf "%s" "$seq_block" | rg -n "spec\\.sequence_execution" >/dev/null; then
  fail "Unified sequence evaluator missing execution-contract dispatch."
fi

if ! printf "%s" "$seq_block" | rg -n "TrialSequenceExecutionKind::LowerLayerDirect" >/dev/null; then
  fail "Unified sequence evaluator missing lower-layer direct execution contract."
fi

if ! printf "%s" "$seq_block" | rg -n "sequence_prefix_density_resolved" >/dev/null; then
  fail "Unified sequence evaluator is not routing multi-step sequence density through the shared transition path."
fi

if ! printf "%s" "$seq_block" | rg -n "TrialSequenceExecutionKind::SequenceState" >/dev/null; then
  fail "Unified sequence evaluator missing sequence-state execution contract."
fi

if printf "%s" "$seq_block" | rg -n "accumulate_component_guess_density_idx|step_guess_shortcut|allow_guess_shortcut" >/dev/null; then
  fail "Unified sequence evaluator still contains guess-shortcut specialization."
fi

if rg -n "step_guess_shortcut|allow_guess_shortcut" src/distributions.cpp >/dev/null; then
  fail "Guess-shortcut payload flags are still present in distributions.cpp."
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

exact_hash_block="$(awk '
/inline std::uint64_t ranked_hash_exact_times\(/ { ++n }
n >= 1 { print }
n >= 1 && /^}/ { exit }
' src/distributions.cpp)"

if [[ -z "$exact_hash_block" ]]; then
  fail "Unable to locate ranked_hash_exact_times definition block."
fi

if printf "%s" "$exact_hash_block" | rg -n "std::vector|std::sort" >/dev/null; then
  fail "ranked_hash_exact_times still materializes/sorts temporary vectors."
fi

bounds_hash_block="$(awk '
/ranked_hash_bounds\(const SourceTimeBoundsMap/ { ++n }
n >= 1 { print }
n >= 1 && /^}/ { exit }
' src/distributions.cpp)"

if [[ -z "$bounds_hash_block" ]]; then
  fail "Unable to locate ranked_hash_bounds definition block."
fi

if printf "%s" "$bounds_hash_block" | rg -n "std::vector|std::sort" >/dev/null; then
  fail "ranked_hash_bounds still materializes/sorts temporary vectors."
fi

if rg -n "std::vector<ExactSourceTimeMap> observed_prefix_exact|observed_prefix_exact\\[rank_idx \\+ 1\\] = observed_prefix_exact\\[rank_idx\\]" src/distributions.cpp >/dev/null; then
  fail "Ranked prefix exact times still use copied prefix-map staging."
fi

echo "[OK] Step 4 structural checks passed."
