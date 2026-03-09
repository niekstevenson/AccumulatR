# Full Plan 3 Sub-Plan: Step 1 Bitset-Only Forced State Core

Date: 2026-03-06
Status: COMPLETE (`FP1.0`..`FP1.4` achieved)
Parent: `dev/full_plan3.MD`
Batch ID range: `FP1.*`

## 1) Objective

Make hot evaluator internals bitset-only for forced state, with no set-based
forced-state plumbing in hot paths.

## 2) Current Anchors (Where Work Applies)

1. Event eval hot chain:
   1. `eval_event_ref_idx` at `src/distributions.cpp:2308`
2. Forced-state carrier in runtime state:
   1. `NodeEvalState` at `src/distributions.cpp:2879`
3. Guard input forced-state plumbing:
   1. `GuardEvalInput` at `src/distributions.cpp:3004`
4. Forced-node evaluation entry:
   1. `eval_node_with_forced_dense_bits` at `src/distributions.cpp:3021`,
      `src/distributions.cpp:4416`
5. Native outcome probability path (bitset boundary conversion):
   1. `native_outcome_probability_bits_impl_idx` at
      `src/distributions.cpp:6203`

## 3) Scope Boundaries

In scope:

1. C++ evaluator internals and callback call chains.
2. Internal state carriers and function signatures.
3. Bitset validity/scoping semantics.

Out of scope:

1. R API signatures.
2. Data-frame input schema.

## 4) Execution Batches

### FP1.0 Baseline + Probe Freeze

1. Capture baseline grep snapshots for forced-state patterns.
2. Capture benchmark + golden baseline under `full_plan3_FP1_...` prefix.

Must record:

1. `rg -n "unordered_set<int>|forced_contains_scoped|forced_set_intersects_scope" src/distributions.cpp`
2. `rg -n "forced_complete_bits|forced_survive_bits|GuardEvalInput|NodeEvalState" src/distributions.cpp`

### FP1.1 Introduce Single Forced-State Internal Carrier

1. Introduce one explicit internal forced-state carrier/view used across
   `NodeEvalState`, `GuardEvalInput`, and event/node evaluators.
2. Ensure carrier supports:
   1. complete bitset
   2. survive bitset
   3. validity flags
   4. scope filter metadata
3. Do not change R-facing API contracts.

### FP1.2 Remove Mixed Signature Plumbing In Hot Chain

1. Remove multi-argument forced-state signature duplication in hot chain where
   carrier can be passed once.
2. Target calls in:
   1. `eval_event_ref_idx`
   2. `eval_node_with_forced_dense_bits`
   3. guard callback plumbing (`make_kernel_guard_eval` path)
3. Keep only boundary conversion from API vectors to bitsets at edge wrappers.

### FP1.3 Remove Set-Based Fallback Helpers From Hot Path

1. Ensure no set-based forced-state helper can execute in hot evaluator loops.
2. If helper exists for compatibility, keep it isolated to boundary/setup code,
   not runtime integrands.

### FP1.4 Verification + Evidence

1. Run full verification loop from parent contract.
2. Capture split-point and forced-state grep evidence.

## 5) Acceptance (What Must Be True)

1. `rg "unordered_set<int>" src/distributions.cpp` matches only setup/
   boundary-compatible code (or zero results).
2. Hot evaluator call chain passes forced state through one unified internal
   carrier, not duplicated argument packs.
3. Kernel event and guard callbacks operate on bitset-backed forced state only.
4. Golden passes with `max_abs_diff <= 1e-4`.
5. No performance gate violation introduced by Step 1 changes.

## 6) What Counts

1. Signature simplification in hot chain with real callsite updates.
2. Deletion or hard isolation of set-based fallback behavior from runtime loops.
3. Evidence proving hot path behavior, not just compile success.

## 7) What Does Not Count

1. Adding wrapper adapters while keeping old mixed signature path active.
2. Leaving set-based helpers reachable from hot evaluator callbacks.
3. Claiming completion from grep-only evidence without correctness/perf checks.

## 8) Hard Stops

1. Any correctness drift above tolerance.
2. Any new hot-path allocation/check bloat introduced by forced-state rewrite.
3. Any mode-specific branching added to compensate for Step 1 refactor.

## 9) Required Evidence Bundle (`FP1.*`)

1. `full_plan3_FP1_*_grep_forcedstate.log`
2. `full_plan3_FP1_*_step3_structure.log`
3. `full_plan3_FP1_*_step4_structure.log`
4. `full_plan3_FP1_*_step5_structure.log`
5. `full_plan3_FP1_*_step6_structure.log`
6. `full_plan3_FP1_*_golden.log`
7. `full_plan3_FP1_*_bench.log`
8. `full_plan3_FP1_*_diff.csv`
9. `full_plan3_FP1_*_summary.txt`

## 10) Exit Rule

Step 1 is complete only when Sections 5 and 9 are satisfied on one commit
lineage.

## 11) Evidence Block: FP1.0 Baseline + Probe Freeze

1. Timestamp (UTC): `20260306T130204Z`
2. Commit SHA: `54437eb`
3. Batch: `FP1`, subbatch `FP1.0`
4. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
6. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. worst ratio: `0.9318182` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. Step 1 probe facts:
   1. `unordered_set<int>`/set-fallback grep produced no matches in
      `src/distributions.cpp` (`required_pattern_1` empty).
   2. Forced-state bitset plumbing snapshot captured for baseline in
      `NodeEvalState`, `GuardEvalInput`, and hot-chain call sites.
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_grep_forcedstate.log`
   4. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_step1_hotchain_probe.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T130204Z_summary.txt`

## 12) Evidence Block: FP1.1 Unified Forced-State Carrier

1. Timestamp (UTC): `20260306T131150Z`
2. Commit SHA: `54437eb`
3. Batch: `FP1`, subbatch `FP1.1`
4. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
6. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. worst ratio: `0.9204545` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP1.1 carrier proof:
   1. `ForcedStateView` introduced and used in `NodeEvalState`,
      `GuardEvalInput`, `eval_event_ref_idx`, and kernel event/guard
      callbacks.
   2. Set-fallback grep remains empty in `src/distributions.cpp`.
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_grep_forcedstate.log`
   4. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_carrier_probe.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T131150Z_summary.txt`

## 13) Evidence Block: FP1.2 Mixed Hot-Chain Signature Removal

1. Timestamp (UTC): `20260306T132012Z`
2. Commit SHA: `54437eb`
3. Batch: `FP1`, subbatch `FP1.2`
4. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
6. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. worst ratio: `0.9431818` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP1.2 proof:
   1. Ranked transition hot loop now calls
      `eval_node_with_forced_state_view` and
      `make_guard_input_forced_state` instead of threading separate
      forced-bit arguments.
   2. Unified coupling generic NodeRef fast branch now uses carrier-based
      node-eval helper.
   3. Legacy wrapper `eval_node_with_forced_dense_bits` remains as
      compatibility shim; hot callsites route via carrier helper.
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_grep_forcedstate.log`
   4. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_mixed_signature_probe.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132012Z_summary.txt`

## 14) Evidence Block: FP1.3 Set-Fallback Deletion/Isolation Check

1. Timestamp (UTC): `20260306T132542Z`
2. Commit SHA: `54437eb`
3. Batch: `FP1`, subbatch `FP1.3`
4. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
6. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. worst ratio: `0.9469697` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP1.3 proof:
   1. No `unordered_set<int>` or set-based forced helper symbols in
      `src/distributions.cpp` (`required_pattern_1` empty).
   2. Hot call chain remains on carrier-based helpers
      (`eval_node_with_forced_state_view`,
      `make_guard_input_forced_state`).
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_grep_forcedstate.log`
   4. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_set_fallback_probe.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T132542Z_summary.txt`

## 15) Evidence Block: FP1.4 Step 1 Closure Decision + Rollup

1. Timestamp (UTC): `20260306T133234Z`
2. Commit SHA: `54437eb`
3. Batch: `FP1`, subbatch `FP1.4`
4. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
6. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. worst ratio: `0.9280303` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. Closure decision:
   1. decision: `ACHIEVED`
   2. criteria:
      1. `unordered_set<int>` hotpath absence: pass
      2. unified carrier in hot chain: pass
      3. kernel callbacks bitset-only forced-state handling: pass
      4. golden gate: pass
      5. perf gate: pass
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_grep_forcedstate.log`
   4. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_step1_closure_probe.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_step1_closure_decision.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_split_inventory.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_step3_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_step4_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_step5_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_step6_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_golden.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_bench.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_benchmark_centralized.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_benchmark_centralized_diff.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_benchmark_step3_comparison.csv`
   16. `dev/scripts/scratch_outputs/full_plan3_FP1_20260306T133234Z_summary.txt`
