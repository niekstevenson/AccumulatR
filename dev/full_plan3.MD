---
editor_options:
  markdown:
    wrap: 72
---

# Phase 3 + 6 Hard-Cutover Execution Contract (No Legacy Fallback)

Last updated: 2026-03-08T16:45:41Z (UTC)
Status: OPEN (RESET UNDER THIS CONTRACT)

## 1) Brutal Current State

1. This contract is not complete.
2. Prior Step 6 evidence exists, but it does not close the full hard-cutover
   criteria defined here.
3. Active structural split points still exist in runtime execution:
   1. Unified coupling dispatcher still routes into pair/nway helper engines:
      `src/distributions.cpp:4238`, `src/distributions.cpp:4243`
   2. Provider-kind branch still controls runtime coupling behavior:
      `src/distributions.cpp:4120`, `src/distributions.cpp:4136`,
      `src/distributions.cpp:4146`
   3. Non-coupling path still evaluates base node density through per-time
      recursive node calls:
      `src/distributions.cpp:5986`, `src/distributions.cpp:6051`
   4. Competitor path still has incremental-batch vs per-node runtime split:
      `src/distributions.cpp:5808`, `src/distributions.cpp:5819`
4. Honest summary: this is partial consolidation, not full single-engine
   closure.

## 2) Scope (Exact)

This plan must deliver all six goals below on one commit lineage:

1. Step 1: bitset-only forced-state internals in hot evaluator chains.
2. Step 2: one shared kernel execution surface for observed/ranked/nonresponse.
3. Step 3: guard/competitor transition compilation path as runtime default.
4. Step 4: ranked flow fully on fused kernel machinery with non-string hot key
   path.
5. Step 5: vectorized distribution execution generalized beyond narrow
   linear-chain branch.
6. Step 6: legacy/fallback branch deletion with one active runtime engine.

## 3) Invariants (Must Hold Throughout)

1. No `std::unordered_set<int>` forced-state operations in hot evaluator paths.
2. No recursive/scenario fallback in hot loops for observed/ranked likelihood
   once the relevant step is declared complete.
3. No string serialization/lookups in integrand inner loops.
4. One runtime kernel execution surface across observed/ranked/nonresponse
   internals.
5. No new mode-specific runtime path added for speed recovery.

## 4) Completion Definition (No Interpretation)

Plan is complete only when all are true:

1. Steps 1 through 6 are individually marked complete under this contract.
2. Structural checks pass.
3. Golden correctness passes with `max_abs_diff <= 1e-4`.
4. Centralized benchmark passes:
   1. no workload ratio > `1.02`
   2. no P1 workload ratio > `1.00`
5. Stability passes 3 consecutive fixed-setting runs.
6. Deletion proof exists for legacy fallback/specialized branches.
7. Required evidence bundle exists for every batch.

## 5) What Counts vs What Fails

### Counts

1. Actual code-path deletion or hard-unreachable proof.
2. One runtime engine with provider/mode differences only as data inputs.
3. Measured correctness + performance + stability on same commit lineage.
4. Explicit split-point inventory shrinking to zero active hot-path splits.

### Automatic Fail (Does Not Count)

1. Wrapper around old specialized engines while keeping them reachable.
2. Labeling one of two active paths as "primary".
3. Vector claims while hot loops still do scalar fallback on target workloads.
4. Single benchmark run without stability series.
5. Completion claims without artifacts listed in Section 9.

## 6) Ordered Execution (Hard Sequence)

1. FP0 baseline freeze + instrumentation lock.
2. FP1 Step 1 implementation and proof.
3. FP2 Step 2 implementation and proof.
4. FP3 Step 3 implementation and proof.
5. FP4 Step 4 implementation and proof.
6. FP5 Step 5 implementation and proof.
7. FP6 Step 6 deletion and proof.
8. FP7 final closure gates.

No skipping. No out-of-order completion claims.

## 7) Sub-Plan Artifacts (Exhaustive Detail)

1. Step 1: `dev/full_plan3_step1_bitset_only_core.MD`
2. Step 2: `dev/full_plan3_step2_unified_kernel_entry.MD`
3. Step 3-4: `dev/full_plan3_step3_4_transitions_ranked_cutover.MD`
4. Step 4 lower-layer: `dev/full_plan3_step4_lower_layer_convergence.MD`
5. Step 5-6: `dev/full_plan3_step5_6_vectorization_legacy_deletion.MD`

These sub-plans are normative; this master file is the contract spine.

## 8) Verification Loop (Run After Every Batch)

1. Build:
   1. `R CMD INSTALL .`
2. Structure:
   1. `bash dev/scripts/check_step3_structure.sh`
   2. `bash dev/scripts/check_step4_structure.sh`
   3. `bash dev/scripts/check_step5_structure.sh`
   4. `bash dev/scripts/check_step6_structure.sh`
3. Correctness:
   1. `Rscript dev/scripts/check_loglik_golden.R`
4. Performance:
   1. `Rscript dev/scripts/bench_centralized.R`
5. Optional profile snapshot for slowdowns:
   1. `ACCUMULATR_PROFILE_RUN_SEC=18 ACCUMULATR_PROFILE_DURATION=18 ACCUMULATR_PROFILE_INTERVAL=1 ACCUMULATR_PROFILE_FILE=dev/scratch_profile_current.txt bash dev/scripts/profile_cpp_simple.sh`

## 9) Evidence Contract (Mandatory Per Batch)

Each batch must produce:

1. Commit SHA.
2. UTC timestamp.
3. Batch ID (`FP0`..`FP7`).
4. Structure logs.
5. Golden log + `max_abs_diff`.
6. Benchmark log + `benchmark_centralized.csv` + diff csv.
7. Stability logs (if FP7 or when performance risk requires).
8. Split-point inventory snapshot (`file:line`).
9. Deletion list (if branch deletion occurred).
10. Explicit declaration:
    `No wrapper-only consolidation introduced = yes/no`

Naming prefix:

1. `dev/scripts/scratch_outputs/full_plan3_<batch>_<timestamp>_*`

Missing artifact means batch invalid.

## 10) P1 Workload Set (Performance Gate)

1. `example_3_stop_na`
2. `example_21_simple_q`
3. `example_22_shared_q`
4. `example_6_dual_path`
5. `depth3_guard_competitor`
6. `shared_gate_nway_shared_triggers`

## 11) Hard Stop Conditions

1. Golden fails (`max_abs_diff > 1e-4`) -> stop and fix before continuing.
2. Any new mode/workload-specific runtime path appears -> reject batch.
3. Missing evidence bundle -> batch invalid.
4. Completion claim with non-empty active split list -> rejected.

## 12) Additional Consolidation/Speed Ideas (Allowed, Same Engine Only)

These are allowed only if they do not re-split architecture:

1. Hoist invariant validation/checking out of hot loops into one-time
   prevalidation.
2. Component/trial SoA scratch reuse to remove per-iteration allocations.
3. Pre-resolve node/label/provider references and reuse compact indices.
4. Adaptive batched integration policy inside one evaluator stack.
5. Fingerprint/cache-key narrowing to remove repeated expensive hashing.

## 13) Status Board (Under This Contract)

1. FP0 baseline freeze: COMPLETE (`20260306T124524Z`, commit `54437eb`)
2. FP1 Step 1: COMPLETE (`FP1.0` + `FP1.1` + `FP1.2` + `FP1.3` + `FP1.4`,
   closure decision `ACHIEVED`, commit `54437eb`)
3. FP2 Step 2: COMPLETE (`FP2.0` + `FP2.1` + `FP2.2` + `FP2.3` + `FP2.4` +
   `FP2.5`, closure decision `ACHIEVED`, commit `54437eb`)
4. FP3 Step 3: COMPLETE (`FP3.0` + `FP3.1` + `FP3.2` + `FP3.3`, closure
   decision `ACHIEVED`, commit `54437eb`)
5. FP4 Step 4: IN PROGRESS (`FP4.0` baseline/probe freeze complete; top-level
   `FP4.1` cut rejected; `FP4L.1` profiling complete; `FP4L.2` + `FP4L.3`
   keepable lower-layer convergence cut complete; `FP4L.4` planner
   convergence landed but Step 4 closure still blocked by residual
   `example_22_shared_q` P1 overhead, commit `54437eb`)
6. FP5 Step 5: NOT STARTED
7. FP6 Step 6: NOT STARTED
8. FP7 final closure: NOT STARTED

Global completion claim: NO

## 14) Immediate Next Action

1. Trim the residual `example_22_shared_q` single-step planner overhead so the
   strict P1 gate returns to `<= 1.00`.
2. Then execute `FP4L.5` closure/deletion-proof decision in
   `dev/full_plan3_step4_lower_layer_convergence.MD`.
3. Keep the failed `FP4.1` artifacts, the `FP4L.1` profile bundle, the
   keepable `FP4L.2`/`FP4L.3` bundle, and the `FP4L.4` planner-convergence
   bundle frozen as proof:
   `full_plan3_FP4_20260308T160703Z_*`,
   `full_plan3_FP4_20260308T161426Z_*`,
   `full_plan3_FP4L_20260308T164541Z_*`,
   `full_plan3_FP4L_20260308T185700Z_*`,
   `full_plan3_FP4L_20260309T091834Z_*`.

## 15) Evidence Block: FP0 Baseline Freeze

1. Timestamp (UTC): `20260306T124524Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9075342` (`example_2_stop_mixture`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. Required artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_split_inventory.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_grep_forcedstate.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_step3_structure.log`
   6. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_step4_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_step5_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_step6_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_golden.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_bench.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_benchmark_centralized.csv`
   12. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_benchmark_centralized_diff.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_benchmark_step3_comparison.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP0_20260306T124524Z_summary.txt`

## 16) Evidence Block: FP1.0 Baseline + Probe Freeze

1. Timestamp (UTC): `20260306T130204Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9318182` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. Step 1 baseline probe:
   1. required set-fallback grep produced no matches in
      `src/distributions.cpp`
   2. forced-state bitset and hot-chain baseline snapshots captured
7. Artifacts:
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

## 17) Evidence Block: FP1.1 Unified Forced-State Carrier

1. Timestamp (UTC): `20260306T131150Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9204545` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. Carrier proof summary:
   1. `ForcedStateView` added and integrated through `NodeEvalState`,
      `GuardEvalInput`, `eval_event_ref_idx`, kernel event callback, and
      kernel guard callback.
   2. required set-fallback grep remains empty in
      `src/distributions.cpp`.
7. Artifacts:
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

## 18) Evidence Block: FP1.2 Mixed Hot-Chain Signature Removal

1. Timestamp (UTC): `20260306T132012Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9431818` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP1.2 proof summary:
   1. Carrier helpers added:
      `eval_node_with_forced_state_view`,
      `make_guard_input_forced_state`.
   2. Ranked and kernel guard hot callsites now route through carrier helper
      entrypoints.
   3. Legacy wrapper remains as compatibility shim only.
7. Artifacts:
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

## 19) Evidence Block: FP1.3 Set-Fallback Deletion/Isolation Check

1. Timestamp (UTC): `20260306T132542Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9469697` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP1.3 proof summary:
   1. `unordered_set<int>` and forced set-fallback helper patterns remain
      absent from `src/distributions.cpp`.
   2. Hot-chain markers confirm carrier-based entrypoints are active.
7. Artifacts:
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

## 20) Evidence Block: FP1.4 Step 1 Closure Decision + Rollup

1. Timestamp (UTC): `20260306T133234Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9280303` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. Closure decision:
   1. decision: `ACHIEVED`
   2. criteria:
      1. `unordered_set<int>` hotpath absence: pass
      2. unified carrier in hot chain: pass
      3. kernel callbacks bitset-only forced-state handling: pass
      4. golden gate: pass
      5. perf gate: pass
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
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

## 21) Evidence Block: FP2.0 Baseline + Call Graph Snapshot

1. Timestamp (UTC): `20260306T134122Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9422348` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP2.0 call-graph freeze proof:
   1. Shared trial runtime kernel entrypoint exists at
      `evaluate_trial_probability_kernel_idx` and is called from `cpp_loglik`.
   2. Remaining Step 2 split points are explicitly captured in
      `mode_entryproof` and `split_inventory`.
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_summary.txt`

## 22) Evidence Block: FP2.1 Unified Internal Execution Contract

1. Timestamp (UTC): `20260306T140655Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9382716` (`example_23_ranked_chain`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP2.1 contract proof:
   1. Unified contribution contract implemented via `TrialKernelIntent` and
      `evaluate_trial_contribution_kernel_idx`.
   2. Trial result transform normalized via `TrialProbabilityTransform`
      (legacy `eval_input.nonresponse_mode` removed).
   3. FP2.0 split predicates removed:
      `spec.metric == KernelContributionMetric::SequenceDensity`,
      `spec.sequence_length == 1`,
      `spec.guess_shortcut`,
      `eval_input.nonresponse_mode`.
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_summary.txt`

## 23) Evidence Block: FP2.2 Observed Path Cutover On Unified Contract

1. Timestamp (UTC): `20260306T150007Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9822920` (`example_2_stop_mixture`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP2.2 cutover proof:
   1. Observed contributions now share one observed intent:
      `TrialKernelIntent::SequenceDensity`.
   2. Legacy observed intents removed:
      `TrialKernelIntent::NodeDensity`,
      `TrialKernelIntent::GuessDensity`,
      `TrialKernelIntent::SequencePrefixDensity`.
   3. One observed evaluator surface now handles both single-step and
      multi-step paths (`evaluate_sequence_density_kernel_idx`).
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_summary.txt`

## 24) Evidence Block: FP2.3 Ranked Path Cutover On Unified Contract

1. Timestamp (UTC): `20260306T151253Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9431818` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP2.3 cutover proof:
   1. Ranked-only execution toggle removed from `TrialEvalInput`.
   2. Ranked and observed paths share one sequence evaluator surface
      (`evaluate_sequence_density_kernel_idx`).
   3. Ranked behavior in execution is determined by `sequence_length`, not a
      separate ranked execution flag.
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_summary.txt`

## 25) Evidence Block: FP2.4 Nonresponse Path Cutover Cleanup

1. Timestamp (UTC): `20260306T151939Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9657534` (`example_2_stop_mixture`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP2.4 cutover proof:
   1. Runtime nonresponse mode-flag branches removed (`TrialEvalMode`,
      `eval_input.mode`, `eval_input.nonresponse_mode`).
   2. Nonresponse path now relies on normalized input shape
      (`sequence_length == 0`) plus complement probability transform.
   3. Caller-side nonresponse coupling-program binding is unconditional.
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_summary.txt`

## 26) Evidence Block: FP2.5 Step 2 Closure Decision + Deletion-Proof Rollup

1. Timestamp (UTC): `20260306T152414Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9602273` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP2.5 closure proof:
   1. `cpp_loglik` and `cpp_loglik_multiple` converge on one trial kernel
      entry (`evaluate_trial_probability_kernel_idx`).
   2. Observed/ranked/nonresponse contribution execution converges on one
      evaluator surface (`evaluate_trial_contribution_kernel_idx`).
   3. Mode differences are normalized as data (`TrialProbabilityTransform`,
      sequence shape) rather than separate evaluator engines.
7. Closure decision:
   1. `decision = ACHIEVED` for Step 2.
   2. Residual split points are within one engine and are carried into FP3+
      contraction work.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step2_closure_probe.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step2_closure_decision.txt`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step3_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step4_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step5_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step6_structure.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_golden.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_bench.log`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_benchmark_centralized.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_benchmark_centralized_diff.csv`
   16. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_benchmark_step3_comparison.csv`
   17. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_summary.txt`

## 27) Evidence Block: FP3.0 Baseline + Transition Split Snapshot

1. Timestamp (UTC): `20260306T153237Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9488636` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP3.0 baseline probe findings:
   1. Competitor path still has active runtime split:
      `eval_kernel_nodes_incremental` vs per-node `eval_kernel_node`.
   2. Compiled transition-mask application is active in runtime path
      (`apply_transition_mask_words` evidence frozen).
   3. Guard/competitor benchmark snapshot frozen in
      `guard_competitor_probe.csv`.
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_split_inventory.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_guard_competitor_probe.csv`
   5. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_step3_structure.log`
   6. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_step4_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_step5_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_step6_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_golden.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_bench.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_benchmark_centralized.csv`
   12. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_benchmark_centralized_diff.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_benchmark_step3_comparison.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T153237Z_summary.txt`

## 28) Evidence Block: FP3.1 Guard Transition Compilation Tightening

1. Timestamp (UTC): `20260306T154212Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9318182` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP3.1 tightening proof:
   1. Guard transition metadata is pre-resolved at `make_guard_input` and
      reused by guard evaluators (`guard_transition_idx` +
      `guard_transition_required`).
   2. Ranked transition source updates skip redundant per-source loops when a
      compiled source mask already covers source bits (`source_mask_covers_ids`).
   3. Ranked transition mask invalidation now uses per-node compiled slot
      mapping (`ranked_invalidate_slot_for_node`) instead of fixed slot `0`.
7. Brutal truth after FP3.1:
   1. Competitor runtime split remains active and is still the FP3.2 closure
      target.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_split_inventory.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_guard_competitor_probe.csv`
   5. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_step3_structure.log`
   6. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_step4_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_step5_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_step6_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_golden.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_bench.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_benchmark_centralized.csv`
   12. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_benchmark_centralized_diff.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_benchmark_step3_comparison.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T154212Z_summary.txt`

## 29) Evidence Block: FP3.2 Competitor Runtime Split Collapse

1. Timestamp (UTC): `20260306T160637Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9431818` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP3.2 collapse proof:
   1. `run_competitor_batch_survival_op` now routes through one evaluator API:
      `eval_kernel_nodes_incremental`.
   2. Per-node fallback loop in this path is deleted.
   3. Runtime differences are policy-only (`state` runtime vs fallback runtime),
      not alternate evaluator engines.
7. Brutal truth after FP3.2:
   1. Competitor API split target is removed.
   2. Step 3 closure decision is still pending FP3.3 evidence rollup.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_split_inventory.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_guard_competitor_probe.csv`
   5. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_step3_structure.log`
   6. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_step4_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_step5_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_step6_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_golden.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_bench.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_benchmark_centralized.csv`
   12. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_benchmark_centralized_diff.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_benchmark_step3_comparison.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP3_20260306T160637Z_summary.txt`

## 30) Evidence Block: FP3.3 Step 3 Closure Decision + Deletion-Proof Rollup

1. Timestamp (UTC): `20260308T154144Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9204545` (`example_16_guard_tie_simple`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. FP3.3 closure proof:
   1. Guard-heavy path is compiled-default (`LinearChainODE` enforced; guard
      transition metadata required).
   2. Competitor hot path split is removed:
      one evaluator API (`eval_kernel_nodes_incremental`) in
      `run_competitor_batch_survival_op`.
   3. Transition updates remain mask+slot-driven; closure probe found no full
      forced-state rebuild assignments inside competitor op block.
7. Closure decision:
   1. `decision = ACHIEVED` for Step 3.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_split_inventory.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_step3_closure_probe.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_step3_closure_decision.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_guard_competitor_probe.csv`
   7. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_step3_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_step4_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_step5_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_step6_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_golden.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_bench.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_benchmark_centralized.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_benchmark_centralized_diff.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_benchmark_step3_comparison.csv`
   16. `dev/scripts/scratch_outputs/full_plan3_FP3_20260308T154144Z_summary.txt`

## 31) Evidence Block: FP4.0 Ranked Baseline + Hot Loop Snapshot

1. Timestamp (UTC): `20260308T154626Z`
2. Commit SHA: `54437eb`
3. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Golden result:
   1. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
5. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. Worst ratio: `0.9383562` (`example_2_stop_mixture`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
   4. Ranked reference workload: `example_23_ranked_chain = 0.9250000`
6. FP4.0 baseline probe findings:
   1. Ranked hot path is frozen in
      `SequenceStateKey -> sequence_state_key -> for_each_sequence_node_transition -> sequence_prefix_density_resolved -> evaluate_sequence_density_kernel_idx`.
   2. Ranked hot-path keys are non-string in `SequenceStateKey`.
   3. Ranked-specific shaping remains outside hot loop in normalization and
      caller-side ranked sequence preparation.
   4. Sequence evaluator still has internal single-step vs multi-step dispatch,
      now explicitly frozen as Step 4 split inventory.
7. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
8. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_ranked_hotpath_inventory.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_ranked_keypath_grep.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_step3_structure.log`
   6. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_step4_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_step5_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_step6_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_golden.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_bench.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_benchmark_centralized.csv`
   12. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_benchmark_centralized_diff.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_benchmark_step3_comparison.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T154626Z_summary.txt`
