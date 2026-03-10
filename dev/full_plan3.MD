---
editor_options:
  markdown:
    wrap: 72
---

# Phase 3 + 6 Hard-Cutover Execution Contract (No Legacy Fallback)

Last updated: 2026-03-10T09:07:31Z (UTC)
Status: OPEN (RESET UNDER THIS CONTRACT)

## 1) Brutal Current State

1. This contract is not complete.
2. Prior Step 6 evidence exists, but it does not close the full hard-cutover
   criteria defined here.
3. Step 6 structural blockers have now been removed on the current commit
   lineage:
   1. no pair/nway helper-engine dispatch remains
   2. no pair/nway helper-engine symbols remain
   3. no provider-kind enum/plan/runtime branching remains in
      `src/distributions.cpp`
4. Honest summary: Step 5 is closed; Step 6 code cut is now structurally in
   place, but closure is not claimed until `FP6.3` audit and `FP7` stability
   are frozen.
5. Honest summary: this is not final closure yet.

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
5. Step 4 deletion map: `dev/full_plan3_step4_deletion_map.MD`
6. Step 4 post-closure cleanup:
   `dev/full_plan3_step4_post_closure_cleanup.MD`
7. Step 5-6: `dev/full_plan3_step5_6_vectorization_legacy_deletion.MD`
8. Step 5 kernel-level cutover:
   `dev/full_plan3_step5_kernel_level_batch_cutover.MD`

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
5. FP4 Step 4: COMPLETE (`FP4.0` baseline/probe freeze complete; top-level
   `FP4.1` cut rejected; `FP4L.1` profiling complete; `FP4L.2` + `FP4L.3`
   keepable lower-layer convergence cut complete; `FP4L.4` planner
   convergence landed; `FP4D.1` + `FP4D.2` removed the internal prep split,
   binder glue, and ranked staging buffers; `FP4D.3` removed the ranked-named
   SoA resolver and unified SoA resolution; `FP4D.4` removed raw
   sequence-length dispatch in favor of a planner-assigned execution contract;
   `FP4D.5` deleted the remaining guess specialization and `FP4L.5` / `FP4D.6`
   closure bundle `20260309T102538Z` passed with `count_ratio_gt_1.02 = 0`
   and `count_p1_ratio_gt_1.00 = 0`,
   commit `1f7f158`)
6. FP5 Step 5: ACHIEVED (`FP5.0` baseline/probe freeze
   `20260309T135122Z` complete; first `FP5.1` runtime attempt
   `20260309T142857Z` rejected and reverted because a direct guarded generic
   `NodeRef` target case regressed by `1.925x` versus `1f7f158`; standard
   suite remained green but does not exercise `generic_noderef_batch_calls`;
   `FP5K.0` target freeze + harness lock `20260309T152055Z` is complete;
   `FP5K.1` kernel-level batch executor contract `20260309T152855Z` is
   complete; `FP5K.2` direct generic `NodeRef` runtime port
   `20260309T153954Z` is complete; `FP5K.3` competitor-bearing `NodeRef`
   batch cutover `20260309T160625Z` is complete; `FP5K.4`
   helper-orchestration deletion cut `20260309T162200Z` is complete; `FP5K.5`
   executor-only stabilization `20260309T165456Z` is complete; `FP5K.6`
   closure bundle `20260309T171110Z` passed and closes Step 5;
   controlling Step 5 execution artifact is now
   `dev/full_plan3_step5_kernel_level_batch_cutover.MD`,
   commit `1f7f158`)
7. FP6 Step 6: IN PROGRESS (`FP6.0` inventory freeze complete; `FP6.1`
   deleted pair/nway helper engines and removed active helper dispatch;
   `FP6.2` removed provider-kind runtime branching; `FP6.3` audit pending,
   commit `1f7f158`)
8. FP7 final closure: NOT STARTED

Global completion claim: NO

## 14) Immediate Next Action

1. Keep the failed `FP4.1` artifacts as negative proof of the wrong
   unification boundary.
2. Keep the keepable Step 4 evidence frozen:
   `full_plan3_FP4L_20260309T091834Z_*`,
   `full_plan3_FP4D_20260309T094513Z_*`,
   `full_plan3_FP4D_20260309T095525Z_*`,
   `full_plan3_FP4D_20260309T100610Z_*`,
   `full_plan3_FP4D_20260309T101011Z_*`,
   `full_plan3_FP4D_20260309T102538Z_*`.
3. Optional follow-on cleanup artifact now exists at
   `dev/full_plan3_step4_post_closure_cleanup.MD`; `P4C.0` / `P4C.1` are
   frozen in `full_plan3_P4C_20260309T124905Z_*`, `P4C.2` is frozen in
   `full_plan3_P4C_20260309T130450Z_*`, `P4C.3` is frozen in
   `full_plan3_P4C_20260309T132214Z_*`, `P4C.4` is frozen in
   `full_plan3_P4C_20260309T132846Z_*`, `P4C.5` closure is frozen in
   `full_plan3_P4C_20260309T133248Z_*`, and the artifact is closed with
   recorded non-blocking residuals.
4. Keep the frozen `FP5.0` evidence bundle:
   `full_plan3_FP5_20260309T135122Z_*`.
5. Keep the rejected `FP5.1` evidence bundle:
   `full_plan3_FP5_20260309T142857Z_*`.
6. Step 5 execution is now governed by
   `dev/full_plan3_step5_kernel_level_batch_cutover.MD`.
7. Keep the frozen `FP5K.0` evidence bundle:
   `full_plan3_FP5K_20260309T152055Z_*`.
8. Keep the frozen `FP5K.1` evidence bundle:
   `full_plan3_FP5K_20260309T152855Z_*`.
9. Keep the frozen `FP5K.2` evidence bundle:
   `full_plan3_FP5K_20260309T153954Z_*`.
10. Keep the frozen `FP5K.3` evidence bundle:
   `full_plan3_FP5K_20260309T160625Z_*`.
11. Keep the frozen `FP5K.4` evidence bundle:
   `full_plan3_FP5K_20260309T162200Z_*`.
12. Keep the frozen `FP5K.5` evidence bundle:
   `full_plan3_FP5K_20260309T165456Z_*`.
13. Keep the frozen `FP5K.6` evidence bundle:
   `full_plan3_FP5K_20260309T171110Z_*`.
14. Keep the frozen `FP6.0` evidence bundle:
   `full_plan3_FP6_20260309T171751Z_*`.
15. Keep the frozen `FP6.1` evidence bundle:
   `full_plan3_FP6_20260310T085328Z_*`.
16. Keep the frozen `FP6.2` evidence bundle:
   `full_plan3_FP6_20260310T090731Z_*`.
17. Next strict action is `FP6.3`:
   freeze the final Step 6 deletion audit, residual split report, and closure
   decision on this commit lineage before moving to `FP7.0`.

## 41) Evidence Block: FP6.0 Legacy Branch Inventory Freeze

1. Timestamp (UTC): `20260309T171751Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. no
4. What landed:
   1. first post-Step-5 full gate freeze for Step 6
   2. explicit residual Step 6 split inventory
   3. explicit call graph from public outcome-probability entry to the live
      specialized/helper surface
   4. explicit deleted-symbol list proving this batch did not claim deletion
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.988636363636364`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
8. Brutal inventory result:
   1. Step 5 remains closed
   2. Step 6 is now started, but active split surface is still non-empty
   3. live Step 6 debt explicitly frozen:
      1. unified coupling still dispatches to specialized pair/nway helpers
      2. specialized pair/nway helper engines are still present and reachable
      3. provider-kind runtime branching and direct fast paths are still live
   4. this does not count as deletion
   5. `FP6.1` is the next strict batch
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_summary.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_deleted_symbols.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_residual_splitpoints.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_callgraph_inventory.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_step3_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_step4_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_step5_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_step6_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_golden.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_bench.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP6_20260309T171751Z_diff.csv`

## 42) Evidence Block: FP6.1 Remove Specialized Coupling Dispatch

1. Timestamp (UTC): `20260310T085328Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. yes
4. What landed:
   1. deleted `evaluate_coupling_mass_nway`
   2. deleted `evaluate_coupling_mass_pair`
   3. pair/nway program kinds now execute inline inside
      `evaluate_outcome_coupling_unified(...)`
   4. `check_step6_structure.sh` now rejects reintroduction of the deleted
      helper symbols
   5. `check_step5_structure.sh` now uses a real function-boundary extraction
      for `evaluate_outcome_coupling_unified(...)` so Step 5 gating remains
      valid after real deletion
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.971590909090903`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
8. Brutal result:
   1. `FP6.1` is achieved
   2. pair/nway helper-engine dispatch is gone
   3. pair/nway helper symbols are deleted
   4. Step 6 is still not achieved because provider-kind fast paths and
      provider execution branching remain live
   5. `FP6.2` is the next strict batch
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_summary.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_deleted_symbols.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_residual_splitpoints.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_callgraph_inventory.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_step3_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_step4_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_step5_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_step6_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_golden.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_bench.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T085328Z_diff.csv`

## 43) Evidence Block: FP6.2 Remove Provider-Kind Alternate Execution Behavior

1. Timestamp (UTC): `20260310T090731Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. yes
4. What landed:
   1. deleted provider-kind enum/plan/runtime-branch machinery
   2. replaced provider-kind branching with once-per-call provider runtime
      resolution using function pointers/data
   3. retained EventRef-backed specialization only as resolved provider data,
      not as live runtime branching in `src/distributions.cpp`
   4. tightened `check_step6_structure.sh` to reject provider-kind alternate
      execution machinery
   5. updated `check_step5_structure.sh` so Step 5 still proves direct NodeRef
      CDF routing through the batch executor adapter under the new resolved
      provider shape
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. targeted provider probe: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.960227272727273`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
8. Targeted provider probe result:
   1. `example_21_simple_q`:
      1. `generic_labelref_batch_fastpath_calls = 1296`
      2. `generic_noderef_batch_calls = 0`
      3. `generic_scalar_fallback_calls = 0`
   2. `example_22_shared_q`:
      1. `generic_labelref_batch_fastpath_calls = 1280`
      2. `generic_noderef_batch_calls = 0`
      3. `generic_scalar_fallback_calls = 0`
   3. `example_3_stop_na`:
      1. `generic_labelref_batch_fastpath_calls = 1600`
      2. `generic_noderef_batch_calls = 0`
      3. `generic_scalar_fallback_calls = 0`
   4. interpretation:
      1. old EventRef-heavy workloads recover
      2. provider specialization remains as resolved data only
      3. runtime provider branching in `src/distributions.cpp` is gone
9. Brutal result:
   1. `FP6.2` is achieved
   2. active Step 6 structural split surface is empty under the current
      contract
   3. Step 6 is still not claimed in this batch because `FP6.3` audit and
      closure proof are still required
   4. `FP6.3` is the next strict batch
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_summary.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_deleted_symbols.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_residual_splitpoints.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_callgraph_inventory.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_targeted_unified_stats.csv`
   8. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_targeted_probe.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_step3_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_step4_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_step5_structure.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_step6_structure.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_golden.log`
   14. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_bench.log`
   15. `dev/scripts/scratch_outputs/full_plan3_FP6_20260310T090731Z_diff.csv`
15. Do not reopen Step 4 unless a later step proves this closure false under
   the same gates.

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

## 32) Evidence Block: FP5.0 Vectorization Baseline + Probe Freeze

1. Timestamp (UTC): `20260309T135122Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. `no`
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
   1. Worst ratio: `0.975609756097567` (`example_2_stop_mixture`)
   2. Count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP5.0 probe findings:
   1. The sampled coupling workload stays on the unified batch-native stack:
      `evaluate_outcome_coupling_unified ->
      integrate_coupling_mass_batch_adaptive ->
      evaluate_labelref_batch/evaluate_coupling_acc_batch -> eval_*_vec`.
   2. The sampled non-coupling workload is still dominated by the direct
      competitor/kernel path:
      `node_density_with_competitors_from_state_dense ->
      eval_node_recursive_dense -> eval_kernel_node_incremental`.
   3. Specialized pair/nway helper engines remain structurally reachable from
      unified coupling and therefore still count as Step 6 debt even though
      they were not sampled in the chosen coupling probe.
   4. Brutal status: Step 5 is not achieved; Step 6 has not started.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_vector_coverage.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_hotloop_probe.csv`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_shared_gate_nway_shared_triggers_profile.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_example_2_stop_mixture_profile.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_step3_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_step4_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_step5_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_step6_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_golden.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_bench.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_benchmark_centralized.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_diff.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_summary.txt`

## 33) Evidence Block: FP5.1 First Attempt Rejected + Reverted

1. Timestamp (UTC): `20260309T142857Z`
2. Commit SHA: `1f7f158`
3. Runtime status:
   1. attempted
   2. rejected
   3. reverted
4. Why this does not count:
   1. the standard example and nested suites do not exercise
      `generic_noderef_batch_calls`
   2. a direct guarded generic `NodeRef` target case was required and was
      explicitly measured
5. Direct target-case result (`pair_guard_generic_noderef`):
   1. attempted current:
      `generic_noderef_batch_calls = 1`,
      `generic_scalar_fallback_calls = 0`,
      `median_per_eval_sec = 0.0018625`
   2. detached worktree baseline at `1f7f158`:
      `generic_noderef_batch_calls = 1`,
      `generic_scalar_fallback_calls = 0`,
      `median_per_eval_sec = 0.0009675`
   3. regression ratio:
      `1.925064599483204`
6. Revert verification:
   1. reverted current:
      `generic_noderef_batch_calls = 1`,
      `generic_scalar_fallback_calls = 0`,
      `median_per_eval_sec = 0.0009475`
   2. reverted/current over baseline ratio:
      `0.9793281653746769`
7. Lightweight gates after revert:
   1. `step5_structure`: pass
   2. `golden`: pass
   3. `max_abs_diff = 1.23378e-05`
8. Closure decision:
   1. `FP5.1` is not achieved
   2. first implementation is rejected as non-keepable
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_summary.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_manifest.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_step5_structure.log`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_golden.log`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_attempt_pair_guard_current_regressed.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_attempt_pair_guard_baseline_1f7f158.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_pair_guard_current_reverted.txt`

## 34) Evidence Block: FP5K.0 Target Freeze + Harness Lock

1. Timestamp (UTC): `20260309T152055Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. no
4. What landed:
   1. repo-owned target harness:
      `dev/scripts/bench_step5_kernel_target_cases.R`
   2. frozen direct target surfaces:
      1. `pair_guard_generic_noderef`
      2. `nested_guard_generic_noderef`
      3. `competitor_guard_generic_noderef`
5. Gate results:
   1. `install`: pass
   2. `step5_structure`: pass
   3. `golden`: pass
   4. `bench`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `1.007948588`
   2. P1 above `1.00`: `0`
8. Target-case freeze:
   1. `pair_guard_generic_noderef`:
      `competitor_count = 0`,
      `median_per_eval_sec = 8.33333333333334e-06`,
      `generic_noderef_batch_per_eval = 1`,
      `generic_scalar_fallback_per_eval = 0`
   2. `nested_guard_generic_noderef`:
      `competitor_count = 0`,
      `median_per_eval_sec = 1.5e-05`,
      `generic_noderef_batch_per_eval = 1`,
      `generic_scalar_fallback_per_eval = 0`
   3. `competitor_guard_generic_noderef`:
      `competitor_count = 2`,
      `median_per_eval_sec = 0.001355`,
      `generic_noderef_batch_per_eval = 1`,
      `generic_scalar_fallback_per_eval = 0`
9. Closure decision:
   1. `FP5K.0` is achieved
   2. Step 5 is still not achieved
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_split_inventory.txt`

## 35) Evidence Block: FP5K.1 Kernel-Level Batch Executor Contract

1. Timestamp (UTC): `20260309T152855Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. yes
4. What landed:
   1. `KernelBatchRuntimeState`
   2. `KernelEventBatchEvalFn`
   3. `KernelGuardBatchEvalFn`
   4. `eval_kernel_node_batch_incremental`
   5. `eval_kernel_nodes_batch_incremental`
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.979066023`
   2. P1 above `1.00`: `0`
8. Direct target-case comparison versus `FP5K.0`:
   1. `pair_guard_generic_noderef`: ratio `1.0`
   2. `nested_guard_generic_noderef`: ratio `1.0`
   3. `competitor_guard_generic_noderef`: ratio `0.9876999`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
9. Honest status:
   1. `FP5K.1` is achieved
   2. this is contract consolidation, not runtime cutover
   3. direct `NodeRef` runtime is not yet ported to the new executor contract
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_split_inventory.txt`

## 36) Evidence Block: FP5K.2 Direct Generic `NodeRef` Runtime Port

1. Timestamp (UTC): `20260309T153954Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. yes
4. What landed:
   1. direct generic `NodeRef` no-competitor runtime in
      `evaluate_outcome_coupling_unified` now routes to
      `eval_node_with_forced_state_view_batch`
   2. that helper now calls `eval_kernel_node_batch_incremental`
   3. Step 5 structure gate now rejects reintroduction of the old direct
      scalar helper call at that site
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `1`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
8. Direct target-case comparison versus `FP5K.0`:
   1. `pair_guard_generic_noderef`: ratio `1`
   2. `nested_guard_generic_noderef`: ratio `1`
   3. `competitor_guard_generic_noderef`: ratio `0.9852399`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
9. Honest status:
   1. `FP5K.2` is achieved
   2. this is the first real runtime cutover onto the kernel executor contract
   3. this does not close Step 5
   4. `FP5K.0` / `FP5K.1` direct target probability values for the two
      no-competitor surfaces were wrong because the old direct branch treated
      `generic.node_id` as a dense node index
   5. `full_plan3_FP5K_20260309T153954Z_direct_port_probe.log` proves the
      corrected resolved-node scalar value matches the new batch value
   6. competitor-bearing `NodeRef` orchestration still remains in
      `node_density_with_competitors_batch_internal`
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_deleted_symbols.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_direct_port_probe.log`

## 37) Evidence Block: FP5K.3 Competitor-Bearing `NodeRef` Batch Cutover

1. Timestamp (UTC): `20260309T160625Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. yes
4. What landed:
   1. `node_density_with_competitors_batch_internal` now does batched
      base-node density plus batched compiled competitor ops
   2. `eval_kernel_nodes_batch_incremental` is now exercised on the
      competitor-bearing target surface
   3. Step 5 structure gate now rejects reintroduction of scalar base-node or
      scalar competitor-op execution inside that helper
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.963576159`
   2. worst workload: `example_10_exclusion`
   3. P1 above `1.00`: `0`
8. Direct target-case comparison versus `FP5K.2`:
   1. `pair_guard_generic_noderef`: ratio `1`
   2. `nested_guard_generic_noderef`: ratio `1`
   3. `competitor_guard_generic_noderef`: ratio `0.9912609`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
9. Honest status:
   1. `FP5K.3` is achieved
   2. direct and competitor-bearing generic `NodeRef` surfaces now share the
      same batch executor boundary
   3. this does not close Step 5
   4. helper-layer orchestration still remains in
      `node_density_with_competitors_batch_internal`
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_deleted_symbols.txt`

## 38) Evidence Block: FP5K.4 Helper-Orchestration Deletion Cut

1. Timestamp (UTC): `20260309T162200Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. yes
4. What landed:
   1. deleted `eval_nodes_batch_with_state_dense`
   2. deleted `run_competitor_batch_survival_op_batch`
   3. deleted `run_competitor_compiled_ops_product_batch`
   4. added executor-side
      `eval_kernel_competitor_product_batch_incremental`
   5. `node_density_with_competitors_batch_internal` now binds state/callbacks,
      runs one batched base-node density call, and delegates the compiled
      competitor-op batch loop to the executor
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.980113636363642`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
8. Direct target-case comparison versus `FP5K.3`:
   1. `pair_guard_generic_noderef`: ratio `1`
   2. `nested_guard_generic_noderef`: ratio `1`
   3. `competitor_guard_generic_noderef`: ratio `1.0025189`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
9. Honest status:
   1. `FP5K.4` is achieved
   2. the helper-local compiled competitor batch engine is deleted
   3. this does not close Step 5
   4. `FP5K.5` is the next strict batch
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_deleted_symbols.txt`

## 39) Evidence Block: FP5K.5 Executor-Only Stabilization

1. Timestamp (UTC): `20260309T165456Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. yes
4. What landed:
   1. executor-side compiled competitor-op batching now resolves target slots
      once and reads runtime slots directly instead of materializing copied
      batch-value vectors
   2. no helper-layer engine or wrapper was added
   3. no runtime split was added
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.977272727272733`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
8. Direct target-case comparison versus `FP5K.4`:
   1. frozen run:
      1. `pair_guard_generic_noderef`: ratio `1`
      2. `nested_guard_generic_noderef`: ratio `1.1111111`
      3. `competitor_guard_generic_noderef`: ratio `1.0018832`
   2. supplemental rerun:
      1. `pair_guard_generic_noderef`: ratio `1`
      2. `nested_guard_generic_noderef`: ratio `1`
      3. `competitor_guard_generic_noderef`: ratio `0.9987437`
   3. all remain `generic_noderef_batch_per_eval = 1`
   4. all remain `generic_scalar_fallback_per_eval = 0`
9. Honest status:
   1. `FP5K.5` is achieved
   2. the only apparent slowdown signal was timer-quantum noise on a tiny
      target surface and disappeared on the supplemental rerun
   3. this does not close Step 5
   4. `FP5K.6` is the next strict batch
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_target_cases_rerun.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_vector_coverage.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_split_inventory.txt`

## 40) Evidence Block: FP5K.6 Step 5 Closure

1. Timestamp (UTC): `20260309T171110Z`
2. Commit SHA: `1f7f158`
3. Runtime code changed:
   1. no
4. What landed:
   1. formal Step 5 closure bundle on the stabilized executor boundary
   2. current direct-target hot-loop probe
   3. current deletion proof and closure decision
5. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
6. Correctness:
   1. `max_abs_diff = 1.23378e-05`
7. Centralized benchmark:
   1. worst ratio: `0.97727272727275`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
8. Direct target-case status:
   1. closure run:
      1. `pair_guard_generic_noderef`: ratio vs `FP5K.5` = `1`
      2. `nested_guard_generic_noderef`: ratio vs `FP5K.5` = `0.9`
      3. `competitor_guard_generic_noderef`: ratio vs `FP5K.5` = `0.9931034`
   2. closure rerun:
      1. `pair_guard_generic_noderef`: ratio vs `FP5K.5` = `1`
      2. `nested_guard_generic_noderef`: ratio vs `FP5K.5` = `0.9`
      3. `competitor_guard_generic_noderef`: ratio vs `FP5K.5` = `0.9956113`
   3. all remain `generic_noderef_batch_per_eval = 1`
   4. all remain `generic_scalar_fallback_per_eval = 0`
9. Closure decision:
   1. `FP5K.6` is achieved
   2. Step 5 is achieved
   3. residual split surface is explicit Step 6 debt
   4. `FP6.0` is the next strict batch
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_target_cases_rerun.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_hotloop_probe.csv`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_vector_coverage.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_split_inventory.txt`
   8. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_deleted_symbols.txt`
   9. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_deletion_proof.txt`
   10. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_step5_closure_decision.txt`
