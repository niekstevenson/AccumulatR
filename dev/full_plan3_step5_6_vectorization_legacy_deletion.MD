# Full Plan 3 Sub-Plan: Step 5-6 Vectorization Rollout + Legacy/Fallback Deletion

Date: 2026-03-09
Status: IN PROGRESS (Step 5 achieved; Step 6 structural deletion complete;
FP6.3 audit/closure pending)
Parent: `dev/full_plan3.MD`
Batch ID range: `FP5.*`, `FP6.*`, `FP7.*`
Controlling Step 5 implementation artifact:
`dev/full_plan3_step5_kernel_level_batch_cutover.MD`

## 1) Objective

1. Step 5: expand vectorized distribution execution across guard/competitor/
   integration paths under one evaluator engine.
2. Step 6: delete legacy/specialized fallback runtime branches and close with
   one active engine.

## 2) Current Anchors (Where Work Applies)

Step 5 anchors:

1. Existing vector primitives:
   1. `uuber::eval_pdf_vec` at `src/distributions.cpp:749`,
      `src/distributions.cpp:3817`, `src/distributions.cpp:5149`,
      `src/distributions.cpp:5404`
   2. `uuber::eval_cdf_vec` at `src/distributions.cpp:755`,
      `src/distributions.cpp:3826`, `src/distributions.cpp:5185`,
      `src/distributions.cpp:5439`
2. Candidate scalar-heavy integration loops:
   1. coupling mass helpers and generic coupling eval around
      `src/distributions.cpp:4122` onward
   2. non-coupling per-time eval in
      `node_density_with_competitors_batch_internal` (`6238` onward)

Step 6 anchors:

1. Final public entry and unified engine surface to audit for closure:
   1. `src/distributions.cpp:6559`, `src/distributions.cpp:6602`
2. Step 6 structural guards:
   1. `dev/scripts/check_step6_structure.sh`
3. Current targeted provider proof surface:
   1. `dev/scripts/bench_examples_models_likelihood.R`

## 3) Scope Boundaries

In scope:

1. Vectorization inside unified engine only.
2. Legacy fallback/specialized helper deletion.
3. Performance stabilization on consolidated architecture.

Out of scope:

1. Reintroducing mode/workload-specific alternate engines.
2. Distribution-specific special-case evaluator stacks.

## 4) Execution Batches

### FP5.0 Vectorization Baseline + Probe

1. Freeze where vector primitives are already used.
2. Identify scalar-heavy loops still inside integration hot paths.

### FP5.1 Expand Vectorized Primitive Usage

1. Move repeated scalar pdf/cdf evaluation in guard/competitor integration
   loops to batched primitive calls.
2. Reuse batched node-time vectors across outputs/components where safe.
3. Keep distribution dispatch centralized in primitive layer.

Brutal rule for this batch:

1. Direct generic `NodeRef` target cases must be measured explicitly.
2. Example/nested suite green status does not count if
   `generic_noderef_batch_calls = 0` on those suites.
3. Any implementation that regresses the direct generic `NodeRef` target case
   is rejected even if centralized benchmark remains green.

### FP5.2 Hot-Loop Bloat Removal (Same Engine)

1. Hoist invariant checks and parameter decoding out of inner loops.
2. Reuse SoA scratch buffers and batched temporary arrays.
3. Keep one evaluator stack (no split fast engine).

### FP5.3 Verification + Evidence

1. Run full verification loop.
2. Capture vectorization coverage and runtime ratio effects.

### FP6.0 Legacy Branch Inventory Freeze

1. Freeze explicit list of specialized/fallback runtime branches.
2. Freeze call graph proving reachability before deletion.

### FP6.1 Remove Specialized Coupling Dispatch

1. Remove active dispatch from unified coupling entry to pair/nway helper
   engines.
2. Ensure generic/unified evaluator path handles those shapes without alternate
   engine branches.

### FP6.2 Remove Residual Fallback Branches

1. Remove provider-specific alternate integration behavior.
2. Remove fallback helper scaffolding that is no longer needed.
3. Keep only one active runtime evaluator stack.

### FP6.3 Deletion Audit + Verification

1. Produce explicit deleted symbol/function list.
2. Produce residual split-point list (must be empty for completion claim).
3. Run full verification loop.

### FP7.0 Final Closure (3-Run Stability)

1. Execute final gates on same commit lineage.
2. Run 3 fixed-setting benchmark passes:
   1. `ACCUMULATR_BENCH_N_REP=6`
   2. `ACCUMULATR_BENCH_TARGET_SEC=0.10`
   3. `ACCUMULATR_BENCH_MIN_INNER_REPS=40`
3. Close only if all thresholds pass.

## 5) Acceptance (What Must Be True)

Step 5 acceptance:

1. Scalar repeated pdf/cdf calls are materially reduced in RK/integration hot
   loops.
2. Distribution-specific branching is confined to primitive dispatch layer.
3. Vectorized evaluation is used beyond narrow linear-chain-only branch.

Step 6 acceptance:

1. No active runtime dispatch to pair/nway specialized engines.
2. No active fallback branch for provider-specific alternate integration engine.
3. One active runtime evaluator engine remains.
4. Deletion list is explicit and reviewable.

Final closure acceptance:

1. Golden passes (`max_abs_diff <= 1e-4`).
2. Performance passes:
   1. no workload > `1.02`
   2. no P1 > `1.00`
3. Stability 3/3 passes with fixed settings.

## 6) What Counts

1. Real vectorization and bloat removal inside one engine.
2. Real deletion of specialized/fallback branches.
3. Full gate pass on one commit lineage with complete artifacts.

## 7) What Does Not Count

1. Keeping old helper engine and calling it behind a generic wrapper.
2. Adding a second engine called "fast path" for selected workloads.
3. Declaring deletion complete without symbol list and residual split report.

## 8) Hard Stops

1. Any re-splitting introduced for speed.
2. Any correctness failure.
3. Any missing final stability evidence.

## 9) Required Evidence Bundle (`FP5.*` / `FP6.*` / `FP7.*`)

1. `full_plan3_FP5_*_vector_coverage.txt`
2. `full_plan3_FP5_*_hotloop_probe.csv`
3. `full_plan3_FP6_*_deleted_symbols.txt`
4. `full_plan3_FP6_*_residual_splitpoints.txt`
5. `full_plan3_FP7_*_stability_run1.log`
6. `full_plan3_FP7_*_stability_run2.log`
7. `full_plan3_FP7_*_stability_run3.log`
8. `full_plan3_FP[567]_*_step3_structure.log`
9. `full_plan3_FP[567]_*_step4_structure.log`
10. `full_plan3_FP[567]_*_step5_structure.log`
11. `full_plan3_FP[567]_*_step6_structure.log`
12. `full_plan3_FP[567]_*_golden.log`
13. `full_plan3_FP[567]_*_bench.log`
14. `full_plan3_FP[567]_*_diff.csv`
15. `full_plan3_FP7_*_closure_summary.txt`

## 10) Exit Rule

Step 5, Step 6, and final closure are complete only when Sections 5 and 9 are
satisfied on one commit lineage.

## 11) FP5.0 Evidence Freeze (`20260309T135122Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `no`
3. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
4. Correctness:
   1. `max_abs_diff = 1.23378e-05`
   2. `tolerance = 1e-4`
5. Performance:
   1. worst ratio: `0.975609756097567` (`example_2_stop_mixture`)
   2. workloads above `1.02`: `0`
   3. P1 above `1.00`: `0`
6. Probe findings that actually count:
   1. `shared_gate_nway_shared_triggers` samples the unified coupling stack
      directly:
      `evaluate_outcome_coupling_unified -> integrate_coupling_mass_batch_adaptive -> evaluate_labelref_batch/evaluate_coupling_acc_batch -> eval_*_vec`
   2. `example_2_stop_mixture` still spends sampled hot work in the direct
      competitor/kernel path:
      `node_density_with_competitors_from_state_dense -> eval_node_recursive_dense -> eval_kernel_node_incremental`
   3. Specialized pair/nway helper engines were not sampled in the chosen
      coupling probe, but they remain structurally reachable and therefore
      still count as live Step 6 debt.
7. Brutal status:
   1. Step 5 is not achieved.
   2. Step 6 has not started in any honest sense.
   3. `FP5.1` must attack repeated per-time guard/competitor work inside the
      same engine, not add a new fast path.
8. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_vector_coverage.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_hotloop_probe.csv`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_summary.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T135122Z_manifest.txt`

## 12) FP5.1 First Attempt (`20260309T142857Z`)

1. Commit: `1f7f158`
2. Runtime status:
   1. attempted
   2. rejected
   3. reverted
3. What was attempted:
   1. precompute event/guard values across the time batch inside
      `node_density_with_competitors_batch_internal`
   2. reuse the existing kernel composition path for generic `NodeRef`
      provider workloads
4. Why this does not count:
   1. current example and nested suites do not exercise
      `generic_noderef_batch_calls`
   2. direct guarded generic `NodeRef` target case
      (`pair_guard_generic_noderef`) regressed badly
5. Direct target-case evidence:
   1. attempted current:
      `generic_noderef_batch_calls = 1`,
      `generic_scalar_fallback_calls = 0`,
      `median_per_eval_sec = 0.0018625`
   2. baseline `1f7f158` worktree:
      `generic_noderef_batch_calls = 1`,
      `generic_scalar_fallback_calls = 0`,
      `median_per_eval_sec = 0.0009675`
   3. regression ratio:
      `1.925064599483204`
6. Revert proof:
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
   2. the first implementation is rejected as non-keepable
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_summary.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_manifest.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_attempt_pair_guard_current_regressed.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_attempt_pair_guard_baseline_1f7f158.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5_20260309T142857Z_pair_guard_current_reverted.txt`
10. Controlling next-step interpretation:
   1. `FP5.1` is redefined by
      `dev/full_plan3_step5_kernel_level_batch_cutover.MD`
   2. helper-local batch caching in `src/distributions.cpp` does not count as
      Step 5 progress

## 13) FP5K.0 Target Freeze + Harness Lock (`20260309T152055Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `no`
3. What landed:
   1. repo-owned target harness:
      `dev/scripts/bench_step5_kernel_target_cases.R`
   2. current direct target cases are now frozen in-repo instead of `/tmp`
4. Gates:
   1. `install`: pass
   2. `step5_structure`: pass
   3. `golden`: pass
   4. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `1.007948588`
   2. P1 above `1.00`: `0`
7. Target harness results that now count:
   1. `pair_guard_generic_noderef`
      1. `generic_noderef_batch_per_eval = 1`
      2. `generic_scalar_fallback_per_eval = 0`
      3. `competitor_count = 0`
   2. `nested_guard_generic_noderef`
      1. `generic_noderef_batch_per_eval = 1`
      2. `generic_scalar_fallback_per_eval = 0`
      3. `competitor_count = 0`
   3. `competitor_guard_generic_noderef`
      1. `generic_noderef_batch_per_eval = 1`
      2. `generic_scalar_fallback_per_eval = 0`
      3. `competitor_count = 2`
8. Closure decision:
   1. `FP5K.0` is achieved
   2. Step 5 is still not achieved
   3. `FP5K.1` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_split_inventory.txt`

## 14) FP5K.1 Kernel-Level Batch Executor Contract (`20260309T152855Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. `KernelBatchRuntimeState`
   2. `KernelEventBatchEvalFn`
   3. `KernelGuardBatchEvalFn`
   4. `eval_kernel_node_batch_incremental`
   5. `eval_kernel_nodes_batch_incremental`
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.979066023`
   2. P1 above `1.00`: `0`
7. Direct target case status:
   1. all three frozen target cases still show
      `generic_noderef_batch_per_eval = 1`
   2. all three still show
      `generic_scalar_fallback_per_eval = 0`
8. Closure decision:
   1. `FP5K.1` is achieved
   2. this does not close Step 5
   3. runtime is not yet ported to the new executor contract
   4. `FP5K.2` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_split_inventory.txt`

## 15) FP5K.2 Direct Generic `NodeRef` Runtime Port (`20260309T153954Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. direct generic `NodeRef` no-competitor runtime now calls the kernel
      executor contract instead of the old scalar helper path
   2. the direct unified entrypoint is now structurally guarded against
      reintroducing that old helper call
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `1`
   2. P1 above `1.00`: `0`
7. Direct target case status:
   1. `pair_guard_generic_noderef`: ratio `1`
   2. `nested_guard_generic_noderef`: ratio `1`
   3. `competitor_guard_generic_noderef`: ratio `0.9852399`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
8. Closure decision:
   1. `FP5K.2` is achieved
   2. this does not close Step 5
   3. `FP5K.0` / `FP5K.1` direct target probability values for the two
      no-competitor surfaces were wrong because the old direct branch treated
      `generic.node_id` as a dense node index
   4. competitor-bearing `NodeRef` orchestration is still not ported to the
      same executor boundary
   5. `FP5K.3` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_direct_port_probe.log`

## 16) FP5K.3 Competitor-Bearing `NodeRef` Batch Cutover (`20260309T160625Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. competitor-bearing generic `NodeRef` batch execution now uses one
      batched base-node density call plus one batched compiled-op loop
   2. `eval_kernel_nodes_batch_incremental` is now exercised on the
      competitor-bearing target surface
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.963576159`
   2. P1 above `1.00`: `0`
7. Direct target case status:
   1. `pair_guard_generic_noderef`: ratio vs `FP5K.2` = `1`
   2. `nested_guard_generic_noderef`: ratio vs `FP5K.2` = `1`
   3. `competitor_guard_generic_noderef`: ratio vs `FP5K.2` = `0.9912609`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
8. Closure decision:
   1. `FP5K.3` is achieved
   2. this does not close Step 5
   3. direct and competitor-bearing generic `NodeRef` surfaces now share the
      same batch executor boundary
   4. helper-layer orchestration still remains and is now the `FP5K.4` target
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_deleted_symbols.txt`

## 17) FP5K.4 Helper-Orchestration Deletion Cut (`20260309T162200Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. deleted the helper-local compiled competitor-op batch engine from
      `src/distributions.cpp`
   2. moved that compiled-op batch loop to
      `src/kernel_executor.cpp`
   3. left `node_density_with_competitors_batch_internal` as thin state/callback
      wiring rather than a second engine
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.980113636363642`
   2. P1 above `1.00`: `0`
7. Direct target case status:
   1. `pair_guard_generic_noderef`: ratio vs `FP5K.3` = `1`
   2. `nested_guard_generic_noderef`: ratio vs `FP5K.3` = `1`
   3. `competitor_guard_generic_noderef`: ratio vs `FP5K.3` = `1.0025189`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
8. Closure decision:
   1. `FP5K.4` is achieved
   2. this does not close Step 5
   3. dead helper-local batch orchestration is deleted
   4. `FP5K.5` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_deleted_symbols.txt`

## 18) FP5K.5 Executor-Only Stabilization (`20260309T165456Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. executor-side compiled competitor-op batching now reads runtime slots
      directly instead of materializing copied batch-value vectors
   2. no helper-layer engine or wrapper was added
   3. no runtime split was added
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.977272727272733`
   2. P1 above `1.00`: `0`
7. Direct target case status:
   1. frozen run:
      1. `pair_guard_generic_noderef`: ratio vs `FP5K.4` = `1`
      2. `nested_guard_generic_noderef`: ratio vs `FP5K.4` = `1.1111111`
      3. `competitor_guard_generic_noderef`: ratio vs `FP5K.4` = `1.0018832`
   2. supplemental rerun:
      1. `pair_guard_generic_noderef`: ratio vs `FP5K.4` = `1`
      2. `nested_guard_generic_noderef`: ratio vs `FP5K.4` = `1`
      3. `competitor_guard_generic_noderef`: ratio vs `FP5K.4` = `0.9987437`
   3. all remain `generic_noderef_batch_per_eval = 1`
   4. all remain `generic_scalar_fallback_per_eval = 0`
8. Closure decision:
   1. `FP5K.5` is achieved
   2. this does not close Step 5
   3. the executor remains on one boundary with no re-split
   4. `FP5K.6` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_target_cases_rerun.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_vector_coverage.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_split_inventory.txt`

## 19) FP5K.6 Step 5 Closure (`20260309T171110Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `no`
3. What landed:
   1. formal Step 5 closure bundle on the stabilized executor boundary
   2. current hot-loop probe for direct target surfaces
   3. current deletion proof and closure decision
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.97727272727275`
   2. P1 above `1.00`: `0`
7. Direct target case status:
   1. all three remain `generic_noderef_batch_per_eval = 1`
   2. all three remain `generic_scalar_fallback_per_eval = 0`
8. Closure decision:
   1. `FP5K.6` is achieved
   2. Step 5 is achieved
   3. Step 6 is still not achieved
   4. `FP6.0` is the next strict batch
9. Evidence bundle:
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

## 20) FP6.0 Legacy Branch Inventory Freeze (`20260309T171751Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `no`
3. What landed:
   1. full gate freeze on the first post-Step-5 baseline
   2. explicit residual Step 6 split inventory
   3. explicit call graph from public outcome probability entry to the live
      specialized/helper surface
   4. explicit deleted-symbol list proving this batch did not pretend to be a
      deletion batch
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.988636363636364`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
7. Brutal inventory result:
   1. Step 5 remains closed
   2. Step 6 is now started, but active split surface is still non-empty
   3. live Step 6 debt explicitly frozen:
      1. unified coupling still dispatches to specialized pair/nway helpers
      2. specialized helper engines are still present and reachable
      3. provider-kind runtime branching and direct fast paths are still live
   4. this batch does not count as deletion
   5. `FP6.1` is the next strict batch
8. Evidence bundle:
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

## 21) FP6.1 Remove Specialized Coupling Dispatch (`20260310T085328Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. deleted `evaluate_coupling_mass_nway`
   2. deleted `evaluate_coupling_mass_pair`
   3. inlined pair/nway execution into `evaluate_outcome_coupling_unified(...)`
   4. tightened `check_step6_structure.sh` so those symbols cannot return
   5. tightened `check_step5_structure.sh` function-boundary extraction so the
      Step 5 gate remains valid after real deletion
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.971590909090903`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
7. Brutal result:
   1. `FP6.1` is achieved
   2. pair/nway helper-engine dispatch is removed
   3. pair/nway helper symbols are deleted
   4. Step 6 is still not achieved because provider-kind fast paths and
      provider execution branching remain live
   5. `FP6.2` is the next strict batch
8. Evidence bundle:
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

## 22) FP6.2 Remove Provider-Kind Alternate Execution Behavior (`20260310T090731Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. deleted provider-kind enum/plan/runtime-branch machinery
   2. replaced provider-kind branching with once-per-call provider runtime
      resolution using function pointers/data
   3. retained EventRef-backed specialization only as resolved provider data,
      not as a live runtime branch in `src/distributions.cpp`
   4. tightened `check_step6_structure.sh` to reject provider-kind alternate
      execution machinery
   5. updated `check_step5_structure.sh` so Step 5 still proves direct NodeRef
      CDF routing through the batch executor adapter under the new resolved
      provider shape
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. targeted provider probe: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.960227272727273`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
7. Targeted provider probe result:
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
8. Brutal result:
   1. `FP6.2` is achieved
   2. active Step 6 structural split surface is empty under the current
      contract
   3. Step 6 is still not claimed in this batch because `FP6.3` audit and
      closure proof are still required
   4. `FP6.3` is the next strict batch
9. Evidence bundle:
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
