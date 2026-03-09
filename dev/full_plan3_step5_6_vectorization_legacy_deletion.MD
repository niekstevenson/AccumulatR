# Full Plan 3 Sub-Plan: Step 5-6 Vectorization Rollout + Legacy/Fallback Deletion

Date: 2026-03-06
Status: NOT STARTED
Parent: `dev/full_plan3.MD`
Batch ID range: `FP5.*`, `FP6.*`, `FP7.*`

## 1) Objective

1. Step 5: expand vectorized distribution execution across guard/competitor/
   integration paths under one evaluator engine.
2. Step 6: delete legacy/specialized fallback runtime branches and close with
   one active engine.

## 2) Current Anchors (Where Work Applies)

Step 5 anchors:

1. Existing vector primitives:
   1. `uuber::eval_pdf_vec` at `src/distributions.cpp:745`,
      `src/distributions.cpp:3657`, `src/distributions.cpp:4888`,
      `src/distributions.cpp:5143`
   2. `uuber::eval_cdf_vec` at `src/distributions.cpp:751`,
      `src/distributions.cpp:3666`, `src/distributions.cpp:4924`,
      `src/distributions.cpp:5178`
2. Candidate scalar-heavy integration loops:
   1. coupling mass helpers and generic coupling eval around
      `src/distributions.cpp:4089` onward
   2. non-coupling per-time eval in
      `node_density_with_competitors_batch_internal` (`5986` onward)

Step 6 anchors:

1. Specialized coupling dispatch still active:
   1. `evaluate_outcome_coupling_unified` dispatch calls at
      `src/distributions.cpp:4238`, `src/distributions.cpp:4243`
2. Specialized helper engines still present:
   1. `evaluate_coupling_mass_nway` at `src/distributions.cpp:4250`
   2. `evaluate_coupling_mass_pair` at `src/distributions.cpp:4323`
3. Provider-mode runtime branching in unified coupling:
   1. `src/distributions.cpp:4120`, `src/distributions.cpp:4136`,
      `src/distributions.cpp:4146`

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
