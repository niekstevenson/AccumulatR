# Phase 3 + 6 No-Shortcut Execution Contract
Last updated: 2026-03-03

## Brutal Status (As Of 2026-03-03)
1. Step 6 is **not done** yet under strict criteria.
2. Correctness is currently passing:
   1. `max_abs_diff=1.2372e-05` and tolerance `1e-04`.
3. Major regressions that blocked Step 6 earlier (`example_6_dual_path`, `example_21_simple_q`, `example_22_shared_q`, `example_3_stop_na`) were fixed.
4. Remaining speed drift still exists in several workloads and must be driven to gate thresholds before Step 6 is declared complete.
   1. Example: `example_21_simple_q` ratio `1.0370`.
   2. Example: `example_2_stop_mixture` ratio `1.0274`.
   3. Example: `example_10_exclusion` ratio `1.0265`.

## What Counts As Real Consolidation
1. Duplicate execution engines are deleted, not wrapped.
2. Hot-path behavior uses one internal execution model with fewer branches and less indirection.
3. Dead code is removed after cutover, not left behind as compatibility shims.
4. Performance improves or is neutral after consolidation under fixed benchmark settings.

## What Does Not Count
1. Moving logic behind helper wrappers while preserving old code paths.
2. Adding dispatch layers in hot loops without deleting the old branch.
3. “Parity only” when speed regresses and the plan goal says “while improving speed”.
4. Marking Step 6 complete with unresolved benchmark regressions.

## Non-Negotiable Invariants
1. No set-based forced-state operations in hot evaluator loops.
2. No recursive fallback execution in observed/ranked hot loops.
3. No string lookup in integrand loops.
4. One internal runtime kernel execution path for observed, ranked, and nonresponse likelihood internals.

## Hard Gates (All Must Pass In Same Commit)
1. Structural gate:
   1. `dev/scripts/check_step3_structure.sh`
   2. `dev/scripts/check_step4_structure.sh`
   3. `dev/scripts/check_step5_structure.sh`
   4. `dev/scripts/check_step6_structure.sh`
2. Correctness gate:
   1. `Rscript dev/scripts/check_loglik_golden.R`
   2. Required: `max_abs_diff <= 1e-4`.
3. Performance gate:
   1. `Rscript dev/scripts/bench_centralized.R`
   2. Gate metric: `median_per_eval_sec` ratio in `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`.
   3. Required:
      1. No workload above `1.02` ratio.
      2. No P1 workload above `1.00` ratio.
      3. P1 workloads: `example_3_stop_na`, `example_21_simple_q`, `example_22_shared_q`, `example_6_dual_path`, `depth3_guard_competitor`, `shared_gate_nway_shared_triggers`.
4. Stability gate:
   1. Run benchmark gate 3 times back-to-back with fixed settings:
      1. `ACCUMULATR_BENCH_N_REP=6`
      2. `ACCUMULATR_BENCH_TARGET_SEC=0.10`
      3. `ACCUMULATR_BENCH_MIN_INNER_REPS=40`
   2. Required: all 3 runs pass performance gate.

## Stop-Ship Rules (No Exceptions)
1. If any gate fails, Step 6 status remains **Blocked**.
2. No merge/claim of completion while any P1 workload regresses.
3. No adding fallback wrappers to “make tests pass” while preserving old engine internals.
4. No benchmark cherry-picking. Only centralized benchmark outputs are accepted as evidence.

## Step-by-Step Remaining Work
### Step A: Lock Evidence Baseline
1. Freeze baseline file and hash for comparison:
   1. `dev/scripts/scratch_outputs/benchmark_centralized_baseline.csv`
2. Done when:
   1. Baseline hash and timestamp are recorded in this file.

### Step B: Finish Real Cutover Cleanup
1. Remove remaining dead compatibility code that is no longer executed by the hot path.
2. Verify no duplicate “legacy evaluator” path remains callable in observed/ranked/nonresponse likelihood internals.
3. Done when:
   1. Structural scripts all pass.
   2. `rg` checks show no forbidden fallback patterns in hot-path internals.

### Step C: Close Performance Drift
1. Target remaining >1.02 workloads until performance gate passes.
2. Prioritize by worst ratio first each run.
3. Done when:
   1. Performance gate passes once.

### Step D: Prove Stability
1. Run the 3-run stability gate with fixed settings.
2. Done when:
   1. All 3 runs pass.

### Step E: Declare Step 6 Achieved
1. Step 6 is achieved only after A+B+C+D all pass on the same commit.
2. Required record in this file:
   1. Commit SHA
   2. UTC timestamp
   3. Golden `max_abs_diff`
   4. Worst benchmark ratio
   5. Path to centralized output artifacts

## Current Step 6 State
1. Status: **Blocked** (pending full hard-gate pass and 3-run stability proof).
2. Last benchmark run timestamp: `2026-03-02T19:01:36Z`.
3. Last known blocker class: residual speed drift above strict `1.02` gate in some examples workloads.

## Completion Record (Fill Only When Gates Pass)
1. Commit SHA: `TBD`
2. Achieved at UTC: `TBD`
3. Golden max_abs_diff: `TBD`
4. Worst ratio (median_per_eval_sec): `TBD`
5. Artifacts:
   1. `dev/scripts/scratch_outputs/benchmark_centralized.csv`
   2. `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`
   3. `dev/scripts/scratch_outputs/benchmark_step3_comparison.csv`
