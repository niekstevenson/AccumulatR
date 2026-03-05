# Phase 3 Step 6 Consolidation Contract (No Spin)

Last updated: 2026-03-04

## Brutal Status Right Now

1.  Step 6 is **not complete**.
2.  Correctness is currently acceptable (`max_abs_diff=1.2372e-05`, tolerance `1e-04`), but architecture is still not fully unified.
3.  Major regressions were fixed (`example_6_dual_path`, `example_21_simple_q`, `example_22_shared_q`, `example_3_stop_na`), but this does not prove consolidation.
4.  Core issue: we still have mode-specific branching and mode-shaped handling in evaluator internals. That is partial convergence, not a single engine.
5.  Honest summary: we improved parity and removed some regressions, but we have not yet delivered a clean, single-path architecture.

## Primary Objective And Tradeoff Policy

1.  True priority is consolidation, not cosmetic parity.
2.  During refactor phases, temporary speed regressions are acceptable.
3.  Correctness is never optional at any stage.
4.  Final state must be one real architecture (no wrappers hiding old paths), then speed is recovered on that architecture.

## Definition Of Real Consolidation (Hard Requirements)

1.  One canonical trial execution model drives observed, ranked, and nonresponse evaluation.
2.  One tree-walk kernel executes all branches; mode differences are input data, not different engines.
3.  Shared logic lives in shared functions, not duplicated by mode.
4.  Legacy/specialized evaluator paths are deleted after cutover.
5.  Hot loops avoid redundant checks that can be prevalidated once per trial/component.

## Automatic Fail (Fake Consolidation)

1.  Wrapping old paths behind a new facade while keeping both active.
2.  Adding dispatch layers inside hot loops instead of deleting legacy branches.
3.  Keeping separate observed/ranked/nonresponse traversal logic with similar-but-not-identical code.
4.  Claiming completion based on tests passing while structural duplication remains.

## Step 6 Is Split Into Two Explicit Milestones

### Step 6A: Structural Consolidation (Primary Goal)

1.  Scope:
    1.  Eliminate duplicate evaluator architecture.
    2.  Unify branch traversal into a single runtime tree walk.
    3.  Move mode differences to normalized input/spec creation.
2.  Intermediate performance rule:
    1.  Slowdown is allowed while migrating.
    2.  Benchmark results are tracked for visibility, not used to block 6A.
3.  6A exit gates (all mandatory):
    1.  Structural scripts pass:
        1.  `dev/scripts/check_step3_structure.sh`
        2.  `dev/scripts/check_step4_structure.sh`
        3.  `dev/scripts/check_step5_structure.sh`
        4.  `dev/scripts/check_step6_structure.sh`
    2.  Correctness passes:
        1.  `Rscript dev/scripts/check_loglik_golden.R`
        2.  Required: `max_abs_diff <= 1e-4`
    3.  Deletion proof:
        1.  Legacy evaluator code paths are removed, not just unused.
        2.  Hot-path call graph shows one tree-walk executor for all modes.

### Step 6B: Optimization On The Consolidated Engine (No Re-Splitting)

1.  Scope:
    1.  Recover/improve speed only after 6A unification is complete.
    2.  No reintroduction of separate mode engines for performance.
2.  6B exit gates:
    1.  `Rscript dev/scripts/bench_centralized.R`
    2.  Metric: `median_per_eval_sec` ratio in `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`
    3.  Required:
        1.  No workload above `1.02`
        2.  No P1 workload above `1.00`
        3.  P1 set: `example_3_stop_na`, `example_21_simple_q`, `example_22_shared_q`, `example_6_dual_path`, `depth3_guard_competitor`, `shared_gate_nway_shared_triggers`
    4.  Stability:
        1.  Run 3 consecutive benchmark passes with fixed settings:
            1.  `ACCUMULATR_BENCH_N_REP=6`
            2.  `ACCUMULATR_BENCH_TARGET_SEC=0.10`
            3.  `ACCUMULATR_BENCH_MIN_INNER_REPS=40`
        2.  All 3 passes must satisfy the performance rules.

## What Must Happen Next (Concrete Remaining Work)

### A. Remove Bloat And Repeated Checks

1.  Audit hot loops for repeated validation/state checks.
2.  Hoist invariant checks to preprocessing/spec construction.
3.  Precompute masks/maps once per component/trial instead of per step.

### B. Unify Similar Functions Into Shared Core Primitives

1.  Collapse near-duplicate observed/ranked/nonresponse helper logic into shared evaluators.
2.  Keep only one transition accumulation routine shape.
3.  Keep only one outcome probability/density accumulation skeleton with pluggable leaf kernels.

### C. Enforce One Tree-Walk Architecture

1.  Build one canonical execution plan format consumed by all modes.
2.  Ensure all branches use the same traversal order and state update mechanism.
3.  Restrict mode-specific code to event normalization and output reduction.

### D. Component-Major / Vector Evaluation Pass

1.  Prototype component-first execution batches where trials share structural fingerprints.
2.  Evaluate vectors of trials per component step (SoA scratch buffers), then reduce per trial.
3.  Keep fallback to scalar execution only as a temporary migration aid; remove it before declaring full consolidation.

## Speed-Up Opportunities To Include (Post-Unification First)

1.  Structural fingerprint caching:
    1.  Group trials with identical active subtree signatures and reuse traversal products.
2.  SoA scratch and memory reuse:
    1.  Replace per-trial temporary allocations with reusable contiguous buffers.
3.  Branch reduction:
    1.  Convert frequent runtime conditionals in hot loops into precomputed flags/masks.
4.  Precomputed node/outcome resolution:
    1.  Remove repeated label-to-node and map lookups from evaluator inner loops.
5.  Batched integration windows:
    1.  Reuse integration grids/time windows for trials sharing timing shape.
6.  Lightweight memoization of deterministic subtree segments:
    1.  Cache immutable segment contributions keyed by component + bounds signature.
7.  Parallelism only after consolidation:
    1.  Parallelize over component batches/fingerprint groups, not over divergent per-mode paths.

## Completion Rules (No Early Victory)

1.  Step 6A is achieved only when structural + correctness + deletion proof pass on one commit.
2.  Step 6B is achieved only when performance + stability gates pass on that consolidated commit lineage.
3.  Phase 3 is complete only when both 6A and 6B are complete.
4.  Any failure keeps status as **Blocked**. No exceptions, no optimistic wording.

## Reporting Format (Required Every Significant Change)

1.  Commit SHA
2.  UTC timestamp
3.  Step status: `6A` or `6B`
4.  Golden result (`max_abs_diff`)
5.  Worst benchmark ratio
6.  Top 5 slowest workloads and their ratios
7.  Explicit list of deleted legacy paths

## Current Status Board

1.  Step 6A (Structural Consolidation): **Blocked**
2.  Step 6B (Optimization On Unified Engine): **Not started by definition** until 6A is complete
3.  Last known correctness status: passing
4.  Last known architecture status: partial consolidation with remaining mode-specific internals

## Active 6A Execution Artifacts

1.  Checklist: `dev/step6a_execution_checklist.md`
2.  Audit + deletion map: `dev/step6a_audit_2026-03-04.md`

## Completion Record (Fill Only When Truly Done)

1.  6A Commit SHA: `TBD`
2.  6A Achieved at UTC: `TBD`
3.  6A Golden max_abs_diff: `TBD`
4.  6A Deleted path summary: `TBD`
5.  6B Commit SHA: `TBD`
6.  6B Achieved at UTC: `TBD`
7.  6B Worst ratio (median_per_eval_sec): `TBD`
8.  Artifacts:
    1.  `dev/scripts/scratch_outputs/benchmark_centralized.csv`
    2.  `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`
    3.  `dev/scripts/scratch_outputs/benchmark_step3_comparison.csv`
