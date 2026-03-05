# Step 6A Execution Checklist (Consolidation First, No Wrapper Exit)
Last updated: 2026-03-04
Phase status: IN PROGRESS

## Mission
Deliver a single internal trial-evaluation architecture. Mode-specific behavior is allowed only at input normalization and output reduction. Legacy internal evaluators must be deleted, not hidden.

## Current Gate Snapshot (2026-03-04)
1. `dev/scripts/check_step3_structure.sh .` -> PASS
2. `dev/scripts/check_step4_structure.sh .` -> PASS
3. `dev/scripts/check_step5_structure.sh .` -> PASS
4. `dev/scripts/check_step6_structure.sh .` -> PASS
5. `Rscript dev/scripts/check_loglik_golden.R` -> PASS

## 6A-A: Hot-Loop Bloat/Check Audit + Deletion Map
1. Audit artifact created:
   1. `dev/step6a_audit_2026-03-04.md`
2. Checklist:
   1. [x] Identify repeated runtime checks in trial hot path with line anchors.
   2. [x] Identify mode-specific duplicated control flow with line anchors.
   3. [x] Produce deletion/merge target map (function-level).
   4. [x] Convert repeated validity checks into prevalidated `TrialEvalInput` (caller pre-normalizes trial inputs).
   5. [x] Remove targeted per-trial ad-hoc allocations in inner path (`contributions`, forced component cache entries, `contribution_runtimes`, `kernel_runtimes`) via reusable caches/scratch.
   6. [x] Remove per-mode contribution builders in `evaluate_trial_probability_kernel_idx`; route via one planner (`build_trial_contributions_unified`).

## 6A-B: Shared Core Primitive Cutover
1. [x] Introduce one canonical contribution planning API that accepts normalized events.
2. [x] Move observed/ranked/nonresponse translation into normalization layer before evaluator call (`cpp_loglik` -> `TrialEvalInput`).
3. [ ] Keep exactly one transition accumulation skeleton (density/mass as leaf strategy only).
4. [x] Delete `*_shared` twins that differ only by trigger plumbing after unified trigger context lands (targeted helper scope complete, including density entrypoint merge).

## 6A-C: One Tree-Walk Runtime
1. [ ] Replace mode branch execution with one tree-walk executor.
2. [ ] Ensure ranked sequence handling uses the same executor contracts as observed/nonresponse.
3. [ ] Restrict mode enums to normalization boundary; no mode branch in hot evaluator loops.
4. [ ] Prove one hot call graph root for trial evaluation.

## 6A-D: Deletion Proof
1. [ ] Remove old evaluators made obsolete by unified planner/executor.
2. [ ] Remove fallback shims and dead compatibility utilities.
3. [ ] Run `rg` checks to show removed symbols are absent (not merely unused).

## 6A-E: Gate Closure
1. [ ] `dev/scripts/check_step3_structure.sh .`
2. [ ] `dev/scripts/check_step4_structure.sh .`
3. [ ] `dev/scripts/check_step5_structure.sh .`
4. [ ] `dev/scripts/check_step6_structure.sh .`
5. [ ] `Rscript dev/scripts/check_loglik_golden.R`
6. [ ] Record 6A completion block in `dev/plan_phase3.MD` with commit SHA + UTC timestamp.

## Stop Conditions
1. No claim of 6A completion while any legacy path in the deletion map still exists.
2. No wrapper-only migrations accepted.
3. Correctness failures block further optimization work.
