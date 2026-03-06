# Step 6A Execution Checklist (Consolidation First, No Wrapper Exit)
Last updated: 2026-03-04
Phase status: 6A COMPLETE ON WORKING TREE (pending commit)

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
3. [x] Keep exactly one transition accumulation skeleton (density/mass as leaf strategy only).
4. [x] Delete `*_shared` twins that differ only by trigger plumbing after unified trigger context lands (targeted helper scope complete, including density entrypoint merge).

## 6A-C: One Tree-Walk Runtime
1. [x] Replace mode branch execution with one tree-walk executor.
2. [x] Ensure ranked sequence handling uses the same executor contracts as observed/nonresponse.
3. [x] Restrict mode enums to normalization boundary; no mode enum branch in hot evaluator loops.
4. [x] Prove one hot call graph root for trial evaluation.

## 6A-D: Deletion Proof
1. [x] Remove old evaluators made obsolete by unified planner/executor (removed legacy coupling-mass dispatcher wrappers).
2. [x] Remove fallback shims and dead compatibility utilities (removed runtime mode enum/input builder shim in trial evaluator path).
3. [x] Run `rg` checks to show removed symbols are absent (not merely unused).

## 6A-E: Gate Closure
1. [x] `dev/scripts/check_step3_structure.sh .`
2. [x] `dev/scripts/check_step4_structure.sh .`
3. [x] `dev/scripts/check_step5_structure.sh .`
4. [x] `dev/scripts/check_step6_structure.sh .`
5. [x] `Rscript dev/scripts/check_loglik_golden.R`
6. [x] Record 6A completion block in `dev/plan_phase3.MD` with commit SHA + UTC timestamp.

## 6A-F: GenericNodeIntegral Unified Runtime Cutover
1. Artifact:
   1. `dev/step6a_generic_nodeintegral_cutover_2026-03-05.md`
2. Checklist:
   1. [x] Add runtime `GenericNodeIntegral` coupling op support in unified program path.
   2. [x] Route `native_outcome_probability_bits_impl_idx` competitor handling through `evaluate_outcome_coupling_unified` only.
   3. [x] Remove scalar fallback integration branch from `native_outcome_probability_bits_impl_idx`.
   4. [x] Implement adaptive batched integration policy inside unified evaluator stack (no second engine).
   5. [x] Re-run structural + golden gates and archive evidence in this checklist + plan completion record.

## Stop Conditions
1. No claim of 6A completion while any legacy path in the deletion map still exists.
2. No wrapper-only migrations accepted.
3. Correctness failures block further optimization work.

## 6A-C/6A-D Deletion Proof Block (2026-03-05T17:05:25Z)
1. Commit context:
   1. Working tree commit base: `8a7c6d1`
   2. Scope: trial-evaluator normalization boundary cleanup + dead coupling-mass dispatcher deletion.
2. Symbols removed from `src/distributions.cpp`:
   1. `TrialKernelMode`
   2. `make_trial_eval_input(...)`
   3. `OutcomeCouplingMassEvaluator`
   4. `evaluate_coupling_mass_pair_program(...)`
   5. `evaluate_coupling_mass_nway_program(...)`
   6. `resolve_outcome_coupling_mass_evaluator(...)`
   7. `evaluate_outcome_coupling_mass(...)`
3. `rg` absence proof (all expected no-match):
   1. `rg -n "TrialKernelMode|make_trial_eval_input\\(|evaluate_outcome_coupling_mass\\(|resolve_outcome_coupling_mass_evaluator\\(|evaluate_coupling_mass_pair_program\\(|evaluate_coupling_mass_nway_program\\(" src/distributions.cpp`
4. Fresh structural/correctness gates after deletion:
   1. `bash dev/scripts/check_step3_structure.sh` -> PASS (`dev/scripts/scratch_outputs/step6ad_delproof2_step3_20260305T170525Z.log`)
   2. `bash dev/scripts/check_step4_structure.sh` -> PASS (`dev/scripts/scratch_outputs/step6ad_delproof2_step4_20260305T170525Z.log`)
   3. `bash dev/scripts/check_step5_structure.sh` -> PASS (`dev/scripts/scratch_outputs/step6ad_delproof2_step5_20260305T170525Z.log`)
   4. `bash dev/scripts/check_step6_structure.sh` -> PASS (`dev/scripts/scratch_outputs/step6ad_delproof2_step6_20260305T170525Z.log`)
   5. `Rscript dev/scripts/check_loglik_golden.R` -> PASS (`max_abs_diff=1.2372e-05`; `dev/scripts/scratch_outputs/step6ad_delproof2_golden_20260305T170525Z.log`)
5. Honest status:
   1. 6A-D deletion proof for the removed legacy symbols above is now in place.
   2. Superseded by the `6A-C Call-Graph Proof Block` below; 6A-C closure evidence is now recorded.

## 6A-C Call-Graph Proof Block (2026-03-05T17:10:07Z)
1. Single-root hot call graph evidence (`src/distributions.cpp`):
   1. `cpp_loglik` -> `evaluate_trial_probability_kernel_idx`:
      1. definition: `8547`
      2. callsite count in hot path: one (`9226`)
   2. `evaluate_trial_probability_kernel_idx` -> `evaluate_trial_contribution_kernel_idx`:
      1. definition: `8519`
      2. callsite count: one (`8630`)
   3. `evaluate_trial_contribution_kernel_idx` -> `evaluate_sequence_density_kernel_idx`:
      1. definition: `8447`
      2. callsite count: one (`8528`)
   4. `evaluate_trial_contribution_kernel_idx` -> `mix_outcome_mass_idx`:
      1. definition: `7056`
      2. callsite count: one (`8541`)
2. Structural split-point closure completed in runtime normalization:
   1. `TrialKernelMode` removed.
   2. `TrialEvalInput` now carries one normalized sequence contract (`sequence_label_data`, `sequence_time_data`, `sequence_length`) for observed and ranked.
3. Gate bundle after closure:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6ac_final2_step3_20260305T170956Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6ac_final2_step4_20260305T170956Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6ac_final2_step5_20260305T170956Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6ac_final2_step6_20260305T170956Z.log`)
   5. golden: PASS (`max_abs_diff=1.2372e-05`; `dev/scripts/scratch_outputs/step6ac_final2_golden_20260305T170956Z.log`)
