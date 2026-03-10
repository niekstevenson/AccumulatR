---
editor_options:
  markdown:
    wrap: 72
---

# Full Plan 3 Sub-Plan: Step 4 Deletion Map

Date: 2026-03-09
Status: COMPLETE (`FP4D.5` deletion + `FP4D.6` closure proof in bundle
`20260309T102538Z`)
Parent: `dev/full_plan3.MD`
Related Step: Step 4
Controlling execution artifact: `dev/full_plan3_step4_lower_layer_convergence.MD`

## 1) Brutal Current State

1. Step 4 is now complete under this deletion contract.
2. Lower-layer convergence improved the real architecture, and the Step 4
   area has now been cut down to:
   1. ranked-only semantics that must remain somewhere
   2. no live ranked-only internal leftovers in the `DELETE OR MERGE` bucket
   3. no live non-ranked specializations in the `JUSTIFY OR DELETE` bucket
3. `FP4D.1` removed the first block of prep leftovers:
   1. internal `ranked_mode` normalization split
   2. inline single-sequence binder glue
   3. separate ranked staging buffers inside `cpp_loglik`
4. Live Step 4 structural leftovers after `FP4D.5` are now:
   1. none
5. Strict Step 4 closure is no longer blocked:
   1. `J1` and `J2` were deleted in `FP4D.5`
   2. frozen closure bundle `20260309T102538Z` has `count_ratio_gt_1.02 = 0`
   3. frozen closure bundle `20260309T102538Z` has `count_p1_ratio_gt_1.00 = 0`
6. Until those categories were separated explicitly, Step 4 was vulnerable to
   fake progress:
   1. moving a branch without deleting it
   2. wrapping old ranked prep in a new shared helper
   3. claiming convergence while internal mode splits remain reachable

## 2) Purpose of This Artifact

This artifact is the Step 4 deletion contract.

It answers three questions with no ambiguity:

1. Which ranked-related structures are semantically required and may remain.
2. Which current ranked-specific branches are structural leftovers and must be
   deleted or merged away.
3. Which remaining specializations are tolerated only with explicit
   justification and evidence.

If a Step 4 change is not reducing the `DELETE` or `JUSTIFY/DELETE` buckets
defined here, it does not count as closure progress.

## 3) Classification Rules

### KEEP

Code stays only if all are true:

1. It represents ranked-only semantics that observed single-step likelihood
   does not have.
2. It cannot be expressed as ordinary lower-layer node/guard/competitor work
   without sequence-state loss.
3. It is not merely normalization, staging, or adapter glue.

### DELETE OR MERGE

Code must go if any are true:

1. It branches on `ranked_mode` inside internal preparation/evaluation.
2. It exists only to reshape ranked input into a second internal format.
3. It duplicates generic logic under a ranked-specific name or helper.
4. It is adapter glue needed only because canonical trial input is still not
   truly canonical.

### JUSTIFY OR DELETE

Code is tolerated temporarily only if:

1. it is not ranked-only semantics
2. it is a specialization kept for measured speed
3. the same semantics could still run on the shared engine
4. explicit evidence shows the specialization is worth keeping

If those conditions are not met, delete it.

## 4) KEEP: Ranked-Only Semantic Layer

These are allowed to remain after Step 4 closure.

### K1. Sequence-state carrier and non-string collapse key

1. [src/distributions.cpp:7380](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7380)
2. [src/distributions.cpp:7390](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7390)
3. [src/distributions.cpp:7496](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7496)

Why it stays:

1. Ranked multi-step likelihood must carry compatible prefix states.
2. Ranked state collapse is a real semantic requirement.
3. The hot key is already non-string, which is the correct direction.

What does not count:

1. Replacing this with string serialization.
2. Pushing single-step observed work through this layer.

### K2. Ranked transition-template compiler

1. [src/distributions.cpp:7539](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7539)
2. [src/distributions.cpp:7680](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7680)

Why it stays:

1. Ranked order semantics require enumerating compatible next transitions.
2. This is not a fake wrapper; it is the ranked scheduler layer.

What does not count:

1. Re-introducing recursive generic fallback underneath it.
2. Letting it own lower-layer node/guard/competitor logic that should be
   shared.

### K3. Ranked transition execution over shared lower-layer ops

1. [src/distributions.cpp:7895](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7895)

Why it stays:

1. Ranked needs a sequence-state executor.
2. The correct consolidation boundary is that this layer should call shared
   lower-layer runtime ops, not that this layer disappears.

Hard limit:

1. This layer may schedule work.
2. It may not grow its own duplicate lower-layer engine.

### K4. Ranked boundary validation semantics

Current boundary code:

1. [src/distributions.cpp:9510](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9510)

What is actually kept:

1. validation for truncation consistency
2. strict monotone time checks
3. duplicate-label rejection

Important distinction:

1. ranked validation semantics are `KEEP`
2. the old staging implementation at this site was `DELETE OR MERGE`
3. that staging implementation was removed in `FP4D.1`

## 5) DELETE OR MERGE: Structural Leftovers

These are current branches/helpers that should not survive Step 4 closure in
their present form.

### D5. Current explicit sequence-length dispatch shape

1. [src/distributions.cpp:8902](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8902)

What is nuanced here:

1. Some dispatch between "needs sequence-state" and "does not need
   sequence-state" is semantically necessary.
2. The current branch form is still a structural Step 4 target because it is
   wired directly to planner shape rather than an explicit canonical execution
   contract.

Required end state:

1. A clear data-driven execution contract such as
   `requires_sequence_state`.
2. Observed single-step stays out of sequence-state machinery.
3. Ranked multi-step keeps using the sequence-state scheduler.

Does not count:

1. Forcing everything through ranked scheduler to delete the branch.
2. Keeping the same branch but only renaming it.

## 6) JUSTIFY OR DELETE: Remaining Specializations

These were not ranked semantics. They are now closed.

### J1. Guess shortcut specialization

1. [src/distributions.cpp:8910](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8910)

Closure result:

1. Deleted in `FP4D.5`.
2. Step 4 closure bundle `20260309T102538Z` proved acceptable performance
   after deletion.

### J2. `allow_guess_shortcut` propagation in internal input

1. [src/distributions.cpp:8445](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8445)
2. [src/distributions.cpp:8493](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8493)

Closure result:

1. Deleted in `FP4D.5`.
2. No surviving internal guess-shortcut payload flag remains in
   `src/distributions.cpp`.

## 7) Step 4 Closure Criteria Under This Deletion Map

Step 4 is not achieved until all are true:

1. `KEEP` items are the only remaining ranked-only structures.
2. Every `DELETE OR MERGE` item is removed or transformed so the old form is
   unreachable.
3. Every `JUSTIFY OR DELETE` item is either deleted or defended with benchmark
   evidence.
4. Observed single-step still does not execute through sequence-state hashing
   or transition enumeration.
5. Ranked transition execution still uses shared lower-layer node/guard/
   competitor runtime ops.
6. Golden passes with `max_abs_diff <= 1e-4`.
7. Centralized benchmark gates pass:
   1. no workload ratio above `1.02`
   2. no P1 workload ratio above `1.00`

## 8) What Counts vs What Fails

### Counts

1. Deleting a listed `D*` item outright.
2. Replacing a listed `D*` item with one canonical internal representation.
3. Shrinking the reachable ranked-only surface until it matches `K*` plus any
   justified `J*`.

### Fails Automatically

1. Adding a new helper layer around a `D*` item.
2. Moving a `D*` item earlier in the pipeline and calling it "preparation."
3. Keeping a ranked-specific name for generic logic.
4. Claiming closure while any `D*` item remains live.

## 9) Ordered Execution Batches

### FP4D.0 Freeze

1. Freeze the keep/delete inventory in this artifact.
2. Use this file as the authoritative Step 4 classification.

### FP4D.1 Canonical input prep cutover

1. Delete `D1`.
2. Delete `D2`.
3. Replace them with one canonical prepared `TrialEvalInput`.
4. Status: `ACHIEVED` in evidence block `FP4D.1`.

### FP4D.2 Ranked staging collapse

1. Delete `D4`.
2. Make ranked validation emit canonical prepared trial inputs directly.
3. Status: `ACHIEVED` in evidence block `FP4D.1`.

### FP4D.3 Generic parameter-resolution naming cutover

1. Delete `D3`.
2. Ensure ranked and non-ranked callers use one generic SoA resolver.
3. Status: `ACHIEVED` in evidence block `FP4D.3`.

### FP4D.4 Execution-contract cleanup

1. Resolve `D5` into explicit data-driven execution semantics.
2. Do not widen the sequence-state boundary.
3. Status: structural target `ACHIEVED`, Step 4 closure still blocked by performance evidence.

### FP4D.5 Specialization decision

1. Delete or justify `J1`.
2. Delete or justify `J2`.

### FP4D.6 Closure proof

1. Re-run structure, golden, and benchmark gates.
2. Produce a remaining-ranked-surface inventory.
3. Prove that only `K*` and justified `J*` items remain.

## 10) Required Evidence for Closure

1. Updated ranked-surface inventory with exact `file:line` anchors.
2. Deletion list showing each `D*` item as removed or replaced.
3. Benchmark evidence for any surviving `J*` item.
4. Explicit declaration:
   `Remaining ranked-only code matches KEEP bucket = yes/no`
5. Explicit declaration:
   `No wrapper-only consolidation introduced = yes/no`

## 11) Immediate Next Action

1. Keep `full_plan3_FP4D_20260309T102538Z_*` frozen as the Step 4 closure
   bundle.
2. Do not reopen this artifact unless a later step breaks one of its closure
   criteria.
3. Non-blocking follow-on cleanup is captured in
   `dev/full_plan3_step4_post_closure_cleanup.MD`.
4. Next work belongs to Step 5.

## 12) Evidence Block: FP4D.1 Canonical Input Prep Cutover

1. Timestamp (UTC): `20260309T094513Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. deleted internal `ranked_mode` normalization helper
   2. deleted inline single-sequence binder helper
   3. deleted separate ranked staging buffers in `cpp_loglik`
   4. replaced them with in-place prepared `TrialEvalInput` builders
4. Code evidence:
   1. canonical in-place observed builder at
      [src/distributions.cpp:8473](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8473)
   2. canonical in-place ranked builder at
      [src/distributions.cpp:8495](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8495)
   3. direct prepared-input assembly in `cpp_loglik` at
      [src/distributions.cpp:9544](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9544)
5. Deletion proof:
   1. `normalize_trial_eval_input_for_trial` removed
   2. `bind_trial_eval_input_inline_storage` removed
   3. `ranked_sequence_label_ids` removed
   4. `ranked_sequence_times` removed
   5. `ranked_sequence_valid` removed
6. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
7. Golden result:
   1. `max_abs_diff = 1.23378e-05`
8. Performance snapshot (`median_per_eval_sec`, current/baseline):
   1. worst ratio: `1.0086957` (`example_22_shared_q`)
   2. `example_21_simple_q`: `0.9259259`
   3. `example_22_shared_q`: `1.0086957`
   4. `example_6_dual_path`: `0.9275362`
   5. `depth3_guard_competitor`: `0.9368018`
   6. `example_23_ranked_chain`: `0.9918699`
   7. `example_2_stop_mixture`: `0.9673601`
   8. count above `1.02`: `0`
   9. P1 above `1.00`: `1`
9. Closure decision:
   1. `FP4D.1`: achieved
   2. `FP4D.2`: achieved
   3. Step 4 overall: not achieved
10. Reason Step 4 remains open:
    1. `D3` still live
    2. `D5` still live
    3. `J1` and `J2` still undecided
    4. strict P1 gate still fails on `example_22_shared_q`
11. Artifacts:
    1. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_manifest.txt`
    2. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_status.txt`
    3. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_install.log`
    4. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_step3_structure.log`
    5. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_step4_structure.log`
    6. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_step5_structure.log`
    7. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_step6_structure.log`
    8. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_golden.log`
    9. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_bench.log`
    10. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_benchmark_centralized.csv`
    11. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_benchmark_centralized_diff.csv`
    12. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_benchmark_step3_comparison.csv`
    13. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T094513Z_summary.txt`

## 13) Removed Leftovers (Already Closed)

### D3. Ranked-named trial-parameter SoA resolver

1. Old site:
   [src/distributions.cpp:7855](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7855)
2. Removal batch:
   `FP4D.3`
3. Replacement:
   one generic context-aware SoA resolver at
   [src/distributions.cpp:1910](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:1910)

### D5. Explicit sequence-length dispatch shape

1. Old site:
   [src/distributions.cpp:8923](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8923)
2. Removal batch:
   `FP4D.4`
3. Replacement:
   planner-assigned execution contract carried by
   [src/distributions.cpp:8403](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8403)
   and consumed by
   [src/distributions.cpp:8935](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8935)

## 14) Evidence Block: FP4D.3 Generic SoA Resolver Cutover

1. Timestamp (UTC): `20260309T095525Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. deleted ranked-named trial-parameter SoA resolver
   2. introduced one generic context-aware resolver
   3. collapsed duplicate manual SoA-resolution logic onto that helper
4. Code evidence:
   1. generic resolver at
      [src/distributions.cpp:1910](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:1910)
   2. ranked transition path now uses that resolver at
      [src/distributions.cpp:7900](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7900)
   3. NodeEvalState now uses that resolver at
      [src/distributions.cpp:3018](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:3018)
5. Deletion proof:
   1. `resolve_trial_params_soa_for_ranked` removed
   2. cache-only `resolve_trial_params_soa(const TrialParamSet*)` removed
   3. duplicate manual SoA resolution removed from `NodeEvalState`
6. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
7. Golden result:
   1. `max_abs_diff = 1.23378e-05`
8. Performance snapshot (`median_per_eval_sec`, current/baseline):
   1. worst ratio: `1.0086957` (`example_22_shared_q`)
   2. `example_21_simple_q`: `0.9259259`
   3. `example_22_shared_q`: `1.0086957`
   4. `example_6_dual_path`: `0.9637681`
   5. `depth3_guard_competitor`: `0.9546760`
   6. `example_23_ranked_chain`: `0.9916667`
   7. `example_2_stop_mixture`: `0.9863014`
   8. count above `1.02`: `0`
   9. P1 above `1.00`: `1`
9. Closure decision:
   1. `FP4D.3`: achieved
   2. Step 4 overall: not achieved
10. Reason Step 4 remains open:
    1. `D5` still live
    2. `J1` and `J2` still undecided
    3. strict P1 gate still fails on `example_22_shared_q`
11. Artifacts:
    1. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_manifest.txt`
    2. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_status.txt`
    3. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_install.log`
    4. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_step3_structure.log`
    5. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_step4_structure.log`
    6. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_step5_structure.log`
    7. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_step6_structure.log`
    8. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_golden.log`
    9. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_bench.log`
    10. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_benchmark_centralized.csv`
    11. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_benchmark_centralized_diff.csv`
    12. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_benchmark_step3_comparison.csv`
    13. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T095525Z_summary.txt`

## 15) Evidence Block: FP4D.4 Execution-Contract Cleanup

1. Timestamp (UTC): `20260309T100610Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. removed raw sequence-length dispatch from `evaluate_sequence_density_kernel_idx`
   2. replaced it with planner-assigned `TrialSequenceExecutionKind`
   3. aligned runtime-slot and sequence-compiler decisions to the same execution contract
   4. updated `check_step4_structure.sh` to reject raw `sequence_length` dispatch and require execution-contract dispatch
4. Code evidence:
   1. execution contract enum and field at
      [src/distributions.cpp:8391](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8391)
      and
      [src/distributions.cpp:8403](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8403)
   2. planner assignment at
      [src/distributions.cpp:8690](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8690)
   3. evaluator switch at
      [src/distributions.cpp:8935](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8935)
   4. structure gate contract check at
      [check_step4_structure.sh](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_step4_structure.sh:30)
5. Deletion proof:
   1. raw `sequence_length == 1u` evaluator dispatch removed
   2. runtime/compiler selection no longer keys off `sequence_outcome_indices.size()`
6. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass at command level, but strict perf gate did not clear
7. Golden result:
   1. `max_abs_diff = 1.23378e-05`
8. Frozen performance evidence:
   1. bundle `20260309T100610Z`: worst `1.0260870` at `example_22_shared_q`
   2. bundle `20260309T101011Z`: worst `1.1304348` at `example_22_shared_q`
   3. non-frozen probe rerun: `example_22_shared_q = 1.0197137`
9. Closure decision:
   1. structural `D5`: achieved
   2. Step 4 overall: not achieved
10. Reason Step 4 remains open:
    1. `J1` and `J2` still undecided
    2. strict performance evidence for `example_22_shared_q` remains unstable/red in this environment
11. Artifacts:
    1. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T100610Z_manifest.txt`
    2. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T100610Z_status.txt`
    3. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T100610Z_summary.txt`
    4. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T101011Z_manifest.txt`
    5. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T101011Z_status.txt`
    6. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T101011Z_summary.txt`

## 16) Evidence Block: FP4D.5 + FP4D.6 Closure Proof

1. Timestamp (UTC): `20260309T102538Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. deleted `J1` guess shortcut specialization from the unified sequence
      evaluator
   2. deleted `J2` internal guess-shortcut payload propagation
   3. strengthened Step 4 structure gate so this specialization cannot return
      unnoticed
4. Deletion proof:
   1. `step_guess_shortcut` removed from `TrialContributionSpec`
   2. `allow_guess_shortcut` removed from `TrialEvalInput`
   3. unified sequence evaluator no longer calls
      `accumulate_component_guess_density_idx`
5. Remaining-ranked-surface proof:
   1. frozen in
      `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T102538Z_remaining_ranked_surface.txt`
   2. only `KEEP`-bucket sequence-state structures remain
6. Gate results:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
7. Golden result:
   1. `max_abs_diff = 1.23378e-05`
8. Performance snapshot (`median_per_eval_sec`, current/baseline):
   1. worst ratio: `0.9917695` (`example_23_ranked_chain`)
   2. `example_5_timeout_guess`: `0.8021996`
   3. `example_22_shared_q`: `0.9826087`
   4. count above `1.02`: `0`
   5. P1 above `1.00`: `0`
9. Closure decision:
   1. `FP4D.5`: achieved
   2. `FP4D.6`: achieved
   3. Step 4 overall: achieved
10. Required declarations:
    1. `Remaining ranked-only code matches KEEP bucket = yes`
    2. `No wrapper-only consolidation introduced = yes`
11. Artifacts:
    1. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T102538Z_manifest.txt`
    2. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T102538Z_status.txt`
    3. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T102538Z_summary.txt`
    4. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T102538Z_step4_closure_decision.txt`
    5. `dev/scripts/scratch_outputs/full_plan3_FP4D_20260309T102538Z_remaining_ranked_surface.txt`
