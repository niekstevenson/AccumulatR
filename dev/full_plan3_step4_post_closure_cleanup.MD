---
editor_options:
  markdown:
    wrap: 72
---

# Full Plan 3 Sub-Plan: Step 4 Post-Closure Cleanup

Date: 2026-03-09
Status: ACHIEVED WITH RECORDED RESIDUALS (`P4C.0` + `P4C.1` + `P4C.2` +
`P4C.3` + `P4C.4` + `P4C.5` frozen in bundles `20260309T124905Z`,
`20260309T130450Z`, `20260309T132214Z`, `20260309T132846Z`, and
`20260309T133248Z`)
Parent: `dev/full_plan3.MD`
Related Step: Step 4 follow-on cleanup after closure
Depends on:
1. `dev/full_plan3_step4_lower_layer_convergence.MD`
2. `dev/full_plan3_step4_deletion_map.MD`

## 1) Purpose

This artifact exists for one reason:

1. Step 4 is closed, but some remaining code in the Step 4 area is still
   larger, more duplicated, or more allocation-heavy than it should be.

This artifact is not a license to reopen Step 4 semantics.

It is a cleanup contract for:

1. representation simplification
2. prep simplification
3. ranked-state representation tightening

while preserving the Step 4 closure boundary:

1. observed single-step stays out of ranked sequence-state scheduling
2. ranked sequence-state machinery stays only where ranked semantics require it

## 2) Brutal Current State

1. Step 4 closure is real and stays closed under this artifact.
2. This optional cleanup artifact is now complete under its own contract.
3. Residual cleanup surface still exists, but it is explicitly recorded as
   non-blocking debt:
   1. payload ownership duplication in `TrialContributionSpec`
   2. repeated map copying inside ranked per-state and per-next-state merge
4. So yes, more cleanup is still possible.
5. No, more cleanup is not required to claim this artifact achieved.

## 3) What Counts

Work under this artifact counts only if at least one is true:

1. a whole live representation variant is deleted
2. a whole live prep/storage buffer family is deleted
3. a repeated hot-loop conversion/copy/sort/hash step is removed
4. a ptr-or-owned payload shape is replaced by one canonical shape
5. a ranked-state representation becomes materially more canonical and smaller

## 4) What Does Not Count

The following do not count:

1. renaming ranked helpers without deleting structure
2. moving ranked prep earlier in the pipeline without deleting the extra form
3. wrapping ptr/storage duality in a new adapter helper
4. reopening top-level observed-vs-ranked engine unification
5. adding another internal mode flag just to describe existing behavior

## 5) Hard Rejection Rules

Any batch under this artifact fails immediately if it does any of the
following:

1. routes observed single-step work into `SequenceStateKey`,
   `sequence_state_key`, `for_each_sequence_node_transition`, or
   `sequence_prefix_density_resolved`
2. reintroduces a new Step 4 specialization that was already deleted
3. adds new wrapper layers while leaving the old storage shape alive
4. claims consolidation without deleting at least one live variant or hot-loop
   operation

## 6) Targets

### C1. Caller-Side Ranked Staging Collapse

Status:

1. achieved in `P4C.1`

Historical surface:

1. [src/distributions.cpp:9192](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9192)
2. [src/distributions.cpp:9230](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9230)
3. [src/distributions.cpp:9344](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9344)
4. [src/distributions.cpp:9500](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9500)

Problem:

1. `cpp_loglik` previously carried a separate ranked staging representation:
   `ranked_label_ids_by_trial` and `ranked_rt_by_trial`.
2. That was prep-side bloat, not semantic ranked scheduling.

What must happen:

1. collapse trial extraction to one canonical prepared-trial assembly path
2. stop carrying a dedicated ranked staging buffer family if direct assembly is
   possible

What counts:

1. deleting `ranked_label_ids_by_trial`
2. deleting `ranked_rt_by_trial`
3. reducing the ranked/observed split in the `cpp_loglik` trial assembly loop

What does not count:

1. building the same ranked vectors through a new helper

### C2. Payload Ownership Flattening

Status:

1. partially achieved in `P4C.2`

Current surface:

1. [src/distributions.cpp:8402](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8402)
2. [src/distributions.cpp:8413](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8413)
3. [src/distributions.cpp:8417](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8417)
4. [src/distributions.cpp:8789](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8789)

Problem:

1. `P4C.2` removed the sequence-time ptr-plus-owned-storage pair.
2. `TrialContributionSpec` still carries ptr-plus-owned-storage variants for:
   1. competitor ids
   2. nonresponse coupling programs
3. This is cleaner than before, but still bloated.

What must happen:

1. move toward one canonical payload ownership model
2. prefer ranges/spans/offsets or one owned flat arena over ptr-or-storage
   duality

What counts:

1. deleting at least one ptr-or-storage pair from `TrialContributionSpec`
2. removing the runtime branch that checks which representation is active

What does not count:

1. adding a new helper that hides the ptr-or-storage branch

### C3. Ranked State Key Canonicalization

Status:

1. achieved in `P4C.3`

Current surface:

1. [src/distributions.cpp:7397](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7397)
2. [src/distributions.cpp:7460](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7460)
3. [src/distributions.cpp:7481](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7481)
4. [src/distributions.cpp:8174](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8174)

Problem:

1. before `P4C.3`, ranked collapse hashed map-backed exact-times and bounds
   state by:
   1. copying map entries to vectors
   2. sorting them
   3. hashing the sorted vectors
2. that materialization tax is now removed by storing ranked source-time state
   in ordered maps and hashing that canonical order directly.

What must happen:

1. make ranked state/keying use a more canonical representation
2. avoid rebuilding sortable vectors every collapse step if state can be stored
   canonically already

What counts:

1. deleting `ranked_hash_exact_times` vector materialization
2. deleting `ranked_hash_bounds` vector materialization
3. replacing map-backed collapse identity with a denser or already-ordered
   form

What does not count:

1. keeping the same copy-sort-hash path behind a helper with a new name

### C4. Ranked State Merge/Copy Tightening

Status:

1. partially achieved in `P4C.4`

Current surface:

1. [src/distributions.cpp:8124](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8124)
2. [src/distributions.cpp:8184](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8184)
3. [src/distributions.cpp:8271](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8271)

Problem:

1. `P4C.4` removed the copied `observed_prefix_exact` prefix-map staging.
2. ranked prefix execution still copies and merges `ExactSourceTimeMap` and
   `SourceTimeBoundsMap` per state and per emitted next-state.
3. This is not a wrong semantic layer, but it is still removable bloat.

What must happen:

1. reduce repeated full-map copying where persistent or overlay-style state is
   sufficient
2. keep the same ranked semantics

What counts:

1. deleting one full-map copy site
2. reducing ranked-state merge to overlay/update operations rather than full
   copy-then-prune

What does not count:

1. moving the copy into a helper and claiming cleanup

## 7) Ordered Batches

### P4C.0 Freeze

1. Freeze current post-Step-4 cleanup targets.
2. Confirm this artifact does not reopen Step 4 closure.

### P4C.1 Caller Prep Collapse

1. Attack `C1`.
2. Delete dedicated ranked staging buffers if direct prepared assembly is
   possible.

### P4C.2 Payload Flattening

1. Attack `C2`.
2. Remove at least one active ptr-or-storage representation branch.

### P4C.3 Ranked State Canonicalization

1. Attack `C3`.
2. Replace copy-sort-hash state identity with a more canonical representation.

### P4C.4 Ranked Merge Tightening

1. Attack `C4`.
2. Reduce repeated map copying/merge work inside ranked collapse.

### P4C.5 Closure Proof

1. Re-run the same Step 4 closure gates.
2. Produce a post-cleanup deletion and simplification proof.
3. Record what still remains and why.

## 8) Acceptance

This artifact is achieved only when all are true:

1. Step 4 closure still holds
2. no new Step 4 specialization or scheduler widening was introduced
3. at least one live target from `C1` through `C4` was materially removed
4. no wrapper-only consolidation was introduced
5. `R CMD INSTALL .` passes
6. `bash dev/scripts/check_step4_structure.sh` passes
7. `Rscript dev/scripts/check_loglik_golden.R` passes
8. centralized benchmark still satisfies:
   1. no workload ratio above `1.02`
   2. no P1 workload ratio above `1.00`

## 9) Required Evidence

1. `full_plan3_P4C_*_manifest.txt`
2. `full_plan3_P4C_*_summary.txt`
3. `full_plan3_P4C_*_status.txt`
4. `full_plan3_P4C_*_deletion_proof.txt`
5. `full_plan3_P4C_*_remaining_surface.txt`
6. explicit declaration:
   `Step 4 closure boundary preserved = yes/no`
7. explicit declaration:
   `No wrapper-only consolidation introduced = yes/no`

## 10) Immediate Next Action

1. `P4C.0` and `P4C.1` are now frozen in
   `full_plan3_P4C_20260309T124905Z_*`.
2. `P4C.2` is now frozen in `full_plan3_P4C_20260309T130450Z_*`.
3. `P4C.3` is now frozen in `full_plan3_P4C_20260309T132214Z_*`.
4. `P4C.4` is now frozen in `full_plan3_P4C_20260309T132846Z_*`.
5. `P4C.5` is now frozen in `full_plan3_P4C_20260309T133248Z_*`.
6. This artifact is now closed.
7. Keep `FP5.0` as the mainline next step in `dev/full_plan3.MD`.
8. If more cleanup is explicitly desired later, start a new artifact rather
   than silently reopening this one.
9. Do not let this closed cleanup artifact delay Step 5/6.

## 11) Evidence Block: P4C.0 Freeze + P4C.1 Caller Prep Collapse

1. Timestamp (UTC): `20260309T124905Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. froze the post-Step-4 cleanup surface under this artifact
   2. deleted the dedicated ranked caller-side staging buffers in `cpp_loglik`
   3. replaced them with staged use of canonical `TrialEvalInput`
4. Code evidence:
   1. ranked staging now initializes directly on `TrialEvalInput` at
      [src/distributions.cpp:8464](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8464)
      and
      [src/distributions.cpp:8472](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8472)
   2. `cpp_loglik` trial creation now stages ranked inputs directly on
      `trial_eval_inputs` at
      [src/distributions.cpp:9225](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9225)
      and
      [src/distributions.cpp:9299](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9299)
   3. final ranked preparation now finalizes those staged inputs in place at
      [src/distributions.cpp:9544](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:9544)
5. Deletion proof:
   1. `ranked_label_ids_by_trial` removed
   2. `ranked_rt_by_trial` removed
   3. net diff in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp) is negative for this cut
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
   1. worst ratio: `0.9772727` (`example_16_guard_tie_simple`)
   2. `example_22_shared_q`: `0.9739130`
   3. `example_21_simple_q`: `0.8518519`
   4. `example_23_ranked_chain`: `0.9708333`
   5. `example_5_timeout_guess`: `0.7860262`
   6. `example_6_dual_path`: `0.9001710`
   7. `depth3_guard_competitor`: `0.9294606`
   8. `example_2_stop_mixture`: `0.9555630`
   9. count above `1.02`: `0`
   10. P1 above `1.00`: `0`
9. Decision:
   1. `P4C.0`: achieved
   2. `P4C.1`: achieved
   3. artifact overall: still open for `C2`, `C3`, `C4`
10. Required declarations:
   1. `Step 4 closure boundary preserved = yes`
   2. `No wrapper-only consolidation introduced = yes`
11. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T124905Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T124905Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T124905Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T124905Z_deletion_proof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T124905Z_remaining_surface.txt`

## 12) Evidence Block: P4C.2 Payload Flattening

1. Timestamp (UTC): `20260309T130450Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. deleted the `TrialContributionSpec` sequence-time ptr-or-owned-storage
      split
   2. deleted the matching `TrialEvalInput` vector wrapper used only to feed
      that split
   3. kept the Step 4 lower-layer boundary unchanged
4. Code evidence:
   1. `TrialContributionSpec` now carries one canonical time payload at
      [src/distributions.cpp:8402](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8402)
      and
      [src/distributions.cpp:8412](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8412)
   2. `TrialEvalInput` now exposes only the canonical time pointer used by the
      evaluator at
      [src/distributions.cpp:8432](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8432)
      and
      [src/distributions.cpp:8442](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8442)
   3. contribution planning now binds directly to that canonical payload at
      [src/distributions.cpp:8789](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8789)
   4. multi-step ranked execution now consumes the canonical time payload
      directly at
      [src/distributions.cpp:8073](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8073)
      and
      [src/distributions.cpp:8936](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8936)
   5. Step 4 structure gate now rejects reintroduction of the deleted
      time-payload split at
      [check_step4_structure.sh:16](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_step4_structure.sh:16)
5. Deletion proof:
   1. `sequence_times_ptr` removed
   2. `sequence_times_storage` removed
   3. `sequence_times_vec` removed
   4. runtime representation check for owned-vs-borrowed time storage removed
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
   1. worst ratio: `0.945205479452066` (`example_2_stop_mixture`)
   2. `example_22_shared_q`: `0.893183038110564`
   3. `example_21_simple_q`: `0.822222222222196`
   4. `example_23_ranked_chain`: `0.933333333333347`
   5. `example_5_timeout_guess`: `0.803493449781648`
   6. `example_6_dual_path`: `0.840579710144874`
   7. `depth3_guard_competitor`: `0.918927545483562`
   8. `example_2_stop_mixture`: `0.945205479452066`
   9. count above `1.02`: `0`
   10. P1 above `1.00`: `0`
9. Decision:
   1. `P4C.2`: achieved
   2. `C2`: partially closed, with only competitor-id and nonresponse
      coupling-program payload duality still live
   3. artifact overall: still open for `C3`, `C4`, and residual `C2`
10. Required declarations:
   1. `Step 4 closure boundary preserved = yes`
   2. `No wrapper-only consolidation introduced = yes`
11. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T130450Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T130450Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T130450Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T130450Z_deletion_proof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T130450Z_remaining_surface.txt`

## 13) Evidence Block: P4C.3 Ranked State Canonicalization

1. Timestamp (UTC): `20260309T132214Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. replaced ranked `ExactSourceTimeMap` and `SourceTimeBoundsMap` storage
      with ordered maps
   2. removed copy-sort-hash materialization from the ranked collapse key path
   3. tightened the Step 4 structure gate so the old temporary-vector hashing
      path cannot quietly return
4. Code evidence:
   1. ranked source-time state is now canonical ordered storage at
      [src/distributions.cpp:47](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:47)
      and
      [src/distributions.cpp:48](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:48)
   2. `SequenceState` still owns the ranked source-time state, but now in the
      canonical ordered form, at
      [src/distributions.cpp:7394](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7394)
      and
      [src/distributions.cpp:7395](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7395)
   3. `ranked_hash_exact_times` now hashes the canonical ordered state
      directly at
      [src/distributions.cpp:7461](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7461)
   4. `ranked_hash_bounds` now hashes the canonical ordered state directly at
      [src/distributions.cpp:7473](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7473)
   5. `sequence_state_key` still builds the compact collapse key, but now from
      canonical ordered source-time state instead of copy-sort materialization,
      at
      [src/distributions.cpp:7486](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:7486)
   6. Step 4 structure gate now rejects temporary-vector hashing in the ranked
      hash helpers at
      [check_step4_structure.sh:80](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_step4_structure.sh:80)
      and
      [check_step4_structure.sh:94](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_step4_structure.sh:94)
5. Deletion proof:
   1. vector materialization removed from `ranked_hash_exact_times`
   2. vector materialization removed from `ranked_hash_bounds`
   3. `std::sort` removed from both ranked hash helpers
   4. ranked source-time state now has canonical ordered storage
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
   1. worst ratio: `0.968927497494148` (`example_2_stop_mixture`)
   2. `example_22_shared_q`: `0.933297326492035`
   3. `example_21_simple_q`: `0.844444444444445`
   4. `example_23_ranked_chain`: `0.91666666666666`
   5. `example_5_timeout_guess`: `0.793573777562123`
   6. `example_6_dual_path`: `0.901771336553929`
   7. `depth3_guard_competitor`: `0.946696457069901`
   8. `example_2_stop_mixture`: `0.968927497494148`
   9. count above `1.02`: `0`
   10. P1 above `1.00`: `0`
9. Decision:
   1. `P4C.3`: achieved
   2. `C3`: achieved
   3. artifact overall: still open for `C4` and residual `C2`
10. Required declarations:
   1. `Step 4 closure boundary preserved = yes`
   2. `No wrapper-only consolidation introduced = yes`
11. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132214Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132214Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132214Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132214Z_deletion_proof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132214Z_remaining_surface.txt`

## 14) Evidence Block: P4C.4 Ranked Merge Tightening

1. Timestamp (UTC): `20260309T132846Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. deleted copied prefix-map staging for ranked observed exact times
   2. replaced that staging with one incremental `prefix_exact` map per rank
   3. tightened the Step 4 structure gate so the deleted copied-prefix pattern
      cannot quietly return
4. Code evidence:
   1. ranked prefix exact times now use one incremental map at
      [src/distributions.cpp:8124](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8124)
   2. per-rank onset-source updates are now collected once and applied to that
      incremental map at
      [src/distributions.cpp:8136](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8136)
      and
      [src/distributions.cpp:8356](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8356)
   3. ranked merge still copies state-owned maps for `merged_exact` and
      `merged_bounds`, which is why `C4` is only partial, at
      [src/distributions.cpp:8184](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8184)
      and
      [src/distributions.cpp:8192](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8192)
   4. Step 4 structure gate now rejects the deleted copied-prefix-map staging
      at
      [check_step4_structure.sh:108](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_step4_structure.sh:108)
5. Deletion proof:
   1. `observed_prefix_exact` vector removed
   2. copied-prefix assignment removed
   3. ranked prefix exact storage now uses one incremental map
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
   1. worst ratio: `0.975609756097567` (`example_2_stop_mixture`)
   2. `example_22_shared_q`: `0.947826086956536`
   3. `example_21_simple_q`: `0.844444444444445`
   4. `example_23_ranked_chain`: `0.920833333333333`
   5. `example_5_timeout_guess`: `0.803493449781648`
   6. `example_6_dual_path`: `0.898550724637646`
   7. `depth3_guard_competitor`: `0.932014044047239`
   8. `example_2_stop_mixture`: `0.975609756097567`
   9. count above `1.02`: `0`
   10. P1 above `1.00`: `0`
9. Decision:
   1. `P4C.4`: achieved
   2. `C4`: partially achieved
   3. artifact overall: still open for residual `C2` and residual `C4`
10. Required declarations:
   1. `Step 4 closure boundary preserved = yes`
   2. `No wrapper-only consolidation introduced = yes`
11. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132846Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132846Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132846Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132846Z_deletion_proof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T132846Z_remaining_surface.txt`

## 15) Evidence Block: P4C.5 Closure Proof

1. Timestamp (UTC): `20260309T133248Z`
2. Commit SHA: `1f7f158`
3. Scope:
   1. no new runtime changes
   2. froze the final post-cleanup state of this artifact
   3. made the closure decision under the artifact's own acceptance rules
4. Closure decision:
   1. `decision = ACHIEVED_WITH_RECORDED_RESIDUALS`
   2. reason:
      1. Step 4 closure boundary is still preserved
      2. real cleanup deletions are present from `P4C.1` through `P4C.4`
      3. no wrapper-only consolidation was introduced
      4. all gate checks still pass
      5. residual `C2/C4` surface is explicitly recorded as non-blocking debt
5. Residual surface kept explicitly visible:
   1. payload duality still live:
      1. `step_competitor_ids_ptrs / step_competitor_ids_storage`
      2. `mass_coupling_programs_ptr / mass_coupling_programs`
   2. ranked merge copies still live:
      1. [src/distributions.cpp:8184](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8184)
      2. [src/distributions.cpp:8192](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8192)
      3. [src/distributions.cpp:8271](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8271)
      4. [src/distributions.cpp:8272](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:8272)
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
   1. worst ratio: `0.994195078080297` (`example_10_exclusion`)
   2. `example_22_shared_q`: `0.959222252227926`
   3. `example_21_simple_q`: `0.833333333333328`
   4. `example_23_ranked_chain`: `0.91666666666668`
   5. `example_5_timeout_guess`: `0.819451183352162`
   6. `example_6_dual_path`: `0.900171032496174`
   7. `depth3_guard_competitor`: `0.968400893712097`
   8. `example_2_stop_mixture`: `0.962245238890756`
   9. count above `1.02`: `0`
   10. P1 above `1.00`: `0`
9. Decision:
   1. artifact status: achieved with recorded residuals
   2. residual `C2`/`C4` surface: explicitly not blocking
   3. recommended next step: `FP5.0`
10. Required declarations:
   1. `Step 4 closure boundary preserved = yes`
   2. `No wrapper-only consolidation introduced = yes`
11. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T133248Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T133248Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T133248Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T133248Z_deletion_proof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T133248Z_remaining_surface.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_P4C_20260309T133248Z_closure_decision.txt`
