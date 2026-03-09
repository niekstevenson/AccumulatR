---
editor_options:
  markdown:
    wrap: 72
---

# Full Plan 3 Sub-Plan: Step 4 Lower-Layer Convergence

Date: 2026-03-08
Status: IN PROGRESS
Parent: `dev/full_plan3.MD`
Related Step: Step 4
Batch ID range: `FP4L.*`

## 1) Brutal Current State

1. The original `FP4.1` top-level cutover attempt failed.
2. Two transient evidence runs already prove that routing single-step observed
   work through ranked sequence-state machinery is not acceptable:
   1. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T160703Z_summary.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP4_20260308T161426Z_summary.txt`
3. Those runs passed correctness but failed the real goal:
   1. `depth3_guard_competitor` slowed by about `19.5x`
   2. `example_2_stop_mixture` slowed by about `17.3x`
   3. `example_16_guard_tie_simple` slowed by about `5.23x`
   4. `example_21_simple_q` slowed by about `1.22x`
   5. `example_22_shared_q` slowed by about `1.13x`
   6. `example_6_dual_path` slowed by about `1.80x`
4. Honest diagnosis:
   1. These workloads share the same model IR tree.
   2. They do not share the same sequence-state needs.
   3. The failed cutover unified the wrong layer.

## 2) Correct Convergence Boundary

The correct unification target is:

1. Shared lower-layer IR/kernel execution:
   1. event/node eval
   2. guard eval
   3. competitor survival/density ops
   4. forced-bit transitions
   5. runtime invalidation/caching
2. Shared planning/input normalization where semantics are the same.
3. Sequence-state scheduling only when sequence semantics require it.

The incorrect unification target is:

1. Forcing observed single-step work to execute through ranked
   sequence-state/scenario expansion just to claim one path.

## 3) What Counts

1. Observed single-step and ranked multi-step share the same lower-layer
   node/guard/competitor primitives.
2. Ranked-specific logic is limited to sequence-state scheduling and state
   collapse.
3. Shared runtime/caching/invalidations move down to the lower layer.
4. Deletions or hard-unreachable proof exist for duplicated lower-layer
   helpers once convergence is complete.

## 4) What Does Not Count

1. Routing single-step observed likelihood through
   `sequence_prefix_density_resolved` or
   `for_each_sequence_node_transition`.
2. Keeping separate lower-layer engines and only wrapping them in one
   top-level API.
3. Claiming convergence because both paths eventually touch IR nodes somewhere.
4. Leaving duplicated guard/competitor/node execution helpers alive and
   reachable.

## 5) Explicit Rejection Criteria

Any batch fails immediately if any of the following are true:

1. Single-step observed work enters sequence-state hashing/collapse
   (`SequenceStateKey`, `sequence_state_key`) in steady-state execution.
2. Single-step observed work enters ranked transition-template enumeration in
   steady-state execution.
3. A new mode-specific fast path is added instead of converging lower-layer
   ops.
4. Performance regression exceeds plan gates on the P1 set.

## 6) Execution Batches

### FP4L.0 Failure Freeze + Profiling Scope

1. Freeze the failed `FP4.1` evidence as negative proof.
2. Define exact workloads for hotspot diagnosis:
   1. `depth3_guard_competitor`
   2. `example_2_stop_mixture`

### FP4L.1 Transient Failed-Path Profiling

1. Reapply the failed `FP4.1` cut transiently for measurement only.
2. Capture C++ profiles for the two target workloads.
3. Measure where `for_each_sequence_node_transition`-driven execution pays
   extra cost:
   1. transition-template enumeration
   2. forced-bit copying
   3. state hashing/collapse
   4. competitor-survival calls
   5. guard eval
4. Restore the pre-profile source state after evidence capture.

### FP4L.2 Shared Lower-Layer Inventory

1. Compare current good observed single-step path with ranked path.
2. Enumerate the genuinely shared lower-layer ops.
3. Enumerate the ranked-only sequence-state ops.
4. Produce a deletion target list for duplicate lower-layer helpers.

### FP4L.3 Lower-Layer Extraction

1. Converge duplicated lower-layer node/guard/competitor execution into shared
   internal ops.
2. Keep sequence-state code outside this layer.
3. Do not route single-step observed work through sequence-state scaffolding.

### FP4L.4 Planner Convergence

1. Make observed and ranked planners emit compatible lower-layer work payloads.
2. Keep ranked state scheduling only where multiple-step ordering semantics
   exist.

### FP4L.5 Verification + Closure Decision

1. Run full gates.
2. Prove that convergence happened at the lower layer, not only at the
   top-level wrapper surface.
3. Record remaining Step 4 blockers explicitly if closure is not achieved.

## 7) Acceptance

Step 4 under this artifact is achieved only when all are true:

1. Ranked hot keys remain non-string in hot loops.
2. Observed single-step does not execute through sequence-state machinery.
3. Ranked and observed share lower-layer node/guard/competitor runtime ops.
4. No duplicated reachable lower-layer engine remains where semantics are the
   same.
5. Golden passes with `max_abs_diff <= 1e-4`.
6. Centralized benchmark gates pass:
   1. no workload ratio above `1.02`
   2. no P1 workload ratio above `1.00`

## 8) Required Evidence Bundle (`FP4L.*`)

1. `full_plan3_FP4L_*_manifest.txt`
2. `full_plan3_FP4L_*_summary.txt`
3. `full_plan3_FP4L_*_status.txt`
4. `full_plan3_FP4L_*_workload_scope.txt`
5. `full_plan3_FP4L_*_depth3_guard_competitor_profile.txt`
6. `full_plan3_FP4L_*_example_2_stop_mixture_profile.txt`
7. `full_plan3_FP4L_*_profile_hotspot_counts.csv`
8. `full_plan3_FP4L_*_split_inventory.txt`
9. `full_plan3_FP4L_*_deletion_targets.txt`
10. Explicit declaration:
    `No wrapper-only consolidation introduced = yes/no`

## 9) Immediate Next Action

1. Trim the residual `example_22_shared_q` single-step planner overhead so the
   strict P1 gate returns to `<= 1.00`.
2. Then execute `FP4L.5`.
3. Do not widen the sequence-state boundary again; observed single-step must
   stay out of `sequence_state_key` / `for_each_sequence_node_transition`.

## 10) Evidence Block: FP4L.1 Transient Failed-Path Profiling

1. Timestamp (UTC): `20260308T164541Z`
2. Commit SHA: `54437eb`
3. Batch: `FP4L`, subbatch `FP4L.1`
4. Mode:
   1. transient reapplication of the failed `FP4.1` cut for measurement only
   2. source restored after capture
5. Workloads:
   1. `depth3_guard_competitor`
   2. `example_2_stop_mixture`
6. Primary findings:
   1. Both sampled workloads spend time inside
      `sequence_prefix_density_resolved -> for_each_sequence_node_transition`.
   2. The main repeated costs are lower-layer runtime/state setup, not
      sequence-state hashing.
   3. `sequence_state_key` / `ranked_hash_*` do not appear as sampled hotspots
      in these single-step failed-path runs.
7. Hotspot evidence:
   1. `depth3_guard_competitor`:
      `eval_kernel_node_incremental=54`,
      `NodeEvalState::NodeEvalState=21`,
      `reset_kernel_runtime=23`,
      `std::vector<double>::assign=43`
   2. `example_2_stop_mixture`:
      `eval_kernel_node_incremental=133`,
      `NodeEvalState::NodeEvalState=57`,
      `reset_kernel_runtime=73`,
      `std::vector<double>::assign=138`
8. Closure implication:
   1. Failed `FP4.1` slowdown comes from pushing single-step observed work into
      transition-scenario execution that repeatedly rebuilds runtime state and
      re-enters competitor kernel eval.
   2. Next keepable cut must target lower-layer runtime reuse, not top-level
      sequence-scheduler unification.
9. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_workload_scope.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_depth3_guard_competitor_profile.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_depth3_guard_competitor_profile.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_example_2_stop_mixture_profile.txt`
   8. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_example_2_stop_mixture_profile.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_profile_hotspot_counts.csv`
   10. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_split_inventory.txt`
   11. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T164541Z_deletion_targets.txt`

## 11) Evidence Block: FP4L.2 + FP4L.3 Keepable Lower-Layer Convergence Cut

1. Timestamp (UTC): `20260308T185700Z`
2. Commit SHA: `54437eb` (working tree)
3. Batch: `FP4L`, subbatch `FP4L.2` + `FP4L.3`
4. What landed:
   1. Ranked transition execution now reuses one `NodeEvalState` per
      transition instead of rebuilding state per density/cdf step.
   2. Observed single-step sequence density now uses the same lower-layer
      state-based dense node/competitor evaluator.
   3. Runtime-slot reuse now applies to all real sequence-density
      contributions, not only multi-step ranked ones.
5. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
6. Golden result:
   1. `max_abs_diff = 1.23378e-05`
7. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. worst ratio: `0.9829545` (`example_16_guard_tie_simple`)
   2. `depth3_guard_competitor`: `0.9214810`
   3. `example_2_stop_mixture`: `0.9369187`
   4. `example_21_simple_q`: `0.8148148`
   5. `example_22_shared_q`: `0.8869565`
   6. `example_6_dual_path`: `0.8260870`
   7. `example_23_ranked_chain`: `0.9512195`
   8. count above `1.02`: `0`
   9. P1 above `1.00`: `0`
8. What this proves:
   1. Lower-layer convergence can be keepable when the shared boundary is
      `NodeEvalState` + compiled competitor ops, not the top-level
      sequence-state scheduler.
   2. Observed single-step still stays out of ranked sequence-state machinery.
   3. The landed cut is not wrapper-only consolidation.
9. What this does not prove:
   1. Step 4 is not closed.
   2. Planner convergence is still pending.
   3. Duplicate lower-layer helper deletion/isolation proof is still pending.
10. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
11. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_workload_scope.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_deletion_targets.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_step3_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_step4_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_step5_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_step6_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_golden.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_bench.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_benchmark_centralized.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_benchmark_centralized_diff.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260308T185700Z_benchmark_step3_comparison.csv`

## 12) Evidence Block: FP4L.4 Planner Convergence

1. Timestamp (UTC): `20260309T091834Z`
2. Commit SHA: `54437eb` (working tree)
3. Batch: `FP4L`, subbatch `FP4L.4`
4. What landed:
   1. `TrialContributionSpec` now carries one per-step payload shape for
      sequence-density planning across observed single-step and ranked
      sequence contributions.
   2. Ranked competitor/persistent-source shaping moved out of
      `sequence_prefix_density_resolved` and into
      `build_trial_contributions_unified`.
   3. Single-step competitor payload uses zero-copy references when the model
      already owns the vector.
5. Gate command status:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
6. Golden result:
   1. `max_abs_diff = 1.23378e-05`
7. Performance snapshot (`median_per_eval_sec` ratio current/baseline):
   1. worst ratio: `1.0086957` (`example_22_shared_q`)
   2. `depth3_guard_competitor`: `0.9297798`
   3. `example_2_stop_mixture`: `0.9707424`
   4. `example_16_guard_tie_simple`: `0.9545455`
   5. `example_21_simple_q`: `0.9259259`
   6. `example_22_shared_q`: `1.0086957`
   7. `example_6_dual_path`: `0.9217751`
   8. `example_23_ranked_chain`: `0.9708333`
   9. count above `1.02`: `0`
   10. P1 above `1.00`: `1`
8. Honest decision:
   1. `FP4L.4` is keepable progress.
   2. `FP4L.4` is not Step 4 closure because the strict P1 gate still fails on
      `example_22_shared_q`.
9. What this proves:
   1. Planner convergence can be made real without routing observed
      single-step through sequence-state machinery.
   2. The canonical planner payload is now per-step, not mode-specific
      single-step fields plus ranked runtime shaping.
10. What remains:
   1. shave the residual `example_22_shared_q` overhead
   2. run `FP4L.5` closure/deletion proof
11. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
12. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_status.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_workload_scope.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_deletion_targets.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_step3_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_step4_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_step5_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_step6_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_golden.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_bench.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_benchmark_centralized.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_benchmark_centralized_diff.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP4L_20260309T091834Z_benchmark_step3_comparison.csv`
