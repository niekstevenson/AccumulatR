# Full Plan 3 Sub-Plan: Step 3-4 Guard/Competitor Compilation + Ranked Cutover

Date: 2026-03-06
Status: IN PROGRESS (Step 3 COMPLETE via `FP3.3`; Step 4 `FP4.0` complete;
top-level `FP4.1` cut rejected; `FP4L.2` + `FP4L.3` keepable lower-layer cut
complete; `FP4L.4` planner convergence landed but Step 4 still not closed;
lower-layer Step 4 work continues in
`dev/full_plan3_step4_lower_layer_convergence.MD`)
Parent: `dev/full_plan3.MD`
Batch ID range: `FP3.*`, `FP4.*`

## 1) Objective

1. Step 3: make compiled dense guard/competitor transitions the runtime default
   behavior.
2. Step 4: complete ranked cutover so ranked execution uses same fused kernel
   machinery and non-string hot key path.

## 2) Current Anchors (Where Work Applies)

Step 3 anchors:

1. Guard transition plumbing:
   1. `guard_effective_survival_internal` at `src/distributions.cpp:4614`
   2. `guard_cdf_internal` at `src/distributions.cpp:5315`
2. Competitor compiled ops:
   1. `fetch_competitor_cluster_cache` around `src/distributions.cpp:5485`
   2. `run_competitor_batch_survival_op` at `src/distributions.cpp:5795`
   3. split point `eval_kernel_nodes_incremental` vs per-node fallback at
      `src/distributions.cpp:5808`, `src/distributions.cpp:5819`
3. Non-coupling competitor density batch entry:
   1. `node_density_with_competitors_batch_internal` at
      `src/distributions.cpp:5986`

Step 4 anchors:

1. Ranked hash keys (non-string path already present):
   1. `ranked_hash_bitset` at `src/distributions.cpp:7190`
   2. `ranked_hash_exact_times` at `src/distributions.cpp:7209`
   3. `ranked_hash_bounds` at `src/distributions.cpp:7230`
2. Ranked transition executor:
   1. `for_each_sequence_node_transition` at `src/distributions.cpp:7585`
3. Mode/plumbing branch points:
   1. evaluation input mode shaping around `src/distributions.cpp:8170`

## 3) Scope Boundaries

In scope:

1. Guard/competitor transition runtime semantics.
2. Ranked transition machinery and key generation in hot loops.
3. Runtime split-point elimination where specialized fallback still active.

Out of scope:

1. Model-level ranked semantics changes.
2. New user-level ranked API behavior.

## 4) Execution Batches

### FP3.0 Baseline + Transition Split Snapshot

1. Freeze current split inventory for guard/competitor execution.
2. Freeze benchmark on guard/competitor-heavy workloads.

### FP3.1 Guard Transition Compilation Tightening

1. Move runtime guard source collection and repeated mask construction into
   precompiled transition ops where possible.
2. Ensure transition ops invalidate only required dirty slots.

### FP3.2 Competitor Runtime Split Collapse

1. Remove active runtime branch between incremental-batch and per-node
   evaluation in hot competitor path.
2. Keep one runtime mechanism with internal policy only (no alternate engine).

### FP3.3 Verification + Evidence

1. Run full verification loop.
2. Record guard/competitor split-point reduction evidence.

### FP4.0 Ranked Baseline + Hot Loop Snapshot

1. Freeze ranked-specific hot-loop call graph and key-path evidence.

### FP4.1 Ranked Fused-Kernel Convergence

1. Ensure ranked node/cdf evaluations are executed through same fused kernel
   machinery used by observed path.
2. Restrict ranked-specific logic to sequence planning and transition
   scheduling.
3. Brutal update:
   1. top-level sequence-scheduler unification for single-step observed work
      failed and is now rejected as the wrong boundary
   2. keepable Step 4 work continues under
      `dev/full_plan3_step4_lower_layer_convergence.MD`

### FP4.2 Remove Residual Ranked Fallback Logic

1. Remove any ranked-special fallback branches from ranked inner loops.
2. Ensure no string key path is used in ranked hot merge loops.

### FP4.3 Verification + Evidence

1. Run full verification loop.
2. Record ranked hot-loop structure and key-path evidence.

## 5) Acceptance (What Must Be True)

Step 3 acceptance:

1. Guard-heavy workloads run compiled transition path as runtime default.
2. Competitor hot path no longer uses active incremental-vs-per-node fallback
   split.
3. Trigger transitions update dirty slots only; no full forced-state rebuild.

Step 4 acceptance:

1. Ranked path uses same fused kernel runtime machinery as observed.
2. Ranked hot keys are non-string in hot loops.
3. Residual ranked-special fallback logic is deleted or hard-unreachable.

Shared acceptance:

1. Golden passes with `max_abs_diff <= 1e-4`.
2. Performance gates remain within thresholds.

## 6) What Counts

1. Runtime split elimination with code-path proof.
2. Ranked convergence demonstrated by call graph + grep evidence.
3. Guard/competitor transition work moved from hot loop to compiled metadata.

## 7) What Does Not Count

1. Keeping both competitor engines and selecting by runtime availability.
2. Ranked wrappers that still call separate evaluator engine internals.
3. String-key elimination outside hot loop only while hot loop still serializes
   strings.

## 8) Hard Stops

1. Any new fallback branch in competitor/guard/ranked hot loops.
2. Any correctness drift beyond tolerance.
3. Any P1 regression above gate without immediate rollback or fix batch.

## 9) Required Evidence Bundle (`FP3.*` / `FP4.*`)

1. `full_plan3_FP3_*_split_inventory.txt`
2. `full_plan3_FP3_*_guard_competitor_probe.csv`
3. `full_plan3_FP4_*_ranked_hotpath_inventory.txt`
4. `full_plan3_FP4_*_ranked_keypath_grep.log`
5. `full_plan3_FP[34]_*_step3_structure.log`
6. `full_plan3_FP[34]_*_step4_structure.log`
7. `full_plan3_FP[34]_*_step5_structure.log`
8. `full_plan3_FP[34]_*_step6_structure.log`
9. `full_plan3_FP[34]_*_golden.log`
10. `full_plan3_FP[34]_*_bench.log`
11. `full_plan3_FP[34]_*_diff.csv`

## 10) Exit Rule

Step 3 and Step 4 are complete only when Sections 5 and 9 are satisfied on one
commit lineage.

## 11) Evidence Block: FP3.0 Baseline + Transition Split Snapshot

1. Timestamp (UTC): `20260306T153237Z`
2. Commit SHA: `54437eb`
3. Batch: `FP3`, subbatch `FP3.0`
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
   1. worst ratio: `0.9488636` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP3.0 baseline probe findings:
   1. Step 3 runtime split is active in competitor path:
      `eval_kernel_nodes_incremental` vs per-node `eval_kernel_node` fallback
      in `run_competitor_batch_survival_op`.
   2. Compiled transition-mask application is present in hot chain
      (`apply_transition_mask_words` path is live).
   3. Guard/competitor benchmark snapshot frozen in
      `guard_competitor_probe.csv`.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
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

## 12) Evidence Block: FP3.1 Guard Transition Compilation Tightening

1. Timestamp (UTC): `20260306T154212Z`
2. Commit SHA: `54437eb`
3. Batch: `FP3`, subbatch `FP3.1`
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
   1. worst ratio: `0.9318182` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP3.1 tightening proof:
   1. Guard transition metadata is pre-resolved at guard-input construction
      (`guard_transition_idx` + `guard_transition` pointer) and reused by
      `guard_transition_required`.
   2. Ranked source-step updates now use compiled mask coverage checks
      (`source_mask_covers_ids`) to skip redundant per-source bit loops when
      masks already cover source bits.
   3. Ranked transition mask invalidation now uses per-node compiled slot
      mapping (`ranked_invalidate_slot_for_node`) instead of hardcoded slot
      `0`.
8. Brutal truth after FP3.1:
   1. Competitor runtime split is still active (`eval_kernel_nodes_incremental`
      vs per-node `eval_kernel_node`) and remains FP3.2 target.
9. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
10. Artifacts:
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

## 13) Evidence Block: FP3.2 Competitor Runtime Split Collapse

1. Timestamp (UTC): `20260306T160637Z`
2. Commit SHA: `54437eb`
3. Batch: `FP3`, subbatch `FP3.2`
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
   1. worst ratio: `0.9431818` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP3.2 split-collapse proof:
   1. `run_competitor_batch_survival_op` now uses one evaluator API path:
      `eval_kernel_nodes_incremental`.
   2. Legacy per-node loop fallback in this path is removed.
   3. Runtime policy is now only state-runtime vs fallback-runtime selection
      inside the same evaluator engine.
8. Brutal truth after FP3.2:
   1. Step 3 split target is removed, but Step 3 closure still requires FP3.3
      explicit completion decision block.
9. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
10. Artifacts:
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

## 14) Evidence Block: FP3.3 Step 3 Closure Decision + Deletion-Proof Rollup

1. Timestamp (UTC): `20260308T154144Z`
2. Commit SHA: `54437eb`
3. Batch: `FP3`, subbatch `FP3.3`
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
   1. worst ratio: `0.9204545` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP3.3 closure proof:
   1. Guard-heavy path remains compiled-default (`LinearChainODE` enforced in
      `guard_cdf_internal`; compiled guard transition metadata required).
   2. Competitor hot path split is absent:
      `run_competitor_batch_survival_op` uses one evaluator API
      (`eval_kernel_nodes_incremental`), no per-node fallback loop in this
      path.
   3. Transition updates are mask+slot-driven; closure probe shows no full
      forced-state rebuild assignments inside competitor op block.
8. Closure decision:
   1. `decision = ACHIEVED` for Step 3.
9. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
10. Artifacts:
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

## 15) Evidence Block: FP4.0 Ranked Baseline + Hot Loop Snapshot

1. Timestamp (UTC): `20260308T154626Z`
2. Commit SHA: `54437eb`
3. Batch: `FP4`, subbatch `FP4.0`
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
   1. worst ratio: `0.9383562` (`example_2_stop_mixture`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
   4. ranked reference workload: `example_23_ranked_chain = 0.9250000`
7. FP4.0 baseline probe findings:
   1. Ranked hot path is centered on `SequenceStateKey`,
      `sequence_state_key`, `for_each_sequence_node_transition`,
      `sequence_prefix_density_resolved`, and
      `evaluate_sequence_density_kernel_idx`.
   2. Hot-path keying is non-string in `SequenceStateKey`; key-block grep has
      no string fields.
   3. Ranked-specific shaping still exists outside hot loop in
      `normalize_trial_eval_input_for_trial` and caller-side
      `ranked_sequence_*` preparation; this is the explicit FP4.1/FP4.2 target.
   4. Sequence evaluator still has internal single-step vs multi-step dispatch
      (`sequence_length == 1`, `use_guess_shortcut`), frozen as Step 4
      split inventory.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
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
