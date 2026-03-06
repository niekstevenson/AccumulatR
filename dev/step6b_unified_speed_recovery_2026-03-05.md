# Step 6B Execution Artifact: Unified Speed Recovery (No Re-Splitting)
Date: 2026-03-05
Status: IN PROGRESS (B3/B4 + B5 + Pair/NWay scalar fallback removal landed; latest B6 3-run stability gate attempt failed)
Owner: Codex + reviewer

## Why This Exists
1. Structural unification has been enforced in runtime outcome-mass evaluation.
2. Correctness is passing, but performance is materially regressed in core example workloads.
3. This artifact defines a strict 6B recovery path that preserves one unified architecture.

## Non-Negotiable Constraints
1. No reintroduction of scalar fallback path in `native_outcome_probability_bits_impl_idx`.
2. No parallel legacy engines hidden behind wrappers.
3. All optimizations must stay inside unified evaluator stack.
4. Correctness (`max_abs_diff <= 1e-4`) must hold after every meaningful optimization batch.

## Current Baseline (Blocking 6B, Post Pair/NWay Scalar-Fallback Removal)
Source: `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv` (`median_per_eval_sec`)
1. Worst: `example_2_stop_mixture = 1.05753424657533`
2. `example_10_exclusion = 1.02843786521231`
3. `example_6_dual_path = 1.02608695652174`
4. `example_7_mixture = 1.01993761019939`
5. `shared_gate_nway_shared_triggers = 1.00546821599453`
6. P1 set status:
   1. `example_3_stop_na = 1.00000000000000`
   2. `example_21_simple_q = 0.925925925925894`
   3. `example_22_shared_q = 1.00173913043479`
   4. `example_6_dual_path = 1.02608695652174`
   5. `depth3_guard_competitor = 1.00478774337695`
   6. `shared_gate_nway_shared_triggers = 1.00546821599453`

## 6B Exit Target (Copied From Phase Contract)
1. No workload above `1.02`.
2. No P1 workload above `1.00`.
3. Pass 3 consecutive benchmark runs with fixed settings:
   1. `ACCUMULATR_BENCH_N_REP=6`
   2. `ACCUMULATR_BENCH_TARGET_SEC=0.10`
   3. `ACCUMULATR_BENCH_MIN_INNER_REPS=40`

## Execution Batches (Strict Order)

### B0. Freeze Baseline + Safety Rails
1. Keep unified structural gates green before and after each batch:
   1. `check_step3_structure.sh`
   2. `check_step4_structure.sh`
   3. `check_step5_structure.sh`
   4. `check_step6_structure.sh`
2. Capture baseline benchmark artifact once per attempt cycle.

### B1. Hotspot Accounting (Mandatory Before More Tuning)
1. Add lightweight counters/timing in unified outcome path:
   1. number of `evaluate_outcome_coupling_unified` invocations per workload
   2. number of `kernel_node_density_entry_idx` evaluations
   3. adaptive segments accepted/split
2. Produce workload-level report for:
   1. `example_21_simple_q`
   2. `example_22_shared_q`
   3. `example_3_stop_na`
3. Do not change math yet; this batch is measurement only.

### B2. Remove Redundant Work In Unified Path
1. Eliminate repeated per-interval allocations in adaptive loop.
2. Reuse scratch vectors across intervals and contributions.
3. Hoist invariant pieces (competitor filtering, resolver outputs, trigger-state invariants).
4. Avoid repeated gate-CDF computations within same state where unchanged.

### B3. Adaptive Policy Redesign (Unified Engine Only)
1. Replace expensive split decision path with low-overhead estimator that avoids doubling evaluations.
2. Use conservative split budget and early-accept criteria tuned to RT scale.
3. Keep finite and infinite handling inside same unified integrator contract.
4. No fallback to scalar branch allowed.

### B4. Unified Fast Paths As Operator Modes (Not Separate Engines)
1. Add zero-competitor generic mode inside unified evaluator (still same entrypoint).
2. Add direct event-node shortcut mode when generic node maps to simple event form.
3. Gate each mode with exact behavior checks and golden validation.

### B5. Regression Closure Pass
1. Run `bench_centralized.R`.
2. Confirm top regressions recovered:
   1. `example_21_simple_q`
   2. `example_22_shared_q`
   3. `example_3_stop_na`
3. If still failing, iterate B2-B4 with evidence, no architecture changes.

### B6. 6B Stability Gate
1. Run 3 consecutive centralized benchmark passes under required env settings.
2. All runs must satisfy 6B thresholds.
3. Record completion block in `dev/plan_phase3.MD`.

## Hard Stop Conditions
1. Any optimization that reintroduces alternate scalar engine is rejected.
2. Any correctness failure blocks merge of optimization batch.
3. Any "improvement" claim without updated benchmark artifacts is rejected.

## Required Evidence Per Batch
1. Commit SHA
2. UTC timestamp
3. Batch ID (`B0`..`B6`)
4. Structural gate statuses
5. Golden status + `max_abs_diff`
6. Benchmark summary:
   1. worst ratio
   2. P1 ratios
7. Explicit statement: `No re-splitting introduced` (`yes/no`)

## Completion Criteria
1. Unified structural gates remain PASS.
2. Golden remains PASS (`max_abs_diff <= 1e-4`).
3. 6B thresholds and stability runs all PASS.
4. Only then mark Phase 3 complete in `dev/plan_phase3.MD`.

## Active Sub-Artifacts
1. GenericNodeIntegral NodeRef batch cutover (no scalar fallback): `dev/step6b_generic_noderef_batch_cutover_2026-03-05.md`

## Evidence Log

### Evidence Block: B0 Baseline Freeze
1. Commit SHA: `8a7c6d1` (working tree dirty)
2. UTC timestamp: `20260305T131424Z`
3. Batch ID: `B0`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_b0_step3_structure_20260305T131424Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_b0_step4_structure_20260305T131424Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_b0_step5_structure_20260305T131424Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_b0_step6_structure_20260305T131424Z.log`)
5. Golden:
   1. status: PASS (`dev/scripts/scratch_outputs/step6b_b0_golden_20260305T131424Z.log`)
   2. `max_abs_diff = 1.2372e-05`
6. Benchmark summary (`median_per_eval_sec`, centralized diff):
   1. worst ratio: `2.17088375504221` (`example_21_simple_q`)
   2. P1 ratios:
      1. `example_3_stop_na = 1.74117647058819`
      2. `example_21_simple_q = 2.17088375504221`
      3. `example_22_shared_q = 1.72521739130433`
      4. `example_6_dual_path = 0.979710144927539`
      5. `depth3_guard_competitor = 0.972231088413661`
      6. `shared_gate_nway_shared_triggers = 1.05434039644566`
7. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: B1 Hotspot Accounting
1. Commit SHA: `8a7c6d1` (working tree dirty, instrumentation-only changes)
2. UTC timestamps:
   1. structure+initial golden/probe run: `20260305T131827Z`
   2. isolated golden rerun: `20260305T132416Z`
   3. targeted probe with emitted stats: `20260305T132329Z`
3. Batch ID: `B1`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_b1_step3_structure_20260305T131827Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_b1_step4_structure_20260305T131827Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_b1_step5_structure_20260305T131827Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_b1_step6_structure_20260305T131827Z.log`)
5. Golden:
   1. initial run status: FAIL (`dev/scripts/scratch_outputs/step6b_b1_golden_20260305T131827Z.log`) due environment build-tools probe, not numerical diff assertion
   2. isolated run status: PASS (`dev/scripts/scratch_outputs/step6b_b1_golden_isolated_20260305T132416Z.log`)
   3. `max_abs_diff = 1.2372e-05`
6. B1 required workload report (targeted examples probe):
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_b1_examples_probe_20260305T132329Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_b1_unified_stats_20260305T132329Z.csv`
   3. `example_21_simple_q`:
      1. `evaluate_outcome_coupling_unified_calls_per_eval = 4`
      2. `kernel_node_density_entry_idx_calls_per_eval = 3154`
      3. `adaptive_segments_accepted_per_eval = 0`
      4. `adaptive_segments_split_per_eval = 0`
   4. `example_22_shared_q`:
      1. `evaluate_outcome_coupling_unified_calls_per_eval = 2`
      2. `kernel_node_density_entry_idx_calls_per_eval = 1599`
      3. `adaptive_segments_accepted_per_eval = 0`
      4. `adaptive_segments_split_per_eval = 0`
   5. `example_3_stop_na`:
      1. `evaluate_outcome_coupling_unified_calls_per_eval = 2`
      2. `kernel_node_density_entry_idx_calls_per_eval = 1519`
      3. `adaptive_segments_accepted_per_eval = 0`
      4. `adaptive_segments_split_per_eval = 0`
7. B1 interpretation:
   1. On the three regressed workloads, adaptive splitter activity is zero.
   2. The dominant operation count is repeated `kernel_node_density_entry_idx` evaluation volume.
   3. This points B2 optimization toward redundant node-density work elimination, not adaptive split-policy tuning.
8. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: B2 Unified Generic Batch Fast Path
1. Commit SHA: `8a7c6d1` (working tree dirty, optimization in unified generic path)
2. UTC timestamps:
   1. structure+golden+targeted probe: `20260305T134026Z`
   2. centralized benchmark: `20260305T134127Z`
   3. targeted probe with explicit fastpath counters: `20260305T134259Z`
3. Batch ID: `B2`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_b2_step3_structure_20260305T134026Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_b2_step4_structure_20260305T134026Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_b2_step5_structure_20260305T134026Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_b2_step6_structure_20260305T134026Z.log`)
5. Golden:
   1. install log: `dev/scripts/scratch_outputs/step6b_b2_install_golden_20260305T134026Z.log`
   2. golden log: `dev/scripts/scratch_outputs/step6b_b2_golden_20260305T134026Z.log`
   3. status: PASS
   4. `max_abs_diff = 1.2372e-05`
6. Centralized benchmark summary (`median_per_eval_sec`, post-B2):
   1. worst ratio: `example_21_simple_q = 1.111111`
   2. P1 ratios:
      1. `example_3_stop_na = 1.047059`
      2. `example_21_simple_q = 1.111111`
      3. `example_22_shared_q = 1.099130`
      4. `example_6_dual_path = 0.973913`
      5. `depth3_guard_competitor = 1.011810`
      6. `shared_gate_nway_shared_triggers = 1.051948`
7. Targeted workload counters (post-B2) vs B1:
   1. `example_21_simple_q`
      1. kernel density calls/eval: `3154 -> 454`
      2. generic fastpath/eval: `0 -> 4` (equal to coupling calls/eval)
   2. `example_22_shared_q`
      1. kernel density calls/eval: `1599 -> 249`
      2. generic fastpath/eval: `0 -> 2`
   3. `example_3_stop_na`
      1. kernel density calls/eval: `1519 -> 169`
      2. generic fastpath/eval: `0 -> 2`
8. Interpretation:
   1. Unified generic coupling now executes through batched labelref path for event-backed shapes.
   2. Major regression source reduced, but 6B threshold still not met (worst `1.111111`, target `<= 1.02`).
9. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: B3/B4 Adaptive Infinite Policy + Unified Operator Modes
1. Commit SHA: `8a7c6d1` (working tree dirty, unified B3/B4 optimizations)
2. UTC timestamps:
   1. structure+golden: `20260305T144430Z`
   2. targeted+shape+centralized+NodeRef probe: `20260305T144507Z`
3. Batch ID: `B3/B4`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_b34_step3_structure_20260305T144430Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_b34_step4_structure_20260305T144430Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_b34_step5_structure_20260305T144430Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_b34_step6_structure_20260305T144430Z.log`)
5. Golden:
   1. install log: `dev/scripts/scratch_outputs/step6b_b34_install_golden_20260305T144430Z.log`
   2. golden log: `dev/scripts/scratch_outputs/step6b_b34_golden_20260305T144430Z.log`
   3. status: PASS
   4. `max_abs_diff = 1.2372e-05`
6. Unified-engine changes implemented:
   1. Infinite-bound coupling integration moved from fixed 660-node batch to adaptive transformed-domain segmentation.
   2. Generic and NWay unified paths now always use `integrate_coupling_mass_batch_adaptive` for finite and infinite bounds.
   3. Added operator-mode shortcuts inside unified path (zero-competitor EventRef/NodeRef mass via direct CDF).
   4. Added EventRef specialized accumulator batch mode inside generic provider to reduce per-call overhead.
7. Targeted regressed generic workloads (post-B3/B4):
   1. probe log: `dev/scripts/scratch_outputs/step6b_b34_targeted_probe_20260305T144507Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_b34_targeted_stats_20260305T144507Z.csv`
   3. `example_21_simple_q`: ratio `0.9144947`, adaptive `accept=20/eval`, `split=4/eval`
   4. `example_22_shared_q`: ratio `0.9739130`, adaptive `accept=14/eval`, `split=6/eval`
   5. `example_3_stop_na`: ratio `0.9791743`, adaptive `accept=12/eval`, `split=4/eval`
8. NodeRef cutover safety check (unchanged guarantees):
   1. probe log: `dev/scripts/scratch_outputs/step6b_b34_noderef_probe_20260305T144507Z.log`
   2. `generic_noderef_batch_calls=400`, `generic_scalar_fallback_calls=0`, `kernel_node_density_entry_idx_calls=0`
9. Centralized benchmark summary (`median_per_eval_sec`):
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_b34_bench_centralized_20260305T144507Z.log`
   2. worst ratio: `example_2_stop_mixture = 1.0491146`
   3. top regressions:
      1. `example_2_stop_mixture = 1.0491146`
      2. `example_10_exclusion = 1.0466807`
      3. `shared_gate_nway_shared_triggers = 1.0246070`
      4. `depth3_guard_competitor = 1.0146824`
10. Interpretation (brutal):
   1. B3/B4 recovered prior P1 generic coupling regressions strongly.
   2. Remaining worst regressions are now mostly on workloads with `coupling_calls_per_eval = 0`, i.e. outside unified coupling mass path.
   3. 6B is still not complete because `max ratio > 1.02` and P1 still has `>1.00` on shared-gate/depth cases.
11. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: B5 Non-Coupling Hot-Loop Reduction Attempt (State Reuse + Cache Reuse)
1. Commit SHA: `8a7c6d1` (working tree dirty, B5 hot-loop reductions in progress)
2. UTC timestamps:
   1. first closure pass (default centralized): `20260305T150149Z`
   2. fixed-gate probe runs (required env): `20260305T151250Z`, `20260305T151524Z`
3. Batch ID: `B5`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_b6t_step3_structure_20260305T151524Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_b6t_step4_structure_20260305T151524Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_b6t_step5_structure_20260305T151524Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_b6t_step6_structure_20260305T151524Z.log`)
5. Golden:
   1. status: PASS (`dev/scripts/scratch_outputs/step6b_b6t_golden_20260305T151524Z.log`)
   2. `max_abs_diff = 1.2372e-05`
6. B5 implementation details:
   1. Runtime unified stats counters are now disabled by default in hot loops, enabled only when explicit stats reset API is called.
   2. `node_density_with_competitors_from_state` now reuses the current `NodeEvalState` for competitor survival (avoids per-call competitor state reconstruction).
   3. Competitor cache lookup and transition-mutation metadata are reused at call sites instead of re-fetch/scanning inside each hot invocation.
   4. Attempted fused scalar fast mode (`target + competitors` single kernel batch) was benchmarked and rejected; it produced broad regressions and was removed.
7. Centralized benchmark summary (required fixed env: `N_REP=6`, `TARGET_SEC=0.10`, `MIN_INNER_REPS=40`):
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_b6t_bench_centralized_20260305T151524Z.log`
   2. worst ratio: `example_2_stop_mixture = 1.05205479452055`
   3. top regressions:
      1. `example_2_stop_mixture = 1.05205479452055`
      2. `example_10_exclusion = 1.03882612647710`
      3. `example_22_shared_q = 1.02956521739127`
      4. `depth3_guard_competitor = 1.02553463134376`
      5. `shared_gate_nway_shared_triggers = 1.02529049897471`
8. P1 status (same fixed env run):
   1. `example_3_stop_na = 1.01176470588234`
   2. `example_21_simple_q = 0.933333333333363`
   3. `example_22_shared_q = 1.02956521739127`
   4. `example_6_dual_path = 0.985507246376863`
   5. `depth3_guard_competitor = 1.02553463134376`
   6. `shared_gate_nway_shared_triggers = 1.02529049897471`
9. Interpretation (brutal):
   1. B5 removed measurable hot-loop bloat and kept correctness/structure green.
   2. 6B exit thresholds are still failing; remaining gap is concentrated in non-coupling/shared-trigger heavy workloads.
   3. Additional speed work is still required before B6 stability gate can be attempted.
10. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: B6 Pair/NWay Scalar Fallback Integral Removal
1. Commit SHA: `8a7c6d1` (working tree dirty, Pair/NWay scalar fallback removal in progress)
2. UTC timestamp: `20260305T152248Z`
3. Batch ID: `B6`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6c_pn_step3_structure_20260305T152248Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6c_pn_step4_structure_20260305T152248Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6c_pn_step5_structure_20260305T152248Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6c_pn_step6_structure_20260305T152248Z.log`)
5. Golden:
   1. status: PASS (`dev/scripts/scratch_outputs/step6c_pn_golden_20260305T152248Z.log`)
   2. `max_abs_diff = 1.2372e-05`
6. Consolidation change:
   1. Removed scalar `integrate_boost_fn_0_to_upper` fallback integration from `evaluate_coupling_mass_nway`.
   2. Removed scalar `integrate_boost_fn_0_to_upper` fallback integration from `evaluate_coupling_mass_pair`.
   3. Pair/NWay coupling mass now routes through batched adaptive integration only.
7. Centralized benchmark summary (required fixed env: `N_REP=6`, `TARGET_SEC=0.10`, `MIN_INNER_REPS=40`):
   1. benchmark log: `dev/scripts/scratch_outputs/step6c_pn_bench_centralized_20260305T152248Z.log`
   2. worst ratio: `example_2_stop_mixture = 1.05753424657533`
   3. top regressions:
      1. `example_2_stop_mixture = 1.05753424657533`
      2. `example_10_exclusion = 1.02843786521231`
      3. `example_6_dual_path = 1.02608695652174`
      4. `example_7_mixture = 1.01993761019939`
      5. `shared_gate_nway_shared_triggers = 1.00546821599453`
8. P1 status:
   1. `example_3_stop_na = 1.00000000000000`
   2. `example_21_simple_q = 0.925925925925894`
   3. `example_22_shared_q = 1.00173913043479`
   4. `example_6_dual_path = 1.02608695652174`
   5. `depth3_guard_competitor = 1.00478774337695`
   6. `shared_gate_nway_shared_triggers = 1.00546821599453`
9. Interpretation (brutal):
   1. Structural scalar Pair/NWay fallback integration is removed.
   2. Correctness and structure remain green.
   3. 6B thresholds are still not met (`max > 1.02` and one P1 case remains >1.00).
10. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: B6 Stability Gate Attempt (3 Consecutive Fixed-Setting Runs)
1. Commit SHA: `8a7c6d1` (working tree dirty, 6A closure + unified hot-path normalization active)
2. UTC bundle timestamp: `20260305T171543Z`
3. Batch ID: `B6 Stability`
4. Required settings used (all three runs):
   1. `ACCUMULATR_BENCH_N_REP=6`
   2. `ACCUMULATR_BENCH_TARGET_SEC=0.10`
   3. `ACCUMULATR_BENCH_MIN_INNER_REPS=40`
5. Run artifacts:
   1. summary: `dev/scripts/scratch_outputs/step6b_stability_summary_20260305T171543Z.txt`
   2. run1 log: `dev/scripts/scratch_outputs/step6b_stability_run1_20260305T171543Z.log`
   3. run1 diff: `dev/scripts/scratch_outputs/benchmark_centralized_diff_run1_20260305T171543Z.csv`
   4. run2 log: `dev/scripts/scratch_outputs/step6b_stability_run2_20260305T171658Z.log`
   5. run2 diff: `dev/scripts/scratch_outputs/benchmark_centralized_diff_run2_20260305T171658Z.csv`
   6. run3 log: `dev/scripts/scratch_outputs/step6b_stability_run3_20260305T171810Z.log`
   7. run3 diff: `dev/scripts/scratch_outputs/benchmark_centralized_diff_run3_20260305T171810Z.csv`
6. Per-run gate results (`median_per_eval_sec`):
   1. Run 1: FAIL
      1. worst: `nested/shared_gate_pair = 1.181506849315`
      2. workloads `> 1.02`: `3`
      3. P1 `> 1.00`: `1` (`shared_gate_nway_shared_triggers = 1.112440191388`)
   2. Run 2: FAIL
      1. worst: `examples_models/example_2_stop_mixture = 1.026054257319`
      2. workloads `> 1.02`: `2`
      3. P1 `> 1.00`: `1` (`depth3_guard_competitor = 1.007660389403`)
   3. Run 3: FAIL
      1. worst: `examples_models/example_2_stop_mixture = 1.020682245501`
      2. workloads `> 1.02`: `1`
      3. P1 `> 1.00`: `1` (`depth3_guard_competitor = 1.007660389403`)
7. Final stability verdict:
   1. `STEP6B_STABILITY_STATUS=FAIL`
8. Interpretation (brutal):
   1. 6B closure gate is not achieved.
   2. Failures now concentrate in non-coupling workloads (`shared_gate_pair/nway`, `example_2_stop_mixture`) plus residual P1 drift (`depth3_guard_competitor` or `shared_gate_nway_shared_triggers` depending on run).
9. Explicit statement: `No re-splitting introduced = yes`
