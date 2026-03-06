# Step 6A-F Cutover Artifact: Unified GenericNodeIntegral Runtime (No Scalar Fallback)
Date: 2026-03-05
Status: IN PROGRESS (`F0`..`F4` complete; `F5` benchmark shows major regressions)
Owner: Codex + reviewer

## Why This Exists
1. Step 6A is currently blocked by `dev/scripts/check_step6_structure.sh` because `native_outcome_probability_bits_impl_idx` still contains scalar adaptive fallback (`integrate_boost_fn_0_to_upper`).
2. Recent speed recovery came from reintroducing that fallback, which improves throughput but violates the consolidation contract.
3. This artifact defines the strict cutover plan to recover a single engine without sacrificing correctness discipline.

## Non-Negotiable Contract
1. One runtime evaluator path for outcome mass with competitors.
2. `native_outcome_probability_bits_impl_idx` must always route competitor handling through `evaluate_outcome_coupling_unified`.
3. No scalar fallback integrator branch in `native_outcome_probability_bits_impl_idx`.
4. Adaptive behavior must exist inside the unified batched engine, not as a second execution engine.
5. Correctness gate (`max_abs_diff <= 1e-4`) is mandatory after every significant batch.

## Scope Of This Cutover
1. Add real runtime `GenericNodeIntegral` coupling operation in unified program model.
2. Extend `evaluate_outcome_coupling_unified` to execute Pair, NWay, and GenericNodeIntegral using one evaluator stack:
   1. `evaluate_labelref_batch`
   2. shared batch integration helper(s)
3. Replace fixed finite batch policy with adaptive batched refinement:
   1. coarse batch
   2. local error estimate
   3. split only high-error intervals
   4. recompute batched nodes
4. Delete scalar fallback path from `native_outcome_probability_bits_impl_idx`.

## Execution Batches (Strict Order)

### F0. Baseline Capture (Required Before Edits)
1. Capture current gate failure artifact:
   1. `bash dev/scripts/check_step6_structure.sh .`
2. Capture correctness baseline:
   1. `Rscript dev/scripts/check_loglik_golden.R`
3. Capture benchmark baseline:
   1. `Rscript dev/scripts/bench_centralized.R`
4. Save outputs under `dev/scripts/scratch_outputs/`.

### F1. Runtime Program Extension (No Behavior Change Yet)
1. Extend runtime coupling op model to include `GenericNodeIntegral`.
2. Add payload needed for generic target + competitor evaluation in unified path.
3. Ensure resolver can emit GenericNodeIntegral program when Pair/NWay specializations do not apply.
4. Do not delete old branch in this batch.
5. Gate after F1:
   1. `Rscript dev/scripts/check_loglik_golden.R` must pass.

### F2. Unified Evaluator Support For GenericNodeIntegral
1. Add GenericNodeIntegral branch in `evaluate_outcome_coupling_unified`.
2. Reuse existing batch label evaluation helpers; no bespoke scalar evaluator.
3. Keep Pair/NWay behavior unchanged.
4. Gate after F2:
   1. `bash dev/scripts/check_step5_structure.sh .`
   2. `Rscript dev/scripts/check_loglik_golden.R`

### F3. Adaptive Batched Integration Policy
1. Implement adaptive finite-range batching in shared integrator utilities:
   1. start from coarse finite segmentation
   2. estimate local quadrature error per interval
   3. split only intervals above tolerance
   4. enforce max refinement depth/interval cap to avoid runaway cost
2. Keep API compatible with unified evaluator callers.
3. No reintroduction of scalar fallback integration.
4. Gate after F3:
   1. `Rscript dev/scripts/check_loglik_golden.R`
   2. spot benchmark on `example_21_simple_q` and `example_22_shared_q`

### F4. Remove Scalar Fallback And Force Unified Routing
1. Delete scalar fallback branch in `native_outcome_probability_bits_impl_idx`.
2. Route all competitor outcome-mass evaluation through `evaluate_outcome_coupling_unified`.
3. Remove now-dead fallback plumbing.
4. Mandatory gates after F4:
   1. `bash dev/scripts/check_step3_structure.sh .`
   2. `bash dev/scripts/check_step4_structure.sh .`
   3. `bash dev/scripts/check_step5_structure.sh .`
   4. `bash dev/scripts/check_step6_structure.sh .`
   5. `Rscript dev/scripts/check_loglik_golden.R`

### F5. Consolidation Proof + Performance Pass
1. Run centralized benchmark:
   1. `Rscript dev/scripts/bench_centralized.R`
2. Verify 6B constraints on `median_per_eval_sec` ratios in:
   1. `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`
3. If any workload regresses:
   1. optimize only inside unified evaluator stack
   2. no branch resurrection

## Hard Stop Conditions
1. Any correctness failure blocks next batch.
2. Any attempt to reintroduce alternate scalar/generic engines blocks merge.
3. Any claim of Step 6A completion while `check_step6_structure.sh` fails is invalid.
4. Any "wrapper" migration that leaves legacy branch logic alive is invalid.

## Required Evidence Record For Each Batch
1. Commit SHA
2. UTC timestamp
3. Batch ID (`F0`..`F5`)
4. Gate command outputs (pass/fail)
5. Golden `max_abs_diff`
6. If benchmarked: worst ratio + P1 ratios
7. Deleted/added symbols summary

## Completion Criteria
1. Same commit (or contiguous final commit set) must satisfy:
   1. `check_step3_structure.sh` PASS
   2. `check_step4_structure.sh` PASS
   3. `check_step5_structure.sh` PASS
   4. `check_step6_structure.sh` PASS
   5. golden PASS with `max_abs_diff <= 1e-4`
2. 6B only starts after the above is true.

## Evidence Log

### F0 Baseline Capture
1. UTC timestamp: `2026-03-05T10:32:55Z`
2. Commit: `8a7c6d1`
3. Command: `bash dev/scripts/check_step6_structure.sh .`
   1. Result: FAIL
   2. Message: `Outcome probability function still contains scalar adaptive integration fallback.`
   3. Log: `dev/scripts/scratch_outputs/step6a_f0_step6_structure_20260305T103255Z.log`
4. Command: `Rscript dev/scripts/check_loglik_golden.R`
   1. Result: PASS
   2. `max_abs_diff=1.2372e-05` (tolerance `1e-4`)
   3. Log: `dev/scripts/scratch_outputs/step6a_f0_golden_20260305T103255Z.log`
5. Command: `Rscript dev/scripts/bench_centralized.R`
   1. Result: PASS
   2. Log: `dev/scripts/scratch_outputs/step6a_f0_bench_centralized_20260305T103255Z.log`
   3. Output artifacts refreshed:
      1. `dev/scripts/scratch_outputs/benchmark_centralized.csv`
      2. `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`
      3. `dev/scripts/scratch_outputs/benchmark_step3_comparison.csv`
6. Benchmark snapshot (`median_per_eval_sec`, current/baseline):
   1. Worst ratio: `1.029733424470` (`nested/shared_gate_nway_shared_triggers`)
   2. P1 ratios:
      1. `example_3_stop_na`: `1.000000000000`
      2. `example_21_simple_q`: `0.925925925926`
      3. `example_22_shared_q`: `0.956521739130`
      4. `example_6_dual_path`: `0.973340490249`
      5. `depth3_guard_competitor`: `0.980210660709`
      6. `shared_gate_nway_shared_triggers`: `1.029733424470`

### F1 Runtime Program Extension
1. UTC timestamp: `2026-03-05T10:41:24Z`
2. Commit: `8a7c6d1`
3. Scope completed:
   1. Added runtime `OutcomeCouplingGenericNodeIntegralPayload`.
   2. Extended `OutcomeCouplingOpKind` with `GenericNodeIntegral`.
   3. Added IR lookup helper: `ir_outcome_coupling_generic_lookup(...)`.
   4. Added resolver split:
      1. `resolve_outcome_coupling_program_impl(..., include_generic_runtime)`
      2. `resolve_outcome_coupling_program(...)` (current behavior, no generic at runtime yet)
      3. `resolve_outcome_coupling_program_with_generic(...)` (ready for F2 cutover)
4. Gate for F1:
   1. Command: `Rscript dev/scripts/check_loglik_golden.R`
   2. Result: PASS
   3. `max_abs_diff=1.2372e-05` (tolerance `1e-4`)
   4. Log: `dev/scripts/scratch_outputs/step6a_f1_golden_20260305T104124Z.log`

### F2 Unified Evaluator Support For GenericNodeIntegral
1. UTC timestamp: `2026-03-05T10:48:18Z`
2. Commit: `8a7c6d1`
3. Scope completed:
   1. Added `OutcomeCouplingOpKind::GenericNodeIntegral` execution branch in `evaluate_outcome_coupling_unified(...)`.
   2. Generic branch evaluates batched node densities over `build_coupling_time_batch_mass(...)` and integrates via `integrate_node_density_batch(...)`.
   3. Enabled generic resolver routing at runtime call sites:
      1. `native_outcome_probability_bits_impl_idx` now resolves with `resolve_outcome_coupling_program_with_generic(...)`.
      2. Nonresponse mass precompile path now resolves with `resolve_outcome_coupling_program_with_generic(...)`.
4. F2 gates:
   1. Command: `bash dev/scripts/check_step5_structure.sh .`
   2. Result: PASS
   3. Log: `dev/scripts/scratch_outputs/step6a_f2_step5_structure_20260305T104818Z.log`
5. Correctness:
   1. Command: `Rscript dev/scripts/check_loglik_golden.R`
   2. Result: PASS
   3. `max_abs_diff=1.2372e-05` (tolerance `1e-4`)
   4. Log: `dev/scripts/scratch_outputs/step6a_f2_golden_20260305T104818Z.log`
6. Boundary status (expected at F2):
   1. Command: `bash dev/scripts/check_step6_structure.sh .`
   2. Result: FAIL (`scalar adaptive integration fallback` still present; scheduled for F4/F3-F4 sequence)
   3. Log: `dev/scripts/scratch_outputs/step6a_f2_step6_structure_20260305T104818Z.log`

### F3 Adaptive Batched Integration Policy
1. UTC timestamp: `2026-03-05T10:55:21Z`
2. Commit: `8a7c6d1`
3. Scope completed:
   1. Added shared adaptive finite-range batch integrator utility for coupling mass.
   2. Implemented per-interval error estimation (`|fine - coarse|`) with selective interval splitting.
   3. Added refinement caps (max depth + max segment budget) to prevent runaway refinement.
   4. Switched unified `GenericNodeIntegral` finite integration path to adaptive batch policy.
   5. Switched unified `NWay` finite integration path to adaptive batch policy.
   6. Kept infinite-bound path on fixed transformed batch for now.
4. F3 gates:
   1. Command: `bash dev/scripts/check_step5_structure.sh .`
   2. Result: PASS
   3. Log: `dev/scripts/scratch_outputs/step6a_f3_step5_structure_20260305T105521Z.log`
5. Correctness:
   1. Command: `Rscript dev/scripts/check_loglik_golden.R`
   2. Result: PASS
   3. `max_abs_diff=1.2372e-05` (tolerance `1e-4`)
   4. Log: `dev/scripts/scratch_outputs/step6a_f3_golden_20260305T105521Z.log`
6. Boundary status (expected before F4):
   1. Command: `bash dev/scripts/check_step6_structure.sh .`
   2. Result: FAIL (`scalar adaptive integration fallback` still present)
   3. Log: `dev/scripts/scratch_outputs/step6a_f3_step6_structure_20260305T105521Z.log`
7. Spot benchmark readout (no spin):
   1. Source: `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`
   2. Captured: `dev/scripts/scratch_outputs/step6a_f3_spot_ratios_20260305T105521Z.txt`
   3. `example_21_simple_q`: `2.121627800640` (regressed)
   4. `example_22_shared_q`: `1.780178233864` (regressed)
   5. Worst observed ratio in current snapshot: `2.121627800640` (`example_21_simple_q`)

### F4 Remove Scalar Fallback And Force Unified Routing
1. UTC timestamp: `2026-03-05T10:58:49Z`
2. Commit: `8a7c6d1`
3. Scope completed:
   1. Deleted scalar adaptive fallback integration branch from `native_outcome_probability_bits_impl_idx`.
   2. `native_outcome_probability_bits_impl_idx` now routes through `evaluate_outcome_coupling_unified(...)` only.
   3. Resolver now emits generic runtime programs as structural fallback when Pair/NWay do not apply.
   4. Generic unified evaluator branch now supports forced-bit inputs (forced complete/survive propagation).
4. Mandatory F4 gates:
   1. `bash dev/scripts/check_step3_structure.sh .` -> PASS
      1. Log: `dev/scripts/scratch_outputs/step6a_f4_step3_structure_20260305T105849Z.log`
   2. `bash dev/scripts/check_step4_structure.sh .` -> PASS
      1. Log: `dev/scripts/scratch_outputs/step6a_f4_step4_structure_20260305T105849Z.log`
   3. `bash dev/scripts/check_step5_structure.sh .` -> PASS
      1. Log: `dev/scripts/scratch_outputs/step6a_f4_step5_structure_20260305T105849Z.log`
   4. `bash dev/scripts/check_step6_structure.sh .` -> PASS
      1. Log: `dev/scripts/scratch_outputs/step6a_f4_step6_structure_20260305T105849Z.log`
   5. `Rscript dev/scripts/check_loglik_golden.R` -> PASS
      1. `max_abs_diff=1.2372e-05` (tolerance `1e-4`)
      2. Log: `dev/scripts/scratch_outputs/step6a_f4_golden_20260305T105849Z.log`

### F5 Consolidation Proof + Performance Pass
1. UTC timestamp: `2026-03-05T10:58:49Z` (benchmark run immediately after F4 gates)
2. Commit: `8a7c6d1`
3. Command:
   1. `Rscript dev/scripts/bench_centralized.R`
4. Artifacts:
   1. `dev/scripts/scratch_outputs/benchmark_centralized.csv`
   2. `dev/scripts/scratch_outputs/benchmark_centralized_diff.csv`
   3. `dev/scripts/scratch_outputs/benchmark_step3_comparison.csv`
   4. Ratio summary: `dev/scripts/scratch_outputs/step6a_f5_bench_ratios_20260305T105849Z.txt`
5. Result: FAIL (performance gate not met)
6. Worst + P1 ratios (`median_per_eval_sec`, current/baseline):
   1. Worst: `example_21_simple_q = 2.166666666667`
   2. `example_3_stop_na = 1.794117647059`
   3. `example_21_simple_q = 2.166666666667`
   4. `example_22_shared_q = 1.728328382393`
   5. `example_6_dual_path = 0.962318840580`
   6. `depth3_guard_competitor = 0.976699648899`
   7. `shared_gate_nway_shared_triggers = 1.030075187970`
7. Honest readout:
   1. Structural consolidation gate is now passing (`F4` success).
   2. Performance is currently far outside `6B` targets.
   3. Next work must tune adaptive batch policy inside unified engine only (no path re-splitting).

### F5-A Adaptive Policy Retune Attempt (Failed)
1. UTC timestamp: `2026-03-05T11:06:14Z`
2. Commit: `8a7c6d1`
3. Change attempted:
   1. Reduced adaptive base segments and refinement caps.
   2. Switched local error estimator from explicit half-interval recomputation to cheap coarse-shape heuristic.
4. Structural/correctness status after change:
   1. `check_step6_structure.sh` PASS
   2. golden PASS (`max_abs_diff=1.2372e-05`)
5. Performance result: WORSE
   1. Ratio summary: `dev/scripts/scratch_outputs/step6a_f5a_bench_ratios_20260305T110614Z.txt`
   2. Worst: `example_21_simple_q = 2.267946959305`
   3. P1:
      1. `example_3_stop_na = 1.867647058824`
      2. `example_21_simple_q = 2.267946959305`
      3. `example_22_shared_q = 1.808695652174`
      4. `example_6_dual_path = 1.028985507246`
      5. `depth3_guard_competitor = 1.036386849665`
      6. `shared_gate_nway_shared_triggers = 1.104579630895`
6. Decision:
   1. Keep this iteration marked as failed optimization.
   2. Next optimization must target call-count reduction and per-outcome repeated work elimination inside unified path, not just quadrature threshold tweaks.
