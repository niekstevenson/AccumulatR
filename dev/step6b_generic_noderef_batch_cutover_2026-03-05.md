# Step 6B Sub-Artifact: GenericNodeIntegral NodeRef Batch Cutover (No Scalar Fallback)
Date: 2026-03-05
Status: IN PROGRESS (N4 complete, N6 performance gate failing)
Owner: Codex + reviewer
Parent artifact: `dev/step6b_unified_speed_recovery_2026-03-05.md`

## Why This Exists
1. `GenericNodeIntegral` now has an event-backed batch fast path and a scalar per-node fallback.
2. Keeping scalar fallback as a long-term path violates the consolidation goal.
3. This artifact defines strict cutover to one batch evaluator stack with provider modes only.

## Current Gap (Brutal)
1. Provider-only cutover is now in place (`EventRef`/`NodeRef` only) with strict hard-fail on invalid generic provider resolution.
2. Centralized suites still mostly exercise `EventRef` generic paths, so NodeRef batch pressure is concentrated in forced-bit/scenario-sensitive shapes.
3. Performance gate is still failing vs B2 on P1 examples (`example_21`, `example_3`, `example_22`), so consolidation is structurally further along but speed recovery remains incomplete.

## Non-Negotiable Contract
1. Keep one integration engine: `evaluate_outcome_coupling_unified`.
2. Keep one batch integration loop.
3. Allowed provider modes inside that loop only:
   1. `EventRefProvider`
   2. `NodeRefProvider`
4. No separate scalar per-node generic evaluator path after cutover.
5. Correctness (`max_abs_diff <= 1e-4`) must pass at each mergeable batch.

## Definition Of Done (Strict)
1. Generic scalar fallback call count is exactly zero on centralized benchmark runs.
2. Code deletion proof:
   1. scalar per-node generic loop removed from GenericNodeIntegral execution branch
   2. no alternate runtime dispatch to a second generic evaluator engine
3. Structural scripts pass:
   1. `dev/scripts/check_step3_structure.sh`
   2. `dev/scripts/check_step4_structure.sh`
   3. `dev/scripts/check_step5_structure.sh`
   4. `dev/scripts/check_step6_structure.sh`
4. Golden passes: `Rscript dev/scripts/check_loglik_golden.R` with `max_abs_diff <= 1e-4`.
5. Performance check after cutover must not regress versus current B2 state on:
   1. `example_3_stop_na`
   2. `example_21_simple_q`
   3. `example_22_shared_q`

## Required Runtime Counters
1. `evaluate_outcome_coupling_unified_calls`
2. `generic_labelref_batch_fastpath_calls`
3. `generic_noderef_batch_calls` (new)
4. `generic_scalar_fallback_calls` (new, target = 0)
5. `kernel_node_density_entry_idx_calls`

## Execution Batches (Strict Order)

### N0. Instrumentation Baseline
1. Add `generic_noderef_batch_calls` and `generic_scalar_fallback_calls` counters.
2. Emit counters in `unified_outcome_stats_cpp` and benchmark capture CSV.
3. Record centralized + targeted baseline before changing math.

### N1. Provider Contract Extraction
1. Define one generic term-provider contract consumed by one batch loop.
2. Keep EventRef provider as first implementation.
3. Introduce NodeRef provider scaffold (no behavior change yet).

### N2. NodeRef Batch Evaluator (Forced-Bit Safe)
1. Implement NodeRef provider to evaluate density/survival vectors over node batches.
2. Support arbitrary IR node ops (`And`, `Or`, `Not`, `Guard`, event nodes).
3. Preserve forced-complete and forced-survive semantics exactly.

### N3. Wire GenericNodeIntegral To Provider Modes
1. Route generic branch through provider selection only:
   1. EventRef provider when event-backed mapping exists
   2. NodeRef provider otherwise
2. Keep one integration loop; provider supplies terms only.

### N4. Delete Scalar Generic Fallback
1. Remove scalar per-node `kernel_node_density_entry_idx` generic fallback loop.
2. Enforce hard fail if provider resolution is impossible (no silent fallback engine).

### N5. Correctness + Structure Gates
1. Run all structure scripts.
2. Run isolated golden.
3. Record logs and max diff.

### N6. Performance Gate
1. Run targeted examples probe with unified stats capture.
2. Run centralized benchmark.
3. Required for this sub-artifact:
   1. `generic_scalar_fallback_calls == 0`
   2. no new slowdown vs B2 on `example_3_stop_na`, `example_21_simple_q`, `example_22_shared_q`

## Hard Stop Conditions
1. Any reintroduction of alternate generic runtime engine is rejected.
2. Any unresolved correctness drift blocks merge.
3. Any claim of cutover completion with non-zero scalar fallback count is rejected.

## Required Evidence Block Format
1. Commit SHA
2. UTC timestamp
3. Batch ID (`N0`..`N6`)
4. Structural gate status
5. Golden status + `max_abs_diff`
6. Counter snapshot:
   1. `generic_labelref_batch_fastpath_calls`
   2. `generic_noderef_batch_calls`
   3. `generic_scalar_fallback_calls`
7. Bench summary:
   1. worst ratio
   2. P1 ratios
8. Explicit statement: `Scalar generic fallback removed = yes/no`

## Completion Record (Fill Only When True)
1. N0 SHA: `TBD`
2. N4 deletion SHA: `TBD`
3. Zero-fallback verified at UTC: `TBD`
4. Golden max_abs_diff: `TBD`
5. Centralized worst ratio post-cutover: `TBD`

## Evidence Log

### Evidence Block: N0 Instrumentation Baseline
1. Commit SHA: `8a7c6d1` (working tree dirty, instrumentation-only for N0 counters)
2. UTC timestamp: `20260305T135113Z`
3. Batch ID: `N0`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_n0_step3_structure_20260305T135113Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_n0_step4_structure_20260305T135113Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_n0_step5_structure_20260305T135113Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_n0_step6_structure_20260305T135113Z.log`)
5. Golden:
   1. install log: `dev/scripts/scratch_outputs/step6b_n0_install_golden_20260305T135113Z.log`
   2. golden log: `dev/scripts/scratch_outputs/step6b_n0_golden_20260305T135113Z.log`
   3. status: PASS
   4. `max_abs_diff = 1.2372e-05`
6. Centralized benchmark summary (`median_per_eval_sec`):
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_n0_bench_centralized_20260305T135113Z.log`
   2. worst ratio: `example_21_simple_q = 1.155555555556`
   3. P1 ratios:
      1. `example_21_simple_q = 1.15555555555554`
      2. `example_3_stop_na = 1.15294117647058`
      3. `example_22_shared_q = 1.12695652173914`
      4. `shared_gate_nway_shared_triggers = 1.07587149692413`
      5. `depth3_guard_competitor = 1.02393871688478`
      6. `example_6_dual_path = 0.997101449275346`
7. Targeted counter baseline (regressed examples):
   1. probe log: `dev/scripts/scratch_outputs/step6b_n0_targeted_probe_20260305T135113Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_n0_targeted_stats_20260305T135113Z.csv`
   3. `example_21_simple_q`:
      1. `generic_labelref_batch_fastpath_calls_per_eval = 4`
      2. `generic_noderef_batch_calls_per_eval = 0`
      3. `generic_scalar_fallback_calls_per_eval = 0`
   4. `example_22_shared_q`:
      1. `generic_labelref_batch_fastpath_calls_per_eval = 2`
      2. `generic_noderef_batch_calls_per_eval = 0`
      3. `generic_scalar_fallback_calls_per_eval = 0`
   5. `example_3_stop_na`:
      1. `generic_labelref_batch_fastpath_calls_per_eval = 2`
      2. `generic_noderef_batch_calls_per_eval = 0`
      3. `generic_scalar_fallback_calls_per_eval = 0`
8. Shape scan counter baseline (full examples suite):
   1. probe log: `dev/scripts/scratch_outputs/step6b_n0_shape_scan_20260305T135113Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_n0_shape_scan_20260305T135113Z.csv`
   3. `generic_scalar_fallback_calls = 0` for all examples in current scan
9. Explicit statement: `Scalar generic fallback removed = no`

### Evidence Block: N1 Provider Contract Extraction (No Math Change)
1. Commit SHA: `8a7c6d1` (working tree dirty, provider refactor only)
2. UTC timestamps:
   1. structure+golden+targeted+shape scan: `20260305T140312Z`
   2. centralized benchmark pass: `20260305T140457Z`
3. Batch ID: `N1`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_n1_step3_structure_20260305T140312Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_n1_step4_structure_20260305T140312Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_n1_step5_structure_20260305T140312Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_n1_step6_structure_20260305T140312Z.log`)
5. Golden:
   1. install log: `dev/scripts/scratch_outputs/step6b_n1_install_golden_20260305T140312Z.log`
   2. golden log: `dev/scripts/scratch_outputs/step6b_n1_golden_20260305T140312Z.log`
   3. status: PASS
   4. `max_abs_diff = 1.2372e-05`
6. Centralized benchmark summary (`median_per_eval_sec`):
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_n1_bench_centralized_20260305T140457Z.log`
   2. worst ratio: `example_3_stop_na = 1.129411764706`
   3. P1 ratios:
      1. `example_3_stop_na = 1.12941176470588`
      2. `example_21_simple_q = 1.12592592592593`
      3. `example_22_shared_q = 1.11304347826089`
      4. `shared_gate_nway_shared_triggers = 1.05468215994532`
      5. `example_6_dual_path = 1.04462620175062`
      6. `depth3_guard_competitor = 1.02776891158634`
7. Counter parity vs N0 on targeted generic cases:
   1. targeted probe log: `dev/scripts/scratch_outputs/step6b_n1_targeted_probe_20260305T140312Z.log`
   2. targeted stats: `dev/scripts/scratch_outputs/step6b_n1_targeted_stats_20260305T140312Z.csv`
   3. parity result:
      1. `coupling_calls_per_eval`: unchanged (`4/2/2` for `example_21/22/3`)
      2. `kernel_density_calls_per_eval`: unchanged (`454/249/169`)
      3. `generic_fastpath_per_eval`: unchanged (`4/2/2`)
      4. `generic_noderef_batch_per_eval`: unchanged (`0/0/0`)
      5. `generic_scalar_fallback_per_eval`: unchanged (`0/0/0`)
8. Shape scan:
   1. probe log: `dev/scripts/scratch_outputs/step6b_n1_shape_scan_20260305T140312Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_n1_shape_scan_20260305T140312Z.csv`
   3. `generic_scalar_fallback_calls = 0` across scan
9. Explicit statement: `Scalar generic fallback removed = no`

### Evidence Block: N2 NodeRef Batch Evaluator (Forced-Bit Safe)
1. Commit SHA: `8a7c6d1` (working tree dirty, N2 batch evaluator wiring in progress)
2. UTC timestamps:
   1. structure+golden: `20260305T142557Z`
   2. targeted+shape+centralized benchmark: `20260305T142631Z`
   3. forced-bit NodeRef probe: `20260305T143019Z`
3. Batch ID: `N2`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_n2_step3_structure_20260305T142557Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_n2_step4_structure_20260305T142557Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_n2_step5_structure_20260305T142557Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_n2_step6_structure_20260305T142557Z.log`)
5. Golden:
   1. install log: `dev/scripts/scratch_outputs/step6b_n2_install_golden_20260305T142557Z.log`
   2. golden log: `dev/scripts/scratch_outputs/step6b_n2_golden_20260305T142557Z.log`
   3. status: PASS
   4. `max_abs_diff = 1.2372e-05`
6. Targeted regressed examples counters (still EventRef-dominant):
   1. probe log: `dev/scripts/scratch_outputs/step6b_n2_targeted_probe_20260305T142631Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_n2_targeted_stats_20260305T142631Z.csv`
   3. `example_21_simple_q`: `generic_fastpath=4/eval`, `generic_noderef_batch=0/eval`, `generic_scalar_fallback=0/eval`, `kernel_density=454/eval`
   4. `example_22_shared_q`: `generic_fastpath=2/eval`, `generic_noderef_batch=0/eval`, `generic_scalar_fallback=0/eval`, `kernel_density=249/eval`
   5. `example_3_stop_na`: `generic_fastpath=2/eval`, `generic_noderef_batch=0/eval`, `generic_scalar_fallback=0/eval`, `kernel_density=169/eval`
7. Forced-bit NodeRef provider activation proof:
   1. probe log: `dev/scripts/scratch_outputs/step6b_n2_noderef_probe_20260305T143019Z.log`
   2. `forced_complete` run (100 calls): `generic_noderef_batch_calls=400`, `generic_scalar_fallback_calls=0`, `kernel_node_density_entry_idx_calls=0`
   3. `forced_survive` run (100 calls): `generic_noderef_batch_calls=400`, `generic_scalar_fallback_calls=0`, `kernel_node_density_entry_idx_calls=0`
8. Step 6 coupling test file:
   1. log: `dev/scripts/scratch_outputs/step6b_n2_test_step6_unified_coupling_20260305T142631Z.log`
   2. status: PASS (`15` assertions)
9. Centralized benchmark summary (`median_per_eval_sec`):
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_n2_bench_centralized_20260305T142631Z.log`
   2. worst ratio: `example_21_simple_q = 1.0925926`
   3. P1 ratios:
      1. `example_21_simple_q = 1.0925926`
      2. `example_3_stop_na = 1.0748003`
      3. `shared_gate_nway_shared_triggers = 1.0693780`
      4. `example_22_shared_q = 1.0434783`
      5. `depth3_guard_competitor = 1.0175551`
      6. `example_6_dual_path = 0.9710145`
10. Shape scan:
   1. probe log: `dev/scripts/scratch_outputs/step6b_n2_shape_scan_20260305T142631Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_n2_shape_scan_20260305T142631Z.csv`
   3. `generic_scalar_fallback_calls = 0` across current examples scan
11. Explicit statement: `Scalar generic fallback removed = no` (loop removed in generic provider execution path, but cutover/deletion proof not closed until N4)

### Evidence Block: N3/N4 Provider-Only Cutover + Scalar-Fallback Scaffolding Removal
1. Commit SHA: `8a7c6d1` (working tree dirty, N3/N4 hardening)
2. UTC timestamps:
   1. structure+golden: `20260305T143400Z`
   2. targeted+shape+centralized+NodeRef probe: `20260305T143439Z`
3. Batch ID: `N3/N4`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_n34_step3_structure_20260305T143400Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_n34_step4_structure_20260305T143400Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_n34_step5_structure_20260305T143400Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_n34_step6_structure_20260305T143400Z.log`)
5. Golden:
   1. install log: `dev/scripts/scratch_outputs/step6b_n34_install_golden_20260305T143400Z.log`
   2. golden log: `dev/scripts/scratch_outputs/step6b_n34_golden_20260305T143400Z.log`
   3. status: PASS
   4. `max_abs_diff = 1.2372e-05`
6. Provider-only dispatch changes:
   1. removed `GenericCouplingProviderKind::None` and `provider.valid` runtime branching
   2. generic provider resolution now hard-fails on invalid node id (no silent `0.0` fallback)
   3. scalar-fallback dispatch increment path removed from `evaluate_outcome_coupling_unified`
7. Targeted regressed examples counters:
   1. probe log: `dev/scripts/scratch_outputs/step6b_n34_targeted_probe_20260305T143439Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_n34_targeted_stats_20260305T143439Z.csv`
   3. `example_21_simple_q`: `generic_fastpath=4/eval`, `generic_noderef_batch=0/eval`, `generic_scalar_fallback=0/eval`, `kernel_density=454/eval`
   4. `example_22_shared_q`: `generic_fastpath=2/eval`, `generic_noderef_batch=0/eval`, `generic_scalar_fallback=0/eval`, `kernel_density=249/eval`
   5. `example_3_stop_na`: `generic_fastpath=2/eval`, `generic_noderef_batch=0/eval`, `generic_scalar_fallback=0/eval`, `kernel_density=169/eval`
8. Forced-bit NodeRef provider activation proof:
   1. probe log: `dev/scripts/scratch_outputs/step6b_n34_noderef_probe_20260305T143439Z.log`
   2. `forced_complete` run (100 calls): `generic_noderef_batch_calls=400`, `generic_scalar_fallback_calls=0`, `kernel_node_density_entry_idx_calls=0`
   3. `forced_survive` run (100 calls): `generic_noderef_batch_calls=400`, `generic_scalar_fallback_calls=0`, `kernel_node_density_entry_idx_calls=0`
9. Centralized benchmark summary (`median_per_eval_sec`):
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_n34_bench_centralized_20260305T143439Z.log`
   2. worst ratio: `example_3_stop_na = 1.161946`
   3. P1 ratios:
      1. `example_3_stop_na = 1.161946`
      2. `example_21_simple_q = 1.131815`
      3. `example_6_dual_path = 1.057971`
      4. `shared_gate_nway_shared_triggers = 1.056049`
      5. `example_22_shared_q = 1.043478`
      6. `depth3_guard_competitor = 1.017555`
10. Shape scan:
   1. probe log: `dev/scripts/scratch_outputs/step6b_n34_shape_scan_20260305T143439Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_n34_shape_scan_20260305T143439Z.csv`
   3. `generic_scalar_fallback_calls = 0` across current examples scan
11. Explicit statement: `Scalar generic fallback removed = yes` (runtime dispatch + per-node scalar generic loop removed; one unified batch evaluator stack remains)

## Fast Generic Reference Cases (For N1-N4)
1. Source artifacts:
   1. probe log: `dev/scripts/scratch_outputs/step6b_generic_ref_probe_20260305T135825Z.log`
   2. stats csv: `dev/scripts/scratch_outputs/step6b_generic_ref_stats_20260305T135825Z.csv`
   3. merged reference table: `dev/scripts/scratch_outputs/step6b_generic_reference_20260305T135825Z.csv`
2. Fastest generic-targeting workloads (`median_per_eval_sec`):
   1. `example_21_simple_q = 0.000370` (generic coupling calls/eval = 4)
   2. `example_22_shared_q = 0.000370` (generic coupling calls/eval = 2)
   3. `example_3_stop_na = 0.000455` (generic coupling calls/eval = 2)
3. Recommended primary micro-reference:
   1. `example_21_simple_q` (fast and highest generic coupling pressure among fast cases)
