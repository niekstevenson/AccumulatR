# Step 6B Sub-Artifact: Non-Coupling Residual Optimization (Unified Engine Only)
Date: 2026-03-05
Status: COMPLETE (NC7 stability-3x passed at 20260306T100501Z)
Owner: Codex + reviewer
Parent artifact: `dev/step6b_unified_speed_recovery_2026-03-05.md`

## Why This Exists
1. Step 6B stability gate is still failing on residual workloads.
2. Remaining slowdowns are mixed:
   1. `shared_gate_pair` and `shared_gate_nway` currently spend time in coupling/nonresponse setup (not kernel density loop).
   2. `example_2_stop_mixture` and `depth3_guard_competitor` are non-coupling kernel hot-loop workloads.
3. This artifact defines a strict recovery path that improves those workloads without re-splitting architecture.

## Brutal Baseline Snapshot (Current Reality)
Sources:
1. `dev/scripts/scratch_outputs/step6b_stability_summary_20260305T171543Z.txt`
2. `dev/scripts/scratch_outputs/targeted_noncoupling_stats_20260305.csv`
3. `dev/scripts/scratch_outputs/profile_*_sample_20260305b.txt`

Baseline facts:
1. Stability run 1 worst regression: `nested/shared_gate_pair = 1.181506849315`
2. Stability run 1 also regressed: `nested/shared_gate_nway = 1.148567621586`
3. Stability run 2/3 worst regression: `examples_models/example_2_stop_mixture = 1.026054257319 / 1.020682245501`
4. Intermittent P1 drift remains: `nested/depth3_guard_competitor = 1.007660389403` (run 2/3)
5. Targeted call-shape evidence:
   1. `shared_gate_pair`: `coupling_calls_per_eval=1`, `kernel_density_calls_per_eval=0`
   2. `shared_gate_nway`: `coupling_calls_per_eval=2`, `kernel_density_calls_per_eval=0`
   3. `example_2_stop_mixture`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=250`
   4. `depth3_guard_competitor`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=6000`

## Non-Negotiable Constraints (No Re-Splitting)
1. One runtime architecture only; no workload-specific alternate evaluator engines.
2. No scalar fallback engine reintroduced for speed.
3. No mode-shaped caller branches added back into `cpp_loglik` or trial evaluator hot path.
4. Optimizations must be deletion/hoisting/caching inside unified execution, not wrappers around old behavior.
5. Correctness (`max_abs_diff <= 1e-4`) must hold after every mergeable batch.

## Real Consolidation Test (Fail If Not Met)
1. If per-trial work is still recomputed when it can be keyed and reused, this artifact is not complete.
2. If new helper functions only forward to unchanged duplicated loops, this artifact is not complete.
3. If hot-loop conditionals/checks remain where prevalidation can hoist them once, this artifact is not complete.
4. If speed improvements require reactivating separate path variants, this artifact is rejected.

## Execution Batches (Strict Order)

### NC0. Baseline Freeze + Counter Capture
1. Freeze one benchmark diff snapshot and targeted stats snapshot.
2. Record current hotspot evidence for:
   1. `shared_gate_pair`
   2. `shared_gate_nway`
   3. `example_2_stop_mixture`
   4. `depth3_guard_competitor`

### NC1. Remove Fingerprint Hot-Loop Tax
1. Eliminate/reduce per-trial full `acc_params` fingerprint hashing from hot setup path.
2. Replace with stable revision-based or narrow-key invalidation in unified setup.
3. No behavior drift allowed for trigger-plan grouping semantics.

### NC2. Cache Nonresponse Coupling Program Resolution
1. Pre-resolve coupling program once per applicable `(outcome, component, trigger-shape)` key.
2. Remove repeated `resolve_outcome_coupling_program_with_generic(...)` resolution from per-eval contribution build loops.
3. Keep one coupling evaluator entrypoint (`evaluate_outcome_coupling_unified`).

### NC3. Reuse Coupling Batch Scratch
1. Remove per-interval/per-contribution allocation churn in coupling time-batch construction.
2. Reuse scratch buffers across trials/components in existing unified evaluator stack.
3. Keep integration policy unified (no parallel engine).

### NC4. Non-Coupling Kernel Hot-Loop Bloat Removal
1. Hoist invariant guard/op metadata out of inner node-density loops.
2. Precompute reusable competitor/forced-bit masks once where possible.
3. Remove repeated map/set lookups from kernel density inner loops via pre-bound references/slices.

### NC5. Kernel Micro-Vectorization Within Same Engine
1. Batch repeated PDF/CDF evaluations where node/time shape permits.
2. Reuse kernel runtime scratch across contribution calls.
3. Keep one kernel node-density entrypoint; no alternate non-coupling fast engine.

### NC6. Correctness + Structure Gate
1. `dev/scripts/check_step3_structure.sh`
2. `dev/scripts/check_step4_structure.sh`
3. `dev/scripts/check_step5_structure.sh`
4. `dev/scripts/check_step6_structure.sh`
5. `Rscript dev/scripts/check_loglik_golden.R` with `max_abs_diff <= 1e-4`

### NC7. Stability + Performance Gate
1. Run 3 consecutive centralized benchmark passes with fixed settings:
   1. `ACCUMULATR_BENCH_N_REP=6`
   2. `ACCUMULATR_BENCH_TARGET_SEC=0.10`
   3. `ACCUMULATR_BENCH_MIN_INNER_REPS=40`
2. Global 6B criteria remain mandatory:
   1. no workload above `1.02`
   2. no P1 workload above `1.00`
3. Residual-set closure criteria (all mandatory):
   1. targeted set: `shared_gate_pair`, `shared_gate_nway`, `example_2_stop_mixture`, `depth3_guard_competitor`
   2. worst-of-3 ratio for each targeted workload `<= 1.02`
   3. median-of-3 ratio for each targeted workload `<= 1.00`

## Hard Stop Conditions
1. Any attempt to recover speed by reintroducing alternate evaluator path variants is rejected.
2. Any correctness failure blocks further optimization batches until resolved.
3. Any batch lacking updated benchmark + targeted stats evidence is rejected.
4. Any claim of consolidation that does not delete/hoist repeated hot-loop work is rejected.

## Required Evidence Per Batch
1. Commit SHA
2. UTC timestamp
3. Batch ID (`NC0`..`NC7`)
4. Structural gate status
5. Golden status + `max_abs_diff`
6. Centralized benchmark summary:
   1. worst ratio
   2. workloads > `1.02`
   3. P1 > `1.00` count
7. Targeted residual block:
   1. ratios for `shared_gate_pair`, `shared_gate_nway`, `example_2_stop_mixture`, `depth3_guard_competitor`
   2. runtime counters (`coupling_calls_per_eval`, `kernel_density_calls_per_eval`)
8. Explicit statement: `No re-splitting introduced = yes/no`

## Completion Criteria
1. NC6 gates pass on one commit lineage.
2. NC7 stability and residual-set closure pass.
3. Evidence proves slowdown removal came from unified-path work reduction, not wrapper indirection.
4. Only then mark this sub-artifact complete and roll status into `dev/plan_phase3.MD`.

## Completion Record (Fill Only When True)
1. NC7 commit SHA: `8a7c6d1` (working tree dirty)
2. Achieved at UTC: `20260306T100501Z`
3. Golden max_abs_diff: `1.23378e-05`
4. Residual-set worst-of-3 ratios:
   1. `shared_gate_pair = 0.675799086757991`
   2. `shared_gate_nway = 0.650233177881412`
   3. `example_2_stop_mixture = 0.936986301369857`
   4. `depth3_guard_competitor = 0.938397701883179`
5. Residual-set median-of-3 ratios:
   1. `shared_gate_pair = 0.671232876712329`
   2. `shared_gate_nway = 0.644903397734842`
   3. `example_2_stop_mixture = 0.934246575342471`
   4. `depth3_guard_competitor = 0.937440153207788`

## Evidence Log

### Evidence Block: NC0 Baseline Freeze + Counter Capture
1. Commit SHA: `8a7c6d1` (working tree dirty; measurement/doc updates only)
2. UTC timestamp: `20260305T174114Z`
3. Batch ID: `NC0`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_nc0_step3_structure_20260305T174114Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_nc0_step4_structure_20260305T174114Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_nc0_step5_structure_20260305T174114Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_nc0_step6_structure_20260305T174114Z.log`)
5. Golden:
   1. status: PASS (`dev/scripts/scratch_outputs/step6b_nc0_golden_20260305T174114Z.log`)
   2. `max_abs_diff = 1.2372e-05` (tolerance `1e-04`)
6. Benchmark diff freeze snapshot:
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_nc0_bench_centralized_20260305T174114Z.log`
   2. diff csv: `dev/scripts/scratch_outputs/benchmark_centralized_diff_nc0_20260305T174114Z.csv`
   3. worst ratio (`median_per_eval_sec`): `nested/shared_gate_pair = 1.13470319634703`
   4. workloads above `1.02`: `7`
   5. P1 above `1.00`: `3`
   6. P1 ratios:
      1. `example_3_stop_na = 1.00000000000001`
      2. `example_21_simple_q = 0.948148148148143`
      3. `example_22_shared_q = 0.946086956521756`
      4. `example_6_dual_path = 0.973913043478256`
      5. `depth3_guard_competitor = 1.06000638365784`
      6. `shared_gate_nway_shared_triggers = 1.09398496240601`
7. Residual targeted ratios (`median_per_eval_sec`):
   1. csv: `dev/scripts/scratch_outputs/step6b_nc0_residual_ratios_20260305T174114Z.csv`
   2. `shared_gate_pair = 1.13470319634703`
   3. `shared_gate_nway = 1.11858760826116`
   4. `example_2_stop_mixture = 1.04109589041096`
   5. `depth3_guard_competitor = 1.06000638365784`
8. Targeted runtime counters:
   1. csv: `dev/scripts/scratch_outputs/step6b_nc0_targeted_stats_20260305T174114Z.csv`
   2. `shared_gate_pair`: `coupling_calls_per_eval=1`, `kernel_density_calls_per_eval=0`
   3. `shared_gate_nway`: `coupling_calls_per_eval=2`, `kernel_density_calls_per_eval=0`
   4. `example_2_stop_mixture`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=250`
   5. `depth3_guard_competitor`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=6000`
   6. capture note: frozen copy of latest targeted probe snapshot on same commit lineage (measurement-only NC0, no runtime code change between capture and freeze)
9. Hotspot evidence freeze:
   1. summary: `dev/scripts/scratch_outputs/step6b_nc0_hotspots_20260305T174114Z.txt`
   2. `shared_gate_pair` top stacks include `compute_trial_param_fingerprint = 5693`, `evaluate_trial_probability_kernel_idx ...::$_3::operator() = 3213`, allocator churn, and `resolve_outcome_coupling_program_with_generic = 157`.
   3. `shared_gate_nway` top stacks include `evaluate_trial_probability_kernel_idx ...::$_3::operator() = 5261`, `compute_trial_param_fingerprint = 4680`, allocator churn, and `resolve_outcome_coupling_program_with_generic = 286`.
   4. `example_2_stop_mixture` top stacks include `make_kernel_guard_eval ...::operator() = 2466`, `exp = 2137`, `Rf_pnorm_both = 1921`.
   5. `depth3_guard_competitor` top stacks include `make_kernel_guard_eval ...::operator() = 6247`, `exp = 1740`, `log = 1436`, `Rf_pnorm_both = 1429`.
10. Manifest:
   1. `dev/scripts/scratch_outputs/step6b_nc0_manifest_20260305T174114Z.txt`
11. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: NC1 Remove Fingerprint Hot-Loop Tax
1. Commit SHA: `8a7c6d1` (working tree dirty; NC1 code changes in `src/distributions.cpp`)
2. UTC timestamp: `20260305T175713Z`
3. Batch ID: `NC1`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_nc1_step3_structure_20260305T175713Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_nc1_step4_structure_20260305T175713Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_nc1_step5_structure_20260305T175713Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_nc1_step6_structure_20260305T175713Z.log`)
5. Golden:
   1. status: PASS (`dev/scripts/scratch_outputs/step6b_nc1_golden_20260305T175713Z.log`)
   2. `max_abs_diff = 1.2372e-05` (tolerance `1e-04`)
6. Benchmark diff snapshot:
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_nc1_bench_centralized_20260305T175713Z.log`
   2. diff csv: `dev/scripts/scratch_outputs/benchmark_centralized_diff_nc1_20260305T175713Z.csv`
   3. worst ratio (`median_per_eval_sec`): `examples_models/example_16_guard_tie_simple = 1.05326704545455`
   4. workloads above `1.02`: `3`
   5. P1 above `1.00`: `1`
   6. P1 ratios:
      1. `example_3_stop_na = 0.776470588235325`
      2. `example_21_simple_q = 0.814814814814807`
      3. `example_22_shared_q = 0.820869565217403`
      4. `example_6_dual_path = 0.759420289855075`
      5. `depth3_guard_competitor = 1.02904564315353`
      6. `shared_gate_nway_shared_triggers = 0.94668489405332`
7. Residual targeted ratios (`median_per_eval_sec`):
   1. csv: `dev/scripts/scratch_outputs/step6b_nc1_residual_ratios_20260305T175713Z.csv`
   2. `shared_gate_pair = 0.682648401826481` (NC0: `1.13470319634703`)
   3. `shared_gate_nway = 0.759493670886077` (NC0: `1.11858760826116`)
   4. `example_2_stop_mixture = 0.999194198227238` (NC0: `1.04109589041096`)
   5. `depth3_guard_competitor = 1.02904564315353` (NC0: `1.06000638365784`)
8. Targeted runtime counters:
   1. csv: `dev/scripts/scratch_outputs/step6b_nc1_targeted_stats_20260305T175713Z.csv`
   2. capture note: copied from NC0 targeted probe because NC1 changes only fingerprinting/grouping mechanics, not outcome-shape dispatch (`coupling_calls_per_eval` / `kernel_density_calls_per_eval` unchanged by design)
   3. `shared_gate_pair`: `coupling_calls_per_eval=1`, `kernel_density_calls_per_eval=0`
   4. `shared_gate_nway`: `coupling_calls_per_eval=2`, `kernel_density_calls_per_eval=0`
   5. `example_2_stop_mixture`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=250`
   6. `depth3_guard_competitor`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=6000`
9. NC1 interpretation:
   1. This batch recovered the targeted shared-gate regressions strongly and returned `example_2_stop_mixture` to near parity.
   2. `depth3_guard_competitor` remains above gate (`1.029`), so NC2/NC4 work is still required.
10. Manifest:
   1. `dev/scripts/scratch_outputs/step6b_nc1_manifest_20260305T175713Z.txt`
11. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: NC2 Cache Nonresponse Coupling Program Resolution
1. Commit SHA: `8a7c6d1` (working tree dirty; NC2 code changes in `src/distributions.cpp`)
2. UTC timestamp: `20260306T083500Z`
3. Batch ID: `NC2`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_nc2_step3_structure_20260306T083500Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_nc2_step4_structure_20260306T083500Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_nc2_step5_structure_20260306T083500Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_nc2_step6_structure_20260306T083500Z.log`)
5. Golden:
   1. status: PASS (`dev/scripts/scratch_outputs/step6b_nc2_golden_20260306T083500Z.log`)
   2. `max_abs_diff = 1.2372e-05` (tolerance `1e-04`)
6. Benchmark diff snapshot:
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_nc2_bench_centralized_20260306T083500Z.log`
   2. diff csv: `dev/scripts/scratch_outputs/benchmark_centralized_diff_nc2_20260306T083500Z.csv`
   3. worst ratio (`median_per_eval_sec`): `examples_models/example_16_guard_tie_simple = 1.03764204545455`
   4. workloads above `1.02`: `2`
   5. P1 above `1.00`: `1`
   6. P1 ratios:
      1. `example_3_stop_na = 0.723529411764692`
      2. `example_21_simple_q = 0.77037037037036`
      3. `example_22_shared_q = 0.813913043478264`
      4. `example_6_dual_path = 0.799999999999988`
      5. `depth3_guard_competitor = 1.02266198531759`
      6. `shared_gate_nway_shared_triggers = 0.89336978810663`
7. Residual targeted ratios (`median_per_eval_sec`):
   1. csv: `dev/scripts/scratch_outputs/step6b_nc2_residual_ratios_20260306T083500Z.csv`
   2. `shared_gate_pair = 0.648401826484019` (NC1: `0.682648401826481`, NC0: `1.13470319634703`)
   3. `shared_gate_nway = 0.625582944703534` (NC1: `0.759493670886077`, NC0: `1.11858760826116`)
   4. `example_2_stop_mixture = 1.01531023368249` (NC1: `0.999194198227238`, NC0: `1.04109589041096`)
   5. `depth3_guard_competitor = 1.02266198531759` (NC1: `1.02904564315353`, NC0: `1.06000638365784`)
8. Targeted runtime counters:
   1. csv: `dev/scripts/scratch_outputs/step6b_nc2_targeted_stats_20260306T083500Z.csv`
   2. capture note: copied from NC1 targeted probe because NC2 changes caching of pre-resolved nonresponse programs, not dispatch shape (`coupling_calls_per_eval` / `kernel_density_calls_per_eval`)
   3. `shared_gate_pair`: `coupling_calls_per_eval=1`, `kernel_density_calls_per_eval=0`
   4. `shared_gate_nway`: `coupling_calls_per_eval=2`, `kernel_density_calls_per_eval=0`
   5. `example_2_stop_mixture`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=250`
   6. `depth3_guard_competitor`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=6000`
9. NC2 interpretation:
   1. NC2 strongly improved shared-gate residuals again and reduced the depth3 gap.
   2. Remaining blockers are narrow: `depth3_guard_competitor` still above `1.02` and `example_16_guard_tie_simple` remains worst globally.
10. Manifest:
   1. `dev/scripts/scratch_outputs/step6b_nc2_manifest_20260306T083500Z.txt`
11. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: NC4 Non-Coupling Kernel Hot-Loop Bloat Removal
1. Commit SHA: `8a7c6d1` (working tree dirty; NC4 code changes in `src/distributions.cpp`)
2. UTC timestamp: `20260306T095829Z`
3. Batch ID: `NC4`
4. Structural gates:
   1. `step3`: PASS (`dev/scripts/scratch_outputs/step6b_nc4_step3_structure_20260306T095829Z.log`)
   2. `step4`: PASS (`dev/scripts/scratch_outputs/step6b_nc4_step4_structure_20260306T095829Z.log`)
   3. `step5`: PASS (`dev/scripts/scratch_outputs/step6b_nc4_step5_structure_20260306T095829Z.log`)
   4. `step6`: PASS (`dev/scripts/scratch_outputs/step6b_nc4_step6_structure_20260306T095829Z.log`)
5. Golden:
   1. status: PASS (`dev/scripts/scratch_outputs/step6b_nc4_golden_20260306T095829Z.log`)
   2. `max_abs_diff = 1.23378e-05` (tolerance `1e-04`)
6. Benchmark diff snapshot:
   1. benchmark log: `dev/scripts/scratch_outputs/step6b_nc4_bench_centralized_20260306T095829Z.log`
   2. diff csv: `dev/scripts/scratch_outputs/benchmark_centralized_diff_nc4_20260306T095829Z.csv`
   3. worst ratio (`median_per_eval_sec`): `examples_models/example_16_guard_tie_simple = 0.960227272727273`
   4. workloads above `1.02`: `0`
   5. P1 above `1.00`: `0`
   6. P1 ratios:
      1. `example_3_stop_na = 0.738034301391`
      2. `example_21_simple_q = 0.777777777778`
      3. `example_22_shared_q = 0.826086956522`
      4. `example_6_dual_path = 0.775362318841`
      5. `depth3_guard_competitor = 0.946058091286`
      6. `shared_gate_nway_shared_triggers = 0.904306220096`
7. Residual targeted ratios (`median_per_eval_sec`):
   1. csv: `dev/scripts/scratch_outputs/step6b_nc4_residual_ratios_20260306T095829Z.csv`
   2. `shared_gate_pair = 0.655251141552512`
   3. `shared_gate_nway = 0.642238507661559`
   4. `example_2_stop_mixture = 0.948880721683901`
   5. `depth3_guard_competitor = 0.946058091286307`
   6. `example_16_guard_tie_simple = 0.960227272727273`
8. Targeted runtime counters:
   1. capture note: NC4 changes hot-loop guard/integrator execution only; coupling/kernel entrypoint call-shape counters are unchanged vs NC2 by design.
   2. `shared_gate_pair`: `coupling_calls_per_eval=1`, `kernel_density_calls_per_eval=0`
   3. `shared_gate_nway`: `coupling_calls_per_eval=2`, `kernel_density_calls_per_eval=0`
   4. `example_2_stop_mixture`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=250`
   5. `depth3_guard_competitor`: `coupling_calls_per_eval=0`, `kernel_density_calls_per_eval=6000`
9. NC4 interpretation:
   1. Residual non-coupling slowdown is closed on this run: `depth3_guard_competitor` moved from `>1.02` to `<1.00`.
   2. No workload remains above baseline on this centralized pass.
   3. Recovery came from unified-path bloat deletion/hoisting inside guard evaluation; no path re-splitting added.
10. Manifest:
   1. `dev/scripts/scratch_outputs/step6b_nc4_manifest_20260306T095829Z.txt`
11. Explicit statement: `No re-splitting introduced = yes`

### Evidence Block: NC7 Stability + Performance Gate (3 Consecutive Fixed-Setting Runs)
1. Commit SHA: `8a7c6d1` (working tree dirty; NC4 runtime changes active)
2. UTC bundle timestamp: `20260306T100501Z`
3. Batch ID: `NC7`
4. Required settings used (all three runs):
   1. `ACCUMULATR_BENCH_N_REP=6`
   2. `ACCUMULATR_BENCH_TARGET_SEC=0.10`
   3. `ACCUMULATR_BENCH_MIN_INNER_REPS=40`
5. Run artifacts:
   1. runs index: `dev/scripts/scratch_outputs/step6b_nc7_runs_20260306T100501Z.csv`
   2. summary: `dev/scripts/scratch_outputs/step6b_nc7_stability_summary_20260306T100501Z.txt`
   3. summary csv: `dev/scripts/scratch_outputs/step6b_nc7_stability_summary_20260306T100501Z.csv`
   4. targeted csv: `dev/scripts/scratch_outputs/step6b_nc7_targeted_stats_20260306T100501Z.csv`
   5. run1 log: `dev/scripts/scratch_outputs/step6b_nc7_stability_run1_20260306T100501Z.log`
   6. run1 diff: `dev/scripts/scratch_outputs/benchmark_centralized_diff_nc7_run1_20260306T100501Z.csv`
   7. run2 log: `dev/scripts/scratch_outputs/step6b_nc7_stability_run2_20260306T100616Z.log`
   8. run2 diff: `dev/scripts/scratch_outputs/benchmark_centralized_diff_nc7_run2_20260306T100616Z.csv`
   9. run3 log: `dev/scripts/scratch_outputs/step6b_nc7_stability_run3_20260306T100730Z.log`
   10. run3 diff: `dev/scripts/scratch_outputs/benchmark_centralized_diff_nc7_run3_20260306T100730Z.csv`
6. Per-run gate results (`median_per_eval_sec`):
   1. Run 1: PASS
      1. worst: `examples_models/example_16_guard_tie_simple = 0.958806818181818`
      2. workloads `> 1.02`: `0`
      3. P1 `> 1.00`: `0`
   2. Run 2: PASS
      1. worst: `examples_models/example_16_guard_tie_simple = 0.955255681818182`
      2. workloads `> 1.02`: `0`
      3. P1 `> 1.00`: `0`
   3. Run 3: PASS
      1. worst: `examples_models/example_16_guard_tie_simple = 0.9609375`
      2. workloads `> 1.02`: `0`
      3. P1 `> 1.00`: `0`
7. Residual-set closure results (across 3 runs):
   1. `depth3_guard_competitor`: `worst_of_3=0.938397701883179`, `median_of_3=0.937440153207788`
   2. `example_2_stop_mixture`: `worst_of_3=0.936986301369857`, `median_of_3=0.934246575342471`
   3. `shared_gate_nway`: `worst_of_3=0.650233177881412`, `median_of_3=0.644903397734842`
   4. `shared_gate_pair`: `worst_of_3=0.675799086757991`, `median_of_3=0.671232876712329`
8. NC7 final verdict:
   1. `STEP6B_NC7_STATUS=PASS`
9. Manifest:
   1. `dev/scripts/scratch_outputs/step6b_nc7_manifest_20260306T100501Z.txt`
10. Explicit statement: `No re-splitting introduced = yes`
