# Full Plan 3 Sub-Plan: Step 2 Unified Kernel Entry (Observed + Ranked + Nonresponse)

Date: 2026-03-06
Status: COMPLETE (`FP2.0` + `FP2.1` + `FP2.2` + `FP2.3` + `FP2.4` + `FP2.5`
complete; closure decision `ACHIEVED`)
Parent: `dev/full_plan3.MD`
Batch ID range: `FP2.*`

## 1) Objective

Route observed, ranked, and nonresponse evaluation through one shared internal
kernel execution surface.

## 2) Current Anchors (Where Work Applies)

1. Trial kernel probability entry:
   1. `evaluate_trial_probability_kernel_idx` at `src/distributions.cpp:8809`
2. Log-likelihood callers:
   1. `cpp_loglik` at `src/distributions.cpp:8953`
   2. `cpp_loglik_multiple` at `src/distributions.cpp:9543`
3. Outcome/node density entries:
   1. `native_outcome_probability_bits_impl_idx` at
      `src/distributions.cpp:6410`
   2. `node_density_entry_idx` at `src/distributions.cpp:6337`
4. Ranked transition machinery entry:
   1. `for_each_sequence_node_transition` at `src/distributions.cpp:7792`

## 3) Scope Boundaries

In scope:

1. Internal evaluator architecture and execution interfaces.
2. Mode normalization for observed/ranked/nonresponse.
3. Hot-loop call graph simplification.

Out of scope:

1. User-facing R API shape.
2. Statistical model semantics.

## 4) Execution Batches

### FP2.0 Baseline + Call Graph Snapshot

1. Freeze current call graph from:
   1. `cpp_loglik`
   2. `cpp_loglik_multiple`
   3. `evaluate_trial_probability_kernel_idx`
2. Identify path-specific execution forks by mode.

### FP2.1 Define One Internal Kernel Execution Contract

1. Introduce one internal execution interface that takes:
   1. target outcome/node references
   2. time inputs
   3. trial/component params
   4. forced/trigger state carrier
   5. execution intent (density/cdf/survival/mass)
2. Ensure mode differences are normalized before entering this interface.

### FP2.2 Route Observed Path Through Unified Entry

1. Replace observed-specific hot evaluator branch calls with unified interface
   invocation.
2. Preserve existing correctness behavior.

### FP2.3 Route Ranked Path Through Same Entry

1. Ensure ranked transition eval uses same execution entry for node/cdf work.
2. Keep ranked-specific sequencing logic only as planner/input shaping, not
   separate evaluator engine.

### FP2.4 Route Nonresponse Path Through Same Entry

1. Ensure nonresponse path uses same unified evaluator entry and same runtime
   state update machinery.
2. Remove duplicated runtime evaluators where behavior is equivalent.

### FP2.5 Verification + Evidence

1. Run full verification loop.
2. Record mode-to-entrypoint mapping proof.

## 5) Acceptance (What Must Be True)

1. `cpp_loglik` and `cpp_loglik_multiple` dispatch to one shared internal
   kernel execution API for observed/ranked/nonresponse internals.
2. No path-specific recursive evaluator call remains in hot loops for these
   modes.
3. Ranked/nonresponse mode-specific code is limited to normalization/planning,
   not separate evaluator engines.
4. Golden passes with `max_abs_diff <= 1e-4`.
5. Performance gates remain within thresholds.

## 6) What Counts

1. Real call graph convergence to one runtime entry.
2. Deletion of duplicated evaluator branches.
3. Evidence showing mode execution convergence.

## 7) What Does Not Count

1. Wrapper function that still branches into old per-mode evaluators.
2. Separate per-mode loops with same function name.
3. Completion claim without call graph proof and benchmark evidence.

## 8) Hard Stops

1. Any new mode-branch in hot execution for speed reasons.
2. Any correctness drift beyond tolerance.
3. Any nonresponse/ranked regressions hidden by selective workload reporting.

## 9) Required Evidence Bundle (`FP2.*`)

1. `full_plan3_FP2_*_callgraph.txt`
2. `full_plan3_FP2_*_mode_entryproof.txt`
3. `full_plan3_FP2_*_step3_structure.log`
4. `full_plan3_FP2_*_step4_structure.log`
5. `full_plan3_FP2_*_step5_structure.log`
6. `full_plan3_FP2_*_step6_structure.log`
7. `full_plan3_FP2_*_golden.log`
8. `full_plan3_FP2_*_bench.log`
9. `full_plan3_FP2_*_diff.csv`

## 10) Exit Rule

Step 2 is complete only when Sections 5 and 9 are satisfied on one commit
lineage.

## 11) Evidence Block: FP2.0 Baseline + Call Graph Snapshot

1. Timestamp (UTC): `20260306T134122Z`
2. Commit SHA: `54437eb`
3. Batch: `FP2`, subbatch `FP2.0`
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
   1. worst ratio: `0.9422348` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP2.0 call-graph freeze facts:
   1. All trial likelihood evaluations pass through one runtime kernel entry
      callsite: `evaluate_trial_probability_kernel_idx` (called from
      `cpp_loglik`).
   2. Mode split points still active inside trial kernel execution:
      `eval_input.nonresponse_mode`,
      `spec.metric == KernelContributionMetric::SequenceDensity`,
      `spec.sequence_length == 1`, and `spec.guess_shortcut`.
   3. Ranked/nonresponse normalization occurs before entry in
      `normalize_trial_eval_input_for_trial`.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T134122Z_summary.txt`

## 12) Evidence Block: FP2.1 Unified Internal Execution Contract

1. Timestamp (UTC): `20260306T140655Z`
2. Commit SHA: `54437eb`
3. Batch: `FP2`, subbatch `FP2.1`
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
   1. worst ratio: `0.9382716` (`example_23_ranked_chain`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP2.1 contract proof:
   1. Contribution execution now routes through one intent-driven contract:
      `TrialKernelIntent` + `evaluate_trial_contribution_kernel_idx`.
   2. Trial result transformation is normalized by
      `TrialProbabilityTransform`; legacy `eval_input.nonresponse_mode`
      predicate is removed.
   3. Deleted split predicates from FP2.0 inventory are absent:
      `spec.metric == KernelContributionMetric::SequenceDensity`,
      `spec.sequence_length == 1`,
      `spec.guess_shortcut`,
      `eval_input.nonresponse_mode`.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T140655Z_summary.txt`

## 13) Evidence Block: FP2.2 Observed Path Cutover On Unified Contract

1. Timestamp (UTC): `20260306T150007Z`
2. Commit SHA: `54437eb`
3. Batch: `FP2`, subbatch `FP2.2`
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
   1. worst ratio: `0.9822920` (`example_2_stop_mixture`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP2.2 cutover proof:
   1. Observed-path contribution intents are now collapsed to one
      `TrialKernelIntent::SequenceDensity` intent (single + multi-step).
   2. Legacy observed split intents (`NodeDensity`, `GuessDensity`,
      `SequencePrefixDensity`) are deleted.
   3. `evaluate_sequence_density_kernel_idx` handles single-step and
      multi-step sequence execution inside one observed evaluator surface.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T150007Z_summary.txt`

## 14) Evidence Block: FP2.3 Ranked Path Cutover On Unified Contract

1. Timestamp (UTC): `20260306T151253Z`
2. Commit SHA: `54437eb`
3. Batch: `FP2`, subbatch `FP2.3`
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
7. FP2.3 cutover proof:
   1. Ranked execution no longer depends on a ranked-only execution toggle in
      `TrialEvalInput`; sequence behavior is driven by `sequence_length`.
   2. Ranked and observed sequence execution share one evaluator surface:
      `evaluate_sequence_density_kernel_idx` with
      `sequence_length == 1` / `sequence_prefix_density_resolved`.
   3. Legacy ranked execution predicate removed:
      `eval_input.enforce_unique_sequence_nodes`.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151253Z_summary.txt`

## 15) Evidence Block: FP2.4 Nonresponse Path Cutover Cleanup

1. Timestamp (UTC): `20260306T151939Z`
2. Commit SHA: `54437eb`
3. Batch: `FP2`, subbatch `FP2.4`
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
   1. worst ratio: `0.9657534` (`example_2_stop_mixture`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP2.4 cutover proof:
   1. Runtime nonresponse mode flag branches removed from execution path
      (`TrialEvalMode`, `eval_input.mode`, `eval_input.nonresponse_mode`).
   2. Nonresponse now executes through normalized input shape
      (`sequence_length == 0`) + `TrialProbabilityTransform::Complement`.
   3. Caller-side nonresponse program pointer binding is unconditional, not
      mode-branch guarded.
8. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
9. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step3_structure.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step4_structure.log`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step5_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_step6_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_golden.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_bench.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_benchmark_centralized.csv`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_benchmark_centralized_diff.csv`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_benchmark_step3_comparison.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T151939Z_summary.txt`

## 16) Evidence Block: FP2.5 Step 2 Closure Decision + Deletion-Proof Rollup

1. Timestamp (UTC): `20260306T152414Z`
2. Commit SHA: `54437eb`
3. Batch: `FP2`, subbatch `FP2.5`
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
   1. worst ratio: `0.9602273` (`example_16_guard_tie_simple`)
   2. count above `1.02`: `0`
   3. P1 above `1.00`: `0`
7. FP2.5 closure proof:
   1. `cpp_loglik` and `cpp_loglik_multiple` both feed one trial kernel entry:
      `evaluate_trial_probability_kernel_idx`.
   2. Observed/ranked/nonresponse route through one contribution execution
      surface: `evaluate_trial_contribution_kernel_idx` with
      `TrialKernelIntent`.
   3. Mode differences remain in normalization/planning only
      (`TrialProbabilityTransform` + normalized input shape); no parallel
      mode-specific evaluator engine remains reachable.
8. Closure decision:
   1. `decision = ACHIEVED` (Step 2 acceptance met).
   2. Residual split points are intent/shape branches inside the same unified
      engine and are forwarded to FP3+ contraction work.
9. Wrapper-only consolidation check:
   1. `No wrapper-only consolidation introduced = yes`
10. Artifacts:
   1. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_status.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_callgraph.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_mode_entryproof.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step2_closure_probe.log`
   7. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step2_closure_decision.txt`
   8. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step3_structure.log`
   9. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step4_structure.log`
   10. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step5_structure.log`
   11. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_step6_structure.log`
   12. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_golden.log`
   13. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_bench.log`
   14. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_benchmark_centralized.csv`
   15. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_benchmark_centralized_diff.csv`
   16. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_benchmark_step3_comparison.csv`
   17. `dev/scripts/scratch_outputs/full_plan3_FP2_20260306T152414Z_summary.txt`
