# Full Plan 3 Step 5 Artifact: Kernel-Level Batch Cutover

Date: 2026-03-09
Status: ACHIEVED (`FP5K.0` through `FP5K.6` complete; Step 6 next)
Parent: `dev/full_plan3.MD`
Related umbrella plan: `dev/full_plan3_step5_6_vectorization_legacy_deletion.MD`
Batch ID range: `FP5K.*`

## 1) Brutal Current State

1. Step 5 is achieved.
2. The current code has vectorized coupling evaluation in parts of the unified
   coupling stack, but the non-coupling guard/competitor path is still driven
   by per-time kernel evaluation.
3. The first `FP5.1` attempt was rejected because it moved batching logic into
   `src/distributions.cpp` helper code instead of making the kernel executor
   itself batch-native.
4. Honest summary: Step 5 will fail again if work keeps happening at the
   caller/helper layer instead of the executor layer.
5. `FP5K.2` fixed the direct generic `NodeRef` no-competitor boundary:
   that runtime now calls the kernel executor contract.
6. `FP5K.3` moved competitor-bearing generic `NodeRef` execution onto the
   same batch executor boundary.
7. `FP5K.4` deleted the batch-only compiled competitor helper engine from
   `src/distributions.cpp`.
8. `FP5K.5` stabilized the executor-side compiled competitor-op loop without
   re-splitting the runtime boundary.
9. `FP5K.6` closed Step 5 on one commit lineage.
10. The remaining split surface is Step 6 debt, not Step 5 debt:
    specialized coupling dispatch and provider-mode structural branches still
    exist and must be deleted under `FP6.*`.

## 2) Rejected Boundary (Must Not Be Repeated)

The following does not count as Step 5 progress:

1. adding cache objects or precompute tables inside
   `node_density_with_competitors_batch_internal`
2. wrapping repeated scalar calls in a new helper while
   `uuber::eval_kernel_node_incremental` remains the true per-time engine
3. claiming "one path" while the heavy work is still a scalar loop hidden under
   a batch wrapper

Why this is rejected:

1. it is pseudo-consolidation
2. it preserves the real old execution boundary
3. it already produced a `1.925x` regression on a direct guarded generic
   `NodeRef` target case in `FP5.1`

## 3) Required Boundary (What Real Consolidation Means Here)

Step 5 only counts if the batch-native cut happens at the kernel execution
layer.

Required runtime boundary:

1. `src/kernel_executor.cpp`
2. `src/kernel_executor.h`

Required caller shape after cutover:

1. `src/distributions.cpp:2143`
   `node_density_with_competitors_batch_internal` becomes a thin coordinator,
   not a hidden scalar engine
2. `src/distributions.cpp:3766`
   `evaluate_coupling_acc_batch` remains on the shared batch-native stack
3. `src/distributions.cpp:3873`
   `evaluate_labelref_batch` remains on the shared batch-native stack
4. `src/distributions.cpp:4261`
   `evaluate_outcome_coupling_unified` remains the unified coupling entry

Required executor property after cutover:

1. `NodeRef` and competitor/guard evaluation must use the same batch-native
   kernel executor contract
2. provider differences may remain as data/provider mode, not as alternate
   evaluator engines

## 4) Scope

In scope:

1. kernel-level batch execution for repeated time evaluation
2. direct guarded generic `NodeRef` target cases
3. competitor and guard reuse inside the same executor
4. deletion of helper-local scalar orchestration that becomes dead afterward

Out of scope:

1. new "fast path" engines
2. workload-specific alternate runtimes
3. helper-local caching in `src/distributions.cpp` counted as Step 5 progress

## 5) Hard Sequence

### FP5K.0 Target Freeze + Harness Lock

1. Freeze the exact direct target cases that Step 5 must satisfy.
2. Add repo-owned harnesses for at least these surfaces:
   1. direct guarded generic `NodeRef`
   2. direct guarded nested `NodeRef`
   3. competitor-bearing generic `NodeRef`
3. Record baseline values for:
   1. probability
   2. `generic_noderef_batch_calls`
   3. `generic_scalar_fallback_calls`
   4. `median_per_eval_sec`

Acceptance:

1. harnesses live in the repo, not `/tmp`
2. direct target cases are reproducible on current `HEAD`
3. baseline evidence is frozen

### FP5K.1 Define One Batch-Native Kernel Executor Contract

1. Add one explicit batch executor API in `src/kernel_executor.h`
2. Implement it in `src/kernel_executor.cpp`
3. Keep scalar `eval_kernel_node_incremental` only as a lower primitive during
   transition if still required internally
4. Do not add a second public evaluator surface in `src/distributions.cpp`

Acceptance:

1. one executor contract exists for vector-of-times evaluation
2. caller code can request batched node evaluation without per-time callback
   orchestration
3. no helper-local "mini executor" is introduced in `src/distributions.cpp`

### FP5K.2 Port Direct Generic `NodeRef` To The Kernel Executor

1. Route direct generic `NodeRef` target cases through the new executor
   contract
2. Keep provider behavior inside the same executor stack
3. Prove:
   1. `generic_noderef_batch_calls >= 1`
   2. `generic_scalar_fallback_calls = 0`

Acceptance:

1. direct target cases use the kernel-level batch executor
2. scalar fallback is not used on those cases
3. probability parity holds

### FP5K.3 Port Guard/Competitor Reuse Into The Same Executor

1. Move repeated guard/event/competitor time evaluation reuse into the kernel
   executor layer
2. Delete or thin duplicated orchestration in
   `node_density_with_competitors_batch_internal`
3. Keep one executor stack for direct and competitor-bearing `NodeRef` cases

Acceptance:

1. repeated per-time guard/competitor work is batch-native in the executor
2. caller-side orchestration shrinks or is deleted
3. direct and competitor-bearing target cases remain on one executor contract

### FP5K.4 Delete Helper-Local Scalar Orchestration

1. Remove batch-only cache/state helpers in `src/distributions.cpp` that would
   be redundant after `FP5K.2` and `FP5K.3`
2. Make `node_density_with_competitors_batch_internal` a thin entry wrapper or
   merge it away if feasible
3. Produce explicit deletion proof

Acceptance:

1. dead helper-local orchestration is deleted
2. no alternate helper engine remains reachable
3. split inventory for Step 5 is strictly smaller than at `FP5.0`

### FP5K.5 Stabilize Inside The Same Executor

1. Recover any speed loss inside the kernel executor only
2. Allow:
   1. scratch reuse
   2. invariant hoisting
   3. compact need-mask/predecoded provider data
3. Do not allow:
   1. new mode-specific engine
   2. special direct-case evaluator outside the executor

Acceptance:

1. direct target cases are not slower than the frozen baseline linearly enough
   to block closure
2. centralized benchmark remains within global gates
3. no re-split appears

### FP5K.6 Step 5 Closure

1. Run full gate bundle
2. Run direct target harnesses
3. Freeze vector coverage, split inventory, and deletion proof

Acceptance:

1. Step 5 acceptance in the umbrella artifact is satisfied
2. direct target cases prove real batch-native execution
3. closure bundle is complete

## 6) Mandatory Target Surfaces

The standard centralized suite is not enough.

The following surfaces must be covered explicitly:

1. direct guarded generic `NodeRef`
2. direct guarded nested `NodeRef`
3. competitor-bearing generic `NodeRef`

For each surface, the frozen evidence must include:

1. probability result
2. `generic_noderef_batch_calls`
3. `generic_scalar_fallback_calls`
4. `median_per_eval_sec`
5. baseline ratio versus frozen pre-cut commit

## 7) Step 5 Acceptance Under This Artifact

Step 5 is achieved only when all are true:

1. repeated per-time guard/competitor work is batch-native at the kernel
   executor layer
2. direct target cases prove batch-native `NodeRef` execution with
   `generic_scalar_fallback_calls = 0`
3. helper-local scalar orchestration is deleted or reduced to thin wiring
4. `max_abs_diff <= 1e-4`
5. centralized benchmark gates pass:
   1. no workload ratio > `1.02`
   2. no P1 workload ratio > `1.00`
6. direct target cases are not regressed versus the frozen baseline badly
   enough to contradict the cutover claim

## 8) What Counts

1. a kernel-level batch executor that real workloads actually use
2. direct `NodeRef` evidence proving that the batch executor is exercised
3. deletion of helper-local orchestration after cutover
4. same-engine speed recovery inside `src/kernel_executor.cpp`

## 9) What Does Not Count

1. cache layers inside `src/distributions.cpp`
2. batch wrappers around scalar kernel execution
3. green example/nested suite results when direct `NodeRef` cases are not
   exercised
4. leaving helper-local scalar orchestration in place and calling it
   "consolidated"

## 10) Hard Stops

1. any new runtime split added for speed
2. any direct target case using scalar fallback after cutover
3. any claim of Step 5 progress without repo-owned target harnesses
4. any batch that is structurally cleaner but still keeps the old execution
   boundary in practice

## 11) Evidence Contract

Each `FP5K.*` batch must produce:

1. `full_plan3_FP5K_*_summary.txt`
2. `full_plan3_FP5K_*_manifest.txt`
3. `full_plan3_FP5K_*_target_cases.txt`
4. `full_plan3_FP5K_*_vector_coverage.txt`
5. `full_plan3_FP5K_*_split_inventory.txt`
6. `full_plan3_FP5K_*_deleted_symbols.txt` if deletion occurs
7. `full_plan3_FP5K_*_step5_structure.log`
8. `full_plan3_FP5K_*_golden.log`
9. `full_plan3_FP5K_*_bench.log`
10. `full_plan3_FP5K_*_diff.csv`

## 12) Control Rule

For Step 5 execution, this artifact is the controlling interpretation of
`FP5.1` and its follow-on batches.

If this artifact conflicts with broader Step 5 wording elsewhere:

1. this artifact wins for implementation boundary
2. the umbrella Step 5-6 artifact still controls final Step 5/6 closure gates

## 13) FP5K.0 Evidence Freeze (`20260309T152055Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `no`
3. What landed:
   1. repo-owned target harness:
      `dev/scripts/bench_step5_kernel_target_cases.R`
   2. frozen target surfaces:
      1. `pair_guard_generic_noderef`
      2. `nested_guard_generic_noderef`
      3. `competitor_guard_generic_noderef`
4. Gate results:
   1. `install`: pass
   2. `step5_structure`: pass
   3. `golden`: pass
   4. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark freeze:
   1. worst ratio: `1.007948588`
   2. P1 above `1.00`: `0`
7. Direct target-case freeze:
   1. `pair_guard_generic_noderef`
      1. surface: `direct_guarded_generic_noderef`
      2. `competitor_count = 0`
      3. `probability = 0.574695700275553`
      4. `median_per_eval_sec = 8.33333333333334e-06`
      5. `generic_noderef_batch_per_eval = 1`
      6. `generic_scalar_fallback_per_eval = 0`
   2. `nested_guard_generic_noderef`
      1. surface: `direct_guarded_nested_noderef`
      2. `competitor_count = 0`
      3. `probability = 0.679796002185848`
      4. `median_per_eval_sec = 1.5e-05`
      5. `generic_noderef_batch_per_eval = 1`
      6. `generic_scalar_fallback_per_eval = 0`
   3. `competitor_guard_generic_noderef`
      1. surface: `competitor_bearing_generic_noderef`
      2. `competitor_count = 2`
      3. `probability = 0.0105555139919877`
      4. `median_per_eval_sec = 0.001355`
      5. `generic_noderef_batch_per_eval = 1`
      6. `generic_scalar_fallback_per_eval = 0`
8. Brutal status:
   1. `FP5K.0` is achieved.
   2. Step 5 is still not achieved.
   3. The direct target baseline is now frozen, so future Step 5 work can no
      longer hide behind standard-suite green status alone.
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152055Z_split_inventory.txt`

## 14) FP5K.1 Evidence Freeze (`20260309T152855Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. kernel-level batch contract in `src/kernel_executor.h`
   2. kernel-level batch executor implementation in `src/kernel_executor.cpp`
   3. Step 5 structure gate now requires that contract
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.979066023`
   2. P1 above `1.00`: `0`
7. Direct target-case comparison versus `FP5K.0`:
   1. `pair_guard_generic_noderef`
      1. ratio current / `FP5K.0` = `1.0`
      2. `generic_noderef_batch_per_eval = 1`
      3. `generic_scalar_fallback_per_eval = 0`
   2. `nested_guard_generic_noderef`
      1. ratio current / `FP5K.0` = `1.0`
      2. `generic_noderef_batch_per_eval = 1`
      3. `generic_scalar_fallback_per_eval = 0`
   3. `competitor_guard_generic_noderef`
      1. ratio current / `FP5K.0` = `0.9876999`
      2. `generic_noderef_batch_per_eval = 1`
      3. `generic_scalar_fallback_per_eval = 0`
8. Honest status:
   1. `FP5K.1` is achieved.
   2. This is contract consolidation, not runtime cutover.
   3. The direct `NodeRef` runtime is still not using the new batch executor.
   4. `FP5K.2` is the first batch that must move real execution to the new
      contract.
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T152855Z_split_inventory.txt`

## 15) FP5K.2 Direct Generic `NodeRef` Runtime Port (`20260309T153954Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. direct generic `NodeRef` no-competitor runtime in
      `evaluate_outcome_coupling_unified` now routes through
      `eval_node_with_forced_state_view_batch`
   2. `eval_node_with_forced_state_view_batch` now calls
      `eval_kernel_node_batch_incremental` in `src/kernel_executor.cpp`
   3. Step 5 structure gate now rejects the old direct scalar helper call in
      that unified entrypoint
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `1`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
7. Direct target-case comparison versus `FP5K.0`:
   1. `pair_guard_generic_noderef`: ratio `1`
   2. `nested_guard_generic_noderef`: ratio `1`
   3. `competitor_guard_generic_noderef`: ratio `0.9852399`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
8. Honest status:
   1. `FP5K.2` is achieved
   2. this is the first real runtime cutover onto the kernel executor contract
   3. only direct generic `NodeRef` no-competitor runtime is ported in this
      batch
   4. competitor-bearing `NodeRef` orchestration still remains in
      `node_density_with_competitors_batch_internal`
   5. `FP5K.0` / `FP5K.1` direct target probability values for the two
      no-competitor surfaces were wrong because the old direct branch treated
      `generic.node_id` as a dense node index
   6. `full_plan3_FP5K_20260309T153954Z_direct_port_probe.log` proves the
      corrected resolved-node scalar value matches the new batch value
   7. Step 5 is still not achieved
   8. `FP5K.3` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_deleted_symbols.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T153954Z_direct_port_probe.log`

## 16) FP5K.3 Competitor-Bearing `NodeRef` Batch Cutover (`20260309T160625Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. `node_density_with_competitors_batch_internal` no longer loops
      per-time over scalar base density plus scalar competitor ops
   2. batched base-node density now uses `eval_node_batch_with_state_dense`
   3. batched compiled competitor ops now use
      `run_competitor_compiled_ops_product_batch`
   4. competitor-bearing `NodeRef` execution now exercises
      `eval_kernel_nodes_batch_incremental` in `src/kernel_executor.cpp`
   5. Step 5 structure gate now rejects reintroduction of scalar base-node or
      scalar competitor-op execution in that batch helper
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.963576159`
   2. worst workload: `example_10_exclusion`
   3. P1 above `1.00`: `0`
7. Direct target-case comparison versus `FP5K.2`:
   1. `pair_guard_generic_noderef`: ratio `1`
   2. `nested_guard_generic_noderef`: ratio `1`
   3. `competitor_guard_generic_noderef`: ratio `0.9912609`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
8. Honest status:
   1. `FP5K.3` is achieved
   2. direct and competitor-bearing generic `NodeRef` surfaces now use the
      same batch executor boundary
   3. helper-layer orchestration still exists in
      `node_density_with_competitors_batch_internal`
   4. Step 5 is still not achieved
   5. `FP5K.4` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T160625Z_deleted_symbols.txt`

## 17) FP5K.4 Helper-Orchestration Deletion Cut (`20260309T162200Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. deleted `eval_nodes_batch_with_state_dense`
   2. deleted `run_competitor_batch_survival_op_batch`
   3. deleted `run_competitor_compiled_ops_product_batch`
   4. added executor-side
      `eval_kernel_competitor_product_batch_incremental`
   5. `node_density_with_competitors_batch_internal` now binds state, runs one
      batched base-node density call, and delegates the compiled competitor-op
      batch loop to the executor
   6. Step 5 structure gate now rejects reintroduction of the deleted helper
      batch engine
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.980113636363642`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
7. Direct target-case comparison versus `FP5K.3`:
   1. `pair_guard_generic_noderef`: ratio `1`
   2. `nested_guard_generic_noderef`: ratio `1`
   3. `competitor_guard_generic_noderef`: ratio `1.0025189`
   4. all remain `generic_noderef_batch_per_eval = 1`
   5. all remain `generic_scalar_fallback_per_eval = 0`
8. Honest status:
   1. `FP5K.4` is achieved
   2. the helper-local compiled competitor-op batch engine is deleted
   3. split inventory is smaller than `FP5.0`
   4. Step 5 is still not achieved
   5. `FP5K.5` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_vector_coverage.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_split_inventory.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T162200Z_deleted_symbols.txt`

## 18) FP5K.5 Executor-Only Stabilization (`20260309T165456Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `yes`
3. What landed:
   1. executor-side batch target-slot resolution now reuses runtime slots
      directly instead of materializing copied `KernelNodeBatchValues`
      vectors in the compiled competitor-op loop
   2. no helper-layer code was added
   3. no runtime split was added
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.977272727272733`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
7. Direct target-case comparison versus `FP5K.4`:
   1. frozen run:
      1. `pair_guard_generic_noderef`: ratio `1`
      2. `nested_guard_generic_noderef`: ratio `1.1111111`
      3. `competitor_guard_generic_noderef`: ratio `1.0018832`
   2. supplemental rerun:
      1. `pair_guard_generic_noderef`: ratio `1`
      2. `nested_guard_generic_noderef`: ratio `1`
      3. `competitor_guard_generic_noderef`: ratio `0.9987437`
   3. all remain `generic_noderef_batch_per_eval = 1`
   4. all remain `generic_scalar_fallback_per_eval = 0`
8. Honest status:
   1. `FP5K.5` is achieved
   2. there is no meaningful speed regression signal after the supplemental
      rerun; the one frozen nested-case blip was timer-quantum noise
   3. Step 5 is still not achieved
   4. `FP5K.6` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_target_cases_rerun.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_vector_coverage.txt`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T165456Z_split_inventory.txt`

## 19) FP5K.6 Step 5 Closure (`20260309T171110Z`)

1. Commit: `1f7f158`
2. Runtime code changed: `no`
3. What landed:
   1. formal closure bundle on the stabilized Step 5 executor boundary
   2. direct target harness rerun included in the closure bundle to rule out
      timer-quantum noise
   3. current hot-loop probe, deletion proof, and closure decision frozen
4. Gates:
   1. `install`: pass
   2. `step3_structure`: pass
   3. `step4_structure`: pass
   4. `step5_structure`: pass
   5. `step6_structure`: pass
   6. `golden`: pass
   7. `bench`: pass
   8. `target_cases`: pass
5. Correctness:
   1. `max_abs_diff = 1.23378e-05`
6. Centralized benchmark:
   1. worst ratio: `0.97727272727275`
   2. worst workload: `example_16_guard_tie_simple`
   3. P1 above `1.00`: `0`
7. Direct target-case comparison:
   1. all three direct target surfaces remain
      `generic_noderef_batch_per_eval = 1`
   2. all three remain
      `generic_scalar_fallback_per_eval = 0`
   3. closure run vs `FP5K.5`:
      1. `pair_guard_generic_noderef`: `1`
      2. `nested_guard_generic_noderef`: `0.9`
      3. `competitor_guard_generic_noderef`: `0.9931034`
8. Closure decision:
   1. `FP5K.6` is achieved
   2. Step 5 is achieved
   3. residual split surface is carried forward explicitly as Step 6 debt
   4. `FP6.0` is the next strict batch
9. Evidence bundle:
   1. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_manifest.txt`
   2. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_summary.txt`
   3. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_target_cases.txt`
   4. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_target_cases_rerun.txt`
   5. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_hotloop_probe.csv`
   6. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_vector_coverage.txt`
   7. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_split_inventory.txt`
   8. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_deleted_symbols.txt`
   9. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_deletion_proof.txt`
   10. `dev/scripts/scratch_outputs/full_plan3_FP5K_20260309T171110Z_step5_closure_decision.txt`
