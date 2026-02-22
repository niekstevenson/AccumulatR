# Phase 3 + 6 Hard-Cutover Execution Plan (No Legacy Fallback)

## Goals
1. Complete Phase 3 (`bitset-only forced/guard/trigger state`).
2. Complete Phase 6 (`single fused kernel execution framework`).
3. Preserve correctness (`max |delta loglik| <= 1e-4`) while improving speed.

## Status Snapshot (2026-02-22)
1. Step 1 (bitset-only forced state): **Mostly complete** in hot evaluator paths; set-based forced state is no longer used in `src/distributions.cpp` hot loops.
2. Step 2 (unified kernel entry): **Complete**; observed/ranked/nonresponse now execute through the shared sequence/mass kernel contribution surface, and duplicated sequence-density entry wrappers were consolidated.
3. Step 3 (guard/competitor transition compilation): **Partially complete**; dense guard transition masks are used, but competitor/guard special cases still branch by structure.
4. Step 4 (ranked cutover): **Mostly complete**; ranked state keys are hash-based (non-string hot key path), but residual ranked-special flow still exists.
5. Step 5 (distribution vectorization rollout): **Partially complete**; vectorized pdf/cdf is used in the linear-chain ODE fast path but not fully generalized.
6. Step 6 (legacy/fallback deletion): **Not complete**; structural fallbacks still exist for non-linear guard graphs and specialized paths.

## Invariants (Must Hold)
1. No `std::unordered_set<int>` forced-state operations in hot evaluator paths.
2. No recursive/scenario fallback in hot loops for observed/ranked likelihood.
3. No string lookups in integrand loops.
4. One runtime kernel entrypoint across observed/ranked/nonresponse internals.

## Current Gaps
1. Non-linear guard graphs still fall back to recursive quadrature paths.
2. Kernel execution is shared broadly; remaining divergence is primarily structural fallback logic for non-linear guard topologies (Step 7/Step 6 scope).
3. Competitor/ranked flows still retain structural fallback logic in inner execution routes.
4. Vectorized distribution primitives are deployed in selected fast paths but not fully rolled out across all integration hotspots.

## Execution Steps

### Step 1: Bitset-Only Forced State Core
1. Introduce a bitset-only forced-state carrier for evaluator internals.
2. Remove `unordered_set<int>` forced parameters from:
   - `NodeEvalState`
   - `GuardEvalInput`
   - `eval_node_with_forced_dense*`
   - `eval_event_ref_idx` hot-call chain
3. Keep boundary conversion only at R API edges; internal execution is bitset-only.
4. Remove set-based helpers from hot paths:
   - `forced_contains_scoped(...)` set fallback
   - `forced_set_intersects_scope(...)` set fallback

Acceptance:
1. `rg "unordered_set<int>" src/distributions.cpp` only matches setup/boundary code, not hot kernel loops.
2. Kernel event/guard callbacks run with bitsets only.

### Step 2: Unified Kernel Entry for Observed + Ranked
1. Add one internal execution entrypoint:
   - takes target node/output indices, time, SoA params, bitset state, and trigger state.
2. Route:
   - observed outcome density/probability
   - ranked prefix transitions
   - nonresponse probability accumulation
   through the same kernel executor surface.
3. Remove duplicated per-path runtime evaluators where equivalent.

Acceptance:
1. `cpp_loglik` and `cpp_loglik_multiple` hit one shared internal kernel execution API.
2. No path-specific recursive evaluator call in hot loops.

### Step 3: Guard/Competitor Transition Compilation
1. Compile guard transition effects into dense bitset transition ops.
2. Replace runtime guard source collection + set updates with precompiled transition masks.
3. Competitor survival path uses compiled transitions and dense node slot reuse.

Acceptance:
1. Guard-heavy workloads avoid general fallback recursion.
2. Trigger flips only invalidate dirty slots; no full rebuild of forced state objects.

### Step 4: Ranked Cutover Completion
1. Keep ranked state in dense bitset + compact map keys (no string key serialization in hot merge loops).
2. Reuse same kernel state transition machinery as observed path.
3. Remove residual ranked-special fallback logic from inner loop.

Acceptance:
1. Ranked path uses same fused kernel machinery and bitset transitions.
2. Ranked state collapse key generation is non-string in hot path.

### Step 5: Distribution Vectorization Rollout
1. Expand vectorized distribution primitives (`eval_pdf_vec/eval_cdf_vec`) to guard/competitor integration nodes.
2. Reuse batched node times across outputs/components for each trial.
3. Keep dispatch distribution-general (lognormal/gamma/exgauss under same primitive interface).

Acceptance:
1. Scalar repeated pdf/cdf calls in RK/integration loops are reduced.
2. No distribution-specific branching outside primitive dispatch layer.

### Step 6: Legacy/Fallback Deletion
1. Remove deprecated fallback branches once parity passes.
2. Keep exactly one internal execution engine path.
3. Prune dead helpers and compatibility shims.

Acceptance:
1. Structural checks pass:
   - no recursive fallback in hot path
   - no set-based forced-state in integrands
   - no string-map lookup in kernel loops

### Step 7: Generalized Guard-DAG ODE Consolidation (New)
1. Extend guard solver from linear chains to general acyclic guard DAGs (branched guard compositions).
2. Compile guard DAG dependencies into dense solver state indices and transition operators.
3. Route non-linear guard CDF/survival through the guard-DAG solver instead of quadrature recursion.
4. Keep quadrature only as rollout safety during validation; remove once parity and stability targets pass.

Acceptance:
1. Non-linear guard graphs (branching blocker/reference structure) use the guard-DAG solver path.
2. Equivalence on guard-heavy suite meets `max |delta loglik| <= 1e-4` vs quadrature baseline.
3. After cutover, guard CDF evaluation no longer depends on general recursive quadrature for supported DAG classes.

## Verification Loop (Run After Each Major Step)
1. Build:
   - `R CMD INSTALL .`
2. Correctness:
   - `Rscript dev/scripts/check_loglik_golden.R`
3. Performance:
   - `Rscript dev/scripts/bench_centralized.R`
4. Profile:
   - `ACCUMULATR_PROFILE_RUN_SEC=18 ACCUMULATR_PROFILE_DURATION=18 ACCUMULATR_PROFILE_INTERVAL=1 ACCUMULATR_PROFILE_FILE=dev/scratch_profile_current.txt bash dev/scripts/profile_cpp_simple.sh`

## Immediate Next Execution Slice
1. Execute Step 7 design + implementation for guard-DAG solver on branching guard topologies.
2. Validate Step 7 parity/perf against existing guard-heavy examples and golden checks.
3. Use Step 7 outcomes to drive Step 6 fallback deletion decisions.
