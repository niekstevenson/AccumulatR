# Phase 3 + 6 Hard-Cutover Execution Plan (No Legacy Fallback)

## Goals
1. Complete Phase 3 (`bitset-only forced/guard/trigger state`).
2. Complete Phase 6 (`single fused kernel execution framework`).
3. Preserve correctness (`max |delta loglik| <= 1e-4`) while improving speed.

## Invariants (Must Hold)
1. No `std::unordered_set<int>` forced-state operations in hot evaluator paths.
2. No recursive/scenario fallback in hot loops for observed/ranked likelihood.
3. No string lookups in integrand loops.
4. One runtime kernel entrypoint across observed/ranked/nonresponse internals.

## Current Gaps
1. Mixed forced-state plumbing still present in evaluator signatures (`unordered_set + bitset`).
2. Guard execution still re-enters general recursive integration paths.
3. Competitor and ranked flows still keep partial fallback branching.
4. Vectorized distribution API exists but is only used in a narrow onset branch.

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
1. Implement Step 1 end-to-end in evaluator internals (bitset-only forced state).
2. Then perform first Step 2 cutover by adding a single internal kernel execution entrypoint and routing observed path through it.
