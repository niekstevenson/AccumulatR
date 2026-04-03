# Scalar To Vector Inventory

Last updated: 2026-03-26

This note inventories the remaining scalar or scalar-shaped execution in the
likelihood framework reachable from `cpp_loglik()` and the closely related
native helpers.

It replaces older inventory notes that referenced deleted files such as
`src/coupling_eval.cpp`. The current production call graph lives mainly in:

- `src/vector_runtime.h`
- `src/vector_program.cpp`
- `src/distributions.cpp`
- `src/exact_outcome_density.cpp`
- `src/ranked_transitions.cpp`
- `src/prep_builder.cpp`

## Bottom Line

- The hot path is no longer predominantly scalar.
- The canonical vector foundation is real:
  - `uuber::VectorProgram`
  - `uuber::VectorLaneBatch`
  - compiled tree vector programs required by `cpp_loglik()`
- `cpp_loglik()` now assigns every eligible trial through grouped batch
  evaluation. It no longer keeps a separate per-trial fallback loop for
  eligible work.
- The remaining non-vector work is mostly scheduler/orchestration:
  - ranked still uses a ranked-specific frontier scheduler
  - exact still uses an exact-specific scenario scheduler
  - outcome-mass, mixture, donor/alias accumulation, and parameter-set helpers
    still keep scalar outer loops
- Honest status:
  - scalar fallback deletion is largely done
  - runtime consolidation is not done

## Scope

Included:

- `src/vector_runtime.h`
- `src/vector_program.cpp`
- `src/distributions.cpp`
- `src/exact_outcome_density.cpp`
- `src/ranked_transitions.cpp`
- `src/prep_builder.cpp`
- helper entry points such as:
  - `cpp_loglik_multiple()`
  - `native_outcome_probability_*()`
  - `native_trial_mixture_*()`

Excluded:

- unrelated scalar utilities outside likelihood evaluation
- deleted historical files mentioned in older plans
- debug/export-only paths unless they mirror production execution

## Terms

- `Fully batched`
  - one execution over many lanes/points without dropping to a separate
    per-point engine
- `Batched inner / scalar outer`
  - leaf math is batched, but outer control flow still loops states, outcomes,
    components, masks, terms, donors, or parameter sets
- `Scalar API surface`
  - helper accepts one upper bound / one trial / one time / one param set even
    if it calls batch kernels internally
- `Scheduler shell`
  - mode-specific frontier/scenario orchestration that sits on top of shared
    vector leaves

## 1. Canonical Vector Foundation Is Real

The framework now has a real vector runtime substrate:

- `uuber::VectorProgram` in `src/vector_runtime.h`
- `uuber::VectorLaneBatch` in `src/vector_runtime.h`
- `compile_tree_vector_program(...)` in `src/vector_program.cpp`
- `cpp_loglik()` requires a compiled tree vector program and compiled tree
  guard metadata
- exact has already adopted the shared lane container:
  - `using ExactScenarioBatch = uuber::VectorLaneBatch`

Status:

- `Implemented`
- This is no longer a scalar framework with a few batch helpers bolted on

## 2. Top-Level cpp_loglik Batching Is Broader Than The Old Note Claimed

`cpp_loglik()` now builds and executes grouped batch plans for:

- direct single-step trials
- grouped direct multi-contribution trials
- ranked sequence trials
- nonresponse / outcome-mass trials
- shared-trigger outcome-mass trials

Important current behavior:

- every eligible trial is expected to be assigned by those grouped batch passes
- if an eligible trial is still unassigned afterwards, `cpp_loglik()` stops
  rather than running a separate per-trial fallback loop

Implication:

- there is no longer a separate "leftover per-trial executor" inside
  `cpp_loglik()`
- singleton work still exists, but it is routed through the same grouped batch
  machinery, usually as one-point batching rather than through a different
  scalar engine

Status:

- `Grouped batch orchestration at the trial level`
- `Not scalar-only at the top level anymore`

## 3. Direct Node Density Is Genuinely Batched Now

The production density entry is `node_density_entry_batch_idx()` in
`src/distributions.cpp`.

Current behavior:

- if per-point `TrialParamsSoA` is already available, it evaluates directly
  through the batch path
- shared-trigger density expands `(mask x point)` batches, optionally
  compresses duplicate `(time, TrialParamsSoA*)` points, evaluates once, then
  reduces by mask weights
- competitor-bearing density first tries specialized overlap batch kernels when
  possible
- if overlap requires exact scenario handling, it routes through
  `exact_outcome_density_batch_from_state(...)`
- otherwise it uses batched tree-node evaluation plus compiled competitor
  survival products
- unsupported shared-trigger density cases return `false`; they do not silently
  re-enter an older scalar density engine

Direct fastpaths that are now batched:

- simple accumulator-event density
- simple competing-event density
- shared-trigger density across mask-expanded SoA points

Status:

- `Batched hot path with multiple vectorized specializations`
- `No separate scalar node-density engine remains in production`

## 4. Outcome Probability Integration Is Mixed

There are now two distinct shapes.

### 4.1 Scalar Helper Still Exists

`integrate_outcome_probability_from_density_idx()`:

- takes one upper bound
- builds a one-point batch callback around `node_density_entry_batch_idx()`
- integrates with `integrate_boost_fn_0_to_upper(...)`

Status:

- `Scalar API surface with batched leaf calls`

### 4.2 Batch Helper Exists And Is Useful

`integrate_outcome_probability_from_density_batch_idx()`:

- accepts many trial-param points
- canonicalizes / deduplicates per-point `TrialParamsSoA`
- builds batched time-query plans
- evaluates density for all active query points in one batch call per
  segment/window
- reuses workspace buffers when supplied

Finite upper bounds:

- one batched interval plan over `[0, upper]`

Infinite upper bounds:

- sequential window growth (`[0, 0.5]`, `[0.5, 1]`, ...)
- each window is batched internally, but the window scheduler itself is scalar

Status:

- `Batched inner evaluation`
- `Scalar outer scheduling still exists for infinite-tail probability integration`

## 5. Ranked Is Batched, But Still Has A Ranked-Specific Scheduler Shell

Production ranked work now routes through
`sequence_prefix_density_batch_resolved(...)` in `src/ranked_transitions.cpp`.
The old scalar ranked oracle is gone from the hot path.

What the current ranked executor does:

- builds a frontier of ranked states
- scans active states
- optionally sorts/group-slices frontier states by reusable keys
- materializes lane groups
- evaluates denom/scenario/competitor work in batch for each group
- collapses next-state candidates and repeats by rank

What is improved:

- ranked no longer falls back to a separate scalar engine
- homogeneous groups can carry a uniform `TrialParamsSoA*` down into exact
  evaluation
- shared-trigger ranked groups can batch `(mask x trial)` work in `cpp_loglik()`

What is still scalar-shaped:

- frontier scanning
- frontier sorting/group scheduling
- scenario collection and state-delta application
- next-state collapse / reduce between ranks

Key ranked-specific runtime types that still exist:

- `RankedProgramSliceRef`
- `RankedBatchPlan`
- `RankedBatchPlanner`

Status:

- `Batched leaves under a ranked scheduler shell`
- `Not yet the canonical vector runtime end-to-end`

## 6. Exact Uses The Shared Lane Batch, But Still Has An Exact-Specific Scheduler Shell

The exact path now shares the lane container and no longer has a separate
scalar exact engine in production.

Current exact behavior:

- `exact_collect_scenarios_batch_impl(...)` seeds scenario points from the
  input times into an `ExactScenarioBatch`
- exact scenario collection is then handled through exact-specific
  aligned/deterministic scenario builders
- `exact_outcome_density_batch_from_state(...)` accumulates scenario weights
  back into output density
- exact competitor survival is batched and invariant-checked; failures are not
  recovered by a scalar fallback path

What is already good:

- exact scenario batches use the shared lane container
- exact competitor survival is batched
- there is no scalar competitor recovery path in production
- singleton exact work is routed through the batched engine

What is still scalar-shaped:

- seeding scenario points one input time at a time
- exact scenario scheduling itself
- accumulation back into output density one scenario point at a time
- exact currently reuses `RankedBatchPlanner` for its scenario planning cache,
  which is a strong sign the scheduler layer is still not canonicalized

Status:

- `Batched leaf math inside an exact scheduler shell`
- `No longer a separate scalar exact runtime`

## 7. Outcome-Mass And Mixture Helpers Still Keep Scalar Outer Loops

### 7.1 Outcome-Mass Batching Is Real, But Outer Loops Remain

`mix_outcome_mass_batch_idx()` is now the common nonresponse/outcome-mass batch
helper.

Current behavior:

- shared-trigger outcome-mass uses mask-expanded SoA batches
- non-shared-trigger outcome-mass batches repeated param points and can chunk
  unique points
- repeated parameter-value groups are compressed before evaluation
- inner probability work uses
  `integrate_outcome_probability_from_density_batch_idx()`

What remains scalar-shaped:

- loops over outcomes
- loops over components
- reduces component-use contributions sequentially

Status:

- `Batched probability integration under scalar outcome/component orchestration`

### 7.2 Native Trial Mixture Is Still Scalar Outside The Batch Leaf

`native_trial_mixture_internal_idx()`:

- assembles mixture terms sequentially
- evaluates each term by calling the batch density entry on a one-time vector
  `{t}`
- accumulates the result term-by-term

Status:

- `Scalar outer orchestration with one-point batching inside`

### 7.3 Parameter-Set Batching Is Still Scalar-Only

`cpp_loglik_multiple()` still loops `params_list` and calls `cpp_loglik()`
one-by-one.

Status:

- `Scalar-only at the parameter-set level`

## 8. Shared-Trigger Handling Is Materially Improved, But Not Universal

What now stays on a batched path:

- direct density via `(mask x point)` expansion in `node_density_entry_batch_idx()`
- ranked grouped shared-trigger sequences in `cpp_loglik()`
- shared-trigger outcome-mass through mask-expanded batch probabilities

What still gates shared-trigger batching:

- cases rejected by `shared_trigger_mask_batch_supported(...)`
- outcome contexts that depend on alias-source or guess-donor machinery
- any caller that only exposes a scalar API and therefore uses one-point
  batching around the shared-trigger batch machinery

Status:

- `Shared-trigger batching is real`
- `Still partial, not universal`

## 9. Alias And Guess-Donor Work Is Vector In Time, But Still Serial Over Donors

The current donor/alias accumulation helpers in `src/distributions.cpp` do
batch over `times`, but they still loop donors or source outcomes sequentially:

- `accumulate_component_guess_density_batch_idx(...)`
- `accumulate_outcome_alias_density_batch_idx(...)`

That means:

- no separate scalar time loop
- but still scalar orchestration across donor/source branches

Status:

- `Vector across time`
- `Scalar across donor/source structure`

## 10. Transitional Runtime Residue Still Exists

The runtime foundation is real, but some transitional architecture remains
visible:

- `VectorProgramDomain` still carries `OutcomeCoupling`
- `VectorPayloadRole` still carries pair / n-way / generic-node role families
- ranked hot execution still depends on ranked-specific planner/slice types
  instead of evaluating the canonical `VectorProgram` directly

These are not proof of scalar execution by themselves, but they are proof that
runtime consolidation is not finished.

Status:

- `Architecture still transitional even where batching is already good`

## 11. Validation And Tooling

Current correctness signal:

- `Rscript dev/scripts/check_loglik_golden.R` passes on the current tree

Current tooling caveat:

- older Step 3-6 structure scripts grep for historical symbols that no longer
  exist
- they are not a reliable source of truth for the current runtime shape without
  being updated

Current instrumentation:

- `src/runtime_stats.cpp` tracks batch / lane / adaptive-segment counters
- it does not provide a complete, authoritative scalar-vs-vector census of
  every remaining orchestration path

## 12. Updated Inventory Summary

Fully or mostly vectorized now:

- canonical tree vector program and shared lane batch
- top-level grouped trial execution in `cpp_loglik()`
- direct node density batch execution
- shared-trigger direct density via mask-expanded SoA batching
- batched exact leaf evaluation
- batched exact competitor survival
- ranked leaf evaluation under grouped frontier slices
- batched outcome-mass inner probability integration

Still scalar-shaped today:

- ranked frontier scheduling / collapse
- exact scenario scheduling / accumulation
- infinite-tail probability window scheduling
- outcome-mass outcome/component loops
- native mixture term loops
- alias / guess donor/source loops
- parameter-set loop in `cpp_loglik_multiple()`

## 13. Honest Status

This codebase is no longer in the "scalar with a few vector helpers" phase.

It is now in the "vectorized core with remaining scheduler shells and scalar
outer orchestration" phase.

If the standard is:

- "Did we remove the old obvious scalar fallback engines from production hot
  paths?": mostly yes.
- "Do we have one canonical vector runtime end-to-end?": not yet.
