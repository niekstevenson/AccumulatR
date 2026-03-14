# Scalar To Vector Inventory

This document inventories every current scalar or scalar-shaped execution path I could find in the likelihood evaluation framework reachable from `cpp_loglik()` and the closely related native probability / mixture helpers.

Scope:

- Included:
  - `src/distributions.cpp`
  - `src/exact_outcome_density.cpp`
  - `src/coupling_eval.cpp`
  - `src/ranked_transitions.cpp`
  - `src/trial_params.h`
  - helper entry points that call the same machinery, such as `cpp_loglik_multiple()`, `native_outcome_probability_*()`, and `native_trial_mixture_*()`
- Not included:
  - unrelated scalar utilities outside likelihood evaluation
  - debug-only exporters except where they mirror production execution

Terminology:

- `Fully batched`: one vectorized execution over many trials / times without dropping to per-point scalar evaluators
- `Pointwise batch`: accepts vector inputs, but loops over points and calls scalar evaluation per point
- `Scalar orchestration`: outer control flow is scalar even if inner kernels evaluate vectors
- `Scalar-only`: no vector/batch entry point is used for that case

This file is intentionally blunt: if a path still loops trials, components, trigger masks, sequence states, or points in C++ rather than grouping them into a batch kernel, it is listed here.

## 1. Top-Level `cpp_loglik()` Coverage Gaps

The top-level multi-trial batching in `cpp_loglik()` is currently narrow. It only batches:

- grouped direct single-step sequence-density trials via `single_step_batch_groups`
  - `src/distributions.cpp:4725-4803`
  - executed at `src/distributions.cpp:4834-4892`
- grouped ranked sequence trials via `ranked_sequence_batch_groups`
  - `src/distributions.cpp:4731-4831`
  - executed at `src/distributions.cpp:4894-4957`

Everything not caught by those two groupers falls through the per-trial loop:

- per-trial fallback loop: `src/distributions.cpp:4962-5021`
- status: `Scalar-only at the trial level`

The current top-level non-batched cases are:

- Any grouped direct single-step case with group size `< 2`
  - `src/distributions.cpp:4834-4838`
- Any grouped ranked case with group size `< 2`
  - `src/distributions.cpp:4895-4899`
- Any trial with a shared-trigger plan
  - skipped before batching at `src/distributions.cpp:4765-4777`
- Any `OutcomeMass` contribution
  - `evaluate_trial_contribution_kernel_idx()` routes those to `mix_outcome_mass_idx()` per trial at `src/distributions.cpp:4054-4069`
- Any sequence trial rejected by the direct batch-spec builder
  - `build_single_step_direct_batch_spec()`: `src/distributions.cpp:3449-3501`
  - rejection conditions include:
    - `probability_transform != Identity`
    - `sequence_length != 1`
    - multiple contributions
    - contribution not lowered as `LowerLayerDirect`
    - non-finite weight / time
- Any ranked trial rejected by the ranked batch-spec builder
  - `build_ranked_sequence_batch_spec()`: `src/distributions.cpp:3503-3555`
  - rejection conditions include:
    - `probability_transform != Identity`
    - `sequence_length <= 1`
    - multiple contributions
    - contribution not lowered as `SequenceState`
    - inconsistent ranked metadata sizes

`cpp_loglik_multiple()` is also scalar over parameter sets:

- loops `params_list` and calls `cpp_loglik()` one-by-one
- `src/distributions.cpp:5033-5050`
- status: `Scalar-only at the parameter-set level`

## 2. Shared-Trigger Enumeration Is Scalar Everywhere

Shared-trigger evaluation still enumerates trigger masks sequentially:

- `evaluate_shared_trigger_states()`: `src/trial_params.h:97-157`
- loops `2^k` Gray-code states one-by-one
- status: `Scalar-only over trigger states`

This scalar mask enumeration is used in all of these places:

- trial likelihood evaluation
  - `src/distributions.cpp:4221-4225`
- single-time density evaluation
  - `src/distributions.cpp:2717-2724`
- outcome-mass probability evaluation
  - `src/distributions.cpp:2807-2810`

Implication:

- Any model with shared triggers is currently scalar over trigger states, even when inner node or coupling evaluation reuses kernel runtime state.
- `cpp_loglik()` also skips these trials before its multi-trial batch grouping:
  - `src/distributions.cpp:4774-4776`

## 3. Contribution And Component Loops Are Still Scalar

These are not trial-vectorized today:

- Per-trial contribution accumulation inside `evaluate_trial_probability_kernel_idx()`
  - `eval_single_state()`: `src/distributions.cpp:4175-4183`
  - loops contributions sequentially
  - status: `Scalar orchestration`
- Outcome-mass mixture over components
  - `mix_outcome_mass_idx()`: `src/distributions.cpp:3182-3238`
  - loops components sequentially
  - status: `Scalar orchestration`
- Native trial mixture over components
  - `native_trial_mixture_internal_idx()`: `src/distributions.cpp:3042-3145`
  - loops components sequentially
  - status: `Scalar orchestration`

These loops matter because even if an inner density call is batched over quadrature nodes or times, the framework still serializes across contributions / mixture components.

## 4. Ranked Likelihood Still Has A Separate Scalar Engine

### 4.1 Single-Trial Or Ungrouped Ranked Work

Single-trial ranked evaluation still dispatches to the old scalar sequence engine:

- dispatch site:
  - `src/distributions.cpp:4012-4041`
- scalar engine:
  - `sequence_prefix_density_resolved()`
  - `src/ranked_transitions.cpp:703-892`
- status: `Scalar-only`

This is still used whenever ranked work does not enter the grouped batch path in `cpp_loglik()`.

### 4.2 Batched Ranked Evaluation Is Still Scalar-Shaped Internally

The batched ranked engine exists:

- `sequence_prefix_density_batch_resolved()`
- `src/ranked_transitions.cpp:894-1185`

But it is not a fully vectorized ranked runtime. It still:

- loops `TrialSequenceState`s sequentially
- builds `ExactScenarioPoint`s one state at a time
- collapses states through hash maps
- updates next states one scenario at a time

Relevant code:

- state assembly and denom preparation:
  - `src/ranked_transitions.cpp:988-1061`
- scenario collection:
  - `src/ranked_transitions.cpp:1093-1098`
- scenario application / next-state collapse:
  - `src/ranked_transitions.cpp:1144-1185`

Status:

- `Scalar orchestration with batched node evaluation inside`

### 4.3 Old Scalar Transition Enumeration Still Exists

The old ranked transition enumerator still exists and is active:

- `for_each_sequence_node_transition()`
  - `src/ranked_transitions.cpp:287-492`
- uses scalar recursive node evaluation internally
  - `src/ranked_transitions.cpp:346-358`

It is still used by:

- the scalar ranked engine:
  - `src/ranked_transitions.cpp:802-865`
- batched debug-oracle checks:
  - `src/ranked_transitions.cpp:1114-1131`

Status:

- `Scalar-only`

## 5. Exact Density: Scalar, Pointwise, And Scalar-Shaped Paths Still Exist

### 5.1 Scalar Exact Density Still Exists

Scalar exact density remains a live production path:

- `exact_outcome_density_from_state()`
  - `src/exact_outcome_density.cpp:1760-1770`
- called from `node_density_with_competitors_from_state_dense()`
  - `src/distributions.cpp:2387-2412`

The scalar exact implementation:

- collects scenarios at one time
  - `src/exact_outcome_density.cpp:1673-1683`
- loops scenarios one-by-one
  - `src/exact_outcome_density.cpp:1686-1714`
- applies competitor survival one scenario at a time when needed
  - `src/exact_outcome_density.cpp:1692-1708`

Status:

- `Scalar-only`

### 5.2 Batched Exact Density Still Has Scalar-Shaped Scenario Expansion

The batched exact density path exists:

- `exact_outcome_density_batch_from_state()`
  - `src/exact_outcome_density.cpp:1773-1850`

But it is still scalar-shaped in these places:

- scenario expansion is step-by-step and point-by-point
  - `exact_collect_scenarios_batch_from_points_impl()`
  - `src/exact_outcome_density.cpp:1270-1465`
- density accumulation loops scenario points one-by-one
  - `src/exact_outcome_density.cpp:1793-1845`

Status:

- `Scalar orchestration with batched node evaluation inside`

### 5.3 Batched Exact Density Still Has Scalar Competitor Fallback

If `exact_competitor_survival_batch()` fails or returns a mismatched size:

- `src/exact_outcome_density.cpp:1805-1845`

the code falls back to:

- constructing `NodeEvalState` per scenario
- calling `competitor_survival_from_state_compiled_ops()` per scenario

Status:

- `Scalar-only fallback inside batched exact density`

### 5.4 Exact Node Batch Has An Explicit Pointwise Fallback

Main batch node evaluator:

- `exact_eval_node_batch_with_points()`
- `src/exact_outcome_density.cpp:1164-1267`

Direct event fastpath is tried first:

- `src/exact_outcome_density.cpp:1181-1193`

If that cannot handle the node, the code falls back to pointwise/shared-state execution:

- `pointwise_fallback`: `src/exact_outcome_density.cpp:1204-1257`

Status:

- `Pointwise batch fallback`

### 5.5 Shared-State Check Is All-Or-Nothing

Exact batching only uses the shared-state kernel path if all active points match on the relevant filtered state:

- relevant-state comparison:
  - `src/exact_outcome_density.cpp:1113-1161`
- shared-state gating:
  - `src/exact_outcome_density.cpp:567-620`

If any active point differs in filtered forced bits or time constraints:

- the whole batch drops to fallback

Status:

- `Pointwise batch fallback caused by state heterogeneity`

### 5.6 Event Batch Can Fall Back To Per-Point Scalar Evaluation

Event batch common path:

- `src/exact_outcome_density.cpp:1010-1055`

Per-point fallback implementation:

- `exact_eval_event_pointwise_batch()`
- `src/exact_outcome_density.cpp:909-935`

That function loops points and calls scalar event evaluation per point:

- `evaluator_eval_event_ref_idx(...)` at `src/exact_outcome_density.cpp:923-929`

Status:

- `Pointwise batch`

### 5.7 Guard Batch Can Fall Back To Per-Point Scalar Evaluation

Guard batch common path:

- `src/exact_outcome_density.cpp:1058-1104`

Per-point fallback implementation:

- `exact_eval_guard_pointwise_batch()`
- `src/exact_outcome_density.cpp:974-1008`

This loops points and uses scalar guard evaluation for each point:

- `shared_guard_eval(op, kernel_batch_values_at(...))`
  - `src/exact_outcome_density.cpp:996-1005`

Status:

- `Pointwise batch`

### 5.8 Some "Batched" Exact Event Work Is Still Per-Point Inside The Batch Loop

`exact_eval_simple_acc_event_batch()` is only fully vector-kernel fast for direct shifted-time cases:

- `src/exact_outcome_density.cpp:132-420`

It still does per-point scalar work in these cases:

- non-direct relative-onset base evaluation
  - inside the batch loop at `src/exact_outcome_density.cpp:267-298`
  - helper does per-point source integration / scalar source eval:
    - `src/exact_outcome_density.cpp:714-905`
- bounded self-conditioning on a relative-onset node
  - lower / upper bound calls back into the same per-point helper
  - `src/exact_outcome_density.cpp:345-383`

Status:

- `Pointwise work inside a batched API`

### 5.9 Exact Event Fastpath Still Rejects Some Scalar-Capable Cases

Current direct event fastpath rejections include:

- density for outcomes with alias sources or guess donors
  - `src/exact_outcome_density.cpp:41-49`
  - rejection at `src/exact_outcome_density.cpp:148-150`
- self exact-time conditioning at a different evaluation time
  - classifier rejects at `src/exact_outcome_density.cpp:102-105`

When these happen, the code does not stay in the direct vector event path.

Status:

- `Pointwise batch fallback or scalar fallback depending on caller`

## 6. Coupling / Non-Ranked Paths Are Not Fully Vectorized

### 6.1 `evaluate_labelref_batch()` Has A Pointwise Scalar Fallback

Label-ref batch entry:

- `src/coupling_eval.cpp:237-300`

It is fully vectorized only when `build_coupling_acc_ref()` succeeds:

- `src/coupling_eval.cpp:64-87`

That specialization is rejected for:

- non-accumulator refs
- inactive components
- any chained onset
  - rejection: `ctx.has_chained_onsets && onset_kind != ONSET_ABSOLUTE`
  - `src/coupling_eval.cpp:79-82`

When specialization fails, `evaluate_labelref_batch()` loops times and calls scalar event evaluation per point:

- `eval_event_unforced(...)`
  - `src/coupling_eval.cpp:17-30`
- pointwise loop:
  - `src/coupling_eval.cpp:286-299`

Status:

- `Pointwise batch`

This affects Generic, Pair, and NWay coupling programs, because they all call `evaluate_labelref_batch()`.

### 6.2 Outcome Coupling Is Per-Trial And Per-Integral, Not Trial-Vectorized

`evaluate_outcome_coupling_unified()`:

- `src/coupling_eval.cpp:751-940`

It evaluates one upper bound at a time. There is no multi-trial batch entry for outcome-mass coupling in `cpp_loglik()`.

Status:

- `Scalar-only at the trial level`

### 6.3 Adaptive Coupling Integration Is Scalar Over Segments

Adaptive integration driver:

- `integrate_coupling_mass_batch_adaptive()`
- `src/coupling_eval.cpp:621-747`

This function:

- maintains a scalar frontier of segments
- accepts / splits segments sequentially
- evaluates batched quadrature nodes within each segment

Status:

- `Scalar orchestration with batched quadrature-node evaluation inside`

### 6.4 Generic Coupling Uses The Label-Ref Pointwise Fallback Too

Generic coupling provider:

- `evaluate_generic_terms_resolved_provider()`
- `src/coupling_eval.cpp:435-523`

If it takes the label-ref provider path:

- `src/coupling_eval.cpp:502-521`

then it relies on `evaluate_labelref_batch()`, which is pointwise whenever specialization fails.

Status:

- `Pointwise batch when label-ref specialization fails`

### 6.5 Pair And N-Way Coupling Still Orchestrate Scalar Integrals Per Trial

N-way program:

- `src/coupling_eval.cpp:791-854`

Pair program:

- `src/coupling_eval.cpp:855-939`

Both:

- run one integral at a time
- call `evaluate_labelref_batch()` repeatedly over quadrature-node batches
- are therefore only as vectorized as the label-ref path underneath

Status:

- `Scalar orchestration with batched or pointwise-batch inner evaluation`

## 7. Direct Node-Density Batch Has Hard Stops

Top-level node-density batch entry:

- `node_density_entry_batch_idx()`
- `src/distributions.cpp:2564-2603`

It refuses shared-trigger batch evaluation outright:

- returns `false` when `shared_trigger_count(*plan_ptr) > 0`
- `src/distributions.cpp:2593-2602`

Status:

- `Scalar-only fallback at the caller level`

This is why `cpp_loglik()` skips shared-trigger trials before grouping:

- `src/distributions.cpp:4774-4776`

## 8. Standalone Probability / Mixture Helpers Still Expose Scalar APIs

These helpers are scalar by interface even if they use batched kernels internally:

- `native_outcome_probability_bits_impl_idx()`
  - `src/distributions.cpp:2730-2811`
  - one upper bound at a time
- `native_trial_mixture_internal_idx()`
  - `src/distributions.cpp:3042-3145`
  - one trial / one time at a time

Status:

- `Scalar-only API surface`

## 9. Instrumentation Gaps

The current runtime stats do not yet prove scalar usage exhaustively.

Known gap:

- `generic_scalar_fallback_calls` exists in the stats struct and R export:
  - `src/runtime_stats.h:24-26`
  - `src/runtime_stats.cpp:176-182`
- but I do not see any increment site for it

What is currently observable:

- exact scalar density calls
  - `src/runtime_stats.cpp:30-35`
- exact batch calls
  - `src/runtime_stats.cpp:44-60`
- exact node pointwise fallback calls
  - `src/runtime_stats.cpp:62-70`
- exact competitor batch calls
  - `src/runtime_stats.cpp:72-90`
- generic label-ref batch fastpath calls
  - `src/runtime_stats.cpp:106-110`
- generic node-ref batch calls
  - `src/runtime_stats.cpp:113-117`

What is not currently counted:

- `evaluate_labelref_batch()` pointwise fallback usage
- shared-trigger mask enumeration usage
- top-level `cpp_loglik()` trials skipped from batching because group size is 1
- grouped-batch spec rejections by reason

## 10. Practical Reading Of The Current State

Today the framework is best described as:

- `Fully batched` for some direct single-step and grouped ranked `cpp_loglik()` cases
- `Pointwise batch` for several inner exact and coupling paths
- `Scalar orchestration` for shared triggers, contribution accumulation, component accumulation, adaptive coupling integration, ranked state scheduling, and exact scenario expansion
- `Scalar-only` for single-trial / ungrouped ranked, top-level `OutcomeMass`, `cpp_loglik_multiple()`, and the standalone scalar helper APIs

In short:

- Ranked is not the only remaining scalar area.
- Non-ranked is not fully vectorized.
- The biggest remaining scalar islands are:
  - top-level per-trial fallback in `cpp_loglik()`
  - shared-trigger enumeration
  - scalar ranked engine for single-trial / ungrouped ranked work
  - pointwise exact and label-ref fallback paths
