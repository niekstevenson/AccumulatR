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
- repeated leftover trials with identical eval shape, params, trigger plan, and
  component weights are now grouped and evaluated once before that loop
- status: `Scalar orchestration at the trial-group level; unique leftovers still scalar`

The current top-level non-batched cases are:

- Shared-trigger trials that are not direct single-step sequence-density groups
  - direct single-step shared-trigger trials can now enter grouped `(trial, mask)` batching in `cpp_loglik()`
  - ranked shared-trigger trials can now also batch `(mask x trial)` sequences in `cpp_loglik()`
  - shared-trigger `OutcomeMass` trials can now bypass the outer Gray-code loop when every contribution is batchable
- Non-shared-trigger `OutcomeMass` trials that do not form a usable batch group
  - grouped nonresponse `OutcomeMass` trials can now batch through coupling in `cpp_loglik()`
  - repeated parameter-value groups inside those batches are compressed before evaluation
  - grouped batches now run through the same `mix_outcome_mass_batch_idx()` helper on compressed point sets and bounded chunks
  - singleton no-trigger `OutcomeMass` evaluation also routes through that same batch helper
  - shared-trigger `OutcomeMass` now expands masks and routes through that same batch helper too
  - there is no longer a separate scalar `OutcomeMass` implementation; the remaining fallback is one-point batching, not a different engine
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

## 2. Shared-Trigger Enumeration Is Only Partially Batched Now

There is now a real shared-trigger mask batch path for direct node density:

- `build_shared_trigger_mask_soa_batch()` materializes per-mask `TrialParamsSoA`
- `node_density_entry_batch_idx()` expands `(mask x point)` batches for shared-trigger direct density
- `cpp_loglik()` can group direct single-step shared-trigger trials and batch `(trial, mask)` points together
- `cpp_loglik()` can also group ranked shared-trigger trials and batch `(mask x trial)` ranked sequences together
- exact-scenario-required direct density can also stay on the batched path now
- status: `Batched for direct density and ranked sequence density; not universal`

Scalar shared-trigger enumeration is still used in these places:

- trial likelihood evaluation
  - fallback path in `evaluate_trial_probability_kernel_idx()` for shared-trigger trials that are not covered by the direct / ranked / outcome-mass batch paths
- single-time density evaluation
  - fallback in `node_density_entry_idx()` when shared-trigger batching is unsupported
- outcome-mass probability evaluation
  - non-batchable coupling cases still fall back to scalar shared-trigger enumeration inside `native_outcome_probability_bits_impl_idx()`

Implication:

- Shared-trigger direct single-step density is no longer scalar by default.
- Shared-trigger ranked likelihood is no longer scalar by default in `cpp_loglik()`.
- Shared-trigger `OutcomeMass` is no longer scalar by default when the underlying coupling program is batchable.
- Alias-source and guess-donor `OutcomeMass` paths are now batchable through coupling.
- Shared-trigger enumeration still remains for direct-density and ranked cases rejected by `shared_trigger_mask_batch_supported()`, plus any coupling case that still misses both batch paths.

## 3. Contribution And Component Loops Are Still Scalar

These are not trial-vectorized today:

- Per-trial contribution accumulation inside `evaluate_trial_probability_kernel_idx()`
  - `eval_single_state()`: `src/distributions.cpp:4175-4183`
  - loops contributions sequentially
  - status: `Scalar orchestration`
- Outcome-mass mixture over components
  - `mix_outcome_mass_batch_idx()` / `mix_outcome_mass_idx()` in `src/distributions.cpp`
  - still loops outcomes and components sequentially inside the batch helper
  - status: scalar-shaped orchestration still present, but the separate scalar `OutcomeMass` engine is gone
- Native trial mixture over components
  - `native_trial_mixture_internal_idx()`: `src/distributions.cpp:3042-3145`
  - loops components sequentially
  - status: `Scalar orchestration`

These loops matter because even if an inner density call is batched over quadrature nodes or times, the framework still serializes across contributions / mixture components.

## 4. Ranked Likelihood Still Carries Scalar Internals

### 4.1 Single-Trial Or Ungrouped Ranked Work

Production ranked dispatch now routes single-trial `SequenceState` work through the batched ranked engine:

- dispatch site:
  - `src/distributions.cpp:4012-4049`
- top-level grouped ranked batching no longer skips size-1 groups:
  - `src/distributions.cpp:4903-4954`
- status: `No longer scalar-only in production`

The old scalar ranked engine still exists, but it is no longer the intended production path for single-trial ranked work.

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

The old ranked transition enumerator still exists and is still scalar:

- `for_each_sequence_node_transition()`
  - `src/ranked_transitions.cpp:287-492`
- uses scalar recursive node evaluation internally
  - `src/ranked_transitions.cpp:346-358`

It is still used by:

- the scalar ranked engine kept as an oracle:
  - `src/ranked_transitions.cpp:703-892`
- batched debug-oracle checks:
  - `src/ranked_transitions.cpp:1114-1131`

Status:

- `Scalar-only`, but no longer the primary production path if ranked batching covers the case

## 5. Exact Density: Pointwise And Scalar-Shaped Paths Still Exist

### 5.1 Singleton Exact Density Now Routes Through The Batch Engine

`exact_outcome_density_from_state()` still exists as a helper entry point:

- `src/exact_outcome_density.cpp`
- called from `node_density_with_competitors_from_state_dense()`
  - `src/distributions.cpp:2387-2412`

But it no longer owns a separate scalar exact implementation. Singleton exact
density now delegates to `exact_outcome_density_batch_from_state()` with a
one-point batch.

Status:

- `No separate scalar exact engine in production`

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

If that cannot handle the node, the code now tries shared-state subgroup
partitioning before it gives up to pointwise execution:

- `pointwise_fallback`: `src/exact_outcome_density.cpp:1204-1257`

Status:

- `Pointwise batch fallback still exists, but only after subgroup partitioning
  fails`

### 5.5 Shared-State Check Is All-Or-Nothing

Exact batching no longer requires the full active set to share one filtered
state. It now partitions active points into relevant-state-compatible subgroups
and evaluates those subgroups through the shared-state kernel path first:

- relevant-state comparison:
  - `src/exact_outcome_density.cpp:1113-1161`
- shared-state gating:
  - `src/exact_outcome_density.cpp:567-620`

Fallback still happens if relevant-state partitioning cannot build valid
subgroups or if subgroup shared evaluation fails:

- the batch drops to pointwise fallback
- event / guard shared-state batching also now bails out when partitioning would
  create one subgroup per active point, because that pattern caused severe exact
  benchmark regressions while providing no real batching benefit

Status:

- `Improved: heterogeneity can stay batched when it partitions cleanly`

### 5.6 Event Batch Has A Targeted Per-Point Fallback Again

Event batch common path:

- `src/exact_outcome_density.cpp:1010-1055`

Current behavior:

- accumulator-backed event refs recurse into `exact_eval_node_batch_with_points()`
- pooled event refs stay in the grouped shared-state batch path
- missing event-node metadata is still treated as invalid IR
- but if relevant-state partitioning would degenerate into singleton groups, the
  pooled-event path now falls back to per-point exact event evaluation instead
  of doing one-point "batched" work repeatedly

Status:

- `Targeted pointwise batch fallback`

### 5.7 Guard Batch Has A Targeted Per-Point Fallback Again

Guard batch common path:

- `src/exact_outcome_density.cpp:1058-1104`

Current behavior:

- grouped shared-state guard evaluation is still used
- missing `op.node_idx` metadata is still treated as invalid IR
- but if relevant-state partitioning would degenerate into singleton groups, the
  guard path now falls back to per-point exact guard evaluation instead of doing
  one-point shared-state batches repeatedly

Status:

- `Targeted pointwise batch fallback`

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

### 6.1 `evaluate_labelref_batch()` No Longer Has A Scalar Execution Loop

Label-ref batch entry:

- `src/coupling_eval.cpp:237-300`

It now has three real batch paths:

- accumulator-ref specialization through `build_coupling_acc_ref()`
- pooled label-ref batching through recursive member batches
- node-batch fallback when coupling metadata can resolve a real event node

The accumulator specialization is rejected for:

- non-accumulator refs
- inactive components
- any chained onset
  - rejection: `ctx.has_chained_onsets && onset_kind != ONSET_ABSOLUTE`
  - `src/coupling_eval.cpp:79-82`

Pooled refs no longer need to fall through the node-batch path just to recurse member-by-member through scalar event evaluation.

In normal coupling execution with valid kernel metadata:

- accumulator refs use the accumulator batch specialization
- pooled refs use pooled label-ref batching
- remaining event refs go through node-batch, including donor-density cases

If a coupling ref now misses all three routes, `evaluate_labelref_batch()` returns failure instead of silently re-entering scalar event evaluation. Pair and n-way coupling payload construction was tightened so valid IR must resolve every coupling event ref back to a concrete event node.

Status:

- `Batched for accumulator-backed refs, pooled refs, donor-aware node-backed refs, and shared-trigger per-point params; no scalar loop remains inside coupling label-ref batch evaluation`

This affects Generic, Pair, and NWay coupling programs, because they all call `evaluate_labelref_batch()`.

### 6.2 Outcome Coupling Is Now Batched Across Per-Point Trial Params, But Not Fully Trial-Vectorized

`evaluate_outcome_coupling_unified()`:

- `src/coupling_eval.cpp:751-940`

There is now a per-point-param batch companion for the same coupling programs:

- shared-trigger mask-state batches can evaluate one coupling integral across many `TrialParamsSoA` points
- `cpp_loglik()` can use that for shared-trigger `OutcomeMass`
- direct `native_outcome_probability_bits_impl_idx()` calls can also use it

Status:

- `Batched across shared-trigger mask points; still scalar across outcomes / components / trials at the top level`

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

### 6.4 Generic Coupling Still Depends On Scalar-Shaped Adaptive Integration

Generic coupling provider:

- `evaluate_generic_terms_resolved_provider()`
- `src/coupling_eval.cpp:435-523`

If it takes the label-ref provider path it now stays inside accumulator, pool,
or node-batch evaluation. The remaining scalar shape here is the adaptive
segment driver around those batched inner calls, not a scalar label-ref loop.

Status:

- `Scalar orchestration with batched inner evaluation`

### 6.5 Pair And N-Way Coupling Still Orchestrate Adaptive Segments Scalar-Style

N-way program:

- `src/coupling_eval.cpp:791-854`

Pair program:

- `src/coupling_eval.cpp:855-939`

Both:

- still split / accept adaptive segments sequentially
- now batch quadrature-node work across shared-trigger mask points
- are therefore no longer scalar by default on shared-trigger mask batches, but still have scalar segment orchestration

Status:

- `Scalar orchestration with batched inner evaluation`

## 7. Direct Node-Density Batch Has Hard Stops

Top-level node-density batch entry:

- `node_density_entry_batch_idx()`
- `src/distributions.cpp:2564-2603`

It no longer refuses shared-trigger batch evaluation outright:

- direct shared-trigger density now batches over `(mask x point)` using per-point `TrialParamsSoA`
- duplicate `(time, TrialParamsSoA*)` points are compressed before batched evaluation, so repeated trials can share exact work
- unsupported shared-trigger cases still return `false` and fall back at the caller
- current hard stops include:
  - alias / guess-donor outcomes
  - cases that cannot use the simple direct shared-trigger batch path

Status:

- `Mixed: batched fast path with scalar fallback`

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

What is currently observable:

- exact batch calls
  - `src/runtime_stats.cpp:34-40`
- exact node pointwise fallback calls
  - `src/runtime_stats.cpp:48-56`
- exact competitor batch calls
  - `src/runtime_stats.cpp:58-76`
- generic label-ref batch fastpath calls
  - `src/runtime_stats.cpp:92-96`
- generic node-ref batch calls
  - `src/runtime_stats.cpp:99-103`
- grouped fallback trial executions in `cpp_loglik()`
  - `src/runtime_stats.cpp`
  - increments in the repeated-leftover grouping step inside `cpp_loglik()`
- direct batch-spec attempts and rejection counts by reason
  - `src/runtime_stats.cpp:131-169`
  - increments in `build_single_step_direct_batch_spec()`
    - `src/distributions.cpp:3829-3873`
- ranked batch-spec attempts and rejection counts by reason
  - `src/runtime_stats.cpp:171-191`
  - increments in `build_ranked_sequence_batch_spec()`
    - `src/distributions.cpp:3875-3928`
- accepted direct batch execution counts
  - `cpp_loglik()` direct grouped batch calls / trials / exec failures
    - `src/runtime_stats.cpp`
    - increments in `src/distributions.cpp:5264-5428`
  - direct node-density batch calls / points / param-matrix usage / simple-fastpath hits / kernel-eval hits / competitor-product hits
    - `src/runtime_stats.cpp`
    - increments in `src/distributions.cpp:2652-2763`
  - kernel event-batch calls / points / param-matrix usage
    - `src/runtime_stats.cpp`
    - increments in `src/distributions.cpp:1494-1549`

What is not currently counted:

- the later top-level reasons a built spec still does not batch
  - e.g. unsupported shared-trigger support-gate checks after spec construction

## 10. Practical Reading Of The Current State

Today the framework is best described as:

- `Fully batched` for some direct single-step and grouped ranked `cpp_loglik()` cases
- `Pointwise batch` for several inner exact and coupling paths
- `Scalar orchestration` for contribution accumulation, component accumulation, adaptive coupling integration, ranked state scheduling, exact scenario expansion, and the remaining unsupported shared-trigger direct/ranked paths
- `Scalar-only` for unique top-level fallback trials, `cpp_loglik_multiple()`, and the standalone scalar helper APIs

In short:

- Ranked is not the only remaining scalar area.
- Non-ranked is not fully vectorized.
- The biggest remaining scalar islands are:
  - unique top-level fallback trials in `cpp_loglik()`
  - shared-trigger enumeration for unsupported direct / ranked paths
  - scalar ranked internals still retained as an oracle / fallback implementation
  - pointwise exact fallback paths, especially singleton-partition escape hatches

## 11. Recommended Next Cuts

Based on the current benchmark behavior, the next work should be prioritized like this:

### 11.2 Top-Level Per-Trial Fallback In `cpp_loglik()`

Even after the ranked, direct, and shared-trigger `OutcomeMass` cuts, anything rejected by the batch-spec builders still falls through the scalar trial loop.

Repeated identical leftover trials are now grouped and evaluated once, so the
remaining problem here is mostly unique leftovers rather than blindly repeated
work.

Current instrumentation changed the picture slightly:

- `example_1_simple` enters the direct batch-spec builder cleanly on every trial in a 10-trial probe
- `example_23_ranked_chain` enters the ranked batch-spec builder cleanly on every trial in a 10-trial probe
- `shared_gate_pair` is pure `OutcomeMass`, so these builders are not the relevant bottleneck there
- `example_22_shared_q` showed only occasional direct-spec input rejection in the same probe

Implication:

- this remains a real scalar island structurally
- but it is smaller now because repeated leftovers no longer each pay a full
  trial evaluation
- the remaining work here is either:
  - making more trials batchable up front
  - or shrinking the singleton / failed-batch leftovers that still miss the
    new non-shared-trigger `OutcomeMass` group path

Current execution-path probes sharpen that further:

- `example_1_simple` batches cleanly through `cpp_loglik()`
- `example_22_shared_q` also mostly batches cleanly at the direct spec level
- both now hit a dedicated simple-competing-race fast path in `node_density_with_competitors_batch_internal()`
  - direct batch calls observed
  - simple-competing fastpath calls observed
  - zero kernel-event batch calls in the updated probe
  - zero competitor-product calls in the updated probe
- on the tracked centralized rows this cut improved:
  - `example_1_simple`
  - `example_21_simple_q`
  - `example_22_shared_q`

That means the already-batched simple competing-race path is no longer the main direct-path regression target.
The remaining direct-path work is now mostly about cases that are still too complex for that specialization.

Target areas:

- unique leftover trials that still never enter the new non-shared-trigger
  `OutcomeMass` group path
- `node_density_with_competitors_batch_internal()` for non-simple competitor shapes
- `build_single_step_direct_batch_spec()`
- `evaluate_trial_contribution_kernel_idx()`

### 11.3 Unsupported Shared-Trigger Direct / Ranked Paths

The shared-trigger outcome-mass alias / guess cut landed, but direct-density and ranked shared-trigger batching still use the narrower support gate.

What remains:

- alias-source outcomes in shared-trigger direct density
- guess-donor outcomes in shared-trigger direct density
- any ranked shared-trigger case rejected because one step fails `shared_trigger_mask_batch_supported()`

Target areas:

- `shared_trigger_mask_batch_supported()`
- `node_density_entry_batch_idx()`
- `build_ranked_sequence_batch_spec()`
- `sequence_prefix_density_batch_resolved()`

### 11.4 Batched Exact Still Needs Deeper Inner Cuts

The separate scalar exact engine is gone, but the exact batch engine still has
real scalar-shaped internals:

- scenario expansion / collapse
- shared-state subgroup orchestration still exists inside the batch path
- singleton-partition escape hatches now drop event / guard evaluation to
  targeted per-point exact fallback to avoid catastrophic fake-batch behavior
- pooled-event shared batch evaluation is still scalar-shaped internally
- scalar competitor fallback inside the batch path

Target area:

- `src/exact_outcome_density.cpp`

### 11.5 Instrumentation Before More Broad Refactors

The runtime stats are still not complete enough to prove where all scalar work remains.

Before another large pass, the next useful instrumentation would be:

- label-ref pointwise fallback counts
- batch-spec rejection counts by reason
- exact scenario count / average scenarios per point
- exact competitor compiled-op counts on the remaining composite exact cases

This matters because the recent duplicate-point compression work already showed that the right optimization depends heavily on whether repeated times actually exist.
