# BEEST Exactness: Problem and Implementation Plan

## 1. The general problem, using the `go / stop / change` example

Consider the model:

- `S = guard(go by stop)`
- `X = and(stop, change)`

and let the trigger fail with probability `q`.

Ignoring the trigger for a moment and looking only at the triggered branch:

- observing `S` at time `t` means:

```text
go = t and stop > t
```

- observing `X` at time `t` means:

```text
max(stop, change) = t
and S did not already happen earlier
```

The exact densities are:

### Exact `S` density

```text
L_S^exact(t) = f_go(t) S_stop(t)
```

There is no extra `change` term. If `stop > t`, then `X` is impossible by `t`, because `X` requires `stop <= t`.

### Exact `X` density

There are two distinct ways for `X` to happen at `t`:

1. `stop = t` and `change <= t`
2. `change = t` and `stop = u < t` for some `u`

So

```text
L_X^exact(t)
= f_stop(t) F_change(t) S_go(t)
 + f_change(t) ∫_0^t f_stop(u) S_go(u) du
```

The second term depends on the internal stop time `u`.

## 2. What the framework is doing now

The current fast path reduces each node to a scalar density/CDF/survival and then combines outcomes as:

```text
density(target node at t) * survival(competing nodes at t)
```

For this model that becomes:

### Current native `S`

`guard(go by stop)` is treated as a node with

```text
f_guard(t) = f_go(t) S_stop(t)
```

and then the competitor `and(stop, change)` is treated as a separate node whose survival is

```text
1 - F_stop(t) F_change(t)
```

So native computes

```text
L_S^native(t) = f_go(t) S_stop(t) [1 - F_stop(t) F_change(t)]
```

But this is wrong for the exact race likelihood.

Why?

Because once we already know `stop > t`, we already know `X` has not happened by `t`.

So the correct competitor factor is

```text
P(X not by t | stop > t) = 1
```

not the unconditional

```text
P(X not by t) = 1 - F_stop(t) F_change(t)
```

### Current native `X`

`and(stop, change)` is treated as one node with density

```text
f_and(t) = f_stop(t) F_change(t) + f_change(t) F_stop(t)
```

and then native multiplies by survival of the competitor `guard(go by stop)`:

```text
L_X^native(t)
= [f_stop(t) F_change(t) + f_change(t) F_stop(t)] [1 - G_guard(t)]
```

where `G_guard(t)` is the CDF of `guard(go by stop)`.

This is also not exact.

The problem is the `change = t` branch:

- if `change = t` and `stop = u < t`, then the correct no-`S` condition is:

```text
go > u
```

not

```text
guard(go by stop) has not happened by t
```

So the current framework replaces a time-dependent condition inside the integral with a single coarse factor at time `t`.

## 3. The general mathematical failure

The current framework is exact only when this factorization is valid:

```text
f_exact(target at t and no competitor before t)
= f_target_node(t) * P(no competitor before t)
```

But for composite nodes this is often false.

The correct object is:

```text
f_exact(target at t)
= sum over scenarios s:
   f(target reaches t via scenario s) *
   P(no competitor before t | scenario s)
```

The scenario matters.

For the BEEST example, the relevant scenarios are things like:

- `go = t`, `stop > t`
- `stop = t`, `change <= t`
- `change = t`, `stop = u < t`

These scenarios imply different conditional constraints on competitors.

The current fast path collapses all those scenarios into one node-level density before competitor handling, so the conditional structure is lost.

## 4. What information must be preserved for exactness

To compute exact first-response likelihood, the evaluator must keep track of branch-level state:

- which source labels are known complete by `t`
- which source labels are known to survive past `t`
- which source labels are known to occur at exact times
- which source labels are only bounded in time

The framework already has representations for these:

- `forced_complete_bits`
- `forced_survive_bits`
- `exact_source_times`
- `source_time_bounds`

These already appear in the scalar evaluation machinery and in ranked sequence evaluation.

## 5. Why the right fix is not in the low-level kernel

The low-level kernel in `src/kernel_executor.cpp` is built around node-local algebra:

- `And`
- `Or`
- `Guard`
- scalar density / survival / CDF for each node

That is useful as a fast approximate/compositional engine, but it is the wrong layer for exact first-passage semantics.

Exactness needs:

```text
branch/scenario enumeration
-> conditional state propagation
-> competitor survival evaluated under that conditional state
```

That is a higher-level transition/state problem, not a kernel-local algebra problem.

## 6. Where the framework already has the right ingredients

The codebase already has most of what is needed:

- scalar conditional evaluation with exact times and bounds in `src/distributions.cpp`
- branch/state propagation in `src/ranked_transitions.cpp`
- specialized coupling/integration in `src/coupling_eval.cpp`

So the clean fix is:

- keep the current kernel fast path for cases where it is exact enough
- add an exact scenario/state evaluator above it

## 7. Recommended architecture

## Principle

Implement exactness on the **scalar semantics path** first, then reuse it from the **vector/integration path**.

That means:

- scalar path defines the meaning
- vector path just batches repeated scalar evaluations for quadrature
- reuse the existing ranked-transition/state machinery rather than introducing a parallel scenario system

## Terminology note: what "source = t" means

For continuous-time branches, phrases like

```text
stop = t
```

are shorthand for:

```text
this density branch is pinned to the current evaluation time t
```

This is not claiming positive probability mass at one exact time for an absolutely continuous variable. It is branch-density bookkeeping.

The one special case that does involve literal point mass is the onset atom already present in the runtime/test harness, where pre-onset latent mass is collapsed onto the onset time. That mixed discrete/continuous behavior must be handled explicitly in the oracle and in the implementation.

## Minimal-change architecture

The codebase already has most of the state representation needed:

- forced complete bits
- forced survive bits
- exact source times
- source time bounds
- ranked branch/scenario enumeration

The main missing piece is narrower than the original draft implied:

- target scenario emission currently carries forced bits well
- it does not explicitly bind "these source labels are pinned to the current time `t` in this density branch"
- the coupling API does not yet thread `exact_source_times` and `source_time_bounds`

So the smallest credible fix is:

1. keep using ranked-transition scenario enumeration
2. extend emitted scenario state with current-time bindings for density branches
3. evaluate competitor survival under that full scenario state
4. reuse that scalar evaluator from coupling integrals

## New scalar primitive

Add a scalar evaluator:

```text
exact_outcome_density_at_time(ctx, node_idx, competitor_ids, t, ...)
```

This computes:

```text
sum over target scenarios s:
  target_scenario_density(s, t) *
  competitor_survival(t | scenario state s)
```

It should be defined in terms of the existing stateful evaluator machinery, not as a new kernel-local algebra.

## Extend ranked transitions narrowly

The ranked transition compiler is already close to what is needed, but it must preserve one more thing:

```text
"these sources are pinned to the current evaluation time t in this density branch"
```

That can be implemented either as:

- a dedicated step such as `BindExactCurrentSources`, or
- a narrow flag on the existing completion step

The important semantic rule is:

- only density-evaluated source branches get current-time bindings
- CDF-only completion means `source <= t`, not `source = t`

For example, for an event branch:

1. evaluate event density at `t`
2. mark its sources complete
3. bind those sources into `exact_source_times[source_id] = t`

That distinguishes:

- source pinned at `t`
from
- source only known complete by `t`

which is exactly the distinction the current generic path loses.

## Use scenario state when checking competitors

After emitting one target scenario, evaluate competitors with:

- that scenario's `forced_complete_bits`
- that scenario's `forced_survive_bits`
- that scenario's `exact_source_times`
- that scenario's `source_time_bounds`

The framework already knows how to consume those conditions at the scalar evaluation level. The exact fix is to make sure the target-side scenario emitter actually supplies the missing current-time information.

## 8. Concrete implementation steps

### Step 1: write a true exact BEEST oracle first

Before changing dispatch, write a trusted reference for the specific BEEST model.

Requirements:

- represent the exact `S` and `X` likelihood pieces, not the current native approximation
- handle the mixed discrete/continuous onset behavior explicitly
- make it easy to compare native vs exact at selected parameter settings and RTs

This can live in:

- `tests/testthat/helper-beest.R`, revised to become exact, or
- a new `dev/test_beest_exact.R` if preserving the current approximation is useful for side-by-side comparison

This step should happen first, because otherwise the implementation has no trustworthy regression target.

If ex-Gaussian onset semantics are changed as described below, the oracle should use those revised onset rules from the start.

### Step 2: reuse the existing ranked-transition scenario state

Do **not** introduce a second parallel `ExactScenarioState` abstraction unless later reuse truly requires it.

The current ranked code already carries the essential state shape:

- weight
- forced complete bits
- forced survive bits
- exact source times
- source time bounds

If sharing is needed beyond sequence-prefix code, lift that existing state into a reusable helper instead of inventing a duplicate representation.

### Step 3: add current-time binding to scenario emission

In `src/ranked_transitions.h` / `src/ranked_transitions.cpp`, extend transition execution so that density branches can bind current-time source labels into emitted scenario state.

Implementation options:

- add `BindExactCurrentSources`, or
- augment the relevant completion step with `bind_exact_current_time = true`

Semantic rule:

- event branch evaluated by density at `t` binds exact time
- sibling branch included only through CDF contributes completion/bounds but not exact current time
- guard blocker branches contribute survive state, not exact current time

### Step 4: generalize the scenario iterator rather than replacing it

Refactor `for_each_sequence_node_transition(...)` into a reusable scenario emitter for:

- ranked sequence evaluation
- single-outcome exact density evaluation

Something like:

```text
for_each_exact_node_scenario(...)
```

which emits:

- scenario weight
- complete bits
- survive bits
- exact times
- bounds

The point is to reuse the current ranked-transition compiler/runtime, not to rebuild scenario enumeration from scratch.

### Step 5: add `exact_outcome_density_at_time(...)`

Add a scalar exact evaluator, for example in:

- `src/exact_outcome_density.h`
- `src/exact_outcome_density.cpp`

or colocated with the existing density path if that keeps dependencies simpler.

Algorithm:

1. enumerate target scenarios for the observed node at time `t`
2. for each scenario, evaluate `competitor_survival_internal(...)` under that scenario state
3. sum the weighted competitor-adjusted branch densities

This is where the BEEST exactness should live semantically.

### Step 6: route single observed outcomes through the exact scalar path

In `src/distributions.cpp`, add targeted dispatch for cases that require exact scenario semantics.

Avoid a large execution-mode redesign if a narrow dispatch point is enough.

Goal:

- simple independent event-vs-event races stay on the current dense fast path
- BEEST-style composite/shared-source races use `exact_outcome_density_at_time(...)`

This is the first place where BEEST-style exactness should show up in the actual trial log-likelihood.

### Step 7: widen the coupling API to carry exact maps and bounds

Before coupling can reuse the scalar exact evaluator, the coupling layer must accept the same conditioning state.

So update the unified coupling interface and helpers to thread:

- `exact_source_times`
- `source_time_bounds`

through the generic path.

Only after that should the generic fallback integrand in `src/coupling_eval.cpp` call the scalar exact evaluator.

Keep these optimized shortcuts:

- specialized `Pair`
- specialized `NWay`
- simple no-forced/no-exact provider fast paths where they are still exact

### Step 8: only then add exactness dispatch

Do not send every model through the expensive exact path.

At prep time, mark outcomes/nodes as requiring exact evaluation when any of these hold:

- target node is composite
- competitor node is composite
- target and competitors share sources
- any guard is involved
- onset dependencies are involved
- shared triggers are involved

For plain independent event-vs-event races, keep the current fast path.

Dispatch should be added only after:

- the oracle exists
- the scalar exact evaluator is correct
- the coupling API can actually carry full conditional state

## 9. Scalar path or vector path?

The correct answer is:

- **exact semantics belong on the scalar path**
- **the vector path should reuse the scalar exact evaluator**

So:

- implement the meaning once, at one `t`
- batch it only for numerical integration

Trying to make exactness a vector-path feature first would solve the wrong problem. The problem is not batching; the problem is preserving conditional branch structure.

## 10. Practical rollout strategy

Recommended order:

1. Write the exact BEEST oracle and explicit BEEST comparison tests
2. Add current-time binding to ranked scenario emission
3. Expose a reusable exact-scenario iterator that emits exact times and bounds
4. Implement `exact_outcome_density_at_time(...)`
5. Route single observed outcomes that need exactness through the scalar exact path
6. Widen coupling APIs to carry exact maps/bounds
7. Reuse the scalar exact evaluator in generic coupling integrals
8. Add prep-time dispatch so simple models stay on the existing fast path
9. Optimize only after correctness is demonstrated

## 11. Separate issue: ex-Gaussian onset semantics

This is a separate semantics cleanup from the BEEST exactness fix, but it should be tracked explicitly because the oracle and runtime need to agree.

For ex-Gaussian accumulators, use:

- no onset atom
- onset as a hard support shift
- latent finishing times truncated to `x > 0`
- post-onset density renormalized analytically

That means:

```text
X ~ ex-Gaussian base distribution
X_obs ~ X | X > 0
T_obs = onset_eff + X_obs
```

with normalization constant

```text
Z = P(X > 0) = 1 - F_base(0)
```

So the shifted observed law becomes:

```text
f_obs(t) = 0                               if t < onset_eff
f_obs(t) = f_base(t - onset_eff) / Z       if t >= onset_eff
```

and

```text
F_obs(t) = 0                                                if t < onset_eff
F_obs(t) = (F_base(t - onset_eff) - F_base(0)) / Z          if t >= onset_eff
```

with

```text
S_obs(t) = 1 - F_obs(t)
```

Important consequences:

- `P(T_obs = onset_eff) = 0`
- there is no discrete onset mass to special-case
- this does not require numerical integration; it only requires evaluating `F_base(0)` once for the parameter set

Suggested implementation scope:

- update the ex-Gaussian onset handling in the core scalar distribution helpers
- update the BEEST oracle/reference code to match
- keep this semantics change logically separate from the shared-source exactness work

This issue is related to BEEST because the oracle must match the intended model semantics, but it is not the same bug as the native composite-race approximation failure.

## 12. Bottom line

The framework's current generic fast path is:

```text
node-local algebra first, competitor survival second
```

The exact first-response path should be:

```text
enumerate target scenarios with full conditional state,
bind current-time sources for density branches,
then evaluate competitors under that conditional state
```

The smallest credible fix is therefore not "build a second scenario engine".

It is:

- reuse the existing ranked-transition machinery
- add the missing current-time binding
- thread exact maps/bounds through coupling
- prove correctness against a true BEEST oracle before widening dispatch

Separately, ex-Gaussian onset semantics should be cleaned up so onset is implemented as truncation-plus-renormalization rather than as an onset atom.
