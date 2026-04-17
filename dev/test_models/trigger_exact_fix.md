# Trigger Exact Fix

## Problem

For the fixed-mixture stop-change model

```r
model <- race_spec() |>
  add_accumulator("S", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_accumulator("change", "lognormal") |>
  add_outcome("S", inhibit("S", by = "stop")) |>
  add_outcome("X", all_of("change", "stop")) |>
  add_component("go_only", members = "S", weight = .75) |>
  add_component("go_stop", members = c("S", "stop", "change"), weight = .25) |>
  add_trigger(
    "stop_trigger",
    members = c("stop", "change"),
    q = 0.05,
    param = "stop_trigger"
  ) |>
  set_mixture_options(mode = "fixed") |>
  set_parameters(list(
    m_go = "S.m",
    m_stop = "stop.m",
    m_change = "change.m",
    s_go = "S.s",
    s_stop = "stop.s",
    s_change = "change.s",
    t0_go = "S.t0",
    t0_change = "change.t0",
    q = c("stop.q", "change.q")
  )) |>
  finalize_model()
```

the intended component-conditioned density is

- `p(S at t | go_stop) = q f_S(t) + (1 - q) f_S(t) S_stop(t)`

but the current engine evaluates

- `p_engine(S at t | go_stop) = [q f_S(t) + (1 - q) f_S(t) S_stop(t)] [1 - F_change(t) F_stop(t)]`

At `t = 0.30`:

- intended: `0.6672304`
- engine: `0.6320930`

So the current engine is wrong on `go_stop -> S`.

## What Is Actually Wrong

There are two separate bugs.

### 1. The observed likelihood path is not exact-conditional

Shared-trigger propagation itself is mostly correct.

The code does correctly:

- build the trigger-state plan in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:1468)
- copy trial params into mutable scratch in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:1617)
- flip triggered accumulators between fail/succeed states in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:1679)
- enumerate and weight all trigger states in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:1716)

That part is not the primary failure.

The wrong step is later, inside the observed likelihood route:

- observed trial dispatch in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:6643)
- shared-trigger density route in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:4806)
- fused observed density in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:4486)

That route computes, per trigger state:

- `density(target outcome at t) * survival(competitor outcome at t)`

Concretely, in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:4517):

- read target density
- multiply by competitor survival

This is only valid when competitor survival is independent of the realized target transition.

For this model, that assumption is false.

If `S = inhibit(S, stop)` fires at time `t`, then the event `stop > t` is already known. Once `stop > t` is known, the competitor

- `X = all_of(change, stop)`

cannot have happened by `t`, regardless of `change`.

So the correct competitor factor after observing `S` is `1`, not marginal `survival(X at t)`.

The current observed path never conditions competitor survival on the realized `S` transition. It uses unconditional marginal competitor survival instead.

That is the main mathematical bug.

### 2. Absolute-onset event CDF ignores `q`

There is also a concrete low-level bug in event evaluation.

For absolute-onset events, the code currently does:

- density: `acc_density_from_cfg(t, onset, q, cfg)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:2190)
- survival: `acc_survival_from_cfg(t, onset, q, cfg)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:2193)
- cdf: `acc_cdf_success_from_cfg(t, onset, cfg)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:2196)

The problem is that `acc_cdf_success_from_cfg(...)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:739) ignores `q`.

So for an event with non-start probability `q`, the engine uses:

- density: `(1 - q) f(t)`
- survival: `q + (1 - q) S(t)`
- cdf: `F(t)`

That is inconsistent. The unconditional event cdf should be:

- `(1 - q) F(t)`

This matters immediately for

- `X = all_of(change, stop)`

because `And` uses child CDFs multiplicatively in [src/kernel_executor.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/kernel_executor.cpp:109).

So in trigger-fail states, where both `change` and `stop` should be absent, the engine still gives `X` a positive CDF and therefore less than unit survival.

That is why the current engine lands exactly on

- `[q f_S(t) + (1 - q) f_S(t) S_stop(t)] [1 - F_change(t) F_stop(t)]`

instead of the intended expression.

## Is There Already An Exact Conditional Path?

Yes.

The exact conditional machinery already exists for ranked likelihood.

The important pieces are:

- transition compilation in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:5781)
- transition enumeration in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:5992)
- exact source-time / forced-bit propagation in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:6036)
- competitor survival under the post-transition forced state in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:6370)

That path does the right kind of thing:

- realize one concrete target transition at time `t`
- add the resulting forced-complete / forced-survive source constraints
- evaluate competitor survival under those constraints

That is exactly what the observed `S` density needs here.

So the situation is:

- there already is an exact conditional path
- observed likelihood is not using it
- trigger propagation is not the core problem
- the observed route is using the wrong abstraction
- and there is a second separate `q`/CDF bug at the event leaf level

## Best Fix

The clean fix is:

1. fix event CDF semantics
2. route observed likelihood through the exact transition-conditioned machinery already used for ranked likelihood

Not:

- more special-case algebra
- more trigger-specific fast paths
- patching this one model only

### A. Fix the event CDF semantics first

Introduce an explicit unconditional event CDF helper:

- `acc_cdf_from_cfg(t, onset, q, cfg) = (1 - q) * F(t - onset - t0)`

Then use that in the absolute-onset event path in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:2188).

That change should make the event triple internally consistent:

- density = unconditional event density
- cdf = unconditional event cdf
- survival = unconditional event survival

The current `acc_cdf_success_from_cfg(...)` name is part of the problem. It mixes two semantics:

- success-given-start cdf
- unconditional event cdf

The code should stop pretending those are interchangeable.

Recommended cleanup:

- add `acc_cdf_from_cfg(...)`
- rename `acc_cdf_success_from_cfg(...)` to something explicit like `acc_cdf_given_start_from_cfg(...)` if it still has a legitimate use
- audit chained-onset branches and make sure they are using the intended version, not whichever helper happens to exist

### B. Replace observed likelihood with the exact transition-conditioned path

The observed path should not go through:

- `kernel_node_density_entry_idx(...)`
- then `density(target) * survival(competitor)`

for likelihood.

Instead, observed likelihood should become a one-event specialization of the ranked path.

There are two reasonable implementations.

#### Option 1: Reuse `ranked_prefix_density_resolved(...)` directly

For an observed trial:

- `ranked_outcome_indices = { outcome_idx }`
- `ranked_node_indices = { node_idx }`
- `ranked_times = { rt }`

and call [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:6160) `ranked_prefix_density_resolved(...)`.

This is mathematically clean and reuses the exact existing path.

Pros:

- smallest semantic change
- no new math
- one exact conditional architecture

Cons:

- it pulls in ranked scaffolding even for simple observed trials

#### Option 2: Extract an observed specialization from the ranked path

Add a new helper, for example:

- `observed_transition_density_resolved(...)`

that is literally the rank-1 case of [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:5992) and [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:6160).

It would:

- compile or reuse a `RankedTransitionCompiler`
- call `for_each_ranked_node_transition(...)` on the target outcome at time `t`
- for each emitted scenario, evaluate competitor survival under the emitted forced bits
- sum the scenario weights

This avoids the wrong marginal-competitor multiplication while keeping the implementation focused on observed trials.

Pros:

- exact
- cleaner than forcing every observed call through full ranked-prefix bookkeeping
- still reuses the same transition compiler and transition semantics

Cons:

- one extra wrapper/helper

### Recommendation

Option 2 is the best balance.

It keeps one exact conditional architecture for likelihood:

- ranked uses transition-conditioned evaluation
- observed becomes the rank-1 specialization of the same transition-conditioned evaluation

That is cleaner than keeping the current observed generic route alive for likelihood.

## What Should Change In The Code

### 1. Add unconditional event CDF helper

Add a helper next to the existing leaf helpers in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:703):

- `acc_cdf_from_cfg(...)`

Then replace the absolute-onset cdf branch in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:2195).

### 2. Add exact observed transition-density helper

Add a helper near the ranked machinery that does:

- resolve outcome/node
- compile target transition scenarios via `for_each_ranked_node_transition(...)`
- evaluate competitor survival under the emitted forced bits
- sum contributions

This helper should take:

- `ctx`
- `outcome_idx`
- `node_idx`
- `t`
- `trial_params`
- `trial_type_key`
- `kernel_runtime`
- `RankedTransitionCompiler`
- `IntegrationSettings`

and return the exact observed outcome density.

### 3. Route `TrialKernelMode::Observed` through that helper

Change [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:6528) `evaluate_sequence_density_kernel_idx(...)`.

Today, observed uses:

- `kernel_node_density_entry_idx(...)`

That should be replaced with the exact observed transition helper.

The current `guess_shortcut` branch should be re-evaluated after this change. Most likely it becomes redundant, because the exact transition path already runs through `eval_event_ref_idx(...)`, which already incorporates guess donor density.

### 4. Stop using the wrong generic competitor multiplication for likelihood

`node_density_with_competitors_from_state(...)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp:4486) can stay as a low-level generic helper if needed by other APIs.

But likelihood should no longer rely on it for observed trials.

That is the key architectural boundary:

- generic node density helper: optional low-level tool
- likelihood: exact transition-conditioned sequence density

## Tests

Add explicit regression tests for this model.

### Mathematical regression

At several times, compare component-conditioned densities against the hand-written formulas:

- `go_only -> S`
- `go_stop -> S`
- `go_stop -> X`

At minimum:

- `t = 0.25`
- `t = 0.30`
- `t = 0.40`
- `t = 0.50`

### Event semantics regression

Add a unit test that an absolute-onset event with `q = 1` has:

- density `0`
- cdf `0`
- survival `1`

The current code violates that cdf condition.

### Non-regression

Run the existing shared-trigger tests and ranked likelihood tests, especially:

- [tests/testthat/test_chained_onset_likelihood.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_chained_onset_likelihood.R)
- [tests/testthat/test_multi_outcome_likelihood.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_multi_outcome_likelihood.R)
- [tests/testthat/test_examples_models.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_examples_models.R)

## Expected Outcome

After the fix:

- `go_stop -> S` should match `q f_S(t) + (1 - q) f_S(t) S_stop(t)`
- `go_stop -> X` should continue to match
- absolute-onset event triples will be internally consistent for `q`

The important point is that this is not a trigger-state enumeration bug in isolation.

The correct diagnosis is:

- shared-trigger state enumeration is mostly fine
- absolute-onset event CDF semantics are wrong
- observed likelihood is using the wrong conditional structure

So the best fix is not “push triggers harder”.

The best fix is:

- make event CDF semantics correct
- make observed likelihood use the exact transition-conditioned path that already exists in ranked likelihood
