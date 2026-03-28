# Vectorization Speed Contract

Last updated: 2026-03-28

This note is ordered by speed impact only.

The standard is working code, lower wall-clock time, and less duplicated
runtime. Thin wrappers do not count.

## Non-Negotiables

1. No new adapter layer that preserves the old scalar outer loop under a batch
   API.
2. No milestone is complete if the hot path still re-enters the old execution
   body one param set, one state, one donor, one term, or one outcome at a
   time.
3. If a batched path is added, the old production path must be deleted or
   demoted out of hot execution in the same batch.
4. Every milestone needs both correctness proof and benchmark proof.
5. If the benchmark does not move materially, the milestone is not done.

## Required Proof For Every Milestone

Correctness:

- `Rscript dev/scripts/check_loglik_golden.R`
- relevant `testthat` coverage for the touched path

Speed:

- benchmark script for the target workload checked into `dev/scripts/`
- before/after numbers stored in `dev/scripts/scratch_outputs/`
- report at least wall-clock time and throughput

Failure conditions:

- no measurable throughput gain
- hidden scalar fallback remains in production hot path
- multiple runtimes still do the same work under different names

## Priority Order

1. Parameter-set batching
2. Ranked/exact scheduler collapse
3. Mixture, donor, alias, and outcome-mass outer-loop batching

## Milestone 1: Parameter-Set Batching

Target code:

- `src/distributions.cpp`: `cpp_loglik_multiple(...)`
- `src/distributions.cpp`: `cpp_loglik(...)`
- any prep/runtime code needed to reuse grouped trial structure across parameter
  sets

Current blocker:

- `cpp_loglik_multiple()` is still a scalar loop over full `cpp_loglik()` calls.

Required end state:

1. `cpp_loglik_multiple()` must stop calling full `cpp_loglik()` one-by-one.
2. Trial grouping, prepared trial structure, and batch-plan construction must
   be reused across parameter sets.
3. Parameter sets must enter the runtime as an added batch dimension, not as a
   repeated top-level restart.
4. Shared-trigger prep, grouped direct batches, ranked grouped batches, and
   outcome-mass grouped batches must all remain reachable in the multi-param
   path.

Hard reject if:

- `cpp_loglik_multiple()` still contains a loop whose body is effectively
  "call `cpp_loglik()` again"
- per-parameter-set trial grouping is rebuilt from scratch in production
- the multi-param path is just a new wrapper around the old single-param path

Success metrics:

1. `cpp_loglik_multiple()` no longer dispatches through repeated full
   `cpp_loglik()` calls.
2. Throughput on large `params_list` workloads improves by at least `2x`.
3. Single-parameter performance regression is less than `5%`.
4. Memory growth is bounded and explained; no unbounded per-param duplication of
   trial-group structures.

Benchmark workloads:

- optimization-style repeated direct trials
- repeated ranked trials
- repeated nonresponse/outcome-mass trials

## Milestone 2: Ranked/Exact Scheduler Collapse

Target code:

- `src/ranked_transitions.cpp`
- `src/ranked_transitions.h`
- `src/exact_outcome_density.cpp`
- `src/exact_outcome_density.h`

Current blocker:

- ranked still owns frontier scheduling
- exact still owns scenario scheduling
- exact still depends on `RankedBatchPlanner`

Required end state:

1. Ranked hot execution is reduced to lane grouping, state-delta application,
   and reduction. It must not behave like a separate runtime family.
2. Exact scenario handling must stop using ranked-specific planner machinery as
   its scheduler cache.
3. Scenario collection and reduction must move closer to the canonical vector
   runtime and away from private ranked/exact execution shells.
4. Ranked and exact must stop doing repeated regroup/rematerialize/collapse work
   when equivalent lane-group execution can do it once.

Hard reject if:

- `RankedBatchPlanner` remains the effective scheduler shell for exact
- ranked still walks frontier states in a way that is logically state-by-state
  execution with batch leaves attached
- the change just renames ranked/exact scheduler code without deleting it

Success metrics:

1. Complex ranked workloads improve by at least `30%` wall-clock time.
2. Exact-heavy guarded/shared-trigger workloads improve by at least `25%`
   wall-clock time.
3. Ranked and exact benchmark gains must come from less scheduler overhead, not
   from looser tolerances or reduced work.
4. Runtime-visible dependence of exact on `RankedBatchPlanner` is removed, or
   reduced to a builder-only artifact that does not drive hot execution.

Benchmark workloads:

- multi-rank shared-trigger likelihood
- guarded exact scenario likelihood
- ranked mixture likelihood

## Milestone 3: Mixture, Donor, Alias, And Outcome-Mass Outer-Loop Batching

Target code:

- `src/distributions.cpp`: `native_trial_mixture_internal_idx(...)`
- `src/distributions.cpp`:
  `accumulate_component_guess_density_batch_idx(...)`
- `src/distributions.cpp`:
  `accumulate_outcome_alias_density_batch_idx(...)`
- `src/distributions.cpp`: `mix_outcome_mass_batch_idx(...)`

Current blocker:

- mixture terms are still evaluated sequentially
- donor/source accumulation is still serial across donors
- outcome-mass still loops outcomes and components outside the batch leaf

Required end state:

1. Mixture terms must be grouped by reusable signatures and evaluated in batch.
2. Donor/source accumulation must batch repeated time/param patterns instead of
   looping donor-by-donor in the hot path.
3. Outcome/component reductions must happen after grouped batch evaluation, not
   inside nested scalar orchestration loops.
4. One-point "batching" inside a serial outer term loop does not count.

Hard reject if:

- `native_trial_mixture_internal_idx()` still loops terms and calls the density
  entry on `{t}` one term at a time
- donor/alias code still serializes all source work despite repeated identical
  point structure
- outcome-mass still pays nested outcome x component scalar orchestration in the
  dominant path

Success metrics:

1. Mixture-heavy workloads improve by at least `1.5x`.
2. Outcome-mass / nonresponse-heavy workloads improve by at least `25%`.
3. Native helper calls used by probability/mixture tooling show the same gain,
   not just `cpp_loglik()`.
4. No correctness drift in donor, alias, guess, or shared-trigger edge cases.

Benchmark workloads:

- fixed and free mixture models
- nonresponse-heavy models
- alias/guess-donor models with repeated trial structure

## What Does Not Belong In The Top 3

- direct-density cleanup that does not move benchmarks
- renaming runtime types without deleting runtime duplication
- stats/test vocabulary cleanup without hot-path change

## Exit Condition

The speed-first phase is done only when all of the following are true:

1. `cpp_loglik_multiple()` is not a scalar wrapper around `cpp_loglik()`.
2. Ranked and exact no longer pay dominant private scheduler-shell overhead.
3. Mixture/donor/outcome-mass helpers no longer fake batching with serial outer
   loops.
4. Benchmarks show real gains, and the gains survive correctness checks.
