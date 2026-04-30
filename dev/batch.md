# Batch Exact Evaluation Plan

This is the implementation contract for making exact likelihood evaluation
batch-oriented. The goal is not a thin vector wrapper around the scalar engine.
The goal is to stop redoing the same source-channel and quadrature work for each
observation when the compiled model already knows which math will be requested.

The first target is the hard exact finite-observation path exposed by the
stimulus-selective stop workload. The implementation must stay general: it should
optimize repeated `source_product(channel, time)` evaluation across any compiled
exact probability program that lowers into the existing exact math IR.

## Current Diagnosis

The scalar hot path is:

```text
evaluate_observed_trials_cached()
  -> evaluate_observation_plan_direct()
  -> evaluate_exact_probability_program()
  -> evaluate_exact_step_distribution()
  -> evaluate_compiled_node_span()
  -> evaluate_compiled_source_product_integral_kernel()
  -> compiled_math_source_product_value_for_ops()
  -> source PDF/CDF/survival calls
```

The profile says this is mostly native exact math and repeated distribution
channel evaluation. It is not mainly R overhead.

The important cost is repeated work:

- each observation enters a scalar observation-plan interpreter;
- each finite source-product integral maps the same Gauss-Legendre rule again;
- each integral sample re-enters source-product channel evaluation;
- distribution functions are called one scalar time at a time;
- equivalent channel/time requests are not shared across observations in a batch.

The fix must attack those facts directly.

## Non-Goals

- Do not add a second semantic evaluator.
- Do not add `if (batch_supported) else scalar` in the hot path and call it done.
- Do not wrap scalar `exact_loglik_for_trial()` in an outer vector loop.
- Do not preserve scalar and batch exact engines indefinitely.
- Do not special-case stimulus-selective stop formulas.
- Do not reintroduce product-rule algebra unless a profile proves it beats the
  batch source-channel approach.
- Do not change quadrature order as part of this work. Quadrature order is a
  separate accuracy policy decision.
- Do not optimize by weakening validation tolerances.

Temporary development switches are acceptable only in scratch scripts or local
experiments. Merged evaluator code must have one canonical execution path for a
compiled shape.

## Success Criteria

Batch exact evaluation is successful only if all of these are true:

1. `dev/validation/run_validation.R` passes with the same tolerances.
2. `dev/scripts/benchmark_speed.R` improves `stim_selective_stop` materially
   without regressing the rest of the benchmark suite outside noise.
3. `dev/scripts/profile_workload_mixed.R` shows fewer samples in scalar
   source-channel calls and finite source-product integral interpretation.
4. The profile explains the speedup. A faster wall clock without a changed hot
   path is not evidence.
5. The implementation removes or consolidates old scalar code for the covered
   exact path. The net codebase should become more coherent, not just larger.
6. The simple cases do not pay for batch machinery they do not need. Their
   compiled plans should contain fewer demands, not runtime filtering.
7. There is no runtime semantic discovery after the likelihood boundary.
8. There is no compatibility fallback in the evaluator.
9. Scratch benchmark files either become documented tools under `dev/scripts/`
   or are deleted.

Minimum useful performance target:

- `stim_selective_stop`, 50 trials, installed package, GL31 finite rule:
  reduce from about 35 ms/eval to below 20 ms/eval.

Stretch target:

- below 15 ms/eval while retaining validation and keeping code size controlled.

If the implementation cannot reach a clear speedup and explain it in the
profile, remove it.

## Core Design

Batch execution has three compiled layers and one runtime pass:

```text
context compilation
  -> observation batch plan
  -> exact time-demand plan
  -> source-channel demand plan

likelihood evaluation
  -> group rows by compatible compiled plan
  -> collect demanded times
  -> evaluate source channels into contiguous buffers
  -> execute exact probability programs over arrays
  -> write log-likelihood results
```

The important abstraction is a demand, not a model-specific formula.

```text
SourceChannelDemand {
  source_id
  condition_key
  channel_mask      // pdf, cdf, survival
  time_slot_id
}

TimeDemand {
  time_slot_id
  producer_kind     // observed_rt, bound_current_time, finite_quad_node, state_time
  parent_slot_id
  quadrature_node
}

BatchGroupKey {
  component_code
  exact_variant_id
  observation_plan_id
  outcome/state code
  probability_program_id
  row_layout_id
}
```

The exact field names can change. The invariants cannot:

- equivalent source/channel/time/condition requests share one buffer slot;
- plan compatibility is decided from compiled ids, not semantic names;
- row membership is gathered once per likelihood call;
- distribution channels are evaluated in contiguous arrays;
- exact programs read precomputed channel buffers when the compiled demand plan
  proves the requested value is available.

## Implementation Order

### 1. Freeze The Baseline

Before code changes, record:

```sh
Rscript dev/validation/run_validation.R
ACCUMULATR_BENCH_INSTALLED=true \
  ACCUMULATR_BENCH_OUT=dev/scripts/scratch_outputs/benchmark_speed_batch_baseline.csv \
  Rscript dev/scripts/benchmark_speed.R
ACCUMULATR_PROFILE_R_SCRIPT=$PWD/dev/scripts/profile_workload_mixed.R \
  ACCUMULATR_PROFILE_CASES=stim_selective_stop \
  ACCUMULATR_PROFILE_TRIALS=50 \
  ACCUMULATR_PROFILE_WORKLOAD_SECONDS=16 \
  ACCUMULATR_PROFILE_DURATION=12 \
  ACCUMULATR_PROFILE_INTERVAL=1 \
  ACCUMULATR_PROFILE_TARGET_SIGNALS=1 \
  ACCUMULATR_PROFILE_FILE=$PWD/dev/scripts/scratch_outputs/profile_batch_baseline.txt \
  bash dev/scripts/profile_cpp_simple.sh
```

Write down:

- evals/sec for `stim_selective_stop`;
- benchmark rows for other exact-heavy cases;
- top profile frames and sample counts;
- current dirty commit hash.

No batch implementation is complete without the matching after numbers.

### 2. Add Compile-Time Batch Metadata Only

Add metadata to the existing exact compiled plan. Do not change behavior yet.

Required metadata:

- stable ids for observation plans and probability programs;
- a compatibility key for rows that can execute together;
- a typed list of time demands for each exact probability program;
- a typed list of source-channel demands for each time demand;
- maximum batch workspace sizes derived from the compiled plan;
- a flag that the plan is batch-complete.

This phase must not add runtime fallback. If a plan cannot be batch-complete,
that is a compiler limitation to fix before enabling batch execution for that
shape.

Validation after this phase must be bitwise or near-bitwise identical because no
evaluation behavior changed.

### 3. Build Batch Grouping

Replace per-observation entry into exact finite evaluation with row groups.

Grouping must use compiled integer keys only:

- component code;
- exact variant id;
- observation plan id;
- probability program id or finite outcome/no-response program id;
- trigger/state code;
- parameter row-layout compatibility.

Forbidden:

- string keys;
- semantic-name lookup;
- rebuilding model structure;
- per-row heap allocation of records;
- grouping by asking whether a scalar path supports the observation.

The first version may group conservatively. Conservative grouping is acceptable
if it is compile-time structural and does not create a fallback path.

### 4. Collect Time Demands

For each group, collect all times required by its compiled demand plan:

- observed RT slots;
- current-time slots used by exact probability nodes;
- finite quadrature nodes mapped from compiled lower/upper slots;
- nested integral sample slots if the compiled plan contains them;
- state/ranked-chain time slots when present.

This must be a deterministic append into pre-sized workspace buffers. The plan
owns the maximum number of slots per row. Evaluation should not grow vectors
because it discovered more structure.

For GL31 finite integrals, the mapped nodes should be produced once per row per
integral bound, then reused by every source-channel demand that reads that node.

### 5. Evaluate Source Channels In Batches

Create a source-channel batch evaluator with a structure-of-arrays layout.

Required buffers:

```text
time[point]
param_row[point]
source_id[demand]
channel_mask[demand]
value[demand, point] or value[flat_demand_point]
```

The exact memory layout can differ, but it must be contiguous enough that the
distribution kernels can loop without scalar interpreter dispatch.

Initial implementation can call scalar R math functions inside a C++ loop, but
the loop must be over contiguous points for one source/channel. That creates the
place where later SIMD or approximate-vector math can be tested honestly.

Do not start with SIMD. First make the data layout vectorizable and prove it
cuts repeated work.

### 6. Execute Exact Programs Over Arrays

Teach the exact math evaluator to consume precomputed source-channel buffers
when executing a batch-complete program.

The evaluator should still execute compiled math nodes. It should not call back
into semantic expressions. The difference is that source-channel reads become
array loads instead of scalar distribution calls when the demand table contains
the requested channel/time/condition.

Integral kernels should become loops over precomputed quadrature point values:

```text
for row in group:
  sum = 0
  for node in finite_nodes:
    product = execute_source_product_from_buffers(row, node)
    sum += weight[node] * product
```

This is still a general source-product integral evaluator. It is not a
stimulus-selective formula.

### 7. Remove Covered Scalar Paths

Once batch exact finite evaluation validates and profiles correctly, delete or
fold the corresponding scalar-only code path.

Required cleanup:

- remove duplicate source-channel leaf evaluation branches that are no longer
  reached for batch-complete exact plans;
- remove development flags and local comparison switches;
- collapse duplicated workspaces into one batch-aware workspace;
- keep scalar helpers only when they are used by non-batch plan compilation,
  debug checks, or genuinely scalar model classes still outside the phase scope;
- document any remaining scalar path with a specific unsupported shape and a
  deletion condition.

The final state should not have two equal-priority exact engines.

## Workspace Contract

Batch workspaces are part of the compiled plan contract.

Allowed:

- pre-sized numeric arrays derived from compiled maximums;
- reusable per-call buffers;
- flat spans into shared arrays;
- debug assertions that compiled ids are valid.

Forbidden:

- per-observation maps;
- per-observation source-product vectors;
- demand vectors that grow because the evaluator discovered structure;
- string-keyed caches;
- fallbacks on missing demand entries;
- optional null workspaces in the hot path.

If a required source-channel demand is missing, that is a compiler bug. The hot
path should not recover from it.

## SIMD Position

SIMD is a second-order optimization. The batch design must make SIMD possible,
but SIMD is not the definition of success.

SIMD becomes worth testing only after:

- source-channel times are contiguous;
- one loop evaluates many points for the same source and channel;
- branchy semantic dispatch is outside the inner loop;
- profile samples move from interpreter overhead into math kernels.

If the profile still shows scalar program dispatch dominating, SIMD work is
premature.

## Accuracy Contract

Batch evaluation must preserve the existing numerical method unless a separate
accuracy decision is made.

- Keep GL31 finite quadrature while implementing batch.
- Keep GL47 tail quadrature.
- Do not switch to Simpson/trapezoid as part of batch work.
- Do not loosen validation tolerances.
- Do not compare only log-likelihood totals; compare validation model checks and
  representative per-trial probabilities when debugging.

After batch is correct, quadrature-order policy can be reconsidered separately.
GL25 is a candidate based on profiling experiments, but it is not part of this
plan.

## Benchmark Contract

Every phase that changes evaluation behavior must run:

```sh
Rscript dev/validation/run_validation.R
ACCUMULATR_BENCH_INSTALLED=true Rscript dev/scripts/benchmark_speed.R
```

For hard-case profiling:

```sh
ACCUMULATR_PROFILE_R_SCRIPT=$PWD/dev/scripts/profile_workload_mixed.R \
ACCUMULATR_PROFILE_CASES=stim_selective_stop \
ACCUMULATR_PROFILE_TRIALS=50 \
ACCUMULATR_PROFILE_WORKLOAD_SECONDS=16 \
ACCUMULATR_PROFILE_DURATION=12 \
ACCUMULATR_PROFILE_INTERVAL=1 \
ACCUMULATR_PROFILE_TARGET_SIGNALS=1 \
ACCUMULATR_PROFILE_FILE=$PWD/dev/scripts/scratch_outputs/profile_batch_after.txt \
bash dev/scripts/profile_cpp_simple.sh
```

Required comparison table:

```text
case                       before ms/eval   after ms/eval   ratio
stim_selective_stop
example_16_guard_tie_simple
example_6_dual_path
stop_change_shared_trigger
```

Required profile comparison:

```text
frame/symbol                                      before samples   after samples
evaluate_observation_plan_direct
evaluate_compiled_source_product_integral_kernel
compiled_math_source_product_value_for_ops
compiled_math_source_product_direct_leaf_scalar
Rf_pnorm/Rf_plnorm/Rf_dlnorm/log/exp
```

A benchmark-only win is insufficient. The profile must show why the win happened.

## Kill Criteria

Remove the batch work if any of these remain true after a serious implementation
attempt:

- validation requires weaker tolerances;
- stimulus-selective speedup is below 20%;
- source-channel samples do not fall in the profile;
- code size grows substantially while scalar exact paths remain the default;
- the implementation is stimulus-selective in disguise;
- simple benchmark cases regress for structural reasons;
- the evaluator contains fallback logic to hide incomplete batch plans.

This is intentionally harsh. A failed batch experiment is useful only if it is
deleted cleanly.

## Expected End State

The clean end state is:

```text
compiled exact plan
  owns observation batch keys
  owns time-demand spans
  owns source-channel-demand spans
  owns workspace sizes

likelihood evaluation
  groups compatible rows
  fills demanded time buffers
  fills source-channel buffers
  executes compiled exact programs over rows
  accumulates log-likelihood
```

The scalar mental model should disappear from the exact finite-observation hot
path. What remains is one compiled evaluator that happens to run one row or many
rows depending on the group size.

