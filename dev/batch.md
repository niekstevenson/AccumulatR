# Batch Exact Evaluation Contract

This is the implementation contract for batched exact likelihood evaluation.
The goal is not to make scalar evaluation look vectorized. The goal is to stop
walking the same compiled exact probability program and source-product kernels
one parameter/data combination at a time.

The correct target is:

```text
group compatible evaluation requests
  -> execute compiled probability programs over active lanes
  -> execute source-product/integral kernels over contiguous active lanes
  -> reduce integral child lanes back to parent lanes
  -> write per-trial log-likelihoods
```

The important distinction is this:

- The vector DAG/tree walk is orchestration.
- The source-product/integral executor is the performance target.

A generic vector tree walk that still calls the scalar source-product machinery
is not success. That implementation was already tried and lost.

## Non-Negotiable Outcome

The implementation is worth keeping only if all of these are true:

1. `dev/validation/run_validation.R` passes with the current tolerances.
2. `dev/scripts/benchmark_speed.R` improves `stim_selective_stop` by at least
   20% against the scalar baseline.
3. No benchmark case regresses for a structural reason.
4. Profiles show reduced scalar observation-plan dispatch, scalar exact program
   dispatch, scalar compiled-node dispatch, and scalar source-product dispatch.
5. Remaining hot source-channel samples are real distribution/math work over
   contiguous active lane loops.
6. Integral-heavy cases do not do more source-channel work than scalar unless
   the extra work is mathematically unavoidable and measured.
7. Covered compiled shapes have one canonical hot path. Scalar remains a
   correctness oracle and an unsupported-shape path during development only.
8. Missing batch metadata is a compiler failure, not a runtime fallback.
9. Code size is justified by the benchmark and profile result.

Minimum useful target:

```text
stim_selective_stop, 50 trials, installed package, current quadrature:
  scalar baseline: about 34-36 ms/eval
  keep threshold: below 28 ms/eval
  preferred target: below 20 ms/eval
```

If a serious implementation cannot meet the keep threshold and explain the win
in the profile, delete it.

## Current Scalar Shape

The current exact finite-observation path is:

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

The useful scalar property is laziness. It evaluates only the branch, source
channel, and quadrature sample demanded by the active scalar path. Batch must
preserve that property.

The harmful scalar property is the unit of execution. The current source-product
integral kernel loops one parent request, one quadrature sample, one compiled
source-product op, and one source-channel call at a time. That repeats dispatch,
cache checks, time rebinding, and distribution calls across lanes that share the
same compiled shape.

## Rejected Implementations

### Eager Source Precompute

Do not collect broad source/time demands and fill dense source buffers before
evaluation.

That design is rejected because it:

- precomputed source/time combinations before branch gates knew which lanes
  remained active;
- made sparse nested integrals dense;
- created Cartesian products of parent and child quadrature nodes;
- moved cost into source precompute and distribution math;
- destroyed the scalar evaluator's useful laziness.

Observed evidence:

```text
case                         scalar/head   eager batch   disabled batch
stim_selective_stop          34.5 ms       746 ms        35.5 ms
example_16_guard_tie_simple  4.08 ms       4.08 ms       4.15 ms
example_6_dual_path          1.26 ms       1.32 ms       1.34 ms
stop_change_shared_trigger   0.14 ms       0.15 ms       0.15 ms
```

### Vector Wrapper Around Scalar Kernels

Do not wrap scalar exact evaluation in a lane loop. Do not build a vector tree
walk whose source-product leaf still calls scalar
`compiled_math_source_product_value_for_ops()` per lane.

A literal vector tree-walk implementation was validated and removed because it
failed the keep criteria:

```text
shape                                      stim_selective_stop
scalar/head                                34-36 ms/eval
full vector tree walk                      about 376 ms/eval
vector tree walk plus RT-free dedup         about 54 ms/eval
source-product op vector rewrite            about 92 ms/eval
```

The lesson is precise: correctness is not enough. A lane-vector wrapper around
the current scalar source-product machinery does not reduce the mathematical
work and adds setup/reduction overhead.

## Required Execution Model

### Lane Definition

A root lane is one evaluation request, not just one trial.

At minimum, a lane identifies:

```text
trial_index
component_code
variant_index
observation_state_code
observation_plan_id
probability_program_id
target_outcome_id when applicable
observed_rt when applicable
first_param_row or row_map layout
mixture/component weight slot when applicable
```

The executor may pack this more tightly, but it must not rediscover the same
information semantically at runtime.

### Grouping

Batch only compatible lanes together. Group keys must be compiled integer keys.

Required grouping dimensions:

- component code;
- exact variant id;
- observation plan id;
- probability program id/root id;
- outcome/state/rank shape;
- parameter row-layout compatibility;
- trigger-state enumeration shape;
- ranked/sequence state shape when ranked observations are involved.

Forbidden grouping behavior:

- string keys;
- semantic-name lookup;
- asking the scalar evaluator whether a row is supported;
- per-row heap allocation in the hot path;
- grouping that forces dense work for lanes that will be gated off immediately.

### Active Sets

Represent active work with compressed active-index arrays. Bit masks are allowed
as secondary metadata, but the hot source loops need compact lane lists.

Required:

```text
active_count
active_lane[i]
```

Reason: gates, zero-width integrals, failed bounds, impossible source
conditions, and zero products will make masks sparse. Iterating all lanes and
checking a boolean is not the intended hot shape.

## Vector Context

The vector context is a structure-of-arrays workspace.

Required logical fields:

```text
lane_count
active_lane[]
param_row[lane]
trial_index[lane]
component_weight[lane]
time_slot[slot][lane]
trigger_state_slot[slot][lane]
sequence_exact_time[source][lane]
sequence_upper_bound[source][lane]
expr_upper_bound[expr][lane]
expr_upper_normalizer[expr][lane]
used_outcome[outcome][lane]
node_value[node][lane]
parent_lane[child_lane]
quadrature_weight[child_lane]
cache_scope_id
```

Implementation can compress or omit fields for unsupported shapes, but a
covered shape must not use scalar workspace mutation shared across lanes.

Required invariants:

- active lanes are the only lanes that may trigger source math;
- child contexts know their parent lane and quadrature weight;
- nested integrals get distinct cache scopes when the same time slot has a new
  meaning;
- buffers are allocated from compiled workspace size estimates and reused;
- source-channel loops operate over contiguous arrays or compressed active
  lists;
- no per-lane maps exist in the hot path.

Forbidden:

- optional null workspaces in covered hot paths;
- vector growth because runtime discovered new structure;
- scalar `CompiledMathWorkspace` rebinding shared across lanes;
- treating inactive lanes as valid work for convenience;
- fallback to scalar evaluation on cache miss.

## Source-Product Kernel Contract

The source-product kernel is the main target. Its loop shape must be inverted
from scalar execution.

Current scalar shape:

```text
for lane:
  for source_product_op:
    value = scalar_source_value(op, lane)
    product *= value
```

Required batch shape:

```text
product[active_lane] = 1
current_active = active_lane

for source_product_op in compiled order:
  values = vector_source_value(op, current_active)
  for lane in current_active:
    product[lane] *= values[lane]
  compact current_active by dropping zero/impossible lanes
```

The op dispatch happens once per source-product op, not once per lane.

### Vector Source Reads

Source reads are demand-driven:

```text
source_value(source_product_channel, channel_kind, time_slot, active_lanes):
  key = source/channel/condition/time/cache_scope
  fill only missing active lanes
  return values for active lanes
```

Required behavior:

- fill only requested active lanes;
- preserve scalar laziness across gates and nested integrals;
- dispatch by distribution/source kind outside the lane loop;
- use one contiguous or compressed loop for one source/channel/condition;
- expose cache miss counts in development profiles;
- reject missing compiled source metadata before execution.

Initial distribution loops may still call scalar R math functions inside the
lane loop. That is acceptable only if source/op dispatch is already outside the
lane loop. SIMD is a later implementation detail.

Forbidden:

- precomputing every source channel for every possible quadrature time;
- filling inactive lanes;
- using only source/time ids as cache keys when nested integrals can rebind the
  same time slot;
- falling back to scalar source-channel evaluation on cache miss.

## Integral Kernel Contract

Integrals are vector control flow plus reduction. They are not an excuse for
dense precompute.

Required finite-integral shape:

```text
parent_active = active lanes entering integral
child_count = 0

for parent in parent_active:
  lower, upper = bounds(parent)
  if upper <= lower:
    result[parent] = 0
    continue

  for q in GL31:
    child_parent[child_count] = parent
    child_time[bind_time_slot][child_count] = map(q, lower, upper)
    child_weight[child_count] = mapped_weight(q, lower, upper)
    child_count += 1

child_active = 0..child_count-1
child_value = eval_integrand(child_context, child_active)

for child in child_active:
  result[child_parent[child]] += child_weight[child] * child_value[child]
```

Tail/no-response integrals use the existing GL47 tail rule with the same child
lane and reduction model.

Required:

- no child lanes for inactive parents;
- no child lanes for invalid or zero-width bounds;
- one vector integrand/source-product execution over child lanes when compiled
  shape is compatible;
- nested integrals repeat the same rule recursively;
- child contexts carry parent lane and quadrature weight explicitly;
- reduction metadata comes from compiled integral metadata;
- quadrature order is unchanged.

Forbidden:

- global dense parent x child quadrature precompute;
- special formulas for named benchmark models;
- hidden product-rule rewrites inside integral execution;
- changing quadrature method as part of batching.

## Compiled Metadata Contract

Compilation must decide whether a compiled probability program is covered by the
batch executor. It must not prescribe eager source buffers.

Required metadata:

- stable ids for observation plans, probability programs, roots, and source
  product kernels;
- grouping keys for compatible root lanes;
- node arity, value kind, result slot, child span, and schedule span;
- time-slot reads/writes;
- source-view and condition ids;
- trigger-state table shape;
- sequence/ranked state slot counts;
- integral metadata: lower slot, upper slot, quadrature rule id, bind time slot,
  integrand root/span id, child workspace estimate, and reduction target;
- source-product metadata: op span, source channel id, source program id,
  channel kind, fill mask, condition slots, time slot, time cap slot, source view
  id, and cache-scope id;
- maximum root lanes, maximum child lanes per integral level, and workspace
  size estimates;
- `batch_complete` flag for the exact covered shape.

Allowed metadata:

- static flags for dense finite-only shapes that can use simpler contiguous
  loops;
- debug-only scalar comparison hooks outside the hot path;
- development counters for cache misses and active-lane compaction.

Forbidden metadata:

- semantic string names in the hot path;
- global lists of source/time pairs to precompute;
- runtime discovery flags such as "maybe batch";
- metadata that requires scalar fallback after execution starts.

## Observation Plan Contract

Observation plans should be batch-scheduled, not scalar-interpreted per trial.

Required:

- group lanes by observation plan id and state code before execution;
- evaluate identical `LogDensity`, `FiniteOutcomeProbability`, and
  `NoResponseProbability` ops as batch requests;
- combine weighted sums over lane values after child probability programs finish;
- preserve RT-free parameter-block dedup semantics for missing-response plans;
- handle latent component choices as separate lanes with component weights, then
  reduce by trial using log-sum-exp.

Forbidden:

- one `evaluate_observation_plan_direct()` call per lane in the covered path;
- per-trial vectors for structural plan values;
- runtime semantic lookup of components, outcomes, or backends.

## Implementation Order

### 1. Baseline

Record validation, benchmark, and profile before changes:

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

Record:

- evals/sec and ms/eval for `stim_selective_stop`;
- benchmark rows for other exact-heavy cases;
- top profile frames and sample counts;
- git commit and dirty status.

### 2. Remove Dead Batch Paths

Remove or quarantine failed eager-precompute machinery before adding new batch
code.

Required:

- no dense source-buffer hot path remains;
- no disabled runtime switches imply a dormant batch engine;
- validation stays green.

### 3. Add Metadata Without Behavior Change

Add batch-completeness and workspace metadata to the compiled exact plan.

Success:

- scalar results are unchanged;
- unsupported shapes are classified before execution;
- no runtime fallback is added.

### 4. Build Grouping Without Behavior Change

Create deterministic integer-key grouping for root evaluation requests.

Success:

- group sizes are visible in a debug benchmark;
- scalar results are unchanged;
- grouping does not allocate per lane in the hot path.

### 5. Implement Vector Source-Product Kernel

Implement the inverted source-product loop and vector source reads first.

This is the first real performance step. If this still calls scalar
`compiled_math_source_product_value_for_ops()` per lane, the design has failed.

Success:

- op dispatch is outside the lane loop;
- source-channel dispatch is outside the lane loop;
- active lanes compact after zero/impossible factors;
- profiles show fewer samples in scalar source-product dispatch.

### 6. Implement Vector Integral Child Contexts

Implement finite and tail integral child-lane creation and reduction.

Success:

- invalid bounds produce no child lanes;
- nested integrals do not create global dense buffers;
- child contexts preserve scalar time-rebinding semantics exactly;
- representative per-trial probabilities match scalar.

### 7. Add Vector DAG/Observation Wrapper

Only after the source-product/integral kernel is real, widen the wrapper to
batch observation plans and compiled math nodes.

Success:

- scalar observation-plan interpretation falls in profiles;
- scalar compiled-node dispatch falls in profiles;
- simple arithmetic nodes do not dominate the implementation.

### 8. Enable Covered Shapes

Enable the batch executor only for `batch_complete` compiled shapes.

During development:

```text
batch_complete     -> batch executor
not batch_complete -> scalar executor, with documented deletion conditions
```

Final covered shapes must not keep equal-priority scalar and batch hot paths.

### 9. Consolidate

After validation, benchmark, and profile prove a win:

- remove superseded scalar hot-path code for covered shapes;
- collapse duplicate workspaces;
- delete development flags;
- delete failed scratch code;
- keep only maintained benchmark/profile scripts.

No cleanup means no accepted win.

## Accuracy Contract

Batch evaluation must preserve the current numerical method.

- Keep GL31 finite quadrature.
- Keep GL47 tail quadrature.
- Do not switch to Simpson/trapezoid as part of this work.
- Do not loosen validation tolerances.
- Compare representative per-trial probabilities while debugging, not only total
  log-likelihoods.

## SIMD Contract

SIMD is not the architecture.

SIMD becomes worth testing only after:

- source-channel loops operate over contiguous or compressed active lane arrays;
- dispatch is outside the lane loop;
- masks can be compressed or skipped cheaply;
- profiles have moved from scalar evaluator dispatch into distribution math.

If profiles still show scalar tree/node/source-product dispatch, SIMD work is
premature.

## Benchmark And Profile Contract

Every behavior-changing phase must run:

```sh
Rscript dev/validation/run_validation.R
ACCUMULATR_BENCH_INSTALLED=true Rscript dev/scripts/benchmark_speed.R
```

Hard-case profile:

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
evaluate_observed_trials_cached
evaluate_observation_plan_direct
evaluate_exact_probability_program
evaluate_exact_step_distribution
evaluate_compiled_node_span
evaluate_compiled_source_product_integral_kernel
compiled_math_source_product_value_for_ops
compiled_math_source_product_direct_leaf_scalar
Rf_pnorm/Rf_plnorm/Rf_dlnorm/log/exp
```

A benchmark-only win is insufficient. The profile must explain why the win
happened.

## Kill Criteria

Delete the batch work if any of these remain true after a serious implementation
attempt:

- validation requires weaker tolerances;
- `stim_selective_stop` speedup is below 20%;
- profile samples do not move away from scalar observation/exact/node/source
  dispatch;
- source-channel calls increase because inactive or unnecessary lanes are
  filled;
- nested or tail integrals become dense Cartesian precompute;
- code size grows while scalar exact paths remain the default for covered
  shapes;
- the implementation is stimulus-selective in disguise;
- simple benchmark cases regress structurally;
- missing batch metadata falls back at runtime.

Failed batch work is useful only if it is removed cleanly.

## Expected End State

```text
compiled exact plan
  owns grouping keys
  owns batch-completeness metadata
  owns source-product kernel metadata
  owns integral child/reduction metadata
  owns workspace size estimates

likelihood evaluation
  groups compatible root requests
  runs covered observation plans over lanes
  runs covered probability programs over lanes
  executes source-product kernels over active lane lists
  creates integral child lanes only when demanded
  reduces child lanes back to parents
  reduces component-choice lanes back to trial log-likelihoods
```

The final covered path is one coherent batch exact evaluator. Scalar execution
is not the optimized implementation; it is the oracle and the temporary path for
unsupported shapes.
