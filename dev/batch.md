# Batch Exact Evaluation Implementation Contract

This document is the implementation contract for turning exact likelihood
evaluation into a coherent lane-batched executor.

The goal is not "some batch helpers exist." The goal is:

```text
group compatible likelihood lanes
  -> walk each compiled probability/evaluation shape once per group
  -> execute source-product and integral work over active lane arrays
  -> call distribution math through lane-batch kernels
  -> reduce child lanes back to parent lanes
  -> write per-trial likelihoods
```

The final implementation must not be a collection of unrelated wrappers around
the scalar evaluator. Covered shapes get one canonical hot path. Unsupported
shapes may use scalar only through an explicit unsupported-shape boundary.

## Current State

This is the current mixed state after the fail-closed batch consolidation work.
It is materially better than scalar replay, but it is not fully SIMD/vectorized.

### Already Batched Or Grouped

- Observation requests are grouped into `ObservationBatchGroup`s for the
  non-identity observation path.
- Exact probability programs have batch execution entry points:
  `evaluate_exact_program_for_observation_batch()`,
  `evaluate_exact_probability_program_op_batch()`, and related trigger/state
  batch functions.
- Compiled math spans have batch entry points:
  `evaluate_compiled_node_span_batch()`,
  `evaluate_compiled_integral_node_batch()`, and
  `evaluate_compiled_integral_node_batch_lane_native()`.
- Covered compiled source and integral nodes no longer fall back to scalar
  evaluation inside the batch executor. If a source/integral fragment reaches
  the batch path without batch support, it throws instead of silently entering
  the deleted scalar compiled source-node or integral evaluator.
- Regular non-integral compiled nodes now dispatch once per scheduled node and
  then apply that operation across lanes. They no longer choose the node kind
  inside every lane iteration.
- Regular arithmetic compiled nodes now use the same batch-array primitive layer
  for constants, time gates, products, sums, probability clamps, complements,
  and negation.
- Finite integral kernels can execute as lane batches:
  parent lanes are filtered, quadrature child lanes are generated, child lanes
  evaluate the integrand in a grouped path, and child values are reduced back to
  parent lanes.
- Source-product op evaluation is grouped:
  `compiled_math_batch_source_product_value_for_ops()` dispatches once per op
  over a lane group. Product state is held in lane-indexed arrays, factors
  update an active mask, and product compaction is deferred to the term/block
  boundary instead of swapping compacted lane buffers after every op.
- Source-product product updates now route through a small batch-array
  primitive layer. Constant factors, direct leaf factors, cached source-program
  factors, expr-upper factors, integral factors, and parent reductions call
  shared array operations instead of open-coded product/update loops at each
  call site.
- Source-product program fills are grouped:
  conditioned fills, exact-gate fills, onset convolution fills, pool-k-of-n
  fills, and leaf fills run over lane arrays.
- The old hot scalar source-product symbols are no longer dominant in the
  stimulus-selective profile. The scalar source-product op evaluator,
  direct-leaf scalar callback, and scalar leaf-fill program have been removed
  from the active compiled evaluator source.
- Batch scratch and observation batch workspace are reused enough that previous
  allocation hot spots such as `std::vector<double>::__append` and `madvise`
  are no longer central in the stimulus-selective profile.
- Finite unranked identity observations now enter the same observation batch
  grouping and exact probability batch executor as non-identity observations.
- Ranked identity observations now enter grouped ranked observation batches.
  Transition roots and ready-expression normalizers are evaluated through the
  compiled batch root executor.
- `OnsetConvolution` source programs now handle both exact-onset and latent
  convolution cases in the batch path. Latent onset convolution uses grouped
  quadrature samples and calls the batch source-program executor for onset
  density.
- `PoolKOfN` source programs now use a lane-batch dynamic program for member
  CDF/survival/PDF composition instead of the scalar source-product program.

### Already Vectorized Distribution Math

- Lognormal batch leaf paths use vector log/exp and batched normal CDF:
  `compiled_math_batch_lognormal_leaf_values_from_times()`.
- Batched normal CDF uses a rational high-accuracy implementation with vector
  `vvexp()` for larger lane counts:
  `compiled_math_batch_normal_cdf_from_z()`.
- Ex-Gaussian batch leaf paths are present and route the normal-CDF and exp work
  through batched helpers. Lower-CDF/tail composition now carries per-lane
  parameter state through scratch arrays instead of reloading source parameters
  during each composition pass.
- Gamma PDF batch paths use vector log/exp for the density calculation.
- Gamma CDF/survival no longer calls `R::pgamma()` in the covered compiled
  leaf path. It uses an internal regularized-gamma implementation driven from
  the batch leaf executor. This is batch-native, but not a hardware SIMD gamma
  CDF.
- Batch leaf evaluation has one dispatcher:
  `compiled_math_batch_leaf_values_from_times()`. Source-program leaves, direct
  leaves, bind-time integral leaves, and top-1 probability leaves route through
  this dispatcher.
- The old scalar distribution leaf primitives have been removed from source:
  `lognormal_*_fast`, `gamma_*_fast_rate`, `exgauss_raw_*`, `lba_*_fast`,
  `rdm_*`, and `standard_leaf_channels*`.
- LBA batch leaf paths now have distribution-specific batch implementations for
  the PDF/CDF/survival formulas. Normal-CDF calls inside those formulas route
  through `compiled_math_batch_normal_cdf_from_z()`.
- RDM batch leaf paths now have distribution-specific batch implementations for
  the common and near-degenerate branches, so the covered compiled path no
  longer calls the scalar RDM leaf primitive. Normal-CDF and normal-log-CDF
  calls inside those formulas route through the batch normal helpers.

### Still Scalar Or Partially Scalar

- Regular non-integral nodes now write to batch-owned node-value storage rather
  than `CompiledMathWorkspace::values`. Arithmetic nodes still apply simple
  operations over lanes, but they no longer use scalar-shaped storage.
- Source-product bound resolution now resolves sequence bounds and condition
  terms into lane arrays. The old
  `CompiledSourceChannels::source_product_resolved_bound_values_from_time_slots()`
  entry point is gone.
- Expr-upper and timed normalizer handling now computes upper-bound time and
  normalizer arrays for active lanes. The old
  `compiled_math_batch_node_evaluator_for_lane()`,
  `compiled_math_batch_expr_upper_bound_for_node()`, and
  `compiled_math_batch_upper_bound_from_terms()` entry points are gone.
- Source-product control no longer compacts or gathers live product lanes after
  every op. Each source-product op evaluates over the original term lane span,
  updates the lane-indexed product/mask arrays, and defers product compaction to
  the term/block boundary.
- Some scalar C++ math remains in non-batch or unsupported paths. That is
  acceptable only while those shapes are explicitly unsupported by the batch
  coverage contract.

## Current Evidence

Installed-package benchmarks show the current state is consolidated but not a
clear speed win over the previous best batch snapshot:

```text
case                         current measured
stim_selective_stop, 400      50.00 ms/eval
stim_selective_stop, 50       27.00 ms/eval
example_6_dual_path, 400       8.00 ms/eval
stop_change_shared_trigger, 400 1.47 ms/eval
example_23_ranked_chain, 400   0.218 ms/eval
```

The batch-array primitive layer is therefore a structural consolidation, not a
speed win by itself. It creates one place to install dense-lane,
compiler-vectorized, or Accelerate-backed kernels, but the current
implementation still uses lane-index gather/scatter loops inside the
primitives. The full-span mask-first source-product op path is structurally
cleaner, but it is not a speed win on the stimulus-selective benchmark because
it now evaluates later op values for lanes that previous per-op compaction would
have skipped. The remaining speed work is inside the array kernels, masked value
evaluation APIs, and integral lane control.

Leaf-vectorization benchmark, run through production likelihood evaluation with
4000 two-choice trials:

```text
case                         current ms/eval  abs loglik error
gamma_pdf_and_survival            0.4691          9.09e-13
exgauss_pdf_and_survival          0.5000          9.09e-13
lba_pdf_and_survival              0.4691          0
rdm_pdf_and_survival              0.5250          9.09e-13
```

The benchmark is `dev/scripts/benchmark_leaf_vectorization.R`. It also checks
the engine log likelihood against independent R formulas.

Interpretation: removing `R::pgamma()` was a real local win in the earlier
before/after benchmark. The current Ex-Gaussian, LBA, and RDM implementations
are cleaner and remove scalar callback entry from the covered compiled path, but
the latest benchmark is a correctness/current-speed check, not evidence of a
production speedup. Do not cite these distribution cleanups as a production
speedup without a profile showing where the saved work appears.

Profile interpretation for `stim_selective_stop`, 400 trials:

- active source/integral batch execution no longer samples the deleted scalar
  compiled source-node or scalar compiled integral entry points;
- old scalar leaf callbacks are effectively gone from the hot profile;
- `erfc` is no longer the main wall;
- top frames are now batch leaf math and vector math:
  `compiled_math_batch_lognormal_leaf_values_from_times()`,
  `compiled_math_batch_normal_cdf_from_z()`, `VVLOG`, and `VVEXP`;
- remaining non-math cost is mostly source-product control, scratch
  `ensure_size()`, vector assignment/allocation noise, and ordinary batch
  lane-control loops.

This is not the finish line. It is a cleaner, fail-closed grouped evaluator
shape. The next speed work must attack the remaining lane-state and math
storage costs, not reintroduce scalar fallbacks.

## Non-Negotiable Design Rules

1. No model-specific fast paths.
2. No scalar replay hidden behind a batch-looking wrapper.
3. No eager dense source precompute across inactive branches or unused
   quadrature children.
4. No fallback to scalar after a covered batch shape has begun execution.
5. No runtime semantic lookup in hot loops.
6. No per-lane heap allocation in hot loops.
7. No loosened validation tolerances.
8. No new approximate distribution math without an error audit over parameter
   ranges, data ranges, tails, and summed likelihoods.
9. No duplicate long-term paths for the same covered shape.
10. Scalar `lane_count = 1` is allowed only as the same batch executor, or as an
    explicitly unsupported-shape boundary.

## Target Architecture

### Lane

A lane is one likelihood contribution request.

At minimum a lane carries:

```text
trial index
component code
variant id
observation plan id
probability program id
observed state/outcome/rank shape
observed time slots
parameter row
component weight
trigger/sequence/ranked state slots
```

The hot path must use integer ids and precompiled offsets. It must not use
strings or semantic names.

### Active Sets

Every batch executor receives compact active lane arrays:

```text
active_count
active_lane[active_count]
```

Boolean masks may exist as side metadata, but distribution math and
source-product loops operate over compact lane lists.

### Workspace

Workspace is structure-of-arrays and reused:

```text
time_slot[slot][lane]
time_valid[slot][lane]
node_value[node][lane]
source_leaf_value[lane]
product[lane]
weight[lane]
parent_lane[child_lane]
bound_lower[lane]
bound_upper[lane]
scratch_x[lane]
scratch_y[lane]
scratch_lane[active_count]
```

Covered hot paths must not mutate scalar `CompiledMathWorkspace` per lane.

### Distribution Kernels

Distribution-specific code is allowed only behind one batch leaf-kernel
interface. That means distribution kernels are specialized by distribution
family, not by benchmark model or evaluator context.

Required interface shape:

```text
batch_leaf_values(
  dist_kind,
  leaf_index,
  onset,
  active_lanes,
  time_by_lane,
  channel_mask,
  source_channels_by_lane,
  scratch,
  out_by_lane
)
```

Source-product program leaves, direct source-product leaves, and bind-time
integral leaves must all call the same distribution-kernel interface.

## Implementation Plan

### Phase 0: Freeze The Measured Baseline

Before further behavior changes, record these files:

```sh
Rscript dev/validation/run_validation.R

ACCUMULATR_BENCH_INSTALLED=true \
ACCUMULATR_BENCH_TRIALS=400 \
ACCUMULATR_BENCH_OUT=dev/scripts/scratch_outputs/benchmark_batch_contract_n400.csv \
Rscript dev/scripts/benchmark_speed.R

ACCUMULATR_BENCH_INSTALLED=true \
ACCUMULATR_BENCH_TRIALS=50 \
ACCUMULATR_BENCH_OUT=dev/scripts/scratch_outputs/benchmark_batch_contract_n50.csv \
Rscript dev/scripts/benchmark_speed.R

ACCUMULATR_PROFILE_R_SCRIPT=$PWD/dev/scripts/profile_workload_mixed.R \
ACCUMULATR_PROFILE_CASES=stim_selective_stop \
ACCUMULATR_PROFILE_TRIALS=400 \
ACCUMULATR_PROFILE_WORKLOAD_SECONDS=30 \
ACCUMULATR_PROFILE_DURATION=30 \
ACCUMULATR_PROFILE_INTERVAL=1 \
ACCUMULATR_PROFILE_TARGET_SIGNALS=1 \
ACCUMULATR_PROFILE_FILE=$PWD/dev/scripts/scratch_outputs/profile_batch_contract_stim_400.txt \
bash dev/scripts/profile_cpp_simple.sh
```

Record:

- git commit and dirty status;
- installed-package benchmark rows;
- eval-only profile top frames;
- counts for scalar symbols listed in the acceptance table below.

### Phase 1: Declare Batch Coverage Explicitly

Add a compiled coverage report for every exact variant and probability program.

Required categories:

```text
BatchComplete
BatchGroupedButScalarLeafMath
BatchGroupedButScalarBounds
BatchGroupedButScalarExprUpper
ScalarIdentityShortcut
Unsupported
```

`BatchGroupedButScalarBounds` and `BatchGroupedButScalarExprUpper` are retained
as zero-count report categories so stale regressions are visible. Covered paths
must not emit them.

Success criteria:

- coverage is determined before execution;
- logs/debug output can summarize counts by category;
- no runtime path asks the scalar evaluator whether it can proceed;
- validation remains 49/49.

### Phase 2: Consolidate The Evaluator Entry Points

Remove the architectural split where identity observations use a separate
scalar exact-trial shortcut.

Required:

- identity finite observations become lanes in the same observation/probability
  batch executor;
- the scalar `evaluate_exact_trials_cached()` implementation is deleted;
- scalar exact evaluation remains only for explicitly unsupported shapes;
- no semantic behavior changes.

Success criteria:

```text
example_6_dual_path profile:
  scalar source-product integral frames                 absent
  scalar source-product op frames                       absent
  batch executor frames explain the work
```

Benchmarks must not regress `example_6_dual_path`, `example_16_guard_tie`, or
`stim_selective_stop`.

### Phase 3: Consolidate Source-Product Leaf Calls

There must be one source-product leaf value path for covered shapes.

Required:

- source-product program leaf fill uses the batch leaf-kernel interface;
- direct source-product leaf values use the same interface;
- bind-time integral direct leaves use the same interface;
- scalar `batch_source_product_*_leaf_value()` callbacks are not used in
  covered batch hot paths;
- direct leaf availability checks are hoisted or made cheap enough to stay out
  of the math profile.

Success criteria:

```text
batch_source_product_lognormal_leaf_value      0 meaningful samples
batch_source_product_gamma_leaf_value          0 meaningful samples
batch_source_product_exgauss_leaf_value        0 meaningful samples
compiled_math_batch_program_leaf_values_for_kind 0 meaningful samples
```

### Phase 4: Vectorize Remaining Distribution Kernels

This phase is distribution math, not evaluator wrapping.

Required order:

1. Gamma CDF/survival.
2. LBA normal-CDF/PDF subexpressions.
3. RDM normal-CDF/log-CDF subexpressions.
4. Ex-Gaussian tail cases and small-lane thresholds.

Gamma CDF is the hardest correctness risk. Use one of these only:

- a high-accuracy vector implementation with an error audit;
- a stable library vector call if available;
- scalar `R::pgamma()` only while gamma CDF is marked
  `BatchGroupedButScalarLeafMath`, not `BatchComplete`.

Current state: the covered compiled path uses an internal regularized-gamma
implementation instead of `R::pgamma()`. Keep the benchmark/error audit in place
because this code is not delegated to R's gamma implementation anymore.

Error audit requirements:

```text
max absolute error
max relative error where reference is not near zero
signed mean error
sum error across realistic likelihood batches
tail error
monotonicity checks
validation models using the distribution
```

Success criteria:

- no validation tolerance changes;
- no systematic signed error over sweeps;
- profile moves time from scalar distribution symbols into vector math or
  unavoidable distribution arithmetic;
- benchmark does not regress non-target models.

### Phase 5: Batch Source Bound Resolution

Status: completed for covered source-product bound plans.

Deleted scalar-ish symbol:

```text
CompiledSourceChannels::source_product_resolved_bound_values_from_time_slots()
```

Required:

- resolve one bound plan over active lanes;
- use precompiled bound term offsets and time slots;
- write `lower`, `upper`, `exact`, and `has_exact` arrays;
- preserve time-valid semantics;
- no per-lane virtual/semantic source lookup.

Success criteria:

```text
source_product_resolved_bound_values_from_time_slots
  drops from visible hot profile to absent/tiny
```

### Phase 6: Batch Expr-Upper And Timed Normalizers

Status: completed for covered expr-upper factors and compiled expr-upper nodes.
Upper-bound time and normalizer values are produced as lane arrays.

Deleted scalar-ish symbols:

```text
compiled_math_batch_node_evaluator_for_lane()
compiled_math_batch_expr_upper_bound_for_node()
compiled_math_batch_upper_bound_from_terms()
```

Required:

- precompute source-view ids and subject ids in metadata;
- evaluate timed upper terms over active lanes;
- preserve sequence upper-bound precedence;
- write upper time and normalizer arrays;
- compact lanes after impossible upper-bound factors.

Success criteria:

- expr-upper symbols stop appearing as meaningful hot frames;
- guarded/tie validation remains green;
- `example_16_guard_tie_simple` improves or does not regress.

### Phase 7: Reduce Source-Product Control Overhead

Only after phases 3-6, optimize control overhead.

Allowed:

- predecode op spans into compact kernel records;
- split constant-only/product-only op spans;
- keep active lane arrays in scratch without repeated clears;
- avoid repeated direct-leaf availability checks inside quadrature loops when
  metadata proves availability.

Forbidden:

- model-specific rewrites;
- source/time precompute across branches;
- duplicate evaluator paths for the same shape.

Success criteria:

```text
compiled_math_batch_source_product_value_for_ops
  is small relative to distribution math and quadrature work
```

### Phase 8: Remove Transitional Scalar Paths For Covered Shapes

After coverage is explicit and benchmarks/profiles prove the batch path:

- delete or quarantine old scalar hot-path code for covered shapes;
- keep scalar only as unsupported-shape executor and test oracle;
- make covered shape scalar entry a failure in debug builds;
- remove development comparison hooks from hot paths;
- keep a small set of benchmark/profile scripts.

Success criteria:

- covered shapes have one hot path;
- scalar remains reachable only through an explicit unsupported category;
- code is smaller or at least structurally simpler than before consolidation.

## Acceptance Table

Every behavior-changing phase must report this table:

```text
case                          before ms/eval   after ms/eval   ratio
stim_selective_stop, 400
stim_selective_stop, 50
example_6_dual_path, 400
example_16_guard_tie, 400
stop_change_shared_trigger, 400
```

Every behavior-changing phase must also report profile counts for:

```text
evaluate_observation_plan_direct
scalar compiled source-product integral frames
scalar compiled source-product op frames
scalar compiled direct-leaf frames
batch_source_product_lognormal_leaf_value
batch_source_product_gamma_leaf_value
batch_source_product_exgauss_leaf_value
compiled_math_batch_program_leaf_values_for_kind
CompiledSourceChannels::source_product_resolved_bound_values_from_time_slots
compiled_math_batch_node_evaluator_for_lane
compiled_math_batch_expr_upper_bound_for_node
R::pgamma / Rf_pgamma
erfc
VVLOG
VVEXP
std::vector<double>::__append
madvise
```

Benchmark-only wins are not accepted. The profile must explain the win.

## Definition Of Done

The batch evaluator is done only when:

1. Validation passes without tolerance changes.
2. Covered shapes do not enter scalar exact/source-product/integral hot paths.
3. Scalar `lane_count = 1` works through the same batch executor where the shape
   is covered.
4. Distribution math is vectorized for lognormal, ex-Gaussian, gamma, LBA, and
   RDM, or the non-vectorized distribution is explicitly not `BatchComplete`.
5. Identity finite observations no longer bypass the batch evaluator.
6. Bound resolution and expr-upper/timed-normalizer work are lane-batched.
7. Allocation and vector growth are absent from the hot profile.
8. Installed-package benchmarks show no structural regressions.
9. The implementation has no model-specific paths.
10. The code has one maintained execution path per covered shape.

## Kill Criteria

Delete or revert a phase if any of these are true:

- validation needs weaker tolerances;
- the profile does not move time away from scalar dispatch or scalar math;
- inactive lanes are evaluated for convenience;
- nested integrals become dense parent x quadrature x source precompute;
- a change only helps `stim_selective_stop` while regressing broader exact-heavy
  cases;
- code adds another long-term path instead of replacing a path;
- benchmark improvements cannot be explained by the profile.
