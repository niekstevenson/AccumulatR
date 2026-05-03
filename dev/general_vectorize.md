# General Vectorized Evaluation Plan

This document is the implementation contract for replacing the current mixed
exact evaluator with a general precompiled vector executor.

The goal is not to collect more fast paths. The goal is:

```text
model/data prep resolves evaluation shape
  -> compiled plan contains route, cache, gate, bound, and leaf eligibility
  -> runtime walks plan steps once per compatible lane batch
  -> runtime executes compact lane arrays and distribution math kernels
  -> runtime reduces only numeric results
```

Anything else is not done.

## Non-Negotiable Requirements

1. Generality

- Every model/data shape not listed in [not_supported.md](not_supported.md) must
  compile into a vector plan or fail as a bug.
- Adding a new unsupported case is allowed only when the semantics are genuinely
  undefined or intentionally out of scope. It is not allowed because a shape is
  inconvenient for the vector executor.
- Unsupported cases must be declared in [not_supported.md](not_supported.md),
  rejected at the R API boundary when possible, rejected in C++ plan validation,
  and covered by tests.

2. Precomputation

- Runtime may read numeric parameters, observed data, state values, and scratch.
- Runtime must not rediscover semantic shape.
- Runtime must not ask "what kind of node/source/product/factor/gate is this?"
  inside per-lane or per-quadrature-child loops.
- Route kind, factor kind, cache slot use, gate requirements, bound strategy,
  source-program relation, and direct-leaf eligibility are plan fields.

3. Vectorization

- Covered evaluation is lane-batched by construction.
- `lane_count == 1` must use the same vector executor, unless the full shape is
  explicitly unsupported.
- A function that accepts arrays but calls the scalar evaluator per lane is a
  thin wrapper. It does not count.
- A function that dispatches per lane on node kind, source-program kind,
  distribution kind, or factor kind does not count.
- Compact active-lane arrays are the unit of work. Dense contiguous spans are a
  permitted specialization of the same representation, not a separate semantic
  path.

4. Efficiency

- Hot loops must be numeric loops over arrays, not semantic tree walking.
- Heap allocation and vector resizing are not permitted in hot loops.
- String lookup, R object lookup, and name matching are not permitted in hot
  loops.
- Source/channel/leaf availability checks must be hoisted into plan validation
  or parent/op preparation outside the repeated math loop.
- Cache computation must be demand-driven by compiled use masks. Caching
  everything because it is available is not efficient.

5. Fast Paths

- Model-specific fast paths are forbidden.
- Benchmark-specific fast paths are forbidden.
- A specialized route is allowed only when it is a general algebraic lowering
  rule with explicit eligibility and fallback behavior in the plan.
- Distribution-specific math kernels are allowed behind the common batch leaf
  kernel interface. Distribution-specific evaluator control flow is not.

## What Counts As Progress

Progress must satisfy all of these:

- It removes runtime semantic discovery from a measured hot path.
- It preserves or improves correctness under validation and golden tests.
- It has installed-package benchmark evidence, not only `pkgload::load_all()`.
- It has profile evidence showing the removed work actually disappeared.
- It deletes or makes unreachable the old duplicate path.

The following do not count:

- moving an `if` chain into a new helper while still running it in the hot loop;
- adding a second evaluator for one example model;
- adding a cache without a compiled use plan;
- accepting scalar fallback after vector execution has started;
- improving source-mode benchmarks while installed-package benchmarks regress;
- leaving both old and new long-term paths active for the same covered shape.

## Current Problems To Fix

The dirty prototype proves that plan-driven dispatch can help some examples, but
it is not the desired architecture.

- `exact_truth.hpp` is absorbing too much evaluator logic. The file now mixes
  planning assumptions, dispatch, vector storage, source-product evaluation,
  integral evaluation, and distribution kernels.
- The current direct-leaf, direct source-product-sum, flat quadrature, cached
  factor, and conditioned direct-leaf routes are partly general ideas but are
  implemented as parallel evaluator branches.
- Cached factor handling is too broad. Unique factors are cached eagerly, then
  runtime still rediscovers which lanes need each slot.
- Direct-leaf eligibility is detected in several places instead of being one
  compiled fact consumed by the executor.
- Gate and bound handling still contains repeated per-call analysis.
- Installed benchmarks show modest wins for some guard/shared-gate examples,
  but no win for `stim_selective_stop`. Source-mode benchmark wins are not valid
  evidence for production speed.

Before new speed work, the implementation must be reshaped so speed comes from
one general vector plan, not from accumulating more route-specific evaluator
copies.

## Target Architecture

### Files

Prefer a split like this:

```text
src/eval/vector_plan.hpp          immutable compiled vector plan structs
src/eval/vector_workspace.hpp     reusable SoA runtime storage
src/eval/vector_executor.hpp      generic vector plan executor
src/eval/vector_source_product.hpp source-product plan execution
src/eval/vector_integral.hpp      integral/factor plan execution
src/eval/vector_leaf.hpp          batch distribution kernel interface
src/eval/exact_planner.hpp        semantic -> vector plan compilation
src/eval/exact_truth.hpp          legacy bridge while migration is incomplete
```

Do not create files that merely wrap old functions. A new file must own a real
layer boundary.

### Vector Plan

The vector plan is immutable and model-specific. It may depend on compiled
observed-data layout where needed, but it must not contain per-call numeric
values.

Core structs:

```text
VectorProgramPlan
  node_plans
  source_product_plans
  source_program_plans
  integral_plans
  factor_plans
  gate_plans
  bound_plans
  leaf_plans
  schedules

VectorNodePlan
  node_id
  route_kind
  output_slot
  input_spans
  source_view_id
  time_id
  aux ids already lowered to typed plan references

VectorIntegralPlan
  bind_time_id
  domain_kind
  quadrature_rule_id
  integrand_plan_id
  factor_plan_span
  gate_plan_span
  upper_bound_plan_id
  reduction_kind
  clean_signed_source_sum
  route_kind

VectorSourceProductPlan
  op_span
  cache_slot_span
  product_init_kind
  clean_signed_source_sum

VectorSourceProductOpPlan
  op_kind
  value_channel_mask
  fill_channel_mask
  constant_value
  source_program_plan_id
  source_channel_id
  leaf_plan_id
  bound_plan_id
  direct_leaf_eligible
  conditioned_leaf_eligible
  time_mode
  cap_mode
  cache_result
  cache_slot

VectorFactorPlan
  factor_kind
  node_id
  integral_plan_id
  upper_time_id
  cache_slot
  use_plan_id
  clamp_probability

VectorGatePlan
  gate_kind
  outcome_bitset_or_id
  lhs_time_id
  rhs_time_id
  expr_upper_plan_id
  source_view_id
  open_mode

VectorBoundPlan
  source_id
  relation
  lower_terms
  upper_terms
  exact_terms
  static_lower
  static_upper
  static_exact
  uses_bind_time
  fallback_time_mode

VectorLeafPlan
  leaf_index
  dist_kind
  onset
  parameter_layout
  q_layout
  t0_layout
  channel_capability_mask
```

These are examples of the required shape, not exact names. The important rule is
that the executor consumes typed fields, not semantic variants.

### Route Kinds

Route kind is a compiled enum. Runtime switches once per node/plan step, never
once per lane.

Required route categories:

```text
Constant
TimeGate
Product
Sum
UnaryClamp
UnaryComplement
UnaryNegate
SourceProgram
SourceProduct
SourceProductSum
IntegralFinite
IntegralZeroToCurrent
IntegralAnalyticAntiderivative
IntegralQuadrature
ProbabilityProgram
TriggerProgram
RankedProgram
Unsupported
```

`IntegralAnalyticAntiderivative` is allowed only as a general algebraic lowering
rule, for example a source product that is exactly a constant scale times one
PDF integrated over a finite interval. It must not be tied to an example model.

`DirectLeaf` should not be a top-level model fast path. It should be an op/factor
eligibility property consumed by the generic integral/source-product executor.
If a direct-leaf tile kernel remains, it must be generated from
`VectorIntegralPlan + VectorSourceProductPlan` and must support every op shape
that declares `direct_leaf_eligible`.

### Cache Plans

Cache use must be compiled, not discovered.

For every source-product sum or integral plan, compile:

```text
cache_slots
  slot_id
  producer_factor_plan_id
  active_use_plan_id
  consuming_term_ids
  gate_dependency_summary
```

For every term/factor use, compile:

```text
factor_use
  factor_plan_id
  cache_slot
  required
  lane_filter_plan_id
```

Runtime is then:

```text
for cache_slot in plan.cache_slots:
  active = execute_use_plan(cache_slot.active_use_plan_id, parent_active)
  if active not empty:
    execute_factor_plan(cache_slot.producer_factor_plan_id, active, cache_values[slot])

for term in plan.terms:
  active = execute_term_gate_plan(term.gate_plan_id, parent_active)
  multiply compiled cache slots in declared order
  evaluate compiled source-product ops
  accumulate
```

Runtime must not scan all terms repeatedly to find whether a cache slot is used.
That scan is planning work.

### Gate Plans

Gate requirements must be lowered before evaluation.

Gate categories:

```text
AlwaysOpen
OutcomeAllowed
OutcomeUsed
TimeCompare
ExprUpperBefore
ExprUpperAfter
SequenceState
CompositeAnd
Unsupported
```

Every source-product term and integral factor receives a gate plan id. Runtime
executes that gate plan over active lanes and returns a compact lane span.

`CompositeAnd` is acceptable if it is a compiled list of gate plan ids. It is
not acceptable if it re-walks semantic expressions.

### Bound Plans

Conditioned source bounds must be compiled into a `VectorBoundPlan`.

Runtime may evaluate numeric bound values over lanes. Runtime may not discover
which bound terms exist or whether they reference bind time.

Required precomputed fields:

```text
has_lower
has_upper
has_exact
uses_bind_time
lower_term_span
upper_term_span
exact_term_span
static_source_view_relation
fallback_time_mode
```

If `uses_bind_time` is true and the consuming route cannot support dynamic
bind-time bounds, that route is ineligible at plan time. Runtime should not
discover this after entering a direct-leaf evaluator.

### Direct-Leaf Eligibility

Direct-leaf eligibility is a compiled property of each source-product op and
factor.

Required fields:

```text
direct_leaf_eligible
conditioned_direct_leaf_eligible
leaf_plan_id
time_mode: BindTime | ParentTimeSlot | ConstantTime | Unsupported
cap_mode: NoCap | BindTime | ParentTimeSlot | ConstantTime | Unsupported
bound_plan_id
bound_uses_bind_time
```

Rules:

- direct leaf eligibility cannot be recomputed in the executor;
- source-channel availability is validated once per parent/op preparation, not
  inside the quadrature child math loop;
- conditioned direct leaves use the same bound plan mechanism as generic
  conditioned source programs;
- ineligible ops still execute through the generic vector source-program path,
  not through scalar replay.

### Workspace

Runtime workspace is structure-of-arrays and reusable.

Required storage classes:

```text
Lane ids:
  active_lanes
  next_lanes
  parent_lanes
  child_lanes
  dense_range descriptors

Numeric arrays:
  time_slot[slot][lane]
  time_valid[slot][lane]
  node_value[node_slot][lane]
  product[lane]
  value[lane]
  weight[lane]
  bound_lower[lane]
  bound_upper[lane]
  bound_exact[lane]
  bound_has_exact[lane]
  source_value[lane]
  leaf_value[lane]
  cache_value[cache_slot][lane]
  scratch arrays by type

Pointers/ids:
  source_channels_by_lane
  source_view_by_lane
  eval_workspace_by_lane
  used_outcomes_by_lane
  leaf_input_by_lane
```

No vector growth is allowed after `ensure_capacity()` for a batch. If a code path
needs new storage, add it to the workspace and size it before hot execution.

### Batch Leaf Kernel Interface

Distribution-specific math is behind one interface:

```text
batch_leaf_values(
  leaf_plan,
  channel_mask,
  active_lanes,
  time_by_lane,
  lane_state,
  scratch,
  out_by_lane
)
```

Allowed distribution-specific code:

- lognormal PDF/CDF/survival formulas;
- gamma PDF/CDF/survival formulas;
- ex-Gaussian PDF/CDF/survival formulas;
- LBA PDF/CDF/survival formulas;
- RDM PDF/CDF/survival formulas.

Not allowed:

- separate evaluator-level lognormal-only direct-leaf control paths;
- benchmark-specific distribution branches;
- scalar R math callbacks in covered batch paths;
- per-lane distribution dispatch.

Dispatch by `dist_kind` once per op group is acceptable. Dispatch by
`dist_kind` per lane is not.

## Implementation Phases

### Phase 0: Establish Honest Baselines

Before changing code further:

1. Install dirty and clean builds with `R CMD INSTALL --preclean`.
2. Run installed benchmarks only.
3. Save source-mode benchmarks only as developer diagnostics.
4. Profile installed evaluation for hard cases.

Commands:

```sh
dev/scripts/benchmark_installed_compare.sh
```

Manual equivalent:

```sh
R CMD INSTALL --preclean -l /tmp/acc_dirty_lib .

ACCUMULATR_BENCH_INSTALLED=true \
R_LIBS=/tmp/acc_dirty_lib \
ACCUMULATR_BENCH_TRIALS=50 \
ACCUMULATR_BENCH_N_REP=5 \
ACCUMULATR_BENCH_OUT=dev/scripts/scratch_outputs/vector_baseline_dirty_n50.csv \
Rscript dev/scripts/benchmark_speed.R

ACCUMULATR_BENCH_INSTALLED=true \
R_LIBS=/tmp/acc_dirty_lib \
ACCUMULATR_BENCH_TRIALS=400 \
ACCUMULATR_BENCH_N_REP=3 \
ACCUMULATR_BENCH_OUT=dev/scripts/scratch_outputs/vector_baseline_dirty_n400.csv \
Rscript dev/scripts/benchmark_speed.R
```

Required record:

```text
git commit
dirty files
install library path
compiler flags
benchmark CSV
top profile frames
test results
```

Do not cite speedups without this record.

### Phase 1: Fix Correctness And Test Hygiene

Tasks:

- Replace or delete `dev/scripts/check_loglik_golden.R` if it keeps comparing
  against stale scratch fixtures.
- Move still-useful custom cases into the maintained golden registry in
  `tests/testthat/helper-loglik-golden.R`.
- Add a regression test for slicing prepared data or explicitly reject sliced
  prepared data with a clear error. A segfault is never acceptable.
- Add tests that every unsupported shape in [not_supported.md](not_supported.md)
  fails at the intended boundary.

Done means:

- `testthat::test_dir("tests/testthat")` passes;
- no known R-level call can crash by passing a data frame shape that should
  produce an R error;
- stale scratch fixtures are not part of correctness claims.

### Phase 2: Add Plan Introspection Before Rewriting

Add a debug/reporting API that dumps compiled evaluator coverage.

Minimum report fields:

```text
variant id
probability program id
node count by route kind
integral count by route kind
source-product op count by op kind
factor count by factor kind
cache slot count
cache slot use counts
gate plan count by gate kind
bound plan count by bound kind
direct-leaf eligible op count
conditioned-direct-leaf eligible op count
unsupported route count
```

This can be internal at first, but it must be callable from a dev script.

Done means:

- every benchmark example emits a coverage report;
- report totals are stable across repeated evaluations;
- no runtime evaluation is needed to compute the report.

### Phase 3: Define The Vector Plan Types

Create the vector plan structs and enums. Do not wire them into evaluation yet.

Tasks:

- Add route kind enums.
- Add factor kind enums.
- Add gate kind enums.
- Add bound strategy enums.
- Add source-product op plan structs.
- Add integral plan structs.
- Add plan validation.

Validation must check:

```text
all ids are in range
all spans are in range
cache slots refer to existing factor plans
cache slot consumers match term factor use
direct-leaf op has valid leaf plan
conditioned direct leaf has valid bound plan
unsupported route has a declared unsupported reason
no invalid semantic index reaches a hot plan field unless that field is optional
```

Done means:

- plan structs compile;
- plan validation runs during context creation;
- no evaluator behavior changes yet.

### Phase 4: Compile Source-Product Plans

Move source-product dispatch decisions from runtime into plan compilation.

For every source-product op compile:

```text
op_kind
channel masks
constant value
source channel id
source program plan id
leaf plan id
bound plan id
time mode
cap mode
direct leaf eligibility
conditioned direct leaf eligibility
cache result flag
cache slot id
```

For every source-product term compile:

```text
source-product plan id
gate plan id
factor use span
expr upper plan span
sign
cleaning mode
```

Remove duplicated runtime helpers that inspect source-program graphs to decide
whether an op is direct-leaf integrable. The compiler may use helpers; the
executor may not.

Done means:

- runtime source-product evaluation switches on compiled `op_kind`;
- direct-leaf eligibility is read from the op plan;
- no runtime code calls semantic graph inspection helpers for eligibility.

### Phase 5: Compile Cache Use Plans

Replace runtime term scans with compiled cache use.

Tasks:

- Compile unique factor slots only when at least one term benefits.
- Compile consuming term ids per slot.
- Compile a lane-use gate plan per slot.
- Compile factor producer route per slot.
- Compile term factor-use order.

Rules:

- A factor used once in an ungated term should usually be inline, not cached.
- A factor used by multiple terms, or by a term where the producer is expensive
  and lane reuse is expected, may be cached.
- A gated cache slot must compute only lanes consumed by open terms.
- Cache policy must be explainable from plan fields. No hidden heuristic in the
  executor.

Done means:

- `compiled_math_batch_factor_active_lanes_for_cache_slot()`-style term scans
  are gone from hot execution;
- cache slot computation is driven by compiled use plans;
- installed benchmarks do not regress hard cases.

### Phase 6: Compile Gate And Bound Plans

Move remaining gate and bound discovery into plan compilation.

Tasks:

- Compile outcome gates into bitset/id checks.
- Compile time gates into typed time comparisons.
- Compile expr-upper gates into expr-upper plan ids and mode.
- Compile sequence-state gates into explicit state ids.
- Compile conditioned source bounds into `VectorBoundPlan`.
- Precompute `uses_bind_time`.

Done means:

- runtime gate execution only evaluates numeric comparisons or compiled state
  checks over active lanes;
- runtime bound execution only fills numeric arrays from compiled term spans;
- no runtime gate/bound code walks source-product term structures to discover
  requirements.

### Phase 7: Replace Parallel Integral Branches With One Planned Executor

Current direct-leaf, direct source-product-sum, flat quadrature, cached-factor,
and generic quadrature branches should converge into one integral executor
driven by `VectorIntegralPlan`.

Permitted route lowering:

```text
IntegralAnalyticAntiderivative
IntegralQuadrature
IntegralQuadratureWithCachedFactors
IntegralQuadratureWithDirectLeafOps
IntegralQuadratureWithSourceProductSum
```

These are plan routes, not model fast paths. Eligibility is determined once by
the compiler.

Executor shape:

```text
prepare parent lanes
prepare child quadrature lanes or analytic interval lanes
execute compiled cache producers
for each compiled term:
  execute gate plan
  multiply compiled factor uses
  execute source-product plan
  reduce into child values
reduce child values to parent lanes
```

Done means:

- route selection happens once per integral plan;
- no branch revalidates direct-leaf support after entering execution;
- source-product-sum and source-product integrals share factor/cache/gate
  machinery;
- duplicate cached-factor dispatch blocks are removed.

### Phase 8: Normalize Direct-Leaf Handling

Direct-leaf handling must be a general property of source-product ops and
integral plans.

Tasks:

- Replace ad hoc conditioned-direct-leaf checks with compiled op fields.
- Use the same bound plan for conditioned direct leaves and generic conditioned
  source programs.
- Pre-resolve parent/op state once per parent lane group:
  leaf input pointers, source times, caps, static bounds, lower/upper CDFs where
  required.
- Keep distribution math behind the batch leaf kernel interface.

If a fused direct-leaf tile evaluator survives, it must consume only:

```text
VectorIntegralPlan
VectorSourceProductPlan
VectorSourceProductOpPlan[]
VectorLeafPlan[]
VectorWorkspace
```

It must not inspect semantic source-program structure.

Done means:

- direct leaf support is general across leaf distributions with batch kernels;
- unsupported direct-leaf subshapes fall back to generic vector evaluation, not
  scalar replay;
- no direct-leaf route is named after or coupled to benchmark examples.

### Phase 9: Remove Legacy Scalar And Duplicate Paths

After the planned executor covers a shape, delete the old path.

Delete or make unreachable:

```text
scalar source-product evaluators for covered shapes
scalar compiled integral evaluators for covered shapes
duplicate direct-leaf eligibility helpers in runtime files
duplicate cached-factor dispatch blocks
runtime term scans for cache slot usage
runtime semantic source-program graph inspection
```

Do not keep legacy code "just in case" unless it is behind an explicit
unsupported boundary and covered by tests.

Done means:

- a source grep for deleted scalar symbols confirms they are absent or only used
  in unsupported/error paths;
- coverage reports show covered examples do not touch scalar categories;
- validation and installed benchmarks pass.

### Phase 10: SIMD-Ready Array Kernels

Only after the plan/executor shape is clean, work on SIMD.

Tasks:

- Align scratch arrays where feasible.
- Track dense spans separately from arbitrary compact lane arrays.
- Add array primitive variants for dense spans:
  multiply, add, clamp, complement, weighted reduce, finite mask.
- Keep gather/scatter variants for sparse lanes.
- Add compiler-vectorized loops first.
- Use Accelerate/vForce only where benchmark and accuracy evidence justify it.

Rules:

- SIMD work must not change semantics.
- SIMD work must not introduce model-specific routes.
- SIMD work must not be mixed with semantic-dispatch cleanup in the same patch.

Done means:

- profiles show time moved into distribution math or vector primitives;
- array primitive tests cover dense and sparse active lanes;
- installed benchmarks improve or remain neutral.

## Acceptance Gates

### Correctness

Required before any merge:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE, helpers=FALSE); testthat::test_dir("tests/testthat", reporter="summary")'
Rscript dev/validation/run_validation.R
Rscript dev/scripts/benchmark_leaf_vectorization.R
```

Golden tests must use maintained fixtures. Scratch-output fixtures are not
correctness sources.

### Installed Benchmarks

Every meaningful speed claim must include installed-package runs.

Minimum benchmark matrix:

```text
n_trials: 50, 400
n_rep: 5 for 50, 3 for 400
cases: all benchmark_speed.R cases
builds: clean baseline, candidate
mode: ACCUMULATR_BENCH_INSTALLED=true
install: R CMD INSTALL --preclean
```

Required report columns:

```text
label
baseline_ms
candidate_ms
ratio
inner_reps
git commit
git dirty
install path
```

Source-mode benchmarks may be kept for iteration, but they cannot justify a
claim.

### Profiles

Hard-case profiles required:

```text
example_6_dual_path, 400 trials
example_16_guard_tie_simple, 400 trials
stop_change_shared_trigger, 400 trials
stim_selective_stop, 400 trials
```

Acceptance targets:

- no scalar evaluator frames for covered shapes;
- no runtime semantic graph inspection frames in hot loops;
- no repeated cache-slot term-scan frame;
- top frames are distribution math, quadrature numeric loops, vector array
  primitives, or unavoidable state transitions;
- non-math control frames must be explainable by compiled plan execution.

### Coverage Report

Every benchmark case must produce:

```text
unsupported_count = 0
scalar_covered_count = 0
runtime_semantic_lookup_count = 0
source_product_ops_planned = source_product_ops_executed
integral_factors_planned = integral_factors_executed
cache_slots_planned = cache_slots_executed_or_skipped_by_empty_use_plan
```

If these counters do not exist yet, the implementation is not mature enough to
claim general vectorization.

## Cleanup Policy For The Current Dirty Prototype

Keep only ideas that can be expressed as plan data:

- cached integral factor kind;
- cache slot plans;
- direct-leaf eligibility;
- conditioned direct-leaf relation;
- analytic PDF antiderivative as a general integral route;
- flat quadrature as a general route;
- parent/op preparation outside quadrature child loops.

Reject or rewrite:

- repeated runtime direct-leaf eligibility checks;
- repeated runtime cache-slot term scans;
- duplicate factor dispatch blocks;
- evaluator branches that cannot be described by `VectorIntegralPlan`;
- distribution-specific evaluator control paths not behind `batch_leaf_values`;
- broad "cache every unique factor" policy without use/cost planning.

The cleanup is not finished when the code is faster on one benchmark. It is
finished when the old duplicated paths are gone and the vector plan explains the
runtime work.

## Final Definition Of Done

The framework is done when all of these are true:

- All supported models compile to a vector plan.
- Unsupported models are exactly those in [not_supported.md](not_supported.md).
- Runtime evaluation never discovers semantic shape in hot loops.
- Runtime dispatch is per compiled plan step, not per lane.
- Source-product ops, factor uses, gates, bounds, and direct-leaf eligibility
  are plan data.
- Covered shapes have no scalar replay.
- Distribution math is reached through one batch leaf-kernel interface.
- Installed benchmarks are neutral or faster across the benchmark matrix.
- Hard-case profiles show time dominated by math and planned numeric control,
  not semantic dispatch.
- The code has one canonical path per covered shape.

Anything short of that is intermediate work, not completion.
