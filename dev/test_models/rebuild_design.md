# Rebuild Design

## Purpose

This document fixes the target architecture before implementation starts.

The old native code failed for two reasons:

1. model meaning was split across too many layers
2. runtime hot paths carried around complexity that should have been compiled away

So this is not a "rewrite the evaluator in cleaner C++" plan.
It is a compiler and semantics redesign.

The acceptance criterion is not that the code looks cleaner.
The acceptance criterion is:

- simple models run like simple models
- complex models are still exact
- there is one coherent design instead of accumulated exceptions

## Existing Surface We Must Preserve

The user-facing R DSL already commits us to the following features:

- accumulators with supported leaf distributions:
  - `lognormal`
  - `gamma`
  - `exgauss`
  - `LBA`
  - `RDM`
- leaf parameters:
  - distribution parameters
  - `q`
  - `t0`
- absolute onsets
- chained onsets via `after(accumulator)` and `after(pool)`
- pools, including nested `k-of-n`
- outcome expressions:
  - direct event
  - `first_of()`
  - `all_of()`
  - `none_of()`
  - `inhibit()`
- components:
  - fixed
  - sampled
  - component-specific `n_outcomes`
- triggers:
  - independent
  - shared
- observation modes:
  - top-1 observed outcome
  - ranked observed outcomes
  - missing later ranks as truncation
- outcome options used by top-1 models:
  - `map_outcome_to`
  - `guess`
  - censor/special outcomes

This surface is visible in:

- [dev/examples/new_API.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/examples/new_API.R)
- [vignettes](/Users/nstevenson/Documents/2025/AccumulatR/vignettes)
- [tests/testthat](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat)

The rebuild may tighten internals.
It should not silently weaken these semantics.

## Hard Rules

These are not suggestions.

1. There is one semantic definition of an outcome event.
2. Observed top-1 and ranked likelihood must use the same transition semantics.
3. Components are compile-time projections, not runtime masks.
4. Leaf helpers must separate unconditional event quantities from conditional-on-start quantities.
5. Simple projected variants must never pay for guards, ties, or trigger state if they do not contain them.
6. There may be at most two evaluator classes in the hot path:
   - direct evaluator
   - exact-transition evaluator
7. No fallback evaluator that keeps the old general complexity alive.
8. Hot numeric code must run on lowered runtime programs, not semantic trees.
9. Trial evaluation must batch over as many trials that share the same compiled path as practical.

If any implementation step violates one of these rules, the design has drifted.

## Core Semantic Model

### 1. Primitive event

Each accumulator defines a primitive event time `T_a`.

That event has:

- start probability
- onset rule
- completion-time distribution conditional on starting

The primitive leaf API must expose three unconditional functions:

- `pdf_event(t)`
- `cdf_event(t)`
- `survival_event(t)`

and optionally separate conditional-on-start functions:

- `pdf_given_start(t)`
- `cdf_given_start(t)`
- `survival_given_start(t)`

These are different objects.
They must not share names or be interchangeable.

This specifically prevents the `q`/CDF ambiguity described in
[dev/test_models/trigger_exact_fix.md](/Users/nstevenson/Documents/2025/AccumulatR/dev/test_models/trigger_exact_fix.md).

### 2. Pool event

A pool defines a new primitive event from child primitive events.

Examples:

- `1-of-n` pool: minimum child completion
- `k-of-n` pool: `k`-th order statistic

A pool is not an observation rule.
It is another event node in the semantic graph.

### 3. Outcome rule

An outcome is a boolean-temporal expression over primitive events.

The expression language is:

- `event(source)`
- `and(args...)`
- `or(args...)`
- `not(arg)`
- `guard(reference, blocker, unless = ...)`

The R DSL helpers map directly to this:

- `all_of(x, y)` -> `and`
- `first_of(x, y)` -> `or`
- `none_of(x)` -> `not`
- `inhibit(x, by = y)` -> `guard(reference = x, blocker = y)`

### 4. Observation

Observation is not part of outcome semantics.
It is a separate layer that asks questions like:

- top-1: which labeled outcome is observed first at time `t`?
- ranked: which ordered sequence of labeled outcomes is observed at times `t1 < t2 < ...`?
- NA/censor/guess/remap: how is the semantic outcome mapped to observed data?

This separation is required.
The old code mixed outcome meaning with observation shortcuts.

## Compiler Architecture

The compiler pipeline should be explicit and staged.

### Stage A: Normalize DSL

Input:

- `race_spec`

Output:

- canonical model graph

Responsibilities:

- normalize expression trees
- validate labels and references
- validate observation declarations
- validate onset dependencies
- validate parameter namespace

No performance work happens here.

### Stage B: Build semantic graph

Input:

- canonical model graph

Output:

- semantic model:
  - accumulators
  - pools
  - outcome expressions
  - trigger declarations
  - component declarations
  - observation spec

This graph is still component-agnostic and evaluator-agnostic.

### Stage C: Project per component

Input:

- semantic model
- component id

Output:

- projected semantic variant

This stage is mandatory.

Projection rules:

- inactive accumulator event -> impossible
- inactive pool members are dropped
- empty pool -> impossible
- singleton `1-of-1` pool -> child
- `and(..., impossible, ...)` -> impossible
- `and(x)` -> `x`
- `or(impossible, x, ...)` -> `or(x, ...)`
- `or(x)` -> `x`
- `guard(reference, impossible blocker)` -> `reference`
- `guard(impossible reference, blocker)` -> impossible
- `not(impossible)` -> always-true structural marker

After simplification:

- rebuild reachable outcomes
- rebuild competitor sets
- rebuild onset reachability
- drop dead accumulators and dead pools

This is the single most important structural step.

Without it, simple components will remain slow.

### Stage D: Classify variant

Each projected variant is classified into one of two evaluator classes.

#### Direct evaluator

Allowed when all of the following hold:

- no guards
- no positive-mass ties induced by shared completion structure
- no shared-trigger exact state expansion
- no ranked observation exact-conditioning requirement beyond direct order-statistic composition
- no observation rule that requires realized-transition conditioning

Typical examples:

- direct races
- direct pools
- direct ranked races without logical outcomes
- independent trigger start probabilities folded into leaves

#### Exact-transition evaluator

Required when any of the following hold:

- guard semantics
- same-time tie mass matters
- shared trigger states affect exact conditioning
- observed outcome density depends on realized transition source
- ranked observation on non-direct outcome expressions
- chained onset needs realized source completion time in a nontrivial logical graph

This classification happens at compile time from the projected variant.

### Stage E: Compile evaluator-specific plan

Direct variants compile to a minimal direct execution plan.
Exact variants compile to a transition enumeration plan.

There is no third architecture.

### Stage F: Lower to runtime program

This stage is separate from semantic projection and separate from evaluator choice.

Input:

- evaluator-specific plan

Output:

- runtime program objects with dense numeric layout

Responsibilities:

- replace string ids with integer slots
- replace semantic references with flat contiguous arrays
- split hot numeric fields into structure-of-arrays layouts
- precompute outcome/source lookup tables
- precompute parameter slot layout
- define the batch key used to group trials onto one kernel path

This is where speed is won or lost.

If evaluators still walk semantic trees or look up ids by string per trial, the design has failed even if the code is otherwise clean.

## Native Module Layout

There is no `src/` tree right now.
The rebuild should introduce one.

Recommended layout:

- `src/leaf/`
  - leaf distributions
  - unconditional and conditional-on-start leaf channels
  - onset helpers
- `src/semantic/`
  - semantic graph types
  - expression node types
  - observation types
- `src/compile/`
  - normalization
  - component projection
  - simplification
  - variant classification
  - per-variant plan builders
- `src/eval/direct/`
  - direct evaluator for leaf/pool/simple race variants
- `src/eval/exact/`
  - exact-transition evaluator
  - tie handling
  - trigger-state enumeration
  - ranked conditioning
- `src/runtime/`
  - context object
  - parameter binding
  - trial input binding
- `src/bridge/`
  - Rcpp entry points only

The bridge layer should stay thin.
No model logic belongs there.

## Core Native Data Structures

These are not final syntax.
They are the minimum conceptual objects the code should implement.

### Semantic graph

```cpp
enum class SourceKind : uint8_t {
  Accumulator,
  Pool
};

enum class ExprKind : uint8_t {
  Event,
  And,
  Or,
  Not,
  Guard,
  Impossible,
  TrueExpr
};

struct OnsetSpec {
  enum class Kind : uint8_t { Absolute, AfterAccumulator, AfterPool };
  Kind kind;
  int source_index;
  double lag;
  double absolute_value;
};

struct LeafSpec {
  std::string id;
  DistKind dist;
  OnsetSpec onset;
  int trigger_index;
  ParamBinding params;
};

struct PoolSpec {
  std::string id;
  int k;
  std::vector<int> member_source_indices;
};

struct ExprNode {
  ExprKind kind;
  int source_index;
  std::vector<int> children;
  int reference_child;
  int blocker_child;
  std::vector<int> unless_children;
};

struct OutcomeSpec {
  std::string label;
  int expr_root;
  OutcomeObservationMap obs_map;
};

struct SemanticModel {
  std::vector<LeafSpec> leaves;
  std::vector<PoolSpec> pools;
  std::vector<ExprNode> expr_nodes;
  std::vector<OutcomeSpec> outcomes;
  std::vector<ComponentSpec> components;
  std::vector<TriggerSpec> triggers;
  ObservationSpec observation;
};
```

### Projected variant

```cpp
struct VariantFeatures {
  bool has_guards;
  bool has_shared_triggers;
  bool has_ranked_observation;
  bool has_positive_mass_ties;
  bool needs_exact_conditioning;
};

enum class EvaluatorClass : uint8_t {
  Direct,
  ExactTransition
};

struct CompiledVariant {
  std::string component_id;
  SemanticModel projected_model;
  VariantFeatures features;
  EvaluatorClass evaluator_class;
  std::unique_ptr<DirectPlan> direct_plan;
  std::unique_ptr<ExactPlan> exact_plan;
};
```

### Runtime context

```cpp
struct CompiledModel {
  SemanticModel base_model;
  std::vector<CompiledVariant> variants;
  MixturePlan mixture_plan;
  ParameterLayout parameter_layout;
};

struct TrialData {
  int trial_id;
  int component_index; // -1 if latent
  ObservedTopK observed;
};

struct BoundTrial {
  const CompiledVariant* variant;
  BoundLeafParams leaf_params;
  BoundTriggerParams trigger_params;
  TrialData trial;
};
```

The hot evaluators should receive `BoundTrial`.
They should not receive raw R objects, raw DSL, or component masks.

### Lowered runtime program

The objects above are still too semantic for the hottest path.
Each evaluator class should lower further into a numeric program.

```cpp
struct RuntimeLayout {
  int n_leaves;
  int n_pools;
  int n_outcomes;
  int n_params;
};

struct TrialBlock {
  int variant_index;
  int start_row;
  int row_count;
};

struct DirectProgram {
  RuntimeLayout layout;

  std::vector<uint8_t> leaf_dist_kind;

  std::vector<int> onset_kind;
  std::vector<int> onset_source_index;
  std::vector<double> onset_lag;
  std::vector<double> onset_abs_value;

  std::vector<int> trigger_index;

  std::vector<int> pool_k;
  std::vector<int> pool_member_offsets;
  std::vector<int> pool_member_indices;
  std::vector<uint8_t> pool_member_kind;

  std::vector<int> outcome_source_index;
  std::vector<uint8_t> outcome_source_kind;
  std::vector<int> observed_label_index;

  ParameterLayout parameter_layout;
};

struct ExactProgram {
  RuntimeLayout layout;

  std::vector<uint8_t> leaf_dist_kind;
  std::vector<int> onset_kind;
  std::vector<int> onset_source_index;
  std::vector<double> onset_lag;
  std::vector<double> onset_abs_value;

  std::vector<int> trigger_kind;
  std::vector<int> trigger_member_offsets;
  std::vector<int> trigger_member_indices;

  std::vector<uint8_t> expr_kind;
  std::vector<int> expr_arg_offsets;
  std::vector<int> expr_args;
  std::vector<int> expr_ref_child;
  std::vector<int> expr_blocker_child;
  std::vector<int> expr_source_index;
  std::vector<uint8_t> expr_source_kind;

  std::vector<int> outcome_expr_root;
  std::vector<int> observed_label_index;

  ParameterLayout parameter_layout;
};
```

This is the level where SIMD and general throughput become feasible.
The evaluator should be fed flat arrays and integer indices, not nested objects.

### Trial scheduling

Trials should be grouped by the runtime path they actually use.

Primary grouping key:

- compiled variant index

Secondary grouping keys when needed:

- observed outcome pattern
- rank count
- censor/remap wrapper pattern
- fixed-vs-latent component status

The important point is that batching should happen above the numeric kernel, not inside it.
The kernel should see one program and one homogeneous block of trials.

That does not mean dozens of model-specific fast paths.
It means one generic batching mechanism feeding one direct kernel or one exact kernel.

## Leaf API

All leaf distributions must satisfy the same contract.

```cpp
struct EventChannels {
  double pdf;
  double cdf;
  double survival;
};

struct LeafDistribution {
  EventChannels event_channels(double t, const BoundLeafParams& pars) const;
  EventChannels given_start_channels(double t, const BoundLeafParams& pars) const;
  double sample_completion_time(const BoundLeafParams& pars, RNG& rng) const;
};
```

Rules:

- `event_channels` are unconditional with `q` included
- onset is applied outside the distribution-specific math by shared helpers
- `sample_completion_time` returns the completion time conditional on starting
- trigger handling is not duplicated inside each distribution

Independent triggers can be compiled into the leaf start probability.
Shared triggers cannot.

## Performance Principles

The rebuild should optimize for broad throughput, not for a zoo of clever exceptions.

Good principles:

- compile away structure before runtime
- batch trials by compiled path
- prefer flat integer-indexed arrays over pointer-heavy object graphs
- use structure-of-arrays for hot numeric fields
- keep string handling, maps, and expression interpretation out of the inner loop
- let the direct backend operate on one stable kernel shape across many simple models
- let the exact backend operate on one stable transition kernel shape across many complex models

Bad principles:

- one huge evaluator that interprets the model per trial
- special closed forms added one by one as model-specific branches
- runtime component filtering
- runtime id lookup
- outcome-by-outcome dispatch inside the hot kernel
- mixing observation wrappers into core numeric loops when they can be applied after

Honest tradeoff:

- maximum elegance is not "one evaluator and one runtime object"
- maximum elegance is "one semantic system, one compile pipeline, two backend families, and low-level programs with no semantic clutter"

If we do not lower far enough, we will get code that is abstract and reusable but slower than it should be.
If we lower too aggressively into bespoke cases, we will get fast-looking code that becomes unmaintainable.
The correct middle is:

- semantic generality at compile time
- dense regular kernels at runtime

## Evaluator APIs

### Direct evaluator API

```cpp
struct DirectEvalResult {
  double loglik;
  bool valid;
};

DirectEvalResult eval_top1_direct(
    const DirectPlan& plan,
    const BoundTrial& trial,
    const NumericOptions& opts);

DirectEvalResult eval_ranked_direct(
    const DirectPlan& plan,
    const BoundTrial& trial,
    const NumericOptions& opts);

ProbabilityTable eval_response_probabilities_direct(
    const DirectPlan& plan,
    const BoundLeafParams& params,
    const NumericOptions& opts);
```

This evaluator is allowed to use compact closed forms, pool order-statistic templates, and direct race formulas.
It is not allowed to know about guards or transition forcing.
It also should not walk semantic trees or inspect component structure per trial.

### Exact-transition evaluator API

```cpp
struct ForcedState {
  BitSet forced_complete;
  BitSet forced_survive;
  std::vector<RealizedTime> realized_times;
};

struct TransitionChoice {
  int outcome_index;
  int transition_node;
  double time;
  ForcedState post_state;
};

struct ExactEvalResult {
  double loglik;
  bool valid;
};

ExactEvalResult eval_top1_exact(
    const ExactPlan& plan,
    const BoundTrial& trial,
    const NumericOptions& opts);

ExactEvalResult eval_ranked_exact(
    const ExactPlan& plan,
    const BoundTrial& trial,
    const NumericOptions& opts);
```

The exact evaluator owns:

- transition enumeration
- post-transition forced-state propagation
- same-time tie handling
- competitor survival under realized conditions
- shared-trigger state expansion

This is the only place where those semantics live.
It should still operate on a lowered `ExactProgram`, not on semantic objects.

## Exact Semantics Strategy

This is the key design choice.

Top-1 and ranked likelihood should both reduce to repeated application of:

1. choose one realized transition at observed time `t`
2. derive the post-transition forced state
3. evaluate non-winning competitors under that realized state

That is the common semantic core.

Top-1 observed is rank-1 exact conditioning.
Ranked observed is repeated exact conditioning with updated state after each rank.

This is stricter than the old design, and that is the point.

## Tie Handling

Ties are not an afterthought.

The compiler must mark a variant as `has_positive_mass_ties` when the structure can produce them.

Examples:

- shared gate in `all_of(F, G)` and `all_of(S, G)`
- pooled completions with shared terminal source
- logically equivalent same-time completions through a shared event

The exact evaluator must then carry tie-resolution rules explicitly.

That means:

- identify tie-generating transition families
- enumerate tie-resolving order conditions
- integrate over the tied pre-times where necessary

This is non-negotiable because [tie.md](/Users/nstevenson/Documents/2025/AccumulatR/tie.md) is a concrete counterexample to marginal competitor survival.

## Trigger Strategy

### Independent trigger

Independent trigger is a leaf-level start probability.

It belongs in:

- leaf unconditional event channels
- simulator leaf start draw

Nothing more complicated is needed.

### Shared trigger

Shared trigger is not a leaf property.
It is a discrete latent state shared across multiple leaves.

It belongs above the leaves as a trigger-state expansion layer.

Compiler responsibilities:

- build trigger groups
- map affected leaves
- determine whether expansion is relevant for the projected variant

Exact evaluator responsibilities:

- enumerate trigger success/fail states
- bind leaf start availability under each state
- weight exact likelihood contributions by trigger-state probability

Direct evaluator may support shared trigger only when the projected structure remains simple enough for a closed direct treatment.
Otherwise classification must force the exact evaluator.

## Observation Wrappers

Observation mapping should be above evaluator output.

That includes:

- `map_outcome_to`
- guess redistribution
- censor/special outcome mapping
- missing later ranks as truncation

The evaluator should operate on semantic outcome labels.
Only then should wrappers convert semantic outcomes into observed outputs.

This prevents observation hacks from contaminating core event semantics.

## R Interface Contract

The public R API should stay essentially the same:

- `finalize_model()`
- `build_param_matrix()`
- `simulate()`
- `prepare_data()`
- `make_context()`
- `log_likelihood()`
- `response_probabilities()`

But the native boundary should become simpler.

Recommended native entry points:

```cpp
SEXP cpp_compile_model(SEXP prep);
SEXP cpp_loglik(SEXP compiled_ctx, SEXP params_mat, SEXP data_df, ...);
SEXP cpp_loglik_multiple(SEXP compiled_ctx, SEXP params_list, SEXP data_df, ...);
SEXP cpp_simulate(SEXP compiled_ctx, SEXP params_mat, SEXP trial_df, ...);
SEXP cpp_response_probabilities(SEXP compiled_ctx, SEXP params_mat, SEXP trial_df, ...);
```

One compiled context object.
No other native architecture visible to R.

## Implementation Order

This is the order I would actually follow.

### Phase 0: Freeze semantics

Deliverable:

- this design doc accepted
- decision on exact meaning of:
  - guard
  - none_of
  - tie resolution
  - shared trigger failure
  - NA/censor/guess/remap behavior

Exit gate:

- no unresolved semantic ambiguity remains for current tests/examples

### Phase 1: Leaf foundation

Scope:

- new `src/` tree
- leaf distributions
- unconditional event channels
- onset handling
- pool event helpers

Must support:

- pdf/cdf/survival for direct leaves
- simulation draws
- no triggers beyond leaf-local start probability

Tests:

- leaf distributions vs reference R formulas
- `q` semantics:
  - `cdf = (1-q) * cdf_given_start`
  - `survival = q + (1-q) * survival_given_start`
- absolute onset shift
- pool order-statistic checks

Benchmark gate:

- simple two-lognormal race probability and likelihood must be within 10% of the old fast path baseline once that baseline is re-established

### Phase 2: Component projection compiler

Scope:

- semantic graph builder
- component projection
- simplification
- variant classification

Tests:

- projected `go` component equals true simple model probabilities and likelihood
- inactive chained-onset source makes dependent accumulator impossible
- no dead nodes remain in projected variants

Benchmark gate:

- a known-component simple projected variant must be within 5% of the equivalent standalone simple model

This is the first hard performance gate.

### Phase 3: Direct evaluator

Scope:

- lowered direct runtime program
- trial blocking by compiled direct variant
- top-1 direct evaluation
- ranked direct evaluation for direct outcomes only
- direct response probabilities
- simulator using compiled variants

Tests:

- simple race
- pools
- nested pools where structure stays direct
- ranked direct outcomes
- independent trigger as leaf start probability

Benchmark gate:

- top-1 simple race faster than or equal to the old architecture under the same compressed workload
- direct kernels accept trial blocks without per-trial structural dispatch

### Phase 4: Exact-transition evaluator for top-1

Scope:

- lowered exact runtime program
- trial blocking by compiled exact variant
- transition enumeration
- forced-state propagation
- guard semantics
- exact competitor conditioning

Tests:

- [dev/test_models/trigger_exact_fix.md](/Users/nstevenson/Documents/2025/AccumulatR/dev/test_models/trigger_exact_fix.md)
- [tie.md](/Users/nstevenson/Documents/2025/AccumulatR/tie.md)
- guarded stop/change examples
- inhibitor-with-protector examples

Benchmark gate:

- exact variants can be slower than direct variants
- but top-1 guarded models must materially beat the old generic guarded path once rebenchmarked

### Phase 5: Shared triggers

Scope:

- shared trigger group compilation
- state expansion
- integration into exact evaluator

Tests:

- shared vs independent trigger examples
- shared trigger with top-1 direct models
- shared trigger with guarded exact models

Correctness gate:

- no separate algebraic special case for the known bug model
- generic exact semantics must make it pass

### Phase 6: Ranked exact conditioning

Scope:

- ranked exact transition conditioning
- missing later ranks as truncation
- multi-rank state propagation

Tests:

- existing ranked direct tests
- ranked with chained onset
- ranked with pool + trigger
- ranked with logical outcomes where supported

Benchmark gate:

- direct ranked models still use direct evaluator
- exact ranked only activates when structurally required

### Phase 7: Observation wrappers

Scope:

- `map_outcome_to`
- guess redistribution
- censor/special outcomes
- component-specific `n_outcomes`

Tests:

- existing example models with NA and guess semantics
- compression of repeated NA trials

Exit gate:

- current public vignette/examples/test semantics are restored on top of the new core

## Test Matrix

Each phase should add or preserve tests in three layers.

### 1. Semantic unit tests

Small exact checks against hand-written formulas.

Examples:

- simple race
- pool order statistic
- guard tie example
- trigger exact counterexample

### 2. Projection/compiler tests

Structural checks.

Examples:

- projected node counts
- no guards in simple projected variants
- no dead accumulators after projection
- classifier selects `Direct` vs `ExactTransition` correctly

### 3. End-to-end user tests

Examples from:

- [dev/examples/new_API.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/examples/new_API.R)
- [vignettes](/Users/nstevenson/Documents/2025/AccumulatR/vignettes)
- [tests/testthat](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat)

The point is to stop relying on only end-to-end tests for semantic regressions.

## Benchmark Policy

Benchmarks need to be explicit and attached to phases.

Recommended benchmark families:

- simple two-choice direct race
- simple `k-of-n` pool
- projected known-component simple variant
- guarded tie model
- stop/change model with shared trigger
- repeated NA/censor compressed workload

Rules:

- compare simple projected variant against equivalent standalone simple model
- compare exact variants against the best old baseline available, but never at the expense of simple variants
- fail the phase if a direct-eligible model is routed to exact evaluation
- fail the phase if hot loops still require semantic-tree walking or runtime id lookup

## What Must Not Be Reintroduced

- runtime component masking
- full `pdf/cdf/survival` computation for every node by default
- different mathematics for top-1 and ranked
- trigger-specific observed-density shortcuts
- tie-specific one-off routes outside the exact evaluator
- R-facing code deciding semantics that belong in native compilation

## Open Decisions

These still need explicit answers before coding the affected phase.

1. Exact support boundary for ranked logical outcomes:
   - all logical outcomes
   - or direct outcomes first, logical ranked later?
2. Whether `none_of()` remains a first-class semantic form or is always lowered to `not(event)` during normalization.
3. Whether shared-trigger direct variants deserve a dedicated closed-form direct plan when structure is otherwise trivial.
4. Whether simulation and likelihood should share the same projected variant plans exactly, or only the same semantic graph.

My recommendation:

- ranked logical outcomes can wait until after ranked direct outcomes are stable
- normalize `none_of()` early
- do not add a third evaluator for shared-trigger direct models unless profiling proves it is necessary
- simulation and likelihood should share the same projected variants, but not necessarily the same low-level plan object

## Recommended First Implementation Slice

The first slice should be:

1. create `src/` skeleton
2. implement leaf channel API
3. implement semantic graph types
4. implement component projection and simplification
5. prove that a projected simple component compiles to a true simple direct variant

Do not start with guards.
Do not start with ranked.
Do not start by porting the old kernel machinery.

If the projection/compiler foundation is wrong, everything built afterward will rot the same way the old code did.
