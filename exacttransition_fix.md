# ExactTransition Fix Plan

## Purpose

The current `ExactTransition` evaluator is not general enough for models
where an event that happened before the observed response time remains
relevant as a time boundary for a competing outcome.

The failure is now explicit in `dev/validation/cases.R`:

- `stop_change_shared_trigger`
- `stim_selective_stop`

These are not trigger-specific or stop-model-specific failures. They
expose a missing layer in the exact evaluator.

The fix must be:

1.  general across
    [`all_of()`](https://niekstevenson.github.io/AccumulatR/reference/all_of.md),
    [`first_of()`](https://niekstevenson.github.io/AccumulatR/reference/first_of.md),
    [`inhibit()`](https://niekstevenson.github.io/AccumulatR/reference/inhibit.md),
    pools, chained onsets, shared triggers, and ranked observations
2.  efficient for simple and moderately complex models
3.  implemented as a clean extension of `ExactTransition`, not as a
    wrapper or a set of model-specific formulas

## Current Problem

The current exact evaluator mostly represents scenario facts as:

- source happened before observed time
- source happened at observed time
- source survived past observed time

In code this is the `ExactRelation` / `RelationView` layer:

``` cpp
enum class ExactRelation {
  Unknown,
  Before,
  At,
  After
};
```

That is too coarse.

`Before(t)` discards the actual completion time. For many models that is
fine. For the new failing cases it is not.

### Stop-Change Counterexample

For:

``` r
S = inhibit(S, by = stop)
X = all_of(change, stop)
```

observing `X` at time `t` has two transition paths:

``` text
1. change = t, stop = u < t
2. stop = t, change < t
```

Path 1 cannot be reduced to:

``` text
stop Before t
```

because the competing `S` outcome wins if:

``` text
S < stop
```

The exact value of `stop = u` is the boundary. Collapsing it to
`Before(t)` loses the information needed to compute competitor non-win
probability.

The correct density contains:

``` text
f_change(t) * integral_0^t f_stop(u) * S_S(u) du
```

not:

``` text
f_change(t) * F_stop(t)
```

### Stim-Selective Counterexample

The stimulus-selective stop model has visible outcomes:

``` r
A = first_of(
  inhibit(A, by = S1),
  all_of(A, S1, inhibit(IS, by = S2))
)
```

and similarly for `B`.

The wrong marginal target is:

``` text
f_A_out(t) * S_B_out(t)
```

That treats `A_out` and `B_out` as independent marginal outcome
processes. They are not. They share `S1`, and the late branch also
shares the ignore-gate state. Correct evaluation must keep the branch
path and its latent ordering.

## Design Principle

The exact evaluator must stop treating composite outcomes as opaque
`pdf/cdf/survival` triples when competitor conditioning matters.

Instead, `ExactTransition` should compile outcome expressions into
temporal transition scenarios:

``` text
observed time t
latent time variables u1, u2, ...
density factors
cdf/survival factors
inequality constraints between times
tie/readiness constraints
trigger-state constraints
```

Then competitor non-win is computed by merging target and competitor
temporal scenarios, not by multiplying target density by marginal
competitor survival.

This is the missing framework layer.

## Keep The Existing Architecture Split

Do not replace the direct evaluator.

The desired architecture remains:

- `Direct`: for simple projected variants that can be evaluated with
  closed products of independent channels
- `ExactTransition`: for variants whose semantics require
  path-conditioned transition evaluation

The fix is inside `ExactTransition`.

Simple variants must not pay for temporal scenario machinery.

## New ExactTransition IR

Add a small temporal constraint IR. Suggested location:

``` text
src/eval/exact_temporal.hpp
```

The names below are illustrative, but the concepts should be explicit.

### Time Variables

``` cpp
enum class ExactTimeVarKind : std::uint8_t {
  Observed,
  PreviousRank,
  Latent
};

struct ExactTimeVar {
  ExactTimeVarKind kind;
  semantic::Index index;
};
```

Required built-in variables:

- observed response time `t`
- previous ranked response times, when evaluating ranked observations
- latent variables introduced by target or competitor transition
  scenarios

### Time Expressions

Keep these intentionally small at first:

``` cpp
enum class ExactTimeExprKind : std::uint8_t {
  Zero,
  Infinity,
  Observed,
  Variable
};

struct ExactTimeExpr {
  ExactTimeExprKind kind;
  semantic::Index var_index;
};
```

Do not add symbolic algebra unless a real case requires it. Current
needs are mostly `0`, `t`, `Inf`, and latent variables.

### Constraints

``` cpp
enum class ExactTimeRelation : std::uint8_t {
  Less,
  LessEqual,
  Equal,
  Greater,
  GreaterEqual
};

struct ExactTimeConstraint {
  semantic::Index lhs_var;
  ExactTimeRelation relation;
  ExactTimeExpr rhs;
};
```

Use strict inequalities conceptually, but recognize that most continuous
boundaries have zero mass. Ties with positive mass should not rely on
strict floating comparisons; they should be represented by explicit
tie/readiness logic.

### Source Factors

``` cpp
enum class ExactTemporalFactorKind : std::uint8_t {
  Pdf,
  Cdf,
  Survival
};

struct ExactTemporalFactor {
  ExactSourceKey source;
  ExactTemporalFactorKind kind;
  ExactTimeExpr time;
};
```

Examples:

``` text
f_stop(u)
S_S(u)
F_change(t)
```

### Expression Factors

Some factors are composite and should initially remain callable through
a compiled sub-plan:

``` cpp
struct ExactTemporalExprFactor {
  semantic::Index expr_idx;
  ExactTemporalFactorKind kind;
  ExactTimeExpr time;
};
```

The compiler should only leave an expression factor in place when it is
safe to evaluate marginally. If the expression’s internal times can
affect a competitor, it must be expanded to temporal terms instead.

### Temporal Terms

``` cpp
struct ExactTemporalTerm {
  std::vector<ExactTimeVar> vars;
  std::vector<ExactTimeConstraint> constraints;
  std::vector<ExactTemporalFactor> source_factors;
  std::vector<ExactTemporalExprFactor> expr_factors;

  semantic::Index active_source_id;
  ExactTimeExpr transition_time;

  bool has_readiness;
  ExactTimeExpr readiness_time;

  double fixed_weight;
};
```

A term represents one disjoint path contribution.

For example, the stop-change `X` path where `change = t` and
`stop = u < t` should compile roughly to:

``` text
vars:
  u_stop

constraints:
  0 < u_stop < t

factors:
  f_change(t)
  f_stop(u_stop)
  S_S(u_stop)

fixed_weight:
  1 - q_trigger
```

The important part is that `S_S(u_stop)` is possible. The old relation
overlay cannot express that.

## Compiler Changes

### 1. Keep Existing Projection

Do not change component projection as part of this fix unless a bug is
found there.

The pipeline remains:

``` text
race_spec
-> semantic model
-> projected semantic variant
-> direct or exact runtime plan
```

The temporal IR is built only for exact variants.

### 2. Add A Temporal Scenario Compiler

Add a compiler that can build temporal terms for these query forms:

``` text
density(expr at t)
cdf(expr by t)
survival(expr past t)
```

The implementation should not walk R-side structures. It should consume
the lowered exact runtime program.

Suggested APIs:

``` cpp
ExactTemporalPlan compile_temporal_plan(const ExactVariantPlan &plan);

std::vector<ExactTemporalTerm> compile_transition_terms(
    const ExactTemporalPlan &plan,
    semantic::Index expr_idx,
    ExactTimeExpr transition_time);

std::vector<ExactTemporalTerm> compile_cdf_terms(
    const ExactTemporalPlan &plan,
    semantic::Index expr_idx,
    ExactTimeExpr upper_time);
```

`survival(expr, t)` should be represented as complement of CDF through
the competitor union machinery, not by prematurely expanding `1 - ...`
into black-box marginal survival when competitor conditioning is active.

### 3. Event Source Rules

For a leaf/source event:

``` text
density(source at x) -> f_source(x)
cdf(source by x)     -> F_source(x)
survival(source x)   -> S_source(x)
```

When source time must be retained, `cdf(source by x)` can instead expand
to:

``` text
integral over u:
  f_source(u)
  0 < u < x
```

The compiler should decide whether to keep `F_source(x)` collapsed or
introduce `u`. That decision is the central optimization.

### 4. Escape Analysis For Latent Times

Do not expand every CDF into a latent integral. That would be correct
but slow.

Add compile-time escape analysis:

A latent source time must be retained when:

- the same source appears as a blocker boundary in another outcome
- the source participates in a shared active tie where readiness order
  matters
- a chained onset depends on that source
- a pool/order statistic needs that child time as a boundary
- the source is referenced by a competitor scenario after merging
  constraints

Otherwise it may be collapsed to CDF/survival.

This gives the desired efficiency:

- simple races stay product formulas
- ordinary
  [`all_of()`](https://niekstevenson.github.io/AccumulatR/reference/all_of.md)
  against independent competitors stays product formulas
- stop-change introduces one latent variable
- stim-selective introduces only the path variables that matter

### 5. Logical Expression Rules

#### `all_of()`

`all_of(a, b, c)` transitions when the last required child becomes
available.

Compiler rule:

- choose one child as active at transition time
- other children are required by transition time
- if any required child’s time escapes, introduce a latent variable
- otherwise use its CDF at transition time

Example:

``` text
all_of(change, stop) at t
```

becomes:

``` text
change active at t, stop before t
stop active at t, change before t
```

If `stop before t` is needed by competitor `inhibit(S, stop)`, expand
it:

``` text
integral f_stop(u) * S_S(u) du
```

#### `first_of()`

`first_of(a, b, c)` transitions when the first child becomes available.

Compiler rule:

- choose one child as active at transition time
- other children must survive beyond transition time
- if the active child shares transition source with another outcome,
  keep readiness time for tie handling

#### `inhibit(reference, blocker)`

`inhibit(reference, blocker)` transitions when `reference` transitions
and the blocker has not transitioned by the reference transition time.

Compiler rule:

``` text
density(reference at x) * survival(blocker past x)
```

If `blocker` is composite, expand only as much as needed. If its
marginal survival is safe, keep it collapsed.

#### `none_of()`

`none_of(x)` remains normalized to `not(x)`.

It should be supported inside conjunctions as an absence condition. Do
not add `none_of` as a separate native semantic form.

### 6. Pool Rules

Pools should be treated as source-like events, but the temporal compiler
must be able to expand them when child times escape.

For `k-of-n`:

- collapsed density/CDF/survival can keep using the existing
  order-statistic code
- expanded terms enumerate which child is the active kth event and which
  children are before/after that time

This mirrors the existing source transition scenario expansion, but the
result should be temporal terms instead of `Before/After` relation
templates.

### 7. Chained Onset Rules

Chained onset must consume exact latent times when available.

The existing oracle already has the idea of `exact_times` for ranked
observations. Extend that mechanism so temporal evaluator frames can
bind:

``` text
source stop has exact latent time u_stop
```

Then dependent leaf channels can evaluate relative to `u_stop`.

If a chained onset source is not bound and cannot be collapsed safely,
the temporal compiler should introduce a latent variable for it.

## Competitor Conditioning

This is the core evaluator change.

Current behavior:

``` text
target transition density at t
* competitor non-win using marginal/relation-overlay survival
```

New behavior:

``` text
sum over target temporal terms:
  target term at t
  * probability(no competitor temporal term wins before/ties ahead)
```

### Inclusion-Exclusion

Keep the existing inclusion-exclusion idea for overlapping competitors,
but apply it to temporal terms.

For each target term:

1.  compile competitor win terms for each competitor outcome by observed
    time
2.  group overlapping competitors as today
3.  enumerate inclusion-exclusion subsets
4.  merge the target term with each competitor-subset term
5.  evaluate the merged temporal integral

This is the clean generalization of the existing competitor block plan.

### Term Merging

Merging two temporal terms means:

- concatenate factors
- concatenate constraints
- unify variables that represent the same exact source time
- reject contradictions
- add tie/readiness constraints when transition sources match

The merge result should be simplified before evaluation.

### Tie And Readiness Semantics

Positive-mass ties must remain in `ExactTransition`.

For shared active transition sources, winner selection depends on
readiness ordering, not just transition time.

Example:

``` r
R1 = all_of(x1, gate)
R2 = all_of(x2, gate)
```

When `gate = t` for both outcomes, the outcome with earlier readiness
wins:

``` text
x1 < x2
```

The temporal IR should represent this directly as a constraint between
readiness variables. This replaces the current special `same_active_*`
correction with a general mechanism.

## Temporal Term Simplifier

Before numeric evaluation, simplify every merged term.

The simplifier is what keeps the approach efficient.

### Required Reductions

1.  Contradiction detection:

``` text
u < v and v < u -> zero
u < 0 -> zero
t < u < 0 -> zero
```

2.  Collapsed independent latent variables:

``` text
integral_0^t f_a(u) du -> F_a(t)
integral_t^Inf f_a(u) du -> S_a(t)
```

3.  Product factor hoisting:

Factors that depend only on observed `t` should be pulled outside
integration.

4.  Duplicate factor elimination where mathematically valid.

5.  Variable ordering:

Sort latent variables topologically from constraints so recursive
integration has tight bounds.

6.  Zero/one channels:

Trigger-failed sources and impossible branches should become zero terms
early.

### Evaluation Classes

After simplification, terms should fall into one of these classes:

``` text
ProductTerm       no latent variables
FiniteIntegral1D  one latent variable
FiniteIntegralND  multiple latent variables
TailIntegral      unbounded latent upper bound
```

Most important failing cases are `ProductTerm` plus `FiniteIntegral1D`.

## Numeric Evaluation

Use the existing quadrature infrastructure where possible.

### Product Terms

Evaluate directly.

### One-Dimensional Integrals

Use existing finite quadrature:

``` cpp
quadrature::build_finite_batch(lower, upper)
```

This handles stop-change and most stim-selective paths.

### Multi-Dimensional Integrals

Use recursive finite integration with variable-specific bounds.

The recursion should:

- evaluate bounds from already-bound parent variables
- bind each latent source exact time into the oracle frame
- cache source channel evaluations per variable/time node where
  practical

Do not add Monte Carlo or approximation here. This is an exact
likelihood engine.

### Infinite Bounds

Use existing tail quadrature only when a variable genuinely has an
infinite upper bound. Many competitor calculations are finite because
top-1 observation conditions all wins before observed `t`.

## Runtime Integration

### Files To Add

Likely new files:

``` text
src/eval/exact_temporal.hpp
src/eval/exact_temporal_compile.hpp
src/eval/exact_temporal_eval.hpp
```

Keep these internal to the exact evaluator.

### Files To Change

Likely touched files:

``` text
src/eval/exact_types.hpp
src/eval/exact_planner.hpp
src/eval/exact_competitor_union.hpp
src/eval/exact_truth.hpp
src/eval/exact_oracle.hpp
src/eval/exact_sequence.hpp
```

The eventual goal is to remove or sharply reduce the role of:

``` text
ExactRelation
RelationView
forced_channels()
relation_template
same_active_* correction code
```

Do not delete them in the first patch unless the replacement is
complete. Start by compiling temporal plans alongside the existing plans
and route only selected exact variants through the new path under tests.

## Implementation Phases

### Phase 1: Lock The Correctness Targets

Already started:

- add `stop_change_shared_trigger` validation
- add `stim_selective_stop` validation

Still needed:

- add comments beside the manual formulas explaining the path
  decomposition
- decide whether `shared_gate_pair NA_mass` tolerance is numerical or
  semantic

Acceptance:

- the new validation cases fail against current engine for the expected
  reason
- existing unrelated validation cases still mostly pass

### Phase 2: Add Temporal IR Without Routing

Add temporal structs and a debug printer.

No likelihood result should change.

Acceptance:

- code compiles
- a unit/debug call can print temporal terms for simple expressions

### Phase 3: Compile Source And Simple Expression Terms

Support:

- leaf/source density/CDF/survival
- `all_of`
- `first_of`
- `inhibit`
- `not` inside conjunction-style absence checks

No pools or chained onset expansion yet, but keep APIs source-general so
pools can be added without redesign.

Acceptance:

- temporal compiler prints expected terms for:
  - simple race
  - `all_of(change, stop)`
  - `inhibit(S, stop)`
  - stim-selective A/B expressions

### Phase 4: Add Term Simplifier And Evaluator

Support:

- product terms
- one-dimensional finite integrals
- observed-time factors
- source exact-time binding in oracle frames

Acceptance:

- hand-run temporal evaluator matches manual stop-change `X` density
- product-only cases match existing exact evaluator

### Phase 5: Temporal Competitor Union

Replace relation-overlay competitor non-win for temporal-routed
variants.

Support:

- inclusion-exclusion over competitor blocks
- term merging
- contradiction detection
- shared active/readiness ordering

Acceptance:

- `stop_change_shared_trigger` passes
- `shared_gate_three_way_tie` and `shared_gate_four_way_tie` still pass
- `guarded_overlapping_competitors` still passes

### Phase 6: Stim-Selective Correctness

Extend term expansion enough for:

- nested `first_of`
- nested `all_of`
- nested `inhibit`
- mapped NA probability computed from corrected visible outcomes

Acceptance:

- `stim_selective_stop` passes
- mapped NA mass agrees with corrected visible response probabilities

### Phase 7: Pools And Chained Onsets

Extend temporal expansion for:

- pool source terms when child times escape
- chained onset exact latent binding

Acceptance:

- existing pooled shared-gate validations still pass
- ranked chained onset tests still pass
- add one explicit temporal pool case where a pooled child time is a
  competitor boundary

### Phase 8: Route Ranked Exact Through Same Terms

Ranked top-k observations should use the same temporal transition terms
with additional exact observed-time bindings.

Acceptance:

- no separate ranked mathematical path for logical outcomes
- ranked independent and ranked chained onset validations still pass

### Phase 9: Remove Obsolete ExactRelation Dependencies

After temporal routing is complete:

- delete old relation-overlay competitor code if unused
- keep only minimal relation helpers if still needed for a reduced fast
  path
- ensure no top-1 exact likelihood path uses marginal competitor
  survival when temporal conditioning is required

Acceptance:

- search confirms exact top-1 likelihood does not route through the old
  relation-overlay shortcut for exact variants that need temporal terms

## Efficiency Policy

The clean solution is not “expand everything into integrals”.

The temporal compiler must aggressively reduce.

### Fast Paths That Should Remain

- simple independent race: direct evaluator
- simple pool versus independent competitor: direct or collapsed exact
  product
- `all_of` with independent competitor: collapsed CDF products
- `inhibit` with independent competitor: collapsed guard
  density/survival

### Cases That Should Use 1D Integration

- stop-change `X`
- shared-gate ties
- stim-selective ignore-active branch
- simple guarded overlapping competitors

### Cases That May Need ND Integration

- multiple escaped latent prerequisites
- nested pools with escaped child times
- chained onset plus guarded competitor

These should be correct first. If benchmarks show a real issue, optimize
the term simplifier and caching before adding special model paths.

## What Not To Do

Do not:

- add a stop-change special case
- add a stim-selective special case
- compute `X` with an R-side helper
- route observed likelihood through a new wrapper that bypasses the
  exact plan
- keep multiplying target marginal density by competitor marginal
  survival for exact variants
- add a third evaluator class
- make direct variants pay temporal compilation or temporal evaluation
  cost

## Acceptance Tests

Minimum validation gate after the fix:

``` text
dev/validation/cases.R
```

must pass for:

- `stop_change_shared_trigger`
- `stim_selective_stop`
- `shared_trigger_conditioning`
- `shared_gate_pair`
- `guarded_positive_mass_tie`
- `shared_gate_three_way_tie`
- `pooled_shared_gate_tie`
- `pooled_guarded_shared_gate_tie`
- `overlapping_composite_competitors`
- `guarded_overlapping_competitors`
- `shared_gate_four_way_tie`

Also run:

``` text
tests/testthat/test_exact_eval_prep.R
tests/testthat/test_loglik_golden.R
tests/testthat/test_response_probabilities.R
```

The golden fixture may need regeneration only after verifying that
changed values are mathematically intended, not because the fixture is
inconvenient.

## Benchmark Gate

Add benchmark rows for:

- simple two-choice direct race
- `stop_change_shared_trigger`
- `stim_selective_stop`
- shared-gate three-way tie
- guarded overlapping competitors

Expected performance:

- direct simple race unchanged
- stop-change and stim-selective slower than current wrong formulas, but
  still practical because they should mostly reduce to 1D integrals
- no large regression on existing exact tie cases

## Final Target

After this fix, `ExactTransition` should answer this question directly:

``` text
What is the probability density that this observed outcome transition occurs
at time t and no competing outcome transition occurs earlier or wins the tie,
given the complete path constraints implied by the target transition?
```

That is the correct abstraction.

Anything less will keep reintroducing the same bug under different model
names.
