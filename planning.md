# Exact Evaluation Planning Contract

## Blunt Diagnosis

The exact evaluator is still too much of an interpreter.

The recent fixes moved correctness in the right direction, but they are
not the final architecture. They made the old evaluator smarter about
relevance, temporal bounds, competitor conditioning, and some expression
kernels. That was useful. It was also a bridge.

The framework still spends too much likelihood time asking structural
questions that the model context already knows the answer to:

- What kind of expression is this?
- Which sources are relevant under this condition?
- What does this guard mean under this
  exact-time/source-bound/source-order state?
- Which normalizer applies here?
- Which competitors are possible?
- Which cache entries are invalidated by this temporary conditioning
  overlay?
- Which recursive evaluator path should be entered next?

That is planning work. It does not belong inside likelihood evaluation.

The target architecture is a compiled conditional evaluation DAG: a
finite, typed, condition-keyed math schedule produced once during
context construction. Evaluation should load trial data, fill
parameterized source channels, and execute integer-coded math nodes. It
should not rediscover the model.

## Core Claim

All current benchmark and validation models are reducible to precompiled
paths.

That is not because the models are simple. It is because their shapes
are fixed before likelihood evaluation:

- The semantic graph is fixed.
- The expression tree is fixed.
- The number of sources is fixed.
- The transition graph is fixed.
- The ranked chain structure is fixed.
- The mixture/component structure is fixed.
- The possible conditioning facts are finite.
- The needed PDFs, CDFs, survivals, normalizers, and competitor terms
  are finite.

Only numeric values vary at evaluation time:

- parameters
- observed times
- observed outcomes
- latent integration coordinates
- state slot values for ranked/exact transitions

So the right abstraction is recursive compiled kernels. Every nested
operation requests typed math from the same compiler. The compiler
canonicalizes the request, deduplicates equivalent requests, assigns
cache slots, and emits a topological schedule.

This is the key shift:

> Runtime recursion through semantic objects must become compile-time
> recursion through typed math requests.

## The Right Mental Model

The final evaluator should look like a network of typed math nodes:

``` text
source_pdf(source, condition, t)
source_cdf(source, condition, t)
source_survival(source, condition, t)

expr_pdf(expr, condition, t)
expr_cdf(expr, condition, t)
expr_survival(expr, condition, t)

scenario_readiness_pdf(scenario, condition, t)
scenario_truth_cdf(scenario, condition, t)
competitor_non_win(target, competitor_group, condition, t)

transition_probability(transition, state, observation)
ranked_step_probability(step, state, observation)
```

Internally those nodes can expand into products, sums, complements,
convolutions, inclusion-exclusion terms, dynamic-programming tables, or
quadrature kernels. But that expansion is a compiler responsibility.

Likelihood evaluation should see already compiled nodes with integer
ids, assigned input slots, assigned cache slots, and fixed child spans.
It should not know whether a node came from an `all_of`, `first_of`,
guard, stop-change model, stimulus-selective model, ranked chain, or
mixture branch except through the node opcodes it executes.

## Compiled Request Identity

Every request for math must have a canonical identity.

``` text
CompiledRequest {
  value_kind      // pdf, cdf, survival, mass, log_mass, expectation
  subject_kind    // source, expression, scenario, transition, competitor, ranked_step
  subject_id
  condition_key
  time_var
  state_key
}
```

The exact fields can be implemented differently, but this identity must
exist conceptually.

Two requests with the same identity must map to the same node. If two
independently discovered paths need
`expr_cdf(expr_7, condition_12, t_readiness)`, they must share one
compiled node and one cache slot. They must not build two evaluator
subtrees.

## Node Kinds

The node vocabulary should be small and typed. The goal is not one
special node per model feature. The goal is a compact algebra that all
model features lower into.

Primitive distribution nodes:

- `SourcePDF`
- `SourceCDF`
- `SourceSurvival`
- `SourceMass`

Expression nodes:

- `ExprPDF`
- `ExprCDF`
- `ExprSurvival`
- `ExprMass`

Algebra nodes:

- `Product`
- `Sum`
- `Difference`
- `Complement`
- `Ratio`
- `Log`
- `Exp`
- `Min`
- `Max`

Integration nodes:

- `Integral`
- `Convolution`
- `BoundedIntegral`
- `QuadratureKernel`

Condition nodes:

- `ConditionNormalizer`
- `ConditionedPDF`
- `ConditionedCDF`
- `ConditionedSurvival`

Scenario and transition nodes:

- `ScenarioReadinessPDF`
- `ScenarioReadinessCDF`
- `ScenarioTruthCDF`
- `ScenarioTailSurvival`
- `TransitionProbability`
- `RankedStepProbability`

Competitor nodes:

- `CompetitorSubsetWinCDF`
- `CompetitorSubsetWinMass`
- `CompetitorNonWin`
- `CompetitorInclusionTerm`

State nodes:

- `StateRead`
- `StateWrite`
- `StateCondition`
- `MixtureSelect`
- `TriggerStateSelect`

This list can shrink during implementation. It should not grow into
feature-specific fast paths. If a proposed node only exists to patch one
model, it is suspect.

## Recursive Compilation

Every semantic operation should compile through a common request
interface.

Expression compilation examples:

``` text
compile(expr, PDF, condition, t)
compile(expr, CDF, condition, t)
compile(expr, SURVIVAL, condition, t)
```

Transition compilation examples:

``` text
compile(scenario, TRUTH_CDF, condition, t)
compile(scenario, READINESS_PDF, condition, t)
compile(target_scenario, COMPETITOR_NON_WIN, condition, t)
compile(transition, PROBABILITY, state_key, observation_key)
```

Ranked-chain compilation examples:

``` text
compile(ranked_step, PROBABILITY, ranked_state_key, observation_key)
compile(next_ranked_state, STATE_WRITE, ranked_state_key, observation_key)
```

The compiler is allowed to recurse. Evaluation is not.

The same model part should compile the same way whether it appears at
top level, inside a guard, inside a competitor term, inside a scenario
readiness calculation, or inside a ranked chain state. That is the main
reason this approach is worth doing.

## Expression Lowering Rules

The exact formulas belong in implementation, but the lowering
responsibilities are clear.

`event/source`:

- PDF/CDF/survival requests lower to primitive source channel nodes.
- Conditions specialize the source channel by condition key.

`all_of`:

- CDF lowers to product of child CDFs under the same condition.
- PDF lowers to sum over active child density times other child CDFs.
- Survival lowers to complement of CDF unless a more stable equivalent
  is compiled.

`first_of`:

- Survival lowers to product of child survivals.
- CDF lowers to complement of survival.
- PDF lowers to sum over active child density times other child
  survivals.

`inhibit(ref, blocker)`:

- Density lowers to
  `ref_pdf(condition, t) * blocker_survival(specialized_condition, t)`.
- CDF lowers to an integral of that density unless a closed
  form/source-specialized channel exists.
- Exact-time, source-order, and upper-bound facts specialize the
  condition key. They do not create runtime branches.

`none_of`/negation:

- Lowers to complement forms with explicit value kind and condition key.

Pools and k-of-n:

- Lower to source-like compiled channels backed by a finite
  dynamic-programming schedule.
- The DP schedule is compiled once. Runtime fills numeric source values
  and executes the table.

Onset/after relations:

- Lower to source-like or expression-like channels with explicit
  latent-time slots.
- Convolutions, shifts, and lower-bound constraints are compiled
  kernels, not runtime expression interpretation.

## Condition System

The condition system is the core of the refactor.

Conditions must be immutable, canonical, and finite. Runtime must select
a condition key, not mutate an oracle overlay.

A condition key should represent facts such as:

- source has exact time `t_i`
- source is known before time `u`
- source is known after time `u`
- expression is known before time `u`
- guard blocker has not occurred before time `u`
- source order relation holds
- ranked-chain lower bound is `u`
- outcome/source has already been consumed by a ranked step
- trigger state is active

The key point is that the fact shape is compiled. The numeric values are
slots.

Example:

``` text
condition_17:
  source_3 exact_at slot ranked_t_0
  source_8 after slot readiness_u
  expr_4 before slot observation_t
```

During evaluation, `ranked_t_0`, `readiness_u`, and `observation_t` get
numeric values. The evaluator does not decide what those facts mean. It
executes the nodes compiled for `condition_17`.

## Normalizers

Normalizers must be explicit nodes.

Runtime must not discover that a normalizer is needed by inspecting a
condition. Runtime must not ask an oracle to recompute a conditional
distribution after temporary overlay mutation.

If a conditioned expression needs:

``` text
P(expr <= upper_bound | condition)
```

then the compiler must emit:

``` text
ConditionNormalizer(expr, condition_key, upper_bound_slot)
```

All dependent conditioned PDFs/CDFs/survivals reference that node by id.

This is non-negotiable. Hidden normalizer discovery is exactly the kind
of repeated planning that makes the current evaluator hard to reason
about and expensive in models like stimulus-selective.

## Time Variables

Time variables must be named compile-time slots.

The evaluator can assign numeric values to slots. It must not invent new
time variables.

Expected slots include:

- observed response time
- scenario readiness integration time
- scenario truth integration time
- guard integration time
- ranked exact time for each observed step
- ranked lower-bound time
- onset latent time
- source-bound time
- component/trigger-specific time if required

The compiler owns the time-variable table. Nodes refer to integer time
ids.

## Scenario and Transition Compilation

A transition probability should lower to a compiled scenario schedule.

For each scenario, the compiler should produce the needed pieces:

- active outcome/source density or mass
- readiness PDF/CDF
- scenario truth CDF
- tail survival
- competitor non-win probability
- normalizers required by conditions
- observation mapping terms

The transition node should combine these pieces with fixed child ids.

Evaluation should not loop over semantic scenarios and decide how each
scenario works. It should loop over compiled scenario node ids.

## Competitor Compilation

Competitor conditioning is one of the current weak points.

The compiler must precompute:

- competitor groups
- relevant competitor scenarios
- same-active-source corrections
- impossible competitor masks
- inclusion-exclusion subsets
- signs
- shared-condition keys
- readiness-time slots
- required competitor CDF/mass nodes

Runtime must not filter competitors based on expression relevance.
Runtime must not build subsets. Runtime must not decide whether an `Or`
is absorbed by a condition. Those are compile-time facts.

The runtime can skip a compiled competitor term only if the compiled
node itself has a numeric impossible flag or zero contribution. That is
numeric pruning, not semantic planning.

## Ranked Chains

Ranked chains are not a separate evaluator architecture.

They are repeated applications of the same transition compiler under a
state key.

A ranked state key should encode finite structural facts:

- which outcomes/sources have already been consumed
- lower-bound slot for the next event
- exact-time slots for previous ranked events
- active trigger/component state
- source-order facts implied by previous events

Each ranked step should select a compiled transition schedule for the
current state key and observation key. State updates should be compiled
writes into the next state key.

Runtime should not reinterpret the ranked chain. It should walk the
precompiled state transition path for the observed ranking.

## Mixtures, Components, and Triggers

Mixtures and trigger states should be outer compiled schedules, not
runtime semantic branches.

The compiler should emit:

- component weight nodes
- component-specific source channel nodes
- trigger-state condition keys
- per-component transition schedules
- mixture sum nodes

Runtime can loop over components because that is math. Runtime cannot
decide the component semantics.

## Reducibility of Current Models

The benchmark and validation families are reducible to this
architecture.

Simple independent races:

- Source PDFs/CDFs/survivals and transition products.
- Almost no overhead beyond primitive channels.

Pools and k-of-n models:

- Finite DP tables over fixed source sets.
- Compiled once, executed as array math.

`all_of` and `first_of` models:

- Recursive expression PDF/CDF/survival lowering.
- Products and sums with fixed child spans.

Guards and nested guards:

- Condition-specialized blocker survival and normalizer nodes.
- Guard semantics are compiled by condition key.

Shared gates and shared sources:

- Shared node identities deduplicate the repeated CDF/PDF/survival
  requests.
- Correctness depends on condition-key equality, not runtime cache
  invalidation.

NA/mapped outcomes:

- Observation mapping compiles to fixed scenario/output masks.
- Runtime selects the observed mapping id.

Mixtures/components/triggers/q:

- Component schedules and trigger-state condition keys are finite.
- Runtime sums weighted compiled paths.

Chained onset:

- Onset constraints become time-slot conditions and convolution/shift
  kernels.
- No need for special runtime interpretation.

Ranked chains:

- Finite sequence of transition schedules over ranked state keys.
- Exact previous times are numeric slots in compiled conditions.

`stop_change`:

- Requires exact transition conditions and competitor non-win terms.
- Fits the scenario/competitor/condition-key model directly.

`stim_selective`:

- Requires repeated conditioned guard/scenario/competitor terms.
- It is slow today because these terms are interpreted and repeatedly
  conditioned at evaluation time.
- It should become a dense schedule of deduplicated condition-keyed
  kernels.

If a validation model fails under this architecture, that is evidence of
a missing lowering rule or wrong math. It is not evidence that the model
is irreducible.

## Why Stimulus-Selective Is Slow

Stimulus-selective models are expensive because they combine the worst
current runtime behaviors:

- repeated conditioned expression evaluation
- guard specialization under different source/order facts
- competitor non-win calculations
- recursive normalizer discovery
- overlay mutation and broad invalidation
- repeated equivalent CDF/PDF requests under conditions that should be
  canonical

Some of the work is real math:

- evaluating source distributions
- multiplying PDF/CDF/survival terms
- integrating guarded densities
- summing competitor inclusion terms
- summing scenario probabilities

That cost will remain.

The reducible overhead is everything else:

- walking expression trees repeatedly
- allocating temporary condition frames
- mutating and restoring oracle overlays
- recomputing which facts matter
- invalidating caches broadly
- rebuilding competitor masks/subsets
- rediscovering normalizers
- selecting recursive evaluator paths dynamically

The target refactor should remove the second category from likelihood
evaluation.

## Evaluation Is Allowed To Do

Likelihood evaluation may:

- load parameter values
- load trial observations
- fill primitive source distribution channels
- write numeric values into preallocated time/state slots
- execute a topological node schedule
- execute compiled quadrature kernels
- read and write assigned cache/workspace slots
- loop over precompiled child spans
- skip nodes using compiled numeric impossible/zero flags
- sum compiled scenario/component/ranked-step terms

That is evaluation.

## Evaluation Must Not Do

Likelihood evaluation must not:

- switch on semantic expression kind
- recurse through semantic expression objects
- allocate condition frames
- mutate oracle overlays
- restore oracle overlays
- broadly invalidate caches because an overlay changed
- discover whether a condition is relevant
- decide which normalizer is needed
- build competitor groups
- build competitor subsets
- decide competitor relevance
- decide whether an `Or` has been absorbed by a condition
- specialize guard semantics
- create source-order facts
- create exact-time facts
- traverse transition semantic graphs
- allocate temporary child vectors for model structure
- use runtime maps keyed by semantic objects in inner loops

If any of these happen in the final exact likelihood path, the refactor
is incomplete.

## Data Structures

The implementation should move toward compact arrays.

Node table:

``` text
node_kind[]
value_kind[]
subject_id[]
condition_id[]
time_id[]
state_id[]
child_begin[]
child_count[]
op_begin[]
op_count[]
cache_slot[]
workspace_slot[]
```

Condition table:

``` text
condition_fact_begin[]
condition_fact_count[]
fact_kind[]
fact_subject_id[]
fact_time_slot[]
fact_relation[]
normalizer_node_id[]
```

Schedule:

``` text
topological_node_ids[]
integral_kernel_ids[]
scenario_node_spans[]
transition_node_spans[]
ranked_state_paths[]
```

Workspace:

``` text
source_channel_values[]
time_slot_values[]
state_slot_values[]
node_cache_values[]
node_cache_epoch_or_valid_bits[]
quadrature_workspace[]
```

The exact layout can change. The direction should not. Inner evaluation
wants arrays and integer ids, not semantic object traversal.

## Implementation Plan

### Phase 1: Freeze The Semantic Contract

Define the complete set of value kinds, subject kinds, condition fact
kinds, and state key kinds.

Success criteria:

- Every current exact expression operation has a compile mapping.
- Every current transition/scenario operation has a compile mapping.
- Ranked-chain and mixture/component state are represented as state
  keys.
- There are no runtime-only semantic cases left undocumented.

### Phase 2: Build Canonical Node And Condition Tables

Implement canonical request interning.

Success criteria:

- Equivalent requests dedupe to one node id.
- Equivalent conditions dedupe to one condition id.
- Normalizers are explicit nodes.
- Cache slots are assigned during planning.
- No likelihood-time condition frame allocation is needed for compiled
  nodes.

### Phase 3: Compile Expression Kernels

Lower expression PDF/CDF/survival/mass requests into typed nodes
recursively.

Success criteria:

- The forced/exact expression evaluator no longer switches over
  expression kind during likelihood evaluation.
- Expression children are fixed spans in the node table.
- `all_of`, `first_of`, guards, negation, pools, and onset expressions
  execute through the same schedule machinery.

### Phase 4: Compile Source Channels

Move source conditioning into compiled source-channel requests.

Success criteria:

- Source exact-time, lower-bound, upper-bound, and order facts are
  condition-keyed slots.
- There is no oracle overlay mutation in the source evaluation path.
- Primitive distribution evaluation is called through precompiled
  channel ids.

### Phase 5: Compile Scenario Truth And Readiness

Lower scenario readiness, truth, and tail terms into nodes.

Success criteria:

- Scenario truth CDF/density is a schedule, not recursive interpreter
  logic.
- Scenario tail survival uses compiled expression/source nodes.
- Scenario node spans are known before likelihood evaluation.

### Phase 6: Compile Competitor Non-Win

Lower competitor conditioning into inclusion/exclusion schedules.

Success criteria:

- Competitor groups are precomputed.
- Competitor subsets and signs are precomputed.
- Same-active-source corrections are precomputed.
- Runtime does not filter competitor relevance.
- Runtime does not build competitor masks.

### Phase 7: Compile Ranked, Mixture, And Trigger Paths

Represent ranked-chain progression and mixture/trigger selection as
compiled state schedules.

Success criteria:

- Ranked steps use compiled transition schedules keyed by ranked state.
- State updates are compiled writes.
- Mixture/component loops are fixed math loops over compiled children.
- Trigger states select condition keys rather than semantic branches.

### Phase 8: Delete The Bridge Architecture

Remove the old interpreter behavior from the exact likelihood path.

Success criteria:

- No overlay mutation.
- No broad cache invalidation.
- No dynamic condition frame allocation in inner evaluation.
- No recursive semantic graph traversal in likelihood evaluation.
- Old bridge classes either disappear or become planning-only helpers.

### Phase 9: Verification And Performance Gate

Use validation and benchmark scripts as gates.

Success criteria:

- All validation models pass, including `stop_change` and
  `stim_selective`.
- Manual derivations remain the correctness target, not engine
  self-consistency.
- Benchmark-speed output records current time, commit, and all model
  timings.
- Simple models do not pay for complex model machinery.
- `stim_selective` is materially faster than the bridge implementation
  or, at minimum, no slower while removing forbidden runtime planning.
- Remaining hot spots correspond to actual math: source distributions,
  products/sums, quadrature, inclusion-exclusion, and scenario
  summation.

## Performance Standard

The standard is not merely “faster”.

The standard is:

> Evaluation cost should be mostly math cost.

For simple models, complexity should remain near the primitive
distribution math. They should not pay for competitor infrastructure,
ranked-state machinery, or guard condition machinery they do not use.

For complex models, cost should scale with the compiled schedule:

- number of source channel evaluations
- number of algebra nodes
- number of quadrature evaluations
- number of competitor inclusion terms
- number of scenario terms
- number of ranked steps

It should not scale with repeated planning work.

Cache behavior should be explainable from the node table. If a model is
slow, profiling should point to named node kinds and scheduled kernels,
not to hash maps, overlay restore logic, semantic tree walks, or
condition discovery.

## Correctness Standard

The engine must match the manually derived model, not its previous
output.

Correctness requires:

- exact transition probabilities are derived from the explicit model
  semantics
- source/order/exact-time conditions are applied mathematically
- guards use the correct conditioned blocker survival/CDF/PDF
- competitor non-win terms condition on the same facts as the target
  scenario
- normalizers are explicit and shared
- ranked-chain state carries exact previous event times and used
  outcomes
- mixtures/components/triggers preserve their intended conditioning

Validation tests should compare against handwritten or independently
derived targets where possible. Engine-vs-engine agreement is not
enough.

## Non-Goals

This refactor should not become:

- a new wrapper around the old oracle
- model-specific fast paths for `stop_change` or `stim_selective`
- a pile of caches hiding interpreter behavior
- another dynamic condition-frame layer
- a broad invalidation scheme with better names
- a runtime planner embedded in the evaluator

If the implementation cannot delete old runtime planning paths, it has
not reached the target.

## Practical Constraint

This is a large refactor. It can be staged. It does not have to land in
one edit.

But each stage must move code from runtime evaluation to context
planning. A stage that only adds another cache, wrapper, or special-case
correction without deleting runtime planning responsibility should be
rejected.

The honest stopping point after each stage should state:

- which planning responsibility was moved out of evaluation
- which forbidden runtime operations remain
- which validation models changed
- which benchmark costs moved from overhead to math

That is how this work stays real.
