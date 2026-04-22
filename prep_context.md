Yes. The clean general fix is to push more into **context preparation**, not by adding special paths, but by finishing the exact runtime as a **compiled execution schedule**.

The goal is simple:

- the exact path to walk should already be determined
- the support-local ids to touch should already be determined
- the exact operations to perform should already be determined
- which features are actually present should already be determined

At runtime, we should then mostly be paying for:

- source channel math
- expression math
- quadrature
- simple array reads and writes

Not for:

- semantic interpretation
- repeated tree traversal decisions
- generic feature plumbing that the model does not use
- wrapper layers deciding what to do next

That is the standard.

**Current problem**
Right now the exact runtime is still too interpreted.

Even after source ids and some templates were lowered, the hot path still does too much dynamic work:

- it switches on `expr_kind` at runtime
- it recursively walks expression children at runtime
- it resolves scenario structure at runtime
- it asks which sources and relations matter at runtime
- it decides which oracle path to use at runtime
- it carries feature complexity even for models that do not use those features

That is why simple exact models like `example_6_dual_path` still pay too much overhead.

The runtime is cleaner than it was, but it is still not honest enough.

**What should go into context**
Context should own the full static execution plan.

Not values.
Not trial-specific state.
Not parameter-specific numerics.

But yes, absolutely:

- the exact path to walk
- the support-local ids to touch
- the exact operations to perform
- which features are actually present

That means the context should compile and store all of the following.

## 1. Support-Local Source Space
Each exact step should operate on a compact local source space, not the full semantic/source-key space.

Context should compile:

- `local_source_ids`
- `local_leaf_source_ids`
- `local_pool_source_ids`
- `global_to_local_source_id`
- `local_to_global_source_id`

The hot path should not ask:

- “what kind of source is this?”
- “what is its ordinal?”
- “which source key does this expr refer to?”

Those questions should be answered once during context creation.

## 2. Compiled Source Evaluation Order
Context should compile a topological source-evaluation schedule for each exact step family.

That schedule should encode, per local source id:

- whether the source is a leaf or pool
- whether the source has absolute onset
- whether the source has recursive onset
- which onset source id it depends on
- which pool member ids it depends on
- whether lower-bound conditioning can matter
- whether exact-time collapse can matter

The runtime should then execute sources in compiled order, not recursively discover dependencies on the fly.

## 3. Compiled Expression Programs
Expression evaluation should not be driven by recursive semantic walking in the hot path.

Context should compile:

- `expr_cdf_program`
- `expr_density_program`
- topological evaluation order
- child spans already flattened to local program indices
- feature flags such as:
  - `has_guard`
  - `has_not`
  - `has_or`
  - `has_relation_overlay`

That means runtime should not keep switching on semantic expression kind in recursive helper calls if the same information can be encoded once in a compiled program.

## 4. Compiled Scenario Step Programs
For each target outcome, context should compile a step program that already says:

- active source id
- before source ids
- after source ids
- ready expression program ids
- tail expression program ids
- competitor block ids
- which local sources are touched by this scenario
- which features this scenario actually uses

That should replace generic “iterate scenario objects and discover what applies” behavior in the hot path.

## 5. Compiled Competitor Programs
Competitor evaluation should also be prepared structurally in context.

For each target outcome, compile:

- competitor blocks
- subset inclusion structure
- singleton/simplified cases
- touched local sources
- touched expression ids

The runtime should evaluate compiled competitor formulas, not interpret general set structure repeatedly.

## 6. Compiled Relation Access Shape
If relation overlays exist, context should compile:

- touched local source ids
- local relation tables or fast support-local access
- feature mask indicating whether relation logic is needed at all

The runtime should not carry relation machinery for models that do not use it.

That is the key opt-in rule:

- complexity should be compiled in only when the model needs it
- not carried generically and skipped by branches

## 7. Feature Masks
Each exact variant, and ideally each exact step program, should carry explicit feature masks.

At minimum:

- `has_triggers`
- `has_exact_times`
- `has_recursive_onsets`
- `has_pools`
- `has_guards`
- `has_relation_overlays`
- `has_ready_exprs`
- `has_tail_exprs`
- `has_competitor_blocks`

This is how complexity becomes opt-in rather than opt-out.

For a simple model like `example_6_dual_path`, the compiled program should simply not include machinery for:

- guards
- relation overlays
- exact-time conditioning
- generic scenario readiness if it is not structurally needed

Not “include it and skip it.”
Just do not include it.

**What should stay out of context**
This line must stay strict.

Context must not store:

- source channel values
- expression numeric values
- readiness values
- current quadrature-node state
- current parameter-row values
- current lower bound
- current exact-time state
- trial-specific caches

Those belong in workspace.

The split must be:

- context = static compiled execution plan
- workspace = mutable numeric arrays for one evaluation

If that line gets blurred, the design will bloat and become incoherent again.

**What runtime should do after this**
The ideal exact runtime shape is:

1. Reset workspace for one exact step.
2. Set current frame times.
3. Run compiled source program over local source ids.
4. Run compiled expression program over local expression ids.
5. Run compiled scenario/competitor accumulators.
6. Return probability.

No semantic interpretation in the hot loop.
No recursive tree walking in the hot loop unless the compiled program explicitly requires a quadrature subframe.

That is the real target.

**What this should replace**
This work is not done if the runtime still fundamentally relies on:

- recursive semantic expression evaluation as the default exact path
- generic scenario interpretation as the default exact path
- dynamic source-key/source-kind resolution in the hot path
- generic feature plumbing surviving in simple models that do not use those features

Those are exactly the things context preparation is supposed to remove.

**Success criteria**
This is only a success if all of the following are true.

## Code completion
- Exact context preparation emits support-local source ids for every exact step family.
- Exact context preparation emits compiled source evaluation order.
- Exact context preparation emits compiled expression programs for `cdf` and `density`.
- Exact context preparation emits compiled scenario/competitor step programs.
- Exact context preparation emits feature masks that let runtime avoid unused machinery entirely.

## Runtime behavior
- The hot exact runtime no longer asks semantic questions such as source kind, source ordinal, or expression kind in the main evaluation loop.
- The hot exact runtime no longer recursively walks generic expression trees as its default evaluation mechanism.
- The hot exact runtime executes compiled arrays/programs plus numeric workspace state.
- Simple exact models only instantiate the machinery implied by their feature mask.

## Cohesion
- There is one clear contract:
  - context owns static structure
  - workspace owns mutable numerics
- There is not a parallel “compiled path” and “old interpreted path” both still live for ordinary exact evaluation.
- Source, expression, scenario, and competitor execution all consume the same compiled support-local id space.

## Consolidation
- There is no semantic/source-key decoding left in the hot exact path.
- There is no generic feature machinery being carried through simple exact models unless the compiled feature mask says it is required.
- There is no repeated dynamic discovery of traversal order in the hot exact path.

## Benchmark acceptance
- `example_6_dual_path` must beat the current state materially.
- `example_10_exclusion`, `example_23_ranked_chain`, and `example_2_stop_mixture` must not regress.
- The top benchmark suite must not regress broadly.

## Profile acceptance
Profiles must show that:

- exact runtime time shifts away from generic evaluator/oracle wrapper layers
- exact runtime time shifts toward actual source math, expression math, and quadrature
- simple exact models no longer spend obvious time in generic feature plumbing they do not use

If the benchmark improves only by one special case, this is not done.
If the code gets cleaner but the generic interpreted path is still the real hot path, this is not done.
If support-local compiled schedules exist but runtime still mostly ignores them, this is not done.

**Failure conditions**
This should be called a failure, not a success, if any of these are true:

- the exact hot path still interprets semantic trees by default
- the exact hot path still resolves source structure dynamically by default
- feature complexity is still opt-out instead of opt-in
- the change adds yet another layer of wrappers without deleting the old interpreted behavior
- the codebase ends up with both compiled and generic paths live for the same ordinary exact work
- the benchmark win depends on one narrow model while the general exact path remains structurally unchanged

That is the brutal standard.

The point of this work is not to add another optimization layer.
The point is to make exact evaluation structurally honest:

- context decides the path
- workspace holds the values
- runtime executes the program

Anything weaker than that is a partial cleanup, not the real fix.
