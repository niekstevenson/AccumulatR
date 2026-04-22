Yes. The clean way to finish this is to stop layering ad hoc caches and instead **collapse the exact runtime onto two explicit abstractions**:

1. **`ExactOracleScratch`**
   This owns all time-local numeric source-channel state.

2. **`RelationView`**
   This owns all forced-relation state without full-width vector copies.

That is the cleanest integration plan. Anything else is just piling hacks on the current code.

**Current problem**
Right now the exact path is split across too many partially-overlapping mechanisms:

- [exact_oracle.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_oracle.hpp) has a top-level current-time source cache.
- [exact_truth.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp) still mutates `oracle_time_` manually and still merges/copies relation vectors.
- [exact_oracle.hpp:229](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_oracle.hpp:229) onset convolution still bypasses the useful cache layer and recomputes upstream base channels directly.
- [exact_truth.hpp:341](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp:341) `merge_forced_relations(...)` still works in terms of full-width source vectors.

That is why the code is not fully coherent yet.

**The right architecture**
Do not add “a second cache here” and “a small branch there.”

Do this instead.

## 1. One Oracle Scratch Object
Replace the current loose `oracle_time_` + partial caches with one proper scratch object that has **explicit time frames**.

It should own:

- fixed step inputs:
  - compiled variant plan
  - param view
  - first param row
  - trigger state
  - exact sequence state

- mutable per-step caches:
  - `conditional_current[source_id]`
  - `base_current[source_id]`
  - `base_lower_bound[source_id]`

- mutable per-recursion frame:
  - `base_frame_stack[depth][source_id]`

So the oracle exposes a small API like:

- `begin_step(observed_time)`
- `set_conditional_time(t)`
- `push_base_time(t)`
- `pop_base_time()`
- `conditional_source(source_id)`
- `base_source(source_id)`
- `lower_bound_base(source_id)`

That gives one single place where exact source-channel values live.

### Why this is clean
Because it separates three different numeric questions that the current code mixes together:

- “what is the source at the current step time?”
- “what is the unconditional base source at the current recursion time?”
- “what is the base source at the fixed lower bound?”

Those are different caches. They should be explicit, not implicit.

## 2. Use Frame Guards, Not Manual `oracle_time_` Mutation
Right now [exact_truth.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp) does this repeatedly:

- save `current_time`
- assign `oracle_time_ = t`
- evaluate
- restore `oracle_time_`

That is ugly and brittle.

Replace it with RAII frame guards like:

- `ConditionalTimeGuard`
- `BaseTimeGuard`

So:

- truth evaluation enters a conditional-time frame
- onset convolution enters a base-time frame
- leaving the scope restores the previous frame automatically

That is much harder to get wrong, and it removes scattered time-state handling from helpers.

## 3. Cache the Onset-Convolution Layer Properly
This is the important part the current top-level cache did **not** solve.

At [exact_oracle.hpp:229](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_oracle.hpp:229), the onset-dependent leaf integral still does:

- for each inner quadrature node `u`
- call `base_source_channels(source_kind, source_index, u)`
- which recursively recomputes upstream source channels

That is the exact hotspot the fresh profiles still point at.

The clean fix is:

- when evaluating a base source at some inner `u`
- enter a `BaseTimeGuard(u)`
- all upstream `base_source(...)` requests inside that frame hit the frame-local cache

So the reuse is:
- not across trials
- not across arbitrary params
- but **within one recursive node**, which is exactly where the repetition currently is

That is the right kind of scratch reuse.

## 4. Remove Full-Width Forced Relation Vectors
This is the second cleanup, and it should be done cleanly, not incrementally.

Right now:
- [exact_truth.hpp:341](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp:341) copies/fills a `source_count`-wide vector
- [exact_truth.hpp:381](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp:381), [exact_truth.hpp:411](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp:411), [exact_truth.hpp:487](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp:487) keep rebuilding that state

That is coherent with the old dense-source-id move, but not elegant.

The clean fix is to replace this with a `RelationView` abstraction:

- `ExactRelation relation_for(source_id) const`

Back it with:

- a parent relation view
- a scenario-local override template compiled in context

So at context creation time, each scenario should carry:
- touched source ids
- touched relations
- local support ids if needed

At runtime, the evaluator should **view** forced relations, not materialize full vectors.

If you want this to stay fast and not regress small exact models, do not make `relation_for` a linear scan over touched ids. Compile support-local ordinals once and make it O(1) on the relevant support.

So the end state is:
- no `merge_forced_relations(...)`
- no `std::vector<ExactRelation>` copies in the hot exact path
- no evaluator API that depends on a full-width relation vector

## 5. Consolidate Oracle Query Semantics
Right now the source path still has too many layers:

- `source_channels(kind,index,t)`
- `source_channels(key)`
- `compute_source_channels(...)`
- `base_source_channels(...)`
- `base_leaf_channels(...)`

The clean version should reduce this to:

- `conditional_source(source_id)`
- `base_source(source_id)`

And internally:
- conditional source = exact-time / lower-bound adjusted view over base source
- base source = recursive unconditional source at current base frame

That makes the call graph much simpler and makes the caching points obvious.

## 6. Keep Context Static, Keep Scratch Mutable
The line must stay strict:

**Context owns only static descriptors**
- source ordinals
- leaf descriptors
- scenario templates
- support-local ordinals
- max recursion depth if we want preallocated frame stacks

**Scratch owns only mutable numeric state**
- channel caches
- expr caches
- relation overlays if needed
- frame stack state

Do not put dynamic numeric caches into the context.
Do not put planner-like static templates into workspaces.

That split is what keeps the whole thing coherent.

**Clean finish plan**

### Phase 1: Oracle Unification
Implement:
- one `ExactOracleScratch`
- explicit conditional/base/lower-bound caches
- frame guards
- onset recursion routed through base-time frames

Delete/replace:
- ad hoc `current_time_` logic
- manual `oracle_time_` save/restore patterns
- scattered top-level-only cache semantics

### Phase 2: Relation View
Implement:
- `RelationView`
- compiled scenario relation templates
- support-local fast lookup

Delete:
- `merge_forced_relations(...)`
- full-width relation vector copying
- evaluator APIs that accept raw relation vectors

### Phase 3: Consolidation
Refactor truth/oracle helpers so they only talk to:
- `ExactOracleScratch`
- `RelationView`

At that point:
- no exact helper should set time manually
- no exact helper should build relation vectors
- no exact helper should call uncached recursive base source resolution directly

**Finishing requirements**

I would call this finished only if all of the following are true.

### Code completion
- [exact_oracle.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_oracle.hpp) has one unified scratch object with explicit time-frame semantics.
- Onset convolution no longer directly recomputes upstream base channels outside that scratch mechanism.
- [exact_truth.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_truth.hpp) no longer mutates `oracle_time_` manually.
- `merge_forced_relations(...)` is deleted.
- `ForcedExprEvaluator` no longer depends on `std::vector<ExactRelation>*`.

### Cohesion
- exact source evaluation has exactly one source of truth for cached channel access
- exact truth helpers use one relation abstraction
- context stores only static execution templates
- scratch stores only mutable numeric values

### Consolidation
- there is no parallel “old way” and “new way” of source lookup left in exact
- there is no top-level cache path plus separate uncached recursive path left
- there is no full-width relation-vector path left in hot exact evaluation

### Benchmark acceptance
- `example_10_exclusion`, `example_23_ranked_chain`, and `example_2_stop_mixture` must all beat the current state, not just one of them
- no regression on the top 8 reference benchmark
- profile evidence must show:
  - `compute_source_channels` reduced materially
  - `base_leaf_channels` reduced materially
  - `source_channels` no longer almost collapses to “cache miss every time”
  - no hot full-width relation merging remains

That is the clean plan.

The brutal truth is that the current code is close enough that random local tweaks would probably produce some gains, but they would also bloat the framework. If we want elegance and speed together, this needs to be finished as a **single exact-runtime consolidation**, not another round of isolated micro-optimizations.
