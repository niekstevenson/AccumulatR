# Component Fix Plan

## Decision

The current architecture is wrong for the required performance target.

`component` is handled too late, as a runtime mask inside one large evaluator. That keeps dead branches, dead guards, dead competitors, and dead accumulators alive long enough to still cost real time. That is why `component = "go"` is still far slower than a true 2-accumulator model.

The fix is to treat components as compile-time projections, not runtime masks.

## Hard Constraints

1. `go` component trials must run within 5% of the standard 2-accumulator model.
2. The implementation must be clean and direct. No thin wrappers, no duplicated execution stacks, no compatibility scaffolding.
3. There must be no old masked hot path left in place. One execution architecture only.

## Target Design

Split the model into:

- a cold mixture shell:
  - component ids
  - per-trial component weights
  - marginalization logic when component is latent
- hot compiled variants:
  - one projected execution variant per component
  - each variant has its own simplified graph, IR, kernel program, outcomes, competitors, and route caches

The evaluator should never receive `component_idx`.

It should only receive the already-selected compiled variant.

## Core Principle

For a known component:

- evaluate exactly one projected variant

For an unknown component:

- evaluate all projected variants
- combine with the component weights afterward

That is:

- known component: `L(trial) = L_c(trial)`
- latent component: `L(trial) = sum_c w_c(trial) * L_c(trial)`

This handles both fixed and sampled mixtures cleanly. The only difference between fixed and sampled mixtures is the source of `w_c(trial)`.

## Implementation Plan

### 1. Introduce component-conditioned compiled variants

Refactor the native context so the hot path runs against a compiled variant, not a full context plus `component_idx`.

Each variant should own:

- projected outcomes
- projected competitors
- projected IR
- projected kernel program
- route caches derived from that projected structure

Models without components should still produce exactly one variant.

### 2. Project and simplify the model before IR compilation

Do not try to solve this in the runtime evaluator.

Projection should happen before or during IR construction, using the component membership as a compile-time condition.

Required simplifications:

- inactive accumulator event -> impossible event
- `guard(ref, impossible)` -> `ref`
- `guard(impossible, blocker)` -> impossible
- `all_of(..., impossible, ...)` -> impossible
- `all_of(x)` -> `x`
- `first_of(impossible, x, ...)` -> `first_of(x, ...)`
- `first_of(x)` -> `x`
- pools drop inactive members
- empty pool -> impossible
- singleton pool -> member

After simplification:

- rebuild outcomes from the simplified graph
- rebuild competitors from the simplified graph
- then compile IR and kernel from that simplified structure

### 3. Ensure the `go` projection collapses structurally

For the motivating model, the `go` projection must collapse to a true 2-accumulator race.

That means:

- only `A` and `B` remain active
- no stop accumulators remain
- no guard nodes remain
- no STOP / NA branch remains

If the projected `go` variant still contains guards, the design has failed.

### 4. Remove runtime component masking from the hot evaluator

Delete the current runtime component machinery from the likelihood hot path.

That includes:

- `component_active_idx`
- `filter_competitor_ids`
- IR component masks for density/probability routing
- component-dependent outcome allowance checks in the evaluator
- per-component runtime cache entries that only exist to compensate for the unsimplified graph

The evaluator should be unaware that components exist.

Component handling belongs above the evaluator, at variant selection / marginalization.

### 5. Rewrite top-level likelihood dispatch around variants

At the top level:

- if trial `component` is present, choose that one variant
- if trial `component` is absent, evaluate all variants and marginalize

Fixed and sampled mixtures should share the same evaluation path.

Only the component weights differ.

### 6. Stop discarding `component` in sample mode

The current behavior that drops the prepared `component` column for sampled mixtures is wrong for the new design.

If the trial data provide a component, use it.

If the trial data do not provide a component, marginalize.

No special-case fallback behavior.

### 7. Keep a shared parameter layout

Do not create a separate parameter interface per component.

Projected variants should still refer back to the original accumulator indexing, so the existing trial parameter layout stays valid.

This keeps the refactor smaller and avoids wrapper layers.

### 8. Re-run route detection on each projected variant

Do not carry forward route decisions from the unsimplified model.

The exact/shared-gate route detection should run on the projected structure.

That is one of the reasons this design is better than runtime masking: simplification can make existing fast paths become available naturally.

### 9. Make the new architecture the only architecture

There should be one execution model:

- compile projected variants
- select one variant or marginalize across variants
- evaluate

Do not keep the old masked evaluator around as a fallback branch.

That would be extra code, weaker guarantees, and permanent maintenance debt.

## Acceptance Gates

### Performance gate

For the motivating model:

- `component = "go"` likelihood must be within 5% of the true 2-accumulator model under the same benchmark workload

### Structural gate

For the `go` projected variant:

- no guard nodes
- no stop accumulators
- no ignore accumulators
- no STOP / NA path

### Profile gate

For `component = "go"` likelihood:

- no `kernel_guard_eval_callback` on the hot stack
- no `guard_cdf_internal` on the hot stack

### Cleanup gate

- no runtime component masking logic remains in the hot evaluator
- no compatibility branch to the old masked path

## What Not To Do

Do not:

- add more guard short-circuits and pretend the architecture is fixed
- keep both the old masked path and the new projected path
- preserve the current behavior that ignores explicit `component` in sampled mixtures

Those options are smaller only superficially. They keep the wrong design alive and will miss the target or make the codebase worse.

## First Step

Build the projection pass.

Concretely:

- add a component-conditioned graph simplifier before IR compilation
- make it produce a projected graph for one component
- verify on the motivating model that the `go` projection contains only `A`, `B`, and direct competition

If that step does not produce the right projected structure, nothing downstream matters.
