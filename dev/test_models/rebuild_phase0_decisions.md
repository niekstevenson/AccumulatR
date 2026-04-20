# Phase 0 Decisions

## Purpose

This document closes the main architecture questions before implementation.

These decisions are intentionally narrow.
They should only be reopened for one of two reasons:

1. a correctness counterexample
2. benchmark evidence that a decision is materially wrong

Not for taste.

## Decision 1: Keep Two Evaluator Classes

Decision:

- keep exactly two evaluator classes:
  - `Direct`
  - `ExactTransition`

Reason:

- one exact-only evaluator is simpler on paper, but it violates the explicit design goal that simple models should not pay for supported complexity
- one "single evaluator" with an internal direct fast mode is just two evaluator classes hidden behind one API
- making the split explicit is less bloated than pretending the split does not exist and then reintroducing it as branches

What this does not mean:

- it does not mean two semantic systems
- it does not mean duplicated model compilation
- it does not mean model-specific fast paths

It means:

- one semantic core
- one compiler pipeline
- two execution backends selected by compiled structure

Rejection criterion:

- do not add a third evaluator class
- do not create special-case evaluators for ties, triggers, or guards

## Decision 2: Normalize `none_of()` Early

Decision:

- normalize `none_of(x)` directly to `not(x)` during DSL normalization

Reason:

- that is already the current R-side meaning
- there is no value in keeping a separate semantic node for the same thing

Consequence:

- the semantic graph has `Not`
- it does not have a separate `NoneOf`

## Decision 3: Guard Semantics Are Realized-Time Semantics

Decision:

- `guard(reference, blocker)` is defined at the realized transition time of `reference`

Operationally:

- if `reference` realizes at time `t`, the guard succeeds only if the blocker is absent by `t`, after applying any `unless` logic

Reason:

- this prevents the old failure mode where observed outcome density was multiplied by marginal competitor survival instead of post-transition conditional survival

Consequence:

- any model whose correctness depends on realized transition conditioning must be classified to `ExactTransition`

## Decision 4: Ties Belong Only In `ExactTransition`

Decision:

- positive-mass tie handling lives only in the exact-transition evaluator

Reason:

- ties are not a leaf tweak
- they arise from realized shared structure
- putting tie handling anywhere else will produce new route-specific exceptions

Consequence:

- if the compiler detects positive-mass ties, the variant is exact

## Decision 5: Shared Trigger Is Not A Leaf Property

Decision:

- independent trigger failure is leaf-local
- shared trigger failure is a latent discrete state above leaves

Reason:

- independent trigger changes only one leaf's start probability
- shared trigger couples several leaves and therefore changes conditional semantics

Consequence:

- no leaf helper should try to emulate shared-trigger semantics by itself

## Decision 6: No Dedicated Shared-Trigger Direct Evaluator At The Start

Decision:

- do not introduce a third evaluator or a dedicated shared-trigger direct backend in Phase 1-5

Reason:

- that would be premature complexity
- if a trivial shared-trigger model deserves a closed direct treatment, it can still be compiled into `Direct` later without adding a new evaluator class

Consequence:

- start with:
  - `Direct` for structurally direct variants
  - `ExactTransition` for anything that needs exact shared-trigger conditioning

## Decision 7: Ranked Logical Outcomes Do Not Block The Rebuild

Decision:

- implement ranked direct outcomes first
- implement ranked logical outcomes only after top-1 exact semantics are stable

Reason:

- ranked logical outcomes are real, but they are not the right place to begin
- trying to solve them before top-1 exact semantics are stable will likely recreate the old architecture failure

Consequence:

- early ranked support is intentionally narrower than the final target

## Decision 8: Simulation And Likelihood Share Projected Variants

Decision:

- simulation and likelihood share the same projected semantic variants
- they do not need to share the same low-level plan struct

Reason:

- projection and simplification are semantic facts
- low-level execution details are backend concerns

Consequence:

- there is one compile-time view of each component-conditioned model

## Decision 9: Observation Wrappers Stay Above Core Semantics

Decision:

- `map_outcome_to`
- `guess`
- truncation of missing later ranks

all stay above core event semantics

Reason:

- these are observation-layer transformations
- they should not influence the meaning of primitive events or outcome expressions

Consequence:

- evaluators operate on semantic outcome labels first
- wrappers transform evaluator results into observed data semantics afterward

## Decision 10: Runtime Must Be Lowered And Batched

Decision:

- evaluators must not consume semantic trees directly
- evaluators must consume lowered runtime programs with integer-indexed dense layouts
- trial evaluation should batch over homogeneous compiled paths whenever practical

Reason:

- general speed comes from regular kernels and contiguous data, not from piling on model-specific shortcuts
- if the hot path still does string lookup, semantic dispatch, or component filtering per trial, simple models will pay for supported complexity again
- this is the prerequisite for later SIMD/vectorization work without redesign

What this does not mean:

- it does not mean adding many evaluator classes
- it does not mean inventing a special fast path for each model family
- it does not mean premature micro-optimization before semantics are correct

It means:

- one compile pipeline lowers semantics into stable runtime programs
- one batching layer groups trials by compiled path
- one direct kernel family and one exact kernel family consume those batches

Rejection criterion:

- do not let the evaluator walk semantic expression trees per trial
- do not let runtime component masking return through the back door
- do not postpone low-level runtime layout until after evaluator code exists

## Immediate Implication

The rebuild starts with this shape:

- one semantic compiler
- one projection pass
- one direct backend
- one exact backend

If the implementation drifts toward:

- one huge generic evaluator
- hidden runtime component masking
- route-specific trigger or tie patches

then the rebuild is already failing its own Phase 0 decisions.
