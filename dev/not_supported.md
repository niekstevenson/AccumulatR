# Model and observation limits

This page describes the restrictions enforced by model finalization, context
construction, and data preparation. See the [logical-rule vignette](../vignettes/logical_rules.Rmd)
and [ranked-response vignette](../vignettes/multi_outcome.Rmd) for supported usage.

## Absence conditions need an event

`none_of()` is a condition evaluated when another rule finishes. It cannot
produce an outcome by itself or serve as a standalone `first_of()` branch.
For example, `first_of("go", none_of("stop"))` has no response time for its
absence branch and is rejected during likelihood context construction.

Use `all_of("go", none_of("stop"))` for a response that requires `go` to
finish before `stop`. A complete guarded rule can appear inside `first_of()`,
as in `first_of(all_of("go", none_of("stop")), "alternative")`.

## Ranked responses

With `n_outcomes > 1`, each outcome must directly name an accumulator or pool.
Chained onsets are supported. The following are unsupported:

- logical outcome expressions such as `all_of()`, `first_of()`, or `inhibit()`;
- `guess` and `map_outcome_to` observation options;
- multiple labels for the same deterministic event source, including a
  singleton pool that aliases an accumulator;
- censoring or truncation of ranked observations.

The requested rank count cannot exceed the number of declared outcomes.
Data must contain a first response/time pair. Subsequent pairs may be jointly
missing, but observed ranks must be consecutive, have distinct labels, and
have strictly increasing times. Simultaneous ranked observations are unsupported.

## Missing observations

A finite response time requires a response label. In a single-response model,
`R = NA` and `rt = NA` represent no observed response. A known label with a
missing response time requires a censoring code or a model with an observation
rule such as guessing or remapping. A censored trial can retain its known
response label or use `NA` when the label is unknown.

## Timing dependencies

`after()` can refer to an accumulator or pool. Outcome labels cannot be onset
sources, and onset dependencies must be acyclic. Its additional lag must be
finite and nonnegative.
