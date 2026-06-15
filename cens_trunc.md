# Censoring and Truncation Likelihood Plan

## Goal

Add likelihood support for optional trial-level censoring and truncation columns:

- `LT`: lower truncation bound
- `UT`: upper truncation bound
- `LC`: lower/left censoring cutoff
- `UC`: upper/right censoring cutoff

If none of these columns are present, likelihood evaluation must follow the
current path with no material slowdown. If any are present, they are interpreted
as observation-layer metadata for first-rank, unranked observations only.

This is a likelihood-only change. Simulation, response probabilities, and
ranked censoring are out of scope.

## Non-negotiable Design Constraint

Do not implement censoring/truncation as a special independent-race
product-of-survivors path. That shortcut is only valid for plain independent
first-finish races. The package already supports structured outcomes,
inhibition, remapping, guessing, mixtures, and terminal no-response. The
general implementation must use the existing exact outcome density machinery
and add interval probabilities to the observation layer.

The right abstraction is:

```text
semantic outcome density
  -> interval probability over outcome time
  -> observed-label / unknown-label observation plan
  -> optional truncation normalization
  -> component mixture marginalization
```

## Data Semantics

The columns are optional. Missing columns mean no censoring/truncation:

- absent `LT`: lower truncation defaults to `0`
- absent `UT`: upper truncation defaults to `Inf`
- absent `LC`: no lower/left censoring
- absent `UC`: no upper/right censoring

For present columns, `NA` means the feature is absent on that row. `LT` and
`UT` may therefore be supplied globally with repeated values, or row-wise.
`LC` and `UC` may also be repeated design cutoffs on exact-response rows.
Their presence does not by itself mean that the row is censored.

The columns are trial-level. After `prepare_data()` expands accumulator rows,
their values must be constant within a trial.

### Exact RT Rows

If `rt` is finite, the row is an exact observation. The numerator is the
existing observed density at `rt`, even if `LC` or `UC` columns are present as
design cutoffs.

If `LT`/`UT` are active, this density is conditioned on the truncation window:

```text
L = density(R, rt) / P(finite observable response in [LT, UT])
```

The exact `rt` must lie inside `[LT, UT]`. Reject the prepared data if this is
violated; do not silently assign a tiny likelihood to structurally invalid data.

If censoring cutoffs are present, the exact `rt` must also lie in the exact
observable region:

- if `LC` is finite, exact `rt` must be greater than `LC`;
- if `UC` is finite, exact `rt` must be less than `UC`.

Boundary equality has zero mass for the currently supported continuous-time
likelihood and should be rejected rather than treated as a special case.

### Censored Rows

A row is an interval observation when `rt` is missing and either a censoring
cutoff is non-missing or truncation bounds are active. If `rt` is missing and
neither censoring nor truncation is active, keep the current missing-RT /
missing-all observation semantics.

The observation interval is:

```text
lower = max(LT, UC) when UC is present, otherwise LT
upper = min(UT, LC) when LC is present, otherwise UT
```

Interpretations:

- `UC` only: right/upper censored, `T >= UC`, bounded above by `UT`
- `LC` only: left/lower censored, `T <= LC`, bounded below by `LT`
- both `UC` and `LC`: interval censored, `UC <= T <= LC`, intersected with
  the truncation window
- truncation only with missing `rt`: interval observed over `[LT, UT]`

Reject rows where the resulting interval is empty, except for exact boundary
equality if a later implementation deliberately supports point-mass models.
The current continuous-time likelihood should reject `upper <= lower`.

If `R` is present, the interval probability is restricted to that observed
label after applying the existing observation remapping/guessing semantics.

If `R` is missing, the interval probability is over all semantic finite
outcomes for the component: the latent winner is unknown. This should not be
implemented as the existing "missing all" complement, because that represents
the current missing-response observation policy, not interval censoring.

For an unknown-label right-censored row with `upper = Inf`, include terminal
no-response mass through the existing terminal no-response primitive. This is
the general analogue of intrinsic omissions: `T = Inf` satisfies `T >= UC`.
Do not include terminal no-response for finite upper bounds or label-specific
censored rows.

## Rejected Combinations

Reject these combinations in `prepare_data()` before entering C++:

- ranked observations (`R2`/`rt2` or higher) with any of `LT`, `UT`, `LC`, `UC`
  present;
- non-numeric censor/truncation columns;
- negative finite bounds;
- finite `UT <= LT`;
- both `UC` and `LC` present with `LC <= UC`;
- exact `rt` outside an active `[LT, UT]` truncation window;
- exact `rt` on the censored side of a finite `LC` or `UC` cutoff;
- non-missing `LC`/`UC` that produce an empty interval after intersecting
  truncation bounds.

If a later implementation hits another combination that requires stateful
ranked conditioning or a new observation model, reject it explicitly. Do not
add fallbacks that pretend to support it.

## R Data Layer

Modify `R/likelihood_param_interface.R`.

### Preparation

Add a small helper, for example:

```r
.detect_observation_bounds <- function(data_df)
```

It should:

- detect which of `LT`, `UT`, `LC`, `UC` are present;
- coerce present columns to numeric;
- validate per-row bounds on the trial-level input;
- leave absent columns absent.

Do not add synthetic bound columns when none are present. The absent case must
remain the current prepared-data shape.

### Layout Attributes

Do not add a separate `has_observation_bounds` attribute. It is redundant.

Extend the existing `layout_cols` attribute to include bound columns:

```r
attr(data_df, "layout_cols") <- c(
  component = ...,
  onset = ...,
  LT = ...,
  UT = ...,
  LC = ...,
  UC = ...
)
```

Use `NA_integer_` for absent bound columns. C++ can compute
`layout.bounds.present` from these indices. This avoids one extra attribute and
keeps all native column lookups in the existing layout mechanism.

Older prepared data that lacks these new `layout_cols` names naturally reads as
"no bounds" through the existing named lookup behavior. Do not force a rebuild
only to add absent metadata.

### Trial-Level Validation

Add `LT`, `UT`, `LC`, and `UC` to `trial_level_columns` when present, so
expanded accumulator rows cannot disagree within a trial.

The existing `.validate_first_rank_trials()` currently permits missing `rt`
only under observation wrappers. Extend that logic so missing `rt` is also
allowed when the row has an observation interval from censoring or truncation.

## C++ Prepared Data Layout

Modify `src/eval/trial_data.hpp`.

Extend `PreparedTrialLayout` with a bound layout:

```cpp
struct PreparedBoundColumnView {
  int lt_col{-1};
  int ut_col{-1};
  int lc_col{-1};
  int uc_col{-1};
  bool present{false};
};
```

or equivalent. Keep it compact and direct.

Add helpers to read bound values for a trial row:

```cpp
struct ObservationBounds {
  bool active{false};
  bool censored{false};
  bool truncation_active{false};
  bool right_censored{false};
  bool left_censored{false};
  double trunc_lower{0.0};
  double trunc_upper{R_PosInf};
  double event_lower{0.0};
  double event_upper{R_PosInf};
};
```

If `layout.bounds.present` is false, the evaluator should not touch bound
columns.

## Exact Interval Probability Primitive

Modify `src/eval/exact_sequence.hpp`.

Add:

```cpp
double exact_outcome_probability_between(
    const ExactVariantPlan &plan,
    const ParamView &params,
    int first_param_row,
    semantic::Index target_idx,
    double lower,
    double upper,
    ExactStepWorkspace *workspace);
```

Behavior:

- invalid target or empty interval returns `0`;
- clamp finite lower to `0`;
- finite upper uses `quadrature::integrate_finite_default(lower, upper, ...)`;
- infinite upper uses the existing positive-tail quadrature with a lower shift:
  evaluate density at `lower + tail_node`;
- return a clamped finite probability.

Keep `exact_finite_outcome_probability()` as the special case `[0, Inf]`, or
rewrite it to call the new helper. The latter is cleaner if it does not change
current behavior.

Add a helper for all-outcome finite mass:

```cpp
double exact_all_outcome_probability_between(...);
```

This should sum semantic outcome probabilities once per outcome, not once per
observed-label branch, to avoid double-counting guess/remap branches when the
winner is unknown.

## Observation Interval Evaluation

Modify `src/eval/observation_model.hpp` and
`src/eval/observation_plan_eval.hpp`.

Do not replace the current scalar `observed_rt` point-likelihood evaluator.
That would put interval branching in the normal no-bound path for no benefit.
Keep `evaluate_observation_plan_at_row()` unchanged for point RTs and existing
missing-RT observations.

Add a separate interval runtime type and evaluator used only when
`layout.bounds.present` and the current row is censored:

```cpp
struct ObservationInterval {
  double lower{0.0};
  double upper{R_PosInf};
  bool include_terminal_no_response{false};
};
```

The interval evaluator can be a small direct helper instead of a new compiled
plan-op family:

```cpp
double evaluate_observed_label_interval_probability(...);
double evaluate_unknown_label_interval_probability(...);
```

This keeps point observations on the current evaluator and avoids adding plan
ops that only exist for row-dependent censoring.

### Branch Planning

For label-present interval rows, reuse the same observed-label branch buckets
used by missing-RT plans:

```text
keep_by_code[label] + missing_rt_by_code[label]
```

but evaluate each branch over `[lower, upper]`, not `[0, Inf]`.

For label-missing interval rows, use a separate list of semantic finite
outcomes for the component. Add this to `ComponentObservationPlan`; do not
reuse `missing_all_branches` or complement logic.

The existing `NoResponseProbability` primitive is sufficient for terminal
no-response mass. Do not add a second op with a broader name unless the
implementation actually needs a separate operation.

## Truncation Normalization

Truncation is a denominator on the observation likelihood.

For an observed component:

```text
log L = log(numerator_c) - log(window_mass_c)
```

For a latent component mixture:

```text
log L =
  log(sum_c weight_c * numerator_c)
  - log(sum_c weight_c * window_mass_c)
```

Do not normalize each component independently before mixture summation. That
would compute:

```text
sum_c weight_c * numerator_c / window_mass_c
```

which is the wrong conditional likelihood when component membership is latent.

The window mass is the probability of the observable sample space under the
component:

- finite truncation window: all finite semantic outcomes in `[LT, UT]`;
- `UT = Inf`: all finite semantic outcomes in `[LT, Inf]`;
- terminal no-response is not part of this finite-response truncation
  denominator.

Unknown-label upper censoring to infinity may include terminal no-response in
the numerator when no truncation normalization is active. If truncation is also
active and terminal no-response would need to be included, reject the
combination until the package has an explicit extended-real sample-space
definition. Do not hide this ambiguity in a fallback.

This denominator must be computed only when truncation is active. If no bound
columns are present, the current likelihood path must remain unchanged.

## Trial Loop Changes

Modify `src/eval/observation_trial_loop.hpp`.

The current fast identity branch is valid only when:

- observation is identity;
- all selected trials have observed components;
- all selected trials have finite first-rank `R`/`rt`;
- no observation bound columns are present.

Add the bound-layout check there. This preserves current speed in the common
case and prevents exact finite-RT evaluation from ignoring truncation.

In the generic observation loop:

1. Read bounds only if `layout.bounds.present`.
2. Classify the row as point or interval.
3. Resolve component choices as today.
4. For each component:
   - evaluate numerator using the point plan or direct interval evaluator;
   - evaluate truncation denominator only if needed;
   - accumulate weighted numerator and weighted denominator separately when
     component membership is latent.
5. Convert to `min_ll` only after the mathematically relevant probability is
   non-positive or non-finite.

Do not route no-bound rows through interval machinery.

## Caching and Efficiency

No-censor/no-truncation evaluation must not pay for interval logic:

- no synthetic columns;
- a single `layout.bounds.present` branch near the trial loop;
- identity fast path unchanged when no bounds are present;
- existing RT-free cache unchanged for bound-free missing-RT observations.

For bound-dependent rows, either disable the existing RT-free cache or extend
its key with:

- event lower;
- event upper;
- truncation lower;
- truncation upper;
- include-terminal flag.

Disabling only for bound-dependent rows is acceptable for the first
implementation. It is honest and avoids corrupt cache hits.

The interval primitive should reuse one `ExactStepWorkspace` per component
through the existing workspace pool. Do not allocate a workspace per quadrature
node.

## Tests Worth Adding

This change is risky enough to justify focused tests.

Add a small test file for first-rank censoring/truncation with a two-lognormal
plain race, comparing against direct R formulas:

- no bound columns exactly preserves existing likelihood;
- finite exact `rt` rows with repeated `LC`/`UC` cutoffs remain exact-density
  rows when the RT falls in the uncensored region;
- exact RT with `LT`/`UT` equals density divided by finite race mass;
- upper-censored unknown label equals race mass over `[UC, Inf]`;
- lower-censored unknown label equals race mass over `[0, LC]`;
- label-specific interval censoring equals the matching defective outcome
  integral;
- latent mixture truncation uses ratio of weighted sums, not sum of ratios;
- ranked observations with any bound column are rejected.

Keep these tests behavioral. Do not test internal op ordering or cache details.

## Implementation Order

1. Add R-side detection, coercion, validation, and layout attrs.
2. Extend `PreparedTrialLayout` and C++ bound readers.
3. Add exact interval probability primitives.
4. Add direct interval observation evaluators and the minimal branch metadata
   they need.
5. Modify the generic observation trial loop for interval numerator and
   truncation denominator handling.
6. Guard the identity fast path so it is used only when no bound columns exist.
7. Add focused tests for the supported cases and rejection cases.
8. Run package tests and compare a representative no-bound benchmark to confirm
   the absent-feature path did not regress materially.

## Expected Limitations After First Implementation

The first implementation should support unranked first-response likelihoods
only. It should reject ranked censoring/truncation. If another combination
requires conditioning on advanced sequence state or a different observation
sample space, reject it explicitly and document the limitation rather than
adding a partial fallback.
