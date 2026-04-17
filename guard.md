# Guarded Shared-Gate Example

This note pins down the likelihood for the guarded/shared-gate model used in `example_16_guard_tie_simple`.

## Model

```r
model <- race_spec() |>
  add_accumulator("go_fast", "lognormal") |>
  add_accumulator("go_slow", "lognormal") |>
  add_accumulator("gate_shared", "lognormal") |>
  add_accumulator("stop_control", "lognormal") |>
  add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
  add_outcome("Slow", all_of("go_slow", "gate_shared")) |>
  finalize_model()
```

Write

- `F = go_fast`
- `S = go_slow`
- `G = gate_shared`
- `C = stop_control`

and let `f_X` / `F_X` denote the density / cdf of `X`.

The observed outcomes are

- `Fast = inhibit(all_of(F, G), by = C)`
- `Slow = all_of(S, G)`

For `Fast` to be observed at time `t`, three things must hold:

1. `C > t`
2. `all_of(F, G)` finishes at `t`
3. `Slow` does not beat `Fast`

The nontrivial part is the third condition. If `G = t` and both `F < t` and `S < t`, then both top-level outcomes occur at the same observed time `t`. In that tied case, `Fast` wins only if `F < S`.

## Hand-Written Density

For this model, the exact observed density for `Fast` at time `t` is

```text
f_Fast(t) =
S_C(t) [
  f_F(t) F_G(t) S_S(t)
  + f_G(t) F_F(t) S_S(t)
  + f_G(t) ∫_0^t f_F(u) (F_S(t) - F_S(u)) du
]
```

The three terms are:

- `f_F(t) F_G(t) S_S(t)`  
  `F` finishes last at `t`, `G < t`, and `Slow` has not finished by `t`.

- `f_G(t) F_F(t) S_S(t)`  
  `G` finishes last at `t`, `F < t`, and `Slow` has not finished by `t`.

- `f_G(t) ∫_0^t f_F(u) (F_S(t) - F_S(u)) du`  
  `G` finishes last at `t`, both `F` and `S` are already done before `t`, and `Fast` wins the same-time tie because `F < S`.

The blocker contributes only through `S_C(t)`, because `Fast` requires `C > t`.

## Numeric Check

These are the parameters from `example_16_guard_tie_simple`:

```r
pars <- c(
  go_fast.m = log(0.28),
  go_fast.s = 0.18,
  go_slow.m = log(0.34),
  go_slow.s = 0.18,
  gate_shared.m = log(0.30),
  gate_shared.s = 0.16,
  stop_control.m = log(0.27),
  stop_control.s = 0.15
)
```

At `t = 0.30`, the following reproduces the hand-written density and compares it to the framework likelihood.

```r
library(AccumulatR)

model <- race_spec() |>
  add_accumulator("go_fast", "lognormal") |>
  add_accumulator("go_slow", "lognormal") |>
  add_accumulator("gate_shared", "lognormal") |>
  add_accumulator("stop_control", "lognormal") |>
  add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
  add_outcome("Slow", all_of("go_slow", "gate_shared")) |>
  finalize_model()

pars <- c(
  go_fast.m = log(0.28),
  go_fast.s = 0.18,
  go_slow.m = log(0.34),
  go_slow.s = 0.18,
  gate_shared.m = log(0.30),
  gate_shared.s = 0.16,
  stop_control.m = log(0.27),
  stop_control.s = 0.15
)

t <- 0.30

f_fast <- function(x) dlnorm(x, pars["go_fast.m"], pars["go_fast.s"])
F_fast <- function(x) plnorm(x, pars["go_fast.m"], pars["go_fast.s"])
f_slow <- function(x) dlnorm(x, pars["go_slow.m"], pars["go_slow.s"])
F_slow <- function(x) plnorm(x, pars["go_slow.m"], pars["go_slow.s"])
f_gate <- function(x) dlnorm(x, pars["gate_shared.m"], pars["gate_shared.s"])
F_gate <- function(x) plnorm(x, pars["gate_shared.m"], pars["gate_shared.s"])
S_stop <- function(x) 1 - plnorm(x, pars["stop_control.m"], pars["stop_control.s"])

hand <- S_stop(t) * (
  f_fast(t) * F_gate(t) * (1 - F_slow(t)) +
  f_gate(t) * F_fast(t) * (1 - F_slow(t)) +
  f_gate(t) * integrate(
    function(u) f_fast(u) * (F_slow(t) - F_slow(u)),
    lower = 0,
    upper = t,
    rel.tol = 1e-10,
    abs.tol = 0
  )$value
)

prepared <- prepare_data(model, data.frame(R = "Fast", rt = t))
ctx <- make_context(model)
params_df <- build_param_matrix(model, pars, trial_df = prepared)
framework <- exp(log_likelihood(ctx, prepared, params_df)[[1]])

comparison <- data.frame(
  method = c("hand_written", "framework"),
  density = c(hand, framework)
)
comparison$diff_from_hand <- comparison$density - hand
comparison
```

With the current exact guarded-pair route, this gives:

```text
       method   density diff_from_hand
1 hand_written 2.0929288      0.0000000
2    framework 2.0929288      0.0000000
```

## Issue

The generic observed-likelihood route is not exact for this shape.

The reason is mathematical, not numerical: the generic route treats the competitor through a marginal non-winning calculation, but this model has positive-mass same-time ties created by the shared gate `G`. The third term above is exactly that tie mass. If it is omitted or handled only through strict survival, the `Fast` density is wrong.

For this model and these parameters at `t = 0.30`:

- hand-written density: `2.0929287888`
- generic observed route before the fix: `1.8598091277`
- current exact guarded-pair route: `2.0929287888`

So this example is a concrete counterexample showing why a generic “target density times competitor non-occurrence” route is not enough once same-time ties have positive mass.
