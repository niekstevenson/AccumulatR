# Prepare behavioral data for likelihood evaluation

Check response labels, times, and observation conditions, and arrange
the data into one row per accumulator per trial for
[`log_likelihood()`](https://niekstevenson.github.io/AccumulatR/reference/log_likelihood.md).

## Usage

``` r
prepare_data(structure, data_df, compress = FALSE)
```

## Arguments

- structure:

  Finalized model structure.

- data_df:

  Data frame with response labels `R` and numeric response times `rt`,
  usually one row per trial. Labels must match the model's outcomes.
  Optional columns specify `component`, `onset`, ranked responses (`R2`,
  `rt2`, and so on), or censoring and truncation; see Details.

- compress:

  If `TRUE`, collapse repeated prepared trials and attach an `expand`
  index mapping original trials to retained trials. Use only when
  repeated observations also share parameter values. Supply parameters
  and `ok` for the retained trials;
  [`log_likelihood()`](https://niekstevenson.github.io/AccumulatR/reference/log_likelihood.md)
  expands the returned values to the original trial order.

## Value

A data frame of class `accumulatr_data`, with accumulator rows grouped
by trial and ordered as in the model. Response and component labels are
factors with model-defined levels; attributes store the likelihood
layout.

## Details

**Trial layout.** Without a `racer` column, each row is one trial.
Trials are numbered consecutively in input order. To supply
accumulator-specific onsets, include `trials` and `racer` columns with
one complete accumulator block per trial in model order. Repeat each
trial's response and component values across its block. An `onset`
replaces a fixed onset or adds an offset to a chained onset. Omit it to
use model defaults.

**Mixtures.** A nonmissing `component` label conditions on that
component. An omitted or `NA` label averages over the model's mixture
probabilities.

**Missing observations.** For single-response models, `R = NA` and
`rt = NA` denote no observed response. A finite `rt` requires a response
label. A known response with `rt = NA` requires a censoring code or a
model with an observation rule such as `guess` or `map_outcome_to`.

**Ranked responses.** Supply paired columns `R2`/`rt2`, `R3`/`rt3`, and
so on. The first pair must be observed. Labels cannot repeat, and
observed times must increase strictly. Later pairs may both be `NA`; all
subsequent ranks must then also be missing. Ranked observations support
neither censoring/truncation nor guessing/remapping rules.

**Censoring and truncation.** `LT` and `UT` define the observation
window. An active truncation window conditions the likelihood on an
observable response with \\LT \le rt \le UT\\. For a censored trial set
`rt = NA` and use:

- `missingness = 1` for \\LT \le rt \< LC\\;

- `missingness = 2` for \\UC \< rt \le UT\\;

- `missingness = 3` for the union of those intervals.

`R` may retain a known response or be `NA`. Uncensored trials use
`missingness = NA`. Missing bounds default to `LT = LC = 0` and
`UC = UT = Inf`. These bounds describe the recorded response time. Times
exactly at `LC` or `UC` remain uncensored.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params_df <- build_param_matrix(
  structure,
  c(m = 0, s = 0.1),
  n_trials = 2
)
data_df <- simulate(structure, params_df, seed = 1)
prepare_data(structure, data_df)
#>   trials     R        rt racer onset   component
#> 1      1 A_win 0.9679031     A     0 __default__
#> 2      2 A_win 0.9198333     A     0 __default__
```
