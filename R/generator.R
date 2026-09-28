# ---- Public API ----------------------------------------------------------------

#' Simulate behavioral data from a model
#'
#' Generate one trial for each accumulator block in `params_df`, using the
#' model's response rules, timing dependencies, triggers, and mixture weights.
#'
#' @param structure Finalized model structure.
#' @param params_df Parameter matrix from [build_param_matrix()], with one row
#'   per accumulator per trial, in the order accumulators were added to the model.
#' @param trial_df Optional trial-by-accumulator conditioning table containing
#'   complete `trials`/`racer` blocks in model order. It may supply per-trial
#'   `component` and per-accumulator `onset` values. Repeat a component label
#'   across all rows of its trial; `NA` draws a component from the mixture.
#'   An onset value replaces a fixed onset or adds a delay to a chained onset;
#'   `NA` uses the onset specified in the model.
#' @param seed Optional seed passed to [set.seed()]. If `NULL`, use the current
#'   state of R's random-number generator.
#' @param keep_detail If `TRUE`, attach a `details` list containing latent
#'   source times and outcome candidates for each trial. Inactive or unreachable
#'   sources have infinite completion times.
#' @param keep_component Whether to keep the chosen mixture component in the
#'   output when the model has multiple components. If `NULL`, fixed mixtures
#'   keep the component label and sampled mixtures drop it.
#' @return A data frame with one row per trial and columns `trials`, `R`, and
#'   `rt`. A trial with no observed response has `R = NA` and `rt = NA`.
#'   If `n_outcomes > 1`, additional ordered response pairs such as `R2`/`rt2`
#'   are included; unobserved later ranks are `NA`.
#' @details Use the same model to build `params_df` and simulate the data.
#'   Supply valid parameter values and keep the matrix and conditioning table
#'   in matching row order; simulation does not check their layout or domains.
#' @seealso [prepare_data()], [log_likelihood()], [set_mixture()]
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params <- c(m = 0, s = 0.1)
#' df <- build_param_matrix(structure, params, n_trials = 3)
#' simulate(structure, df, seed = 123)
#' @export
simulate <- function(structure,
                     params_df,
                     trial_df = NULL,
                     seed = NULL,
                     keep_detail = FALSE,
                     keep_component = NULL) {
  prep <- structure$prep
  component <- if (is.null(trial_df$component)) NULL else {
    first_rows <- seq.int(1L, nrow(params_df), by = length(prep$accumulators))
    match(trial_df$component[first_rows], prep$components$ids)
  }
  onset <- if (is.null(trial_df$onset)) NULL else as.numeric(trial_df$onset)
  if (!is.null(seed)) set.seed(seed)
  keep_component <- keep_component %||% (prep$components$mode != "sample")
  simulate_cpp(prep, structure$simulation, params_df,
               component, onset,
               isTRUE(keep_detail),
               isTRUE(keep_component) && length(prep$components$ids) > 1L)
}
