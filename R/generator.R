# ---- Public API ----------------------------------------------------------------

#' Simulate behavioral data from a model
#'
#' @param structure Finalized model structure.
#' @param params_df Canonical parameter matrix from `build_param_matrix()`.
#' @param trial_df Optional trial-by-accumulator conditioning table containing
#'   complete `trials`/`racer` blocks in model order. It may supply per-trial
#'   `component` and per-accumulator `onset` values.
#' @param seed Optional random-number seed.
#' @param keep_detail If `TRUE`, include latent source times and outcome candidates.
#'   Inactive or unreachable sources have infinite completion times.
#' @param keep_component Whether to keep the chosen mixture component in the
#'   output when the model has multiple components. If `NULL`, fixed mixtures
#'   keep the component label and sampled mixtures drop it.
#' @return A data frame of simulated behavioral data. For standard models this
#'   includes `trials`, `R`, and `rt`. If `n_outcomes > 1`, additional ordered
#'   response columns such as `R2`/`rt2` are included.
#' @details Simulation uses R's random-number generator from C++. Its compiled
#'   plan is cached on the finalized model and rebuilt after serialization.
#'   Execution trusts the canonical parameter matrix and complete conditioning
#'   table supplied by the caller; it does not validate them.
#'   Seeds are reproducible within this implementation, but do not reproduce
#'   the former R simulator's random-number sequence.
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
