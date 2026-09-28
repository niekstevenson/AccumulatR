.normalize_trial_df <- function(trial_df, n_trials, acc_ids, comp_ids) {
  if (is.null(trial_df)) {
    return(list(
      component = rep.int(NA_character_, n_trials),
      onset = NULL
    ))
  }

  trial_df <- as.data.frame(trial_df)
  if (!all(c("trials", "racer") %in% names(trial_df))) {
    stop("trial_df must contain trials and racer columns", call. = FALSE)
  }
  n_acc <- length(acc_ids)
  expected_trials <- rep.int(seq_len(n_trials), rep.int(n_acc, n_trials))
  if (!is.numeric(trial_df$trials) ||
      length(trial_df$trials) != length(expected_trials) ||
      any(!is.finite(trial_df$trials)) ||
      any(trial_df$trials != expected_trials) ||
      !identical(as.character(trial_df$racer), rep(acc_ids, times = n_trials))) {
    stop(
      "trial_df must contain one complete trials-by-accumulator block in model order",
      call. = FALSE
    )
  }

  component <- rep.int(NA_character_, n_trials)
  if ("component" %in% names(trial_df)) {
    values <- as.character(trial_df$component)
    component <- values[seq.int(1L, length(values), by = n_acc)]
    if (!identical(values, rep.int(component, rep.int(n_acc, n_trials)))) {
      stop("trial_df$component must be constant within each trial", call. = FALSE)
    }
    unknown <- unique(component[!is.na(component) & !component %in% comp_ids])
    if (length(unknown)) {
      stop("Unknown trial_df component: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }

  onset <- NULL
  if ("onset" %in% names(trial_df)) {
    if (!is.numeric(trial_df$onset) || any(!is.na(trial_df$onset) & !is.finite(trial_df$onset))) {
      stop("trial_df$onset must be numeric", call. = FALSE)
    }
    onset <- as.numeric(trial_df$onset)
  }

  list(component = component, onset = onset)
}

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
  acc_ids <- names(prep$accumulators)
  n_acc <- length(acc_ids)
  n_params <- max(vapply(prep$accumulators, function(acc) {
    length(dist_param_names(acc$dist))
  }, integer(1)))
  columns <- c("q", "t0", paste0("p", seq_len(n_params)))
  if (!is.matrix(params_df) || !is.numeric(params_df) ||
      nrow(params_df) %% n_acc != 0L ||
      !identical(colnames(params_df)[seq_along(columns)], columns)) {
    stop("params_df must be a canonical trials-by-accumulator parameter matrix from build_param_matrix()", call. = FALSE)
  }
  storage.mode(params_df) <- "double"
  trial <- .normalize_trial_df(trial_df, nrow(params_df) %/% n_acc,
                               acc_ids, prep$components$ids)
  if (!is.null(seed)) set.seed(seed)
  keep_component <- keep_component %||% (prep$components$mode != "sample")
  simulate_cpp(prep, structure$simulation, params_df,
               match(trial$component, prep$components$ids), trial$onset,
               isTRUE(keep_detail),
               isTRUE(keep_component) && length(prep$components$ids) > 1L)
}
