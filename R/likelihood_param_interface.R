.validate_ranked_observation_columns <- function(data_df) {
  nm <- names(data_df)
  max_rank <- 1L
  rank <- 2L
  repeat {
    r_col <- paste0("R", rank)
    rt_col <- paste0("rt", rank)
    has_r <- r_col %in% nm
    has_rt <- rt_col %in% nm
    if (!has_r && !has_rt) {
      break
    }
    if (xor(has_r, has_rt)) {
      stop(sprintf("Ranked observations must provide paired columns '%s' and '%s'", r_col, rt_col), call. = FALSE)
    }
    max_rank <- rank
    rank <- rank + 1L
  }

  rank_r <- grep("^R[0-9]+$", nm, value = TRUE)
  if (length(rank_r) > 0L) {
    idx_r <- suppressWarnings(as.integer(sub("^R", "", rank_r)))
    idx_r <- idx_r[is.finite(idx_r)]
    if (length(idx_r) > 0L) {
      bad <- sort(unique(idx_r[idx_r > max_rank]))
      if (length(bad) > 0L) {
        stop(
          "Ranked observation columns must be contiguous from R2/rt2 with no gaps. Unexpected columns detected at ranks: ",
          paste(bad, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  rank_rt <- grep("^rt[0-9]+$", nm, value = TRUE)
  if (length(rank_rt) > 0L) {
    idx_rt <- suppressWarnings(as.integer(sub("^rt", "", rank_rt)))
    idx_rt <- idx_rt[is.finite(idx_rt)]
    if (length(idx_rt) > 0L) {
      bad <- sort(unique(idx_rt[idx_rt > max_rank]))
      if (length(bad) > 0L) {
        stop(
          "Ranked observation columns must be contiguous from R2/rt2 with no gaps. Unexpected rt columns detected at ranks: ",
          paste(bad, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  list(max_rank = max_rank)
}

.has_observation_wrappers <- function(prep) {
  outcome_defs <- prep$outcomes %||% list()
  any(vapply(outcome_defs, function(outcome) {
    options <- outcome$options %||% list()
    !is.null(options$guess) || !is.null(options$map_outcome_to)
  }, logical(1)))
}

.observation_bound_columns <- c("LT", "UT", "LC", "UC")

.prepare_observation_bounds <- function(data_df, max_rank) {
  columns <- intersect(c(.observation_bound_columns, "missingness"), names(data_df))
  if (length(columns) == 0L) {
    return(data_df)
  }
  for (column in intersect(.observation_bound_columns, columns)) {
    value <- data_df[[column]]
    if ((!is.numeric(value) && !(is.logical(value) && all(is.na(value)))) ||
        any(is.nan(value))) {
      stop(sprintf("Observation bound column '%s' must contain numbers or NA", column), call. = FALSE)
    }
    data_df[[column]] <- as.numeric(value)
  }
  if ("missingness" %in% columns) {
    missingness <- data_df$missingness
    if ((!is.numeric(missingness) &&
         !(is.logical(missingness) && all(is.na(missingness)))) ||
        any(is.nan(missingness)) ||
        any(!is.na(missingness) & !missingness %in% 1:3)) {
      stop("Observation column 'missingness' must contain 1, 2, 3, or NA", call. = FALSE)
    }
    missingness <- as.integer(missingness)
  } else {
    missingness <- rep.int(NA_integer_, nrow(data_df))
  }

  bound <- function(column, default) {
    if (!column %in% names(data_df)) {
      return(rep.int(default, nrow(data_df)))
    }
    out <- data_df[[column]]
    out[is.na(out)] <- default
    out
  }
  LT <- bound("LT", 0)
  UT <- bound("UT", Inf)
  LC <- bound("LC", 0)
  UC <- bound("UC", Inf)
  if (any(LT < 0 | UT < 0 | LC < 0 | UC < 0)) {
    stop("Observation bounds must be non-negative numbers", call. = FALSE)
  }
  if (any(!is.finite(LT)) || any(UT <= LT)) {
    stop("Truncation bounds require finite LT and UT > LT", call. = FALSE)
  }
  active_bounds <- !is.na(missingness) | LT > 0 | is.finite(UT)
  if (max_rank > 1L && any(active_bounds)) {
    stop("ranked observations do not support censoring or truncation", call. = FALSE)
  }

  rt <- data_df$rt
  censored <- !is.na(missingness)
  if (any(censored & !is.na(rt))) {
    stop("Censored observations must have rt = NA", call. = FALSE)
  }
  lower <- missingness %in% c(1L, 3L)
  upper <- missingness %in% c(2L, 3L)
  if (any(lower & (LC < LT | LC > UT))) {
    stop("Lower-censored observations require LT <= LC <= UT", call. = FALSE)
  }
  if (any(upper & (UC < LT | UC > UT))) {
    stop("Upper-censored observations require LT <= UC <= UT", call. = FALSE)
  }
  if (any(missingness == 3L & LC > UC, na.rm = TRUE)) {
    stop("Both-censored observations require LC <= UC", call. = FALSE)
  }

  exact <- !censored & !is.na(rt)
  if (any(exact & (rt < LT | rt > UT))) {
    stop("Exact RT values must lie inside the truncation window [LT, UT]", call. = FALSE)
  }
  data_df[.observation_bound_columns] <- list(LT, UT, LC, UC)
  data_df$missingness <- missingness
  data_df
}

.validate_trial_level_columns <- function(data_df, columns) {
  columns <- intersect(columns, names(data_df))
  trial <- data_df$trials
  first <- match(trial, trial)
  for (col in columns) {
    x <- data_df[[col]]
    reference <- x[first]
    same <- (is.na(x) & is.na(reference)) |
      (!is.na(x) & !is.na(reference) & x == reference)
    if (!all(same)) {
      stop(
        "Prepared data must keep trial-level column '", col,
        "' constant within each trial",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

.validate_ranked_trials <- function(data_df, max_rank) {
  if (max_rank <= 1L || nrow(data_df) == 0L) {
    return(invisible(NULL))
  }
  trial <- as.integer(data_df$trials)
  starts <- c(1L, which(trial[-1L] != trial[-nrow(data_df)]) + 1L)
  for (start in starts) {
    label1 <- as.character(data_df$R[[start]])
    rt1 <- data_df$rt[[start]]
    if (is.na(label1) || is.na(rt1) || !is.finite(rt1)) {
      stop(
        "Ranked observations must provide a finite first-rank R/rt pair for every trial",
        call. = FALSE
      )
    }
    seen_labels <- label1
    prev_rt <- rt1
    terminated <- FALSE
    for (rank in 2:max_rank) {
      r_col <- paste0("R", rank)
      rt_col <- paste0("rt", rank)
      label <- as.character(data_df[[r_col]][[start]])
      rt <- data_df[[rt_col]][[start]]
      label_missing <- is.na(label)
      time_missing <- is.na(rt)
      if (label_missing && time_missing) {
        terminated <- TRUE
        next
      }
      if (terminated) {
        stop(
          "Ranked observations must stop cleanly after the last observed rank",
          call. = FALSE
        )
      }
      if (label_missing || time_missing) {
        stop(
          sprintf("Ranked observations must provide paired values in '%s'/'%s' within each trial", r_col, rt_col),
          call. = FALSE
        )
      }
      if (!is.finite(rt)) {
        stop("Ranked observations must provide finite RT values for observed ranks", call. = FALSE)
      }
      if (!(rt > prev_rt)) {
        stop("Ranked observation times must be strictly increasing within each trial", call. = FALSE)
      }
      if (label %in% seen_labels) {
        stop("Ranked observations must not repeat the same outcome label within a trial", call. = FALSE)
      }
      seen_labels <- c(seen_labels, label)
      prev_rt <- rt
    }
  }
  invisible(NULL)
}

.validate_first_rank_trials <- function(data_df,
                                        allow_missing_all = FALSE,
                                        allow_missing_rt = FALSE) {
  if (nrow(data_df) == 0L) {
    return(invisible(NULL))
  }
  trial <- as.integer(data_df$trials)
  starts <- c(1L, which(trial[-1L] != trial[-nrow(data_df)]) + 1L)
  for (start in starts) {
    label_missing <- is.na(as.character(data_df$R[[start]]))
    rt <- data_df$rt[[start]]
    time_missing <- is.na(rt)
    censored <- "missingness" %in% names(data_df) &&
      !is.na(data_df$missingness[[start]])
    if (label_missing && !time_missing) {
      stop("finite RT with missing response label is not supported", call. = FALSE)
    }
    if (label_missing && time_missing) {
      if (!allow_missing_all) {
        stop(
          "Identity observations require a finite first-rank R/rt pair for every trial",
          call. = FALSE
        )
      }
      next
    }
    if (time_missing && !allow_missing_rt && !censored) {
      stop(
        "Identity observations require a finite first-rank R/rt pair for every trial",
        call. = FALSE
      )
    }
    if (!time_missing && !is.finite(rt)) {
      stop("Observed RT values must be finite", call. = FALSE)
    }
  }
  invisible(NULL)
}

.compress_prepared_trials <- function(data_df, n_accumulators) {
  n_rows <- nrow(data_df)
  trial_starts <- seq.int(1L, n_rows, by = n_accumulators)
  n_trials <- length(trial_starts)
  if (n_trials <= 1L) {
    return(data_df)
  }
  trial_ends <- trial_starts + n_accumulators - 1L
  sig_cols <- setdiff(names(data_df), "trials")
  signatures <- character(n_trials)
  for (i in seq_len(n_trials)) {
    block <- data_df[trial_starts[[i]]:trial_ends[[i]], sig_cols, drop = FALSE]
    rownames(block) <- NULL
    signatures[[i]] <- rawToChar(serialize(block, NULL, ascii = TRUE))
  }
  keep_trials <- which(!duplicated(signatures))
  if (length(keep_trials) == n_trials) {
    return(data_df)
  }
  keep_rows <- unlist(Map(seq.int, trial_starts[keep_trials], trial_ends[keep_trials]), use.names = FALSE)
  out <- data_df[keep_rows, , drop = FALSE]
  out$trials <- rep.int(
    seq_along(keep_trials),
    rep.int(n_accumulators, length(keep_trials))
  )
  rownames(out) <- NULL
  attr(out, "expand") <- as.integer(match(signatures, signatures[keep_trials]))
  class(out) <- class(data_df)
  out
}

.attach_prepared_layout_attrs <- function(data_df, max_rank) {
  rank_names <- c("R", if (max_rank > 1L) paste0("R", seq.int(2L, max_rank)))
  time_names <- c("rt", if (max_rank > 1L) paste0("rt", seq.int(2L, max_rank)))

  attr(data_df, "layout_cols") <- setNames(
    as.integer(match(c("component", "onset", .observation_bound_columns, "missingness"), names(data_df))),
    c("component", "onset", .observation_bound_columns, "missingness")
  )
  attr(data_df, "label_cols") <- as.integer(match(rank_names, names(data_df)))
  attr(data_df, "time_cols") <- as.integer(match(time_names, names(data_df)))
  attr(data_df, "max_rank") <- as.integer(max_rank)
  data_df
}

.prepare_data_structure <- function(structure, data_df, compress = FALSE) {
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain R/rt per trial", call. = FALSE)
  }
  prep_eval_base <- structure$prep
  data_df <- as.data.frame(data_df)
  required_cols <- c("R", "rt")
  missing_cols <- setdiff(required_cols, names(data_df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Data frame must include columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!"trials" %in% names(data_df)) {
    if ("racer" %in% names(data_df)) {
      stop("Racer-level data must include a 'trials' column", call. = FALSE)
    }
    data_df$trials <- seq_len(nrow(data_df))
  }
  if (!is.numeric(data_df$rt)) {
    stop("Data column 'rt' must be numeric", call. = FALSE)
  }
  data_df$rt <- as.numeric(data_df$rt)
  rank_info <- .validate_ranked_observation_columns(data_df)
  data_df <- .prepare_observation_bounds(data_df, rank_info$max_rank)
  if (!"racer" %in% names(data_df)) {
    data_df <- .expand_accumulator_rows(structure, data_df)
  } else {
    acc_ids <- names(prep_eval_base$accumulators)
    n_accumulators <- length(acc_ids)
    if (nrow(data_df) %% n_accumulators != 0L) {
      stop(
        "Racer-level data must contain one complete accumulator block per trial in model order",
        call. = FALSE
      )
    }
    n_trials <- nrow(data_df) %/% n_accumulators
    starts <- 1L + seq.int(0L, n_trials - 1L) * n_accumulators
    trial_ids <- as.character(data_df$trials)
    block_ids <- trial_ids[starts]
    if (anyNA(block_ids) || anyDuplicated(block_ids) ||
        !identical(trial_ids, rep(block_ids, each = n_accumulators)) ||
        !identical(as.character(data_df$racer), rep(acc_ids, times = n_trials))) {
      stop(
        "Racer-level data must contain one complete accumulator block per trial in model order",
        call. = FALSE
      )
    }
    data_df$trials <- rep.int(seq_len(n_trials), rep.int(n_accumulators, n_trials))
  }
  if (!"onset" %in% names(data_df)) {
    acc_defs <- prep_eval_base$accumulators %||% list()
    acc_onset <- vapply(acc_defs, function(a) a$onset %||% 0, numeric(1))
    acc_ids <- names(acc_defs)
    onset_map <- setNames(acc_onset, acc_ids)
    data_df$onset <- vapply(as.character(data_df$racer), function(acc) {
      onset_map[[acc]] %||% 0
    }, numeric(1))
  } else if (!is.numeric(data_df$onset) || any(!is.finite(data_df$onset))) {
    stop("Data column 'onset' must contain finite numbers", call. = FALSE)
  }
  outcome_levels <- unique(names(prep_eval_base$outcomes %||% list()))
  if (length(outcome_levels) == 0L) {
    stop("Model must define outcomes", call. = FALSE)
  }
  for (rank in seq_len(rank_info$max_rank)) {
    r_col <- if (rank == 1L) "R" else paste0("R", rank)
    if (!r_col %in% names(data_df)) {
      next
    }
    data_df[[r_col]] <- .normalize_prepared_index_column(
      data_df[[r_col]],
      outcome_levels,
      r_col
    )
  }
  if (rank_info$max_rank > 1L) {
    for (rank in 2:rank_info$max_rank) {
      column <- paste0("rt", rank)
      if (!is.numeric(data_df[[column]])) {
        stop("Data column '", column, "' must be numeric", call. = FALSE)
      }
      data_df[[column]] <- as.numeric(data_df[[column]])
    }
  }
  component_levels <- prep_eval_base$components$ids
  if (!"component" %in% names(data_df)) {
    data_df$component <- if (length(component_levels) <= 1L) {
      "__default__"
    } else {
      NA_character_
    }
  }
  data_df$component <- .normalize_prepared_index_column(
    data_df$component,
    component_levels,
    "component"
  )
  if (rank_info$max_rank > 1L && .has_observation_wrappers(prep_eval_base)) {
    stop("ranked observations do not support observation wrappers", call. = FALSE)
  }
  allowed_by_component <- prep_eval_base$outcomes_by_component
  component_chr <- as.character(data_df$component)
  for (rank in seq_len(rank_info$max_rank)) {
    r_col <- if (rank == 1L) "R" else paste0("R", rank)
    if (!r_col %in% names(data_df)) {
      next
    }
    label_chr <- as.character(data_df[[r_col]])
    bad <- !is.na(label_chr) & !is.na(component_chr) & !mapply(
      function(lbl, cid) lbl %in% (allowed_by_component[[cid]] %||% character(0)),
      label_chr,
      component_chr,
      USE.NAMES = FALSE
    )
    if (any(bad)) {
      bad_idx <- which(bad)[1L]
      stop(
        sprintf(
          "Outcome '%s' is not allowed for component '%s' in column '%s'",
          label_chr[[bad_idx]],
          component_chr[[bad_idx]],
          r_col
        ),
        call. = FALSE
      )
    }
  }
  trial_level_columns <- c(
    "component",
    "R",
    "rt",
    .observation_bound_columns,
    "missingness",
    unlist(lapply(seq.int(2L, rank_info$max_rank), function(rank) c(paste0("R", rank), paste0("rt", rank))), use.names = FALSE)
  )
  .validate_trial_level_columns(data_df, trial_level_columns)
  .validate_first_rank_trials(
    data_df,
    allow_missing_all = rank_info$max_rank == 1L,
    allow_missing_rt = rank_info$max_rank == 1L && .has_observation_wrappers(prep_eval_base)
  )
  .validate_ranked_trials(data_df, rank_info$max_rank)
  class(data_df) <- unique(c("accumulatr_data", class(data_df)))
  if (isTRUE(compress)) {
    data_df <- .compress_prepared_trials(
      data_df,
      length(prep_eval_base$accumulators)
    )
  }
  .attach_prepared_layout_attrs(data_df, rank_info$max_rank)
}

#' Prepare behavioral data for likelihood evaluation
#'
#' `prepare_data()` expands trial-level observations to the accumulator layout
#' expected by the compiled likelihood code and tags the result as trusted
#' likelihood input.
#'
#' @param structure Finalized model structure.
#' @param data_df Behavioral data. In the simplest case this contains `trials`,
#'   `R`, and `rt`; for multi-outcome models it can also contain `R2`, `rt2`,
#'   and so on. Optional `LT`/`UT` columns define a truncation window. For a
#'   censored trial, set `rt = NA` and use `missingness = 1` for `[LT, LC)`,
#'   `2` for `(UC, UT]`, or `3` for their union; `R` may retain the known
#'   response or be `NA`. Uncensored trials use `missingness = NA`. Missing
#'   bounds default to `LT = LC = 0` and `UC = UT = Inf`. Codes 1--3 apply to
#'   the response observation, not to individual accumulators. Censoring and
#'   truncation are not supported for ranked observations.
#' @param compress If `TRUE`, collapse repeated prepared trials and attach
#'   an `expand` index so `log_likelihood()` can return trial-level values on
#'   the original trial scale.
#'   Defaults to `FALSE`.
#' @return An `accumulatr_data` object.
#' @details The likelihood for an active truncation window is conditioned on a
#'   observable response in `[LT, UT]`. Censoring comparisons are strict, so
#'   observations exactly at `LC` or `UC` remain uncensored.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_param_matrix(
#'   structure,
#'   c(m = 0, s = 0.1),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' prepare_data(structure, data_df)
#' @export
prepare_data <- function(structure,
                         data_df,
                         compress = FALSE) {
  .prepare_data_structure(
    structure = structure,
    data_df = data_df,
    compress = compress
  )
}

#' Build a compiled likelihood context from a model
#'
#' A context stores compiled model/runtime state only. Behavioral data are
#' prepared separately with `prepare_data()` and supplied to
#' `log_likelihood()`.
#'
#' @param structure Finalized model structure.
#' @param diagnostics If `TRUE`, collect symbolic/compiled complexity metrics.
#' @return An `accumulatr_context` object.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' make_context(structure)
#' @export
make_context <- function(structure, diagnostics = FALSE) {
  prep <- structure$prep
  structure(list(
    cpp = semantic_make_likelihood_context_prep_cpp(prep, isTRUE(diagnostics)),
    n_accumulators = length(prep$accumulators),
    outcome_labels = unique(names(prep$outcomes)),
    observed_outcome_labels = Reduce(
      union,
      prep$observed_outcomes_by_component,
      init = character(0)
    )
  ), class = "accumulatr_context")
}

#' Return compiled exact complexity metrics
#'
#' @param context Context created with `make_context(diagnostics = TRUE)`.
#' @return A list with per-variant and total symbolic/compiled metrics.
#' @export
complexity_metrics <- function(context) {
  if (!isTRUE(context$cpp$has_complexity_metrics)) {
    stop(
      "complexity metrics were not collected; create the context with diagnostics = TRUE",
      call. = FALSE
    )
  }
  semantic_complexity_metrics_context_cpp(context$cpp$native)
}

.normalize_prepared_index_column <- function(x, levels, column_name) {
  if (!is.character(x) && !is.factor(x)) {
    stop("Prepared data column '", column_name, "' must contain labels", call. = FALSE)
  }
  labels <- as.character(x)
  unknown <- unique(labels[!is.na(labels) & !labels %in% levels])
  if (length(unknown)) {
    stop(
      "Unknown labels in prepared data column '", column_name, "': ",
      paste(utils::head(unknown, 5L), collapse = ", "),
      call. = FALSE
    )
  }
  factor(labels, levels = levels)
}

#' Evaluate marginal response probabilities
#'
#' `response_probabilities()` evaluates the model-implied marginal probability
#' of each observed response label for a compiled context and canonical
#' parameter matrix. Mixture components are marginalized according to the model.
#'
#' @param context Context created with `make_context()`.
#' @param parameters A canonical parameter matrix created by
#'   `build_param_matrix()`.
#' @param include_na If `TRUE`, include residual mass as `"NA"`.
#' @return A named numeric vector of marginal response probabilities. Names are
#'   observed outcome labels. When `include_na = TRUE`, a residual `"NA"` entry
#'   is included if the model assigns probability mass to unobserved or
#'   `NA`-mapped outcomes.
#' @examples
#' spec <- race_spec() |>
#'   add_accumulator("left", "lognormal") |>
#'   add_accumulator("right", "lognormal") |>
#'   add_outcome("left", "left") |>
#'   add_outcome("right", "right") |>
#'   set_parameters(separate = list(m = TRUE))
#'
#' model <- finalize_model(spec)
#' params <- build_param_matrix(
#'   model,
#'   c(left.m = log(0.25), right.m = log(0.40), s = 0.20),
#'   n_trials = 1
#' )
#'
#' response_probabilities(make_context(model), params)
#' @export
response_probabilities <- function(context, parameters, include_na = TRUE) {
  probability <- as.numeric(semantic_response_probabilities_context_cpp(
    context$cpp$native,
    parameters
  ))
  result <- stats::setNames(probability, context$outcome_labels)
  result <- result[intersect(names(result), context$observed_outcome_labels)]
  residual <- 1.0 - sum(result)
  if (include_na && residual > .Machine$double.eps) {
    result["NA"] <- residual
  }
  result
}

#' Evaluate log-likelihoods of behavioral data
#'
#' Compute the summed log-likelihood by default, or trial-wise log-likelihoods
#' when `sum = FALSE`.
#'
#' @param context Context created with `make_context()`.
#' @param data Prepared data created with `prepare_data()`.
#' @param parameters A canonical numeric parameter matrix created by
#'   `build_param_matrix()`.
#' @param ok Logical vector marking which trials should contribute to the
#'   likelihood. Trials marked `FALSE` are assigned `min_ll`.
#' @param sum If `TRUE`, return the summed log-likelihood. If `FALSE`, return
#'   trial-wise log-likelihood values.
#' @param min_ll Minimum log-likelihood value used for excluded or impossible
#'   trials.
#' @return A summed log-likelihood by default, or a numeric vector of
#'   trial-wise log-likelihood values when `sum = FALSE`.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_param_matrix(
#'   structure,
#'   c(m = 0, s = 0.1),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' prepared <- prepare_data(structure, data_df)
#' ctx <- make_context(structure)
#' log_likelihood(ctx, prepared, params_df)
#' @export
log_likelihood <- function(context,
                           data,
                           parameters,
                           ok = NULL,
                           sum = TRUE,
                           min_ll = log(1e-10)) {
  cpp_ctx <- context$cpp
  value <- semantic_loglik_context_cpp(
    cpp_ctx$native,
    parameters,
    data,
    ok,
    min_ll
  )
  if (sum) {
    base::sum(value)
  } else {
    value
  }
}
