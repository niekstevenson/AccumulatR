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

  ranked_columns <- grep("^(R|rt)[0-9]+$", nm, value = TRUE)
  indices <- suppressWarnings(as.integer(sub("^(R|rt)", "", ranked_columns)))
  bad <- sort(unique(indices[is.finite(indices) & indices > max_rank]))
  if (length(bad)) {
    stop(
      "Ranked observation columns must be contiguous from R2/rt2 with no gaps. Unexpected columns detected at ranks: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  list(max_rank = max_rank)
}

.has_observation_wrappers <- function(prep) {
  any(vapply(prep$outcomes, function(outcome) {
    options <- outcome$options
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
  if (max_rank <= 1L) {
    return(invisible(NULL))
  }
  for (start in seq_len(nrow(data_df))) {
    label1 <- as.character(data_df$R[[start]])
    rt1 <- data_df$rt[[start]]
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
  label_missing <- is.na(data_df$R)
  time_missing <- is.na(data_df$rt)
  censored <- if ("missingness" %in% names(data_df)) {
    !is.na(data_df$missingness)
  } else {
    FALSE
  }
  if (any(label_missing & !time_missing)) {
    stop("finite RT with missing response label is not supported", call. = FALSE)
  }
  if ((!allow_missing_all && any(label_missing)) ||
      (!allow_missing_rt && any(!label_missing & time_missing & !censored))) {
    stop(
      "Identity observations require a finite first-rank R/rt pair for every trial",
      call. = FALSE
    )
  }
  if (any(!time_missing & !is.finite(data_df$rt))) {
    stop("Observed RT values must be finite", call. = FALSE)
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

#' Prepare behavioral data for likelihood evaluation
#'
#' Check response labels, times, and observation conditions, and arrange the
#' data into one row per accumulator per trial for [log_likelihood()].
#'
#' @param structure Finalized model structure.
#' @param data_df Data frame with response labels `R` and numeric response times
#'   `rt`, usually one row per trial. Labels must match the model's outcomes.
#'   Optional columns specify `component`, `onset`, ranked responses (`R2`,
#'   `rt2`, and so on), or censoring and truncation; see Details.
#' @param compress If `TRUE`, collapse repeated prepared trials and attach
#'   an `expand` index mapping original trials to retained trials. Use only
#'   when repeated observations also share parameter values. Supply parameters
#'   and `ok` for the retained trials; [log_likelihood()] expands the returned
#'   values to the original trial order.
#' @return A data frame of class `accumulatr_data`, with accumulator rows grouped
#'   by trial and ordered as in the model. Response and component labels are
#'   factors with model-defined levels; attributes store the likelihood layout.
#' @details
#' **Trial layout.** Without a `racer` column, each row is one trial. Trials
#' are numbered consecutively in input order. To supply accumulator-specific
#' onsets, include `trials` and `racer` columns with one complete accumulator
#' block per trial in model order. Repeat each trial's response and component
#' values across its block. An `onset` replaces a fixed onset or adds an offset
#' to a chained onset. Omit it to use model defaults.
#'
#' **Mixtures.** A nonmissing `component` label conditions on that component.
#' An omitted or `NA` label averages over the model's mixture probabilities.
#'
#' **Missing observations.** For single-response models, `R = NA` and `rt = NA`
#' denote no observed response. A finite `rt` requires a response label.
#' A known response with `rt = NA` requires a censoring code or a model with
#' an observation rule such as `guess` or `map_outcome_to`.
#'
#' **Ranked responses.** Supply paired columns `R2`/`rt2`, `R3`/`rt3`, and so
#' on. The first pair must be observed. Labels cannot repeat, and observed
#' times must increase strictly. Later pairs may both be `NA`; all subsequent
#' ranks must then also be missing. Ranked observations support neither
#' censoring/truncation nor guessing/remapping rules.
#'
#' **Censoring and truncation.** `LT` and `UT` define the observation window.
#' An active truncation window conditions the likelihood on an observable
#' response with \eqn{LT \le rt \le UT}. For a censored trial set `rt = NA`
#' and use:
#' - `missingness = 1` for \eqn{LT \le rt < LC};
#' - `missingness = 2` for \eqn{UC < rt \le UT};
#' - `missingness = 3` for the union of those intervals.
#'
#' `R` may retain a known response or be `NA`. Uncensored trials use
#' `missingness = NA`. Missing bounds default to `LT = LC = 0` and
#' `UC = UT = Inf`. These bounds describe the recorded response time.
#' Times exactly at `LC` or `UC` remain uncensored.
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
prepare_data <- function(structure, data_df, compress = FALSE) {
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain R/rt per trial", call. = FALSE)
  }
  prep_eval_base <- structure$prep
  data_df <- as.data.frame(data_df)
  racer_level <- "racer" %in% names(data_df)
  required_cols <- c("R", "rt")
  missing_cols <- setdiff(required_cols, names(data_df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Data frame must include columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!"trials" %in% names(data_df)) {
    if (racer_level) {
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
  acc_ids <- names(prep_eval_base$accumulators)
  n_accumulators <- length(acc_ids)
  if (!racer_level) {
    n_trials <- nrow(data_df)
    data_df$trials <- seq_len(n_trials)
    data_df <- data_df[rep(seq_len(n_trials), each = n_accumulators), , drop = FALSE]
    data_df$racer <- rep(acc_ids, times = n_trials)
    rownames(data_df) <- NULL
  } else {
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
    acc_onset <- vapply(prep_eval_base$accumulators, `[[`, numeric(1), "onset")
    data_df$onset <- rep(unname(acc_onset), times = n_trials)
  } else if (!is.numeric(data_df$onset) || any(!is.finite(data_df$onset))) {
    stop("Data column 'onset' must contain finite numbers", call. = FALSE)
  }
  outcome_levels <- unique(names(prep_eval_base$outcomes))
  for (rank in seq_len(rank_info$max_rank)) {
    r_col <- if (rank == 1L) "R" else paste0("R", rank)
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
      component_levels[[1L]]
    } else {
      NA_character_
    }
  }
  data_df$component <- .normalize_prepared_index_column(
    data_df$component,
    component_levels,
    "component"
  )
  has_wrappers <- .has_observation_wrappers(prep_eval_base)
  if (rank_info$max_rank > 1L && has_wrappers) {
    stop("ranked observations do not support observation wrappers", call. = FALSE)
  }
  trial_level_columns <- c(
    "component", "R", "rt", .observation_bound_columns, "missingness",
    if (rank_info$max_rank > 1L) {
      c(paste0("R", 2:rank_info$max_rank), paste0("rt", 2:rank_info$max_rank))
    }
  )
  if (racer_level) {
    .validate_trial_level_columns(data_df, trial_level_columns)
  }
  trial_data <- data_df[
    seq.int(1L, nrow(data_df), by = n_accumulators),
    intersect(trial_level_columns, names(data_df)), drop = FALSE
  ]
  allowed_by_component <- prep_eval_base$outcomes_by_component
  component_chr <- as.character(trial_data$component)
  for (rank in seq_len(rank_info$max_rank)) {
    r_col <- if (rank == 1L) "R" else paste0("R", rank)
    label_chr <- as.character(trial_data[[r_col]])
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
  .validate_first_rank_trials(
    trial_data,
    allow_missing_all = rank_info$max_rank == 1L,
    allow_missing_rt = rank_info$max_rank == 1L && has_wrappers
  )
  .validate_ranked_trials(trial_data, rank_info$max_rank)
  class(data_df) <- unique(c("accumulatr_data", class(data_df)))
  if (isTRUE(compress)) {
    data_df <- .compress_prepared_trials(
      data_df,
      n_accumulators
    )
  }
  .attach_prepared_layout_attrs(data_df, rank_info$max_rank)
}

#' Build a compiled likelihood context from a model
#'
#' Compile the model's response rules and dependencies for likelihood
#' evaluation. Reuse the context across candidate parameter values and datasets
#' for the same model. Prepare each dataset with [prepare_data()].
#'
#' @param structure Finalized model structure.
#' @param diagnostics If `TRUE`, collect model compilation statistics for
#'   [complexity_metrics()].
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
    outcome_labels = unique(names(prep$outcomes)),
    observed_outcome_labels = Reduce(
      union,
      prep$observed_outcomes_by_component,
      init = character(0)
    )
  ), class = "accumulatr_context")
}

#' Inspect the size of a compiled likelihood plan
#'
#' Report how many symbolic regions, numerical operations, and integration
#' kernels the model requires. Use these counts to investigate models whose
#' context construction or likelihood evaluation is expensive.
#'
#' @param context Context created with `make_context(diagnostics = TRUE)`.
#' @return A list containing a `variants` data frame and a `total` list.
#'   Each variant is a compiled component plan. Columns count symbolic regions
#'   and cells, compiled roots and nodes, and integral kernels. Fields beginning
#'   with `max_` report maxima; other total fields sum across variants.
#' @export
complexity_metrics <- function(context) {
  semantic_complexity_metrics_context_cpp(context$cpp)
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
#' Calculate the probability of each first response, integrated over response
#' time and averaged over mixture components. With multiple trial blocks in
#' `parameters`, return the mean probabilities across those blocks.
#'
#' @param context Context created with `make_context()`.
#' @param parameters Parameter matrix from [build_param_matrix()]. Use
#'   `n_trials = 1` for one set of response probabilities.
#' @param include_na If `TRUE`, include a `"NA"` entry when there is residual
#'   probability of no observed response.
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
    context$cpp,
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
#' Build the context and prepare the observations once, then reuse them when
#' evaluating candidate parameter matrices for the same model.
#'
#' @param context Context created with `make_context()`.
#' @param data Prepared data created with `prepare_data()`.
#' @param parameters Numeric matrix from [build_param_matrix()], with one
#'   accumulator block per prepared trial in matching order.
#' @param ok Optional logical vector with one value per prepared trial.
#'   `TRUE` evaluates that trial; `FALSE` assigns `min_ll`. These assigned
#'   values are included in the sum. For compressed data, use the retained
#'   trial order.
#' @param sum If `TRUE`, return the summed log-likelihood. If `FALSE`, return
#'   trial-wise log-likelihood values.
#' @param min_ll Minimum log-likelihood value used for excluded or impossible
#'   trials.
#' @return A summed log-likelihood by default, or a numeric vector of
#'   trial-wise log-likelihood values when `sum = FALSE`.
#' @details Response-time observations contribute densities, so a
#'   log-likelihood can be positive. Missing responses and censoring contribute
#'   probability masses according to the observation rules.
#'
#'   Use [prepare_data()] and [build_param_matrix()] for the same model.
#'   Evaluation assumes matching layouts and valid parameters and does not
#'   repeat preparation checks.
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
  value <- semantic_loglik_context_cpp(
    context$cpp,
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
