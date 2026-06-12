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

.prepare_likelihood_prep <- function(structure, prep = NULL) {
  structure <- .as_model_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("model structure must include model_spec; rebuild with finalize_model")
  }
  prep_eval_base <- prep %||% structure$prep %||% NULL
  if (is.null(prep_eval_base)) {
    prep_eval_base <- prepare_model(structure$model_spec)
  }
  list(structure = structure, prep = prep_eval_base)
}

.required_p_slots_from_prep <- function(prep) {
  acc_defs <- prep$accumulators %||% list()
  if (length(acc_defs) == 0L) {
    return(0L)
  }
  max(vapply(acc_defs, function(acc) {
    length(dist_param_names(acc$dist %||% ""))
  }, integer(1)), 0L)
}

.outcome_allowed_components <- function(prep) {
  outcome_defs <- prep$outcomes %||% list()
  outcome_labels <- names(outcome_defs)
  component_ids <- prep$components$ids %||% "__default__"
  out <- setNames(vector("list", length(component_ids)), component_ids)
  for (cid in component_ids) {
    out[[cid]] <- character(0)
  }
  for (label in outcome_labels) {
    options <- outcome_defs[[label]]$options %||% list()
    allowed <- options$component %||% component_ids
    if (length(allowed) == 0L) {
      allowed <- component_ids
    }
    for (cid in allowed) {
      out[[cid]] <- c(out[[cid]], label)
    }
  }
  lapply(out, unique)
}

.observed_outcome_allowed_components <- function(prep) {
  outcome_defs <- prep$outcomes %||% list()
  component_ids <- prep$components$ids %||% "__default__"
  out <- setNames(vector("list", length(component_ids)), component_ids)
  for (cid in component_ids) {
    out[[cid]] <- character(0)
  }
  for (label in names(outcome_defs)) {
    options <- outcome_defs[[label]]$options %||% list()
    allowed <- options$component %||% component_ids
    if (length(allowed) == 0L) {
      allowed <- component_ids
    }
    observed_labels <- if (!is.null(options$guess)) {
      as.character(options$guess$labels %||% character(0))
    } else {
      label
    }
    if (!is.null(options$map_outcome_to)) {
      map_to <- options$map_outcome_to
      observed_labels <- if (length(map_to) == 0L || is.na(map_to)) {
        character(0)
      } else {
        as.character(map_to[[1L]])
      }
    }
    observed_labels <- observed_labels[!is.na(observed_labels) & nzchar(observed_labels)]
    for (cid in allowed) {
      out[[cid]] <- c(out[[cid]], observed_labels)
    }
  }
  lapply(out, unique)
}

.has_observation_wrappers <- function(prep) {
  outcome_defs <- prep$outcomes %||% list()
  any(vapply(outcome_defs, function(outcome) {
    options <- outcome$options %||% list()
    !is.null(options$guess) || !is.null(options$map_outcome_to)
  }, logical(1)))
}

.validate_trial_level_columns <- function(data_df, columns) {
  columns <- intersect(columns, names(data_df))
  if (length(columns) == 0L || nrow(data_df) <= 1L) {
    return(invisible(NULL))
  }
  trial <- as.integer(data_df$trials)
  starts <- c(1L, which(trial[-1L] != trial[-nrow(data_df)]) + 1L)
  ends <- c(starts[-1L] - 1L, nrow(data_df))
  if (length(starts) == nrow(data_df)) {
    return(invisible(NULL))
  }
  same_value <- function(x, rows) {
    ref <- x[[rows[[1L]]]]
    vals <- x[rows]
    same_na <- is.na(vals) & is.na(ref)
    same_val <- !is.na(vals) & !is.na(ref) & vals == ref
    all(same_na | same_val)
  }
  for (col in columns) {
    x <- data_df[[col]]
    for (i in seq_along(starts)) {
      rows <- starts[[i]]:ends[[i]]
      if (!same_value(x, rows)) {
        stop(
          sprintf(
            "Prepared data must keep trial-level column '%s' constant within each trial",
            col
          ),
          call. = FALSE
        )
      }
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
    if (time_missing && !allow_missing_rt) {
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

.compress_prepared_trials <- function(data_df) {
  trial <- as.integer(data_df$trials)
  n_rows <- length(trial)
  if (n_rows == 0L) {
    return(data_df)
  }
  trial_starts <- c(1L, which(trial[-1L] != trial[-n_rows]) + 1L)
  n_trials <- length(trial_starts)
  if (n_trials <= 1L) {
    return(data_df)
  }
  trial_ends <- c(trial_starts[-1L] - 1L, n_rows)
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
    trial_ends[keep_trials] - trial_starts[keep_trials] + 1L
  )
  rownames(out) <- NULL
  attr(out, "expand") <- as.integer(match(signatures, signatures[keep_trials]))
  class(out) <- class(data_df)
  out
}

.attach_prepared_layout_attrs <- function(data_df, max_rank) {
  trial <- as.integer(data_df$trials)
  n_rows <- length(trial)
  trial_starts <- c(1L, which(trial[-1L] != trial[-n_rows]) + 1L)

  rank_names <- c("R", if (max_rank > 1L) paste0("R", seq.int(2L, max_rank)))
  time_names <- c("rt", if (max_rank > 1L) paste0("rt", seq.int(2L, max_rank)))

  attr(data_df, "trials_start_rows") <- as.integer(trial_starts)
  attr(data_df, "layout_cols") <- setNames(
    as.integer(match(c("component", "onset"), names(data_df))),
    c("component", "onset")
  )
  attr(data_df, "label_cols") <- as.integer(match(rank_names, names(data_df)))
  attr(data_df, "time_cols") <- as.integer(match(time_names, names(data_df)))
  attr(data_df, "max_rank") <- as.integer(max_rank)
  data_df
}

.prepare_data_structure <- function(structure, data_df, prep = NULL, compress = FALSE) {
  if (inherits(data_df, "accumulatr_data")) {
    if (is.null(attr(data_df, "max_rank", exact = TRUE))) {
      attr(data_df, "max_rank") <- .validate_ranked_observation_columns(data_df)$max_rank
    }
    return(.attach_prepared_layout_attrs(data_df, attr(data_df, "max_rank", exact = TRUE)))
  }
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain R/rt per trial", call. = FALSE)
  }
  prep_info <- .prepare_likelihood_prep(structure, prep = prep)
  structure <- prep_info$structure
  prep_eval_base <- prep_info$prep
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
  rank_info <- .validate_ranked_observation_columns(data_df)
  if (!"racer" %in% names(data_df)) {
    data_df <- .expand_accumulator_rows(structure, data_df)
  }
  if (!"onset" %in% names(data_df)) {
    acc_defs <- prep_eval_base$accumulators %||% list()
    acc_onset <- vapply(acc_defs, function(a) a$onset %||% 0, numeric(1))
    acc_ids <- names(acc_defs)
    onset_map <- setNames(acc_onset, acc_ids)
    data_df$onset <- vapply(as.character(data_df$racer), function(acc) {
      onset_map[[acc]] %||% 0
    }, numeric(1))
  }
  data_df$trials <- {
    trial_in <- as.integer(data_df$trials)
    trial_out <- integer(length(trial_in))
    trial_idx <- 0L
    last_trial <- NA_integer_
    for (i in seq_along(trial_in)) {
      if (i == 1L || trial_in[[i]] != last_trial) {
        trial_idx <- trial_idx + 1L
        last_trial <- trial_in[[i]]
      }
      trial_out[[i]] <- trial_idx
    }
    trial_out
  }
  data_df$rt <- as.numeric(data_df$rt)
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
      data_df[[paste0("rt", rank)]] <- as.numeric(data_df[[paste0("rt", rank)]])
    }
  }
  component_levels <- structure$components$component_id %||% "__default__"
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
  allowed_by_component <- .outcome_allowed_components(prep_eval_base)
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
    data_df <- .compress_prepared_trials(data_df)
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
#'   and so on.
#' @param compress If `TRUE`, collapse repeated prepared trials and attach
#'   an `expand` index so `log_likelihood()` can return trial-level values on
#'   the original trial scale.
#'   Defaults to `FALSE`.
#' @param prep Optional preprocessed model bundle.
#' @return An `accumulatr_data` object.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_param_matrix(
#'   spec,
#'   c(m = 0, s = 0.1),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' prepare_data(structure, data_df)
#' @export
prepare_data <- function(structure,
                         data_df,
                         compress = FALSE,
                         prep = NULL) {
  .prepare_data_structure(
    structure = structure,
    data_df = data_df,
    compress = compress,
    prep = prep
  )
}

.make_context_structure <- function(structure, prep = NULL, diagnostics = FALSE) {
  prep_info <- .prepare_likelihood_prep(structure, prep = prep)
  prep_eval_base <- prep_info$prep
  structure(list(
    cpp = .make_likelihood_context_prep(prep_eval_base, diagnostics = diagnostics),
    required_p_slots = .required_p_slots_from_prep(prep_eval_base)
  ), class = "accumulatr_context")
}

#' Build a compiled likelihood context from a model
#'
#' A context stores compiled model/runtime state only. Behavioral data are
#' prepared separately with `prepare_data()` and supplied to
#' `log_likelihood()`.
#'
#' @param structure Finalized model structure.
#' @param prep Optional preprocessed model bundle.
#' @param diagnostics If `TRUE`, collect symbolic/compiled complexity metrics.
#' @return An `accumulatr_context` object.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' make_context(structure)
#' @export
make_context <- function(structure, prep = NULL, diagnostics = FALSE) {
  .make_context_structure(
    structure = structure,
    prep = prep,
    diagnostics = diagnostics
  )
}

#' Return compiled exact complexity metrics
#'
#' @param context Context created with `make_context(diagnostics = TRUE)`.
#' @return A list with per-variant and total symbolic/compiled metrics.
#' @export
complexity_metrics <- function(context) {
  context <- .validate_context(context)
  if (!isTRUE(context$cpp$has_complexity_metrics)) {
    stop(
      "complexity metrics were not collected; create the context with diagnostics = TRUE",
      call. = FALSE
    )
  }
  .complexity_metrics_context(context$cpp$native)
}

.validate_context <- function(context) {
  if (inherits(context, "accumulatr_context")) {
    return(context)
  }
  stop("context must be created via make_context()", call. = FALSE)
}

.validate_prepared_data <- function(data) {
  if (inherits(data, "accumulatr_data")) {
    required_attrs <- c("trials_start_rows", "layout_cols", "label_cols", "time_cols", "max_rank")
    missing_attrs <- required_attrs[vapply(required_attrs, function(attr_name) {
      is.null(attr(data, attr_name, exact = TRUE))
    }, logical(1))]
    if (length(missing_attrs) > 0L) {
      stop("prepared data are missing native layout metadata; rebuild with prepare_data()", call. = FALSE)
    }
    return(data)
  }
  stop("data must be created via prepare_data()", call. = FALSE)
}

.normalize_prepared_index_column <- function(x, levels, column_name, allow_empty = FALSE) {
  if (length(levels) == 0L) {
    if (allow_empty) {
      return(x)
    }
    stop(sprintf("Prepared data column '%s' has no valid levels in the model", column_name), call. = FALSE)
  }
  if (is.factor(x)) {
    raw_lbl <- as.character(x)
    out <- factor(raw_lbl, levels = levels, ordered = is.ordered(x))
    bad <- !is.na(raw_lbl) & nzchar(raw_lbl) & is.na(out)
    if (any(bad)) {
      stop(
        sprintf(
          "Prepared data encountered unknown labels in '%s': %s",
          column_name,
          paste(utils::head(unique(raw_lbl[bad]), 5L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(out)
  }
  if (is.character(x)) {
    out <- factor(x, levels = levels)
    bad <- !is.na(x) & nzchar(x) & is.na(out)
    if (any(bad)) {
      stop(
        sprintf(
          "Prepared data encountered unknown labels in '%s': %s",
          column_name,
          paste(utils::head(unique(x[bad]), 5L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(out)
  }
  if (is.integer(x)) {
    bad <- !is.na(x) & (x < 1L | x > length(levels))
    if (any(bad)) {
      stop(
        sprintf(
          "Prepared data integer '%s' values must fall in [1, %d]",
          column_name,
          length(levels)
        ),
        call. = FALSE
      )
    }
    raw_lbl <- rep.int(NA_character_, length(x))
    keep <- !is.na(x)
    raw_lbl[keep] <- levels[x[keep]]
    return(factor(raw_lbl, levels = levels))
  }
  stop(
    sprintf(
      "Prepared data column '%s' must be a factor, character, or integer-coded vector",
      column_name
    ),
    call. = FALSE
  )
}

.param_matrix_to_rows <- function(structure, params_mat) {
  prep <- structure$prep %||% stop("structure missing prep")
  acc_defs <- prep$accumulators %||% list()
  acc_ids <- names(acc_defs)
  n_acc <- length(acc_ids)
  if (is.null(colnames(params_mat))) {
    stop("Parameter matrix must have column names")
  }
  if (nrow(params_mat) %% n_acc != 0) {
    stop(sprintf("Parameter rows (%d) not divisible by accumulators (%d)", nrow(params_mat), n_acc))
  }
  n_trials <- nrow(params_mat) / n_acc
  p_cols <- grep("^p[0-9]+$", colnames(params_mat), value = TRUE)
  if (length(p_cols) == 0L) stop("Parameter matrix must include p1.. columns")
  p_cols <- p_cols[order(suppressWarnings(as.integer(sub("^p", "", p_cols))))]
  p_cols <- p_cols[seq_len(min(length(p_cols), 8L))]
  dist_param_list <- lapply(acc_defs, function(acc) dist_param_names(acc$dist))

  # Component leaders for weight params
  comp_ids <- structure$components$component_id %||% character(0)
  comp_leader_idx <- setNames(rep(NA_integer_, length(comp_ids)), comp_ids)
  for (i in seq_along(acc_ids)) {
    comps <- acc_defs[[acc_ids[[i]]]]$components %||% character(0)
    for (c_id in comps) {
      if (is.na(comp_leader_idx[[c_id]])) comp_leader_idx[[c_id]] <- i
    }
  }
  comp_attrs <- structure$components$attrs %||% vector("list", length(comp_ids))
  comp_weight_param <- setNames(rep(NA_character_, length(comp_ids)), comp_ids)
  for (i in seq_along(comp_ids)) {
    attrs <- comp_attrs[[i]] %||% list()
    if (!is.null(attrs$weight_param) && nzchar(attrs$weight_param)) {
      comp_weight_param[[i]] <- attrs$weight_param
    }
  }
  # All column names we may emit (dist params can differ by accumulator)
  global_dist_params <- unique(unlist(dist_param_list))
  weight_param_names <- unique(comp_weight_param[!is.na(comp_weight_param) & nzchar(comp_weight_param)])
  all_cols <- c(
    "trials", "racer", "q", "t0", "onset",
    "shared_trigger_id", "shared_trigger_q",
    global_dist_params,
    weight_param_names
  )

  rows <- vector("list", nrow(params_mat))
  row_idx <- 1L
  for (t in seq_len(n_trials)) {
    weight_vals <- list()
    for (i in seq_along(comp_ids)) {
      wp <- comp_weight_param[[i]]
      if (is.na(wp) || !nzchar(wp)) next
      leader_idx <- comp_leader_idx[[comp_ids[[i]]]]
      if (is.na(leader_idx)) next
      row_num <- (t - 1L) * n_acc + leader_idx
      weight_vals[[wp]] <- params_mat[row_num, "w"]
    }
    for (i in seq_len(n_acc)) {
      acc <- acc_defs[[acc_ids[[i]]]] %||% list()
      dist_params <- dist_param_list[[i]] %||% character(0)
      if (length(dist_params) > length(p_cols)) {
        stop(
          sprintf(
            "Accumulator '%s' requires %d distribution params, but matrix has only %d p-slots",
            acc_ids[[i]],
            length(dist_params),
            length(p_cols)
          ),
          call. = FALSE
        )
      }
      p_vals <- params_mat[row_idx, p_cols][seq_along(dist_params)]
      row <- list(
        trials = t,
        racer = i,
        q = params_mat[row_idx, "q"],
        t0 = params_mat[row_idx, "t0"],
        onset = acc$onset %||% 0,
        shared_trigger_id = acc$shared_trigger_id %||% NA_character_,
        shared_trigger_q = params_mat[row_idx, "q"]
      )
      for (k in seq_along(dist_params)) {
        row[[dist_params[[k]]]] <- as.numeric(p_vals[[k]])
      }
      # pad any missing columns with NA and order consistently
      missing_cols <- setdiff(all_cols, names(row))
      if (length(missing_cols)) {
        for (nm in missing_cols) row[[nm]] <- NA
      }
      row <- row[all_cols]
      if (length(weight_vals) > 0) {
        for (nm in names(weight_vals)) {
          row[[nm]] <- weight_vals[[nm]]
        }
      }
      rows[[row_idx]] <- row
      row_idx <- row_idx + 1L
    }
  }
  df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(df) <- NULL
  df
}

.params_df_to_matrix <- function(df) {
  if (!"trials" %in% names(df)) {
    stop("Parameter table must include 'trials'", call. = FALSE)
  }

  n <- nrow(df)
  if (is.null(n) || n == 0L) {
    stop("Parameter table must contain rows", call. = FALSE)
  }

  # Core columns
  trial <- as.numeric(df$trials)
  acc_idx <- if ("racer_index" %in% names(df)) {
    as.numeric(df$racer_index)
  } else if ("racer" %in% names(df)) {
    as.numeric(df$racer)
  } else {
    stop("Parameter table must include 'racer_index' (or numeric 'racer')", call. = FALSE)
  }
  onset <- as.numeric(df$onset %||% 0)
  q <- as.numeric(df$q %||% 0)
  t0 <- as.numeric(df$t0 %||% 0)
  shared_q <- if ("shared_trigger_q" %in% names(df)) as.numeric(df$shared_trigger_q) else rep(NA_real_, n)
  p_mat <- matrix(NA_real_, nrow = n, ncol = 8L)

  # Prefer explicit p-slot columns first.
  for (slot in seq_len(8L)) {
    nm <- paste0("p", slot)
    if (!nm %in% names(df)) next
    p_mat[, slot] <- as.numeric(df[[nm]])
  }

  # If dist is available, map named params to canonical p-slot order per row.
  if ("dist" %in% names(df)) {
    dvals <- tolower(as.character(df$dist))
    dvals[is.na(dvals)] <- ""
    for (i in seq_len(n)) {
      d <- dvals[[i]]
      if (!nzchar(d)) next
      d_params <- dist_param_names(d)
      if (length(d_params) == 0L) next
      d_params <- d_params[seq_len(min(length(d_params), 8L))]
      for (slot in seq_along(d_params)) {
        if (!is.na(p_mat[i, slot])) next
        nm <- d_params[[slot]]
        if (!nm %in% names(df)) next
        val <- suppressWarnings(as.numeric(df[[nm]][[i]]))
        if (!is.na(val)) {
          p_mat[i, slot] <- val
        }
      }
    }
  }

  fill_slot <- function(slot, candidates) {
    if (slot < 1L || slot > ncol(p_mat)) return(invisible(NULL))
    for (nm in candidates) {
      if (!nm %in% names(df)) next
      col <- as.numeric(df[[nm]])
      mask <- is.na(p_mat[, slot]) & !is.na(col)
      if (any(mask)) {
        p_mat[mask, slot] <<- col[mask]
      }
    }
    invisible(NULL)
  }
  fill_slot(1L, c("m", "shape", "mu", "v"))
  fill_slot(2L, c("s", "rate", "sigma", "B"))
  fill_slot(3L, c("tau", "A"))
  fill_slot(4L, c("sv", "s"))

  has_slot <- vapply(seq_len(8L), function(slot) any(!is.na(p_mat[, slot])), logical(1))
  max_params <- if (any(has_slot)) max(which(has_slot)) else 0L
  max_params <- as.integer(min(max_params, 8L))

  cols <- list(
    trials = trial,
    racer_index = acc_idx,
    onset = onset,
    q = q,
    t0 = t0,
    shared_trigger_q = shared_q
  )
  if (max_params >= 1L) {
    for (slot in seq_len(max_params)) {
      cols[[paste0("p", slot)]] <- p_mat[, slot]
    }
  }

  base_cols <- c(
    "trials",
    "racer_index",
    "racer",
    "onset",
    "q",
    "t0",
    "shared_trigger_q",
    paste0("p", seq_len(8L))
  )
  extra_numeric <- setdiff(names(df), base_cols)
  for (nm in extra_numeric) {
    col <- df[[nm]]
    if (is.numeric(col) || is.integer(col)) {
      cols[[nm]] <- as.numeric(col)
    }
  }

  mat <- do.call(cbind, cols)
  storage.mode(mat) <- "double"
  attr(mat, "trials") <- trial
  colnames(mat) <- names(cols)
  mat
}

.likelihood_params_by_trial <- function(params_df) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    return(list())
  }
  trials <- params_df$trials %||% seq_len(nrow(params_df))
  if (is.factor(trials)) trials <- as.character(trials)
  split(params_df, as.character(trials))
}

.coerce_loglik_param_matrix <- function(parameters) {
  if (inherits(parameters, "param_matrix") || is.matrix(parameters)) {
    return(parameters)
  }
  if (is.data.frame(parameters)) {
    return(.params_df_to_matrix(parameters))
  }
  stop("parameters must be a param_matrix, matrix, or data frame", call. = FALSE)
}

.canonicalize_loglik_param_matrix <- function(param_mat, required_p_slots) {
  cn <- colnames(param_mat)
  if (is.null(cn) || anyNA(cn) || any(!nzchar(cn))) {
    stop("Parameter matrix must have named columns", call. = FALSE)
  }
  required_cols <- c("q", "t0", if (required_p_slots > 0L) paste0("p", seq_len(required_p_slots)))
  missing_cols <- setdiff(required_cols, cn)
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Parameter matrix is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  canonical_prefix <- identical(cn[seq_along(required_cols)], required_cols)
  if (canonical_prefix && is.double(param_mat) && is.matrix(param_mat)) {
    return(param_mat)
  }
  extra_cols <- setdiff(cn, required_cols)
  out <- unclass(as.matrix(param_mat[, c(required_cols, extra_cols), drop = FALSE]))
  storage.mode(out) <- "double"
  colnames(out) <- c(required_cols, extra_cols)
  out
}

.coerce_response_probability_params <- function(structure, params_df) {
  if (inherits(params_df, "param_matrix") || is.matrix(params_df)) {
    n_acc <- length(structure$prep$accumulators %||% list())
    if (n_acc <= 0L || nrow(params_df) %% n_acc != 0L) {
      stop(
        "response_probabilities() expects a rectangular param_matrix; use a row-form parameter data frame for prepared-row-aligned inputs",
        call. = FALSE
      )
    }
    params_df <- .param_matrix_to_rows(structure, params_df)
  } else if (is.data.frame(params_df)) {
    params_df <- as.data.frame(params_df, stringsAsFactors = FALSE)
  } else {
    stop("params_df must be a param_matrix, matrix, or data frame", call. = FALSE)
  }

  if (!"trials" %in% names(params_df)) {
    params_df$trials <- 1L
  }
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }
  params_df
}

.response_component_accumulator_indices <- function(structure) {
  acc_ids <- names(structure$prep$accumulators %||% list())
  comp_ids <- structure$components$component_id %||% "__default__"
  out <- setNames(vector("list", length(comp_ids)), comp_ids)
  if (identical(comp_ids, "__default__")) {
    out[["__default__"]] <- seq_along(acc_ids)
    return(out)
  }
  for (comp_id in comp_ids) {
    active <- vapply(structure$prep$accumulators %||% list(), function(acc) {
      comp_id %in% (acc$components %||% character(0))
    }, logical(1))
    out[[comp_id]] <- which(active)
  }
  out
}

.response_param_accumulator_indices <- function(structure, rows_df) {
  acc_ids <- names(structure$prep$accumulators %||% list())
  if ("racer_index" %in% names(rows_df)) {
    idx <- as.integer(rows_df$racer_index)
  } else if ("racer" %in% names(rows_df)) {
    acc <- rows_df$racer
    if (is.factor(acc)) {
      idx <- match(as.character(acc), acc_ids)
    } else if (is.numeric(acc) || is.integer(acc)) {
      idx <- as.integer(acc)
    } else {
      idx <- match(as.character(acc), acc_ids)
    }
  } else {
    stop(
      "response_probabilities() parameter rows must include racer information",
      call. = FALSE
    )
  }
  if (anyNA(idx) || any(idx < 1L) || any(idx > length(acc_ids))) {
    stop(
      "response_probabilities() parameter rows use invalid racer indices",
      call. = FALSE
    )
  }
  idx
}

.align_response_probability_param_rows <- function(structure,
                                                   source_rows,
                                                   active_accumulator_indices,
                                                   trial_index,
                                                   component = NULL) {
  source_rows <- as.data.frame(source_rows, stringsAsFactors = FALSE)
  if (!is.null(component) && "component" %in% names(source_rows)) {
    comp_chr <- as.character(source_rows$component)
    source_rows <- source_rows[is.na(comp_chr) | comp_chr == component, , drop = FALSE]
  }
  source_acc <- .response_param_accumulator_indices(structure, source_rows)
  match_idx <- match(active_accumulator_indices, source_acc)
  if (anyNA(match_idx)) {
    stop(
      sprintf(
        "response_probabilities() could not align parameter rows for component '%s'",
        component %||% "__default__"
      ),
      call. = FALSE
    )
  }
  out <- source_rows[match_idx, , drop = FALSE]
  out$trials <- trial_index
  out$racer_index <- active_accumulator_indices
  out$racer <- active_accumulator_indices
  out
}

.response_probability_param_layout <- function(structure, params_rows, required_p_slots) {
  params_by_trial <- .likelihood_params_by_trial(params_rows)
  if (length(params_by_trial) == 0L) {
    stop("Parameter data frame must contain at least one row", call. = FALSE)
  }
  acc_ids <- names(structure$prep$accumulators %||% list())
  component_ids <- structure$components$component_id %||% "__default__"
  component_accumulators <- .response_component_accumulator_indices(structure)

  full_start_rows <- rep.int(NA_integer_, length(params_by_trial))
  component_start_rows <- matrix(
    NA_integer_,
    nrow = length(params_by_trial),
    ncol = length(component_ids),
    dimnames = list(NULL, component_ids)
  )
  aligned_rows <- list()
  next_row <- 1L

  append_block <- function(block, trial_index, component_id = NULL) {
    if (nrow(block) == 0L) {
      return(invisible(NULL))
    }
    aligned_rows[[length(aligned_rows) + 1L]] <<- block
    if (is.null(component_id)) {
      full_start_rows[[trial_index]] <<- next_row
    } else {
      component_start_rows[trial_index, component_id] <<- next_row
    }
    next_row <<- next_row + nrow(block)
    invisible(NULL)
  }

  for (trial_index in seq_along(params_by_trial)) {
    source_rows <- as.data.frame(params_by_trial[[trial_index]], stringsAsFactors = FALSE)
    component_values <- if ("component" %in% names(source_rows)) {
      as.character(source_rows$component)
    } else {
      rep(NA_character_, nrow(source_rows))
    }
    listed_components <- unique(component_values[!is.na(component_values) & nzchar(component_values)])
    if (length(listed_components) == 0L) {
      block <- .align_response_probability_param_rows(
        structure,
        source_rows,
        seq_along(acc_ids),
        trial_index
      )
      append_block(block, trial_index)
      next
    }

    unknown_components <- setdiff(listed_components, component_ids)
    if (length(unknown_components) > 0L) {
      stop(
        sprintf(
          "response_probabilities() parameter rows reference unknown component(s): %s",
          paste(unknown_components, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    for (component_id in component_ids[component_ids %in% listed_components]) {
      active_accumulators <- component_accumulators[[component_id]] %||% integer(0)
      if (length(active_accumulators) == 0L) {
        next
      }
      block <- .align_response_probability_param_rows(
        structure,
        source_rows,
        active_accumulators,
        trial_index,
        component_id
      )
      append_block(block, trial_index, component_id)
    }
  }

  if (length(aligned_rows) == 0L) {
    stop("response_probabilities() could not build parameter rows", call. = FALSE)
  }
  param_mat <- .params_df_to_matrix(do.call(rbind, aligned_rows))
  param_mat <- .canonicalize_loglik_param_matrix(param_mat, required_p_slots)
  list(
    params = param_mat,
    layout = list(
      full_start_rows = as.integer(full_start_rows),
      component_start_rows = component_start_rows
    )
  )
}

#' Evaluate marginal response probabilities
#'
#' `response_probabilities()` evaluates the model-implied marginal probability
#' of each observed response label for a finalized model and parameter set.
#' Component labels in row-form parameters condition the calculation on those
#' observed components; latent or `NA` components are marginalized according to
#' the model's mixture specification.
#'
#' @param structure Finalized model structure.
#' @param params_df A parameter data frame or rectangular parameter matrix.
#' @param include_na If `TRUE`, include residual mass as `"NA"`.
#' @param ... Unused; for S3 compatibility.
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
#'   spec,
#'   c(left.m = log(0.25), right.m = log(0.40), s = 0.20),
#'   n_trials = 1
#' )
#'
#' response_probabilities(model, params)
#' @export
response_probabilities <- function(structure, params_df, include_na = TRUE, ...) {
  UseMethod("response_probabilities")
}

#' @rdname response_probabilities
#' @export
response_probabilities.model_structure <- function(structure,
                                                   params_df,
                                                   include_na = TRUE,
                                                   ...) {
  if (is.null(params_df) || NROW(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row", call. = FALSE)
  }

  prep_info <- .prepare_likelihood_prep(structure)
  structure <- prep_info$structure
  params_rows <- .coerce_response_probability_params(structure, params_df)
  context <- .make_context_structure(structure, prep = prep_info$prep)
  query <- .response_probability_param_layout(
    structure,
    params_rows,
    context$required_p_slots
  )
  labels <- names(prep_info$prep$outcomes %||% list())
  probability <- as.numeric(.response_probabilities_context(
    context$cpp$native,
    query$params,
    query$layout
  ))
  if (length(probability) < length(labels)) {
    probability <- c(probability, numeric(length(labels) - length(probability)))
  } else if (length(probability) > length(labels)) {
    probability <- probability[seq_along(labels)]
  }
  result <- stats::setNames(probability, labels)
  observed_by_component <- .observed_outcome_allowed_components(prep_info$prep)
  observed_labels <- Reduce(union, observed_by_component, init = character(0))
  result <- result[intersect(names(result), observed_labels)]
  residual <- 1.0 - sum(result)
  if (!is.finite(residual) || residual < 0) {
    residual <- 0.0
  }
  if (include_na && residual > .Machine$double.eps) {
    result["NA"] <- residual
  }
  if (!include_na) {
    result <- result[setdiff(names(result), "NA")]
  }
  result
}

#' @rdname response_probabilities
#' @export
response_probabilities.default <- function(structure,
                                           params_df,
                                           include_na = TRUE,
                                           ...) {
  stop(
    "response_probabilities() expects a finalized model structure",
    call. = FALSE
  )
}

#' Evaluate log-likelihoods of behavioral data
#'
#' Compute the summed log-likelihood by default, or trial-wise log-likelihoods
#' when `sum = FALSE`.
#'
#' @param context Context created with `make_context()`.
#' @param data Prepared data created with `prepare_data()`.
#' @param parameters A parameter data frame, or a list of parameter data frames.
#' @param ok Logical vector marking which trials should contribute to the
#'   likelihood. Trials marked `FALSE` are assigned `min_ll`.
#' @param sum If `TRUE`, return the summed log-likelihood. If `FALSE`, return
#'   trial-wise log-likelihood values.
#' @param min_ll Minimum log-likelihood value used for excluded or impossible
#'   trials.
#' @param ... Unused; for S3 compatibility.
#' @return A summed log-likelihood by default, or a numeric vector of
#'   trial-wise log-likelihood values when `sum = FALSE`.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_param_matrix(
#'   spec,
#'   c(m = 0, s = 0.1),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' prepared <- prepare_data(structure, data_df)
#' ctx <- make_context(structure)
#' log_likelihood(ctx, prepared, params_df)
#' @export
log_likelihood <- function(context, data, parameters, ok = NULL, sum = TRUE, min_ll = log(1e-10), ...) {
  UseMethod("log_likelihood")
}

#' @rdname log_likelihood
#' @export
log_likelihood.accumulatr_context <- function(context,
                                              data,
                                              parameters,
                                              ok = NULL,
                                              sum = TRUE,
                                              min_ll = log(1e-10),
                                              ...) {
  if (is.null(ok) || length(ok) == 0L) {
    ok <- NULL
  }

  cpp_ctx <- context$cpp
  value <- .loglik_context(
    cpp_ctx$native,
    parameters,
    data,
    ok = ok,
    min_ll = min_ll
  )
  if (isTRUE(sum)) {
    base::sum(value)
  } else {
    value
  }
}

#' @rdname log_likelihood
#' @export
log_likelihood.default <- function(context, ...) {
  stop("log_likelihood() expects a context created with make_context()", call. = FALSE)
}
