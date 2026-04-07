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

.prepare_likelihood_prep <- function(structure, prep = NULL, native_bundle = NULL) {
  structure <- .as_model_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("model structure must include model_spec; rebuild with finalize_model")
  }
  if (!is.null(native_bundle)) {
    prep <- native_bundle$prep %||% prep
  }
  prep_eval_base <- prep %||% structure$prep %||% NULL
  if (is.null(prep_eval_base)) {
    prep_eval_base <- prepare_model(structure$model_spec)
  }
  if (!is.null(prep_eval_base[[".runtime"]]) && !is.null(prep_eval_base[[".id_index"]])) {
    return(list(structure = structure, prep = prep_eval_base))
  }
  outcome_ids <- names(prep_eval_base[["outcomes"]] %||% list())
  acc_ids <- names(prep_eval_base[["accumulators"]] %||% list())
  pool_ids <- names(prep_eval_base[["pools"]] %||% list())
  all_ids <- unique(c(outcome_ids, acc_ids, pool_ids))
  prep_eval_base[[".id_index"]] <- setNames(seq_along(all_ids), all_ids)
  prep_eval_base[[".label_cache"]] <- new.env(parent = emptyenv(), hash = TRUE)
  prep_eval_base <- .precompile_likelihood_expressions(prep_eval_base)
  prep_eval_base[[".competitors"]] <- .prepare_competitor_map(prep_eval_base)
  prep_eval_base[[".runtime"]] <- list(
    expr_compiled = prep_eval_base[[".expr_compiled"]],
    label_cache = prep_eval_base[[".label_cache"]],
    competitor_map = prep_eval_base[[".competitors"]],
    id_index = prep_eval_base[[".id_index"]],
    pool_members_cache = new.env(parent = emptyenv(), hash = TRUE),
    cache_bundle = .build_likelihood_cache_bundle(prep_eval_base)
  )
  .refresh_compiled_prep_refs(prep_eval_base)
  list(structure = structure, prep = prep_eval_base)
}

.prepare_data_structure <- function(structure, data_df, prep = NULL, native_bundle = NULL) {
  if (inherits(data_df, "accumulatr_data")) {
    return(data_df)
  }
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain R/rt per trial", call. = FALSE)
  }
  prep_info <- .prepare_likelihood_prep(structure, prep = prep, native_bundle = native_bundle)
  structure <- prep_info$structure
  prep_eval_base <- prep_info$prep
  data_df <- as.data.frame(data_df)
  required_cols <- c("R", "rt")
  missing_cols <- setdiff(required_cols, names(data_df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Data frame must include columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!"trial" %in% names(data_df)) {
    if ("trials" %in% names(data_df)) {
      data_df$trial <- as.integer(data_df$trials)
    } else {
      data_df$trial <- seq_len(nrow(data_df))
    }
  }
  rank_info <- .validate_ranked_observation_columns(data_df)
  if (!"accumulator" %in% names(data_df)) {
    data_df <- .expand_accumulator_rows(structure, data_df)
  }
  if (!"onset" %in% names(data_df)) {
    acc_defs <- prep_eval_base$accumulators %||% list()
    acc_onset <- vapply(acc_defs, function(a) a$onset %||% 0, numeric(1))
    acc_ids <- names(acc_defs)
    onset_map <- setNames(acc_onset, acc_ids)
    data_df$onset <- vapply(as.character(data_df$accumulator), function(acc) {
      onset_map[[acc]] %||% 0
    }, numeric(1))
  }
  data_df$trial <- {
    trial_in <- as.integer(data_df$trial)
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
  if ("component" %in% names(data_df)) {
    component_levels <- structure$components$component_id %||% character(0)
    data_df$component <- .normalize_prepared_index_column(
      data_df$component,
      component_levels,
      "component"
    )
  }
  class(data_df) <- unique(c("accumulatr_data", class(data_df)))
  data_df
}

#' Prepare behavioral data for likelihood evaluation
#'
#' `prepare_data()` expands trial-level observations to the accumulator layout
#' expected by the compiled likelihood code and tags the result as trusted
#' likelihood input.
#'
#' @param structure Finalized model structure.
#' @param data_df Behavioral data. In the simplest case this contains `trial`,
#'   `R`, and `rt`; for multi-outcome models it can also contain `R2`, `rt2`,
#'   and so on.
#' @param prep Optional preprocessed model bundle.
#' @param native_bundle Optional serialized native bundle.
#' @return An `accumulatr_data` object.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_param_matrix(
#'   spec,
#'   c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' prepare_data(structure, data_df)
#' @export
prepare_data <- function(structure,
                         data_df,
                         prep = NULL,
                         native_bundle = NULL) {
  .prepare_data_structure(
    structure = structure,
    data_df = data_df,
    prep = prep,
    native_bundle = native_bundle
  )
}

.make_context_structure <- function(structure, prep = NULL, native_bundle = NULL) {
  prep_info <- .prepare_likelihood_prep(structure, prep = prep, native_bundle = native_bundle)
  prep_eval_base <- prep_info$prep
  native_ctx <- .prep_native_context(prep_eval_base)
  if (!inherits(native_ctx, "externalptr")) {
    native_ctx <- NULL
  }
  structure(list(
    native_ctx = native_ctx
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
#' @param native_bundle Optional serialized native bundle.
#' @return An `accumulatr_context` object.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' make_context(structure)
#' @export
make_context <- function(structure, prep = NULL, native_bundle = NULL) {
  .make_context_structure(
    structure = structure,
    prep = prep,
    native_bundle = native_bundle
  )
}

.validate_context <- function(context) {
  if (inherits(context, "accumulatr_context")) {
    return(context)
  }
  stop("context must be created via make_context()", call. = FALSE)
}

.validate_prepared_data <- function(data) {
  if (inherits(data, "accumulatr_data")) {
    return(data)
  }
  stop("data must be created via prepare_data()", call. = FALSE)
}

.component_index_from_prep <- function(prep, component) {
  comp_label <- as.character(component %||% "__default__")
  comp_label <- if (length(comp_label) > 0L) comp_label[[1L]] else "__default__"
  if (is.na(comp_label) || !nzchar(comp_label) || identical(comp_label, "__default__")) {
    return(-1L)
  }
  comp_ids <- as.character((prep[["components"]] %||% list())[["ids"]] %||% "__default__")
  idx <- match(comp_label, comp_ids)
  if (is.na(idx)) {
    stop(sprintf("Unknown component '%s' for IR outcome probability", comp_label), call. = FALSE)
  }
  as.integer(idx - 1L)
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
    return(x)
  }
  stop(
    sprintf(
      "Prepared data column '%s' must be a factor, character, or integer vector",
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
    "trial", "accumulator", "q", "t0", "onset",
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
        trial = t,
        accumulator = i,
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

.likelihood_component_rows <- function(trial_rows, component) {
  if (is.null(trial_rows) || nrow(trial_rows) == 0L) {
    return(trial_rows)
  }
  if (!"component" %in% names(trial_rows)) {
    return(trial_rows)
  }
  comp_label <- component %||% "__default__"
  mask <- is.na(trial_rows$component)
  if (!identical(comp_label, "__default__")) {
    mask <- mask | trial_rows$component == comp_label
  }
  trial_rows[mask, , drop = FALSE]
}

.model_spec_with_params <- function(model_spec, params_df) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    return(model_spec)
  }
  spec_copy <- unserialize(serialize(model_spec, NULL))
  prep_tmp <- .prepare_model_for_likelihood(model_spec)
  param_state <- .generator_param_state_from_rows(prep_tmp, params_df)
  acc_overrides <- param_state$accs %||% list()
  if (length(acc_overrides) > 0L) {
    acc_defs <- spec_copy$accumulators %||% list()
    acc_index <- setNames(seq_along(acc_defs), vapply(acc_defs, `[[`, character(1), "id"))
    for (acc_id in names(acc_overrides)) {
      idx <- acc_index[[acc_id]]
      if (is.na(idx)) next
      override <- acc_overrides[[acc_id]]
      acc <- acc_defs[[idx]]
      acc$params <- override$params
      acc$onset <- override$onset
      acc$q <- override$q
      acc_defs[[idx]] <- acc
    }
    spec_copy$accumulators <- acc_defs
  }
  shared_overrides <- param_state$shared %||% list()
  if (length(shared_overrides) > 0L && length(spec_copy$groups) > 0L) {
    for (g in seq_along(spec_copy$groups)) {
      grp <- spec_copy$groups[[g]] %||% list()
      trig <- grp$attrs$shared_trigger %||% NULL
      if (is.null(trig)) next
      trig_id <- trig$id %||% grp$id
      if (is.null(trig_id) || !nzchar(trig_id)) next
      override <- shared_overrides[[trig_id]] %||% NULL
      if (is.null(override) || is.null(override$prob)) next
      grp$attrs$shared_trigger$q <- override$prob
      spec_copy$groups[[g]] <- grp
    }
  }
  spec_copy
}

.likelihood_build_trial_plan <- function(structure, params_df, prep) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    return(list(order = character(0), trials = list()))
  }
  params_df <- as.data.frame(params_df)
  if (is.null(rownames(params_df))) {
    rownames(params_df) <- as.character(seq_len(nrow(params_df)))
  }
  if (!"trial" %in% names(params_df)) params_df$trial <- 1L
  params_df$trial <- params_df$trial
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }
  trial_rows_map <- .likelihood_params_by_trial(params_df)
  component_ids <- structure$components$component_id
  comp_labels <- component_ids
  if (length(comp_labels) == 0L) comp_labels <- "__default__"
  comp_labels <- unique(c(comp_labels, "__default__"))

  trials <- list()
  trial_order <- character(0)
  for (trial_key in names(trial_rows_map)) {
    trial_df <- as.data.frame(trial_rows_map[[trial_key]])
    trial_id <- trial_df$trial[[1]] %||% NA
    comp_entries <- list()
    for (comp_id in comp_labels) {
      rows_df <- .likelihood_component_rows(trial_df, comp_id)
      comp_entries[[comp_id]] <- list(rows = rows_df)
    }
    trials[[trial_key]] <- list(
      trial_id = trial_id,
      rows = trial_df,
      components = comp_entries
    )
    trial_order <- c(trial_order, trial_key)
  }
  order_idx <- order(suppressWarnings(as.numeric(trial_order)))
  trial_order <- trial_order[order_idx]
  trials <- trials[trial_order]
  list(
    order = trial_order,
    trials = trials,
    .native_cache = new.env(parent = emptyenv())
  )
}

.likelihood_trial_metadata <- function(structure, prep) {
  structure_hash <- .structure_hash_value(prep)
  outcome_defs <- prep[["outcomes"]] %||% list()
  competitor_map <- .prep_competitors(prep) %||% list()
  outcome_labels <- names(outcome_defs)
  use_indexed_competitors <- isTRUE(attr(competitor_map, "by_index")) &&
    length(competitor_map) == length(outcome_defs)
  compiled_nodes <- lapply(outcome_defs, function(def) {
    expr <- def[["expr"]]
    if (is.null(expr)) {
      return(NULL)
    }
    .expr_lookup_compiled(expr, prep)
  })
  names(compiled_nodes) <- outcome_labels
  compiled_competitors <- lapply(seq_along(outcome_defs), function(idx) {
    comp_exprs <- if (use_indexed_competitors) {
      competitor_map[[idx]] %||% list()
    } else {
      competitor_map[[outcome_labels[[idx]]]] %||% list()
    }
    if (length(comp_exprs) == 0) {
      return(integer(0))
    }
    comp_nodes <- lapply(comp_exprs, function(ex) .expr_lookup_compiled(ex, prep))
    if (any(vapply(comp_nodes, is.null, logical(1)))) {
      return(NULL)
    }
    ids <- vapply(comp_nodes, function(node) as.integer(node$id %||% NA_integer_), integer(1))
    if (any(is.na(ids))) {
      return(NULL)
    }
    ids
  })
  names(compiled_competitors) <- outcome_labels
  na_source_idx <- Filter(function(idx) {
    def <- outcome_defs[[idx]]
    map_to <- def[["options"]][["map_outcome_to"]]
    if (is.null(map_to)) {
      return(FALSE)
    }
    if (is.na(map_to)) {
      return(TRUE)
    }
    identical(map_to, "NA")
  }, seq_along(outcome_defs))
  na_source_specs <- NULL
  if (length(na_source_idx) > 0L) {
    na_source_specs <- lapply(na_source_idx, function(idx) {
      node <- compiled_nodes[[idx]]
      def <- outcome_defs[[idx]] %||% list()
      comp_ids <- compiled_competitors[[idx]]
      if (is.null(node) || is.null(comp_ids)) {
        return(NULL)
      }
      list(
        node_id = as.integer(node$id),
        competitor_ids = comp_ids %||% integer(0),
        source_label = outcome_labels[[idx]]
      )
    })
    if (any(vapply(na_source_specs, is.null, logical(1)))) na_source_specs <- NULL
  }
  guess_target_specs <- list()
  for (idx in seq_along(outcome_defs)) {
    def <- outcome_defs[[idx]]
    opts <- def[["options"]] %||% list()
    guess_opts <- opts[["guess"]]
    if (is.null(guess_opts)) next
    donor_node <- compiled_nodes[[idx]]
    if (is.null(donor_node)) next
    donor_comp <- compiled_competitors[[idx]] %||% integer(0)
    labels <- guess_opts[["labels"]] %||% character(0)
    weights <- guess_opts[["weights"]] %||% numeric(0)
    if (!length(labels) || length(labels) != length(weights)) next
    rt_policy <- guess_opts[["rt_policy"]] %||% "keep"
    for (j in seq_along(labels)) {
      tgt <- labels[[j]]
      tgt_key <- if (is.na(tgt)) "NA" else as.character(tgt)
      release_val <- as.numeric(weights[[j]])
      release_val <- if (is.finite(release_val)) release_val else 0
      if (release_val < 0) release_val <- 0
      if (release_val > 1) release_val <- 1
      donor_rec <- list(
        node_id = as.integer(donor_node$id),
        competitor_ids = donor_comp %||% integer(0),
        source_label = outcome_labels[[idx]],
        release = release_val,
        rt_policy = rt_policy
      )
      guess_target_specs[[tgt_key]] <- c(guess_target_specs[[tgt_key]] %||% list(), list(donor_rec))
    }
  }
  alias_specs <- lapply(seq_along(outcome_defs), function(idx) {
    def <- outcome_defs[[idx]]
    alias_refs <- def[["options"]][["alias_of"]]
    if (is.null(alias_refs)) {
      return(NULL)
    }
    refs <- as.character(alias_refs)
    alias_sources <- lapply(refs, function(ref_lbl) {
      ref_idx <- which(outcome_labels == ref_lbl)
      if (length(ref_idx) == 0) {
        return(NULL)
      }
      refs_out <- lapply(ref_idx, function(j) {
        node <- compiled_nodes[[j]]
        comp_ids <- compiled_competitors[[j]]
        if (is.null(node) || is.null(comp_ids)) {
          return(NULL)
        }
        list(
          node_id = as.integer(node$id),
          competitor_ids = comp_ids %||% integer(0),
          source_label = ref_lbl
        )
      })
      refs_out
    })
    alias_sources <- unlist(alias_sources, recursive = FALSE)
    if (any(vapply(alias_sources, is.null, logical(1)))) {
      return(NULL)
    }
    alias_sources
  })
  names(alias_specs) <- outcome_labels
  rel_tol <- .integrate_rel_tol()
  abs_tol <- .integrate_abs_tol()
  max_depth <- 12L
  list(
    structure_hash = structure_hash,
    compiled_nodes = compiled_nodes,
    compiled_competitors = compiled_competitors,
    na_source_specs = na_source_specs,
    guess_target_specs = guess_target_specs,
    alias_specs = alias_specs,
    rel_tol = as.numeric(rel_tol),
    abs_tol = as.numeric(abs_tol),
    max_depth = as.integer(max_depth)
  )
}

.params_df_to_matrix <- function(df) {
  if (!"trial" %in% names(df)) {
    stop("Parameter table must include 'trial'", call. = FALSE)
  }

  n <- nrow(df)
  if (is.null(n) || n == 0L) {
    stop("Parameter table must contain rows", call. = FALSE)
  }

  # Core columns
  trial <- as.numeric(df$trial)
  acc_idx <- if ("accumulator_index" %in% names(df)) {
    as.numeric(df$accumulator_index)
  } else if ("accumulator" %in% names(df)) {
    as.numeric(df$accumulator)
  } else {
    stop("Parameter table must include 'accumulator_index' (or numeric 'accumulator')", call. = FALSE)
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
    trial = trial,
    accumulator_index = acc_idx,
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

  mat <- do.call(cbind, cols)
  storage.mode(mat) <- "double"
  attr(mat, "trial") <- trial
  colnames(mat) <- names(cols)
  mat
}

.likelihood_params_by_trial <- function(params_df) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    return(list())
  }
  trials <- params_df$trial %||% seq_len(nrow(params_df))
  if (is.factor(trials)) trials <- as.character(trials)
  split(params_df, as.character(trials))
}

.likelihood_component_weights <- function(structure, trial_id, available_components,
                                          base_weights,
                                          trial_rows = NULL) {
  mode <- structure$components$mode[[1]] %||% "fixed"
  comp_table <- structure$components
  if (!identical(mode, "sample")) {
    weights <- base_weights[match(available_components, comp_table$component_id)]
    if (length(weights) != length(available_components) || all(is.na(weights))) {
      weights <- rep(1 / length(available_components), length(available_components))
    } else {
      if (any(!is.finite(weights) | weights < 0)) {
        stop("Component weights must be non-negative finite numbers")
      }
      total <- sum(weights)
      if (!is.finite(total) || total <= 0) {
        weights <- rep(1 / length(available_components), length(available_components))
      } else {
        weights <- weights / total
      }
    }
    return(weights)
  }
  ref_id <- comp_table$reference[[1]] %||% available_components[[1]]
  comp_attrs <- comp_table$attrs
  weights <- numeric(length(available_components))
  names(weights) <- available_components
  sum_nonref <- 0
  for (i in seq_along(available_components)) {
    cid <- available_components[[i]]
    if (identical(cid, ref_id)) next
    attrs <- comp_attrs[[cid]] %||% list()
    wp <- attrs$weight_param
    if (is.null(wp) || !nzchar(wp)) {
      stop(sprintf("Component '%s' missing weight_param for sampled mixture", cid))
    }
    val <- NA_real_
    if (!is.null(trial_rows) && wp %in% names(trial_rows)) {
      vals <- trial_rows[[wp]]
      if (!is.null(vals) && length(vals) > 0) val <- as.numeric(vals[[1]])
    }
    if (is.na(val)) {
      idx <- which(comp_table$component_id == cid)
      if (length(idx) == 1) val <- comp_table$weight[[idx]]
    }
    if (is.na(val) || !is.finite(val) || val < 0 || val > 1) {
      stop(sprintf("Mixture weight '%s' must be a probability in [0,1]", wp))
    }
    weights[[i]] <- val
    sum_nonref <- sum_nonref + val
  }
  if (ref_id %in% available_components) {
    ref_idx <- match(ref_id, available_components)
    ref_weight <- 1 - sum_nonref
    if (ref_weight < -1e-8) stop("Mixture weights sum to >1; check non-reference weight params")
    if (ref_weight < 0) ref_weight <- 0
    weights[[ref_idx]] <- ref_weight
  }
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    weights <- rep(1 / length(available_components), length(available_components))
  } else {
    weights <- weights / total
  }
  weights
}

.resolve_outcome_def_index <- function(prep, outcome_label, component) {
  all_names <- names(prep[["outcomes"]])
  idxs <- which(all_names == outcome_label)
  if (length(idxs) == 0) {
    return(NA_integer_)
  }
  comp_label <- component
  if (is.factor(comp_label)) comp_label <- as.character(comp_label)
  if (length(comp_label) == 0L || is.null(comp_label) || is.na(comp_label) || !nzchar(comp_label)) {
    return(idxs[[1]])
  }
  # First pass: explicit component match
  for (i in idxs) {
    opts <- prep[["outcomes"]][[i]][["options"]] %||% list()
    comps <- opts[["component"]]
    if (!is.null(comps) && comp_label %in% comps) {
      return(i)
    }
  }
  # Second pass: unrestricted outcome
  for (i in idxs) {
    opts <- prep[["outcomes"]][[i]][["options"]] %||% list()
    if (is.null(opts[["component"]]) || length(opts[["component"]]) == 0L) {
      return(i)
    }
  }
  NA_integer_
}

.likelihood_response_prob_component <- function(prep, outcome_label, component,
                                                trial_rows = NULL,
                                                trial_state = NULL) {
  selected_idx <- .resolve_outcome_def_index(prep, outcome_label, component)
  if (is.na(selected_idx)) {
    return(0.0)
  }
  out_def <- prep[["outcomes"]][[selected_idx]]
  if (!is.null(out_def[["options"]][["alias_of"]])) {
    refs <- as.character(out_def[["options"]][["alias_of"]])
    vals <- vapply(refs, function(lbl) {
      .likelihood_response_prob_component(
        prep = prep,
        outcome_label = lbl,
        component = component,
        trial_rows = trial_rows,
        trial_state = trial_state
      )
    }, numeric(1))
    return(sum(vals, na.rm = TRUE))
  }

  # Native-only probability (upper = Inf) so shared triggers and guards follow the
  # same path as the main likelihood.
  native_ctx <- .prep_native_context(prep)
  if (!inherits(native_ctx, "externalptr")) {
    stop("Native context required for response probabilities; R fallback removed", call. = FALSE)
  }

  compiled <- .expr_lookup_compiled(out_def[["expr"]], prep)
  if (is.null(compiled)) {
    stop(sprintf("No compiled node for outcome '%s'", outcome_label), call. = FALSE)
  }

  competitor_map <- .prep_competitors(prep) %||% list()
  use_indexed_competitors <- isTRUE(attr(competitor_map, "by_index")) &&
    length(competitor_map) == length(prep[["outcomes"]])
  comp_exprs <- if (use_indexed_competitors) {
    competitor_map[[selected_idx]] %||% list()
  } else {
    competitor_map[[outcome_label]] %||% list()
  }
  comp_ids <- integer(0)
  if (length(comp_exprs) > 0L) {
    comp_nodes <- lapply(comp_exprs, function(ex) .expr_lookup_compiled(ex, prep))
    if (!any(vapply(comp_nodes, is.null, logical(1)))) {
      comp_ids <- vapply(comp_nodes, function(node) as.integer(node$id %||% NA_integer_), integer(1))
    }
  }

  trial_rows_df <- if (is.null(trial_rows)) data.frame() else as.data.frame(trial_rows)
  component_idx <- .component_index_from_prep(prep, component)
  prob_native <- tryCatch(
    native_outcome_probability_params_cpp_idx(
      native_ctx,
      as.integer(compiled$id),
      Inf,
      as.integer(component_idx),
      integer(0),
      integer(0),
      as.integer(comp_ids),
      .integrate_rel_tol(),
      .integrate_abs_tol(),
      12L,
      trial_rows_df
    ),
    error = function(e) NA_real_
  )
  if (!is.finite(prob_native) || prob_native < 0) {
    stop(sprintf("Native outcome probability failed for outcome '%s'", outcome_label), call. = FALSE)
  }
  base <- as.numeric(prob_native)

  if (!identical(outcome_label, "GUESS")) {
    gp <- .get_component_attr(prep, component, "guess")
    if (!is.null(gp) && !is.null(gp[["weights"]])) {
      keep_key <- if (is.na(outcome_label)) "<NA>" else as.character(outcome_label)
      keep <- gp[["weights"]][[outcome_label]] %||% gp[["weights"]][[keep_key]] %||% 1.0
      base <- base * as.numeric(keep)
    }
  }
  base
}

.aggregate_observed_probs <- function(prep, probs, include_na = TRUE, component = NULL) {
  labels <- unique(names(prep[["outcomes"]]))
  labels <- Filter(function(lbl) {
    idx <- .resolve_outcome_def_index(prep, lbl, component)
    if (is.na(idx)) return(FALSE)
    is.null(prep[["outcomes"]][[idx]][["options"]][["alias_of"]])
  }, labels)
  obs <- numeric(0)
  na_sum <- 0.0
  for (lbl in labels) {
    idx <- .resolve_outcome_def_index(prep, lbl, component)
    if (is.na(idx)) next
    out_def <- prep[["outcomes"]][[idx]]
    map_to <- out_def[["options"]][["map_outcome_to"]] %||% NULL
    prob <- probs[[lbl]] %||% 0.0
    if (is.null(map_to)) {
      current <- obs[lbl]
      if (length(current) == 0L || is.na(current)) current <- 0.0
      obs[lbl] <- current + prob
    } else {
      map_label <- as.character(map_to)[1]
      if (!nzchar(map_label) || identical(map_label, "NA") || is.na(map_label)) {
        na_sum <- na_sum + prob
      } else {
        current <- obs[map_label]
        if (length(current) == 0L || is.na(current)) current <- 0.0
        obs[map_label] <- current + prob
      }
    }
  }
  if (include_na && (na_sum > 0 || length(obs) == 0L)) {
    obs <- c(obs, "NA" = na_sum)
  }
  if (include_na) {
    total <- sum(obs)
    resid <- 1.0 - total
    if (!is.finite(resid)) resid <- 0.0
    if (resid > .Machine$double.eps) {
      if ("NA" %in% names(obs)) {
        obs[["NA"]] <- obs[["NA"]] + resid
      } else {
        obs <- c(obs, "NA" = resid)
      }
    }
  }
  obs
}

#' @rdname response_probabilities
#' @export
response_probabilities.default <- function(structure, ...) {
  stop("response_probabilities() expects a model_structure", call. = FALSE)
}

#' Compute predicted response probabilities
#'
#' This function returns the model-implied response probabilities for a given
#' parameter set. It is useful when you want the predicted response distribution
#' without evaluating a full trial-by-trial likelihood.
#'
#' @param structure Finalized model structure.
#' @param params_df Parameter data frame for one or more trials.
#' @param include_na If `TRUE`, include any leftover probability mass assigned
#'   to missing responses.
#' @param ... Unused; for S3 compatibility.
#' @return A named numeric vector of response probabilities.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_param_matrix(
#'   spec,
#'   c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0),
#'   n_trials = 1
#' )
#' response_probabilities(structure, params_df)
#' @export
response_probabilities <- function(structure,
                                   params_df,
                                   include_na = TRUE,
                                   ...) {
  UseMethod("response_probabilities")
}

#' @rdname response_probabilities
#' @export
response_probabilities.model_structure <- function(structure,
                                                   params_df,
                                                   include_na = TRUE,
                                                   ...) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row")
  }
  structure <- .as_model_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("model structure must include model_spec; rebuild with finalize_model")
  }
  if (is.matrix(params_df) || inherits(params_df, "param_matrix")) {
    params_df <- .param_matrix_to_rows(structure, params_df)
  }
  model_spec <- .model_spec_with_params(structure$model_spec, params_df)
  prep_eval_base <- .prepare_model_for_likelihood(model_spec)
  comp_ids <- structure$components$component_id

  if (!"trial" %in% names(params_df)) params_df$trial <- 1L
  params_df$trial <- params_df$trial
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }

  trial_ids <- sort(unique(params_df$trial))
  accum <- list()
  weights_per_trial <- numeric(length(trial_ids))

  for (i in seq_along(trial_ids)) {
    tid <- trial_ids[[i]]
    trial_rows <- params_df[params_df$trial == tid, , drop = FALSE]
    trial_state <- .eval_state_create()
    comps <- comp_ids
    if ("component" %in% names(trial_rows)) {
      listed <- unique(trial_rows$component)
      listed <- listed[!is.na(listed)]
      if (length(listed) > 0L) comps <- intersect(comp_ids, listed)
    }
    if (length(comps) == 0L) comps <- comp_ids
    weights <- .likelihood_component_weights(
      structure,
      tid,
      comps,
      structure$components$weight,
      trial_rows = trial_rows
    )
    trial_probs <- numeric(0)
    for (idx in seq_along(comps)) {
      comp_id <- comps[[idx]]
      comp_rows <- .likelihood_component_rows(trial_rows, comp_id)
      labels <- unique(names(prep_eval_base$outcomes))
      base_probs <- setNames(vapply(labels, function(lbl) {
        .likelihood_response_prob_component(
          prep_eval_base,
          lbl,
          comp_id,
          trial_rows = comp_rows,
          trial_state = trial_state
        )
      }, numeric(1)), labels)
      comp_probs <- .aggregate_observed_probs(
        prep_eval_base,
        base_probs,
        include_na = include_na,
        component = comp_id
      )
      factor <- weights[[idx]]
      weighted <- factor * comp_probs
      if (length(trial_probs) == 0L) {
        trial_probs <- weighted
      } else {
        for (lbl in names(weighted)) {
          current <- trial_probs[lbl]
          if (length(current) == 0L || is.na(current)) current <- 0
          trial_probs[lbl] <- current + weighted[[lbl]]
        }
      }
    }
    if (include_na) {
      total <- sum(trial_probs)
      resid <- 1.0 - total
      if (!is.finite(resid)) resid <- 0.0
      if (resid > .Machine$double.eps) {
        current <- trial_probs[["NA"]]
        if (is.null(current) || is.na(current)) current <- 0.0
        trial_probs[["NA"]] <- current + resid
      }
    }
    accum[[i]] <- trial_probs
    weights_per_trial[[i]] <- 1.0
  }
  total_weight <- sum(weights_per_trial)
  combined_labels <- Reduce(union, lapply(accum, names))
  result <- setNames(numeric(length(combined_labels)), combined_labels)
  for (i in seq_along(accum)) {
    vec <- accum[[i]]
    w <- weights_per_trial[[i]] / total_weight
    for (lbl in names(vec)) {
      current <- result[lbl]
      if (length(current) == 0L || is.na(current)) current <- 0
      result[lbl] <- current + vec[[lbl]] * w
    }
  }
  result
}


#' Evaluate the log-likelihood of behavioral data
#'
#' Compute the log-likelihood of prepared behavioral data under one or more
#' candidate parameter sets.
#'
#' @param context Context created with `make_context()`.
#' @param data Prepared data created with `prepare_data()`.
#' @param parameters A parameter data frame, or a list of parameter data frames.
#' @param ok Logical vector marking which trials should contribute to the
#'   likelihood. Trials marked `FALSE` are assigned `min_ll`.
#' @param expand Optional index vector used to expand compressed trial-level
#'   results back to the original trial count.
#' @param min_ll Minimum log-likelihood value used for excluded or impossible
#'   trials.
#' @param ... Unused; for S3 compatibility.
#' @return A numeric vector of log-likelihood values.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_param_matrix(
#'   spec,
#'   c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' prepared <- prepare_data(structure, data_df)
#' ctx <- make_context(structure)
#' log_likelihood(ctx, prepared, list(params_df))
#' @export
log_likelihood <- function(context, data, parameters, ok = NULL, expand = NULL, min_ll = log(1e-10), ...) {
  UseMethod("log_likelihood")
}

#' @rdname log_likelihood
#' @export
log_likelihood.accumulatr_context <- function(context,
                                              data,
                                              parameters,
                                              ok = NULL,
                                              expand = NULL,
                                              min_ll = log(1e-10),
                                              ...) {
  ctx <- .validate_context(context)
  data_df <- .validate_prepared_data(data)
  native_ctx <- ctx$native_ctx
  if (is.null(native_ctx) || !inherits(native_ctx, "externalptr")) {
    stop("Native context required for log_likelihood", call. = FALSE)
  }
  params_list <- if (is.data.frame(parameters) || is.matrix(parameters) || inherits(parameters, "param_matrix")) {
    list(parameters)
  } else {
    parameters
  }
  params_list <- lapply(params_list, function(pm) {
    if (inherits(pm, "param_matrix") || is.matrix(pm)) {
      pm
    } else {
      as.matrix(pm)
    }
  })
  trial <- data_df$trial
  n_trials <- if (length(trial) == 0L) {
    0L
  } else {
    1L + sum(trial[-1L] != trial[-length(trial)])
  }
  if (is.null(expand) || length(expand) == 0L) {
    expand <- seq_len(n_trials)
  }
  expand <- as.integer(expand)
  if (is.null(ok) || length(ok) == 0L) {
    ok <- rep_len(TRUE, n_trials)
  }
  ok <- as.logical(ok)
  ok[is.na(ok)] <- FALSE
  cpp_loglik_multiple(
    native_ctx,
    params_list,
    data_df,
    ok,
    expand,
    min_ll,
    .integrate_rel_tol(),
    .integrate_abs_tol(),
    getOption("uuber.integrate.max.depth", 12L)
  )
}

#' @rdname log_likelihood
#' @export
log_likelihood.default <- function(context, ...) {
  stop("log_likelihood() expects a context created with make_context()", call. = FALSE)
}
