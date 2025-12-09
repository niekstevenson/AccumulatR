.acc_numeric_scalar <- function(value, default, name) {
  if (is.null(value) || length(value) == 0L) return(default)
  val <- as.numeric(value)[[1]]
  if (!is.finite(val)) {
    stop(sprintf("Accumulator parameter '%s' must be a finite numeric value", name),
         call. = FALSE)
  }
  val
}

.acc_extract_args <- function(acc_def) {
  dist_raw <- acc_def[["dist"]]
  if (is.null(dist_raw) || length(dist_raw) == 0L) {
    stop("Accumulator definition missing 'dist'", call. = FALSE)
  }
  dist_chr <- as.character(dist_raw)[[1]]
  if (!nzchar(dist_chr)) {
    stop("Accumulator definition has empty distribution name", call. = FALSE)
  }
  dist_chr <- tolower(dist_chr)

  params <- acc_def[["params"]]
  if (is.null(params) || length(params) == 0L) {
    stop(sprintf("Accumulator '%s' missing parameter list", dist_chr), call. = FALSE)
  }
  if (!is.list(params)) {
    params <- as.list(params)
  }

  list(
    dist = dist_chr,
    params = params,
    onset = .acc_numeric_scalar(acc_def[["onset"]], 0, "onset"),
    q = .acc_numeric_scalar(acc_def[["q"]], 0, "q")
  )
}

.acc_eval <- function(fn, acc_def, t) {
  t_vals <- as.numeric(t)
  if (length(t_vals) == 0L) return(numeric(0))
  args <- .acc_extract_args(acc_def)
  res <- vapply(
    t_vals,
    function(tt) fn(tt, args$onset, args$q, args$dist, args$params),
    numeric(1),
    USE.NAMES = FALSE
  )
  if (length(res) == 1L) res[[1]] else res
}

.acc_density <- function(acc_def, t) {
  .acc_eval(acc_density_cpp, acc_def, t)
}

.acc_density_success <- function(acc_def, t) {
  .acc_eval(acc_density_success_cpp, acc_def, t)
}

.acc_survival <- function(acc_def, t) {
  .acc_eval(acc_survival_cpp, acc_def, t)
}

.acc_cdf_success <- function(acc_def, t) {
  .acc_eval(acc_cdf_success_cpp, acc_def, t)
}

.acc_shared_trigger_id <- function(acc_def) {
  sid <- acc_def[['shared_trigger_id']]
  if (is.null(sid) || is.na(sid) || identical(sid, "")) return(NA_character_)
  as.character(sid)[[1]]
}

.pool_active_members <- function(prep, pool_id, component) {
  pool_defs <- prep[["pools"]]
  acc_defs <- prep[["accumulators"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))

  runtime_cache <- .prep_runtime_get(prep, "pool_members_cache")
  component_key <- if (is.null(component) || length(component) == 0L) {
    "__default__"
  } else {
    as.character(component)[[1]]
  }
  cached_members <- NULL
  if (is.environment(runtime_cache)) {
    pool_cache <- runtime_cache[[pool_id]]
    if (is.environment(pool_cache)) {
      cached_members <- pool_cache[[component_key]]
    }
  }
  if (!is.null(cached_members)) return(cached_members)

  members <- pool_def[["members"]]
  active <- character(0)
  for (m in members) {
    if (!is.null(acc_defs[[m]])) {
      comps <- acc_defs[[m]][["components"]]
      if (length(comps) == 0 || is.null(component) ||
          identical(component, "__default__") || component %in% comps) {
        active <- c(active, m)
      }
    } else if (!is.null(pool_defs[[m]])) {
      active <- c(active, m)
    }
  }
  if (is.environment(runtime_cache)) {
    pool_cache <- runtime_cache[[pool_id]]
    if (!is.environment(pool_cache)) {
      pool_cache <- new.env(parent = emptyenv(), hash = TRUE)
      runtime_cache[[pool_id]] <- pool_cache
    }
    pool_cache[[component_key]] <- active
  }
  active
}

.pool_density_fast_value <- function(prep, pool_id, component, t) {
  if (is.na(t) || !is.finite(t) || t < 0) return(0.0)
  pool_defs <- prep[["pools"]]
  acc_defs <- prep[["accumulators"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) return(0.0)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) return(0.0)
  k <- as.integer(pool_def[["k"]] %||% 1L)
  if (k < 1L || k > n) return(0.0)
  dens_vec <- numeric(n)
  surv_vec <- numeric(n)
  for (i in seq_len(n)) {
    mid <- members[[i]]
    if (!is.null(acc_defs[[mid]])) {
      acc_def <- acc_defs[[mid]]
      dens_vec[[i]] <- .acc_density(acc_def, t)
      surv_vec[[i]] <- .acc_survival(acc_def, t)
    } else if (!is.null(pool_defs[[mid]])) {
      dens_vec[[i]] <- .pool_density_fast_value(prep, mid, component, t)
      surv_vec[[i]] <- .pool_survival_fast_value(prep, mid, component, t)
    } else {
      dens_vec[[i]] <- 0.0
      surv_vec[[i]] <- 1.0
    }
  }
  total_density <- pool_density_fast_cpp(dens_vec, surv_vec, as.integer(k))
  if (!is.finite(total_density) || total_density <= 0) 0.0 else as.numeric(total_density)
}

.pool_survival_fast_value <- function(prep, pool_id, component, t) {
  if (is.na(t)) return(1.0)
  pool_defs <- prep[["pools"]]
  acc_defs <- prep[["accumulators"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) return(1.0)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) return(1.0)
  k <- as.integer(pool_def[["k"]] %||% 1L)
  if (k > n) return(1.0)
  if (k < 1L) return(0.0)

  if (!is.finite(t)) {
    if (t < 0) return(1.0)
    t <- Inf
  }

  Svec <- numeric(n)
  for (i in seq_len(n)) {
    mid <- members[[i]]
    surv <- if (!is.null(acc_defs[[mid]])) {
      .acc_survival(acc_defs[[mid]], t)
    } else if (!is.null(pool_defs[[mid]])) {
      .pool_survival_fast_value(prep, mid, component, t)
    } else {
      1.0
    }
    Svec[[i]] <- surv
  }
  val <- pool_survival_fast_cpp(Svec, as.integer(k))
  if (!is.finite(val)) 0.0 else as.numeric(val)
}

.build_pool_templates <- function(pool_id, members, member_ids, pool_idx, k) {
  member_ids_int <- as.integer(member_ids)
  pool_idx_int <- if (is.na(pool_idx)) NA_integer_ else as.integer(pool_idx)
  pool_build_templates_cpp(length(members), member_ids_int, pool_idx_int, as.integer(k))
}

.event_survival_at <- function(prep, id, component, t,
                               forced_complete = integer(0),
                               forced_survive = integer(0)) {
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) return(1.0)
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) return(0.0)
  }
  is_pool <- !is.null(prep[["pools"]][[id]])
  if (is_pool) {
    # Fast-path when no conditional forcing is present
    if ((length(forced_complete) == 0L) && (length(forced_survive) == 0L)) {
      return(.pool_survival_fast_value(prep, id, component, t))
    }
    return(.pool_survival(prep, id, component, t,
                          forced_complete = forced_complete,
                          forced_survive = forced_survive))
  }
  acc_def <- prep[["accumulators"]][[id]]
  if (!is.null(acc_def)) return(.acc_survival(acc_def, t))
  0.0
}

.event_cdf_at <- function(prep, id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0)) {
  1.0 - .event_survival_at(prep, id, component, t,
                           forced_complete = forced_complete,
                           forced_survive = forced_survive)
}

.event_density_at <- function(prep, id, component, t,
                              forced_complete = integer(0),
                              forced_survive = integer(0)) {
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
  }

  pool_def <- prep[["pools"]][[id]]
  if (!is.null(pool_def)) {
    # Fast-path when no conditional forcing is present
    if ((length(forced_complete) == 0L) && (length(forced_survive) == 0L)) {
      dens <- .pool_density_fast_value(prep, id, component, t)
      attr(dens, "scenarios") <- list()
      return(dens)
    }
    return(.pool_density(
      prep, id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    ))
  }

  acc_def <- prep[["accumulators"]][[id]]
  if (!is.null(acc_def)) {
    dens <- .acc_density(acc_def, t)
    attr(dens, "scenarios") <- list()
    return(dens)
  }

  val <- 0.0
  attr(val, "scenarios") <- list()
  val
}

.pool_density <- function(prep, pool_id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0)) {
  pool_defs <- prep[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  k <- as.integer(pool_def[["k"]] %||% 1L)
  if (k < 1L || k > n) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }

  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)

  member_ids <- .labels_to_ids(prep, members)
  pool_idx <- .label_to_id(prep, pool_id)
  bundle_info <- .cache_bundle_ensure(prep)
  prep <- bundle_info$prep
  bundle <- bundle_info$bundle
  template_key <- .pool_template_cache_key(pool_id, component, k)
  templates_info <- .bundle_pool_templates_get(bundle, template_key)
  if (is.null(templates_info)) {
    templates_info <- .build_pool_templates(pool_id, members, member_ids, pool_idx, k)
    .bundle_pool_templates_set(bundle, template_key, templates_info)
  }
  templates <- templates_info$templates
  finisher_map <- templates_info$finisher_map

  dens_vec <- numeric(n)
  cdf_vec <- numeric(n)
  surv_vec <- numeric(n)
  acc_defs <- prep[["accumulators"]]
  pool_defs <- prep[["pools"]]
  for (i in seq_len(n)) {
    mid <- members[[i]]
    dens_vec[[i]] <- if (!is.null(acc_defs[[mid]])) {
      acc_def <- acc_defs[[mid]]
      shared_id <- .acc_shared_trigger_id(acc_def)
      if (!is.na(shared_id)) {
        q_val <- acc_def[['q']] %||% 0
        (1 - q_val) * .acc_density_success(acc_def, t)
      } else {
        .acc_density(acc_def, t)
      }
    } else if (!is.null(pool_defs[[mid]])) {
      .pool_density(prep, mid, component, t,
                    forced_complete = forced_complete,
                    forced_survive = forced_survive)
    } else {
      0.0
    }
    cdf_vec[[i]] <- if (!is.null(acc_defs[[mid]])) {
      acc_def <- acc_defs[[mid]]
      shared_id <- .acc_shared_trigger_id(acc_def)
      if (!is.na(shared_id)) {
        .acc_cdf_success(acc_def, t)
      } else {
        1.0 - .acc_survival(acc_def, t)
      }
    } else if (!is.null(pool_defs[[mid]])) {
      .event_cdf_at(prep, mid, component, t,
                    forced_complete = forced_complete,
                    forced_survive = forced_survive)
    } else {
      0.0
    }
    surv_vec[[i]] <- 1.0 - cdf_vec[[i]]
  }

  shared_ids <- vapply(seq_len(n), function(i) {
    mid <- members[[i]]
    acc_def <- acc_defs[[mid]]
    if (is.null(acc_def)) NA_character_ else .acc_shared_trigger_id(acc_def)
  }, character(1))

  cdf_success_vec <- numeric(n)
  surv_success_vec <- numeric(n)
  for (i in seq_len(n)) {
    mid <- members[[i]]
    acc_def <- acc_defs[[mid]]
    if (!is.null(acc_def) && !is.na(shared_ids[[i]]) && nzchar(shared_ids[[i]])) {
      cdf_success_vec[[i]] <- .acc_cdf_success(acc_def, t)
      surv_success_vec[[i]] <- 1.0 - cdf_success_vec[[i]]
    } else {
      cdf_success_vec[[i]] <- cdf_vec[[i]]
      surv_success_vec[[i]] <- surv_vec[[i]]
    }
  }

  shared_groups <- integer(n)
  if (n > 0L) {
    non_na <- which(!is.na(shared_ids) & nzchar(shared_ids))
    if (length(non_na) > 0L) {
      unique_vals <- unique(shared_ids[non_na])
      shared_groups[non_na] <- match(shared_ids[non_na], unique_vals)
    }
  }

  res <- pool_density_combine_cpp(
    dens_vec,
    cdf_vec,
    surv_vec,
    cdf_success_vec,
    surv_success_vec,
    as.integer(shared_groups),
    templates,
    as.integer(forced_complete),
    as.integer(forced_survive)
  )

  total_density <- as.numeric(res$value)
  if (!is.finite(total_density) || total_density < 0) total_density <- 0.0

  scenarios_raw <- res$scenarios %||% list()
  scenarios <- list()
  if (length(scenarios_raw) > 0L) {
    for (sr in scenarios_raw) {
      sc <- .make_scenario_record(
        prep,
        sr$weight,
        sr$forced_complete,
        sr$forced_survive
      )
      if (!is.null(sc)) scenarios[[length(scenarios) + 1L]] <- sc
    }
  }
  attr(total_density, "scenarios") <- scenarios
  total_density
}

.pool_survival <- function(prep, pool_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0)) {
  if (!is.finite(t)) return(0.0)
  pool_defs <- prep[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  k <- as.integer(pool_def[["k"]] %||% 1L)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) return(1.0)
  if (k > n) return(1.0)
  if (k < 1L) return(0.0)

  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)

  Fvec <- numeric(n)
  for (i in seq_len(n)) {
    id <- members[[i]]
    Fvec[[i]] <- .event_cdf_at(prep, id, component, t,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive)
  }
  as.numeric(pool_survival_general_cpp(Fvec, as.integer(k)))
}
