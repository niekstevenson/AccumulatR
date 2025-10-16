`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

is_model_tables <- function(x) {
  is.list(x) &&
    all(c("struct_nodes", "component_nodes", "param_values") %in% names(x))
}

.make_param_id <- function(struct_id, param_name) {
  sprintf("%s::%s", struct_id, param_name)
}

model_to_tables <- function(model) {
  if (inherits(model, "race_spec") && exists("build_model", mode = "function")) {
    model <- build_model(model)
  }
  if (!is.list(model) || is.null(model$accumulators) || is.null(model$outcomes)) {
    stop("model_to_tables expects a race model specification")
  }

  components_info <- list(ids = "__default__", weights = 1, attrs = list(list()))
  if (!is.null(model$metadata) && !is.null(model$metadata$mixture) &&
      !is.null(model$metadata$mixture$components)) {
    comps <- model$metadata$mixture$components
    if (length(comps) > 0) {
      components_info$ids <- vapply(comps, `[[`, character(1), "id")
      components_info$weights <- vapply(comps, function(cmp) cmp$weight %||% 1, numeric(1))
      components_info$attrs <- lapply(comps, function(cmp) cmp$attrs %||% list())
      weight_params <- vapply(comps, function(cmp) cmp$attrs$weight_param %||% NA_character_, character(1))
    } else {
      weight_params <- NA_character_
    }
  } else {
    weight_params <- NA_character_
  }

  struct_rows <- list()
  param_rows <- list()

  add_struct_row <- function(struct_id, label, node_type,
                             dist_family = NA_character_,
                             calc_mode = NA_character_,
                             param_slots = NULL,
                             payload = NULL) {
    struct_rows[[length(struct_rows) + 1L]] <<- data.frame(
      struct_id = struct_id,
      label = label,
      node_type = node_type,
      dist_family = dist_family,
      calc_mode = calc_mode,
      param_slots = I(list(param_slots %||% list())),
      payload = I(list(payload %||% list())),
      stringsAsFactors = FALSE
    )
  }

  shared_param_map <- list()
  shared_param_values <- list()
  shared_group_counter <- 0L
  if (!is.null(model$groups) && length(model$groups) > 0) {
    for (grp in model$groups) {
      sp <- grp$attrs$shared_params %||% NULL
      if (is.null(sp) || length(sp) == 0) next
      shared_group_counter <- shared_group_counter + 1L
      gid <- grp$id %||% sprintf("shared_group_%d", shared_group_counter)
      members <- grp$members %||% character(0)
      if (length(members) == 0) next
      for (member in members) {
        if (is.null(shared_param_map[[member]])) shared_param_map[[member]] <- list()
        for (pname in names(sp)) {
          shared_id <- sprintf("shared:%s:%s", gid, pname)
          shared_param_map[[member]][[pname]] <- shared_id
          if (!is.null(sp[[pname]])) {
            if (is.null(shared_param_values[[shared_id]])) {
              shared_param_values[[shared_id]] <- sp[[pname]]
            }
          }
        }
      }
    }
  }

  added_param_ids <- character(0)

  for (acc in model$accumulators) {
    struct_id <- sprintf("acc:%s", acc$id)
    param_slots <- list()
    if (!is.null(acc$params) && length(acc$params) > 0) {
      for (pname in names(acc$params)) {
        shared_mapping <- shared_param_map[[acc$id]] %||% list()
        param_id <- shared_mapping[[pname]] %||% .make_param_id(struct_id, pname)
        param_slots[[pname]] <- param_id
        if (!(param_id %in% added_param_ids)) {
          param_value <- shared_param_values[[param_id]]
          if (is.null(param_value)) {
            param_value <- acc$params[[pname]]
            shared_param_values[[param_id]] <- param_value
          }
          param_rows[[length(param_rows) + 1L]] <- data.frame(
            param_id = param_id,
            value = param_value,
            stringsAsFactors = FALSE
          )
          added_param_ids <- c(added_param_ids, param_id)
        }
      }
    }
    payload <- list(
      onset = acc$onset %||% 0,
      q = acc$q %||% 0,
      components = acc$components %||% NULL
    )
    add_struct_row(struct_id, acc$id, "accumulator", acc$dist, "density",
                   param_slots = param_slots, payload = payload)
  }

  if (!is.null(model$pools) && length(model$pools) > 0) {
    for (pool in model$pools) {
      struct_id <- sprintf("pool:%s", pool$id)
      payload <- list(
        members = pool$members %||% character(0),
        rule = pool$rule %||% list()
      )
      add_struct_row(struct_id, pool$id, "pool", NA_character_, "density",
                     payload = payload)
    }
  }

  for (outcome in model$outcomes) {
    struct_id <- sprintf("outcome:%s", outcome$label)
    payload <- list(
      label = outcome$label,
      expr = outcome$expr,
      options = outcome$options %||% list()
    )
    add_struct_row(struct_id, outcome$label, "outcome", NA_character_, "probability",
                   payload = payload)
  }

  struct_nodes <- if (length(struct_rows) > 0) do.call(rbind, struct_rows) else
    data.frame(struct_id = character(0), label = character(0),
               node_type = character(0), dist_family = character(0),
               calc_mode = character(0), param_slots = I(list()),
               payload = I(list()), stringsAsFactors = FALSE)

  param_values <- if (length(param_rows) > 0) do.call(rbind, param_rows) else
    data.frame(param_id = character(0), value = numeric(0),
               stringsAsFactors = FALSE)

  struct_edges_df <- data.frame(
    edge_id = character(0),
    parent_struct_id = character(0),
    child_struct_id = character(0),
    role = character(0),
    weight_default = numeric(0),
    order_idx = integer(0),
    stringsAsFactors = FALSE
  )

  component_ids <- components_info$ids
  component_nodes <- expand.grid(
    component_id = component_ids,
    struct_id = struct_nodes$struct_id,
    stringsAsFactors = FALSE
  )

  component_nodes$active <- 1L
  component_nodes$payload_override <- I(vector("list", nrow(component_nodes)))

  for (i in seq_len(nrow(struct_nodes))) {
    payload <- struct_nodes$payload[[i]] %||% list()
    allowed <- payload$components %||% NULL
    if (!is.null(allowed) && length(allowed) > 0) {
      struct_id <- struct_nodes$struct_id[[i]]
      idx <- component_nodes$struct_id == struct_id
      component_nodes$active[idx] <- ifelse(component_nodes$component_id[idx] %in% allowed, 1L, 0L)
    }
  }

  components_df <- data.frame(
    component_id = component_ids,
    weight = components_info$weights,
    weight_param = if (length(weight_params) == length(component_ids)) weight_params else rep(NA_character_, length(component_ids)),
    attrs = I(components_info$attrs),
    stringsAsFactors = FALSE
  )

  tables <- list(
    struct_nodes = struct_nodes,
    struct_edges = struct_edges_df,
    component_nodes = component_nodes,
    component_edges = data.frame(
      component_id = character(0),
      edge_id = character(0),
      active = integer(0),
      weight_override = numeric(0),
      stringsAsFactors = FALSE
    ),
    param_values = param_values,
    components = components_df,
    metadata = model$metadata %||% list(),
    groups = model$groups %||% list()
  )

  class(tables) <- c("model_tables", class(tables))
  tables
}

tables_to_model <- function(tables) {
  if (!is_model_tables(tables)) stop("tables_to_model expects output from model_to_tables")

  struct_nodes <- tables$struct_nodes
  param_lookup <- setNames(tables$param_values$value, tables$param_values$param_id)

  get_param <- function(pid) {
    if (is.null(pid) || is.na(pid) || pid == "") return(NA_real_)
    if (!pid %in% names(param_lookup)) {
      stop(sprintf("Parameter '%s' not found in param_values", pid))
    }
    param_lookup[[pid]]
  }

  accumulator_rows <- struct_nodes[struct_nodes$node_type == "accumulator", , drop = FALSE]
  accumulators <- lapply(seq_len(nrow(accumulator_rows)), function(i) {
    row <- accumulator_rows[i, ]
    slots <- row$param_slots[[1]] %||% list()
    params <- if (length(slots) > 0) {
      out <- lapply(slots, get_param)
      names(out) <- names(slots)
      out
    } else list()
    payload <- row$payload[[1]] %||% list()
    list(
      id = row$label,
      dist = row$dist_family,
      onset = payload$onset %||% 0,
      q = payload$q %||% 0,
      params = params,
      components = payload$components %||% NULL
    )
  })

  pool_rows <- struct_nodes[struct_nodes$node_type == "pool", , drop = FALSE]
  pools <- lapply(seq_len(nrow(pool_rows)), function(i) {
    row <- pool_rows[i, ]
    payload <- row$payload[[1]] %||% list()
    members <- payload$members %||% character(0)
    rule <- payload$rule %||% list()
    list(
      id = row$label,
      members = members,
      rule = rule
    )
  })

  outcome_rows <- struct_nodes[struct_nodes$node_type == "outcome", , drop = FALSE]
  outcomes <- lapply(seq_len(nrow(outcome_rows)), function(i) {
    row <- outcome_rows[i, ]
    payload <- row$payload[[1]] %||% list()
    list(
      label = payload$label %||% row$label,
      expr = payload$expr,
      options = payload$options %||% list()
    )
  })

  metadata <- tables$metadata %||% list()
  groups <- tables$groups %||% list()

  if (!is.null(tables$components) && nrow(tables$components) > 0) {
    comps <- lapply(seq_len(nrow(tables$components)), function(i) {
      list(
        id = tables$components$component_id[[i]],
        weight = tables$components$weight[[i]],
        attrs = tables$components$attrs[[i]] %||% list()
      )
    })
    metadata$mixture <- list(components = comps)
  }

  model <- list(
    accumulators = accumulators,
    pools = pools,
    outcomes = outcomes,
    metadata = metadata,
    groups = groups
  )
  class(model) <- "race_model_spec"
  model
}

simulate_model_from_tables <- function(model_tables, n_trials, seed = NULL,
                                       keep_detail = FALSE) {
  if (!missing(seed) && !is.null(seed)) set.seed(seed)
  simulate_model(tables_to_model(model_tables), n_trials = n_trials,
                 seed = NULL, keep_detail = keep_detail)
}
