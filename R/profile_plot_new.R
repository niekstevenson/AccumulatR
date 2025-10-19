# Profile plot system for new API likelihood
# Tests parameter identifiability by varying each parameter around its true value

#' Extract all parameters from a model
#'
#' @param model Model specification
#' @return List with parameter paths and values
.extract_model_params <- function(model, include_component_weights = NULL) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  prep <- prepare_model(model)

  params <- list()
  param_info <- list()

  # Detect shared parameter groups from the original model definition
  group_list <- model$groups %||% list()
  shared_groups <- list()
  skip_pairs <- list()  # (acc_id, par_name) combinations to skip extracting per-acc
  if (length(group_list) > 0) {
    for (grp in group_list) {
      attrs <- grp$attrs %||% list()
      sp <- attrs$shared_params %||% NULL
      if (!is.null(sp) && length(sp) > 0) {
        gid <- grp$id
        members <- as.character(grp$members)
        for (par_name in names(sp)) {
          par_val <- sp[[par_name]]
          if (is.numeric(par_val) && length(par_val) == 1 && is.finite(par_val)) {
            key <- sprintf("shared:%s:%s", gid, par_name)
            params[[key]] <- par_val
            param_info[[key]] <- list(
              type = "shared_group_param",
              group_id = gid,
              members = members,
              par_name = par_name,
              label = sprintf("%s[%s]", par_name, gid)
            )
            # mark member params to skip individual extraction
            for (acc_id in members) {
              skip_pairs[[paste(acc_id, par_name, sep = "::")]] <- TRUE
            }
          }
        }
      }
    }
  }

  # Extract accumulator parameters, skipping those covered by shared groups
  acc_defs <- prep$accumulators
  for (acc_id in names(acc_defs)) {
    acc_def <- acc_defs[[acc_id]]
    acc_params <- acc_def$params

    for (par_name in names(acc_params)) {
      if (isTRUE(skip_pairs[[paste(acc_id, par_name, sep = "::")]])) next
      par_val <- acc_params[[par_name]]
      if (is.numeric(par_val) && length(par_val) == 1 && is.finite(par_val)) {
        key <- sprintf("acc:%s:%s", acc_id, par_name)
        params[[key]] <- par_val
        param_info[[key]] <- list(
          type = "accumulator",
          acc_id = acc_id,
          par_name = par_name,
          label = sprintf("%s[%s]", par_name, acc_id)
        )
      }
    }
  }

  # Extract component weights (optional)
  comp_info <- prep$components
  if (length(comp_info$ids) > 1) {
    # Default behavior: auto-detect; include weights only if a weight_param is set in metadata
    policy <- include_component_weights %||% "auto"
    weights <- comp_info$weights
    comps_md <- NULL
    if (!is.null(model$metadata) && !is.null(model$metadata$mixture) && !is.null(model$metadata$mixture$components)) {
      comps_md <- model$metadata$mixture$components
    }
    for (i in seq_along(comp_info$ids)[-length(comp_info$ids)]) {
      comp_id <- comp_info$ids[[i]]
      has_param <- FALSE
      if (!is.null(comps_md)) {
        # find matching comp by id
        match_idx <- which(vapply(comps_md, function(x) identical(x$id, comp_id), logical(1)))
        if (length(match_idx) == 1) {
          attrs <- comps_md[[match_idx]]$attrs %||% list()
          has_param <- !is.null(attrs$weight_param)
        }
      }
      include_this <- FALSE
      if (isTRUE(policy)) include_this <- TRUE else if (identical(policy, "auto")) include_this <- has_param else include_this <- FALSE
      if (!include_this) next
      key <- sprintf("comp_weight:%s", comp_id)
      w <- weights[[i]]
      logit_w <- log(w / (1 - w))
      params[[key]] <- logit_w
      label_txt <- if (has_param && !is.null(comps_md)) {
        match_idx <- which(vapply(comps_md, function(x) identical(x$id, comp_id), logical(1)))
        wp <- if (length(match_idx) == 1) (comps_md[[match_idx]]$attrs$weight_param %||% NULL) else NULL
        if (!is.null(wp)) sprintf("%s (w[%s])", wp, comp_id) else sprintf("logit(w[%s])", comp_id)
      } else sprintf("logit(w[%s])", comp_id)
      param_info[[key]] <- list(
        type = "component_weight",
        comp_id = comp_id,
        index = i,
        label = label_txt
      )
    }
  }

  list(
    values = params,
    info = param_info,
    prep = prep
  )
}

#' Modify model parameters
#'
#' @param model Original model
#' @param param_key Parameter key from .extract_model_params
#' @param new_value New parameter value
#' @return Modified model
.modify_model_param <- function(model, param_info, param_key, new_value) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  info <- param_info[[param_key]]

  if (identical(info$type, "accumulator")) {
    # Modify accumulator parameter
    acc_id <- info$acc_id
    par_name <- info$par_name

    # Find and modify the accumulator
    for (i in seq_along(model$accumulators)) {
      if (model$accumulators[[i]]$id == acc_id) {
        model$accumulators[[i]]$params[[par_name]] <- new_value
        break
      }
    }
  } else if (identical(info$type, "component_weight")) {
    # Modify component weight (from logit back to probability)
    comp_id <- info$comp_id
    comp_idx <- info$index

    # Convert logit back to weight
    # This is tricky - we need to renormalize all weights
    # For now, just modify the target weight and renormalize
    if (!is.null(model$metadata$mixture$components)) {
      # Extract current weights
      comps <- model$metadata$mixture$components
      n_comps <- length(comps)

      # Simple approach: set this weight and distribute remainder
      w_new <- 1 / (1 + exp(-new_value))

      # Renormalize
      old_weight <- comps[[comp_idx]]$weight %||% (1 / n_comps)
      scale <- (1 - w_new) / (1 - old_weight)

      for (i in seq_along(comps)) {
        if (i == comp_idx) {
          comps[[i]]$weight <- w_new
        } else {
          old_w <- comps[[i]]$weight %||% (1 / n_comps)
          comps[[i]]$weight <- old_w * scale
        }
      }

      model$metadata$mixture$components <- comps
    }
  } else if (identical(info$type, "shared_group_param")) {
    # Update the group's shared parameter so it applies to all members
    gid <- info$group_id
    par_name <- info$par_name
    groups <- model$groups %||% list()
    if (length(groups) > 0) {
      for (i in seq_along(groups)) {
        if (identical(groups[[i]]$id, gid)) {
          if (is.null(groups[[i]]$attrs)) groups[[i]]$attrs <- list()
          if (is.null(groups[[i]]$attrs$shared_params)) groups[[i]]$attrs$shared_params <- list()
          groups[[i]]$attrs$shared_params[[par_name]] <- new_value
        }
      }
      model$groups <- groups
    }
  }

  model
}

#' Compute profile likelihood for a single parameter
#'
#' @param model Model specification with true parameter values
#' @param data Simulated data
#' @param param_key Parameter to profile
#' @param param_info Parameter information from .extract_model_params
#' @param percent Percentage range to explore (default 0.2 = ±20%)
#' @param n_points Number of points in profile (default 20)
#' @return Data frame with parameter values and log-likelihoods
.profile_single_param <- function(model, data, param_key, param_info,
                                  percent = 0.2, n_points = 20) {
  info <- param_info[[param_key]]
  true_value <- param_info[[param_key]]$true_value

  # Create grid around true value
  if (true_value == 0) {
    span <- percent
  } else {
    span <- abs(true_value) * percent
  }

  grid <- seq(true_value - span, true_value + span, length.out = n_points)

  # Compute likelihood for each point
  log_liks <- numeric(length(grid))

  for (i in seq_along(grid)) {
    # Modify model with new parameter value
    mod_model <- .modify_model_param(model, param_info, param_key, grid[[i]])

    # Compute likelihood
    log_liks[[i]] <- tryCatch({
      compute_loglik(mod_model, data)
    }, error = function(e) {
      -Inf
    })
  }

  data.frame(
    parameter = param_key,
    label = info$label,
    value = grid,
    true_value = true_value,
    loglik = log_liks,
    stringsAsFactors = FALSE
  )
}

#' Compute profile likelihoods for all parameters in a model
#'
#' @param model Model specification with true parameter values
#' @param data Simulated data from this model
#' @param percent Percentage range to explore (default 0.2 = ±20%)
#' @param n_points Number of points per parameter (default 20)
#' @param n_cores Number of cores for parallel (default 1)
#' @return Data frame with profile results
#' @export
profile_likelihood <- function(model, data,
                              percent = 0.1,
                              n_points = 10,
                              n_cores = 1,
                              include_component_weights = "auto") {

  cat("Extracting parameters...\n")
  param_extract <- .extract_model_params(model, include_component_weights = include_component_weights)
  params <- param_extract$values
  param_info <- param_extract$info

  # Add true values to param_info
  for (key in names(params)) {
    param_info[[key]]$true_value <- params[[key]]
  }

  cat(sprintf("Found %d parameters to profile\n", length(params)))

  # Compute true log-likelihood
  cat("Computing log-likelihood at true parameters...\n")
  true_ll <- compute_loglik(model, data)
  cat(sprintf("True log-likelihood: %.4f\n\n", true_ll))

  # Profile each parameter
  if (n_cores > 1) {
    cat(sprintf("Profiling parameters in parallel (%d cores)...\n", n_cores))

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      # Export necessary objects
      parallel::clusterExport(cl, c("model", "data", "param_info",
                                   "percent", "n_points",
                                   ".profile_single_param",
                                   ".modify_model_param",
                                   "compute_loglik"),
                            envir = environment())
      profiles <- parallel::parLapply(cl, names(params), function(key) {
        .profile_single_param(model, data, key, param_info, percent, n_points)
      })
    } else {
      profiles <- parallel::mclapply(names(params), function(key) {
        .profile_single_param(model, data, key, param_info, percent, n_points)
      }, mc.cores = n_cores)
    }
  } else {
    cat("Profiling parameters sequentially...\n")
    profiles <- vector("list", length(params))

    for (i in seq_along(params)) {
      key <- names(params)[[i]]
      cat(sprintf("  [%d/%d] %s\n", i, length(params), param_info[[key]]$label))
      profiles[[i]] <- .profile_single_param(model, data, key, param_info,
                                             percent, n_points)
    }
  }

  # Combine results
  result <- do.call(rbind, profiles)

  # Add metadata
  attr(result, "true_loglik") <- true_ll
  attr(result, "n_params") <- length(params)
  attr(result, "n_trials") <- nrow(data)

  cat("\nProfile likelihood computation complete!\n")
  result
}

#' Plot profile likelihoods
#'
#' @param profile_result Result from profile_likelihood()
#' @param max_plots Maximum number of plots per page (default 9)
#' @export
plot_profile <- function(profile_result, max_plots = 9) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  true_ll <- attr(profile_result, "true_loglik")
  params <- unique(profile_result$parameter)
  n_params <- length(params)

  # Determine grid layout
  n_rows <- ceiling(sqrt(min(n_params, max_plots)))
  n_cols <- ceiling(min(n_params, max_plots) / n_rows)

  # Create plots
  n_pages <- ceiling(n_params / max_plots)

  for (page in seq_len(n_pages)) {
    start_idx <- (page - 1) * max_plots + 1
    end_idx <- min(page * max_plots, n_params)
    page_params <- params[start_idx:end_idx]

    plot_list <- lapply(page_params, function(param) {
      subset_data <- profile_result[profile_result$parameter == param, ]
      label <- subset_data$label[[1]]
      true_val <- subset_data$true_value[[1]]

      ggplot2::ggplot(subset_data, ggplot2::aes(x = value, y = loglik)) +
        ggplot2::geom_line(color = "blue", size = 1) +
        ggplot2::geom_point(color = "blue", size = 2) +
        ggplot2::geom_vline(xintercept = true_val,
                           linetype = "dashed", color = "red") +
        ggplot2::geom_hline(yintercept = true_ll,
                           linetype = "dotted", color = "gray50") +
        ggplot2::labs(
          title = label,
          x = "Parameter value",
          y = "Log-likelihood"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 10),
          axis.title = ggplot2::element_text(size = 9)
        )
    })

    # Arrange plots
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(grobs = plot_list,
                              nrow = n_rows, ncol = n_cols)
    } else {
      # Just print them sequentially
      for (p in plot_list) {
        print(p)
      }
    }

    if (page < n_pages) {
      readline("Press [Enter] for next page...")
    }
  }

  invisible(profile_result)
}
