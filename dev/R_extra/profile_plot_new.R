`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.profile_param_curve <- function(structure,
                                 model_spec,
                                 base_params,
                                 data,
                                 param_name,
                                 label,
                                 percent,
                                 n_points,
                                 n_trials,
                                 ctx) {
  base_value <- base_params[[param_name]]
  span <- if (!is.finite(base_value) || base_value == 0) {
    percent
  } else {
    abs(base_value) * percent
  }
  grid <- seq(base_value - span, base_value + span, length.out = n_points)
  param_tables <- lapply(grid, function(val) {
    candidate <- base_params
    candidate[[param_name]] <- val
    build_params_df(model_spec, candidate, n_trials = n_trials)
  })
  log_liks <- native_loglikelihood_param_repeat(ctx, param_tables)
  data.frame(
    parameter = param_name,
    label = label,
    value = grid,
    true_value = base_value,
    loglik = log_liks,
    stringsAsFactors = FALSE
  )
}

#' Compute profile likelihoods for all parameters in a model
#'
#' This version operates directly on the new param-table API:
#' you supply the generator structure, the model specification,
#' and the named vector of core parameter values that are passed to
#' `build_params_df()`. Each grid point rebuilds the params table and
#' evaluates `log_likelihood_from_params()`.
#'
#' @param structure Generator structure from `build_generator_structure()`
#' @param model_spec Model definition (race_spec or race_model_spec)
#' @param base_params Named numeric vector of core parameter values
#' @param data Data frame of simulated trials
#' @param percent Fractional range +/- around each parameter value
#' @param n_points Number of grid points per parameter
#' @param n_cores Optional parallel core count
#' @param param_labels Optional named character vector for plotting labels
#' @export
profile_likelihood <- function(structure,
                               model_spec,
                               base_params,
                               data,
                               percent = 0.1,
                               n_points = 10,
                               n_cores = 1,
                               param_labels = NULL) {
  if (is.null(names(base_params))) {
    stop("base_params must be a named numeric vector")
  }
  n_trials <- length(unique(data$trial %||% seq_len(nrow(data))))
  base_table <- build_params_df(model_spec, base_params, n_trials = n_trials)
  ctx <- build_likelihood_context(
    structure = structure,
    params_df = base_table,
    data_df = data
  )
  true_ll <- native_loglikelihood_param_repeat(ctx, list(base_table))[[1]]
  label_map <- param_labels %||% stats::setNames(names(base_params), names(base_params))
  param_names <- names(base_params)
  compute_curve <- function(param_name) {
    .profile_param_curve(
      structure = structure,
      model_spec = model_spec,
      base_params = base_params,
      data = data,
      param_name = param_name,
      label = label_map[[param_name]] %||% param_name,
      percent = percent,
      n_points = n_points,
      n_trials = n_trials,
      ctx = ctx
    )
  }
  if (length(param_names) == 0) {
    stop("No parameters found in base_params")
  }
  if (n_cores > 1) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      profiles <- parallel::parLapply(cl, param_names, compute_curve)
    } else {
      profiles <- parallel::mclapply(param_names, compute_curve, mc.cores = n_cores)
    }
  } else {
    profiles <- lapply(param_names, compute_curve)
  }
  result <- do.call(rbind, profiles)
  attr(result, "true_loglik") <- true_ll
  attr(result, "n_params") <- length(param_names)
  attr(result, "n_trials") <- nrow(data)
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
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package required for plotting")
  }
  true_ll <- attr(profile_result, "true_loglik")
  params <- unique(profile_result$parameter)
  n_params <- length(params)
  n_rows <- ceiling(sqrt(min(n_params, max_plots)))
  n_cols <- ceiling(min(n_params, max_plots) / n_rows)
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
        ggplot2::geom_line(color = "blue", linewidth = 1) +
        ggplot2::geom_point(color = "blue", size = 2) +
        ggplot2::geom_vline(xintercept = true_val, linetype = "dashed", color = "red") +
        ggplot2::geom_hline(yintercept = true_ll, linetype = "dotted", color = "gray50") +
        ggplot2::labs(
          title = label,
          x = "Parameter value",
          y = "Log-likelihood"
        ) +
        ggplot2::theme_minimal()
    })
    grid <- gridExtra::arrangeGrob(grobs = plot_list, nrow = n_rows, ncol = n_cols)
    grid::grid.newpage()
    grid::grid.draw(grid)
  }
  invisible(profile_result)
}
