.test_model_complete_parameters <- function(structure, param_values) {
  required <- sampled_pars(structure)
  extras <- setdiff(names(param_values), required)
  if (length(extras) > 0L) {
    stop(
      "Unknown parameter values supplied: ",
      paste(extras, collapse = ", "),
      call. = FALSE
    )
  }

  missing <- setdiff(required, names(param_values))
  if (length(missing) > 0L) {
    warning(
      "Assuming 0 for missing parameter(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
    param_values <- c(param_values, stats::setNames(rep(0, length(missing)), missing))
  }

  param_values[required]
}

.test_model_domain <- function(structure, external_name) {
  lookup <- .parameter_name_lookup(race_model(structure))
  targets <- names(lookup)[unname(lookup) == external_name]
  suffixes <- sub("^.*\\.", "", targets)

  positive_suffixes <- c("s", "sdlog", "sigma", "tau", "shape", "rate", "B", "A", "sv")
  nonnegative_suffixes <- c("q", "t0")

  if (length(suffixes) > 0L && all(suffixes %in% positive_suffixes)) {
    return("positive")
  }
  if (length(suffixes) > 0L && all(suffixes %in% nonnegative_suffixes)) {
    return("nonnegative")
  }
  "unrestricted"
}

.test_model_profile_grid <- function(center,
                                     domain,
                                     points = 21L,
                                     span = 0.5) {
  points <- as.integer(points)
  if (!is.finite(center)) {
    return(rep(center, points))
  }

  if (identical(domain, "positive")) {
    if (center > 0) {
      lower <- max(center * (1 - span), .Machine$double.eps)
      upper <- center * (1 + span)
    } else {
      lower <- .Machine$double.eps
      upper <- max(span, 0.25)
    }
    return(seq(lower, upper, length.out = points))
  }

  if (identical(domain, "nonnegative")) {
    if (center > 0) {
      lower <- max(0, center * (1 - span))
      upper <- center * (1 + span)
    } else {
      lower <- 0
      upper <- max(span, 0.25)
    }
    return(seq(lower, upper, length.out = points))
  }

  width <- max(abs(center) * span, span)
  seq(center - width, center + width, length.out = points)
}

.test_model_simulated_probabilities <- function(sim_df, labels) {
  sim_labels <- as.character(sim_df$R)
  sim_labels[is.na(sim_labels)] <- "NA"
  probs <- prop.table(table(sim_labels))
  out <- stats::setNames(numeric(length(labels)), labels)
  shared <- intersect(labels, names(probs))
  out[shared] <- as.numeric(probs[shared])
  out
}

.test_model_comparison_table <- function(analytical, simulated) {
  labels <- names(analytical)
  comparison <- data.frame(
    response = labels,
    analytical = round(as.numeric(analytical[labels]), 3),
    simulated = round(as.numeric(simulated[labels]), 3),
    abs_diff = round(abs(as.numeric(analytical[labels]) - as.numeric(simulated[labels])), 3),
    stringsAsFactors = FALSE
  )

  keep <- !(comparison$analytical == 0 & comparison$simulated == 0)
  comparison[keep, , drop = FALSE]
}

.test_model_profile_one <- function(structure,
                                    context,
                                    data,
                                    baseline_values,
                                    parameter,
                                    points,
                                    span) {
  center <- as.numeric(baseline_values[[parameter]])[1]
  grid <- .test_model_profile_grid(
    center = center,
    domain = .test_model_domain(structure, parameter),
    points = points,
    span = span
  )

  ll <- vapply(grid, function(value) {
    candidate <- baseline_values
    candidate[[parameter]] <- value
    n_trials <- if (nrow(data) == 0L) {
      0L
    } else {
      1L + sum(data$trial[-1L] != data$trial[-nrow(data)])
    }
    pm <- tryCatch(
      build_param_matrix(
        structure,
        candidate,
        n_trials = n_trials
      ),
      error = function(e) NULL
    )
    if (is.null(pm)) {
      return(NA_real_)
    }
    out <- tryCatch(
      as.numeric(log_likelihood(context, data, pm)),
      error = function(e) NA_real_
    )
    out[[1]]
  }, numeric(1))

  data.frame(
    parameter = parameter,
    value = grid,
    log_likelihood = ll,
    stringsAsFactors = FALSE
  )
}

.plot_test_model_profiles <- function(profiles, baseline_values) {
  if (length(profiles) == 0L) {
    return(invisible(NULL))
  }

  n <- length(profiles)
  n_col <- ceiling(sqrt(n))
  n_row <- ceiling(n / n_col)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(n_row, n_col))

  for (nm in names(profiles)) {
    prof <- profiles[[nm]]
    ok <- is.finite(prof$log_likelihood)
    if (!any(ok)) {
      graphics::plot.new()
      graphics::title(main = nm)
      graphics::text(0.5, 0.5, "No finite likelihood values")
      next
    }

    graphics::plot(
      prof$value[ok],
      prof$log_likelihood[ok],
      type = "l",
      lwd = 2,
      xlab = nm,
      ylab = "Log-likelihood",
      main = nm
    )
    graphics::abline(v = baseline_values[[nm]], col = "firebrick", lty = 2)
    base_idx <- which.min(abs(prof$value - baseline_values[[nm]]))
    if (length(base_idx) == 1L && is.finite(prof$log_likelihood[[base_idx]])) {
      graphics::points(
        prof$value[[base_idx]],
        prof$log_likelihood[[base_idx]],
        pch = 19,
        col = "firebrick"
      )
    }
  }

  invisible(NULL)
}

#' Simulate from a model and inspect whether it behaves as expected
#'
#' `test_model()` is a quick diagnostic for checking whether a model and a
#' parameter vector behave sensibly. It does three things:
#'
#' 1. fills in any missing parameters with `0` and warns which ones were
#'    assumed to be zero
#' 2. compares analytical response probabilities with response proportions from
#'    simulated behavioral data
#' 3. profiles the log-likelihood for each parameter you supplied, holding the
#'    others fixed
#'
#' This is especially useful when you are building a new model, checking a
#' custom `set_parameters()` mapping, or verifying that a parameter vector gives
#' sensible behavior before fitting real behavioral data.
#'
#' @param model Model specification or finalized model.
#' @param param_values Named numeric vector of parameter values. Parameters not
#'   supplied are assumed to be `0`, with a warning. Names should follow
#'   `sampled_pars(model)`, so custom names from `set_parameters()` are
#'   supported.
#' @param n_trials Number of trials to simulate.
#' @param seed Random seed used for simulation.
#' @param include_na If `TRUE`, include missing-response probability mass in the
#'   analytical and simulated comparison.
#' @param profile_points Number of grid points per parameter profile.
#' @param profile_span Relative profiling span. Positive and non-negative
#'   parameters are profiled on a multiplicative scale around their supplied
#'   value; unrestricted parameters are profiled on a symmetric additive scale.
#' @param plot If `TRUE`, draw one profile plot per supplied parameter.
#' @return Invisibly returns a list with the completed parameter vector, the
#'   analytical and simulated probability comparison, the simulated data, the
#'   prepared data, the model context, and the profile tables.
#' @examples
#' spec <- race_spec() |>
#'   add_accumulator("go", "lognormal") |>
#'   add_accumulator("stop", "lognormal") |>
#'   add_outcome("go", "go") |>
#'   add_outcome("stop", "stop") |>
#'   set_parameters(list(
#'     drift = c("go.m", "stop.m")
#'   ))
#'
#' res <- test_model(
#'   spec,
#'   c(drift = log(0.30), go.s = 0.16, stop.s = 0.18),
#'   n_trials = 50,
#'   plot = FALSE
#' )
#'
#' res$comparison
#' @export
test_model <- function(model,
                       param_values,
                       n_trials,
                       seed = 123,
                       include_na = TRUE,
                       profile_points = 21L,
                       profile_span = 0.5,
                       plot = TRUE) {
  structure <- .as_model_structure(model)

  if (missing(param_values) || is.null(param_values) || length(param_values) == 0L) {
    stop("test_model() requires a named numeric vector of parameter values", call. = FALSE)
  }
  if (is.null(names(param_values)) || any(!nzchar(names(param_values)))) {
    stop("test_model() requires a named numeric vector of parameter values", call. = FALSE)
  }
  if (anyDuplicated(names(param_values))) {
    stop("Parameter names supplied to test_model() must be unique", call. = FALSE)
  }
  if (!is.numeric(param_values)) {
    stop("test_model() requires numeric parameter values", call. = FALSE)
  }
  if (!is.numeric(n_trials) || length(n_trials) != 1L || !is.finite(n_trials) ||
      n_trials < 1L || !isTRUE(all.equal(n_trials, round(n_trials)))) {
    stop("n_trials must be a single positive integer", call. = FALSE)
  }
  if (!is.numeric(profile_points) || length(profile_points) != 1L ||
      profile_points < 3L || !isTRUE(all.equal(profile_points, round(profile_points)))) {
    stop("profile_points must be an integer >= 3", call. = FALSE)
  }
  if (!is.numeric(profile_span) || length(profile_span) != 1L ||
      !is.finite(profile_span) || profile_span <= 0) {
    stop("profile_span must be a single positive numeric value", call. = FALSE)
  }

  supplied_names <- names(param_values)
  completed_values <- .test_model_complete_parameters(structure, param_values)

  params_df_sim <- build_param_matrix(
    structure,
    completed_values,
    n_trials = as.integer(n_trials)
  )
  analytical <- response_probabilities(
    structure,
    build_param_matrix(structure, completed_values, n_trials = 1L),
    include_na = include_na
  )

  sim_df <- simulate(structure, params_df_sim, seed = seed)
  prepared_data <- prepare_data(structure, sim_df)
  context <- make_context(structure)

  labels <- names(analytical)
  simulated <- .test_model_simulated_probabilities(sim_df, labels)
  comparison <- .test_model_comparison_table(analytical, simulated)

  print(comparison, row.names = FALSE)

  profile_names <- intersect(supplied_names, names(completed_values))
  profiles <- lapply(profile_names, function(nm) {
    .test_model_profile_one(
      structure = structure,
      context = context,
      data = prepared_data,
      baseline_values = completed_values,
      parameter = nm,
      points = as.integer(profile_points),
      span = profile_span
    )
  })
  names(profiles) <- profile_names

  if (isTRUE(plot) && length(profiles) > 0L) {
    .plot_test_model_profiles(profiles, completed_values[profile_names])
  }

  invisible(list(
    parameters = completed_values,
    comparison = comparison,
    simulated_data = sim_df,
    context = context,
    profiles = profiles
  ))
}
