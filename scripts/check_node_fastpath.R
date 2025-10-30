#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
example_name <- if (length(args) >= 1) args[[1]] else "example_12_inhibitor_with_protector"

message("Loading likelihood stack…")
source("R/generator_new.R")
source("examples/new_API.R")
source("R/super_large_likelihood.R")

if (!example_name %in% names(new_api_examples)) {
  stop(sprintf("Unknown example '%s'. Available names: %s",
               example_name,
               paste(sort(names(new_api_examples)), collapse = ", ")))
}

model <- new_api_examples[[example_name]]
prep <- .prepare_model_for_likelihood(model)
compiled <- .prep_expr_compiled(prep)
nodes <- compiled$nodes %||% list()

if (length(nodes) == 0) {
  stop("Compiled node table is empty; nothing to inspect.")
}

# choose a handful of evaluation times
t_grid <- c(0.1, 0.2, 0.4, 0.6)

summaries <- lapply(seq_along(nodes), function(idx) {
  node <- nodes[[idx]]
  out <- data.frame(
    node_id = idx,
    kind = node$kind %||% NA_character_,
    has_fast_density = is.function(node$density_fast_fn),
    has_fast_surv = is.function(node$surv_fast_fn),
    has_fast_cdf = is.function(node$cdf_fast_fn),
    uses_fast_density = NA,
    uses_fast_surv = NA,
    max_scenarios = NA_integer_,
    stringsAsFactors = FALSE
  )
  dens_fast_flag <- FALSE
  surv_fast_flag <- FALSE
  scenario_max <- 0L
  for (tt in t_grid) {
    dens_val <- .node_density(
      node, tt, prep, component = "__default__",
      forced_complete = integer(0),
      forced_survive = integer(0),
      cache = NULL,
      prefer_fast = TRUE
    )
    dens_fast_flag <- dens_fast_flag || isTRUE(attr(dens_val, "fast_used"))
    surv_val <- .node_survival_cond(
      node, tt, prep, component = "__default__",
      forced_complete = integer(0),
      forced_survive = integer(0),
      cache = NULL,
      prefer_fast = TRUE
    )
    surv_fast_flag <- surv_fast_flag || isTRUE(attr(surv_val, "fast_used"))
    scen <- .node_scenarios_at(
      node, tt, prep, component = "__default__",
      forced_complete = integer(0),
      forced_survive = integer(0),
      cache = NULL
    )
    scenario_max <- max(scenario_max, length(scen))
  }
  out$uses_fast_density <- dens_fast_flag
  out$uses_fast_surv <- surv_fast_flag
  out$max_scenarios <- scenario_max
  out
})

summary_df <- do.call(rbind, summaries)

message(sprintf("Summary for %s:", example_name))
print(summary_df)

message("\nLegend:")
message("  has_fast_*      -> node compiled a fast closure")
message("  uses_fast_*     -> at least one t in grid triggered the fast closure")
message("  max_scenarios   -> maximum number of scenario records returned over the grid")
