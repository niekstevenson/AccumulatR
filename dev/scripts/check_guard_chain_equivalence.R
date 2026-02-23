#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(AccumulatR))

`%||%` <- function(x, y) if (is.null(x)) y else x

run_with_guard_mode <- function(mode, fn) {
  old <- Sys.getenv("ACCUMULATR_GUARD_CDF_MODE", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("ACCUMULATR_GUARD_CDF_MODE")
    } else {
      Sys.setenv(ACCUMULATR_GUARD_CDF_MODE = old)
    }
  }, add = TRUE)
  if (identical(mode, "auto")) {
    Sys.unsetenv("ACCUMULATR_GUARD_CDF_MODE")
  } else {
    Sys.setenv(ACCUMULATR_GUARD_CDF_MODE = mode)
  }
  fn()
}

normalize_probs <- function(p) {
  vals <- as.numeric(p)
  nms <- names(p)
  if (is.null(nms)) {
    nms <- rep("", length(vals))
  }
  nms[is.na(nms)] <- "NA"
  names(vals) <- nms
  vals
}

max_named_abs_diff <- function(a, b) {
  keys <- sort(unique(c(names(a), names(b))))
  av <- setNames(rep(0.0, length(keys)), keys)
  bv <- av
  av[names(a)] <- a
  bv[names(b)] <- b
  max(abs(av - bv))
}

make_chain_expr <- function(ids) {
  if (length(ids) == 1L) {
    return(ids[[1L]])
  }
  inhibit(ids[[1L]], by = make_chain_expr(ids[-1L]))
}

build_chain_case <- function(depth, rep_idx) {
  ids <- sprintf("a%02d", seq_len(depth + 1L))
  spec_obj <- race_spec()
  for (id in ids) {
    spec_obj <- add_accumulator(spec_obj, id, "lognormal")
  }
  spec_obj <- add_accumulator(spec_obj, "comp", "lognormal")
  spec_obj <- add_outcome(spec_obj, "CHAIN", make_chain_expr(ids))
  spec_obj <- add_outcome(spec_obj, "COMP", "comp")
  structure <- finalize_model(spec_obj)

  all_ids <- c(ids, "comp")
  params <- numeric(2L * length(all_ids))
  names(params) <- as.vector(rbind(
    paste0(all_ids, ".meanlog"),
    paste0(all_ids, ".sdlog")
  ))
  for (i in seq_along(all_ids)) {
    mean_val <- 0.22 + 0.028 * i + 0.003 * rep_idx
    sd_val <- 0.11 + 0.004 * ((i + rep_idx) %% 7L)
    params[paste0(all_ids[[i]], ".meanlog")] <- log(mean_val)
    params[paste0(all_ids[[i]], ".sdlog")] <- sd_val
  }

  list(
    depth = depth,
    rep_idx = rep_idx,
    spec_obj = spec_obj,
    structure = structure,
    params = params
  )
}

depth_min <- as.integer(Sys.getenv("ACCUMULATR_EQ_DEPTH_MIN", "2"))
depth_max <- as.integer(Sys.getenv("ACCUMULATR_EQ_DEPTH_MAX", "8"))
reps_per_depth <- as.integer(Sys.getenv("ACCUMULATR_EQ_REPS", "4"))
n_trials <- as.integer(Sys.getenv("ACCUMULATR_EQ_TRIALS", "500"))
run_ll <- Sys.getenv("ACCUMULATR_EQ_RUN_LL", "1") != "0"
run_prob <- Sys.getenv("ACCUMULATR_EQ_RUN_PROB", "1") != "0"

depth_min <- if (is.na(depth_min)) 2L else depth_min
depth_max <- if (is.na(depth_max)) 8L else depth_max
if (depth_max < depth_min) {
  stop("depth_max must be >= depth_min")
}
reps_per_depth <- if (is.na(reps_per_depth) || reps_per_depth < 1L) 1L else reps_per_depth
n_trials <- if (is.na(n_trials) || n_trials < 1L) 1L else n_trials
depths <- depth_min:depth_max

cat(sprintf(
  "Running guard equivalence audit depths=%d..%d reps=%d n_trials=%d run_prob=%s run_ll=%s\n",
  depth_min, depth_max, reps_per_depth, n_trials,
  if (run_prob) "yes" else "no",
  if (run_ll) "yes" else "no"
))

rows <- vector("list", length(depths) * reps_per_depth)
row_i <- 0L

for (depth in depths) {
  for (rep_idx in seq_len(reps_per_depth)) {
    row_i <- row_i + 1L
    cs <- build_chain_case(depth, rep_idx)

    prob_diff <- NA_real_
    if (run_prob) {
      prob_params <- build_param_matrix(cs$spec_obj, cs$params, n_trials = 1L)
      p_auto <- run_with_guard_mode("auto", function() {
        normalize_probs(response_probabilities(cs$structure, prob_params, include_na = TRUE))
      })
      p_quad <- run_with_guard_mode("quadrature", function() {
        normalize_probs(response_probabilities(cs$structure, prob_params, include_na = TRUE))
      })
      prob_diff <- max_named_abs_diff(p_auto, p_quad)
    }

    ll_auto <- NA_real_
    ll_quad <- NA_real_
    ll_diff <- NA_real_
    if (run_ll) {
      sim_params <- build_param_matrix(cs$spec_obj, cs$params, n_trials = n_trials)
      sim_seed <- 10000L + depth * 100L + rep_idx
      data_df <- simulate(cs$structure, sim_params, seed = sim_seed, keep_component = TRUE)
      ctx <- build_likelihood_context(cs$structure, data_df)
      slim_params <- build_param_matrix(
        cs$spec_obj,
        cs$params,
        n_trials = max(data_df$trial),
        layout = ctx$param_layout
      )

      ll_auto <- run_with_guard_mode("auto", function() {
        as.numeric(log_likelihood(ctx, slim_params))
      })
      ll_quad <- run_with_guard_mode("quadrature", function() {
        as.numeric(log_likelihood(ctx, slim_params))
      })
      ll_diff <- abs(ll_auto - ll_quad)
    }

    rows[[row_i]] <- data.frame(
      depth = depth,
      rep = rep_idx,
      prob_diff = prob_diff,
      ll_diff = ll_diff,
      ll_auto = ll_auto,
      ll_quad = ll_quad,
      stringsAsFactors = FALSE
    )
    cat(sprintf(
      "depth=%d rep=%d prob_diff=%.9g ll_diff=%.9g\n",
      depth, rep_idx, prob_diff, ll_diff
    ))
    flush.console()
  }
}

res <- do.call(rbind, rows)
summary_by_depth <- do.call(
  rbind,
  lapply(split(res, res$depth), function(df) {
    data.frame(
      depth = df$depth[[1L]],
      prob_diff = if (all(is.na(df$prob_diff))) NA_real_ else max(df$prob_diff, na.rm = TRUE),
      ll_diff = if (all(is.na(df$ll_diff))) NA_real_ else max(df$ll_diff, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

cat("\nMax diffs by depth:\n")
print(summary_by_depth, row.names = FALSE)

global_prob <- if (all(is.na(res$prob_diff))) NA_real_ else max(res$prob_diff %||% 0, na.rm = TRUE)
global_ll <- if (all(is.na(res$ll_diff))) NA_real_ else max(res$ll_diff %||% 0, na.rm = TRUE)
cat(sprintf(
  "\nGlobal max prob diff: %.9g\nGlobal max ll diff: %.9g\n",
  global_prob, global_ll
))

out_dir <- file.path("dev", "scripts", "scratch_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "guard_chain_equivalence_depth2_10.csv")
utils::write.csv(res, out_file, row.names = FALSE)
cat(sprintf("Wrote results to %s\n", out_file))
