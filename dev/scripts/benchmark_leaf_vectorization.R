#!/usr/bin/env Rscript

script_file <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) == 0L) {
    stop("run this script with Rscript dev/scripts/benchmark_leaf_vectorization.R")
  }
  sub("^--file=", "", file_arg[[1]])
}

repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

suppressPackageStartupMessages({
  if (identical(Sys.getenv("ACCUMULATR_BENCH_INSTALLED", "true"), "true")) {
    library(AccumulatR)
  } else {
    library(pkgload)
    load_all(repo_root, quiet = TRUE, helpers = FALSE)
  }
})

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- Sys.getenv(
  "ACCUMULATR_BENCH_OUT",
  file.path(out_dir, "benchmark_leaf_vectorization.csv")
)

n_trials <- as.integer(Sys.getenv("ACCUMULATR_BENCH_TRIALS", "400"))
n_rep <- as.integer(Sys.getenv("ACCUMULATR_BENCH_N_REP", "5"))
target_sec <- as.numeric(Sys.getenv("ACCUMULATR_BENCH_TARGET_SEC", "0.10"))
min_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BENCH_MIN_INNER_REPS", "1"))
max_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BENCH_MAX_INNER_REPS", "500"))
if (!is.finite(n_trials) || n_trials < 2L) n_trials <- 400L
if (!is.finite(n_rep) || n_rep < 1L) n_rep <- 5L
if (!is.finite(target_sec) || target_sec <= 0) target_sec <- 0.10
if (!is.finite(min_inner_reps) || min_inner_reps < 1L) min_inner_reps <- 1L
if (!is.finite(max_inner_reps) || max_inner_reps < min_inner_reps) {
  max_inner_reps <- min_inner_reps
}

clamp01 <- function(x) pmin(pmax(x, 0), 1)
safe_density <- function(x) ifelse(is.finite(x) & x > 0, x, 0)

lba_denom <- function(v, sv) {
  denom <- stats::pnorm(v / sv)
  if (!is.finite(denom) || denom < 1e-10) denom <- 1e-10
  denom
}

lba_pdf <- function(x, v, B, A, sv) {
  if (!is.finite(x) || x <= 0 || !is.finite(v) || !is.finite(B) ||
      !is.finite(A) || !is.finite(sv) || sv <= 0) {
    return(0)
  }
  denom <- lba_denom(v, sv)
  if (A > 1e-10) {
    zs <- x * sv
    cmz <- B - x * v
    cz <- cmz / zs
    cz_max <- (cmz - A) / zs
    out <- (v * (stats::pnorm(cz) - stats::pnorm(cz_max)) +
      sv * (stats::dnorm(cz_max) - stats::dnorm(cz))) / (A * denom)
  } else {
    out <- stats::dnorm(B / x, mean = v, sd = sv) * B / (x * x * denom)
  }
  safe_density(out)
}

lba_cdf <- function(x, v, B, A, sv) {
  if (!is.finite(x) || x <= 0 || !is.finite(v) || !is.finite(B) ||
      !is.finite(A) || !is.finite(sv) || sv <= 0) {
    return(0)
  }
  denom <- lba_denom(v, sv)
  if (A > 1e-10) {
    zs <- x * sv
    cmz <- B - x * v
    xx <- cmz - A
    cz <- cmz / zs
    cz_max <- xx / zs
    out <- (1 + (zs * (stats::dnorm(cz_max) - stats::dnorm(cz)) +
      xx * stats::pnorm(cz_max) -
      cmz * stats::pnorm(cz)) / A) / denom
  } else {
    out <- stats::pnorm(B / x, mean = v, sd = sv, lower.tail = FALSE) / denom
  }
  clamp01(out)
}

rdm_pigt0 <- function(x, k, l) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l)) return(0)
  if (abs(l) < 1e-12) {
    return(clamp01(2 * stats::pnorm(k / sqrt(x), lower.tail = FALSE)))
  }
  mu <- k / l
  lambda <- k * k
  p1 <- stats::pnorm(sqrt(lambda / x) * (1 + x / mu), lower.tail = FALSE)
  p2 <- stats::pnorm(sqrt(lambda / x) * (1 - x / mu), lower.tail = FALSE)
  clamp01(exp(exp(log(2 * lambda) - log(mu)) + log(max(1e-300, p1))) + p2)
}

rdm_digt0 <- function(x, k, l) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l)) return(0)
  lambda <- k * k
  exponent <- if (l == 0) {
    -0.5 * lambda / x
  } else {
    mu <- k / l
    -(lambda / (2 * x)) * ((x * x) / (mu * mu) - 2 * x / mu + 1)
  }
  exp(exponent + 0.5 * log(lambda) - 0.5 * log(2 * x * x * x * pi))
}

rdm_pigt <- function(x, k, l, a, threshold = 1e-10) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l) ||
      !is.finite(a)) {
    return(0)
  }
  if (a < threshold) return(rdm_pigt0(x, k, l))
  sqt <- sqrt(x)
  lgt <- log(x)
  if (l < threshold) {
    t5a <- 2 * stats::pnorm((k + a) / sqt) - 1
    t5b <- 2 * stats::pnorm((-k - a) / sqt) - 1
    t6a <- -0.5 * ((k + a)^2 / x - log(2) - log(pi) + lgt) - log(a)
    t6b <- -0.5 * ((k - a)^2 / x - log(2) - log(pi) + lgt) - log(a)
    out <- 1 + exp(t6a) - exp(t6b) +
      ((-k + a) * t5a - (k - a) * t5b) / (2 * a)
  } else {
    t1a <- exp(-0.5 * (k - a - x * l)^2 / x)
    t1b <- exp(-0.5 * (a + k - x * l)^2 / x)
    t1 <- exp(0.5 * (lgt - log(2) - log(pi))) * (t1a - t1b)
    t2a <- exp(2 * l * (k - a) +
      stats::pnorm(-(k - a + x * l) / sqt, log.p = TRUE))
    t2b <- exp(2 * l * (k + a) +
      stats::pnorm(-(k + a + x * l) / sqt, log.p = TRUE))
    t2 <- a + (t2b - t2a) / (2 * l)
    t4a <- 2 * stats::pnorm((k + a) / sqt - sqt * l) - 1
    t4b <- 2 * stats::pnorm((k - a) / sqt - sqt * l) - 1
    t4 <- 0.5 * (x * l - a - k + 0.5 / l) * t4a +
      0.5 * (k - a - x * l - 0.5 / l) * t4b
    out <- 0.5 * (t4 + t2 + t1) / a
  }
  clamp01(out)
}

rdm_digt <- function(x, k, l, a, threshold = 1e-10) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l) ||
      !is.finite(a)) {
    return(0)
  }
  if (a < threshold) return(safe_density(rdm_digt0(x, k, l)))
  if (l < threshold) {
    term <- exp(-(k - a)^2 / (2 * x)) - exp(-(k + a)^2 / (2 * x))
    out <- exp(-0.5 * (log(2) + log(pi) + log(x)) +
      log(max(1e-300, term)) - log(2) - log(a))
  } else {
    sqt <- sqrt(x)
    t1a <- -(a - k + x * l)^2 / (2 * x)
    t1b <- -(a + k - x * l)^2 / (2 * x)
    t1 <- (1 / sqrt(2)) * (exp(t1a) - exp(t1b)) / (sqrt(pi) * sqt)
    t2a <- 2 * stats::pnorm((-k + a) / sqt + sqt * l) - 1
    t2b <- 2 * stats::pnorm((k + a) / sqt - sqt * l) - 1
    t2 <- exp(log(0.5) + log(l)) * (t2a + t2b)
    out <- exp(log(max(1e-300, t1 + t2)) - log(2) - log(a))
  }
  safe_density(out)
}

exgauss_raw_pdf <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || !is.finite(mu) || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(tau) || tau <= 0) {
    return(0)
  }
  inv_tau <- 1 / tau
  z <- (x - mu) / sigma
  exponent <- sigma * sigma * inv_tau * inv_tau * 0.5 - (x - mu) * inv_tau
  out <- inv_tau * exp(exponent) * stats::pnorm(z - sigma * inv_tau)
  safe_density(out)
}

exgauss_raw_cdf <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || !is.finite(mu) || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(tau) || tau <= 0) {
    return(0)
  }
  inv_tau <- 1 / tau
  z <- (x - mu) / sigma
  exponent <- sigma * sigma * inv_tau * inv_tau * 0.5 - (x - mu) * inv_tau
  clamp01(stats::pnorm(z) -
    exp(exponent) * stats::pnorm(z - sigma * inv_tau))
}

exgauss_pdf <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || x <= 0) return(0)
  lower_cdf <- exgauss_raw_cdf(0, mu, sigma, tau)
  lower_survival <- 1 - lower_cdf
  if (!is.finite(lower_survival) || lower_survival <= 0) return(0)
  exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival
}

exgauss_cdf <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || x <= 0) return(0)
  lower_cdf <- exgauss_raw_cdf(0, mu, sigma, tau)
  lower_survival <- 1 - lower_cdf
  if (!is.finite(lower_survival) || lower_survival <= 0) return(0)
  clamp01((exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) / lower_survival)
}

leaf_base <- function(dist, x, params, channel) {
  if (!is.finite(x) || x <= 0) {
    return(if (identical(channel, "survival")) 1 else 0)
  }
  q <- params[["q"]]
  if (identical(dist, "gamma")) {
    pdf <- stats::dgamma(x, shape = params[["shape"]], rate = params[["rate"]])
    cdf <- stats::pgamma(x, shape = params[["shape"]], rate = params[["rate"]])
  } else if (identical(dist, "exgauss")) {
    pdf <- exgauss_pdf(x, params[["mu"]], params[["sigma"]], params[["tau"]])
    cdf <- exgauss_cdf(x, params[["mu"]], params[["sigma"]], params[["tau"]])
  } else if (identical(dist, "LBA")) {
    pdf <- lba_pdf(x, params[["v"]], params[["B"]], params[["A"]], params[["sv"]])
    cdf <- lba_cdf(x, params[["v"]], params[["B"]], params[["A"]], params[["sv"]])
  } else if (identical(dist, "RDM")) {
    v_sc <- params[["v"]] / params[["s"]]
    B_sc <- params[["B"]] / params[["s"]]
    A_sc <- params[["A"]] / params[["s"]]
    a <- 0.5 * A_sc
    pdf <- rdm_digt(x, B_sc + a, v_sc, a)
    cdf <- rdm_pigt(x, B_sc + a, v_sc, a)
  } else {
    pdf <- stats::dlnorm(x, meanlog = params[["m"]], sdlog = params[["s"]])
    cdf <- stats::plnorm(x, meanlog = params[["m"]], sdlog = params[["s"]])
  }
  if (identical(channel, "pdf")) {
    return((1 - q) * safe_density(pdf))
  }
  cdf <- clamp01((1 - q) * clamp01(cdf))
  if (identical(channel, "cdf")) cdf else clamp01(1 - cdf)
}

make_case <- function(label, dist, a_params) {
  structure <- race_spec() |>
    add_accumulator("a", dist) |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    finalize_model()
  trial_df <- data.frame(
    trial = seq_len(n_trials),
    R = rep(c("A", "B"), length.out = n_trials),
    rt = seq(0.24, 0.92, length.out = n_trials),
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, trial_df)
  b_params <- c(m = log(0.54), s = 0.20, q = 0.03, t0 = 0.015)
  params <- c(setNames(a_params, paste0("a.", names(a_params))),
              setNames(b_params, paste0("b.", names(b_params))))
  params_mat <- build_param_matrix(structure, params, trial_df = prepared)
  ctx <- make_context(structure)
  eval_fn <- function() log_likelihood(ctx, prepared, params_mat)
  engine <- as.numeric(eval_fn())
  a_ref <- as.list(a_params)
  b_ref <- as.list(b_params)
  manual_terms <- vapply(seq_len(nrow(trial_df)), function(i) {
    rt <- trial_df$rt[[i]]
    if (identical(trial_df$R[[i]], "A")) {
      p <- leaf_base(dist, rt - a_ref$t0, a_ref, "pdf") *
        leaf_base("lognormal", rt - b_ref$t0, b_ref, "survival")
    } else {
      p <- leaf_base("lognormal", rt - b_ref$t0, b_ref, "pdf") *
        leaf_base(dist, rt - a_ref$t0, a_ref, "survival")
    }
    log(max(1e-10, p))
  }, numeric(1))
  list(
    label = label,
    eval_fn = eval_fn,
    engine = engine,
    manual = sum(manual_terms),
    abs_error = abs(engine - sum(manual_terms))
  )
}

choose_inner_reps <- function(fn) {
  calib_reps <- 1L
  elapsed <- 0
  while ((!is.finite(elapsed) || elapsed <= 0) && calib_reps <= 4096L) {
    elapsed <- system.time({
      for (i in seq_len(calib_reps)) fn()
    })[["elapsed"]]
    if (!is.finite(elapsed) || elapsed <= 0) {
      calib_reps <- calib_reps * 2L
    }
  }
  if (!is.finite(elapsed) || elapsed <= 0) return(min_inner_reps)
  per_eval <- elapsed / calib_reps
  if (!is.finite(per_eval) || per_eval <= 0) return(min_inner_reps)
  as.integer(max(min_inner_reps, min(max_inner_reps, ceiling(target_sec / per_eval))))
}

run_case <- function(case) {
  inner_reps <- choose_inner_reps(case$eval_fn)
  timings <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    gc(FALSE)
    timings[[i]] <- system.time({
      for (j in seq_len(inner_reps)) case$eval_fn()
    })[["elapsed"]]
  }
  data.frame(
    label = case$label,
    n_trials = n_trials,
    n_rep = n_rep,
    inner_reps = inner_reps,
    median_per_eval_ms = 1000 * stats::median(timings) / inner_reps,
    min_per_eval_ms = 1000 * min(timings) / inner_reps,
    max_per_eval_ms = 1000 * max(timings) / inner_reps,
    engine_loglik = case$engine,
    manual_loglik = case$manual,
    abs_loglik_error = case$abs_error,
    stringsAsFactors = FALSE
  )
}

cases <- list(
  make_case(
    "gamma_pdf_and_survival",
    "gamma",
    c(shape = 3.0, rate = 8.0, q = 0.04, t0 = 0.02)
  ),
  make_case(
    "exgauss_pdf_and_survival",
    "exgauss",
    c(mu = 0.28, sigma = 0.07, tau = 0.12, q = 0.03, t0 = 0.03)
  ),
  make_case(
    "lba_pdf_and_survival",
    "LBA",
    c(v = 2.0, B = 1.2, A = 0.4, sv = 0.6, q = 0.02, t0 = 0.04)
  ),
  make_case(
    "rdm_pdf_and_survival",
    "RDM",
    c(v = 1.5, B = 1.0, A = 0.3, s = 1.1, q = 0.02, t0 = 0.04)
  )
)

out <- do.call(rbind, lapply(cases, run_case))
utils::write.csv(out, out_file, row.names = FALSE)
print(out, row.names = FALSE)
cat("Wrote leaf vectorization benchmark to", out_file, "\n")
