suppressPackageStartupMessages({
  library(pkgload)
  loaded <- FALSE
  try({
    load_all(".", compile = FALSE, quiet = TRUE)
    loaded <- TRUE
  }, silent = TRUE)
  if (!loaded) {
    load_all(".", compile = TRUE, quiet = TRUE)
  }
})

stop_selective_model <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("S1", "lognormal") |>
    add_accumulator("IS", "lognormal") |>
    add_accumulator("S2", "lognormal") |>
    add_outcome(
      "A",
      first_of(
        inhibit("A", by = "S1"),
        all_of("A", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome(
      "B",
      first_of(
        inhibit("B", by = "S1"),
        all_of("B", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome(
      "STOP",
      all_of("S1", inhibit("S2", by = "IS")),
      options = list(map_outcome_to = NA_character_)
    ) |>
    add_component("go", members = c("A", "B")) |>
    add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
    set_parameters(list(
      m_go = c("A.m", "B.m"),
      s_go = c("A.s", "B.s"),
      t0_go = c("A.t0", "B.t0")
    )) |>
    finalize_model()
}

stop_selective_params <- c(
  m_go = log(0.30), s_go = 0.18, t0_go = 0.05,
  S1.m = log(0.26), S1.s = 0.18, S1.q = 0.00, S1.t0 = 0.00,
  IS.m = log(0.35), IS.s = 0.18, IS.q = 0.00, IS.t0 = 0.00,
  S2.m = log(0.32), S2.s = 0.18, S2.q = 0.00, S2.t0 = 0.00
)

go_parts <- function(params) {
  list(
    m = unname(params[["m_go"]]),
    s = unname(params[["s_go"]]),
    q = 0.0,
    t0 = unname(params[["t0_go"]])
  )
}

acc_parts <- function(prefix, params) {
  list(
    m = unname(params[[paste0(prefix, ".m")]]),
    s = unname(params[[paste0(prefix, ".s")]]),
    q = unname(params[[paste0(prefix, ".q")]]),
    t0 = unname(params[[paste0(prefix, ".t0")]])
  )
}

acc_pdf <- function(t, p) {
  out <- numeric(length(t))
  ok <- is.finite(t) & t >= p$t0
  out[ok] <- (1.0 - p$q) * dlnorm(t[ok] - p$t0, meanlog = p$m, sdlog = p$s)
  out[!is.finite(out) | out < 0.0] <- 0.0
  out
}

acc_cdf <- function(t, p) {
  out <- numeric(length(t))
  out[is.infinite(t) & t > 0.0] <- 1.0 - p$q
  ok <- is.finite(t) & t >= p$t0
  out[ok] <- (1.0 - p$q) * plnorm(t[ok] - p$t0, meanlog = p$m, sdlog = p$s)
  pmin(pmax(out, 0.0), 1.0)
}

acc_survival <- function(t, p) {
  1.0 - acc_cdf(t, p)
}

integrate_vec <- function(fn,
                          lower,
                          upper,
                          rel_tol = 1e-10,
                          abs_tol = 1e-12,
                          subdivisions = 400L) {
  if (is.finite(upper) && upper <= lower) {
    return(0.0)
  }
  integrate(
    function(x) {
      y <- fn(x)
      y[!is.finite(y)] <- 0.0
      y
    },
    lower = lower,
    upper = upper,
    rel.tol = rel_tol,
    abs.tol = abs_tol,
    subdivisions = subdivisions,
    stop.on.error = FALSE
  )$value
}

ignore_pdf <- function(t, is_par, s2_par) {
  acc_pdf(t, is_par) * acc_survival(t, s2_par)
}

ignore_cdf <- function(t, is_par, s2_par) {
  if (!is.finite(t) || t <= 0.0) {
    return(0.0)
  }
  integrate_vec(function(x) ignore_pdf(x, is_par, s2_par), 0.0, t)
}

# Validated stimulus-selective stop reduction:
# - before S1, the stop component behaves like a standard A/B go race
# - after S1, the ignore branch is represented by a latent process H with
#   density h(t) = f_IS(t) * S_S2(t)
# - the candidate A density is
#   f_A(t) * S_B(t) * [S_S1(t) + F_S1(t) * H(t)] +
#   h(t) * ∫_0^t f_A(a) * S_B(a) * F_S1(a) da
# - B is symmetric
#
manual_beest_like_density <- function(component, response, rt, params) {
  go_a <- go_parts(params)
  go_b <- go_parts(params)
  s1 <- acc_parts("S1", params)
  is_par <- acc_parts("IS", params)
  s2 <- acc_parts("S2", params)

  if (component == "go") {
    if (identical(response, "A")) {
      return(acc_pdf(rt, go_a) * acc_survival(rt, go_b))
    }
    if (identical(response, "B")) {
      return(acc_pdf(rt, go_b) * acc_survival(rt, go_a))
    }
    stop("go component only supports observed A/B in this manual helper", call. = FALSE)
  }

  target <- if (identical(response, "A")) go_a else go_b
  competitor <- if (identical(response, "A")) go_b else go_a

  h_t <- ignore_pdf(rt, is_par, s2)
  H_t <- ignore_cdf(rt, is_par, s2)
  base <- acc_pdf(rt, target) * acc_survival(rt, competitor) *
    (acc_survival(rt, s1) + acc_cdf(rt, s1) * H_t)

  if (!is.finite(rt) || rt <= 0.0) {
    return(base)
  }

  inner <- integrate_vec(
    function(t) {
      acc_pdf(t, target) *
        acc_survival(t, competitor) *
        acc_cdf(t, s1)
    },
    0.0,
    rt
  )

  base + h_t * inner
}

manual_beest_like_visible_probability <- function(target, competitor, s1,
                                                  is_par, s2) {
  early <- integrate_vec(
    function(t) {
      acc_pdf(t, target) *
        acc_survival(t, competitor) *
        acc_survival(t, s1)
    },
    0.0,
    Inf
  )
  ignore_mass <- integrate_vec(function(t) ignore_pdf(t, is_par, s2), 0.0, Inf)
  late_target <- integrate_vec(
    function(t) {
      acc_pdf(t, target) *
        acc_survival(t, competitor) *
        acc_cdf(t, s1)
    },
    0.0,
    Inf
  )
  early + ignore_mass * late_target
}

manual_beest_like_mass <- function(component, params) {
  go_a <- go_parts(params)
  go_b <- go_parts(params)
  s1 <- acc_parts("S1", params)
  is_par <- acc_parts("IS", params)
  s2 <- acc_parts("S2", params)

  if (component == "go") {
    p_a <- integrate_vec(
      function(t) acc_pdf(t, go_a) * acc_survival(t, go_b),
      0.0,
      Inf
    )
    p_b <- integrate_vec(
      function(t) acc_pdf(t, go_b) * acc_survival(t, go_a),
      0.0,
      Inf
    )
    return(max(0.0, 1.0 - p_a - p_b))
  }

  p_a <- manual_beest_like_visible_probability(go_a, go_b, s1, is_par, s2)
  p_b <- manual_beest_like_visible_probability(go_b, go_a, s1, is_par, s2)
  max(0.0, 1.0 - p_a - p_b)
}

engine_loglik <- function(model, params, data_df) {
  prepared <- prepare_data(model, data_df)
  ctx <- make_context(model)
  params_df <- build_param_matrix(model$model_spec, params, trial_df = prepared)
  as.numeric(log_likelihood(ctx, prepared, params_df, min_ll = -1e12))
}

compare_manual_beest_like <- function(model, params, data_df) {
  component <- as.character(data_df$component[[1]])
  response <- as.character(data_df$R[[1]])
  rt <- data_df$rt[[1]]

  manual_loglik <- if (is.na(data_df$R[[1]])) {
    log(manual_beest_like_mass(component, params))
  } else {
    log(manual_beest_like_density(component, response, rt, params))
  }

  data.frame(
    component = component,
    response = response,
    rt = rt,
    engine_loglik = engine_loglik(model, params, data_df),
    manual_loglik = manual_loglik,
    diff = engine_loglik(model, params, data_df) - manual_loglik,
    stringsAsFactors = FALSE
  )
}

model <- stop_selective_model()

cases <- list(
  data.frame(trial = 1L, component = "go", R = "A", rt = 0.42, stringsAsFactors = FALSE),
  data.frame(trial = 1L, component = "go", R = "B", rt = 0.47, stringsAsFactors = FALSE),
  data.frame(trial = 1L, component = "stop", R = "A", rt = 0.42, stringsAsFactors = FALSE),
  data.frame(trial = 1L, component = "stop", R = "B", rt = 0.47, stringsAsFactors = FALSE),
  data.frame(trial = 1L, component = "stop", R = NA_character_, rt = NA_real_, stringsAsFactors = FALSE)
)

comparison <- do.call(
  rbind,
  lapply(cases, function(df) compare_manual_beest_like(model, stop_selective_params, df))
)

print(comparison)
