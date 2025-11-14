if (!exists(".precompile_likelihood_expressions", mode = "function")) {
  pkg_root <- testthat::test_path("..", "..")
  r_sources <- file.path(
    pkg_root,
    c(
      "R/dist.R",
      "R/likelihood_common.R",
      "R/likelihood_cache.R",
      "R/pool_math.R",
      "R/likelihood_primitives.R",
      "R/likelihood_prep.R",
      "R/likelihood_integrate.R",
      "R/likelihood_kernels.R"
    )
  )
  invisible(lapply(r_sources, base::source))
}

has_build_tools <- TRUE
if (requireNamespace("pkgbuild", quietly = TRUE)) {
  has_build_tools <- isTRUE(pkgbuild::has_build_tools())
}

build_test_prep <- function() {
  old_dir <- getwd()
  pkg_root <- testthat::test_path("..", "..")
  setwd(pkg_root)
  on.exit(setwd(old_dir), add = TRUE)

  acc_defs <- list(
    A = list(
      dist = "lognormal",
      params = list(meanlog = -0.2, sdlog = 0.45),
      onset = 0.1,
      q = 0.0
    ),
    B = list(
      dist = "gamma",
      params = list(shape = 2.0, rate = 3.0),
      onset = 0.0,
      q = 0.0
    )
  )

  pool_defs <- list(
    P = list(
      members = c("A", "B"),
      k = 1L
    )
  )

  outcomes <- list(
    accA = list(expr = list(kind = "event", source = "A")),
    poolP = list(expr = list(kind = "event", source = "P")),
    andAB = list(expr = list(
      kind = "and",
      args = list(
        list(kind = "event", source = "A"),
        list(kind = "event", source = "B")
      )
    )),
    orAB = list(expr = list(
      kind = "or",
      args = list(
        list(kind = "event", source = "A"),
        list(kind = "event", source = "B")
      )
    )),
    notB = list(expr = list(
      kind = "not",
      arg = list(kind = "event", source = "B")
    )),
    deadlinePseudo = list(expr = list(kind = "event", source = "__DEADLINE__"))
  )

  prep <- list(
    accumulators = acc_defs,
    pools = pool_defs,
    outcomes = outcomes,
    components = list(
      ids = "__default__",
      weights = 1,
      attrs = list(`__default__` = list(deadline = Inf, guess = NULL)),
      has_weight_param = FALSE
    ),
    default_deadline = Inf,
    special_outcomes = list()
  )

  prep <- .precompile_likelihood_expressions(prep)
  prep <- .refresh_compiled_prep_refs(prep)

  id_index <- setNames(
    seq_along(c(names(acc_defs), names(pool_defs))),
    c(names(acc_defs), names(pool_defs))
  )

  prep$.runtime <- list(
    expr_compiled = prep[[".expr_compiled"]],
    label_cache = new.env(parent = emptyenv(), hash = TRUE),
    competitor_map = list(),
    id_index = id_index,
    pool_members_cache = new.env(parent = emptyenv(), hash = TRUE)
  )
  prep[[".id_index"]] <- id_index
  prep[[".label_cache"]] <- prep$.runtime$label_cache
  prep <- .prep_set_cache_bundle(prep, .build_likelihood_cache_bundle(prep))
  prep <- .refresh_compiled_prep_refs(prep)
  prep
}

with_node_option <- function(flag, expr) {
  old <- getOption("uuber.use.native.node.eval")
  on.exit(options(uuber.use.native.node.eval = old), add = TRUE)
  options(uuber.use.native.node.eval = flag)
  force(expr)
}

test_that("native node evaluator matches R kernels when enabled", {
  skip_if_not(has_build_tools, "Native build tools are unavailable")
  prep <- build_test_prep()
  compiled_nodes <- prep[[".expr_compiled"]]$nodes
  node_ids <- vapply(
    prep$outcomes,
    function(outcome) attr(outcome$expr, ".lik_id", exact = TRUE),
    integer(1)
  )
  expect_false(any(is.na(node_ids)))

  test_nodes <- c(
    accA = node_ids[["accA"]],
    poolP = node_ids[["poolP"]],
    andAB = node_ids[["andAB"]],
    orAB = node_ids[["orAB"]],
    notB = node_ids[["notB"]]
  )

  forced_survive_vec <- as.integer(prep[[".id_index"]][["A"]])
  forced_survive_vec <- forced_survive_vec[!is.na(forced_survive_vec)]
  expect_true(length(forced_survive_vec) > 0L)

  eval_entry <- function(node, tt,
                         forced_complete = integer(0),
                         forced_survive = integer(0)) {
    state <- .eval_state_create()
    dens <- .node_density(
      node, tt, prep, NULL,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
    surv <- .node_survival_cond(
      node, tt, prep, NULL,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
    cdf <- .node_cdf_cond(
      node, tt, prep, NULL,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
    c(density = dens, survival = surv, cdf = cdf)
  }

  collect_results <- function(flag) {
    with_node_option(flag, {
      res <- list()
      times <- c(0.2, 0.5, 1.0)
      for (label in names(test_nodes)) {
        node <- compiled_nodes[[test_nodes[[label]]]]
        for (tt in times) {
          key <- sprintf("%s|%0.2f", label, tt)
          res[[key]] <- eval_entry(node, tt)
        }
      }
      forced_node <- compiled_nodes[[test_nodes[["orAB"]]]]
      res[["orAB|forced"]] <- eval_entry(
        forced_node,
        0.65,
        forced_survive = forced_survive_vec
      )
      res
    })
  }

  baseline <- collect_results(FALSE)
  native <- collect_results(TRUE)

  expect_equal(
    unlist(native),
    unlist(baseline),
    tolerance = 1e-10
  )
})

test_that("native node evaluator handles pseudo labels", {
  skip_if_not(has_build_tools, "Native build tools are unavailable")
  prep <- build_test_prep()
  compiled_nodes <- prep[[".expr_compiled"]]$nodes
  pseudo_id <- attr(prep$outcomes$deadlinePseudo$expr, ".lik_id", exact = TRUE)
  pseudo_node <- compiled_nodes[[pseudo_id]]
  expect_false(is.null(pseudo_node))

  times <- c(-0.5, 0, 0.5)
  expected_cdf <- vapply(times, function(tt) .node_cdf_cond(pseudo_node, tt, prep, NULL), numeric(1))
  expected_surv <- vapply(times, function(tt) .node_survival_cond(pseudo_node, tt, prep, NULL), numeric(1))

  native_cdf <- with_node_option(TRUE, {
    vapply(times, function(tt) .node_cdf_cond(pseudo_node, tt, prep, NULL), numeric(1))
  })
  native_surv <- with_node_option(TRUE, {
    vapply(times, function(tt) .node_survival_cond(pseudo_node, tt, prep, NULL), numeric(1))
  })

  expect_equal(native_cdf, expected_cdf)
  expect_equal(native_surv, expected_surv)
})

test_that("native outcome probability matches R fallback for infinite upper", {
  prep <- build_test_prep()
  expr <- prep$outcomes$andAB$expr
  state <- .eval_state_create()
  integrand <- function(t) {
    .integrand_outcome_density(
      t,
      expr = expr,
      prep = prep,
      component = NULL,
      competitor_exprs = list(),
      state = state
    )
  }
  ref <- .native_integrate(
    integrand,
    lower = 0,
    upper = Inf,
    rel.tol = .integrate_rel_tol(),
    abs.tol = .integrate_abs_tol(),
    max.depth = 12L
  )
  native_val <- .integrate_outcome_probability(
    expr,
    prep,
    component = NULL,
    upper_limit = Inf,
    competitor_exprs = list(),
    state = .eval_state_create()
  )
  expect_equal(native_val, ref, tolerance = 1e-8)
})
