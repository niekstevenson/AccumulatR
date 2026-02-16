suppressPackageStartupMessages({
  library(pkgload)
})

loaded <- FALSE
try({
  load_all('.', compile = FALSE, quiet = TRUE)
  loaded <- TRUE
}, silent = TRUE)
if (!loaded) {
  load_all('.', compile = TRUE, quiet = TRUE)
}

n_trials <- as.integer(Sys.getenv('N_TRIALS', '3000'))
if (!is.finite(n_trials) || n_trials < 100L) n_trials <- 3000L
seed <- as.integer(Sys.getenv('SEED', '123456'))

model_spec <- race_spec() |>
  add_accumulator('A', 'lognormal') |>
  add_accumulator('B', 'lognormal') |>
  add_accumulator('C', 'lognormal', onset = after('B')) |>
  add_outcome('A', 'A') |>
  add_outcome('C', 'C')

structure <- finalize_model(model_spec)

true_params <- c(
  A.meanlog = log(0.28),
  A.sdlog = 0.14,
  B.meanlog = log(0.10),
  B.sdlog = 0.10,
  C.meanlog = log(0.15),
  C.sdlog = 0.10
)

set.seed(seed)
params_df <- build_param_matrix(model_spec, true_params, n_trials = n_trials)
sim <- simulate(structure, params_df)

obs <- data.frame(
  trial = sim$trial,
  R = as.character(sim$R),
  rt = sim$rt,
  stringsAsFactors = FALSE
)

ctx <- build_likelihood_context(structure, obs)

ll_at <- function(theta) {
  pm <- build_param_matrix(
    model_spec,
    theta,
    n_trials = max(obs$trial),
    layout = ctx$param_layout
  )
  as.numeric(log_likelihood(ctx, pm))
}

ll_true <- ll_at(true_params)

# 1D fixed-others profiles (all but one parameter fixed at true)
profile_grid <- list(
  A.meanlog = seq(log(0.18), log(0.40), length.out = 41L),
  A.sdlog   = seq(0.06, 0.26, length.out = 41L),
  B.meanlog = seq(log(0.05), log(0.20), length.out = 51L),
  B.sdlog   = seq(0.04, 0.24, length.out = 51L),
  C.meanlog = seq(log(0.08), log(0.30), length.out = 51L),
  C.sdlog   = seq(0.04, 0.24, length.out = 51L)
)

profile_rows <- list()
summary_rows <- list()

cat('Computing 1D profiles...\n')
for (nm in names(profile_grid)) {
  g <- profile_grid[[nm]]
  vals <- numeric(length(g))
  for (i in seq_along(g)) {
    th <- true_params
    th[[nm]] <- g[[i]]
    vals[[i]] <- ll_at(th)
  }
  i_max <- which.max(vals)
  mle_1d <- g[[i_max]]
  max_ll <- vals[[i_max]]
  summary_rows[[nm]] <- data.frame(
    param = nm,
    true_value = true_params[[nm]],
    mle_1d = mle_1d,
    delta = mle_1d - true_params[[nm]],
    ll_true = ll_true,
    ll_max = max_ll,
    ll_drop_at_true = max_ll - ll_true,
    stringsAsFactors = FALSE
  )
  profile_rows[[nm]] <- data.frame(
    param = nm,
    value = g,
    ll = vals,
    delta_ll = vals - max_ll,
    stringsAsFactors = FALSE
  )
}

profile_1d <- do.call(rbind, profile_rows)
summary_1d <- do.call(rbind, summary_rows)
rownames(summary_1d) <- NULL

# 2D grid for B.meanlog/C.meanlog to diagnose identifiability ridge
cat('Computing 2D profile for B.meanlog x C.meanlog...\n')
b_grid <- seq(log(0.05), log(0.20), length.out = 31L)
c_grid <- seq(log(0.08), log(0.30), length.out = 31L)
ll_bc <- matrix(NA_real_, nrow = length(b_grid), ncol = length(c_grid))
for (i in seq_along(b_grid)) {
  for (j in seq_along(c_grid)) {
    th <- true_params
    th[['B.meanlog']] <- b_grid[[i]]
    th[['C.meanlog']] <- c_grid[[j]]
    ll_bc[i, j] <- ll_at(th)
  }
}
idx <- which(ll_bc == max(ll_bc), arr.ind = TRUE)[1, , drop = TRUE]
bc_best <- c(B.meanlog = b_grid[[idx[1]]], C.meanlog = c_grid[[idx[2]]])

# 2D grid for B.sdlog/C.sdlog
cat('Computing 2D profile for B.sdlog x C.sdlog...\n')
bs_grid <- seq(0.04, 0.24, length.out = 31L)
cs_grid <- seq(0.04, 0.24, length.out = 31L)
ll_bcs <- matrix(NA_real_, nrow = length(bs_grid), ncol = length(cs_grid))
for (i in seq_along(bs_grid)) {
  for (j in seq_along(cs_grid)) {
    th <- true_params
    th[['B.sdlog']] <- bs_grid[[i]]
    th[['C.sdlog']] <- cs_grid[[j]]
    ll_bcs[i, j] <- ll_at(th)
  }
}
idx_s <- which(ll_bcs == max(ll_bcs), arr.ind = TRUE)[1, , drop = TRUE]
bcs_best <- c(B.sdlog = bs_grid[[idx_s[1]]], C.sdlog = cs_grid[[idx_s[2]]])

# quick unconstrained fit from vignette-style start to compare with profile picture
cat('Running optim() fit from vignette-like start...\n')
neg_loglik <- function(theta) {
  th <- theta
  th[c('A.sdlog', 'B.sdlog', 'C.sdlog')] <- exp(th[c('A.sdlog', 'B.sdlog', 'C.sdlog')])
  -ll_at(th)
}
start <- c(
  A.meanlog = log(0.22),
  A.sdlog = log(0.10),
  B.meanlog = log(0.28),
  B.sdlog = log(0.10),
  C.meanlog = log(0.28),
  C.sdlog = log(0.10)
)
fit_far <- optim(start, neg_loglik, method = 'Nelder-Mead', control = list(maxit = 800))
fit_far_par <- fit_far$par
fit_far_par[c('A.sdlog', 'B.sdlog', 'C.sdlog')] <- exp(fit_far_par[c('A.sdlog', 'B.sdlog', 'C.sdlog')])

cat('Running optim() fit from near-true start...\n')
start_near <- c(
  A.meanlog = log(0.27),
  A.sdlog = log(0.13),
  B.meanlog = log(0.11),
  B.sdlog = log(0.11),
  C.meanlog = log(0.14),
  C.sdlog = log(0.11)
)
fit_near <- optim(start_near, neg_loglik, method = 'Nelder-Mead', control = list(maxit = 1200))
fit_near_par <- fit_near$par
fit_near_par[c('A.sdlog', 'B.sdlog', 'C.sdlog')] <- exp(fit_near_par[c('A.sdlog', 'B.sdlog', 'C.sdlog')])

out_dir <- file.path('dev', 'scripts', 'scratch_outputs')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(profile_1d, file.path(out_dir, 'profile_unobserved_chain_1d.csv'), row.names = FALSE)
write.csv(summary_1d, file.path(out_dir, 'profile_unobserved_chain_1d_summary.csv'), row.names = FALSE)
write.csv(
  expand.grid(B.meanlog = b_grid, C.meanlog = c_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
    transform(ll = as.vector(ll_bc), delta_ll = as.vector(ll_bc - max(ll_bc))),
  file.path(out_dir, 'profile_unobserved_chain_2d_Bmeanlog_Cmeanlog.csv'),
  row.names = FALSE
)
write.csv(
  expand.grid(B.sdlog = bs_grid, C.sdlog = cs_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
    transform(ll = as.vector(ll_bcs), delta_ll = as.vector(ll_bcs - max(ll_bcs))),
  file.path(out_dir, 'profile_unobserved_chain_2d_Bsdlog_Csdlog.csv'),
  row.names = FALSE
)
write.csv(
  data.frame(
    metric = c(
      'n_trials', 'll_true', 'll_2d_BC_meanlog_max',
      'bc_best_B.meanlog', 'bc_best_C.meanlog',
      'll_2d_BC_sdlog_max', 'bcs_best_B.sdlog', 'bcs_best_C.sdlog',
      'optim_far_negloglik',
      'optim_far_A.meanlog', 'optim_far_A.sdlog',
      'optim_far_B.meanlog', 'optim_far_B.sdlog',
      'optim_far_C.meanlog', 'optim_far_C.sdlog',
      'optim_near_negloglik',
      'optim_near_A.meanlog', 'optim_near_A.sdlog',
      'optim_near_B.meanlog', 'optim_near_B.sdlog',
      'optim_near_C.meanlog', 'optim_near_C.sdlog'
    ),
    value = c(
      n_trials, ll_true, max(ll_bc),
      bc_best[['B.meanlog']], bc_best[['C.meanlog']],
      max(ll_bcs), bcs_best[['B.sdlog']], bcs_best[['C.sdlog']],
      fit_far$value,
      fit_far_par[['A.meanlog']], fit_far_par[['A.sdlog']],
      fit_far_par[['B.meanlog']], fit_far_par[['B.sdlog']],
      fit_far_par[['C.meanlog']], fit_far_par[['C.sdlog']],
      fit_near$value,
      fit_near_par[['A.meanlog']], fit_near_par[['A.sdlog']],
      fit_near_par[['B.meanlog']], fit_near_par[['B.sdlog']],
      fit_near_par[['C.meanlog']], fit_near_par[['C.sdlog']]
    )
  ),
  file.path(out_dir, 'profile_unobserved_chain_metrics.csv'),
  row.names = FALSE
)

cat('\n=== PROFILE SUMMARY ===\n')
print(summary_1d)
cat('\n2D best (B.meanlog, C.meanlog):\n')
print(bc_best)
cat('\n2D best (B.sdlog, C.sdlog):\n')
print(bcs_best)
cat('\noptim recovered params (far start):\n')
print(fit_far_par)
cat('\noptim recovered params (near start):\n')
print(fit_near_par)
cat('\noutputs written to: ', out_dir, '\n', sep = '')
