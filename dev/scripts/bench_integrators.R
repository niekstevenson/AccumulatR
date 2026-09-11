# Run from the AccumulatR repository. Nothing is installed or changed in src/.
library(AccumulatR)
library(EMC2)
library(Rcpp)
Sys.setenv(VECLIB_MAXIMUM_THREADS = '1', OMP_NUM_THREADS = '1')
fits <- '/Users/nstevenson/Documents/2026/StopData/3B_stopchange/WeberResults/results_RDEX'
subjects <- c(105, 115, 116, 118, 121, 126)
draws <- 100L

# Compile the real executor, replacing only its composite-integration call.
build <- tempfile('quadrature-')
dir.create(build)
file.copy('src', build, recursive = TRUE)
file.rename(file.path(build, 'src/semantic_bridge.cpp'), file.path(build, 'src/semantic_bridge.hpp'))
file.copy(paste0('dev/scripts/bench_integrators.', c('cpp', 'hpp')), build)
header <- file.path(build, 'src/eval/exact_adaptive.hpp')
code <- readLines(header)
code <- sub('inline constexpr double kAdaptive', 'inline double kAdaptive', code, fixed = TRUE)
writeLines(code, header)
header <- file.path(build, 'src/eval/exact_compiled_lane_eval.hpp')
code <- readLines(header)
code <- sub('#include "exact_adaptive.hpp"', '#include "bench_integrators.hpp"', code, fixed = TRUE)
code <- sub('  adaptive_integrate_lane_batch(', '  benchmark_integrate_lane_batch(', code, fixed = TRUE)
writeLines(code, header)
compiler_flags <- Sys.getenv(c('PKG_CPPFLAGS', 'PKG_LIBS'))
Sys.setenv(PKG_LIBS = '-framework Accelerate')
system2(file.path(R.home('bin'), 'R'), c('CMD', 'SHLIB', '-o',
  file.path(build, 'lane_math.so'), file.path(build, 'src/lane_math.cpp')))
Sys.setenv(PKG_CPPFLAGS = paste('-I', build),
           PKG_LIBS = paste(file.path(build, 'src/lane_math.o'), '-framework Accelerate'))
sourceCpp(file.path(build, 'bench_integrators.cpp'))
do.call(Sys.setenv, as.list(compiler_flags))

methods <- read.table(header = TRUE, text = '
name                method absolute relative
GL31                    31        0        0
GK15_1e-3_abs1e-12        0    1e-12     1e-3
GK15_1e-4_abs1e-12        0    1e-12     1e-4
GK15_previous            0     1e-8     1e-6
reference                0    1e-12    1e-10
')

run_subject <- function(subject, methods, draws = 100L, repetitions = 3L) {
  load(file.path(fits, paste0('subject_', subject, '.RData')))
  particles <- do.call(rbind, lapply(emc, function(ch) t(ch$samples$alpha[, 1, ])))
  index <- unique(round(seq(1, nrow(particles), length.out = draws)))
  particles <- particles[index, , drop = FALSE]
  data <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model
  context <- make_context(model()$spec)$cpp$native
  parameters <- lapply(seq_len(nrow(particles)), function(i) {
    pars <- EMC2:::get_pars_oo(particles[i, ], data, model)
    EMC2:::.accumulatr_runtime_parameters(pars, attr(data, 'AccumulatR_bridge')$bridge)
  })
  result <- list()
  for (i in sample(seq_len(nrow(methods)))) {
    m <- methods[i, ]
    cat(subject, m$name, '\n'); flush.console()
    result[[m$name]] <- benchmark_integrators(context, data, parameters,
      m$method, m$absolute, m$relative, repetitions)
    cat('  seconds:', result[[m$name]]$seconds, '\n'); flush.console()
  }
  list(index = index, expand = attr(data, 'expand'), methods = result)
}

results <- list()
dir.create('dev/scripts/scratch_outputs', showWarnings = FALSE)
set.seed(19)
for (subject in subjects) {
  results[[as.character(subject)]] <- run_subject(subject, methods, draws)
  saveRDS(results, 'dev/scripts/scratch_outputs/benchmark_integrators.rds')
}

rows <- list()
for (subject in names(results)) {
  result <- results[[subject]]
  reference <- result$methods$reference$loglik
  for (name in names(result$methods)) {
    m <- result$methods[[name]]
    error <- m$loglik - reference
    delta <- colSums(error[result$expand, , drop = FALSE])
    rows[[length(rows) + 1L]] <- data.frame(subject = subject, method = name,
      seconds_100 = 100 * median(m$seconds) / ncol(error),
      max_density_error_percent = 100 * max(abs(expm1(error))),
      max_abs_loglik_error = max(abs(delta)), sd_loglik_error = sd(delta))
  }
}
rows <- do.call(rbind, rows)
print(rows, digits = 4, row.names = FALSE)
write.csv(rows, 'dev/scripts/scratch_outputs/benchmark_integrators_summary.csv', row.names = FALSE)
