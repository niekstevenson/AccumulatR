# Source from the repository root. No package is installed or modified.
# Each named case contains: spec (finalized model), data (prepared observations),
# parameters (list of parameter matrices), expand (original-trial row indices).
# A case can instead contain interleave: a list of such cases to alternate.
# baseline is a git revision or a directory containing its src/ snapshot.
# Both snapshots must use the current native context interface.
benchmark_cumulative <- function(cases, baseline, repetitions = 5, candidate = '.') {
  build <- tempfile('cumulative-')
  dir.create(build)
  dir.create(file.path(build,'baseline'))
  if (dir.exists(baseline)) {
    file.copy(file.path(baseline,'src'),file.path(build,'baseline'),recursive=TRUE)
  } else {
    archive <- file.path(build,'baseline.tar')
    system2('git',c('archive','-o',shQuote(archive),baseline,'src'))
    untar(archive,exdir=file.path(build,'baseline'))
  }
  dir.create(file.path(build,'current'))
  file.copy(file.path(candidate,'src'),file.path(build,'current'),recursive=TRUE)
  flags <- Sys.getenv('PKG_LIBS')
  on.exit(Sys.setenv(PKG_LIBS=flags))
  evaluators <- lapply(c('baseline','current'),function(version) {
    path <- file.path(build,version)
    Sys.setenv(PKG_LIBS='-framework Accelerate')
    math <- file.path(path,'src/lane_math.cpp')
    system2(file.path(R.home('bin'),'R'),c('CMD','SHLIB','--preclean','-o',
      shQuote(file.path(path,'lane_math.so')),shQuote(math)))
    Sys.setenv(PKG_LIBS=paste(shQuote(sub('cpp$','o',math)),'-framework Accelerate'))
    file.rename(file.path(path,'src/semantic_bridge.cpp'),file.path(path,'src/semantic_bridge.hpp'))
    file.copy('dev/scripts/bench_cumulative.cpp',path)
    env <- new.env()
    Rcpp::sourceCpp(file.path(path,'bench_cumulative.cpp'),env=env)
    env
  })
  names(evaluators) <- c('baseline','current')
  results <- lapply(cases,function(case) {
    values <- list()
    for (version in sample(names(evaluators))) {
      env <- evaluators[[version]]
      if (is.null(case$interleave)) {
        context <- env$benchmark_context(case$spec$prep)
        values[[version]] <- env$benchmark_particles(context,case$data,case$parameters,repetitions)
      } else {
        contexts <- lapply(case$interleave,function(x)env$benchmark_context(x$spec$prep))
        values[[version]] <- env$benchmark_interleaved(contexts,
          lapply(case$interleave,`[[`,'data'),lapply(case$interleave,`[[`,'parameters'),repetitions)
      }
    }
    difference <- if (is.null(case$interleave)) {
      colSums((values$current$loglik-values$baseline$loglik)[case$expand,,drop=FALSE])
    } else unlist(lapply(seq_along(case$interleave),function(i)
      colSums((values$current$loglik[[i]]-values$baseline$loglik[[i]])[
        case$interleave[[i]]$expand,,drop=FALSE])))
    list(seconds=vapply(values,function(x)median(x$seconds),numeric(1)),
         total_ll_difference=difference,
         results=values)
  })
  results
}
