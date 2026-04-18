.semantic_bridge_env <- new.env(parent = baseenv())

.load_semantic_bridge <- function(rebuild = FALSE, root = ".") {
  if (!rebuild && isTRUE(.semantic_bridge_env$loaded)) {
    return(invisible(NULL))
  }
  root <- normalizePath(root, mustWork = TRUE)
  Rcpp::sourceCpp(
    file.path(root, "src", "bridge", "semantic_bridge.cpp"),
    env = .semantic_bridge_env,
    rebuild = rebuild,
    showOutput = FALSE,
    verbose = FALSE
  )
  .semantic_bridge_env$loaded <- TRUE
  invisible(NULL)
}

.compile_semantic_prep <- function(prep, rebuild = FALSE, root = ".") {
  .load_semantic_bridge(rebuild = rebuild, root = root)
  .semantic_bridge_env$semantic_compile_prep_cpp(prep)
}

.validate_semantic_prep <- function(prep, rebuild = FALSE, root = ".") {
  .load_semantic_bridge(rebuild = rebuild, root = root)
  .semantic_bridge_env$semantic_validate_prep_cpp(prep)
}

.project_semantic_prep <- function(prep, rebuild = FALSE, root = ".") {
  .load_semantic_bridge(rebuild = rebuild, root = root)
  .semantic_bridge_env$semantic_project_prep_cpp(prep)
}

.lower_direct_prep <- function(prep, rebuild = FALSE, root = ".") {
  .load_semantic_bridge(rebuild = rebuild, root = root)
  .semantic_bridge_env$semantic_lower_direct_prep_cpp(prep)
}

.lower_exact_prep <- function(prep, rebuild = FALSE, root = ".") {
  .load_semantic_bridge(rebuild = rebuild, root = root)
  .semantic_bridge_env$semantic_lower_exact_prep_cpp(prep)
}

.direct_loglik_prep <- function(prep,
                                params,
                                data,
                                min_ll = log(1e-10),
                                rebuild = FALSE,
                                root = ".") {
  .load_semantic_bridge(rebuild = rebuild, root = root)
  .semantic_bridge_env$semantic_direct_loglik_prep_cpp(
    prep,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE),
    as.numeric(min_ll)
  )
}

.exact_loglik_prep <- function(prep,
                               params,
                               data,
                               min_ll = log(1e-10),
                               rebuild = FALSE,
                               root = ".") {
  .load_semantic_bridge(rebuild = rebuild, root = root)
  .semantic_bridge_env$semantic_exact_loglik_prep_cpp(
    prep,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE),
    as.numeric(min_ll)
  )
}
