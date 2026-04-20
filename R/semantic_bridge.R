.compile_semantic_prep <- function(prep) {
  semantic_compile_prep_cpp(prep)
}

.validate_semantic_prep <- function(prep) {
  semantic_validate_prep_cpp(prep)
}

.project_semantic_prep <- function(prep) {
  semantic_project_prep_cpp(prep)
}

.lower_direct_prep <- function(prep) {
  semantic_lower_direct_prep_cpp(prep)
}

.lower_exact_prep <- function(prep) {
  semantic_lower_exact_prep_cpp(prep)
}

.direct_loglik_prep <- function(prep,
                                params,
                                data,
                                min_ll = log(1e-10)) {
  semantic_direct_loglik_prep_cpp(
    prep,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE),
    as.numeric(min_ll)
  )
}

.exact_loglik_prep <- function(prep,
                               params,
                               data,
                               min_ll = log(1e-10)) {
  semantic_exact_loglik_prep_cpp(
    prep,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE),
    as.numeric(min_ll)
  )
}

.direct_prob_prep <- function(prep,
                              params,
                              data) {
  semantic_direct_prob_prep_cpp(
    prep,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE)
  )
}

.exact_prob_prep <- function(prep,
                             params,
                             data) {
  semantic_exact_prob_prep_cpp(
    prep,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE)
  )
}

.observed_loglik_prep <- function(prep,
                                  params,
                                  data,
                                  min_ll = log(1e-10)) {
  semantic_observed_loglik_prep_cpp(
    prep,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE),
    as.numeric(min_ll)
  )
}
