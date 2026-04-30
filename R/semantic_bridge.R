.make_likelihood_context_prep <- function(prep) {
  .Call(`_AccumulatR_semantic_make_likelihood_context_prep_cpp`, prep)
}

.loglik_context <- function(context,
                            params,
                            data,
                            ok = NULL,
                            expand = NULL,
                            min_ll = log(1e-10)) {
  .Call(
    `_AccumulatR_semantic_loglik_context_cpp`,
    context,
    params,
    data,
    ok,
    expand,
    min_ll
  )
}

.probability_context <- function(context,
                                 params,
                                 data) {
  .Call(
    `_AccumulatR_semantic_probability_context_cpp`,
    context,
    params,
    data
  )
}
