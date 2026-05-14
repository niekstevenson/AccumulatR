.make_likelihood_context_prep <- function(prep) {
  .Call(`_AccumulatR_semantic_make_likelihood_context_prep_cpp`, prep)
}

.complexity_metrics_context <- function(context) {
  .Call(`_AccumulatR_semantic_complexity_metrics_context_cpp`, context)
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

.response_probabilities_context <- function(context,
                                            params,
                                            layout) {
  .Call(
    `_AccumulatR_semantic_response_probabilities_context_cpp`,
    context,
    params,
    layout
  )
}
