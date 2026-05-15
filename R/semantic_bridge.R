.make_likelihood_context_prep <- function(prep, diagnostics = FALSE) {
  .Call(
    `_AccumulatR_semantic_make_likelihood_context_prep_cpp`,
    prep,
    isTRUE(diagnostics)
  )
}

.complexity_metrics_context <- function(context) {
  .Call(`_AccumulatR_semantic_complexity_metrics_context_cpp`, context)
}

.loglik_total_context <- function(context,
                                  params,
                                  data,
                                  ok = NULL,
                                  trial_weights = NULL,
                                  min_ll = log(1e-10)) {
  .Call(
    `_AccumulatR_semantic_loglik_total_context_cpp`,
    context,
    params,
    data,
    ok,
    trial_weights,
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
