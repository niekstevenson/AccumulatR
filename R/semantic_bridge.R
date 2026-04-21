.make_likelihood_context_prep <- function(prep) {
  .Call(`_AccumulatR_semantic_make_likelihood_context_prep_cpp`, prep)
}

.prepare_trial_layout <- function(data) {
  .Call(
    `_AccumulatR_semantic_prepare_trial_layout_cpp`,
    as.data.frame(data, stringsAsFactors = FALSE)
  )
}

.loglik_context <- function(context,
                            layout,
                            params,
                            data,
                            min_ll = log(1e-10)) {
  .Call(
    `_AccumulatR_semantic_loglik_context_cpp`,
    context,
    layout,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE),
    as.numeric(min_ll)
  )
}
