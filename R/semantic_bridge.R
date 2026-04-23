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
                            ok = NULL,
                            expand = NULL,
                            min_ll = log(1e-10)) {
  .Call(
    `_AccumulatR_semantic_loglik_context_cpp`,
    context,
    layout,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE),
    if (is.null(ok)) NULL else as.logical(ok),
    if (is.null(expand)) NULL else as.integer(expand),
    as.numeric(min_ll)
  )
}

.probability_context <- function(context,
                                 layout,
                                 params,
                                 data) {
  .Call(
    `_AccumulatR_semantic_probability_context_cpp`,
    context,
    layout,
    as.matrix(params),
    as.data.frame(data, stringsAsFactors = FALSE)
  )
}
