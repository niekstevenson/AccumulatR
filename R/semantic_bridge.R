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

.profiler_available <- function() {
  .Call(`_AccumulatR_semantic_profiler_is_available_cpp`)
}

.profiler_start <- function(path, frequency = 500L) {
  .Call(
    `_AccumulatR_semantic_profiler_start_cpp`,
    as.character(path),
    as.integer(frequency)
  )
}

.profiler_flush <- function() {
  .Call(`_AccumulatR_semantic_profiler_flush_cpp`)
}

.profiler_stop <- function() {
  .Call(`_AccumulatR_semantic_profiler_stop_cpp`)
}
