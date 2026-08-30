.model_view <- function(model) {
  prep <- if (inherits(model, "model_structure")) {
    model$prep
  } else {
    prepare_model(.validate_race_spec_input(model, "model view"))
  }
  list(
    accumulators = prep$accumulators,
    pools = prep$pools,
    outcomes = prep$outcomes
  )
}
