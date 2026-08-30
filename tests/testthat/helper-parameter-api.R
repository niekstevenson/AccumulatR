test_separate_all_parameters <- function(spec) {
  groups <- par_names(spec)
  set_parameters(spec, separate = stats::setNames(rep(list(TRUE), length(groups)), groups))
}
