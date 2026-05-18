test_separate_all_parameters <- function(spec) {
  lookup <- AccumulatR:::.default_parameter_lookup(spec)
  groups <- unique(unname(lookup))
  groups <- setdiff(groups[nzchar(groups)], AccumulatR:::.mixture_weight_parameter_names(spec))
  if (length(groups) == 0L) {
    return(spec)
  }
  set_parameters(spec, separate = stats::setNames(rep(list(TRUE), length(groups)), groups))
}
