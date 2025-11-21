test_that(".prep_native_payload strips R closures", {
  skip_if_not_installed("AccumulatR")

  source("dev/examples/new_API.R")
  spec <- new_api_examples[[1]]
  params <- new_api_example_params[[1]]
  structure <- build_generator_structure(apply_core_params_to_spec(spec, params))
  prep <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  prep <- AccumulatR:::.likelihood_seed_prep_from_params(prep, build_params_df(spec, params))
  prep <- AccumulatR:::.prep_set_cache_bundle(prep, AccumulatR:::.build_likelihood_cache_bundle(prep))

  trimmed <- AccumulatR:::.prep_native_payload(prep)

  expect_type(trimmed, "list")
  expect_true(is.null(trimmed$.expr_compiled$nodes[[1]]$density_fn))
  expect_true(is.null(trimmed$.expr_compiled$nodes[[1]]$ops))
})

