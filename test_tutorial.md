## Running `dev/scripts/test.R` Without `devtools::load_all()`

Some environments block `processx` (and therefore `devtools::load_all()`), which breaks `dev/scripts/test.R`. To exercise the full test suite anyway, install the package into a temporary library and run the script from there:

1. **Install the package to a writable temp library**  
   ```bash
   cd /Users/nstevenson/Documents/2025/UberUbermodel_newercopy
   mkdir -p /tmp/accumlib
   R CMD INSTALL --library=/tmp/accumlib .
   ```

2. **Run the test script in a vanilla R session pointing at that library**  
   ```bash
   R --vanilla --quiet <<'RS'
   rm(list = ls())
   .libPaths(c('/tmp/accumlib', .libPaths()))
   library(AccumulatR)
   source('dev/examples/new_API.R')
   source('dev/scripts/compare_likelihoods.R')

   example_param_tables <- lapply(names(new_api_examples), function(id) {
     build_params_df(new_api_examples[[id]], new_api_example_params[[id]], n_trials = 1L)
   })
   names(example_param_tables) <- names(new_api_examples)
   example_ids <- names(new_api_examples)[c(1:11, 13:16, 18)]
   response_summary <- compare_response_suite(
     new_api_examples,
     example_param_tables,
     example_ids
   )
   print(response_summary[c('model', 'max_abs_diff')])
   names(response_summary$R) <- example_ids
   print(response_summary$R)

   param_example <- 'example_3_stop_na'
   param_results <- run_param_table_benchmark(
     model_spec = new_api_examples[[param_example]],
     core_params = new_api_example_params[[param_example]],
     n_trials = 1000L,
     seed = 2025,
     bench_reps = 3
   )
   cat('\\nLog-likelihood difference (native - R):', param_results$loglik$loglik_diff, '\\n')
   cat('Native batch used:', param_results$loglik$native_batch_used, '\\n')
   print(param_results$speed$timings)
   RS
   ```

This reproduces the regression table (`max_abs_diff = 0` for every example) and the benchmark/consistency checks, without relying on `devtools::load_all()`. Clean up the temporary library afterwards if needed:

```bash
rm -rf /tmp/accumlib
```
