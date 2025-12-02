## Running `dev/scripts/test.R` Without `devtools::load_all()`

Some environments block `processx` (and therefore `devtools::load_all()`), which breaks `dev/scripts/test.R`. To exercise the dev script anyway, install the package into a temporary library and run the script from there:

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
   source('dev/scripts/test.R')
   RS
   ```

This exercises the native probability, simulation, and profiling flow without relying on `devtools::load_all()`. Clean up the temporary library afterwards if needed:

```bash
rm -rf /tmp/accumlib
```

## Profiling Native Code (Examples 1–3)

Use `dev/scripts/profile_examples.R` (invoked by `dev/scripts/profile_native_mac.sh`) to drive a native-only workload for the first three `new_api_examples`. This script loops over examples 1, 2, and 3, builds their generator structures, and runs the native log-likelihood so that external profilers (e.g., `perf`, Instruments) can attribute time inside `src/AccumulatR.so`.

1. Install the package into the temp library if you have not already (same as above).
2. Run the profiling workload:
   ```bash
   cd /Users/nstevenson/Documents/2025/UberUbermodel_newercopy
   R --vanilla --quiet -f dev/scripts/profile_examples.R
   ```
   Set `UUBER_PROFILE_TRIALS` to change the per-example trial count (default `2000`).
3. Wrap the command with your preferred profiler to capture C++ stacks. For example, on Linux:
   ```bash
   perf record -g -- R --vanilla --quiet -f dev/scripts/profile_examples.R
   perf report
   ```
   On macOS you can use Instruments’ “Time Profiler” template and point it at `R --vanilla --quiet -f dev/scripts/profile_examples.R`.
