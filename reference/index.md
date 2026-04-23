# Package index

## Package Overview

Documentation Overview

- [`AccumulatR`](https://niekstevenson.github.io/AccumulatR/reference/AccumulatR-package.md)
  [`AccumulatR-package`](https://niekstevenson.github.io/AccumulatR/reference/AccumulatR-package.md)
  : AccumulatR: Simulate and fit evidence-accumulation models

## Model Specification

Define accumulators, outcomes, components, triggers, and finalize a race
model.

- [`race_spec()`](https://niekstevenson.github.io/AccumulatR/reference/race_spec.md)
  : Start a race-model specification
- [`add_accumulator()`](https://niekstevenson.github.io/AccumulatR/reference/add_accumulator.md)
  : Add an accumulator to a model
- [`add_pool()`](https://niekstevenson.github.io/AccumulatR/reference/add_pool.md)
  : Pool several accumulators under a shared label
- [`add_outcome()`](https://niekstevenson.github.io/AccumulatR/reference/add_outcome.md)
  : Define an observed response
- [`add_component()`](https://niekstevenson.github.io/AccumulatR/reference/add_component.md)
  : Define a mixture component
- [`set_parameters()`](https://niekstevenson.github.io/AccumulatR/reference/set_parameters.md)
  : Define the external parameter names for a model
- [`add_trigger()`](https://niekstevenson.github.io/AccumulatR/reference/add_trigger.md)
  : Add a shared trigger or gate
- [`set_metadata()`](https://niekstevenson.github.io/AccumulatR/reference/set_metadata.md)
  : Store model-level metadata
- [`set_mixture_options()`](https://niekstevenson.github.io/AccumulatR/reference/set_mixture_options.md)
  : Control how mixture components are combined
- [`finalize_model()`](https://niekstevenson.github.io/AccumulatR/reference/finalize_model.md)
  : Compile a model for simulation and fitting

## Outcome Expressions and Timing

Compose outcome rules and onset dependencies.

- [`build_outcome_expr()`](https://niekstevenson.github.io/AccumulatR/reference/build_outcome_expr.md)
  : Turn a response rule into an internal expression
- [`all_of()`](https://niekstevenson.github.io/AccumulatR/reference/all_of.md)
  : Define a response that requires several processes to finish
- [`first_of()`](https://niekstevenson.github.io/AccumulatR/reference/first_of.md)
  : Define a response that occurs when the first listed process finishes
- [`none_of()`](https://niekstevenson.github.io/AccumulatR/reference/none_of.md)
  : Define the absence of an event
- [`inhibit()`](https://niekstevenson.github.io/AccumulatR/reference/inhibit.md)
  : Define a response that is blocked by another process
- [`expr_guard()`](https://niekstevenson.github.io/AccumulatR/reference/expr_guard.md)
  : Build a blocking rule explicitly
- [`after()`](https://niekstevenson.github.io/AccumulatR/reference/after.md)
  : Start one accumulator after another process finishes

## Data Preparation

Inspect model parameters and prepare trial-level parameter matrices.

- [`prepare_data()`](https://niekstevenson.github.io/AccumulatR/reference/prepare_data.md)
  : Prepare behavioral data for likelihood evaluation
- [`sampled_pars()`](https://niekstevenson.github.io/AccumulatR/reference/sampled_pars.md)
  : List the free parameters implied by a model
- [`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md)
  : Create trial-level parameter values

## Likelihood and Simulation

Build native likelihood contexts, simulate observations, and evaluate
probabilities.

- [`make_context()`](https://niekstevenson.github.io/AccumulatR/reference/make_context.md)
  : Build a compiled likelihood context from a model
- [`log_likelihood()`](https://niekstevenson.github.io/AccumulatR/reference/log_likelihood.md)
  : Evaluate the log-likelihood of behavioral data
- [`simulate()`](https://niekstevenson.github.io/AccumulatR/reference/simulate.md)
  : Simulate behavioral data from a model
- [`response_probabilities()`](https://niekstevenson.github.io/AccumulatR/reference/response_probabilities.md)
  : Evaluate marginal response probabilities

## Visualization

Plot and inspect model structure.

- [`processing_tree()`](https://niekstevenson.github.io/AccumulatR/reference/processing_tree.md)
  : Draw a processing tree for the responses in a model
- [`plot_accumulators()`](https://niekstevenson.github.io/AccumulatR/reference/plot_accumulators.md)
  : Plot accumulators and their timing relations
