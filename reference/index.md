# Package index

## Package Overview

Start here for the main package concepts and workflow.

- [`AccumulatR`](AccumulatR-package.md)
  [`AccumulatR-package`](AccumulatR-package.md) : AccumulatR: Simulate
  and fit evidence-accumulation models

## Model Specification

Define accumulators, outcomes, components, triggers, and finalize a race
model.

- [`race_spec()`](race_spec.md) : Start a race-model specification
- [`add_accumulator()`](add_accumulator.md) : Add an accumulator to a
  model
- [`add_pool()`](add_pool.md) : Pool several accumulators under a shared
  label
- [`add_outcome()`](add_outcome.md) : Define an observed response
- [`add_component()`](add_component.md) : Define a mixture component
- [`add_shared_params()`](add_shared_params.md) : Tie parameters across
  accumulators
- [`add_trigger()`](add_trigger.md) : Add a shared trigger or gate
- [`set_metadata()`](set_metadata.md) : Store model-level settings
- [`set_mixture_options()`](set_mixture_options.md) : Control how
  mixture components are combined
- [`nest_accumulators()`](nest_accumulators.md) : Expand behavioral data
  to one row per accumulator
- [`finalize_model()`](finalize_model.md) : Compile a model for
  simulation and fitting

## Outcome Expressions and Timing

Compose outcome rules and onset dependencies.

- [`build_outcome_expr()`](build_outcome_expr.md) : Turn a response rule
  into an internal expression
- [`all_of()`](all_of.md) : Define a response that requires several
  processes to finish
- [`first_of()`](first_of.md) : Define a response that occurs when the
  first listed process finishes
- [`none_of()`](none_of.md) [`exclude()`](none_of.md) : Define the
  absence of an event
- [`inhibit()`](inhibit.md) : Define a response that is blocked by
  another process
- [`expr_guard()`](expr_guard.md) : Build a blocking rule explicitly
- [`after()`](after.md) : Start one accumulator after another process
  finishes

## Data Preparation

Inspect model parameters and prepare trial-level parameter matrices.

- [`sampled_pars()`](sampled_pars.md) : List the free parameters implied
  by a model
- [`build_param_matrix()`](build_param_matrix.md) : Create trial-level
  parameter values

## Likelihood and Simulation

Build native likelihood contexts, simulate observations, and evaluate
probabilities.

- [`build_likelihood_context()`](build_likelihood_context.md) : Prepare
  behavioral data for repeated likelihood evaluation
- [`ensure_native_ctx()`](ensure_native_ctx.md) : Refresh the compiled
  backend inside a context
- [`log_likelihood()`](log_likelihood.md) : Evaluate the log-likelihood
  of behavioral data
- [`response_probabilities()`](response_probabilities.md) : Compute
  predicted response probabilities
- [`simulate()`](simulate.md) : Simulate behavioral data from a model

## Visualization

Plot and inspect model structure.

- [`processing_tree()`](processing_tree.md) : Draw a processing tree for
  the responses in a model
- [`plot_accumulators()`](plot_accumulators.md) : Plot accumulators and
  their timing relations
