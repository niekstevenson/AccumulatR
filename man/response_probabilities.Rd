% Generated manually to document the user-facing S3 generic.
\name{response_probabilities}
\alias{response_probabilities}
\alias{response_probabilities.model_structure}
\alias{response_probabilities.default}
\title{Evaluate marginal response probabilities}
\usage{
response_probabilities(structure, params_df, include_na = TRUE, ...)

\method{response_probabilities}{model_structure}(
  structure,
  params_df,
  include_na = TRUE,
  ...
)

\method{response_probabilities}{default}(
  structure,
  params_df,
  include_na = TRUE,
  ...
)
}
\arguments{
\item{structure}{Finalized model structure.}

\item{params_df}{A parameter data frame or rectangular parameter matrix.}

\item{include_na}{If \code{TRUE}, include residual mass as \code{"NA"}.}

\item{...}{Unused; for S3 compatibility.}
}
\value{
A named numeric vector of marginal response probabilities.
}
\description{
Build query trials for each observable response label, evaluate those
probabilities with the compiled likelihood kernels, and aggregate the results
across mixture components and parameter trials.
}
\examples{
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params_df <- build_param_matrix(
  spec,
  c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0),
  n_trials = 1
)
response_probabilities(structure, params_df)
}
