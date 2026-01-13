# Sample posterior of TMB models in SAMtool

A convenient wrapper function (`posterior`) to sample the posterior
using MCMC in rstan and returns a `stanfit` object for diagnostics. Use
`RCMstan` to update the RCM and the enclosed operating model with MCMC
samples..

## Usage

``` r
posterior(x, ...)

# S4 method for class 'RCModel'
posterior(
  x,
  priors_only = FALSE,
  laplace = FALSE,
  chains = 2,
  iter = 2000,
  warmup = floor(iter/2),
  thin = 5,
  seed = 34,
  init = "last.par.best",
  cores = chains,
  ...
)

# S4 method for class 'Assessment'
posterior(x, priors_only = FALSE, ...)

RCMstan(RCModel, stanfit, sim, cores = 1, silent = FALSE)
```

## Arguments

- x:

  An object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
  or [RCModel](https://samtool.openmse.com/reference/RCModel-class.md).

- ...:

  Additional arguments to pass to
  [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  via
  [`tmbstan::tmbstan`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html).

- priors_only:

  Logical, whether to set the likelihood to zero and sample the priors
  only.

- laplace:

  Logical, whether to do the Laplace approximation for random
  parameters.

- chains:

  The numer of MCMC chains.

- iter:

  The number of iterations for each chain, including warmup.

- warmup:

  The number of burnin iterations

- thin:

  The frequency at which iterations are kept (e.g., `5` saves every
  fifth iteration)

- seed:

  Seed for random number generator during the MCMC.

- init:

  The initial values of parameters for starting the MCMC chain. See
  [`tmbstan::tmbstan`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html).

- cores:

  The number of cores for running in parallel, e.g., one core per MCMC
  chain. Used in `RCMstan` for reconstructing the population.

- RCModel:

  An object of class `RCModel`

- stanfit:

  An object of class `stanfit` returned by `posterior`.

- sim:

  A matrix of `RCModel@OM@nsim` rows and 2 columns that specifies the
  samples used to update the operating model. The first column specifies
  the chain and the second columns specifies the MCMC iteration.

- silent:

  Logical to indicate if progress messages should be printed to console.

## Value

`posterior` returns an object of class `stanfit`. See `class?stanfit`.

`RCMstan` returns an updated `RCModel`.

## Online Documentation

A vignette on the steps to run the MCMC is available on the openMSE
[website](https://openmse.com/tutorial-rcm/4-case-study-mcmc/).

## Author

Q. Huynh
