# Generate simulated data from TMB models in SAMtool

A convenient wrapper function (`simulate`) to simulate data (and process
error) from the likelihood function.

## Usage

``` r
simulate(object, ...)

# S4 method for class 'Assessment'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  process_error = FALSE,
  refit = FALSE,
  cores = 1,
  ...
)

# S4 method for class 'RCModel'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  process_error = FALSE,
  refit = FALSE,
  cores = 1,
  ...
)
```

## Arguments

- object:

  An object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
  or [RCModel](https://samtool.openmse.com/reference/RCModel-class.md)
  containing the fitted model.

- ...:

  Additional arguments

- nsim:

  Number of simulations

- seed:

  Used for the random number generator

- process_error:

  Logical, indicates if process error is re-sampled in the simulation.

- refit:

  Logical, whether to re-fit the model for each simulated dataset.

- cores:

  The number of CPUs for parallel processing for model re-fitting if
  `refit = TRUE`.

## Value

A [sim](https://samtool.openmse.com/reference/sim-class.md) object
returning the original data, simulated data, original parameters,
parameters estimated from simulated data, and process error used to
simulate data. then a nested list of model output (`opt`, `SD`, and
`report`).

## Details

Process error, e.g., recruitment deviations, will be re-sampled in the
simulation.

## Author

Q. Huynh
