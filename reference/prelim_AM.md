# Preliminary Assessments in MSE

Evaluates the likely performance of Assessment models in the operating
model. This function will apply the assessment model for Data generated
during the historical period of the MSE, and report the convergence rate
for the model and total time elapsed in running the assessments.

## Usage

``` r
prelim_AM(x, Assess, ncpus = NULL, ...)
```

## Arguments

- x:

  Either a `Hist`, `Data` or `OM` object.

- Assess:

  An Assess function of class `Assess`.

- ncpus:

  Numeric, the number of CPUs to run the Assessment model (will run in
  parallel if greater than 1).

- ...:

  Arguments to be passed to `Assess`, e.g., model configurations.

## Value

Returns invisibly a list of
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
objects of length `OM@nsim`. Messages via console.

## Author

Q. Huynh

## Examples

``` r
# \donttest{
prelim_AM(MSEtool::SimulatedData, SP)
#> ✔ Running SP with 3 simulations for MSEtool::SimulatedData.
#> ✔ Assessments complete.
#> ✔ Total time to run 3 assessments: 0.1 seconds
#> ✔ 0 of 3 simulations (0%) failed to converge.
# }
```
