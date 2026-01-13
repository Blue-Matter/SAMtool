# Simple Stock Synthesis

A simple age-structured model
([SCA_Pope](https://samtool.openmse.com/reference/SCA.md)) fitted to a
time series of catch going back to unfished conditions. Terminal
depletion (ratio of current total biomass to unfished biomass) is by
default fixed to 0.4. Selectivity is fixed to the maturity ogive,
although it can be overridden with the start argument. The sole
parameter estimated is R0 (unfished recruitment), with no process error.

## Usage

``` r
SSS(
  x = 1,
  Data,
  dep = 0.4,
  SR = c("BH", "Ricker"),
  rescale = "mean1",
  start = NULL,
  prior = list(),
  silent = TRUE,
  opt_hess = FALSE,
  n_restart = ifelse(opt_hess, 0, 1),
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  ...
)
```

## Arguments

- x:

  A position in the Data object (by default, equal to one for
  assessments).

- Data:

  An object of class Data

- dep:

  Depletion value to use in the model. Can be an expression that will be
  evaluated inside the function.

- SR:

  Stock-recruit function (either `"BH"` for Beverton-Holt or
  `"Ricker"`).

- rescale:

  A multiplicative factor that rescales the catch in the assessment
  model, which can improve convergence. By default, `"mean1"` scales the
  catch so that time series mean is 1, otherwise a numeric. Output is
  re-converted back to original units.

- start:

  Optional named list of starting values. Entries can be expressions
  that are evaluated in the function:

  - `R0` Unfished recruitment

  - `vul_par` A length-two vector for the age of 95% and 50% fleet
    selectivity. Fixed to maturity otherwise.

- prior:

  A named list for the parameters of any priors to be added to the
  model. See details in `SCA_Pope`.

- silent:

  Logical, passed to
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html),
  whether TMB will print trace information during optimization. Used for
  diagnostics for model convergence.

- opt_hess:

  Logical, whether the hessian function will be passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) during
  optimization (this generally reduces the number of iterations to
  convergence, but is memory and time intensive and does not guarantee
  an increase in convergence rate).

- n_restart:

  The number of restarts (calls to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)) in the
  optimization procedure, so long as the model hasn't converged. The
  optimization continues from the parameters from the previous
  (re)start.

- control:

  A named list of arguments for optimization to be passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- ...:

  Other arguments to be passed (not currently used).

## Value

An object of class
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md).

## Details

In SAMtool, SSS is an implementation of
[SCA_Pope](https://samtool.openmse.com/reference/SCA.md) with fixed
final depletion (in terms of total biomass, not spawning biomass)
assumption.

## References

Cope, J.M. 2013. Implementing a statistical catch-at-age model (Stock
Synthesis) as a tool for deriving overfishing limits in data-limited
situations. Fisheries Research 142:3-14.

## Author

Q. Huynh

## Examples

``` r
res <- SSS(Data = Red_snapper)

SSS_MP <- make_MP(SSS, HCR40_10, dep = 0.3) # Always assume depletion = 0.3
```
