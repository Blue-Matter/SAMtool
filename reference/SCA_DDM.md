# SCA models with time-varying natural mortality

A modification of `SCA` that incorporates density-dependent effects on M
based on biomass depletion (Forrest et al. 2018). Set the bounds of M in
the `M_bounds` argument, a length-2 vector where the first entry is M0,
the M as B/B0 \>= 1, and the second entry is M1, the M as B/B0
approaches zero. Note that M0 can be greater than M1 (compensatory) or
M0 can be less than M1 (depensatory).

## Usage

``` r
SCA_DDM(
  x = 1,
  Data,
  AddInd = "B",
  SR = c("BH", "Ricker", "none"),
  vulnerability = c("logistic", "dome"),
  catch_eq = c("Baranov", "Pope"),
  CAA_dist = c("multinomial", "lognormal"),
  CAA_multiplier = 50,
  rescale = "mean1",
  max_age = Data@MaxAge,
  start = NULL,
  prior = list(),
  fix_h = TRUE,
  fix_F_equilibrium = TRUE,
  fix_omega = TRUE,
  fix_tau = TRUE,
  LWT = list(),
  early_dev = c("comp_onegen", "comp", "all"),
  late_dev = "comp50",
  M_bounds = NULL,
  integrate = FALSE,
  silent = TRUE,
  opt_hess = FALSE,
  n_restart = ifelse(opt_hess, 0, 1),
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  inner.control = list(),
  ...
)
```

## Arguments

- x:

  A position in the Data object (by default, equal to one for
  assessments).

- Data:

  An object of class Data

- AddInd:

  A vector of integers or character strings indicating the indices to be
  used in the model. Integers assign the index to the corresponding
  index in Data@AddInd, "B" (or 0) represents total biomass in Data@Ind,
  "VB" represents vulnerable biomass in Data@VInd, and "SSB" represents
  spawning stock biomass in Data@SpInd. Vulnerability to the survey is
  fixed in the model.

- SR:

  Stock-recruit function (either `"BH"` for Beverton-Holt, `"Ricker"`,
  or `"none"` for constant mean recruitment).

- vulnerability:

  Whether estimated vulnerability is `"logistic"` or `"dome"`
  (double-normal). See details for parameterization.

- catch_eq:

  Whether to use the Baranov equation or Pope's approximation to
  calculate the predicted catch at age in the model.

- CAA_dist:

  Whether a multinomial or lognormal distribution is used for likelihood
  of the catch-at-age matrix. See details.

- CAA_multiplier:

  Numeric for data weighting of catch-at-age matrix if
  `CAA_hist = "multinomial"`. Otherwise ignored. See details.

- rescale:

  A multiplicative factor that rescales the catch in the assessment
  model, which can improve convergence. By default, `"mean1"` scales the
  catch so that time series mean is 1, otherwise a numeric. Output is
  re-converted back to original units.

- max_age:

  Integer, the maximum age (plus-group) in the model.

- start:

  Optional list of starting values. Entries can be expressions that are
  evaluated in the function. See details.

- prior:

  A named list for the parameters of any priors to be added to the
  model. See below.

- fix_h:

  Logical, whether to fix steepness to value in `Data@steep` in the
  model for `SCA`. This only affects calculation of reference points for
  `SCA2`.

- fix_F_equilibrium:

  Logical, whether the equilibrium fishing mortality prior to the first
  year of the model is estimated. If `TRUE`, `F_equilibrium` is fixed to
  value provided in `start` (if provided), otherwise, equal to zero
  (assumes unfished conditions).

- fix_omega:

  Logical, whether the standard deviation of the catch is fixed. If
  `TRUE`, omega is fixed to value provided in `start` (if provided),
  otherwise, value based on `Data@CV_Cat`.

- fix_tau:

  Logical, the standard deviation of the recruitment deviations is
  fixed. If `TRUE`, tau is fixed to value provided in `start` (if
  provided), otherwise, value based on `Data@sigmaR`.

- LWT:

  A named list (Index, CAA, Catch) of likelihood weights for the data
  components. For the index, a vector of length survey. For CAL and
  Catch, a single value.

- early_dev:

  Numeric or character string describing the years for which recruitment
  deviations are estimated in `SCA`. By default, equal to
  `"comp_onegen"`, where rec devs are estimated one full generation
  prior to the first year when catch-at-age (CAA) data are available.
  With `"comp"`, rec devs are estimated starting in the first year with
  CAA. With `"all"`, rec devs start at the beginning of the model. If
  numeric, the number of years after the first year of the model for
  which to start estimating rec devs. Use negative numbers for years
  prior to the first year.

- late_dev:

  Typically, a numeric for the number of most recent years in which
  recruitment deviations will not be estimated in `SCA` (recruitment in
  these years will be based on the mean predicted by stock-recruit
  relationship). By default, `"comp50"` uses the number of ages (smaller
  than the mode) for which the catch-at-age matrix has less than half
  the abundance than that at the mode.

- M_bounds:

  A numeric vector of length 2 to indicate the M as B/B0 approaches zero
  and one, respectively. By default, set to 75% and 125%, respectively,
  of `Data@Mort[x]`.

- integrate:

  Logical, whether the likelihood of the model integrates over the
  likelihood of the recruitment deviations (thus, treating it as a
  random effects/state-space variable). Otherwise, recruitment
  deviations are penalized parameters.

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
  an increase in convergence rate). Ignored if `integrate = TRUE`.

- n_restart:

  The number of restarts (calls to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)) in the
  optimization procedure, so long as the model hasn't converged. The
  optimization continues from the parameters from the previous
  (re)start.

- control:

  A named list of arguments for optimization to be passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- inner.control:

  A named list of arguments for optimization of the random effects,
  which is passed on to
  [`TMB::newton()`](https://rdrr.io/pkg/TMB/man/newton.html).

- ...:

  Other arguments to be passed, including `yind` (an expression for the
  vector of years to include in the model, useful for debugging for data
  lags), `M_at_age` (set to TRUE to specify a matrix of M by year and
  age from the operating model and the bias parameter), `IAA_hist` (an
  array of index age proportions by year, age, survey), and `IAA_n` (a
  matrix of multinomial sample size by year and survey).

## Value

An object of class
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md).

## Details

See [SCA](https://samtool.openmse.com/reference/SCA.md) for more
information on all arguments.

## Online Documentation

Model description and equations are available on the openMSE
[website](https://openmse.com/features-assessment-models/2-sca/).

## References

Forrest, R.E., Holt, K.R., and Kronlund, A.R. 2018. Performance of
alternative harvest control rules for two Pacific groundfish stocks with
uncertain natural mortality: Bias, robustness and trade-offs. Fisheries
Research 2016: 259-286.

## See also

[SCA](https://samtool.openmse.com/reference/SCA.md)
[SCA_RWM](https://samtool.openmse.com/reference/SCA_RWM.md)
[plot.Assessment](https://samtool.openmse.com/reference/plot.Assessment.md)
[summary.Assessment](https://samtool.openmse.com/reference/summary.Assessment.md)
[retrospective](https://samtool.openmse.com/reference/retrospective.md)
[profile](https://samtool.openmse.com/reference/profile.md)
[make_MP](https://samtool.openmse.com/reference/make_MP.md)

## Author

Q. Huynh

## Examples

``` r
res <- SCA_DDM(Data = MSEtool::SimulatedData)
```
