# Age-structured model using fishery length composition

A single-fleet assessment that fits to catch, indices of abundance, and
fishery length compositions. See
[SCA](https://samtool.openmse.com/reference/SCA.md) for all details.

## Usage

``` r
SCA_CAL(
  x = 1,
  Data,
  AddInd = "B",
  SR = c("BH", "Ricker", "none"),
  vulnerability = c("logistic", "dome"),
  catch_eq = c("Baranov", "Pope"),
  CAL_dist = c("multinomial", "lognormal"),
  CAL_multiplier = 50,
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

- CAL_dist:

  Character, the statistical distribution for the likelihood of the
  catch-at-length.

- CAL_multiplier:

  Numeric for data weighting of catch-at-length matrix if
  `CAL_hist = "multinomial"`. A value smaller than one rescales annual
  sample sizes to this fraction of the original sample size. Values
  greater than one generates a cap of the annual sample size to this
  value.

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

## Online Documentation

Model description and equations are available on the openMSE
[website](https://openmse.com/features-assessment-models/2-sca/).

## Author

Q. Huynh
