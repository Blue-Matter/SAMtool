# Delay - Difference Stock Assessment in TMB

A simple delay-difference assessment model using a time-series of
catches and a relative abundance index and coded in TMB. The model can
be conditioned on either (1) effort and estimates predicted catch or (2)
catch and estimates a predicted index. In the state-space version
`DD_SS`, recruitment deviations from the stock-recruit relationship are
estimated.

## Usage

``` r
DD_TMB(
  x = 1,
  Data,
  condition = c("catch", "effort"),
  AddInd = "B",
  SR = c("BH", "Ricker"),
  rescale = "mean1",
  MW = FALSE,
  start = NULL,
  prior = list(),
  fix_h = TRUE,
  dep = 1,
  LWT = list(),
  n_itF = 3L,
  silent = TRUE,
  opt_hess = FALSE,
  n_restart = ifelse(opt_hess, 0, 1),
  control = list(iter.max = 5000, eval.max = 10000),
  ...
)

DD_SS(
  x = 1,
  Data,
  condition = c("catch", "effort"),
  AddInd = "B",
  SR = c("BH", "Ricker"),
  rescale = "mean1",
  MW = FALSE,
  start = NULL,
  prior = list(),
  fix_h = TRUE,
  fix_sd = FALSE,
  fix_tau = TRUE,
  dep = 1,
  LWT = list(),
  n_itF = 3L,
  integrate = FALSE,
  silent = TRUE,
  opt_hess = FALSE,
  n_restart = ifelse(opt_hess, 0, 1),
  control = list(iter.max = 5000, eval.max = 10000),
  inner.control = list(),
  ...
)
```

## Arguments

- x:

  An index for the objects in `Data` when running in closed loop
  simulation. Otherwise, equals to 1 when running an assessment.

- Data:

  An object of class
  [MSEtool::Data](https://msetool.openmse.com/reference/Data-class.html).

- condition:

  A string to indicate whether to condition the model on catch or effort
  (ratio of catch and index).

- AddInd:

  A vector of integers or character strings indicating the indices to be
  used in the model. Integers assign the index to the corresponding
  index in Data@AddInd, "B" (or 0) represents total biomass in Data@Ind,
  "VB" represents vulnerable biomass in Data@VInd, and "SSB" represents
  spawning stock biomass in Data@SpInd.

- SR:

  Stock-recruit function (either `"BH"` for Beverton-Holt or
  `"Ricker"`).

- rescale:

  A multiplicative factor that rescales the catch in the assessment
  model, which can improve convergence. By default, `"mean1"` scales the
  catch so that time series mean is 1, otherwise a numeric. Output is
  re-converted back to original units.

- MW:

  Logical, whether to fit to mean weight. In closed-loop simulation,
  mean weight will be grabbed from `Data@Misc[[x]]$MW`, otherwise
  calculated from `Data@CAL`.

- start:

  Optional list of starting values. Entries can be expressions that are
  evaluated in the function. See details.

- prior:

  A named list for the parameters of any priors to be added to the
  model. See below.

- fix_h:

  Logical, whether to fix steepness to value in `Data@steep` in the
  assessment model. Automatically false if a prior is used.

- dep:

  The initial depletion in the first year of the model. A tight prior is
  placed on the model objective function to estimate the equilibrium
  fishing mortality rate that corresponds to the initial depletion. Due
  to this tight prior, this F should not be considered to be an
  independent model parameter. Set to zero to eliminate this prior.

- LWT:

  A named list of likelihood weights. For `LWT$Index`, a vector of
  likelihood weights for each survey, while for `LWT$MW` a numeric.

- n_itF:

  Integer, the number of iterations to solve F within an annual time
  step when conditioning on catch.

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

  A named list of parameters regarding optimization to be passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- ...:

  Additional arguments (not currently used).

- fix_sd:

  Logical, whether the standard deviation of the data in the likelihood
  (index for conditioning on catch or catch for conditioning on effort).
  If `TRUE`, the SD is fixed to value provided in `start` (if provided),
  otherwise, value based on either `Data@CV_Cat` or `Data@CV_Ind`.

- fix_tau:

  Logical, the standard deviation of the recruitment deviations is
  fixed. If `TRUE`, tau is fixed to value provided in `start` (if
  provided), otherwise, equal to 1.

- integrate:

  Logical, whether the likelihood of the model integrates over the
  likelihood of the recruitment deviations (thus, treating it as a
  random effects/state-space variable). Otherwise, recruitment
  deviations are penalized parameters.

- inner.control:

  A named list of arguments for optimization of the random effects,
  which is passed on to
  [`TMB::newton()`](https://rdrr.io/pkg/TMB/man/newton.html) via
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

## Value

An object of
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
containing objects and output from TMB.

## Details

For `start` (optional), a named list of starting values of estimates can
be provided for:

- `R0` Unfished recruitment. Otherwise, `Data@OM$R0[x]` is used in
  closed-loop, and 400% of mean catch otherwise.

- `h` Steepness. Otherwise, `Data@steep[x]` is used, or 0.9 if empty.

- `M` Natural mortality. Otherwise, `Data@Mort[x]` is used.

- `k` Age of knife-edge maturity. By default, the age of 50% maturity
  calculated from the slots in the Data object.

- `Rho` Delay-difference rho parameter. Otherwise, calculated from
  biological parameters in the Data object.

- `Alpha` Delay-difference alpha parameter. Otherwise, calculated from
  biological parameters in the Data object.

- `q_effort` Scalar coefficient when conditioning on effort (to scale to
  F). Otherwise, 1 is the default.

- `F_equilibrium` Equilibrium fishing mortality rate leading into first
  year of the model (to determine initial depletion). By default, 0.

- `omega` Lognormal SD of the catch (observation error) when
  conditioning on effort. By default, `Data@CV_Cat[x]`.

- `tau` Lognormal SD of the recruitment deviations (process error) for
  `DD_SS`. By default, `Data@sigmaR[x]`.

- `sigma` Lognormal SD of the index (observation error) when
  conditioning on catch. By default, `Data@CV_Ind[x]`. Not used if
  multiple indices are used.

- `sigma_W` Lognormal SD of the mean weight (observation error). By
  default, 0.1.

Multiple indices are supported in the model. Data@Ind, Data@VInd, and
Data@SpInd are all assumed to be biomass-based. For Data@AddInd,
Data@I_units are used to identify a biomass vs. abundance-based index.

Similar to many other assessment models, the model depends on
assumptions such as stationary productivity and proportionality between
the abundance index and real abundance. Unsurprisingly the extent to
which these assumptions are violated tends to be the biggest driver of
performance for this method.

## Priors

The following priors can be added as a named list, e.g.,
`prior = list(M = c(0.25, 0.15), h = c(0.7, 0.1)`. For each parameter
below, provide a vector of values as described:

- `R0` - A vector of length 3. The first value indicates the
  distribution of the prior: `1` for lognormal, `2` for uniform on
  `log(R0)`, `3` for uniform on R0. If lognormal, the second and third
  values are the prior mean (in normal space) and SD (in log space).
  Otherwise, the second and third values are the lower and upper bounds
  of the uniform distribution (values in normal space).

- `h` - A vector of length 2 for the prior mean and SD, both in normal
  space. Beverton-Holt steepness uses a beta distribution, while Ricker
  steepness uses a normal distribution.

- `M` - A vector of length 2 for the prior mean (in normal space) and SD
  (in log space). Lognormal prior.

- `q` - A matrix for nsurvey rows and 2 columns. The first column is the
  prior mean (in normal space) and the second column for the SD (in log
  space). Use `NA` in rows corresponding to indices without priors.

See online documentation for more details.

## Online Documentation

Model description and equations are available on the openMSE
[website](https://openmse.com/features-assessment-models/1-dd/).

## Required Data

- `DD_TMB`: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge

- `DD_SS`: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge

## Optional Data

- `DD_TMB`: steep

- `DD_SS`: steep, CV_Cat

## References

Carruthers, T, Walters, C.J,, and McAllister, M.K. 2012. Evaluating
methods that classify fisheries stock status using only fisheries catch
data. Fisheries Research 119-120:66-79.

Hilborn, R., and Walters, C., 1992. Quantitative Fisheries Stock
Assessment: Choice, Dynamics and Uncertainty. Chapman and Hall, New
York.

## See also

[plot.Assessment](https://samtool.openmse.com/reference/plot.Assessment.md)
[summary.Assessment](https://samtool.openmse.com/reference/summary.Assessment.md)
[retrospective](https://samtool.openmse.com/reference/retrospective.md)
[profile](https://samtool.openmse.com/reference/profile.md)
[make_MP](https://samtool.openmse.com/reference/make_MP.md)

## Author

T. Carruthers & Z. Siders. Zach Siders coded the TMB function.

## Examples

``` r
# \donttest{
#### Observation-error delay difference model
res <- DD_TMB(x = 3, Data = MSEtool::SimulatedData)

# Provide starting values
start <- list(h = 0.95)
res <- DD_TMB(x = 3, Data = MSEtool::SimulatedData, start = start)

summary(res@SD) # Parameter estimates
#>                Estimate   Std. Error
#> R0x       -1.828102e+00 2.164137e-02
#> log_sigma -1.071746e+00 1.414214e-02
#> R0         1.579920e+02 3.419164e+00
#> h          9.500000e-01 0.000000e+00
#> q          1.675884e-04 5.816162e-06
#> sigma      3.424100e-01 4.842410e-03

### State-space version
### Set recruitment variability SD = 0.3 (since fix_tau = TRUE)
res <- DD_SS(x = 3, Data = MSEtool::SimulatedData, start = list(tau = 0.3))
# }
```
