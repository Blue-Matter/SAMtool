# Continuous Delay-differential assessment model

A catch and index-based assessment model. Compared to the discrete
delay-difference (annual time-step in production and fishing), the
delay-differential model (cDD) is based on continuous recruitment and
fishing mortality within a time-step. The continuous model works much
better for populations with high turnover (e.g. high F or M, continuous
reproduction). This model is conditioned on catch and fits to the
observed index. In the state-space version (cDD_SS), recruitment
deviations from the stock-recruit relationship are estimated.

## Usage

``` r
cDD(
  x = 1,
  Data,
  AddInd = "B",
  SR = c("BH", "Ricker"),
  rescale = "mean1",
  MW = FALSE,
  start = NULL,
  prior = list(),
  fix_h = TRUE,
  dep = 1,
  LWT = list(),
  n_itF = 5L,
  silent = TRUE,
  opt_hess = FALSE,
  n_restart = ifelse(opt_hess, 0, 1),
  control = list(iter.max = 5000, eval.max = 10000),
  ...
)

cDD_SS(
  x = 1,
  Data,
  AddInd = "B",
  SR = c("BH", "Ricker"),
  rescale = "mean1",
  MW = FALSE,
  start = NULL,
  prior = list(),
  fix_h = TRUE,
  fix_sigma = FALSE,
  fix_tau = TRUE,
  dep = 1,
  LWT = list(),
  n_itF = 5L,
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
  assessment model.

- dep:

  The initial depletion in the first year of the model. A tight prior is
  placed on the model objective function to estimate the equilibrium
  fishing mortality corresponding to the initial depletion. Due to this
  tight prior, this F should not be considered to be an independent
  model parameter. Set to zero to eliminate this prior.

- LWT:

  A named list of likelihood weights. For `LWT$Index`, a vector of
  likelihood weights for each survey, while for `LWT$MW` a numeric.

- n_itF:

  Integer, the number of iterations to solve F conditional on the
  observed catch.

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

- fix_sigma:

  Logical, whether the standard deviation of the index is fixed. If
  `TRUE`, sigma is fixed to value provided in `start` (if provided),
  otherwise, value based on `Data@CV_Ind`.

- fix_tau:

  Logical, the standard deviation of the recruitment deviations is
  fixed. If `TRUE`, tau is fixed to value provided in `start` (if
  provided), otherwise, equal to 1.

- integrate:

  Logical, whether the likelihood of the model integrates over the
  likelihood of the recruitment deviations (thus, treating it as a
  state-space variable). Otherwise, recruitment deviations are penalized
  parameters.

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

- `Kappa` Delay-differential Kappa parameter. Otherwise, calculated from
  biological parameters in the Data object.

- `F_equilibrium` Equilibrium fishing mortality leading into first year
  of the model (to determine initial depletion). By default, 0.

- `tau` Lognormal SD of the recruitment deviations (process error) for
  `DD_SS`. By default, `Data@sigmaR[x]`.

- `sigma` Lognormal SD of the index (observation error). By default,
  `Data@CV_Ind[x]`. Not used if multiple indices are used.

- `sigma_W` Lognormal SD of the mean weight (observation error). By
  default, 0.1.

Multiple indices are supported in the model. Data@Ind, Data@VInd, and
Data@SpInd are all assumed to be biomass-based. For Data@AddInd,
Data@I_units are used to identify a biomass vs. abundance-based index.

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

- `cDD`: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge

- `cDD_SS`: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge

## Optional Data

- `cDD`: steep

- `cDD_SS`: steep, CV_Ind, sigmaR

## References

Hilborn, R., and Walters, C., 1992. Quantitative Fisheries Stock
Assessment: Choice, Dynamics and Uncertainty. Chapman and Hall, New
York.

## See also

[DD_TMB](https://samtool.openmse.com/reference/DD_TMB.md)
[plot.Assessment](https://samtool.openmse.com/reference/plot.Assessment.md)
[summary.Assessment](https://samtool.openmse.com/reference/summary.Assessment.md)
[retrospective](https://samtool.openmse.com/reference/retrospective.md)
[profile](https://samtool.openmse.com/reference/profile.md)
[make_MP](https://samtool.openmse.com/reference/make_MP.md)

## Author

Q. Huynh

## Examples

``` r
#### Observation-error delay difference model
res <- cDD(Data = MSEtool::Red_snapper)

### State-space version
### Also set recruitment variability SD = 0.6 (since fix_tau = TRUE)
res <- cDD_SS(Data = MSEtool::Red_snapper, start = list(tau = 0.6))

summary(res@SD) # Parameter estimates
#>                 Estimate Std. Error
#> R0x         -0.516462146        NaN
#> log_sigma   -1.825697311        NaN
#> log_rec_dev  0.007309376        NaN
#> log_rec_dev  0.033324174        NaN
#> log_rec_dev  0.036138491        NaN
#> log_rec_dev  0.037280642        NaN
#> log_rec_dev  0.053679032        NaN
#> log_rec_dev  0.113799092        NaN
#> log_rec_dev  0.148088168        NaN
#> log_rec_dev  0.185460723        NaN
#> log_rec_dev  0.242453173        NaN
#> log_rec_dev  0.310530201        NaN
#> log_rec_dev  0.342690648        NaN
#> log_rec_dev  0.337102550        NaN
#> log_rec_dev  0.330495999        NaN
#> log_rec_dev  0.366257963        NaN
#> log_rec_dev  0.384039989        NaN
#> log_rec_dev  0.402561567        NaN
#> log_rec_dev  0.453567855        NaN
#> log_rec_dev  0.414847098        NaN
#> log_rec_dev  0.355096744        NaN
#> log_rec_dev  0.323088402        NaN
#> log_rec_dev  0.318067480        NaN
#> log_rec_dev  0.327920944        NaN
#> log_rec_dev  0.306529182        NaN
#> log_rec_dev  0.302725964        NaN
#> log_rec_dev  0.353198959        NaN
#> log_rec_dev  0.384240430        NaN
#> log_rec_dev  0.446717608        NaN
#> log_rec_dev  0.539499589        NaN
#> log_rec_dev  0.564339576        NaN
#> log_rec_dev  0.511080900        NaN
#> log_rec_dev  0.563475086        NaN
#> log_rec_dev  0.539197361        NaN
#> log_rec_dev  0.367188417        NaN
#> log_rec_dev  0.221219093        NaN
#> log_rec_dev  0.139529613        NaN
#> log_rec_dev  0.085437945        NaN
#> log_rec_dev  0.065626784        NaN
#> log_rec_dev  0.062538212        NaN
#> log_rec_dev  0.074750056        NaN
#> log_rec_dev  0.041520994        NaN
#> log_rec_dev  0.007606233        NaN
#> log_rec_dev  0.007096973        NaN
#> log_rec_dev  0.053798631        NaN
#> log_rec_dev  0.095699895        NaN
#> log_rec_dev  0.092595003        NaN
#> log_rec_dev  0.070754701        NaN
#> log_rec_dev  0.059231729        NaN
#> log_rec_dev  0.050443507        NaN
#> log_rec_dev  0.049333273        NaN
#> log_rec_dev  0.049128428        NaN
#> log_rec_dev  0.031476104        NaN
#> log_rec_dev  0.028465058        NaN
#> log_rec_dev  0.011997883        NaN
#> log_rec_dev -0.014335908        NaN
#> log_rec_dev -0.013501650        NaN
#> log_rec_dev -0.023189784        NaN
#> log_rec_dev -0.013922131        NaN
#> log_rec_dev -0.037682299        NaN
#> log_rec_dev -0.105658047        NaN
#> log_rec_dev -0.150263442        NaN
#> log_rec_dev -0.157945761        NaN
#> log_rec_dev -0.150706928        NaN
#> log_rec_dev -0.167906572        NaN
#> log_rec_dev -0.229036300        NaN
#> log_rec_dev -0.274486497        NaN
#> log_rec_dev -0.279307235        NaN
#> log_rec_dev -0.248849718        NaN
#> log_rec_dev -0.214099416        NaN
#> log_rec_dev -0.205747970        NaN
#> log_rec_dev -0.231667491        NaN
#> log_rec_dev -0.256883800        NaN
#> log_rec_dev -0.273316451        NaN
#> log_rec_dev -0.270855428        NaN
#> log_rec_dev -0.275867215        NaN
#> log_rec_dev -0.270787993        NaN
#> log_rec_dev -0.275187001        NaN
#> log_rec_dev -0.282803279        NaN
#> log_rec_dev -0.249544800        NaN
#> log_rec_dev -0.229031544        NaN
#> log_rec_dev -0.218633177        NaN
#> log_rec_dev -0.186554106        NaN
#> log_rec_dev -0.150406056        NaN
#> log_rec_dev -0.107573490        NaN
#> log_rec_dev -0.063303275        NaN
#> log_rec_dev -0.065666120        NaN
#> log_rec_dev -0.033815268        NaN
#> log_rec_dev -0.025354704        NaN
#> log_rec_dev -0.019843728        NaN
#> log_rec_dev -0.023144074        NaN
#> log_rec_dev -0.030285521        NaN
#> log_rec_dev -0.061928335        NaN
#> log_rec_dev -0.079617299        NaN
#> log_rec_dev -0.115037875        NaN
#> log_rec_dev -0.157590805        NaN
#> log_rec_dev -0.193396401        NaN
#> log_rec_dev -0.211749821        NaN
#> log_rec_dev -0.206408787        NaN
#> log_rec_dev -0.215012407        NaN
#> log_rec_dev -0.202239284        NaN
#> log_rec_dev -0.164449617        NaN
#> log_rec_dev -0.168729177        NaN
#> log_rec_dev -0.166330932        NaN
#> log_rec_dev -0.123512727        NaN
#> log_rec_dev -0.155587811        NaN
#> log_rec_dev -0.179035004        NaN
#> log_rec_dev -0.271359916        NaN
#> log_rec_dev -0.358704133        NaN
#> log_rec_dev -0.356621232        NaN
#> log_rec_dev -0.411378487        NaN
#> log_rec_dev -0.447722953        NaN
#> log_rec_dev -0.482605564        NaN
#> log_rec_dev -0.533069456        NaN
#> log_rec_dev -0.554356514        NaN
#> log_rec_dev -0.608834077        NaN
#> log_rec_dev -0.658653773        NaN
#> log_rec_dev -0.648402316        NaN
#> log_rec_dev -0.732410697        NaN
#> log_rec_dev -0.792387014        NaN
#> log_rec_dev -0.788267328        NaN
#> log_rec_dev -0.714021933        NaN
#> log_rec_dev -0.477556493        NaN
#> log_rec_dev -0.166597079        NaN
#> log_rec_dev  0.113171464        NaN
#> log_rec_dev  0.776095656        NaN
#> log_rec_dev  0.065704176        NaN
#> log_rec_dev  1.475117497        NaN
#> log_rec_dev -0.469887281        NaN
#> log_rec_dev -0.594670030        NaN
#> log_rec_dev  1.412706083        NaN
#> log_rec_dev  0.923574520        NaN
#> log_rec_dev  1.385879522        NaN
#> log_rec_dev  0.175963525        NaN
#> log_rec_dev  0.077785849        NaN
#> log_rec_dev  0.205758440        NaN
#> log_rec_dev  2.070546485        NaN
#> log_rec_dev  0.593703916        NaN
#> log_rec_dev  1.252685884        NaN
#> log_rec_dev  0.482254069        NaN
#> log_rec_dev  0.125041458        NaN
#> log_rec_dev  0.000000000        NaN
#> R0           0.944441461        NaN
#> h            0.990000000        NaN
#> q            0.008799889        NaN
#> sigma        0.161105264        NaN
```
