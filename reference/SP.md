# Surplus production model with FMSY and MSY as leading parameters

A surplus production model that uses only a time-series of catches and a
relative abundance index and coded in TMB. The base model, `SP`, is
conditioned on catch and estimates a predicted index. Continuous surplus
production and fishing is modeled with sub-annual time steps which
should approximate the behavior of ASPIC (Prager 1994). The Fox model,
`SP_Fox`, fixes BMSY/K = 0.37 (1/e). The state-space version, `SP_SS`
estimates annual deviates in biomass. An option allows for setting a
prior for the intrinsic rate of increase. The function for the `spict`
model (Pedersen and Berg, 2016) is available in
[MSEextra](https://msetool.openmse.com/reference/MSEextra.html).

## Usage

``` r
SP(
  x = 1,
  Data,
  AddInd = "B",
  rescale = "mean1",
  start = NULL,
  prior = list(),
  fix_dep = TRUE,
  fix_n = TRUE,
  LWT = NULL,
  n_seas = 4L,
  n_itF = 3L,
  Euler_Lotka = 0L,
  SR_type = c("BH", "Ricker"),
  silent = TRUE,
  opt_hess = FALSE,
  n_restart = ifelse(opt_hess, 0, 1),
  control = list(iter.max = 5000, eval.max = 10000),
  ...
)

SP_SS(
  x = 1,
  Data,
  AddInd = "B",
  rescale = "mean1",
  start = NULL,
  prior = list(),
  fix_dep = TRUE,
  fix_n = TRUE,
  fix_sigma = TRUE,
  fix_tau = TRUE,
  LWT = NULL,
  early_dev = c("all", "index"),
  n_seas = 4L,
  n_itF = 3L,
  Euler_Lotka = 0L,
  SR_type = c("BH", "Ricker"),
  integrate = FALSE,
  silent = TRUE,
  opt_hess = FALSE,
  n_restart = ifelse(opt_hess, 0, 1),
  control = list(iter.max = 5000, eval.max = 10000),
  inner.control = list(),
  ...
)

SP_Fox(x = 1, Data, ...)
```

## Arguments

- x:

  An index for the objects in `Data` when running in
  [runMSE](https://msetool.openmse.com/reference/runMSE.html).
  Otherwise, equals to 1 When running an assessment interactively.

- Data:

  An object of class Data.

- AddInd:

  A vector of integers or character strings indicating the indices to be
  used in the model. Integers assign the index to the corresponding
  index in Data@AddInd, "B" (or 0) represents total biomass in Data@Ind,
  "VB" represents vulnerable biomass in Data@VInd, and "SSB" represents
  spawning stock biomass in Data@SpInd.

- rescale:

  A multiplicative factor that rescales the catch in the assessment
  model, which can improve convergence. By default, `"mean1"` scales the
  catch so that time series mean is 1, otherwise a numeric. Output is
  re-converted back to original units.

- start:

  Optional list of starting values. Entries can be expressions that are
  evaluated in the function. See details.

- prior:

  A named list for the parameters of any priors to be added to the
  model. See details.

- fix_dep:

  Logical, whether to fix the initial depletion (ratio of biomass to
  carrying capacity in the first year of the model). If `TRUE`, uses the
  value in `start`, otherwise equal to 1 (unfished conditions).

- fix_n:

  Logical, whether to fix the exponent of the production function. If
  `TRUE`, uses the value in `start`, otherwise equal to `n = 2`, where
  the biomass at MSY is half of carrying capacity.

- LWT:

  A vector of likelihood weights for each survey.

- n_seas:

  Integer, the number of seasons in the model for calculating continuous
  surplus production.

- n_itF:

  Integer, the number of iterations to solve F conditional on the
  observed catch given multiple seasons within an annual time step.
  Ignored if `n_seas` = 1.

- Euler_Lotka:

  Integer. If greater than zero, the function will calculate a prior for
  the intrinsic rate of increase to use in the estimation model (in lieu
  of an explicit prior in argument `prior`). The value of this argument
  specifies the number of stochastic samples used to calculate the prior
  SD. See section on priors below.

- SR_type:

  If `use_r_prior = TRUE`, the stock-recruit relationship used to
  calculate the stock-recruit alpha parameter from steepness and
  unfished spawners-per-recruit. Used to develop the r prior.

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

  For `SP_Fox`, additional arguments to pass to `SP`.

- fix_sigma:

  Logical, whether the standard deviation of the index is fixed. If
  `TRUE`, sigma is fixed to value provided in `start` (if provided),
  otherwise, value based on `Data@CV_Ind`.

- fix_tau:

  Logical, the standard deviation of the biomass deviations is fixed. If
  `TRUE`, tau is fixed to value provided in `start` (if provided),
  otherwise, equal to 0.1.

- early_dev:

  Character string describing the years for which biomass deviations are
  estimated in `SP_SS`. By default, deviations are estimated in each
  year of the model (`"all"`), while deviations could also be estimated
  once index data are available (`"index"`).

- integrate:

  Logical, whether the likelihood of the model integrates over the
  likelihood of the biomass deviations (thus, treating it as a
  state-space variable).

- inner.control:

  A named list of arguments for optimization of the random effects,
  which is passed on to
  [newton](https://rdrr.io/pkg/TMB/man/newton.html) via
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

## Value

An object of
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
containing objects and output from TMB.

## Details

For `start` (optional), a named list of starting values of estimates can
be provided for:

- `MSY` Maximum sustainable yield.. Otherwise, 300% of mean catch by
  default.

- `FMSY` Steepness. Otherwise, `Data@Mort[x]` or 0.2 is used.

- `dep` Initial depletion (B/B0) in the first year of the model. By
  default, 1.

- `n` The production function exponent that determines BMSY/B0. By
  default, 2 so that BMSY/B0 = 0.5.

- `sigma` Lognormal SD of the index (observation error). By default,
  0.05. Not used with multiple indices.

- `tau` Lognormal SD of the biomass deviations (process error) in
  `SP_SS`. By default, 0.1.

Multiple indices are supported in the model.

Tip: to create the Fox model (Fox 1970), just fix n = 1. See example.

## Note

The model uses the Fletcher (1978) formulation and is parameterized with
FMSY and MSY as leading parameters. The default conditions assume
unfished conditions in the first year of the time series and a symmetric
production function (n = 2).

## Priors

The following priors can be added as a named list, e.g., prior = list(r
= c(0.25, 0.15), MSY = c(50, 0.1). For each parameter below, provide a
vector of values as described:

- `r` - A vector of length 2 for the lognormal prior mean (normal space)
  and SD (lognormal space).

- `MSY` - A vector of length 2 for the lognormal prior mean (normal
  space) and SD (lognormal space).

In lieu of an explicit r prior provided by the user, set argument
`Euler_Lotka = TRUE` to calculate the prior mean and SD using the
Euler-Lotka method (Equation 15a of McAllister et al. 2001). The
Euler-Lotka method is modified to multiply the left-hand side of
equation 15a by the alpha parameter of the stock-recruit relationship
(Stanley et al. 2009). Natural mortality and steepness are sampled in
order to generate a prior distribution for r. See
`vignette("Surplus_production")` for more details.

## Online Documentation

Model description and equations are available on the openMSE
[website](https://openmse.com/features-assessment-models/3-sp/).

## Required Data

- `SP`: Cat, Ind

- `SP_SS`: Cat, Ind

## Optional Data

`SP_SS`: CV_Ind

## References

Fletcher, R. I. 1978. On the restructuring of the Pella-Tomlinson
system. Fishery Bulletin 76:515:521.

Fox, W.W. 1970. An exponential surplus-yield model for optimizing
exploited fish populations. Transactions of the American Fisheries
Society 99:80-88.

McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
demographic methods to construct Bayesian priors for the intrinsic rate
of increase in the Schaefer model and implications for stock rebuilding.
Can. J. Fish. Aquat. Sci. 58: 1871-1890.

Pedersen, M. W. and Berg, C. W. 2017. A stochastic surplus production
model in continuous time. Fish and Fisheries. 18:226-243.

Pella, J. J. and Tomlinson, P. K. 1969. A generalized stock production
model. Inter-Am. Trop. Tuna Comm., Bull. 13:419-496.

Prager, M. H. 1994. A suite of extensions to a nonequilibrium
surplus-production model. Fishery Bulletin 92:374-389.

Stanley, R.D., M. McAllister, P. Starr and N. Olsen. 2009. Stock
assessment for bocaccio (Sebastes paucispinis) in British Columbia
waters. DFO Can. Sci. Advis. Sec. Res. Doc. 2009/055. xiv + 200 p.

## See also

[SP_production](https://samtool.openmse.com/reference/SP_production.md)
[plot.Assessment](https://samtool.openmse.com/reference/plot.Assessment.md)
[summary.Assessment](https://samtool.openmse.com/reference/summary.Assessment.md)
[retrospective](https://samtool.openmse.com/reference/retrospective.md)
[profile](https://samtool.openmse.com/reference/profile.md)
[make_MP](https://samtool.openmse.com/reference/make_MP.md)

## Author

Q. Huynh

## Examples

``` r
data(swordfish)

#### Observation-error surplus production model
res <- SP(Data = swordfish)

# Provide starting values, assume B/K = 0.875 in first year of model
# and symmetrical production curve (n = 2)
start <- list(dep = 0.875, n = 2)
res <- SP(Data = swordfish, start = start)

# \donttest{
plot(res)
#> ℹ Writing markdown file: /tmp/RtmpRktVJQ/report_SP.Rmd
#> ℹ Rendering markdown file: /tmp/RtmpRktVJQ/report_SP.Rmd
#> ℹ See help(plot.Assessment) to adjust report and file directory.
#> ✔ Rendered file: /tmp/RtmpRktVJQ/report_SP.html
profile(res, FMSY = seq(0.1, 0.4, 0.01))

#> An object of class "prof"
#> Slot "Model":
#> [1] "SP"
#> 
#> Slot "Name":
#> [1] "North Atlantic Swordfish (Source: ASPIC software)"
#> 
#> Slot "Par":
#> [1] "FMSY"
#> 
#> Slot "MLE":
#> [1] 0.2779916
#> 
#> Slot "grid":
#>    FMSY      MSY          nll
#> 1  0.10 14279.16 3.229494e+00
#> 2  0.11 14279.16 2.870962e+00
#> 3  0.12 14279.16 2.530610e+00
#> 4  0.13 14279.16 2.209910e+00
#> 5  0.14 14279.16 1.909978e+00
#> 6  0.15 14279.16 1.631634e+00
#> 7  0.16 14279.16 1.375443e+00
#> 8  0.17 14279.16 1.141752e+00
#> 9  0.18 14279.16 9.307120e-01
#> 10 0.19 14279.16 7.423049e-01
#> 11 0.20 14279.16 5.763609e-01
#> 12 0.21 14279.16 4.325774e-01
#> 13 0.22 14279.16 3.105355e-01
#> 14 0.23 14279.16 2.097149e-01
#> 15 0.24 14279.16 1.295084e-01
#> 16 0.25 14279.16 6.923471e-02
#> 17 0.26 14279.16 2.815079e-02
#> 18 0.27 14279.16 5.463285e-03
#> 19 0.28 14279.16 3.392047e-04
#> 20 0.29 14279.16 1.191583e-02
#> 21 0.30 14279.16 3.930994e-02
#> 22 0.31 14279.16 8.162629e-02
#> 23 0.32 14279.16 1.999949e+01
#> 24 0.33 14279.16 2.074303e-01
#> 25 0.34 14279.16 5.523536e+00
#> 26 0.35 14279.16 3.822012e-01
#> 27 0.36 14279.16 2.003749e+01
#> 28 0.37 14279.16 2.004625e+01
#> 29 0.38 14279.16 2.005471e+01
#> 30 0.39 14279.16 1.466457e+01
#> 31 0.40 14279.16 2.007081e+01
#> 
retrospective(res)





#>                   Mohn's rho
#> Fishing mortality      0.129
#> F/F[MSY]              -0.145
#> Biomass               -0.091
#> B/B[MSY]               0.152
#> B/B[0]                 0.152
# }

#### State-space version
res_SS <- SP_SS(Data = swordfish, start = list(dep = 0.875, sigma = 0.1, tau = 0.1))

# \donttest{
plot(res_SS)
#> ℹ Writing markdown file: /tmp/RtmpRktVJQ/report_SP_SS.Rmd
#> ℹ Rendering markdown file: /tmp/RtmpRktVJQ/report_SP_SS.Rmd
#> ℹ See help(plot.Assessment) to adjust report and file directory.
#> ✔ Rendered file: /tmp/RtmpRktVJQ/report_SP_SS.html
# }

#### Fox model
res_Fox <- SP(Data = swordfish, start = list(n = 1), fix_n = TRUE)
res_Fox2 <- SP_Fox(Data = swordfish)

#### SP with r prior calculated internally (100 stochastic samples to get prior SD)
res_prior <- SP(Data = SimulatedData, Euler_Lotka = 100)

#### Pass an r prior to the model with mean = 0.35, lognormal sd = 0.10
res_prior2 <- SP(Data = SimulatedData, prior = list(r = c(0.35, 0.10)))

#### Pass MSY prior to the model with mean = 1500, lognormal sd = 0.05
res_prior3 <- SP(Data = SimulatedData, prior = list(MSY = c(1500, 0.05)))
```
