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
#>                 Estimate   Std. Error
#> R0x         -0.480921977 0.0547469283
#> log_sigma   -1.847399768 0.0642471864
#> log_rec_dev  0.081585070 0.6049890134
#> log_rec_dev  0.111204459 0.6125896948
#> log_rec_dev  0.107320301 0.6117276030
#> log_rec_dev  0.097502590 0.6093867300
#> log_rec_dev  0.107642611 0.6124732877
#> log_rec_dev  0.141448991 0.6224345033
#> log_rec_dev  0.180771699 0.6348172246
#> log_rec_dev  0.222164672 0.6488328380
#> log_rec_dev  0.269644486 0.6659360832
#> log_rec_dev  0.317762970 0.6845939442
#> log_rec_dev  0.348398947 0.6979846928
#> log_rec_dev  0.349580453 0.7002867979
#> log_rec_dev  0.339556165 0.6976676386
#> log_rec_dev  0.370496407 0.7111749853
#> log_rec_dev  0.383391062 0.7173274611
#> log_rec_dev  0.401657087 0.7251996545
#> log_rec_dev  0.443006160 0.7430141582
#> log_rec_dev  0.405949604 0.7262769119
#> log_rec_dev  0.354624978 0.7038049429
#> log_rec_dev  0.322123690 0.6897534938
#> log_rec_dev  0.323791520 0.6895292464
#> log_rec_dev  0.328405697 0.6910295135
#> log_rec_dev  0.315129195 0.6863550363
#> log_rec_dev  0.324127285 0.6906353866
#> log_rec_dev  0.369394827 0.7103066557
#> log_rec_dev  0.408292392 0.7296413384
#> log_rec_dev  0.452470791 0.7532794228
#> log_rec_dev  0.526303164 0.7939200118
#> log_rec_dev  0.529253552 0.7976696913
#> log_rec_dev  0.490694033 0.7756108253
#> log_rec_dev  0.511394181 0.7802830699
#> log_rec_dev  0.482631971 0.7577464578
#> log_rec_dev  0.368157177 0.7021809328
#> log_rec_dev  0.241152862 0.6519127298
#> log_rec_dev  0.158676910 0.6230205089
#> log_rec_dev  0.105521755 0.6055393119
#> log_rec_dev  0.081408084 0.5970068904
#> log_rec_dev  0.068900324 0.5920876627
#> log_rec_dev  0.068758105 0.5909036485
#> log_rec_dev  0.036589901 0.5830421567
#> log_rec_dev  0.007572692 0.5762423034
#> log_rec_dev  0.009040654 0.5761856474
#> log_rec_dev  0.045865479 0.5839648984
#> log_rec_dev  0.076351132 0.5905813875
#> log_rec_dev  0.092545044 0.5941784237
#> log_rec_dev  0.074278715 0.5899813333
#> log_rec_dev  0.059469496 0.5864783964
#> log_rec_dev  0.053894613 0.5849304078
#> log_rec_dev  0.053599475 0.5844590360
#> log_rec_dev  0.043765636 0.5817722050
#> log_rec_dev  0.040903529 0.5803912726
#> log_rec_dev  0.051544370 0.5817906102
#> log_rec_dev  0.024806508 0.5752340260
#> log_rec_dev -0.004889809 0.5680770979
#> log_rec_dev -0.007056536 0.5664722907
#> log_rec_dev -0.020612409 0.5623878474
#> log_rec_dev -0.022961763 0.5602921459
#> log_rec_dev -0.054161330 0.5528385116
#> log_rec_dev -0.109253247 0.5416116409
#> log_rec_dev -0.152853818 0.5331830037
#> log_rec_dev -0.162744742 0.5304118986
#> log_rec_dev -0.158920497 0.5300289360
#> log_rec_dev -0.177715656 0.5265977789
#> log_rec_dev -0.223381135 0.5192267185
#> log_rec_dev -0.247064846 0.5151353057
#> log_rec_dev -0.245102022 0.5148845035
#> log_rec_dev -0.212012651 0.5193251935
#> log_rec_dev -0.198147789 0.5203086227
#> log_rec_dev -0.204945697 0.5180512683
#> log_rec_dev -0.229859020 0.5136683298
#> log_rec_dev -0.266243970 0.5078817272
#> log_rec_dev -0.296364464 0.5033193424
#> log_rec_dev -0.308549771 0.5014166639
#> log_rec_dev -0.313257897 0.5007592226
#> log_rec_dev -0.317281189 0.5002980196
#> log_rec_dev -0.312028428 0.5009010428
#> log_rec_dev -0.294194043 0.5033231183
#> log_rec_dev -0.270553218 0.5066674393
#> log_rec_dev -0.277470186 0.5060675851
#> log_rec_dev -0.269795717 0.5077810115
#> log_rec_dev -0.237517495 0.5129418042
#> log_rec_dev -0.206445748 0.5180627102
#> log_rec_dev -0.150568850 0.5271094564
#> log_rec_dev -0.067232387 0.5409977420
#> log_rec_dev -0.037356343 0.5473062479
#> log_rec_dev -0.033371112 0.5495362372
#> log_rec_dev -0.016728778 0.5531673048
#> log_rec_dev -0.015000225 0.5525047551
#> log_rec_dev  0.008968312 0.5546974320
#> log_rec_dev -0.037513481 0.5449921520
#> log_rec_dev -0.066624856 0.5378218180
#> log_rec_dev -0.127094550 0.5266562305
#> log_rec_dev -0.175935789 0.5177771594
#> log_rec_dev -0.195989017 0.5122797796
#> log_rec_dev -0.237236537 0.5043396956
#> log_rec_dev -0.267611927 0.4986674074
#> log_rec_dev -0.282639257 0.4956353776
#> log_rec_dev -0.311088342 0.4917940502
#> log_rec_dev -0.297865462 0.4929142260
#> log_rec_dev -0.241687628 0.4992563033
#> log_rec_dev -0.196259646 0.5050644519
#> log_rec_dev -0.177485825 0.5084715775
#> log_rec_dev -0.155908234 0.5114189655
#> log_rec_dev -0.128707615 0.5117234144
#> log_rec_dev -0.215310418 0.4981221900
#> log_rec_dev -0.313304951 0.4831539732
#> log_rec_dev -0.391948348 0.4723955630
#> log_rec_dev -0.421249151 0.4672175147
#> log_rec_dev -0.434554321 0.4637887831
#> log_rec_dev -0.488333279 0.4571589866
#> log_rec_dev -0.560629841 0.4491554235
#> log_rec_dev -0.630111978 0.4402148163
#> log_rec_dev -0.668453474 0.4327777918
#> log_rec_dev -0.696920625 0.4296070486
#> log_rec_dev -0.733067064 0.4248328389
#> log_rec_dev -0.795098689 0.4184514191
#> log_rec_dev -0.810627817 0.4141120997
#> log_rec_dev -0.826021444 0.4109087374
#> log_rec_dev -0.829308739 0.4105798155
#> log_rec_dev -0.728542734 0.4196009803
#> log_rec_dev -0.488689593 0.4381690277
#> log_rec_dev -0.146390874 0.4581568685
#> log_rec_dev  0.101649462 0.4453275698
#> log_rec_dev  0.441372596 0.3590932125
#> log_rec_dev -0.049117429 0.4938668468
#> log_rec_dev  1.435446867 0.1638065774
#> log_rec_dev -0.601908343 0.4651834714
#> log_rec_dev -0.693002663 0.4560710622
#> log_rec_dev  1.463035045 0.2130146328
#> log_rec_dev  0.839238495 0.6837686624
#> log_rec_dev  1.327529994 0.4130839066
#> log_rec_dev  0.128809363 0.5990381429
#> log_rec_dev  0.027016017 0.5540030390
#> log_rec_dev  0.165569813 0.6025117387
#> log_rec_dev  2.067425294 0.2393331131
#> log_rec_dev  0.590287995 0.9240651438
#> log_rec_dev  1.198969896 0.9607615899
#> log_rec_dev  0.488657328 0.8414100549
#> log_rec_dev  0.125367054 0.6371168564
#> log_rec_dev  0.000000000 0.6000000000
#> R0           0.978610663 0.0535759278
#> h            0.990000000 0.0000000000
#> q            0.008233661 0.0005346941
#> sigma        0.157646551 0.0101283474
```
