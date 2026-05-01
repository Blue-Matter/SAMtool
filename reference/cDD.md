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
#>                 Estimate  Std. Error
#> R0x         -0.480922690 0.054755037
#> log_sigma   -1.847400156 0.064246874
#> log_rec_dev  0.081585552 0.604988498
#> log_rec_dev  0.111204258 0.612587292
#> log_rec_dev  0.107317828 0.611723512
#> log_rec_dev  0.097498772 0.609381361
#> log_rec_dev  0.107643998 0.612466570
#> log_rec_dev  0.141450359 0.622426174
#> log_rec_dev  0.180772523 0.634807000
#> log_rec_dev  0.222163701 0.648820348
#> log_rec_dev  0.269649807 0.665920925
#> log_rec_dev  0.317761596 0.684575002
#> log_rec_dev  0.348398633 0.697962544
#> log_rec_dev  0.349581087 0.700262933
#> log_rec_dev  0.339558849 0.697642877
#> log_rec_dev  0.370494049 0.711145666
#> log_rec_dev  0.383392479 0.717295328
#> log_rec_dev  0.401659019 0.725163642
#> log_rec_dev  0.443006415 0.742971004
#> log_rec_dev  0.405950433 0.726235883
#> log_rec_dev  0.354623491 0.703767185
#> log_rec_dev  0.322121692 0.689716945
#> log_rec_dev  0.323792447 0.689490538
#> log_rec_dev  0.328406517 0.690987769
#> log_rec_dev  0.315125904 0.686311693
#> log_rec_dev  0.324127270 0.690587929
#> log_rec_dev  0.369391520 0.710248453
#> log_rec_dev  0.408294437 0.729571683
#> log_rec_dev  0.452465677 0.753192651
#> log_rec_dev  0.526305443 0.793804107
#> log_rec_dev  0.529256339 0.797543327
#> log_rec_dev  0.490693518 0.775490054
#> log_rec_dev  0.511397887 0.780147886
#> log_rec_dev  0.482632852 0.757616010
#> log_rec_dev  0.368159238 0.702081430
#> log_rec_dev  0.241149902 0.651837297
#> log_rec_dev  0.158674162 0.622955214
#> log_rec_dev  0.105521401 0.605477839
#> log_rec_dev  0.081409337 0.596944578
#> log_rec_dev  0.068900511 0.592022787
#> log_rec_dev  0.068757180 0.590834220
#> log_rec_dev  0.036591432 0.582972727
#> log_rec_dev  0.007573868 0.576172448
#> log_rec_dev  0.009039320 0.576110670
#> log_rec_dev  0.045867697 0.583877677
#> log_rec_dev  0.076352243 0.590481686
#> log_rec_dev  0.092550847 0.594067806
#> log_rec_dev  0.074280593 0.589868225
#> log_rec_dev  0.059466663 0.586362320
#> log_rec_dev  0.053896760 0.584808590
#> log_rec_dev  0.053596310 0.584330283
#> log_rec_dev  0.043759633 0.581638717
#> log_rec_dev  0.040902411 0.580250216
#> log_rec_dev  0.051550661 0.581636681
#> log_rec_dev  0.024809164 0.575079304
#> log_rec_dev -0.004892577 0.567922602
#> log_rec_dev -0.007052784 0.566307859
#> log_rec_dev -0.020611503 0.562217145
#> log_rec_dev -0.022963039 0.560110364
#> log_rec_dev -0.054165205 0.552655513
#> log_rec_dev -0.109252507 0.541436024
#> log_rec_dev -0.152850736 0.533011784
#> log_rec_dev -0.162743886 0.530233729
#> log_rec_dev -0.158929396 0.529839190
#> log_rec_dev -0.177709635 0.526402538
#> log_rec_dev -0.223378861 0.519037543
#> log_rec_dev -0.247060264 0.514942721
#> log_rec_dev -0.245104123 0.514679408
#> log_rec_dev -0.212018157 0.519092180
#> log_rec_dev -0.198147510 0.520049851
#> log_rec_dev -0.204947991 0.517776095
#> log_rec_dev -0.229863456 0.513387415
#> log_rec_dev -0.266253538 0.507602630
#> log_rec_dev -0.296355848 0.503038559
#> log_rec_dev -0.308545734 0.501128268
#> log_rec_dev -0.313247839 0.500458982
#> log_rec_dev -0.317297845 0.499988480
#> log_rec_dev -0.312030681 0.500571055
#> log_rec_dev -0.294178095 0.502963272
#> log_rec_dev -0.270546210 0.506272599
#> log_rec_dev -0.277473634 0.505658739
#> log_rec_dev -0.269799659 0.507346704
#> log_rec_dev -0.237520149 0.512457301
#> log_rec_dev -0.206456617 0.517525074
#> log_rec_dev -0.150560382 0.526485451
#> log_rec_dev -0.067215236 0.540240287
#> log_rec_dev -0.037339249 0.546476578
#> log_rec_dev -0.033384560 0.548680914
#> log_rec_dev -0.016737194 0.552261849
#> log_rec_dev -0.015008130 0.551574100
#> log_rec_dev  0.008963531 0.553703826
#> log_rec_dev -0.037507878 0.544064999
#> log_rec_dev -0.066610882 0.536932026
#> log_rec_dev -0.127080856 0.525853550
#> log_rec_dev -0.175933456 0.517038478
#> log_rec_dev -0.196005352 0.511556978
#> log_rec_dev -0.237244134 0.503653061
#> log_rec_dev -0.267633509 0.498004798
#> log_rec_dev -0.282640857 0.494973424
#> log_rec_dev -0.311076024 0.491149028
#> log_rec_dev -0.297844026 0.492242403
#> log_rec_dev -0.241677021 0.498511044
#> log_rec_dev -0.196282332 0.504259829
#> log_rec_dev -0.177485964 0.507641928
#> log_rec_dev -0.155903901 0.510567615
#> log_rec_dev -0.128720656 0.510844655
#> log_rec_dev -0.215280362 0.497358657
#> log_rec_dev -0.313313040 0.482515922
#> log_rec_dev -0.391953204 0.471847450
#> log_rec_dev -0.421239214 0.466709621
#> log_rec_dev -0.434542915 0.463308517
#> log_rec_dev -0.488346647 0.456735580
#> log_rec_dev -0.560649877 0.448791617
#> log_rec_dev -0.630095233 0.439891174
#> log_rec_dev -0.668434222 0.432482347
#> log_rec_dev -0.696914473 0.429339633
#> log_rec_dev -0.733073117 0.424591698
#> log_rec_dev -0.795100355 0.418239676
#> log_rec_dev -0.810627831 0.413912168
#> log_rec_dev -0.826024354 0.410720031
#> log_rec_dev -0.829314221 0.410401743
#> log_rec_dev -0.728546468 0.419413509
#> log_rec_dev -0.488680712 0.437934476
#> log_rec_dev -0.146387657 0.457833884
#> log_rec_dev  0.101647830 0.444968852
#> log_rec_dev  0.441377702 0.358842055
#> log_rec_dev -0.049131629 0.493604853
#> log_rec_dev  1.435452211 0.163802698
#> log_rec_dev -0.601919689 0.465119074
#> log_rec_dev -0.693009262 0.456042201
#> log_rec_dev  1.463035339 0.213009828
#> log_rec_dev  0.839240564 0.683712905
#> log_rec_dev  1.327530117 0.413063724
#> log_rec_dev  0.128812802 0.599026300
#> log_rec_dev  0.027014260 0.553999434
#> log_rec_dev  0.165567682 0.602508315
#> log_rec_dev  2.067426165 0.239331720
#> log_rec_dev  0.590296517 0.924061547
#> log_rec_dev  1.198961579 0.960762363
#> log_rec_dev  0.488662732 0.841409737
#> log_rec_dev  0.125368999 0.637116815
#> log_rec_dev  0.000000000 0.600000000
#> R0           0.978609965 0.053583825
#> h            0.990000000 0.000000000
#> q            0.008233668 0.000534874
#> sigma        0.157646490 0.010128294
```
